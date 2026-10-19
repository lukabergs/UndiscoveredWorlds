#include "io/run_bundle.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>

namespace runbundle
{
std::uint64_t hashbytes(const void* data, std::size_t size)
{
    std::uint64_t hash = 1469598103934665603ULL;
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; index++)
    {
        hash ^= bytes[index];
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string jsonstring(std::string_view text)
{
    std::ostringstream out;
    out << '"';
    for (const char character : text)
    {
        switch (character)
        {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (static_cast<unsigned char>(character) < 0x20)
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(character) << std::dec;
            else
                out << character;
        }
    }
    out << '"';
    return out.str();
}

namespace
{
std::string hexhash(std::uint64_t hash)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

std::string shapetext(const std::vector<std::size_t>& shape, bool json)
{
    std::ostringstream out;
    out << (json ? '[' : '(');
    for (std::size_t index = 0; index < shape.size(); index++)
    {
        if (index)
            out << ", ";
        out << shape[index];
    }
    if (!json && shape.size() == 1)
        out << ',';
    out << (json ? ']' : ')');
    return out.str();
}

std::string utcnow()
{
    const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm utc{};
    gmtime_s(&utc, &now);
    std::ostringstream out;
    out << std::put_time(&utc, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

void writetext(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    out << text;
    if (!out)
        throw std::runtime_error("Failed to write " + path.string());
}
}

void writenpy(const std::filesystem::path& path, const std::string& dtype,
    const std::vector<std::size_t>& shape, const void* data, std::size_t bytes)
{
    std::string header = "{'descr': '" + dtype + "', 'fortran_order': False, 'shape': " + shapetext(shape, false) + ", }";
    // Version 1.0: magic (6) + version (2) + length (2) + header, padded to 64 bytes.
    const std::size_t prefix = 10;
    const std::size_t total = ((prefix + header.size() + 1 + 63) / 64) * 64;
    header.append(total - prefix - header.size() - 1, ' ');
    header += '\n';

    std::ofstream out(path, std::ios::binary);
    out.write("\x93NUMPY", 6);
    const unsigned char version[2] = { 1, 0 };
    out.write(reinterpret_cast<const char*>(version), 2);
    const auto length = static_cast<std::uint16_t>(header.size());
    const unsigned char lengthbytes[2] = { static_cast<unsigned char>(length & 0xff), static_cast<unsigned char>(length >> 8) };
    out.write(reinterpret_cast<const char*>(lengthbytes), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    if (bytes)
        out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out)
        throw std::runtime_error("Failed to write " + path.string());
}

Writer::Writer(std::filesystem::path directory) : directory_(std::move(directory))
{
    if (std::filesystem::exists(directory_) && !std::filesystem::is_empty(directory_))
        throw std::runtime_error("Run bundle directory is not empty: " + directory_.string());
    std::filesystem::create_directories(directory_ / "fields");
    std::filesystem::create_directories(directory_ / "tables");
}

void Writer::addfield(const FieldEntry& field)
{
    const std::string relative = "fields/" + field.component + "/" + field.id.substr(field.id.find('.') + 1) + ".npy";
    std::filesystem::create_directories((directory_ / relative).parent_path());
    writenpy(directory_ / relative, field.dtype, field.shape, field.data, field.bytes);
    WrittenField written{ field, relative, hashbytes(field.data, field.bytes) };
    written.entry.data = nullptr;
    fields_.push_back(std::move(written));
}

void Writer::addexport(const FieldEntry& field, const std::string& rule, bool recomputed)
{
    const std::string relative = "latlon/" + field.component + "/" + field.id.substr(field.id.find('.') + 1) + ".npy";
    std::filesystem::create_directories((directory_ / relative).parent_path());
    writenpy(directory_ / relative, field.dtype, field.shape, field.data, field.bytes);
    WrittenExport written{ field, relative, hashbytes(field.data, field.bytes), rule, recomputed };
    written.entry.data = nullptr;
    exports_.push_back(std::move(written));
}

void Writer::addcells(const ArrayEntry& cells)
{
    const std::string relative = "grid/" + cells.id + ".npy";
    std::filesystem::create_directories((directory_ / relative).parent_path());
    writenpy(directory_ / relative, cells.dtype, cells.shape, cells.data, cells.bytes);
    WrittenArray written{ cells, relative, hashbytes(cells.data, cells.bytes) };
    written.entry.data = nullptr;
    cells_.push_back(std::move(written));
}

void Writer::addtable(const TableEntry& table)
{
    const std::string relative = "tables/" + table.component + "/" + table.id.substr(table.id.find('.') + 1) + table.extension;
    std::filesystem::create_directories((directory_ / relative).parent_path());
    writetext(directory_ / relative, table.text);
    WrittenTable written{ table, relative, hashbytes(table.text.data(), table.text.size()) };
    written.entry.text.clear();
    tables_.push_back(std::move(written));
}

void Writer::finish(const Metadata& metadata)
{
    std::ostringstream out;
    out << std::setprecision(6) << std::fixed;
    out << "{\n";
    out << "  \"schema\": " << jsonstring(kSchema) << ",\n";
    out << "  \"created\": " << jsonstring(utcnow()) << ",\n";
    out << "  \"seed\": " << metadata.seed << ",\n";
    const GridSpec& grid = metadata.grid;
    if (grid.kind == "latlon")
    {
        if (!cells_.empty())
            throw std::logic_error("A lat-lon bundle has no cell table");
        out << "  \"grid\": {\"kind\": \"latlon\", \"columns\": " << metadata.columns << ", \"rows\": " << metadata.rows
            << ", \"registration\": \"cell-centred\", \"layout\": \"row-major, north row first\"";
    }
    else if (grid.kind == "cube")
    {
        if (grid.facecells < 1 || metadata.columns != 4 * grid.facecells || metadata.rows * 2 != metadata.columns)
            throw std::logic_error("A cube bundle needs N >= 1 and its lat-lon export of 4 N x 2 N");
        out << "  \"grid\": {\"kind\": \"cube\", \"face_cells\": " << grid.facecells
            << ", \"orientation\": \"uw-cube-1\", \"layout\": \"face, row (top first), column\"";
    }
    else if (grid.kind == "cube-window")
    {
        const bool ids = std::any_of(cells_.begin(), cells_.end(), [](const WrittenArray& array) { return array.entry.id == "cell_ids"; });
        if (grid.facecells < 1 || grid.basefacecells < 1 || grid.level < 0 || grid.facecells % grid.basefacecells != 0 || !ids ||
            metadata.columns < 1 || metadata.rows < 1)
            throw std::logic_error("A cube-window bundle needs N_L a multiple of N >= 1, a level, its cell_ids and a lat-lon export box");
        out << "  \"grid\": {\"kind\": \"cube-window\", \"face_cells\": " << grid.facecells << ", \"level\": " << grid.level
            << ", \"base_face_cells\": " << grid.basefacecells
            << ", \"orientation\": \"uw-cube-1\", \"layout\": \"cell list, ascending id\", \"cell_count\": " << grid.cellcount;
    }
    else
        throw std::logic_error("Unsupported bundle grid kind: " + grid.kind);
    if (grid.radiuskm > 0.0)
        out << ", \"radius_km\": " << grid.radiuskm;
    if (grid.kind == "cube")
        out << ", \"export_columns\": " << metadata.columns;
    if (grid.kind == "cube-window")
        out << ", \"export_columns\": " << metadata.columns << ", \"export_rows\": " << metadata.rows;
    out << "},\n";
    out << "  \"recipe\": " << jsonstring(metadata.recipe) << ",\n";
    out << "  \"command\": [";
    for (std::size_t index = 0; index < metadata.command.size(); index++)
        out << (index ? ", " : "") << jsonstring(metadata.command[index]);
    out << "],\n";
    out << "  \"properties\": {";
    for (std::size_t index = 0; index < metadata.properties.size(); index++)
        out << (index ? ", " : "") << jsonstring(metadata.properties[index].first) << ": " << jsonstring(metadata.properties[index].second);
    out << "},\n";
    if (!metadata.parameters.empty())
        out << "  \"parameters\": " << metadata.parameters << ",\n";
    for (const auto& [name, value] : metadata.blocks)
        out << "  " << jsonstring(name) << ": " << value << ",\n";
    out << "  \"steps\": [";
    for (std::size_t index = 0; index < metadata.steps.size(); index++)
    {
        const StepRecord& step = metadata.steps[index];
        out << (index ? ",\n" : "\n") << "    {\"index\": " << step.index << ", \"process\": " << jsonstring(step.process)
            << ", \"label\": " << jsonstring(step.label);
        if (!step.overrides.empty())
        {
            out << ", \"overrides\": {";
            for (std::size_t entry = 0; entry < step.overrides.size(); entry++)
                out << (entry ? ", " : "") << jsonstring(step.overrides[entry].first) << ": " << step.overrides[entry].second;
            out << "}";
        }
        out << ", \"seconds\": " << step.seconds << ", \"writes\": [";
        for (std::size_t write = 0; write < step.writes.size(); write++)
            out << (write ? ", " : "") << jsonstring(step.writes[write]);
        out << "]}";
    }
    out << (metadata.steps.empty() ? "],\n" : "\n  ],\n");
    out << "  \"fields\": [";
    for (std::size_t index = 0; index < fields_.size(); index++)
    {
        const FieldEntry& field = fields_[index].entry;
        out << (index ? ",\n" : "\n") << "    {\"id\": " << jsonstring(field.id) << ", \"component\": " << jsonstring(field.component)
            << ", \"layout\": " << jsonstring(field.layout) << ", \"dtype\": " << jsonstring(field.dtype)
            << ", \"shape\": " << shapetext(field.shape, true) << ", \"units\": " << jsonstring(field.units)
            << ", \"file\": " << jsonstring(fields_[index].file) << ", \"hash\": " << jsonstring(hexhash(fields_[index].hash))
            << ", \"written_by\": " << jsonstring(field.writtenby) << ", \"description\": " << jsonstring(field.description) << "}";
    }
    out << (fields_.empty() ? "],\n" : "\n  ],\n");
    if (!exports_.empty())
    {
        out << "  \"exports\": [";
        for (std::size_t index = 0; index < exports_.size(); index++)
        {
            const WrittenExport& written = exports_[index];
            out << (index ? ",\n" : "\n") << "    {\"id\": " << jsonstring(written.entry.id) << ", \"dtype\": " << jsonstring(written.entry.dtype)
                << ", \"shape\": " << shapetext(written.entry.shape, true) << ", \"file\": " << jsonstring(written.file)
                << ", \"hash\": " << jsonstring(hexhash(written.hash)) << ", \"rule\": " << jsonstring(written.rule);
            // Every direction field records whether it was recomputed; another field only when it was.
            if (written.rule == "direction" || written.recomputed)
                out << ", \"recomputed\": " << (written.recomputed ? "true" : "false");
            out << "}";
        }
        out << "\n  ],\n";
    }
    if (!cells_.empty())
    {
        out << "  \"cells\": [";
        for (std::size_t index = 0; index < cells_.size(); index++)
        {
            const WrittenArray& written = cells_[index];
            out << (index ? ",\n" : "\n") << "    {\"id\": " << jsonstring(written.entry.id) << ", \"dtype\": " << jsonstring(written.entry.dtype)
                << ", \"shape\": " << shapetext(written.entry.shape, true) << ", \"file\": " << jsonstring(written.file)
                << ", \"hash\": " << jsonstring(hexhash(written.hash)) << "}";
        }
        out << "\n  ],\n";
    }
    out << "  \"tables\": [";
    for (std::size_t index = 0; index < tables_.size(); index++)
    {
        const TableEntry& table = tables_[index].entry;
        out << (index ? ",\n" : "\n") << "    {\"id\": " << jsonstring(table.id) << ", \"component\": " << jsonstring(table.component)
            << ", \"file\": " << jsonstring(tables_[index].file) << ", \"hash\": " << jsonstring(hexhash(tables_[index].hash))
            << ", \"written_by\": " << jsonstring(table.writtenby) << ", \"description\": " << jsonstring(table.description) << "}";
    }
    out << (tables_.empty() ? "]\n" : "\n  ]\n");
    out << "}\n";
    writetext(directory_ / "manifest.json", out.str());
}

namespace
{
// The JSON a manifest uses, parsed into a tree.
struct JsonValue
{
    double number = 0.0;
    bool boolean = false;
    std::string text;
    std::vector<JsonValue> items;
    std::map<std::string, JsonValue> members;

    const JsonValue& at(const std::string& key) const
    {
        static const JsonValue missing;
        const auto found = members.find(key);
        return found == members.end() ? missing : found->second;
    }
};

class JsonParser
{
public:
    explicit JsonParser(const std::string& text) : text_(text) {}

    JsonValue parse()
    {
        JsonValue value = parsevalue();
        skipspace();
        if (position_ != text_.size())
            fail();
        return value;
    }

private:
    [[noreturn]] void fail() const
    {
        throw std::runtime_error("malformed JSON at offset " + std::to_string(position_));
    }

    void skipspace()
    {
        while (position_ < text_.size() && std::string_view(" \n\r\t").find(text_[position_]) != std::string_view::npos)
            position_++;
    }

    bool consume(char character)
    {
        skipspace();
        if (position_ >= text_.size() || text_[position_] != character)
            return false;
        position_++;
        return true;
    }

    std::string parsestring()
    {
        if (!consume('"'))
            fail();
        std::string out;
        while (position_ < text_.size() && text_[position_] != '"')
        {
            char character = text_[position_++];
            if (character == '\\')
            {
                if (position_ >= text_.size())
                    fail();
                const char escaped = text_[position_++];
                switch (escaped)
                {
                case 'n': character = '\n'; break;
                case 'r': character = '\r'; break;
                case 't': character = '\t'; break;
                case 'u':
                    // jsonstring escapes only control characters this way.
                    if (position_ + 4 > text_.size())
                        fail();
                    character = static_cast<char>(std::strtol(text_.substr(position_, 4).c_str(), nullptr, 16));
                    position_ += 4;
                    break;
                default: character = escaped; break;
                }
            }
            out += character;
        }
        if (position_ >= text_.size())
            fail();
        position_++;
        return out;
    }

    JsonValue parsevalue()
    {
        skipspace();
        if (position_ >= text_.size())
            fail();
        JsonValue value;
        const char first = text_[position_];
        if (first == '{')
        {
            position_++;
            if (consume('}'))
                return value;
            do
            {
                std::string key = parsestring();
                if (!consume(':'))
                    fail();
                value.members[std::move(key)] = parsevalue();
            } while (consume(','));
            if (!consume('}'))
                fail();
        }
        else if (first == '[')
        {
            position_++;
            if (consume(']'))
                return value;
            do
                value.items.push_back(parsevalue());
            while (consume(','));
            if (!consume(']'))
                fail();
        }
        else if (first == '"')
            value.text = parsestring();
        else if (text_.compare(position_, 4, "true") == 0)
        {
            value.boolean = true;
            position_ += 4;
        }
        else if (text_.compare(position_, 4, "null") == 0)
            position_ += 4;
        else if (text_.compare(position_, 5, "false") == 0)
            position_ += 5;
        else
        {
            char* end = nullptr;
            value.number = std::strtod(text_.c_str() + position_, &end);
            if (end == text_.c_str() + position_)
                fail();
            position_ = static_cast<std::size_t>(end - text_.c_str());
        }
        return value;
    }

    const std::string& text_;
    std::size_t position_ = 0;
};

std::string readtext(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Failed to read " + path.string());
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

ManifestEntry manifestentry(const JsonValue& value)
{
    ManifestEntry entry;
    entry.id = value.at("id").text;
    entry.component = value.at("component").text;
    entry.layout = value.at("layout").text;
    entry.dtype = value.at("dtype").text;
    for (const JsonValue& extent : value.at("shape").items)
        entry.shape.push_back(static_cast<std::size_t>(extent.number));
    entry.units = value.at("units").text;
    entry.description = value.at("description").text;
    entry.file = value.at("file").text;
    entry.writtenby = value.at("written_by").text;
    entry.rule = value.at("rule").text;
    entry.recomputed = value.at("recomputed").boolean;
    return entry;
}
}

Reader::Reader(std::filesystem::path directory) : directory_(std::move(directory))
{
    const std::filesystem::path manifestpath = directory_ / "manifest.json";
    const std::string text = readtext(manifestpath);
    JsonValue manifest;
    try
    {
        manifest = JsonParser(text).parse();
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error(manifestpath.string() + ": " + error.what());
    }
    schema_ = manifest.at("schema").text;
    if (schema_ != kSchema && schema_ != kSchemaV1)
        throw std::runtime_error(manifestpath.string() + ": not a " + kSchema + " or " + kSchemaV1 + " manifest");
    seed_ = static_cast<long>(manifest.at("seed").number);
    const JsonValue& grid = manifest.at("grid");
    // v1 has no kind: lat-lon.
    grid_.kind = schema_ == kSchemaV1 || grid.at("kind").text.empty() ? "latlon" : grid.at("kind").text;
    grid_.radiuskm = grid.at("radius_km").number;
    if (grid_.kind == "latlon")
    {
        columns_ = static_cast<int>(grid.at("columns").number);
        rows_ = static_cast<int>(grid.at("rows").number);
    }
    else if (grid_.kind == "cube")
    {
        if (grid.at("orientation").text != "uw-cube-1")
            throw std::runtime_error(manifestpath.string() + ": unknown cube orientation " + grid.at("orientation").text);
        grid_.facecells = static_cast<int>(grid.at("face_cells").number);
        columns_ = static_cast<int>(grid.at("export_columns").number);
        rows_ = columns_ / 2;
    }
    else if (grid_.kind == "cube-window")
    {
        if (grid.at("orientation").text != "uw-cube-1")
            throw std::runtime_error(manifestpath.string() + ": unknown cube orientation " + grid.at("orientation").text);
        grid_.facecells = static_cast<int>(grid.at("face_cells").number);
        grid_.level = static_cast<int>(grid.at("level").number);
        grid_.basefacecells = static_cast<int>(grid.at("base_face_cells").number);
        grid_.cellcount = static_cast<std::size_t>(grid.at("cell_count").number);
        columns_ = static_cast<int>(grid.at("export_columns").number);
        rows_ = static_cast<int>(grid.at("export_rows").number);
    }
    else
        throw std::runtime_error(manifestpath.string() + ": unsupported grid kind " + grid_.kind);
    recipe_ = manifest.at("recipe").text;
    for (const JsonValue& field : manifest.at("fields").items)
        fields_.push_back(manifestentry(field));
    for (const JsonValue& field : manifest.at("exports").items)
        exports_.push_back(manifestentry(field));
    for (const JsonValue& cells : manifest.at("cells").items)
        cells_.push_back(manifestentry(cells));
    for (const JsonValue& table : manifest.at("tables").items)
        tables_.push_back(manifestentry(table));
}

std::vector<unsigned char> Reader::fielddata(const ManifestEntry& field) const
{
    const std::filesystem::path path = directory_ / field.file;
    std::ifstream in(path, std::ios::binary);
    unsigned char prefix[10] = {};
    in.read(reinterpret_cast<char*>(prefix), sizeof(prefix));
    if (!in || std::string(reinterpret_cast<const char*>(prefix), 6) != "\x93NUMPY" || prefix[6] != 1)
        throw std::runtime_error("Not a version 1 .npy file: " + path.string());
    std::string header(static_cast<std::size_t>(prefix[8]) | (static_cast<std::size_t>(prefix[9]) << 8), ' ');
    in.read(header.data(), static_cast<std::streamsize>(header.size()));
    // The header must be the one writenpy produces for the manifest's dtype and shape.
    const std::string expected = "{'descr': '" + field.dtype + "', 'fortran_order': False, 'shape': " + shapetext(field.shape, false) + ", }";
    if (!in || header.compare(0, expected.size(), expected) != 0)
        throw std::runtime_error("Field file does not match its manifest entry: " + path.string());

    // Every dtype of the catalogue ends in its element size: "<i2", "<f4", "|u1".
    std::size_t bytes = field.dtype.empty() ? 0 : static_cast<std::size_t>(field.dtype.back() - '0');
    for (const std::size_t extent : field.shape)
        bytes *= extent;
    std::vector<unsigned char> data(bytes);
    if (bytes)
        in.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(bytes));
    if (!in || in.peek() != std::ifstream::traits_type::eof())
        throw std::runtime_error("Field file does not hold its declared shape: " + path.string());
    return data;
}

std::string Reader::tabletext(const ManifestEntry& table) const
{
    return readtext(directory_ / table.file);
}
}
