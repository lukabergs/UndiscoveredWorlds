#include "io/regional_bundle.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "core/world/world_fields.hpp"

namespace regionalbundle
{
namespace
{
std::string jsonnumber(double value)
{
    if (!std::isfinite(value))
        return "null";
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
}

// A field under its catalogue id: component, element type, layout and units from the
// catalogue; the element type must be the one the data carries.
template<typename T>
void addcataloguefield(runbundle::Writer& writer, std::vector<std::string>& written, const std::string& id,
                       const std::vector<T>& values, const char* dtype, const std::vector<std::size_t>& shape,
                       const std::string& description, const std::string& writtenby)
{
    worldmodel::FieldId field{};
    if (!worldmodel::findfield(id, field))
        throw std::logic_error("writeregionalbundle: " + id + " is not in the world catalogue");
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(field);
    if (std::string(worldmodel::elementname(spec.element)) != dtype || worldmodel::elementsize(spec.element) != sizeof(T))
        throw std::logic_error("writeregionalbundle: " + id + " is " + worldmodel::elementname(spec.element) +
                               " in the catalogue, not " + dtype + " of " + std::to_string(sizeof(T)) + " bytes");
    written.push_back(id);
    runbundle::FieldEntry entry;
    entry.id = id;
    entry.component = worldmodel::componentname(spec.component);
    entry.units = spec.units;
    entry.description = description;
    entry.layout = worldmodel::layoutname(spec.layout);
    entry.dtype = dtype;
    entry.shape = shape;
    entry.data = values.data();
    entry.bytes = values.size() * sizeof(T);
    entry.writtenby = writtenby;
    writer.addfield(entry);
}

// SHA-256 (FIPS 180-4).
class Sha256
{
public:
    void update(const unsigned char* data, std::size_t size)
    {
        for (std::size_t index = 0; index < size; index++)
        {
            block_[used_++] = data[index];
            if (used_ == 64)
            {
                compress();
                used_ = 0;
            }
        }
        length_ += static_cast<std::uint64_t>(size);
    }

    std::string hex()
    {
        const std::uint64_t bits = length_ * 8;
        const unsigned char one = 0x80, zero = 0;
        update(&one, 1);
        while (used_ != 56)
            update(&zero, 1);
        unsigned char tail[8];
        for (int index = 0; index < 8; index++)
            tail[index] = static_cast<unsigned char>(bits >> (56 - 8 * index));
        update(tail, 8);
        std::ostringstream out;
        for (const std::uint32_t word : state_)
            out << std::hex << std::setw(8) << std::setfill('0') << word;
        return out.str();
    }

private:
    static std::uint32_t rotate(std::uint32_t value, int bits) { return (value >> bits) | (value << (32 - bits)); }

    void compress()
    {
        static constexpr std::uint32_t k[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
        std::uint32_t w[64];
        for (int index = 0; index < 16; index++)
            w[index] = (static_cast<std::uint32_t>(block_[4 * index]) << 24) |
                       (static_cast<std::uint32_t>(block_[4 * index + 1]) << 16) |
                       (static_cast<std::uint32_t>(block_[4 * index + 2]) << 8) | static_cast<std::uint32_t>(block_[4 * index + 3]);
        for (int index = 16; index < 64; index++)
        {
            const std::uint32_t s0 = rotate(w[index - 15], 7) ^ rotate(w[index - 15], 18) ^ (w[index - 15] >> 3);
            const std::uint32_t s1 = rotate(w[index - 2], 17) ^ rotate(w[index - 2], 19) ^ (w[index - 2] >> 10);
            w[index] = w[index - 16] + s0 + w[index - 7] + s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int index = 0; index < 64; index++)
        {
            const std::uint32_t t1 = h + (rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25)) + ((e & f) ^ (~e & g)) + k[index] + w[index];
            const std::uint32_t t2 = (rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += h;
    }

    std::array<std::uint32_t, 8> state_ = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                           0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    unsigned char block_[64] = {};
    std::size_t used_ = 0;
    std::uint64_t length_ = 0;
};
}

std::string filesha256(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("Cannot read " + path.string());
    const std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    Sha256 hash;
    hash.update(bytes.data(), bytes.size());
    return hash.hex();
}

std::string windowjson(const landscape::RegionalWindow& window)
{
    const int width = window.worldwidth, height = window.worldheight;
    // The east edge lies one turn further when the box crosses the date line: written
    // within -180 .. 180, so west > east then.
    int endcolumn = window.column0 + window.columns;
    if (endcolumn > width)
        endcolumn -= width;
    const double west = -180.0 + 360.0 * window.column0 / width;
    const double east = -180.0 + 360.0 * endcolumn / width;
    const double north = 90.0 - 180.0 * window.row0 / height;
    const double south = 90.0 - 180.0 * (window.row0 + window.rows) / height;
    std::ostringstream out;
    out << "{\"scale\": " << window.scale << ", \"world_columns\": " << width << ", \"world_rows\": " << height
        << ", \"column0\": " << window.column0 << ", \"columns\": " << window.columns << ", \"row0\": " << window.row0
        << ", \"rows\": " << window.rows
        << ", \"fine_column0\": " << static_cast<std::int64_t>(window.column0) * window.scale
        << ", \"fine_row0\": " << static_cast<std::int64_t>(window.row0) * window.scale
        << ", \"west\": " << jsonnumber(west) << ", \"south\": " << jsonnumber(south) << ", \"east\": " << jsonnumber(east)
        << ", \"north\": " << jsonnumber(north) << ", \"halo\": " << window.halo << "}";
    return out.str();
}

void writeregionalbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                         const Metadata& metadata)
{
    const landscape::RegionalWindow& window = fields.window;
    const bool pixelgrid = metadata.pixelcolumns.has_value();
    if (pixelgrid != metadata.pixelrows.has_value() || pixelgrid != metadata.pixelcells.has_value() ||
        (pixelgrid && (*metadata.pixelcolumns < 1 || *metadata.pixelrows < 1)))
        throw std::invalid_argument("writeregionalbundle: the pixel grid needs its columns, rows and cells, all positive");
    const std::size_t columns = static_cast<std::size_t>(pixelgrid ? *metadata.pixelcolumns : window.finewidth());
    const std::size_t rows = static_cast<std::size_t>(pixelgrid ? *metadata.pixelrows : window.fineheight());
    const std::size_t size = columns * rows;
    const auto checksize = [&](std::size_t count, const char* name)
    {
        if (count != size)
            throw std::invalid_argument(std::string("writeregionalbundle: ") + name +
                                    (pixelgrid ? " is not the pixel grid's size" : " is not the kept window's size"));
    };
    checksize(fields.elevation_m.size(), "elevation");
    checksize(fields.sea.size(), "sea");
    checksize(fields.riverdirection.size(), "river direction");
    checksize(fields.discharge_m3s.size(), "discharge");
    checksize(fields.januaryflow.size(), "January flow");
    checksize(fields.julyflow.size(), "July flow");
    checksize(fields.lakesurface_m.size(), "lake surface");
    checksize(fields.januarytemperature_c.size(), "January temperature");
    checksize(fields.julytemperature_c.size(), "July temperature");
    checksize(fields.januaryprecipitation_mm.size(), "January precipitation");
    checksize(fields.julyprecipitation_mm.size(), "July precipitation");
    checksize(fields.icethickness_m.size(), "ice thickness");
    checksize(fields.climate.size(), "climate class");

    runbundle::Writer writer(directory);
    const std::vector<std::size_t> grid = {rows, columns};
    const std::string& by = metadata.writtenby;
    std::vector<std::string> written;

    addcataloguefield(writer, written, "geosphere.base_elevation", fields.elevation_m, "<f4", grid,
                      "Regional terrain: the upsampled global base elevation with noise, carved by the analytical relief on the pinned channels", by);
    // The sea mask: no catalogue grid holds it (geosphere.sea_cells is a count).
    runbundle::FieldEntry sea;
    sea.id = "geosphere.sea";
    sea.component = "geosphere";
    sea.units = "0/1";
    sea.description = "Sea cells of the regional pass: elevation at or below the sea level and the global cell or one of its eight neighbours sea";
    sea.layout = "grid";
    sea.dtype = "|u1";
    sea.shape = grid;
    sea.data = fields.sea.data();
    sea.bytes = fields.sea.size();
    sea.writtenby = by;
    writer.addfield(sea);
    written.push_back(sea.id);
    const std::vector<std::int32_t> sealevel = {static_cast<std::int32_t>(metadata.sealevel)};
    addcataloguefield(writer, written, "geosphere.sea_level", sealevel, "<i4", {}, "Sea level of the world bundle on the elevation datum", by);

    // Not hydrosphere.river_direction: a world bundle sets that on every land cell, the
    // regional pass only on the pinned channels.
    std::vector<std::int16_t> direction(size);
    for (std::size_t cell = 0; cell < size; cell++)
        direction[cell] = static_cast<std::int16_t>(fields.riverdirection[cell]);
    runbundle::FieldEntry channel;
    channel.id = "hydrosphere.channel_direction";
    channel.component = "hydrosphere";
    channel.units = "direction code";
    channel.description = "Flow direction of the pinned channel cells (1..8 as river_direction), 0 elsewhere";
    channel.layout = "grid";
    channel.dtype = "<i2";
    channel.shape = grid;
    channel.data = direction.data();
    channel.bytes = direction.size() * sizeof(std::int16_t);
    channel.writtenby = by;
    writer.addfield(channel);
    written.push_back(channel.id);
    addcataloguefield(writer, written, "hydrosphere.discharge", fields.discharge_m3s, "<f4", grid,
                      "Discharge: the global discharge on pinned channels, the fine accumulation elsewhere", by);
    addcataloguefield(writer, written, "hydrosphere.river_flow_january", fields.januaryflow, "<i4", grid,
                      "January river flow of the global network on pinned channels, 0 elsewhere", by);
    addcataloguefield(writer, written, "hydrosphere.river_flow_july", fields.julyflow, "<i4", grid,
                      "July river flow of the global network on pinned channels, 0 elsewhere", by);
    addcataloguefield(writer, written, "hydrosphere.lake_surface", fields.lakesurface_m, "<f4", grid,
                      "Lake surface of the global lake where the fine bed lies below it, 0 elsewhere", by);
    addcataloguefield(writer, written, "cryosphere.ice_thickness", fields.icethickness_m, "<f4", grid,
                      "Present-day ice thickness: the global ice upsampled where the global cell or a neighbour "
                      "has ice; with regional.ice_ela_thickness_fraction above 0, at least that fraction of the "
                      "height above the present climate's equilibrium line (1 m to ice_ela_thickness_max_m) on "
                      "fine land and lake cells above it (plus an edge noise); 0 elsewhere",
                      by);
    addcataloguefield(writer, written, "atmosphere.temperature_january", fields.januarytemperature_c, "<f4", grid,
                      "January temperature, interpolated and corrected by the lapse rate for the fine elevation", by);
    addcataloguefield(writer, written, "atmosphere.temperature_july", fields.julytemperature_c, "<f4", grid,
                      "July temperature, interpolated and corrected by the lapse rate for the fine elevation", by);
    addcataloguefield(writer, written, "atmosphere.precipitation_january", fields.januaryprecipitation_mm, "<f4", grid,
                      "January precipitation, interpolated and clamped at zero", by);
    addcataloguefield(writer, written, "atmosphere.precipitation_july", fields.julyprecipitation_mm, "<f4", grid,
                      "July precipitation, interpolated and clamped at zero", by);
    addcataloguefield(writer, written, "atmosphere.climate_class", fields.climate, "<i2", grid,
                      "Climate class of the global cell", by);

    runbundle::Metadata manifest = metadata.run;
    manifest.columns = static_cast<int>(columns);
    manifest.rows = static_cast<int>(rows);
    for (runbundle::StepRecord& step : manifest.steps)
        step.writes = written;
    // The refinement fields, each only when set, before the window block's closing brace.
    std::string windowblock = windowjson(window);
    std::ostringstream refinement;
    if (metadata.level)
        refinement << ", \"level\": " << *metadata.level;
    if (metadata.factor)
        refinement << ", \"factor\": " << *metadata.factor;
    if (metadata.metresperpixel)
        refinement << ", \"metres_per_pixel\": " << jsonnumber(*metadata.metresperpixel);
    if (metadata.resampled)
        refinement << ", \"resampled\": " << (*metadata.resampled ? "true" : "false");
    if (pixelgrid)
        refinement << ", \"pixel_columns\": " << *metadata.pixelcolumns << ", \"pixel_rows\": " << *metadata.pixelrows
                   << ", \"pixel_cells\": " << jsonnumber(*metadata.pixelcells);
    windowblock.insert(windowblock.size() - 1, refinement.str());
    manifest.blocks.push_back({"window", windowblock});
    std::ostringstream source;
    source << "{\"bundle\": " << runbundle::jsonstring(metadata.sourcebundle)
           << ", \"manifest_sha256\": " << runbundle::jsonstring(metadata.sourcemanifestsha256)
           << ", \"seed\": " << metadata.sourceseed << "}";
    manifest.blocks.push_back({"source", source.str()});
    std::ostringstream diagnostics;
    diagnostics << "{\"blockmeanerror_m\": " << jsonnumber(fields.blockmeanerror_m)
                << ", \"reliefconverged\": " << (fields.reliefconverged ? "true" : "false")
                << ", \"changedreceivers\": " << fields.changedreceivers
                << ", \"meanlowering_m\": " << jsonnumber(fields.meanlowering_m)
                << ", \"maxlowering_m\": " << jsonnumber(fields.maxlowering_m)
                << ", \"sharepasses\": " << fields.sharepasses << ", \"sharemisses\": " << fields.sharemisses;
    // A refinement's levels (multiscale plan S2.3), only when there are any, so a single
    // pass's block is the one written before them.
    if (!metadata.levels.empty())
    {
        diagnostics << ", \"levels\": [";
        for (std::size_t index = 0; index < metadata.levels.size(); index++)
        {
            const landscape::LevelDiagnostics& level = metadata.levels[index];
            const landscape::ProjectionReport& projection = level.projection;
            diagnostics << (index > 0 ? ", " : "") << "{\"level\": " << level.level
                        << ", \"cellsize_km\": " << jsonnumber(level.cellsize_km)
                        << ", \"kept_columns\": " << level.window.columns << ", \"kept_rows\": " << level.window.rows
                        << ", \"domain_columns\": " << level.window.domaincolumns
                        << ", \"domain_rows\": " << level.window.domainrows << ", \"halo\": " << level.window.halo
                        << ", \"blockmeanerror_m\": " << jsonnumber(level.blockmeanerror_m)
                        << ", \"parents\": " << projection.parents
                        << ", \"maxparentmeanerror_m\": " << jsonnumber(projection.maxparentmeanerror_m)
                        << ", \"globalmeanerror_m\": " << jsonnumber(level.globalmeanerror_m)
                        << ", \"smoothpasses\": " << level.smoothing.passes
                        << ", \"maxlambdabefore_m\": " << jsonnumber(level.smoothing.maxlambdabefore_m)
                        << ", \"maxlambdastepbefore_m\": " << jsonnumber(level.smoothing.maxlambdastepbefore_m)
                        << ", \"maxsmoothmove_m\": " << jsonnumber(level.smoothing.maxmove_m)
                        << ", \"maxlambda_m\": " << jsonnumber(projection.maxlambda_m)
                        << ", \"maxlambdastep_m\": " << jsonnumber(projection.maxlambdastep_m)
                        << ", \"maxborderstep_m\": " << jsonnumber(projection.maxborderstep_m)
                        << ", \"borderstep99_m\": " << jsonnumber(projection.borderstep99_m)
                        << ", \"meanbound_m\": " << jsonnumber(level.meanbound_m)
                        << ", \"bedsmovedparents\": " << projection.bedsmovedparents
                        << ", \"forcedparents\": " << projection.forcedparents
                        << ", \"maxforcedmove_m\": " << jsonnumber(projection.maxforcedmove_m)
                        << ", \"driedcells\": " << level.driedcells << ", \"shoresbelowlake\": " << level.shoresbelowlake
                        << ", \"landbelowsea\": " << level.landbelowsea
                        << ", \"maxlandbelowsea_m\": " << jsonnumber(level.maxlandbelowsea_m)
                        << ", \"maxlandbelowsea_column\": " << level.maxlandbelowseacolumn
                        << ", \"maxlandbelowsea_row\": " << level.maxlandbelowsearow
                        << ", \"rimraised\": " << level.rimraised
                        << ", \"maxrimraise_m\": " << jsonnumber(level.maxrimraise_m)
                        << ", \"maxborderstep_parent_column\": " << projection.maxborderstepcolumn
                        << ", \"maxborderstep_parent_row\": " << projection.maxbordersteprow
                        << ", \"maxforced_parent_column\": " << projection.maxforcedcolumn
                        << ", \"maxforced_parent_row\": " << projection.maxforcedrow
                        << ", \"seainlandparents\": " << level.seainlandparents
                        << ", \"seabeyondparent\": " << level.seabeyondparent
                        << ", \"drybelowlake\": " << level.drybelowlake
                        // Phase S3's diagnostics (multiscale plan S3.0; refinement_levels.hpp).
                        << ", \"prolongationclamped\": " << level.prolongationclamped
                        << ", \"maxprolongationclamp_m\": " << jsonnumber(level.maxprolongationclamp_m)
                        << ", \"smoothmove99_m\": " << jsonnumber(level.smoothmove99_m)
                        << ", \"maxchannelrise_m\": " << jsonnumber(level.maxchannelrise_m)
                        << ", \"landcells\": " << level.landcells << ", \"channelcells\": " << level.channelcells
                        << ", \"drainagedensity_per_km\": " << jsonnumber(level.drainagedensity_per_km)
                        << ", \"tributaries\": " << level.tributaries << ", \"tributarycells\": " << level.tributarycells
                        << ", \"lakeedges\": " << level.lakeedges
                        << ", \"lakeedgealigned\": " << jsonnumber(level.lakeedgealigned)
                        << ", \"landroughness_m\": " << jsonnumber(level.landroughness_m)
                        << ", \"valleycells\": " << level.valleycells << ", \"troughcells\": " << level.troughcells
                        << ", \"cirquecells\": " << level.cirquecells
                        << ", \"ringbasincells\": " << level.ringbasincells
                        << ", \"unreachedcells\": " << level.unreachedcells
                        << ", \"sinkjunctions\": " << level.sinkjunctions
                        << ", \"diagonallinks\": " << level.diagonallinks
                        << ", \"cornerlinks\": " << level.cornerlinks
                        << ", \"mergedpieces\": " << level.mergedpieces << ", \"drytails\": " << level.drytails
                        << ", \"solvelowering50_m\": " << jsonnumber(level.solvelowering50_m)
                        << ", \"solvelowering99_m\": " << jsonnumber(level.solvelowering99_m)
                        << ", \"reliefconverged\": " << (level.reliefconverged ? "true" : "false")
                        << ", \"changedreceivers\": " << level.changedreceivers
                        << ", \"seconds\": " << jsonnumber(level.seconds) << ", \"steps\": [";
            for (std::size_t step = 0; step < level.steps.size(); step++)
                diagnostics << (step > 0 ? ", " : "") << runbundle::jsonstring(landscape::levelstepname(level.steps[step]));
            diagnostics << "]}";
        }
        diagnostics << "]";
    }
    diagnostics << "}";
    manifest.blocks.push_back({"diagnostics", diagnostics.str()});
    writer.finish(manifest);
}
}
