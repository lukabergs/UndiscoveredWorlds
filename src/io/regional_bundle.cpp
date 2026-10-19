#include "io/regional_bundle.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_remap.hpp"
#include "core/world/world_fields.hpp"

namespace regionalbundle
{
namespace
{
// Sees a field right after it is written, with its remap rule (a cube window's lat-lon export).
using ExportHook = std::function<void(const runbundle::FieldEntry&, worldmodel::RemapRule)>;

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
                       const std::string& description, const std::string& writtenby, const ExportHook* hook = nullptr)
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
    if (hook)
        (*hook)(entry, spec.rule);
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

namespace
{
void checkfieldsizes(const landscape::RegionalFields& fields, std::size_t size, const char* what)
{
    const auto checksize = [&](std::size_t count, const char* name)
    {
        if (count != size)
            throw std::invalid_argument(std::string("writeregionalbundle: ") + name + what);
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
}

// The regional fields on a grid of shape `grid` (size cells), in the order of every regional
// bundle; the hook, when set, sees each grid field right after it is written, with its remap rule.
std::vector<std::string> addregionalfields(runbundle::Writer& writer, const landscape::RegionalFields& fields,
                                           const std::vector<std::size_t>& grid, std::size_t size, const Metadata& metadata,
                                           const ExportHook* hook)
{
    const std::string& by = metadata.writtenby;
    std::vector<std::string> written;

    addcataloguefield(writer, written, "geosphere.base_elevation", fields.elevation_m, "<f4", grid,
                      "Regional terrain: the upsampled global base elevation with noise, carved by the analytical relief on the pinned channels", by, hook);
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
    if (hook)
        (*hook)(sea, worldmodel::RemapRule::nearest);
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
    if (hook)
        (*hook)(channel, worldmodel::RemapRule::direction);
    written.push_back(channel.id);
    addcataloguefield(writer, written, "hydrosphere.discharge", fields.discharge_m3s, "<f4", grid,
                      "Discharge: the global discharge on pinned channels, the fine accumulation elsewhere", by, hook);
    addcataloguefield(writer, written, "hydrosphere.river_flow_january", fields.januaryflow, "<i4", grid,
                      "January river flow of the global network on pinned channels, 0 elsewhere", by, hook);
    addcataloguefield(writer, written, "hydrosphere.river_flow_july", fields.julyflow, "<i4", grid,
                      "July river flow of the global network on pinned channels, 0 elsewhere", by, hook);
    addcataloguefield(writer, written, "hydrosphere.lake_surface", fields.lakesurface_m, "<f4", grid,
                      "Lake surface of the global lake where the fine bed lies below it, 0 elsewhere", by, hook);
    addcataloguefield(writer, written, "cryosphere.ice_thickness", fields.icethickness_m, "<f4", grid,
                      "Present-day ice thickness: the global ice upsampled where the global cell or a neighbour "
                      "has ice; with regional.ice_ela_thickness_fraction above 0, at least that fraction of the "
                      "height above the present climate's equilibrium line (1 m to ice_ela_thickness_max_m) on "
                      "fine land and lake cells above it (plus an edge noise); 0 elsewhere",
                      by, hook);
    addcataloguefield(writer, written, "atmosphere.temperature_january", fields.januarytemperature_c, "<f4", grid,
                      "January temperature, interpolated and corrected by the lapse rate for the fine elevation", by, hook);
    addcataloguefield(writer, written, "atmosphere.temperature_july", fields.julytemperature_c, "<f4", grid,
                      "July temperature, interpolated and corrected by the lapse rate for the fine elevation", by, hook);
    addcataloguefield(writer, written, "atmosphere.precipitation_january", fields.januaryprecipitation_mm, "<f4", grid,
                      "January precipitation, interpolated and clamped at zero", by, hook);
    addcataloguefield(writer, written, "atmosphere.precipitation_july", fields.julyprecipitation_mm, "<f4", grid,
                      "July precipitation, interpolated and clamped at zero", by, hook);
    addcataloguefield(writer, written, "atmosphere.climate_class", fields.climate, "<i2", grid,
                      "Climate class of the global cell", by, hook);

    return written;
}

// The manifest of a regional bundle of columns x rows (the fine grid, the pixel grid or a
// cube window's export box): the steps' writes, the window, source and diagnostics blocks.
runbundle::Metadata regionalmanifest(const landscape::RegionalFields& fields, const Metadata& metadata, std::size_t columns,
                                     std::size_t rows, const std::vector<std::string>& written)
{
    const landscape::RegionalWindow& window = fields.window;
    const bool pixelgrid = metadata.pixelcolumns.has_value();
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
                        << ", \"maxblockerror_column\": " << level.maxblockerrorcolumn
                        << ", \"maxblockerror_row\": " << level.maxblockerrorrow
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
                        << ", \"foldcells\": " << level.foldcells
                        << ", \"foldfloorclamped\": " << level.foldfloorclamped
                        << ", \"foldamplitude_p50_m\": " << jsonnumber(level.foldamplitude_p50_m)
                        << ", \"foldamplitude_p99_m\": " << jsonnumber(level.foldamplitude_p99_m)
                        << ", \"ringbasincells\": " << level.ringbasincells
                        << ", \"unreachedcells\": " << level.unreachedcells
                        << ", \"sinkjunctions\": " << level.sinkjunctions
                        << ", \"diagonallinks\": " << level.diagonallinks
                        << ", \"cornerlinks\": " << level.cornerlinks
                        << ", \"mergedpieces\": " << level.mergedpieces << ", \"drytails\": " << level.drytails
                        << ", \"drychannelinlake\": " << level.drychannelinlake
                        << ", \"channelbetweenlakes\": " << level.channelbetweenlakes
                        << ", \"driedchannels\": " << level.driedchannels
                        << ", \"wetchannelsabove\": " << level.wetchannelsabove
                        << ", \"wetchannelexcess_m\": " << jsonnumber(level.wetchannelexcess_m)
                        << ", \"pockets\": " << level.pockets << ", \"pocketcells\": " << level.pocketcells
                        << ", \"pocketsfilled\": " << level.pocketsfilled << ", \"pocketcellsfilled\": " << level.pocketcellsfilled
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
    return manifest;
}
}

void writeregionalbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                         const Metadata& metadata)
{
    const bool pixelgrid = metadata.pixelcolumns.has_value();
    if (pixelgrid != metadata.pixelrows.has_value() || pixelgrid != metadata.pixelcells.has_value() ||
        (pixelgrid && (*metadata.pixelcolumns < 1 || *metadata.pixelrows < 1)))
        throw std::invalid_argument("writeregionalbundle: the pixel grid needs its columns, rows and cells, all positive");
    const std::size_t columns = static_cast<std::size_t>(pixelgrid ? *metadata.pixelcolumns : fields.window.finewidth());
    const std::size_t rows = static_cast<std::size_t>(pixelgrid ? *metadata.pixelrows : fields.window.fineheight());
    const std::size_t size = columns * rows;
    checkfieldsizes(fields, size, pixelgrid ? " is not the pixel grid's size" : " is not the kept window's size");
    runbundle::Writer writer(directory);
    const std::vector<std::string> written = addregionalfields(writer, fields, {rows, columns}, size, metadata, nullptr);
    writer.finish(regionalmanifest(fields, metadata, columns, rows, written));
}

namespace
{
template<typename T>
void addarray(runbundle::Writer& writer, const char* id, const char* dtype, std::vector<std::size_t> shape, const std::vector<T>& values)
{
    runbundle::ArrayEntry entry;
    entry.id = id;
    entry.dtype = dtype;
    entry.shape = std::move(shape);
    entry.data = values.data();
    entry.bytes = values.size() * sizeof(T);
    writer.addcells(entry);
}

// The cell table of a cube window (planet_bundle.cpp's for a cube world, plus the level keys).
void writewindowcells(runbundle::Writer& writer, const mesh::CellGrid& grid)
{
    const std::size_t count = grid.cellcount();
    constexpr double degrees = 180.0 / 3.14159265358979323846;
    std::vector<std::uint64_t> keys(count);
    std::vector<double> centre(3 * count), latitude(count), longitude(count), area(count);
    std::vector<std::int32_t> neighbours(mesh::slotcount * count);
    std::vector<std::int8_t> reverse(mesh::slotcount * count);
    for (mesh::CellId c = 0; c < count; c++)
    {
        keys[c] = grid.globalkey(c);
        const mesh::Vec3 p = grid.centre(c);
        centre[3 * c] = p.x;
        centre[3 * c + 1] = p.y;
        centre[3 * c + 2] = p.z;
        latitude[c] = grid.latitude_rad(c) * degrees;
        longitude[c] = grid.longitude_rad(c) * degrees;
        area[c] = grid.area_m2(c);
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            neighbours[mesh::slotcount * c + (k - 1)] = grid.neighbour(c, k);
            reverse[mesh::slotcount * c + (k - 1)] = static_cast<std::int8_t>(grid.reverse(c, k));
        }
    }
    const std::size_t slots = static_cast<std::size_t>(mesh::slotcount);
    addarray(writer, "cell_ids", "<u8", {count}, keys);
    addarray(writer, "centre", "<f8", {count, 3}, centre);
    addarray(writer, "latitude_deg", "<f8", {count}, latitude);
    addarray(writer, "longitude_deg", "<f8", {count}, longitude);
    addarray(writer, "area_m2", "<f8", {count}, area);
    addarray(writer, "neighbours", "<i4", {count, slots}, neighbours);
    addarray(writer, "reverse", "|i1", {count, slots}, reverse);
}

// A field's values remapped to the export box by its rule; a direction field as 0. A cell the
// kept cells do not cover (no overlap) is not extrapolated: NaN for a float field, 0 for an
// integer one (the catalogue has no missing value; export_covered marks the cells).
template<typename T>
std::vector<unsigned char> exportvalues(const runbundle::FieldEntry& entry, worldmodel::RemapRule rule, const mesh::RemapWeights& weights,
                                        const mesh::CellGrid& source, const mesh::CellGrid& destination,
                                        const std::vector<std::uint8_t>& covered)
{
    const T* data = static_cast<const T*>(entry.data);
    const std::vector<T> values(data, data + entry.bytes / sizeof(T));
    std::vector<T> out;
    switch (rule)
    {
    case worldmodel::RemapRule::mean: out = mesh::remapmean(weights, source, destination, values); break;
    case worldmodel::RemapRule::total: out = mesh::remaptotal(weights, values); break;
    case worldmodel::RemapRule::nearest: out = mesh::remapnearest(source, destination, values); break;
    case worldmodel::RemapRule::direction: out.assign(destination.cellcount(), T{}); break;
    }
    for (std::size_t d = 0; d < out.size(); d++)
        if (!covered[d])
        {
            if constexpr (std::is_floating_point_v<T>)
                out[d] = std::numeric_limits<T>::quiet_NaN();
            else
                out[d] = T{};
        }
    std::vector<unsigned char> bytes(out.size() * sizeof(T));
    if (!bytes.empty())
        std::memcpy(bytes.data(), out.data(), bytes.size());
    return bytes;
}
}

void writecubewindowbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                           const mesh::CellGrid& kept, int basefacecells, int level, const Metadata& metadata)
{
    const landscape::RegionalWindow& window = fields.window;
    if (kept.kind() != mesh::GridKind::cube || !kept.issubgrid() || basefacecells < 1 ||
        static_cast<std::int64_t>(kept.facecells()) != static_cast<std::int64_t>(basefacecells) * window.scale ||
        window.worldwidth != 4 * basefacecells || level < 0)
        throw std::invalid_argument("writecubewindowbundle: the kept cells must be a window of N_L = N x scale of a world of 4 N columns");
    if (metadata.pixelcolumns || metadata.pixelrows || metadata.pixelcells)
        throw std::invalid_argument("writecubewindowbundle: a cube window has no pixel grid");
    const std::size_t count = kept.cellcount();
    checkfieldsizes(fields, count, " is not the cube window's cell count");

    // The export box: the kept box of the level's lat-lon grid (kW columns).
    const int scale = window.scale;
    mesh::LatLonWindowBox box;
    box.column0 = static_cast<std::int64_t>(window.column0) * scale;
    box.columns = window.finewidth();
    box.row0 = window.row0 * scale;
    box.rows = window.fineheight();
    const std::shared_ptr<const mesh::CellGrid> exportgrid = mesh::latlonwindowgrid(window.worldwidth * scale, kept.radius_m(), box);
    const mesh::RemapWeights weights = mesh::makeremap(kept, *exportgrid);
    const std::size_t columns = static_cast<std::size_t>(box.columns), rows = static_cast<std::size_t>(box.rows);
    // An export cell is covered when it overlaps a kept cell (its remap row has weight).
    std::vector<std::uint8_t> covered(exportgrid->cellcount());
    for (std::size_t d = 0; d < covered.size(); d++)
        covered[d] = weights.rowstart[d + 1] > weights.rowstart[d] ? 1 : 0;

    runbundle::Writer writer(directory);
    const ExportHook hook = [&](const runbundle::FieldEntry& entry, worldmodel::RemapRule rule) {
        std::vector<unsigned char> bytes;
        if (entry.dtype == "<f4")
            bytes = exportvalues<float>(entry, rule, weights, kept, *exportgrid, covered);
        else if (entry.dtype == "|u1")
            bytes = exportvalues<std::uint8_t>(entry, rule, weights, kept, *exportgrid, covered);
        else if (entry.dtype == "<i2")
            bytes = exportvalues<std::int16_t>(entry, rule, weights, kept, *exportgrid, covered);
        else if (entry.dtype == "<i4")
            bytes = exportvalues<std::int32_t>(entry, rule, weights, kept, *exportgrid, covered);
        else
            throw std::logic_error("writecubewindowbundle: no export for " + entry.id + " of " + entry.dtype);
        runbundle::FieldEntry exported = entry;
        exported.shape = {rows, columns};
        exported.data = bytes.data();
        exported.bytes = bytes.size();
        writer.addexport(exported, worldmodel::remaprulename(rule), false);
    };
    const std::vector<std::string> written = addregionalfields(writer, fields, {count}, count, metadata, &hook);
    // The export's coverage mask (export only: no native field).
    runbundle::FieldEntry mask;
    mask.id = "geosphere.export_covered";
    mask.component = "geosphere";
    mask.units = "0/1";
    mask.description = "Export cells that overlap a kept cell of the window: elsewhere float exports are NaN and integer ones 0";
    mask.layout = "grid";
    mask.dtype = "|u1";
    mask.shape = {rows, columns};
    mask.data = covered.data();
    mask.bytes = covered.size();
    mask.writtenby = metadata.writtenby;
    writer.addexport(mask, "coverage", false);
    writewindowcells(writer, kept);
    runbundle::Metadata manifest = regionalmanifest(fields, metadata, columns, rows, written);
    manifest.grid.kind = "cube-window";
    manifest.grid.facecells = kept.facecells();
    manifest.grid.level = level;
    manifest.grid.basefacecells = basefacecells;
    manifest.grid.cellcount = count;
    manifest.grid.radiuskm = kept.radius_m() / 1000.0;
    writer.finish(manifest);
}

std::vector<mesh::LevelKey> cubewindowkeys(const runbundle::Reader& reader)
{
    if (reader.grid().kind != "cube-window")
        throw std::runtime_error("cubewindowkeys: " + reader.directory().string() + " is not a cube-window bundle");
    for (const runbundle::ManifestEntry& entry : reader.cells())
    {
        if (entry.id != "cell_ids")
            continue;
        if (entry.dtype != "<u8" || entry.shape.size() != 1 || entry.shape[0] != reader.grid().cellcount)
            throw std::runtime_error("cubewindowkeys: cell_ids does not hold cell_count u8 keys: " + reader.directory().string());
        const std::vector<unsigned char> bytes = reader.fielddata(entry);
        std::vector<mesh::LevelKey> keys(entry.shape[0]);
        if (!keys.empty())
            std::memcpy(keys.data(), bytes.data(), bytes.size());
        return keys;
    }
    throw std::runtime_error("cubewindowkeys: no cell_ids in " + reader.directory().string());
}
}
