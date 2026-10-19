// Bundle v2 (docs/architecture/sphere-mesh.md sections 6 and 7, T8 W1c): a cube
// world round-trips bit-identically with its cell table, the lat-lon export keeps
// a constant field exact, and a v1 bundle still loads as a lat-lon world.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "io/regional_bundle.hpp"
#include "io/run_bundle.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "world/planet.hpp"
#include "world/planet_bundle.hpp"
#include "world/planet_export.hpp"
#include "world/planet_fields.hpp"

namespace
{
using worldmodel::FieldId;

std::filesystem::path freshdirectory(const std::string& name)
{
    const std::filesystem::path directory = std::filesystem::temp_directory_path() / "uw-bundle-v2-tests" / name;
    std::filesystem::remove_all(directory);
    return directory;
}

runbundle::Metadata testmetadata()
{
    runbundle::Metadata metadata;
    metadata.seed = 1;
    metadata.recipe = "bundle-v2-test";
    return metadata;
}

// Seed 1: values of several element types and layouts on every cell.
void seedfields(planet& world)
{
    std::mt19937 random(1);
    std::uniform_real_distribution<float> height(-4000.0f, 4000.0f);
    std::uniform_int_distribution<int> code(0, 8);
    world.setsealevel(0);
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
    {
        world.setnomfloat(c, height(random));
        world.setmountainheight(c, height(random) * 0.25f);
        world.setriverdir(c, code(random));
        world.setriverjan(c, code(random) * 1000);
        world.settectonicconvergence(c, code(random) * 10);
        world.setisland(c, code(random) > 4);
        world.setjantempfloat(c, height(random) * 0.01f);
        world.setseasonalrainfloat(seasonjuly, c, height(random) * 0.05f + 200.0f);
        world.setownerpolityid(c, code(random) - 1);
    }
}

std::unique_ptr<planet> cubeplanet(int facecells)
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    return world;
}

const runbundle::ManifestEntry* findentry(const std::vector<runbundle::ManifestEntry>& entries, const std::string& id)
{
    for (const runbundle::ManifestEntry& entry : entries)
        if (entry.id == id)
            return &entry;
    return nullptr;
}

template<typename T>
std::vector<T> values(const runbundle::Reader& reader, const runbundle::ManifestEntry& entry)
{
    const std::vector<unsigned char> bytes = reader.fielddata(entry);
    std::vector<T> out(bytes.size() / sizeof(T));
    std::memcpy(out.data(), bytes.data(), out.size() * sizeof(T));
    return out;
}

std::string readtext(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
}

TEST(BundleV2, CubePlanetRoundTripsBitIdentically)
{
    const std::unique_ptr<planet> world = cubeplanet(8);
    world->setradius(3.0e6f);
    seedfields(*world);
    const std::filesystem::path directory = freshdirectory("cube-round-trip");
    writeplanetbundle(*world, directory, testmetadata());

    const runbundle::Reader reader(directory);
    EXPECT_EQ(reader.schema(), runbundle::kSchema);
    EXPECT_EQ(reader.grid().kind, "cube");
    EXPECT_EQ(reader.grid().facecells, 8);
    EXPECT_EQ(reader.columns(), 32);
    EXPECT_EQ(reader.rows(), 16);
    EXPECT_DOUBLE_EQ(reader.grid().radiuskm, 3000.0);
    const runbundle::ManifestEntry* base = findentry(reader.fields(), "geosphere.base_elevation");
    ASSERT_NE(base, nullptr);
    EXPECT_EQ(base->shape, (std::vector<std::size_t>{ 6, 8, 8 }));
    const runbundle::ManifestEntry* rain = findentry(reader.fields(), "atmosphere.precipitation");
    ASSERT_NE(rain, nullptr);
    EXPECT_EQ(rain->shape, (std::vector<std::size_t>{ 4, 6, 8, 8 }));

    auto loaded = std::make_unique<planet>();
    const PlanetBundleSource source = loadplanetbundle(*loaded, directory);
    ASSERT_EQ(loaded->grid().kind(), mesh::GridKind::cube);
    ASSERT_EQ(loaded->grid().facecells(), 8);
    EXPECT_EQ(loaded->grid().radius_m(), world->grid().radius_m());
    EXPECT_EQ(loaded->radius(), world->radius());

    for (const FieldId seeded : { FieldId::base_elevation, FieldId::mountain_height, FieldId::river_direction, FieldId::river_flow_january,
             FieldId::tectonic_convergence, FieldId::island, FieldId::temperature_january, FieldId::precipitation, FieldId::owner_polity })
        EXPECT_NE(std::find(source.loadedfields.begin(), source.loadedfields.end(), seeded), source.loadedfields.end())
            << worldmodel::qualifiedname(seeded);
    for (const FieldId id : source.loadedfields)
    {
        const PlanetFieldData before = readplanetfield(*world, id);
        const PlanetFieldData after = readplanetfield(*loaded, id);
        EXPECT_EQ(before.shape, after.shape) << worldmodel::qualifiedname(id);
        EXPECT_TRUE(before.bytes == after.bytes) << worldmodel::qualifiedname(id);
    }
    for (std::size_t index = 0; index < worldmodel::kTableCount; index++)
    {
        const auto id = static_cast<worldmodel::TableId>(index);
        EXPECT_EQ(planettabletext(*world, id), planettabletext(*loaded, id)) << worldmodel::qualifiedname(id);
    }
    std::filesystem::remove_all(directory);
}

TEST(BundleV2, CellTableMatchesCellGrid)
{
    const std::unique_ptr<planet> world = cubeplanet(8);
    const std::filesystem::path directory = freshdirectory("cell-table");
    writeplanetbundle(*world, directory, testmetadata());
    const runbundle::Reader reader(directory);
    const mesh::CellGrid& grid = world->grid();
    const std::size_t count = grid.cellcount();
    ASSERT_EQ(reader.cells().size(), 6u);

    const auto entry = [&](const char* id)
    {
        const runbundle::ManifestEntry* found = findentry(reader.cells(), id);
        EXPECT_NE(found, nullptr) << id;
        return found;
    };
    const std::vector<double> centre = values<double>(reader, *entry("centre"));
    const std::vector<double> latitude = values<double>(reader, *entry("latitude_deg"));
    const std::vector<double> longitude = values<double>(reader, *entry("longitude_deg"));
    const std::vector<double> area = values<double>(reader, *entry("area_m2"));
    const std::vector<std::int32_t> neighbours = values<std::int32_t>(reader, *entry("neighbours"));
    const std::vector<std::int8_t> reverse = values<std::int8_t>(reader, *entry("reverse"));
    EXPECT_EQ(entry("centre")->shape, (std::vector<std::size_t>{ count, 3 }));
    EXPECT_EQ(entry("neighbours")->shape, (std::vector<std::size_t>{ count, 8 }));
    ASSERT_EQ(centre.size(), 3 * count);
    ASSERT_EQ(neighbours.size(), 8 * count);
    ASSERT_EQ(reverse.size(), 8 * count);
    constexpr double degrees = 180.0 / 3.14159265358979323846;
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        EXPECT_EQ(centre[3 * c], p.x);
        EXPECT_EQ(centre[3 * c + 1], p.y);
        EXPECT_EQ(centre[3 * c + 2], p.z);
        EXPECT_EQ(latitude[c], grid.latitude_rad(c) * degrees);
        EXPECT_EQ(longitude[c], grid.longitude_rad(c) * degrees);
        EXPECT_EQ(area[c], grid.area_m2(c));
        for (int k = 1; k <= 8; k++)
        {
            EXPECT_EQ(neighbours[8 * c + (k - 1)], grid.neighbour(c, k));
            EXPECT_EQ(reverse[8 * c + (k - 1)], grid.reverse(c, k));
        }
    }
    std::filesystem::remove_all(directory);
}

TEST(BundleV2, ExportOfConstantFieldIsExact)
{
    const std::unique_ptr<planet> world = cubeplanet(8);
    world->setsealevel(0);
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        world->setnomfloat(c, 1234.5f);
        world->setjantempfloat(c, -3.25f);
        world->settectonicconvergence(c, 42);
        world->setriverdir(c, 3);
    }

    bool called = false;
    const std::unique_ptr<planet> latlon = exportlatlon(*world, 32, [&](planet&) { called = true; });
    EXPECT_TRUE(called);
    ASSERT_EQ(latlon->grid().kind(), mesh::GridKind::latlon);
    ASSERT_EQ(latlon->grid().columns(), 32);
    for (mesh::CellId c = 0; c < latlon->cellcount(); c++)
    {
        ASSERT_EQ(latlon->nomfloat(c), 1234.5f) << c;
        ASSERT_EQ(latlon->jantempfloat(c), -3.25f) << c;
        ASSERT_EQ(latlon->tectonicconvergence(c), 42) << c;
        ASSERT_EQ(latlon->riverdir(c), 0) << c;   // not recomputed by this callback
    }

    // The import back onto the cube keeps the constant too.
    const std::unique_ptr<planet> back = cubeplanet(8);
    importfields(*back, *latlon, { FieldId::base_elevation, FieldId::temperature_january });
    for (mesh::CellId c = 0; c < back->cellcount(); c++)
    {
        ASSERT_EQ(back->nomfloat(c), 1234.5f) << c;
        ASSERT_EQ(back->jantempfloat(c), -3.25f) << c;
    }
    EXPECT_THROW(importfields(*back, *latlon, { FieldId::river_direction }), std::logic_error);

    // The bundle's export: the same values, direction fields 0 and not recomputed.
    const std::filesystem::path directory = freshdirectory("constant-export");
    writeplanetbundle(*world, directory, testmetadata());
    const runbundle::Reader reader(directory);
    const runbundle::ManifestEntry* base = findentry(reader.exports(), "geosphere.base_elevation");
    ASSERT_NE(base, nullptr);
    EXPECT_EQ(base->rule, "mean");
    EXPECT_EQ(base->file, "latlon/geosphere/base_elevation.npy");
    EXPECT_EQ(base->shape, (std::vector<std::size_t>{ 16, 32 }));
    for (const float value : values<float>(reader, *base))
        ASSERT_EQ(value, 1234.5f);
    const runbundle::ManifestEntry* river = findentry(reader.exports(), "hydrosphere.river_direction");
    ASSERT_NE(river, nullptr);
    EXPECT_EQ(river->rule, "direction");
    EXPECT_FALSE(river->recomputed);
    for (const std::int16_t value : values<std::int16_t>(reader, *river))
        ASSERT_EQ(value, 0);
    const runbundle::ManifestEntry* rain = findentry(reader.exports(), "atmosphere.precipitation");
    ASSERT_NE(rain, nullptr);
    EXPECT_EQ(rain->shape, (std::vector<std::size_t>{ 4, 16, 32 }));
    EXPECT_EQ(findentry(reader.exports(), "geosphere.sea_level"), nullptr);   // scalars have no export
    std::filesystem::remove_all(directory);
}

TEST(BundleV2, V1BundleStillLoads)
{
    auto world = std::make_unique<planet>();
    world->setwidth(31);
    seedfields(*world);
    const std::filesystem::path directory = freshdirectory("v1");
    writeplanetbundle(*world, directory, testmetadata());
    {
        const runbundle::Reader reader(directory);
        EXPECT_EQ(reader.grid().kind, "latlon");
        EXPECT_TRUE(reader.cells().empty());
        EXPECT_TRUE(reader.exports().empty());
    }

    // The v1 manifest: its schema and its grid block without kind and radius.
    std::string manifest = readtext(directory / "manifest.json");
    const std::string v2 = std::string("\"") + runbundle::kSchema + "\"";
    ASSERT_NE(manifest.find(v2), std::string::npos);
    manifest.replace(manifest.find(v2), v2.size(), std::string("\"") + runbundle::kSchemaV1 + "\"");
    const std::size_t begin = manifest.find("\"grid\": {");
    const std::size_t end = manifest.find("},\n", begin);
    ASSERT_NE(begin, std::string::npos);
    ASSERT_NE(end, std::string::npos);
    manifest.replace(begin, end + 1 - begin,
        "\"grid\": {\"columns\": 32, \"rows\": 16, \"registration\": \"cell-centred\", \"layout\": \"row-major, north row first\"}");
    {
        std::ofstream out(directory / "manifest.json", std::ios::binary | std::ios::trunc);
        out << manifest;
    }

    auto loaded = std::make_unique<planet>();
    const PlanetBundleSource source = loadplanetbundle(*loaded, directory);
    EXPECT_EQ(source.reader->schema(), runbundle::kSchemaV1);
    EXPECT_EQ(source.reader->grid().kind, "latlon");
    ASSERT_EQ(loaded->grid().kind(), mesh::GridKind::latlon);
    ASSERT_EQ(loaded->grid().columns(), 32);
    ASSERT_FALSE(source.loadedfields.empty());
    for (const FieldId id : source.loadedfields)
    {
        const PlanetFieldData before = readplanetfield(*world, id);
        const PlanetFieldData after = readplanetfield(*loaded, id);
        EXPECT_EQ(before.shape, after.shape) << worldmodel::qualifiedname(id);
        EXPECT_TRUE(before.bytes == after.bytes) << worldmodel::qualifiedname(id);
    }
    std::filesystem::remove_all(directory);
}

// A cube-window bundle (sphere-mesh.md section 7, T8 W2-W1c): written from synthetic fields on
// the level-1 cells of a box across the face 0|1 edge, read back bit for bit with its keys, its
// cell table and its lat-lon export of the box. Kept for the Python round trip (uwtools.bundle)
// when UW_CUBE_WINDOW_BUNDLE names a fresh directory.
TEST(BundleV2, CubeWindowRoundTrips)
{
    const int n = 8, scale = 4, level = 1, columns = 4 * n;
    const double radius = 6371000.0;
    const landscape::RegionalWindow window = landscape::regionalwindow(columns, scale, 2, 22.5, -11.25, 67.5, 33.75);
    ASSERT_EQ(window.column0, 18);
    ASSERT_EQ(window.row0, 5);
    const auto base = mesh::makecubegrid(n, radius);
    landscape::KeyLevelWindow keyed;
    keyed.kind = mesh::GridKind::cube;
    keyed.basecells = n;
    keyed.factor = scale;
    keyed.level = level;
    keyed.keptglobal = landscape::cellsinpixelbox(*base, columns, window.column0, window.row0,
                                                  window.column0 + window.columns - 1, window.row0 + window.rows - 1);
    const auto kept = mesh::cubewindowgrid(n * scale, radius, keyed.keptkeys(), true);
    const std::size_t count = kept->cellcount();
    ASSERT_GT(count, 0u);
    std::set<int> faces;

    landscape::RegionalFields fields;
    fields.window = window;
    const auto fill = [count](auto& values, auto value) { values.assign(count, value); };
    fill(fields.elevation_m, 0.0f);
    fill(fields.sea, std::uint8_t{ 0 });
    fill(fields.riverdirection, std::uint8_t{ 0 });
    fill(fields.discharge_m3s, 1.5f);
    fill(fields.januaryflow, 0);
    fill(fields.julyflow, 0);
    fill(fields.lakesurface_m, 0.0f);
    fill(fields.januarytemperature_c, 0.0f);
    fill(fields.julytemperature_c, 0.0f);
    fill(fields.januaryprecipitation_mm, 7.0f);
    fill(fields.julyprecipitation_mm, 7.0f);
    fill(fields.icethickness_m, 0.0f);
    fill(fields.climate, std::int16_t{ 0 });
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::Vec3 p = kept->centre(c);
        const mesh::LevelKey key = kept->globalkey(c);
        faces.insert(mesh::cubecell(n * scale, key).face);
        fields.elevation_m[c] = static_cast<float>(1000.0 * p.z + 200.0 * p.x);
        fields.sea[c] = fields.elevation_m[c] < 0.0f ? 1 : 0;
        fields.riverdirection[c] = static_cast<std::uint8_t>(key % 9);
        fields.januaryflow[c] = static_cast<std::int32_t>(key % 1000);
        fields.julyflow[c] = static_cast<std::int32_t>(key % 997);
        fields.januarytemperature_c[c] = static_cast<float>(kept->latitude_rad(c));
        fields.climate[c] = static_cast<std::int16_t>(key % 30);
    }
    EXPECT_EQ(faces, (std::set<int>{ 0, 1 }));

    regionalbundle::Metadata metadata;
    metadata.run.seed = 1;
    metadata.run.recipe = "region";
    metadata.sourcebundle = "synthetic";
    metadata.sourcemanifestsha256 = std::string(64, 'a');
    metadata.sourceseed = 1;
    metadata.level = level;
    metadata.factor = scale;
    const char* keep = std::getenv("UW_CUBE_WINDOW_BUNDLE");
    const std::filesystem::path directory = keep != nullptr && *keep ? std::filesystem::path(keep) : freshdirectory("cube-window");
    regionalbundle::writecubewindowbundle(directory, fields, *kept, n, level, metadata);

    const runbundle::Reader reader(directory);
    EXPECT_EQ(reader.grid().kind, "cube-window");
    EXPECT_EQ(reader.grid().facecells, n * scale);
    EXPECT_EQ(reader.grid().level, level);
    EXPECT_EQ(reader.grid().basefacecells, n);
    EXPECT_EQ(reader.grid().cellcount, count);
    EXPECT_DOUBLE_EQ(reader.grid().radiuskm, 6371.0);
    EXPECT_EQ(reader.columns(), window.finewidth());
    EXPECT_EQ(reader.rows(), window.fineheight());
    const std::vector<mesh::LevelKey> keys = regionalbundle::cubewindowkeys(reader);
    ASSERT_EQ(keys.size(), count);
    for (mesh::CellId c = 0; c < count; c++)
        EXPECT_EQ(keys[c], kept->globalkey(c));

    const auto bytesof = [](const auto& values) {
        std::vector<unsigned char> bytes(values.size() * sizeof(values[0]));
        std::memcpy(bytes.data(), values.data(), bytes.size());
        return bytes;
    };
    std::size_t checked = 0;
    for (const runbundle::ManifestEntry& entry : reader.fields())
    {
        if (entry.id == "geosphere.sea_level")
            continue;
        EXPECT_EQ(entry.shape, (std::vector<std::size_t>{ count })) << entry.id;
        if (entry.id == "geosphere.base_elevation")
        {
            EXPECT_TRUE(reader.fielddata(entry) == bytesof(fields.elevation_m));
            checked++;
        }
        if (entry.id == "hydrosphere.river_flow_january")
        {
            EXPECT_TRUE(reader.fielddata(entry) == bytesof(fields.januaryflow));
            checked++;
        }
    }
    EXPECT_EQ(checked, 2u);
    for (const runbundle::ManifestEntry& entry : reader.cells())
        if (entry.id == "neighbours")
        {
            const std::vector<unsigned char> bytes = reader.fielddata(entry);
            std::vector<std::int32_t> neighbours(bytes.size() / sizeof(std::int32_t));
            std::memcpy(neighbours.data(), bytes.data(), bytes.size());
            for (mesh::CellId c = 0; c < count; c++)
                for (int k = 1; k <= mesh::slotcount; k++)
                    EXPECT_EQ(neighbours[mesh::slotcount * c + (k - 1)], kept->neighbour(c, k));
        }

    // The export: every grid field on the box's level lat-lon grid; a constant mean exact where
    // the window covers the export cell, the elevation near its formula there.
    mesh::LatLonWindowBox box;
    box.column0 = static_cast<std::int64_t>(window.column0) * scale;
    box.columns = window.finewidth();
    box.row0 = window.row0 * scale;
    box.rows = window.fineheight();
    const auto exportgrid = mesh::latlonwindowgrid(columns * scale, radius, box);
    EXPECT_EQ(reader.exports().size(), reader.fields().size() - 1);
    std::size_t covered = 0;
    for (const runbundle::ManifestEntry& entry : reader.exports())
    {
        EXPECT_EQ(entry.shape, (std::vector<std::size_t>{ static_cast<std::size_t>(box.rows), static_cast<std::size_t>(box.columns) }))
            << entry.id;
        if (entry.id != "atmosphere.precipitation_july" && entry.id != "geosphere.base_elevation")
            continue;
        EXPECT_EQ(entry.rule, "mean");
        const std::vector<unsigned char> bytes = reader.fielddata(entry);
        std::vector<float> values(bytes.size() / sizeof(float));
        std::memcpy(values.data(), bytes.data(), bytes.size());
        for (mesh::CellId d = 0; d < exportgrid->cellcount(); d++)
        {
            const mesh::Vec3 p = exportgrid->centre(d);
            if (kept->locate(p) == mesh::nocell)
                continue;
            if (entry.id == "atmosphere.precipitation_july")
            {
                EXPECT_EQ(values[d], 7.0f);
                covered++;
            }
            else
                EXPECT_NEAR(values[d], 1000.0 * p.z + 200.0 * p.x, 60.0);
        }
    }
    EXPECT_GT(covered, exportgrid->cellcount() / 2);
    if (keep == nullptr || !*keep)
        std::filesystem::remove_all(directory);
}
