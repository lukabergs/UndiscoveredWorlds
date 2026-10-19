// The lat-lon export of a cube world for its consumers (docs/architecture/sphere-mesh.md
// sections 6, 7 and 9; T8 W2-IO): the routing recomputed on the export (directions
// downhill on the refilled surface, discharge accumulated and monotone downstream),
// settlements relocated to export cells, the bundle's recomputed and shadow exports,
// and the .uww round trip of the export.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "io/run_bundle.hpp"
#include "pipeline/climate_adapter.hpp"
#include "pipeline/export_routing.hpp"
#include "world/planet.hpp"
#include "world/planet_bundle.hpp"
#include "world/planet_export.hpp"
#include "world/planet_fields.hpp"
#include "world/world_setup.hpp"

namespace
{
using worldmodel::FieldId;

constexpr int facecells = 16;
constexpr int columns = 4 * facecells;
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr double pi = 3.14159265358979323846;

std::filesystem::path freshpath(const std::string& name)
{
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "uw-export-tests" / name;
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path.parent_path());
    return path;
}

// A cube N 16 world: a tilted plane of land over a sea, with a 600 m cone-shaped pit at
// (60 N, 0 E) that the routing has not filled, a wet and warm climate (no transmission
// losses: precipitation exceeds the open-water evaporation), wetter in January than in
// July, a drainage area on land (the world was routed) and river flows.
std::unique_ptr<planet> syntheticcubeworld()
{
    auto world = std::make_unique<planet>();
    initialiseworld(*world);
    world->setseed(1);
    world->setgrid(mesh::GridKind::cube, columns);
    world->setsealevel(0);
    world->setriverfactor(1.0f);
    const mesh::CellGrid& grid = world->grid();
    const double latitude = 60.0 * pi / 180.0;
    const mesh::Vec3 pit{ std::cos(latitude), 0.0, std::sin(latitude) };
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        double height = 2500.0 * p.z + 1200.0 * p.x - 400.0;
        const double distance = mesh::greatcircle_m(p, pit, grid.radius_m());
        if (distance < 500.0e3)
            height -= 600.0 * (1.0 - distance / 500.0e3);
        world->setnomfloat(c, static_cast<float>(height));
        world->filledelevation()[c] = static_cast<float>(height);
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            world->setseasonaltempfloat(season, c, 15.0f);
            world->setseasonalrainfloat(season, c, season == seasonjanuary ? 300.0f : season == seasonjuly ? 100.0f : 200.0f);
        }
        const bool land = height > 0.0;
        world->drainagearea()[c] = land ? 1.0f : 0.0f;
        world->setriverjan(c, land ? 300 : 0);
        world->setriverjul(c, land ? 100 : 0);
    }
    return world;
}

TEST(ExportCube, DirectionsPointDownhillOnTheRefilledSurface)
{
    const auto world = syntheticcubeworld();
    const auto exported = exportworld(*world);
    const planet& view = *exported;
    const mesh::CellGrid& grid = view.grid();
    ASSERT_EQ(grid.kind(), mesh::GridKind::latlon);
    ASSERT_EQ(grid.columns(), columns);

    int land = 0, raised = 0;
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        const float filled = view.filledelevation()[c];
        EXPECT_GE(filled, view.nomfloat(c)) << "cell " << c;
        if (filled > view.nomfloat(c) + 1.0f)
            raised++;
        if (view.basebelowsea(c))
        {
            EXPECT_EQ(view.flowdirection()[c], 0) << "sea cell " << c;
            continue;
        }
        land++;
        const int code = view.flowdirection()[c];
        ASSERT_GE(code, 1) << "land cell " << c;
        ASSERT_LE(code, 8) << "land cell " << c;
        EXPECT_EQ(view.riverdir(c), code) << "land cell " << c;
        const std::int32_t next = grid.neighbour(c, code);
        ASSERT_GE(next, 0) << "land cell " << c;
        const mesh::CellId r = static_cast<mesh::CellId>(next);
        if (!view.basebelowsea(r))
            EXPECT_LE(view.filledelevation()[r], filled) << "land cell " << c << " to " << r;
    }
    EXPECT_GT(land, 100);
    EXPECT_GT(raised, 0) << "the pit is filled on the export";
}

TEST(ExportCube, DischargeIsRecomputedAndGrowsDownstream)
{
    const auto world = syntheticcubeworld();
    const auto exported = exportworld(*world);
    const planet& view = *exported;
    const mesh::CellGrid& grid = view.grid();

    double runofftotal = 0.0, outlettotal = 0.0, outletarea = 0.0, landarea = 0.0;
    mesh::CellId largest = 0;
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        if (view.basebelowsea(c))
            continue;
        landarea += grid.area_m2(c);
        runofftotal += view.runoff()[c] / 1000.0 * grid.area_m2(c) / secondsperyear;
        const mesh::CellId r = static_cast<mesh::CellId>(grid.neighbour(c, view.flowdirection()[c]));
        if (view.basebelowsea(r))
        {
            outlettotal += view.discharge()[c];
            outletarea += view.drainagearea()[c] * 1.0e6;
            if (view.discharge()[c] > view.discharge()[largest])
                largest = c;
            continue;
        }
        EXPECT_GE(view.discharge()[r], view.discharge()[c]) << "land cell " << c << " to " << r;
        EXPECT_GE(view.drainagearea()[r], view.drainagearea()[c]) << "land cell " << c << " to " << r;
    }
    // Every land cell drains to one outlet: the outlets carry all the runoff and the land area.
    EXPECT_NEAR(outlettotal / runofftotal, 1.0, 1.0e-4);
    EXPECT_NEAR(outletarea / landarea, 1.0, 1.0e-4);

    // hydrology.rivers' legacy flows: the wetter January carries more.
    const int january = view.riverjan(largest), july = view.riverjul(largest);
    EXPECT_GT(july, 0);
    EXPECT_GT(january, july);
}

TEST(ExportCube, SettlementsAndDepositsMoveToTheExportCellOfTheirCentre)
{
    auto world = syntheticcubeworld();
    const mesh::CellGrid& cube = world->grid();
    for (const mesh::CellId c : { mesh::CellId{ 0 }, mesh::CellId{ 1000 }, static_cast<mesh::CellId>(cube.cellcount() - 1) })
    {
        Settlement settlement;
        settlement.id = static_cast<int>(world->settlements().size());
        settlement.cell = c;
        world->settlements().push_back(settlement);
        mineralsystems::Deposit deposit;
        deposit.id = static_cast<std::uint32_t>(world->mineraldeposits().size());
        deposit.surface.cell = c;
        world->mineraldeposits().push_back(deposit);
    }
    const auto exported = exportlatlon(*world, columns);
    const mesh::CellGrid& latlon = exported->grid();
    ASSERT_EQ(exported->settlements().size(), 3u);
    const double diagonal = std::sqrt(2.0) * 2.0 * pi * latlon.radius_m() / columns;
    for (std::size_t index = 0; index < 3; index++)
    {
        const Settlement& original = world->settlements()[index];
        const Settlement& moved = exported->settlements()[index];
        EXPECT_EQ(moved.cell, latlon.locate(cube.centre(original.cell)));
        EXPECT_EQ(moved.x, latlon.column(moved.cell));
        EXPECT_EQ(moved.y, latlon.row(moved.cell));
        EXPECT_LT(mesh::greatcircle_m(cube.centre(original.cell), latlon.centre(moved.cell), latlon.radius_m()), diagonal);
        const mineralsystems::SurfaceState& deposit = exported->mineraldeposits()[index].surface;
        EXPECT_EQ(deposit.cell, moved.cell);
        EXPECT_EQ(deposit.worldX, moved.x);
        EXPECT_EQ(deposit.worldY, moved.y);
    }
}

const runbundle::ManifestEntry* findentry(const std::vector<runbundle::ManifestEntry>& entries, FieldId id)
{
    const std::string name = worldmodel::qualifiedname(id);
    for (const runbundle::ManifestEntry& entry : entries)
        if (entry.id == name)
            return &entry;
    return nullptr;
}

TEST(ExportCube, BundleExportsTheRecomputedAndTheShadowFields)
{
    const auto world = syntheticcubeworld();
    // A stand-in for the climate shadow: a lat-lon world of the export's grid with climate's own fields.
    auto shadow = std::make_unique<planet>();
    shadow->setgrid(mesh::GridKind::latlon, columns);
    for (mesh::CellId c = 0; c < shadow->cellcount(); c++)
        shadow->setjanmountainraindir(c, static_cast<int>(c % 8) + 1);
    for (int x = 0; x < columns; x++)
        shadow->sethorse(x, 2, x % 7);

    PlanetBundleExport exports = bundleexport();
    exports.climateshadow = shadow.get();
    exports.shadowowned = climateshadowowned;
    runbundle::Metadata metadata;
    metadata.seed = 1;
    metadata.recipe = "export-test";
    const std::filesystem::path directory = freshpath("bundle");
    writeplanetbundle(*world, directory, metadata, {}, nullptr, {}, {}, exports);

    const runbundle::Reader reader(directory);
    const auto expected = exportworld(*world);
    for (const FieldId id : { FieldId::flow_direction, FieldId::river_direction, FieldId::discharge, FieldId::river_flow_january })
    {
        const runbundle::ManifestEntry* entry = findentry(reader.exports(), id);
        ASSERT_NE(entry, nullptr) << worldmodel::qualifiedname(id);
        EXPECT_TRUE(entry->recomputed) << worldmodel::qualifiedname(id);
        EXPECT_EQ(reader.fielddata(*entry), readplanetfield(*expected, id).bytes) << worldmodel::qualifiedname(id);
    }
    const runbundle::ManifestEntry* ridges = findentry(reader.exports(), FieldId::mountain_ridges);
    ASSERT_NE(ridges, nullptr);
    EXPECT_FALSE(ridges->recomputed);

    const runbundle::ManifestEntry* raindirection = findentry(reader.exports(), FieldId::mountain_precipitation_direction_january);
    ASSERT_NE(raindirection, nullptr);
    EXPECT_TRUE(raindirection->recomputed);
    EXPECT_EQ(reader.fielddata(*raindirection), readplanetfield(*shadow, FieldId::mountain_precipitation_direction_january).bytes);

    const runbundle::ManifestEntry* horse = findentry(reader.fields(), FieldId::subtropical_high_latitudes);
    ASSERT_NE(horse, nullptr);
    EXPECT_EQ(horse->shape, (std::vector<std::size_t>{ 6, static_cast<std::size_t>(columns) }));
    EXPECT_EQ(reader.fielddata(*horse), readplanetfield(*shadow, FieldId::subtropical_high_latitudes).bytes);
    EXPECT_EQ(findentry(reader.exports(), FieldId::subtropical_high_latitudes), nullptr);

    // The bundle still loads as a cube world; the lat-lon profile is carried, not loaded.
    auto loaded = std::make_unique<planet>();
    const PlanetBundleSource source = loadplanetbundle(*loaded, directory);
    EXPECT_EQ(loaded->grid().kind(), mesh::GridKind::cube);
    EXPECT_NE(std::find(source.carriedfields.begin(), source.carriedfields.end(), FieldId::subtropical_high_latitudes),
        source.carriedfields.end());

    // The desktop's world of the bundle: its export with the shadow's rain directions and profile.
    const auto desktop = exportbundleworld(*loaded, source);
    ASSERT_EQ(desktop->grid().kind(), mesh::GridKind::latlon);
    for (mesh::CellId c = 0; c < desktop->cellcount(); c++)
        ASSERT_EQ(desktop->janmountainraindir(c), shadow->janmountainraindir(c)) << "cell " << c;
    for (int x = 0; x < columns; x++)
        ASSERT_EQ(desktop->horse(x, 2), x % 7) << "column " << x;
}

TEST(ExportCube, UwwOfACubeWorldIsItsExport)
{
    auto world = syntheticcubeworld();
    Settlement settlement;
    settlement.id = 0;
    settlement.cell = 1000;
    world->settlements().push_back(settlement);
    const std::filesystem::path path = freshpath("cube.uww");
    saveworldfile(*world, path.string());
    const auto expected = exportworld(*world);

    // Planets are large: on the heap.
    auto holder = std::make_unique<planet>();
    planet& loaded = *holder;
    initialiseworld(loaded);
    ASSERT_TRUE(loaded.loadworld(path.string()));
    const mesh::CellGrid& grid = loaded.grid();
    ASSERT_EQ(grid.kind(), mesh::GridKind::latlon);
    ASSERT_EQ(grid.columns(), columns);
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        ASSERT_NEAR(loaded.nomfloat(c), expected->nomfloat(c), 0.5f) << "cell " << c;
        ASSERT_EQ(loaded.riverdir(c), expected->riverdir(c)) << "cell " << c;
        ASSERT_EQ(loaded.riverjan(c), expected->riverjan(c)) << "cell " << c;
        ASSERT_EQ(loaded.riverjul(c), expected->riverjul(c)) << "cell " << c;
    }
    ASSERT_EQ(loaded.settlements().size(), 1u);
    EXPECT_EQ(loaded.settlements()[0].x, expected->settlements()[0].x);
    EXPECT_EQ(loaded.settlements()[0].y, expected->settlements()[0].y);

    // Saved again, the loaded export is the same file.
    const std::filesystem::path again = freshpath("again.uww");
    saveworldfile(loaded, again.string());
    EXPECT_EQ(std::filesystem::file_size(again), std::filesystem::file_size(path));
}
}

