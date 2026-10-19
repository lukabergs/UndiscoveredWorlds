// The last non-engine processes flagged cube-ready (docs/architecture/sphere-mesh.md section 12,
// T8 W2-F): each runs on a cube planet of N 16 with synthetic inputs through the stage executor,
// leaves finite fields and changes its declared writes where expected. Also the land and sea
// totals and the field storage (read and load by cell id) on a cube. Seed 1.

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "core/classes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/world/world_fields.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "world/planet.hpp"
#include "world/planet_fields.hpp"
#include "world/world_parameters.hpp"
#include "world/world_setup.hpp"

namespace
{
using S = GenerationStageId;
using worldmodel::FieldId;

constexpr int facecells = 16;
constexpr float sealevel = 31000.0f;

// A cube planet with a continent where x + 0.3 y > 0.1 (rising inland), sea elsewhere, a
// convergent oceanic boundary along y = 0 in the sea (boundary distance 0 there, 30 elsewhere).
std::unique_ptr<planet> continentworld()
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(48000);
    world->setsealevel(static_cast<int>(sealevel));
    world->setgravity(1.0f);
    const mesh::CellGrid& grid = world->grid();
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double land = p.x + 0.3 * p.y - 0.1;
        world->settectonicboundarydistance(c, 30);
        if (land > 0.0)
        {
            world->setnomfloat(c, static_cast<float>(sealevel + 50.0 + 2500.0 * land + 200.0 * std::sin(9.0 * p.z)));
            continue;
        }
        world->setnomfloat(c, static_cast<float>(sealevel - 3000.0 + 2500.0 * land + 300.0 * p.z));
        world->settectoniccrustclass(c, CrustClass::oceanic);
        if (std::abs(p.y) < 0.12)
        {
            world->settectonicboundarytype(c, BoundaryType::convergent);
            world->settectonicconvergence(c, 200);
            world->settectonicuplifttendency(c, 0.5f);
            world->settectonicboundarydistance(c, 0);
        }
    }
    return world;
}

bool runstage(GenerationStageId stage, planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context = {})
{
    const WorldProcessDefinition& process = getworldprocess(stage);
    EXPECT_TRUE(process.cubeready) << process.id;
    WorldStageParameters parameters;
    parameters.seaLevel = world.sealevel();
    std::string error;
    const bool ran = executegenerationstage(stage, world, scratch, parameters, context, &error);
    EXPECT_TRUE(ran) << process.id << ": " << error;
    return ran;
}

int nonfinite(const planet& world)
{
    int count = 0;
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        if (!std::isfinite(world.nomfloat(c)) || !std::isfinite(world.volcanofloat(c)) || !std::isfinite(world.roughness(c)))
            count++;
    return count;
}

double landshare(const planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    double land = 0.0;
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        land += world.basebelowsea(c) ? 0.0 : grid.area_m2(c);
    return land / grid.totalarea_m2();
}

// Twelve 5 x 5 discs, the template arrays the wetlands take.
struct Templates
{
    boolshapetemplate shapes[12];
    Templates()
    {
        for (boolshapetemplate& shape : shapes)
        {
            shape.setxsize(4);
            shape.setysize(4);
            for (int i = 0; i <= 4; i++)
                for (int j = 0; j <= 4; j++)
                    shape.setpoint(i, j, (i - 2) * (i - 2) + (j - 2) * (j - 2) <= 4);
        }
    }
};

// Restores a run-time tuning value at the end of a test.
template<typename T>
struct Restore
{
    T& value;
    T saved;
    explicit Restore(T& target) : value(target), saved(target) {}
    ~Restore() { value = saved; }
};
}

TEST(CubeFlags, InlandSeaPolicyFillsAndDrains)
{
    // Land at 500 m with an ocean where x < -0.2 and two closed seas 500 m deep inland.
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(48000);
    world->setsealevel(static_cast<int>(sealevel));
    const mesh::CellGrid& grid = world->grid();
    const mesh::Vec3 first = mesh::normalise(mesh::Vec3{ 0.3, 0.9, 0.3 });
    const mesh::Vec3 second = mesh::normalise(mesh::Vec3{ 0.3, -0.9, -0.3 });
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const bool inland = mesh::dot(p, first) > std::cos(0.2) || mesh::dot(p, second) > std::cos(0.2);
        world->setnomfloat(c, p.x < -0.2 ? sealevel - 3000.0f : inland ? sealevel - 500.0f : sealevel + 500.0f);
    }

    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    detectinlandseas(*world, scratch);
    ASSERT_EQ(scratch.inlandSeaComponents.size(), 2u);
    scratch.inlandSeaComponents[0].policy = GenerationComponentPolicy::fill;
    scratch.inlandSeaComponents[1].policy = GenerationComponentPolicy::drain;
    const std::vector<int> components = scratch.inlandSeaComponentIds;

    ASSERT_TRUE(runstage(S::inland_seas, *world, scratch));
    EXPECT_EQ(nonfinite(*world), 0);
    // The filled sea is land, the drained one joined the ocean by a channel: no closed sea is left.
    int filled = 0, drained = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        filled += components[c] == 1 && !world->basebelowsea(c) ? 1 : 0;
        drained += components[c] == 2 && world->basebelowsea(c) ? 1 : 0;
    }
    EXPECT_GT(filled, 0);
    EXPECT_GT(drained, 0);
    EXPECT_TRUE(scratch.inlandSeaComponents.empty());
    // The declared write: the shelves of the coasts, rebuilt.
    int shelves = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        shelves += world->continentalshelves()[c] ? 1 : 0;
    EXPECT_GT(shelves, 0);
}

TEST(CubeFlags, TectonicVolcanoesRaiseConesAndHotspots)
{
    const std::unique_ptr<planet> world = continentworld();
    std::vector<float> before(world->cellcount());
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        before[c] = world->nomfloat(c);
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    ASSERT_TRUE(runstage(S::tectonic_volcanoes, *world, scratch));
    EXPECT_EQ(nonfinite(*world), 0);
    int volcanoes = 0, strato = 0, raised = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        volcanoes += world->volcanofloat(c) != 0.0f ? 1 : 0;
        strato += world->volcanofloat(c) != 0.0f && world->strato(c) ? 1 : 0;
        raised += world->nomfloat(c) > before[c] ? 1 : 0;
    }
    EXPECT_GT(volcanoes, 0);
    EXPECT_GT(strato, 0);   // the arc cones on the convergent boundary
    EXPECT_GT(raised, 0);   // the uplift kernels
}

TEST(CubeFlags, SeaLevelPlacesTheLandFraction)
{
    Restore<float> fraction(tuning::world::ocean::landFraction);
    Restore<bool> aftercoast(tuning::world::ocean::landFractionAfterCoast);
    tuning::world::ocean::landFraction = 0.3f;
    tuning::world::ocean::landFractionAfterCoast = true;

    const std::unique_ptr<planet> world = continentworld();
    const int oldlevel = world->sealevel();
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    ASSERT_TRUE(runstage(S::sea_level, *world, scratch));
    EXPECT_NE(world->sealevel(), oldlevel);
    // One cell of N 16 is 1/1536 of the surface.
    EXPECT_NEAR(landshare(*world), 0.3, 0.01);
    int shelves = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        shelves += world->continentalshelves()[c] ? 1 : 0;
    EXPECT_GT(shelves, 0);
}

TEST(CubeFlags, DeltasAndWetlandsAfterTheCoarseClimate)
{
    // The ClimateAdapter test's two continents, the coarse climate (32 columns, 2 coupling
    // iterations), roughness and the hydrology, then the deltas stage with distributaries set
    // on some land cells.
    auto world = std::make_unique<planet>();
    initialiseworld(*world);
    world->setseed(1);
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    const mesh::CellGrid& grid = world->grid();
    const mesh::Vec3 first = mesh::unitvector(0.3, 0.5);
    const mesh::Vec3 second = mesh::unitvector(-0.4, 2.8);
    const float level = static_cast<float>(world->sealevel());
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double land = std::max(mesh::dot(p, first) - 0.75, mesh::dot(p, second) - 0.8);
        world->setnomfloat(c, land > 0.0 ? level + 200.0f + static_cast<float>(land) * 8000.0f : level - 3000.0f);
    }
    getlandandseatotals(*world);
    world->climatesimulation().resolution = { 32, 32, 64, 16 };
    world->climatesimulation().couplingIterationLimit = 2;

    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    for (const GenerationStageId stage : { S::terrain_roughness, S::global_temperature, S::pressure, S::winds, S::ocean_currents,
             S::sea_surface_temperatures, S::sea_ice_and_tides, S::rainfall, S::flow_routing, S::rivers, S::lake_levels,
             S::mountain_temperature_lapse, S::climates_and_biomes })
        ASSERT_TRUE(runstage(stage, *world, scratch));

    int marked = 0;
    for (mesh::CellId c = 0; c < grid.cellcount() && marked < 24; c += 7)
    {
        if (world->sea(c) == 0 && world->lakesurfacefloat(c) == 0.0f)
        {
            world->setdeltadir(c, 5);
            marked++;
        }
    }
    ASSERT_GT(marked, 0);
    std::vector<float> roughness(grid.cellcount());
    int classesbefore = 0;
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        roughness[c] = world->roughness(c);
        classesbefore += world->sea(c) == 0 && world->climate(c) > 0 ? 1 : 0;
    }

    Templates templates;
    GenerationExecutionContext context;
    context.smalllake = templates.shapes;
    ASSERT_TRUE(runstage(S::deltas_wetlands_and_roughness, *world, scratch, context));
    EXPECT_EQ(nonfinite(*world), 0);

    int wetlands = 0, rougher = 0, land = 0, classes = 0;
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        wetlands += world->special(c) >= 130 && world->special(c) <= 132 ? 1 : 0;
        rougher += world->roughness(c) != roughness[c] ? 1 : 0;
        land += world->sea(c) == 0 ? 1 : 0;
        classes += world->sea(c) == 0 && world->climate(c) > 0 ? 1 : 0;
    }
    EXPECT_GT(wetlands, 0);
    EXPECT_EQ(rougher, static_cast<int>(grid.cellcount()));   // the regional roughness on every cell
    // The classes recomputed on the shadow and imported: as many land cells classified as after
    // climate.classification (some coastal land cells have none from either, 22 of 353 here).
    ASSERT_GT(land, 0);
    EXPECT_GT(classes, land / 2);
    EXPECT_EQ(classes, classesbefore);
}

TEST(CubeFlags, SubgridStatisticsOnTheCellGrid)
{
    const std::unique_ptr<planet> world = continentworld();
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        world->erodibility()[c] = 1.0f;
        world->runoff()[c] = 500.0f;
        world->erosionrate()[c] = 100.0f;
    }
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    ASSERT_TRUE(runstage(S::subgrid_statistics, *world, scratch));
    const planet& view = *world;
    int land = 0, relief = 0, spacing = 0, bad = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const float r = view.reliefamplitude()[c];
        const float s = view.valleyspacing()[c];
        bad += std::isfinite(r) && std::isfinite(s) ? 0 : 1;
        if (world->basebelowsea(c))
        {
            bad += r == 0.0f && s == 0.0f ? 0 : 1;
            continue;
        }
        land++;
        relief += r > 0.0f ? 1 : 0;
        spacing += s > 0.0f ? 1 : 0;
    }
    EXPECT_EQ(bad, 0);
    ASSERT_GT(land, 0);
    EXPECT_EQ(relief, land);
    EXPECT_EQ(spacing, land);
}

TEST(CubeFlags, LandAndSeaTotalsCountEveryCell)
{
    const std::unique_ptr<planet> world = continentworld();
    getlandandseatotals(*world);
    int land = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        land += world->basebelowsea(c) ? 0 : 1;
    EXPECT_EQ(world->landtotal(), land);
    EXPECT_EQ(world->landtotal() + world->seatotal(), static_cast<int>(world->cellcount()));
}

TEST(CubeFlags, GridFieldsReadAndLoadByCellId)
{
    // Every field stored by cell (shape [cell] or [layer, cell]) reads from a cube world and loads
    // into another one unchanged.
    const std::unique_ptr<planet> world = continentworld();
    const std::size_t cells = world->cellcount();
    int checked = 0;
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
        if (spec.layout == worldmodel::Layout::scalar)
            continue;
        const PlanetFieldData data = readplanetfield(*world, id);
        if (data.shape.empty() || data.shape.back() != cells)
            continue;
        auto copy = std::make_unique<planet>();
        copy->setgrid(mesh::GridKind::cube, 4 * facecells);
        ASSERT_TRUE(loadplanetfield(*copy, data)) << worldmodel::qualifiedname(id);
        const PlanetFieldData again = readplanetfield(*copy, id);
        EXPECT_EQ(again.shape, data.shape) << worldmodel::qualifiedname(id);
        EXPECT_EQ(again.bytes, data.bytes) << worldmodel::qualifiedname(id);
        checked++;
    }
    EXPECT_GT(checked, 50);
}

// The geology processes and layers.finalize on a cube world (T8 E5): the engine on the world's
// own cube (the identity hand-off) for 10 Myr, then the physical layers with the engine's
// deposits, each on the cell that contains its longitude and colatitude.
TEST(CubeFlags, GeologyAndFinalizeOnTheWorldCube)
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setgravity(1.0f);
    GenerationScratch scratch;
    WorldStageParameters parameters;
    parameters.tectonicWorldGrid = true;
    parameters.tectonicDurationMyr = 10;
    GenerationExecutionContext context;
    for (const S stage : { S::geosphere_initialization, S::plate_tectonics, S::finalize_layers })
    {
        const WorldProcessDefinition& process = getworldprocess(stage);
        EXPECT_TRUE(process.cubeready) << process.id;
        std::string error;
        ASSERT_TRUE(executegenerationstage(stage, *world, scratch, parameters, context, &error)) << process.id << ": " << error;
    }
    EXPECT_EQ(nonfinite(*world), 0);
    const double land = landshare(*world);
    EXPECT_GT(land, 0.05);
    EXPECT_LT(land, 0.95);
    int primary = 0;
    for (const mineralsystems::Deposit& deposit : world->mineraldeposits())
    {
        ASSERT_TRUE(deposit.surface.evaluated);
        ASSERT_LT(deposit.surface.cell, world->cellcount());
        if (deposit.system == mineralsystems::DepositSystem::placer_gold)
            continue;
        ++primary;
        const double pi = 3.14159265358979323846;
        EXPECT_EQ(deposit.surface.cell,
            world->grid().locate(mesh::unitvector(pi * (0.5 - deposit.v), pi * (2.0 * deposit.u - 1.0))));
    }
    EXPECT_GT(primary, 0);
}
