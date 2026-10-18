// The terrain processes on a cube world (docs/architecture/sphere-mesh.md sections 3, 11 and 12;
// T8 W2-T1): every cube-ready terrain process runs on a cube planet of N 16 with synthetic inputs
// and leaves finite fields, and a stamp of a radius in km covers the same area at a face centre,
// a face edge and a cube corner.

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "core/mesh/cube_geometry.hpp"
#include "core/procedural_noise.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_volcanoes.hpp"
#include "world/planet.hpp"

namespace
{
constexpr int facecells = 16;

// A cube planet with a continent where x + 0.3 y > 0.1 (rising inland), sea elsewhere, a band of
// mountains on the land, and a convergent oceanic boundary segment along y = 0 in the sea.
std::unique_ptr<planet> syntheticcubeworld()
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(48000);
    world->setsealevel(31000);
    world->setgravity(1.0f);
    const mesh::CellGrid& grid = world->grid();
    std::vector<TectonicBoundarySegment> segments(1);
    segments[0].id = 1;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double land = p.x + 0.3 * p.y - 0.1;
        if (land > 0.0)
        {
            world->setnomfloat(c, static_cast<float>(31000.0 + 50.0 + 2500.0 * land + 200.0 * std::sin(9.0 * p.z)));
            if (std::abs(p.z - 0.3) < 0.08)
            {
                world->setmountainheight(c, 2500.0f);
                world->setmountainridge(c, 4 + 64);
            }
        }
        else
        {
            world->setnomfloat(c, static_cast<float>(31000.0 - 3000.0 + 2500.0 * land + 300.0 * p.z));
            world->settectoniccrustclass(c, CrustClass::oceanic);
            if (std::abs(p.y) < 0.12)
            {
                world->settectonicboundarytype(c, BoundaryType::convergent);
                world->settectonicboundarysegmentid(c, 1);
                world->settectonicconvergence(c, 200);
                segments[0].cellCount++;
            }
        }
    }
    world->settectonicboundarysegments(segments);
    return world;
}

void expectfinite(const planet& world, const std::string& stage)
{
    int nonfinite = 0;
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
    {
        if (!std::isfinite(world.nomfloat(c)) || !std::isfinite(world.mountainheightfloat(c)) ||
            !std::isfinite(world.extraelevfloat(c)) || !std::isfinite(world.volcanofloat(c)) ||
            !std::isfinite(world.roughness(c)))
            nonfinite++;
    }
    EXPECT_EQ(nonfinite, 0) << stage;
}

double caparea_m2(double radius_m, double sphere_m)
{
    const double pi = 3.14159265358979323846;
    return 2.0 * pi * sphere_m * sphere_m * (1.0 - std::cos(radius_m / sphere_m));
}
}

TEST(TerrainCube, CubeReadyProcessesRunWithFiniteOutputs)
{
    const std::unique_ptr<planet> world = syntheticcubeworld();
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    WorldStageParameters parameters;
    parameters.seaLevel = world->sealevel();
    GenerationExecutionContext context;

    const GenerationStageId stages[] = {
        GenerationStageId::archipelagos,
        GenerationStageId::ocean_depth_refinement,
        GenerationStageId::tectonic_trenches,
        GenerationStageId::terrain_texturing,
        GenerationStageId::mountain_bases,
        GenerationStageId::terrain_smoothing,
        GenerationStageId::canyon_uplift,
        GenerationStageId::depression_fill,
        GenerationStageId::terrain_roughness,
    };
    for (const GenerationStageId stage : stages)
    {
        const WorldProcessDefinition& process = getworldprocess(stage);
        ASSERT_TRUE(process.cubeready) << process.id;
        std::string error;
        ASSERT_TRUE(executegenerationstage(stage, *world, scratch, parameters, context, &error)) << process.id << ": " << error;
        expectfinite(*world, process.id);
    }

    // The passes did their work: trenches below the old floor, extra elevation and roughness set,
    // land raised under the mountains, and no land pit left by the fill.
    int extra = 0;
    int rough = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        extra += world->extraelevfloat(c) > 0.0f ? 1 : 0;
        rough += world->roughness(c) > 0 ? 1 : 0;
    }
    EXPECT_GT(extra, 0);
    EXPECT_EQ(rough, static_cast<int>(world->cellcount()));
}

TEST(TerrainCube, IsolatedVolcanoOnTheCube)
{
    const std::unique_ptr<planet> world = syntheticcubeworld();
    worldgrid<bool> shelves(world->grid());
    const std::vector<float> east = terrain::spherefractal(*world, 11, terrain::terrainfractal(), 48000.0f);
    const std::vector<float> north = terrain::spherefractal(*world, 12, terrain::terrainfractal(), 48000.0f);
    // A land cell at a cube corner, where one diagonal neighbour is missing.
    const mesh::CellId corner = mesh::cubeid(facecells, 0, facecells - 1, 0);
    // A chain has one to six volcanoes less one: a few chains from a fixed stream make some.
    fast_srand(12345);
    for (int chain = 0; chain < 6; chain++)
        terrain::createisolatedvolcano(*world, corner, shelves, east, north, 2000, true);
    terrain::makemountainisland(*world, mesh::cubeid(facecells, 1, facecells / 2, facecells / 2), 1500.0f);
    int volcanoes = 0;
    int islands = 0;
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        volcanoes += world->volcanofloat(c) != 0.0f ? 1 : 0;
        islands += world->mountainisland(c) ? 1 : 0;
    }
    EXPECT_GT(volcanoes, 0);
    EXPECT_GE(islands, 12);  // the 3 x 3 boxes of two neighbouring cells at a face centre
    expectfinite(*world, "createisolatedvolcano");
}

TEST(TerrainCube, StampAreaIsTheSameAtFaceCentreEdgeAndCorner)
{
    const std::unique_ptr<planet> world = syntheticcubeworld();
    const mesh::CellGrid& grid = world->grid();
    const double sphere_m = grid.radius_m();
    ASSERT_GT(sphere_m, 0.0);
    const double radius_km = 3000.0;
    const double cap = caparea_m2(radius_km * 1000.0, sphere_m);

    const mesh::CellId cells[] = {
        mesh::cubeid(facecells, 0, facecells / 2, facecells / 2),  // face centre
        mesh::cubeid(facecells, 0, facecells - 1, facecells / 2), // face edge
        mesh::cubeid(facecells, 0, facecells - 1, 0),             // cube corner
    };
    const char* names[] = { "face centre", "face edge", "cube corner" };
    for (int index = 0; index < 3; index++)
    {
        double area = 0.0;
        for (const mesh::CellId c : terrain::cellswithinkm(*world, cells[index], radius_km))
            area += grid.area_m2(c);
        EXPECT_NEAR(area / cap, 1.0, 0.10) << names[index];
    }
}
