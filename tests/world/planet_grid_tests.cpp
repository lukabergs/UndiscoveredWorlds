// The world's cell grid (docs/architecture/sphere-mesh.md section 4, T8 W1b):
// fields by cell id on a cube planet, (x, y) and cell id agreeing on lat-lon,
// and a cube world refusing a process that is not cube-ready.

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "world/planet.hpp"
#include "world/planet_fields.hpp"

namespace
{
std::unique_ptr<planet> cubeplanet(int facecells)
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    return world;
}
}

TEST(PlanetGrid, CubePlanetSetsAndGetsByCell)
{
    const std::unique_ptr<planet> world = cubeplanet(8);
    ASSERT_EQ(world->grid().kind(), mesh::GridKind::cube);
    ASSERT_EQ(world->grid().facecells(), 8);
    ASSERT_EQ(world->cellcount(), 6u * 8u * 8u);
    EXPECT_EQ(world->width(), -1);
    EXPECT_EQ(world->height(), -1);

    world->setsealevel(0);
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        world->setnomfloat(c, 0.5f + static_cast<float>(c));
        world->setmountainheight(c, static_cast<float>(c % 7));
        world->setriverdir(c, static_cast<int>(c % 9));
        world->settectonicconvergence(c, static_cast<int>(c % 101));
        world->setjantempfloat(c, -10.0f + 0.25f * static_cast<float>(c));
        world->setseasonalrainfloat(seasonjuly, c, static_cast<float>(c) * 2.0f);
        world->setownerpolityid(c, static_cast<int>(c) - 3);
    }
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        EXPECT_EQ(world->nomfloat(c), 0.5f + static_cast<float>(c));
        EXPECT_EQ(world->mountainheightfloat(c), static_cast<float>(c % 7));
        EXPECT_EQ(world->mapfloat(c), world->nomfloat(c) + world->mountainheightfloat(c));
        EXPECT_FALSE(world->basebelowsea(c));
        EXPECT_EQ(world->riverdir(c), static_cast<int>(c % 9));
        EXPECT_EQ(world->tectonicconvergence(c), static_cast<int>(c % 101));
        EXPECT_EQ(world->jantempfloat(c), -10.0f + 0.25f * static_cast<float>(c));
        EXPECT_EQ(world->seasonalrainfloat(seasonjuly, c), static_cast<float>(c) * 2.0f);
        EXPECT_EQ(world->ownerpolityid(c), static_cast<int>(c) - 3);
    }

    // The fields are stored in id order: the bundle bytes are the cells in turn.
    const PlanetFieldData base = readplanetfield(*world, worldmodel::FieldId::base_elevation);
    ASSERT_EQ(base.shape, std::vector<std::size_t>{ world->cellcount() });
    const auto* values = reinterpret_cast<const float*>(base.bytes.data());
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
        EXPECT_EQ(values[c], world->nomfloat(c));

    // The (x, y) forms have no cell on the cube.
    EXPECT_EQ(world->nomfloat(0, 0), 0.0f);
}

TEST(PlanetGrid, LatLonXYEqualsByCell)
{
    auto world = std::make_unique<planet>();
    world->setwidth(31);
    const mesh::CellGrid& grid = world->grid();
    ASSERT_EQ(grid.kind(), mesh::GridKind::latlon);
    ASSERT_EQ(grid.columns(), 32);
    ASSERT_EQ(grid.rows(), 16);
    ASSERT_EQ(world->cellcount(), 32u * 16u);

    for (int x = 0; x <= world->width(); x++)
        for (int y = 0; y <= world->height(); y++)
        {
            const float value = static_cast<float>(1000 * x + y) + 0.25f;
            world->setnomfloat(x, y, value);
            world->setlakesurface(x, y, value * 0.5f);
            world->setspecial(x, y, (x + y) % 140);
            world->setjultempfloat(x, y, -value);
            world->setseasonaltempfloat(seasonapril, x, y, value + 1.0f);
            world->setfisheryreserve(x, y, (x * 7 + y) % 101);
        }
    for (int x = 0; x <= world->width(); x++)
        for (int y = 0; y <= world->height(); y++)
        {
            const mesh::CellId c = grid.cellof(x, y);
            ASSERT_EQ(c, static_cast<mesh::CellId>(y * 32 + x));
            EXPECT_EQ(world->nomfloat(c), world->nomfloat(x, y));
            EXPECT_EQ(world->mapfloat(c), world->mapfloat(x, y));
            EXPECT_EQ(world->sea(c), world->sea(x, y));
            EXPECT_EQ(world->truelake(c), world->truelake(x, y));
            EXPECT_EQ(world->lakesurfacefloat(c), world->lakesurfacefloat(x, y));
            EXPECT_EQ(world->special(c), world->special(x, y));
            EXPECT_EQ(world->jultempfloat(c), world->jultempfloat(x, y));
            EXPECT_EQ(world->jultemp(c), world->jultemp(x, y));
            EXPECT_EQ(world->seasonaltempfloat(seasonapril, c), world->seasonaltempfloat(seasonapril, x, y));
            EXPECT_EQ(world->avetempfloat(c), world->avetempfloat(x, y));
            EXPECT_EQ(world->fisheryreserve(c), world->fisheryreserve(x, y));
            EXPECT_EQ(world->climate().nomfloat(c), world->nomfloat(x, y));
        }

    // Row-major storage: the bundle field of shape [row, column] is the grid's memory.
    const PlanetFieldData base = readplanetfield(*world, worldmodel::FieldId::base_elevation);
    ASSERT_EQ(base.shape, (std::vector<std::size_t>{ 16, 32 }));
    const auto* values = reinterpret_cast<const float*>(base.bytes.data());
    EXPECT_EQ(values[3 * 32 + 5], world->nomfloat(5, 3));
}

TEST(PlanetGrid, CubeWorldRefusesANonReadyProcess)
{
    const std::unique_ptr<planet> world = cubeplanet(8);
    const WorldProcessDefinition& process = getworldprocess(GenerationStageId::terraforming);
    ASSERT_FALSE(process.cubeready);

    GenerationScratch scratch;
    WorldStageParameters parameters;
    GenerationExecutionContext context;
    std::string error;
    EXPECT_FALSE(executegenerationstage(GenerationStageId::terraforming, *world, scratch, parameters, context, &error));
    EXPECT_NE(error.find(process.id), std::string::npos) << error;
    EXPECT_NE(error.find("cube"), std::string::npos) << error;
}
