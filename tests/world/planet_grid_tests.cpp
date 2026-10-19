// The world's cell grid (docs/architecture/sphere-mesh.md section 4, T8 W1b):
// fields by cell id on a cube planet, (x, y) and cell id agreeing on lat-lon,
// and every process cube-ready (since T8 W2-FE).

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "core/climate_seasons.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "world/planet.hpp"
#include "world/planet_fields.hpp"
#include "world/world_setup.hpp"

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

// terrain.fluvial_erosion was the last process ported (T8 W2-FE): a cube world now refuses none,
// so the refusal in executegenerationstage has no process left to test it with.
TEST(PlanetGrid, EveryProcessIsCubeReady)
{
    for (const WorldProcessDefinition& process : getworldprocesses())
        EXPECT_TRUE(process.cubeready) << process.id;
}

TEST(PlanetGrid, MaxRiverFlowReadsEveryCubeCell)
{
    // setmaxriverflow loops over cells, not over the lat-lon width x height.
    const std::unique_ptr<planet> world = cubeplanet(8);
    const mesh::CellId last = static_cast<mesh::CellId>(world->cellcount() - 1);
    world->setriverjan(last, 300);
    world->setriverjul(last, 500);
    world->setmaxriverflow();
    EXPECT_EQ(world->maxriverflow(), 400);
    EXPECT_EQ(world->sharedgrid().get(), &world->grid());
}

// A .uww of version 18 (T8 E5): its tectonic records held raster centroids (cell centres at whole
// numbers, y south), lengths in equatorial cells and southward y velocities; loading converts
// them to degrees, km and northward components, the record table and the cell field alike.
TEST(PlanetGrid, Version18SaveConvertsTheTectonicRecords)
{
    constexpr int columns = 64;
    auto holder = std::make_unique<planet>();
    planet& world = *holder;
    initialiseworld(world);
    world.setwidth(columns - 1);
    TectonicBoundarySegment segment;
    segment.id = 1;
    segment.centroidLongitudeDeg = 10.0f;  // the version 18 meaning: raster x
    segment.centroidLatitudeDeg = 5.0f;    // raster y
    segment.lengthKm = 3.0f;               // equatorial cells
    world.settectonicboundarysegments(std::vector<TectonicBoundarySegment>{ segment });
    TectonicDeformingRegion region;
    region.id = 1;
    region.centroidLongitudeDeg = 40.0f;
    region.centroidLatitudeDeg = 20.0f;
    region.averageInterpolatedVelocityX = 0.5f;
    region.averageInterpolatedVelocityY = 0.25f;  // southward
    world.settectonicdeformingregions(std::vector<TectonicDeformingRegion>{ region });
    world.settectonicdeformationvelocityx(7, 3, 0.75f);
    world.settectonicdeformationvelocityy(7, 3, 0.5f);  // southward

    const std::filesystem::path path = std::filesystem::temp_directory_path() / "uw-planet-grid-v18.uww";
    world.saveworld(path.string());
    // The layout of 18 and 19 is the same: the version line makes the file a version 18 save.
    std::string text;
    {
        std::ifstream in(path, std::ios::binary);
        std::ostringstream all;
        all << in.rdbuf();
        text = all.str();
    }
    ASSERT_TRUE(text.compare(0, 2, "19") == 0 && (text[2] == '\n' || text[2] == '\r'));
    text.replace(0, 2, "18");
    {
        std::ofstream out(path, std::ios::binary);
        out << text;
    }

    auto loadedholder = std::make_unique<planet>();
    planet& loaded = *loadedholder;
    initialiseworld(loaded);
    ASSERT_TRUE(loaded.loadworld(path.string()));
    std::filesystem::remove(path);
    ASSERT_EQ(loaded.tectonicboundarysegments().size(), 1u);
    ASSERT_EQ(loaded.tectonicdeformingregions().size(), 1u);
    const TectonicBoundarySegment& s = loaded.tectonicboundarysegments()[0];
    EXPECT_NEAR(s.centroidLongitudeDeg, 10.5f / columns * 360.0f - 180.0f, 1e-4f);
    EXPECT_NEAR(s.centroidLatitudeDeg, 90.0f - 5.5f / (columns / 2) * 180.0f, 1e-4f);
    const double cellkm = 2.0 * 3.14159265358979323846 * loaded.grid().radius_m() / 1000.0 / columns;
    EXPECT_NEAR(s.lengthKm, 3.0 * cellkm, 1e-3 * cellkm);
    const TectonicDeformingRegion& r = loaded.tectonicdeformingregions()[0];
    EXPECT_NEAR(r.centroidLongitudeDeg, 40.5f / columns * 360.0f - 180.0f, 1e-4f);
    EXPECT_NEAR(r.centroidLatitudeDeg, 90.0f - 20.5f / (columns / 2) * 180.0f, 1e-4f);
    EXPECT_FLOAT_EQ(r.averageInterpolatedVelocityX, 0.5f);
    EXPECT_FLOAT_EQ(r.averageInterpolatedVelocityY, -0.25f);
    EXPECT_FLOAT_EQ(loaded.tectonicdeformationvelocityx(7, 3), 0.75f);
    EXPECT_FLOAT_EQ(loaded.tectonicdeformationvelocityy(7, 3), -0.5f);
}
