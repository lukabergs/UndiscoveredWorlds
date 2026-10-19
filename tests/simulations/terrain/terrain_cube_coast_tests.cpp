// The coast, shelf, ocean-feature, arid-landform and import terrain functions on a cube planet
// (T8 W2-T2; docs/architecture/sphere-mesh.md section 12): each runs on a cube of N 16 with
// synthetic inputs and leaves finite fields, and a coastline across a face edge is treated like
// one in a face interior. Seed 1 throughout.

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "core/classes.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/procedural_noise.hpp"
#include "simulations/terrain/climate_landforms.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "simulations/terrain/import_elements.hpp"
#include "simulations/terrain/terrain.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "world/planet.hpp"

namespace
{
using mesh::CellId;
namespace cg = terrain::coastgrid;

constexpr int facecells = 16;
constexpr int sealevel = 1000;
constexpr int maxelevation = 20000;

// Land (sealevel + 500) inside spherical caps of the given angular radius around the axes,
// sea (sealevel - 3000) elsewhere.
std::unique_ptr<planet> cubeplanet(const std::vector<mesh::Vec3>& axes, double radius_rad)
{
    auto world = std::make_unique<planet>();
    world->setradius(6371000.0f);
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(maxelevation);
    world->setsealevel(sealevel);
    fast_srand(1);
    const mesh::CellGrid& grid = world->grid();
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        bool land = false;
        for (const mesh::Vec3& axis : axes)
            if (mesh::dot(grid.centre(c), mesh::normalise(axis)) > std::cos(radius_rad))
                land = true;
        world->setnomfloat(c, static_cast<float>(land ? sealevel + 500 : sealevel - 3000));
    }
    return world;
}

const mesh::Vec3 faceinterior{ 1.0, 0.0, 0.0 };       // the centre of face 0
const mesh::Vec3 faceedge{ 1.0, 1.0, 0.0 };           // the middle of the edge of faces 0 and 1

int coastcells(const planet& world)
{
    int count = 0;
    for (CellId c = 0; c < world.cellcount(); c++)
        if (cg::coast(world, c))
            count++;
    return count;
}

void expectfinite(const planet& world)
{
    for (CellId c = 0; c < world.cellcount(); c++)
    {
        ASSERT_TRUE(std::isfinite(world.nomfloat(c))) << "cell " << c;
        ASSERT_TRUE(std::isfinite(world.mapfloat(c))) << "cell " << c;
        ASSERT_TRUE(std::isfinite(world.roughness(c))) << "cell " << c;
        ASSERT_TRUE(std::isfinite(world.lakesurfacefloat(c))) << "cell " << c;
    }
}

// Twelve 5 x 5 discs, the template arrays the arid landforms take.
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
}

TEST(TerrainCubeCoast, WalkHelpersCrossFaceEdgesStraight)
{
    const std::unique_ptr<planet> world = cubeplanet({ faceinterior }, 0.3);
    const mesh::CellGrid& grid = world->grid();
    // Four walks of 4 N cells east along the equator come back to the start.
    const CellId start = grid.locate(mesh::unitvector(0.01, 0.01));
    EXPECT_EQ(cg::offset(grid, start, 4 * facecells, 0), static_cast<std::int32_t>(start));
    // A walk north from the equator crosses into the north cap and comes back south on the
    // far side: 4 N steps go round a full meridian.
    EXPECT_EQ(cg::offset(grid, start, 0, -4 * facecells), static_cast<std::int32_t>(start));
    // A box of half-width 1 is the cell and its eight neighbours (seven at a cube corner).
    int cells = 0;
    cg::box(grid, start, 1, [&](CellId, int, int) { cells++; });
    EXPECT_EQ(cells, 9);
}

TEST(TerrainCubeCoast, FaceEdgeCoastlineCountsAsInterior)
{
    const double radius = 0.45;
    const std::unique_ptr<planet> interior = cubeplanet({ faceinterior }, radius);
    const std::unique_ptr<planet> edge = cubeplanet({ faceedge }, radius);
    const int interiorcount = coastcells(*interior);
    const int edgecount = coastcells(*edge);
    ASSERT_GT(interiorcount, 0);
    std::cout << "coast cells: interior " << interiorcount << ", across the face edge " << edgecount << "\n";
    EXPECT_NEAR(static_cast<double>(edgecount), static_cast<double>(interiorcount), 0.10 * interiorcount)
        << "interior " << interiorcount << ", across the face edge " << edgecount;

    // The outline (land side) too.
    int interioroutline = 0, edgeoutline = 0;
    for (CellId c = 0; c < interior->cellcount(); c++)
    {
        interioroutline += cg::outline(*interior, c) ? 1 : 0;
        edgeoutline += cg::outline(*edge, c) ? 1 : 0;
    }
    EXPECT_NEAR(static_cast<double>(edgeoutline), static_cast<double>(interioroutline), 0.10 * interioroutline);
}

TEST(TerrainCubeCoast, CoastCleanupsRunAndStayFinite)
{
    const std::unique_ptr<planet> world = cubeplanet({ faceinterior, mesh::Vec3{ 0.0, 0.0, 1.0 } }, 0.55);
    const int before = coastcells(*world);

    removesmallseas(*world, 20, sealevel + 1);
    widenchannels(*world);
    loweroceans(*world);
    removestraights(*world);
    normalisecoasts(*world, 13, 11, 4);
    checkislands(*world);
    extendnoshade(*world);
    expectfinite(*world);

    // The caps keep a coast of the same order.
    const int after = coastcells(*world);
    EXPECT_GT(after, before / 2);
    EXPECT_LT(after, before * 2);

    // The queries.
    const mesh::CellGrid& grid = world->grid();
    const CellId inland = grid.locate(mesh::normalise(faceinterior));
    const CellId ocean = grid.locate(mesh::unitvector(-0.6, 2.5));
    EXPECT_EQ(world->sea(inland), 0);
    EXPECT_EQ(world->sea(ocean), 1);
    EXPECT_FALSE(cg::vaguelycoastal(*world, inland));
    const std::int32_t sea = cg::nearestsea(*world, inland, false, 5000.0);
    ASSERT_GE(sea, 0);
    EXPECT_EQ(world->sea(static_cast<CellId>(sea)), 1);
    EXPECT_GE(cg::landdistance(*world, inland), -1);
}

TEST(TerrainCubeCoast, ShelvesTrenchesAndRidgesOnTheCube)
{
    const std::unique_ptr<planet> world = cubeplanet({ faceinterior, faceedge }, 0.4);
    const mesh::CellGrid& grid = world->grid();

    std::vector<int> voronoi;
    makeshelvesvoronoi(*world, voronoi, 4);
    ASSERT_EQ(voronoi.size(), world->cellcount());
    for (const int id : voronoi)
        EXPECT_GT(id, 0);

    generateimportedshelves(*world);
    const worldgrid<bool>& shelves = world->continentalshelves();
    int shelfcells = 0;
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        if (world->sea(c) == 0)
            EXPECT_TRUE(shelves[c]) << "land is shelf";
        shelfcells += shelves[c] ? 1 : 0;
    }
    EXPECT_GT(shelfcells, 0);
    EXPECT_LT(static_cast<std::size_t>(shelfcells), world->cellcount());
    expectfinite(*world);

    generateimportedoceanridges(*world);
    generateimportedseabed(*world);
    expectfinite(*world);
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        EXPECT_GE(world->oceanridgeangle(c), 0);
        EXPECT_LE(world->oceanridgeangle(c), 360);
    }
    (void)grid;
}

TEST(TerrainCubeCoast, AridLandformsAndRoughnessOnTheCube)
{
    const std::unique_ptr<planet> world = cubeplanet({ faceinterior, faceedge }, 0.6);
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        world->setclimate(c, world->sea(c) ? 0 : 5); // hot desert on all land
        world->setroughness(c, 1000.0f);
    }
    Templates templates;
    createergs(*world, templates.shapes, templates.shapes, templates.shapes);
    createsaltpans(*world, templates.shapes, templates.shapes);
    refineroughnessmap(*world);
    expectfinite(*world);

    int specials = 0;
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        if (world->special(c) == 110 || world->special(c) == 120)
        {
            specials++;
            EXPECT_EQ(world->sea(c), 0);
            EXPECT_GT(world->lakesurfacefloat(c), static_cast<float>(sealevel));
        }
        EXPECT_GE(world->roughness(c), 0.0f);
    }
    EXPECT_GT(specials, 0);
}

TEST(TerrainCubeCoast, ImportedMountainsOnTheCube)
{
    const std::unique_ptr<planet> world = cubeplanet({ faceinterior, faceedge }, 0.6);
    const mesh::CellGrid& grid = world->grid();
    std::vector<std::vector<int>> raw(ARRAYWIDTH, std::vector<int>(ARRAYHEIGHT, 0));
    std::vector<std::vector<bool>> okmountains(ARRAYWIDTH, std::vector<bool>(ARRAYHEIGHT, false));

    // A ridge along the equator across the face edge of faces 0 and 1.
    int ridgecells = 0;
    for (double longitude = -0.4; longitude <= 1.9; longitude += 0.02)
    {
        const CellId c = grid.locate(mesh::unitvector(0.0, longitude));
        int x = 0, y = 0;
        legacyslot(grid, c, x, y);
        if (raw[static_cast<std::size_t>(x)][static_cast<std::size_t>(y)] == 0)
            ridgecells++;
        raw[static_cast<std::size_t>(x)][static_cast<std::size_t>(y)] = 2000;
    }
    ASSERT_GT(ridgecells, 10);

    createmountainsfromraw(*world, raw, okmountains);
    expectfinite(*world);

    int mountains = 0;
    for (CellId c = 0; c < world->cellcount(); c++)
    {
        if (world->mountainridge(c) == 0)
            continue;
        mountains++;
        EXPECT_GE(world->mountainheightfloat(c), 0.0f);
        // Every ridge has its other end at the neighbour.
        for (int slot = 1; slot <= 8; slot++)
        {
            if ((world->mountainridge(c) & (1 << (slot - 1))) == 0)
                continue;
            const std::int32_t n = grid.neighbour(c, slot);
            ASSERT_GE(n, 0);
        }
        int x = 0, y = 0;
        legacyslot(grid, c, x, y);
        EXPECT_TRUE(okmountains[static_cast<std::size_t>(x)][static_cast<std::size_t>(y)]);
    }
    EXPECT_GT(mountains, 0);
}
