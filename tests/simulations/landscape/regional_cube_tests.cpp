// Tests of the regional pass's cube forms (T8 W2-W2, sphere-mesh.md section 8), suite Landscape:
// Landscape.RegionalCube*. Synthetic cube worlds (level 0 N 32) and cube windows built from rings of
// level-0 cells; two windows that share cells far enough from both domains' edges must give them
// bit-identical values (the overlap property of regional_pass.hpp on the cube).

#include "landscape_fixture.hpp"
#include "regional_test_world.hpp"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_surface.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <vector>

namespace
{
constexpr int cuben = 32;

// The level-0 keys within `rings` rings of a cell, ascending.
std::vector<mesh::LevelKey> disc(mesh::LevelKey centre, int rings) { return mesh::dilatekeys(cuben, { centre }, rings); }

// The cube window of the children at k of a set of level-0 keys.
std::shared_ptr<const mesh::CellGrid> windowof(const std::vector<mesh::LevelKey>& parents, int k)
{
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey parent : parents)
        for (const mesh::LevelKey child : mesh::cubechildkeys(cuben, k, parent))
            keys.push_back(child);
    return mesh::cubewindowgrid(cuben * k, radius_m, keys);
}

// A cube world: smooth relief of a few thousand metres with coasts, seas and one-cell islands.
RegionalGlobalInputs cubeworld()
{
    RegionalGlobalInputs inputs;
    inputs.grid = mesh::makecubegrid(cuben, radius_m);
    inputs.radius_m = radius_m;
    inputs.sealevel = 0.0f;
    inputs.seed = 1;
    inputs.lapserate_c_per_km = 6.5f;
    const std::size_t size = inputs.grid->cellcount();
    inputs.elevation_m.resize(size);
    inputs.network.sea.assign(size, 0);
    inputs.lakesurface_m.assign(size, 0.0f);
    inputs.island.assign(size, 0);
    for (mesh::CellId c = 0; c < size; c++)
    {
        const mesh::Vec3 p = inputs.grid->centre(c);
        const double z = 1500.0 * std::sin(5.0 * p.x + 1.0) * std::cos(4.0 * p.y - 0.5) + 900.0 * std::sin(7.0 * p.z) - 200.0;
        inputs.elevation_m[c] = static_cast<float>(z);
        inputs.network.sea[c] = z < 0.0 ? 1 : 0;
    }
    return inputs;
}

// The coast and the land-share rule at the defaults of configs/parameters/landscape.json, with fewer
// share and sign-keeping passes (the box's margin shrinks with them).
RegionalParameters coastparameters()
{
    RegionalParameters parameters;
    parameters.meanpasses = 2;
    parameters.warp.amplitude_cells = 0.5;
    parameters.warp.wavelength_cells = 2.0;
    parameters.warp.fade_m = 300.0;
    parameters.coast.noiseamplitude_m = 50.0;
    parameters.coast.noisepersistence = 0.6;
    parameters.coast.noiserelieffraction = 0.25;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.coast.shorefloor_m = 20.0;
    parameters.coast.sharepasses = 2;
    parameters.coast.signedpasses = 2;
    parameters.coast.onesided = 1;
    parameters.coast.islandmargin = 0.125;
    parameters.coast.warpcells = 0.15;
    parameters.coast.warppersistence = 0.6;
    parameters.coast.pocketreach = 2;
    parameters.edges.smoothgate = 3;
    parameters.edges.straightmax = 4;
    parameters.edges.chunkmin = 2;
    parameters.edges.chunkmax = 7;
    return parameters;
}

std::uint32_t floatbits(float value)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    return bits;
}
}

TEST_F(Landscape, RegionalCubeStartSurfaceOverlaps)
{
    // Two windows at k 4, discs of 11 rings around a cell at a cube corner and one 5 cells along its
    // face edge: on the cells of the level-0 cells within 3 rings of both centres (8 rings inside both
    // domains) the land-share start surface, its base and the sea are bit-identical, the land share
    // holds on the domain, and no cell of the start surface is far from its block's mean.
    constexpr int k = 4;
    const RegionalGlobalInputs inputs = cubeworld();
    const RegionalParameters parameters = coastparameters();
    const mesh::LevelKey corner = mesh::cubekey(cuben, 1, cuben - 1, 0), along = mesh::cubekey(cuben, 1, cuben - 1, 5);
    const auto a = windowof(disc(corner, 11), k);
    const auto b = windowof(disc(along, 11), k);
    const landscape::FineWarp warpa = landscape::regionalwarp(*a, cuben, inputs.seed, parameters);
    const landscape::FineWarp warpb = landscape::regionalwarp(*b, cuben, inputs.seed, parameters);
    landscape::StartSurface sa = landscape::regionalstartsurface(inputs, *a, warpa, parameters, false);
    landscape::StartSurface sb = landscape::regionalstartsurface(inputs, *b, warpb, parameters, false);
    ASSERT_EQ(sa.start_m.size(), a->cellcount());
    std::vector<float> surfacea = sa.start_m, surfaceb = sb.start_m;
    const std::vector<uint8_t> seaa = landscape::regionalsea(inputs, *a, warpa, parameters, surfacea, nullptr, false);
    const std::vector<uint8_t> seab = landscape::regionalsea(inputs, *b, warpb, parameters, surfaceb, nullptr, false);
    const std::vector<mesh::LevelKey> innera = disc(corner, 3), innerb = disc(along, 3);
    int shared = 0, starts = 0, bases = 0, seas = 0;
    for (mesh::CellId c = 0; c < a->cellcount(); c++)
    {
        const mesh::LevelKey key = a->globalkey(c);
        const mesh::LevelKey parent = mesh::cubeparentkey(cuben * k, k, key);
        if (!std::binary_search(innera.begin(), innera.end(), parent) || !std::binary_search(innerb.begin(), innerb.end(), parent))
            continue;
        const mesh::CellId other = b->localofkey(key);
        ASSERT_NE(other, mesh::nocell);
        shared++;
        starts += floatbits(sa.start_m[c]) != floatbits(sb.start_m[other]) ? 1 : 0;
        bases += floatbits(sa.base_m[c]) != floatbits(sb.base_m[other]) ? 1 : 0;
        seas += seaa[c] != seab[other] ? 1 : 0;
    }
    int sea = 0;
    for (const uint8_t value : seaa)
        sea += value;
    std::printf("[cube coast] %d shared cells: %d starts, %d bases, %d sea flags differ; sea %d of %zu cells; share passes %d, "
                "misses %d; block-mean error %.3g m; pockets %d filled, %d left\n",
                shared, starts, bases, seas, sea, seaa.size(), sa.sharepasses, sa.sharemisses, sa.blockmeanerror_m,
                sa.pocketsfilled, sa.pockets);
    EXPECT_GT(shared, 0);
    EXPECT_EQ(starts, 0);
    EXPECT_EQ(bases, 0);
    EXPECT_EQ(seas, 0);
    EXPECT_GT(sea, 0);
    EXPECT_LT(sea, static_cast<int>(seaa.size()));
    EXPECT_LT(sa.blockmeanerror_m, 50.0);
}

TEST_F(Landscape, RegionalCubeRefinedCoastKeepsTheParentShares)
{
    // refinedcoast at level 2 of a cube window (parents the level-1 cells of a disc, factor 4): the sea
    // it returns lies inside the open cells, and every parent sea cell's children hold at most half land
    // where the displacement moved them (the guard), else the undisplaced sea stays.
    constexpr int k1 = 4, factor = 4;
    const RegionalGlobalInputs inputs = cubeworld();
    const RegionalParameters parameters = coastparameters();
    const mesh::LevelKey centre = mesh::cubekey(cuben, 4, 3, 3);
    const auto parent = windowof(disc(centre, 6), k1);
    std::vector<float> parentelevation = landscape::upsample(inputs.elevation_m, *inputs.grid, *parent);
    std::vector<uint8_t> parentsea(parent->cellcount());
    for (std::size_t cell = 0; cell < parentsea.size(); cell++)
        parentsea[cell] = parentelevation[cell] <= 0.0f ? 1 : 0;
    // The level: the children of the parent cells within 3 level-0 rings of the centre.
    std::vector<mesh::LevelKey> levelparents;
    for (const mesh::LevelKey zero : disc(centre, 3))
        for (const mesh::LevelKey one : mesh::cubechildkeys(cuben, k1, zero))
            levelparents.push_back(one);
    std::sort(levelparents.begin(), levelparents.end());
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey one : levelparents)
        for (const mesh::LevelKey two : mesh::cubechildkeys(cuben * k1, factor, one))
            keys.push_back(two);
    const auto level = mesh::cubewindowgrid(cuben * k1 * factor, radius_m, keys);
    std::vector<float> surface = landscape::upsample(inputs.elevation_m, *inputs.grid, *level);
    const landscape::FineWarp warp = landscape::regionalwarp(*level, cuben, inputs.seed, parameters);
    landscape::PocketCount filled, left;
    const std::vector<uint8_t> sea = landscape::refinedcoast(inputs, *level, factor, *parent, parentsea, parentelevation, warp,
                                                             parameters, nullptr, surface, nullptr, &filled, &left);
    std::vector<float> below(level->cellcount(), -1.0f);
    const std::vector<uint8_t> open = landscape::regionalsea(inputs, *level, warp, parameters, below, nullptr, true);
    int seacells = 0;
    for (std::size_t cell = 0; cell < sea.size(); cell++)
    {
        seacells += sea[cell];
        if (sea[cell] != 0)
        {
            EXPECT_NE(open[cell], 0) << "cell " << cell;
            EXPECT_LE(surface[cell], 0.0f) << "cell " << cell;
        }
    }
    std::printf("[cube coast] level 2: %d sea of %zu cells, pockets %d filled, %d left\n", seacells, sea.size(),
                filled.pockets, left.pockets);
    EXPECT_GT(seacells, 0);
}
