// Tests of landscape/level_blocks.hpp and the G gather of upsampling.hpp (T8 W2-W2b, sphere-mesh.md
// section 8, "Cube rules"), suite Landscape: Landscape.LevelBlocks* and Landscape.UpsampleGather*.

#include "landscape_fixture.hpp"
#include "regional_test_world.hpp"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace
{
std::uint32_t floatbits(float value)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    return bits;
}
}

TEST_F(Landscape, LevelBlocksOfALatLonDomainAreItsGlobalCells)
{
    // A window across the date line at k = 4 with a halo of 2: its blocks are the domain's global
    // cells in row-major order, each parent the global cell's id, its children the k x k block in
    // row-major order (ascending keys).
    constexpr int k = 4;
    const landscape::RegionalWindow window = landscape::regionalwindow(32, k, 2, 170.0, -10.0, -170.0, 10.0);
    const auto grid = landscape::regionalwindowgrid(window, radius_m);
    const landscape::WindowBlocks blocks = landscape::windowblocks(*grid, k);
    ASSERT_EQ(blocks.count(), static_cast<std::size_t>(window.domaincolumns) * static_cast<std::size_t>(window.domainrows));
    const int width = window.domainfinewidth();
    for (int dr = 0; dr < window.domainrows; dr++)
        for (int dc = 0; dc < window.domaincolumns; dc++)
        {
            const std::size_t b = static_cast<std::size_t>(dr) * static_cast<std::size_t>(window.domaincolumns) +
                                  static_cast<std::size_t>(dc);
            EXPECT_EQ(blocks.parents[b], window.globalcell(dc * k, dr * k));
            const mesh::CellId* children = blocks.childrenof(b);
            int n = 0;
            for (int y = dr * k; y < (dr + 1) * k; y++)
                for (int x = dc * k; x < (dc + 1) * k; x++)
                    EXPECT_EQ(children[n++], static_cast<mesh::CellId>(y * width + x));
            for (int m = 0; m < k * k; m++)
                EXPECT_EQ(landscape::parentkey(*grid, grid->globalkey(children[m]), k), blocks.parents[b]);
        }
}

TEST_F(Landscape, LevelBlocksOfACubeWindowFollowThePerFaceQuadtree)
{
    // Level 1 (factor 4) of an N 8 cube over the children of level-0 cells on three faces, one of
    // them at a cube corner: each block is a level-0 cell with its 16 children (cubechildkeys,
    // ascending), the parents in the order of their first child's local id (ascending keys here).
    constexpr int n = 8, factor = 4;
    const std::vector<mesh::LevelKey> parents{ mesh::cubekey(n, 0, 3, 4), mesh::cubekey(n, 0, 7, 7),
                                               mesh::cubekey(n, 1, 0, 7), mesh::cubekey(n, 4, 0, 0) };
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey parent : parents)
        for (const mesh::LevelKey child : mesh::cubechildkeys(n, factor, parent))
            keys.push_back(child);
    const auto grid = mesh::cubewindowgrid(n * factor, radius_m, keys);
    const landscape::WindowBlocks blocks = landscape::windowblocks(*grid, factor);
    std::vector<mesh::LevelKey> sorted = parents;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_EQ(blocks.parents, sorted);
    for (std::size_t b = 0; b < blocks.count(); b++)
    {
        const std::vector<mesh::LevelKey> expected = mesh::cubechildkeys(n, factor, blocks.parents[b]);
        EXPECT_EQ(landscape::childkeys(*grid, blocks.parents[b], factor), expected);
        for (std::size_t m = 0; m < blocks.perparent(); m++)
        {
            EXPECT_EQ(grid->globalkey(blocks.childrenof(b)[m]), expected[m]);
            EXPECT_EQ(landscape::parentkey(*grid, expected[m], factor), blocks.parents[b]);
        }
    }
    // A window missing one child of a parent is not whole parents.
    keys.pop_back();
    EXPECT_THROW(landscape::windowblocks(*mesh::cubewindowgrid(n * factor, radius_m, keys), factor), std::invalid_argument);
}

TEST_F(Landscape, UpsampleGatherEqualsBicubicOnLatLon)
{
    // On lat-lon the gather on the window grid is upsample's bicubic bit for bit: a normal, a
    // date-line and a polar window of the synthetic world at k = 4 and 5.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const double boxes[3][4] = { { 20.0, -20.0, 60.0, 10.0 }, { 170.0, -10.0, -170.0, 10.0 }, { -30.0, 70.0, 30.0, 90.0 } };
    for (const int k : { 4, 5 })
        for (const auto& box : boxes)
        {
            const landscape::RegionalWindow window = landscape::regionalwindow(32, k, 3, box[0], box[1], box[2], box[3]);
            const std::vector<float> raster = landscape::upsample(inputs.elevation_m, inputs.width, inputs.height, window);
            const std::vector<float> gathered =
                landscape::upsample(inputs.elevation_m, *inputs.grid, *landscape::regionalwindowgrid(window, radius_m));
            ASSERT_EQ(gathered.size(), raster.size());
            int different = 0;
            for (std::size_t cell = 0; cell < raster.size(); cell++)
                different += floatbits(gathered[cell]) != floatbits(raster[cell]) ? 1 : 0;
            EXPECT_EQ(different, 0) << "k " << k << ", west " << box[0];
        }
}

TEST_F(Landscape, UpsampleGatherOnTheCube)
{
    // On the cube the gather is G (mesh::cubicstencil, true-position ghosts): a constant comes back
    // exactly, and the smooth field z (the unit vector's third component) at the level-1 cells of
    // level-0 cells inside a face, at a face edge and at a cube corner is within the stencil's error.
    constexpr int n = 16, factor = 4;
    const auto levelzero = mesh::makecubegrid(n, radius_m);
    std::vector<float> constant(levelzero->cellcount(), 812.5f), smooth(levelzero->cellcount());
    for (mesh::CellId c = 0; c < levelzero->cellcount(); c++)
        smooth[c] = static_cast<float>(levelzero->centre(c).z);
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey parent : { mesh::cubekey(n, 2, 7, 8), mesh::cubekey(n, 2, 15, 8), mesh::cubekey(n, 2, 15, 15),
                                         mesh::cubekey(n, 4, 0, 0) })
        for (const mesh::LevelKey child : mesh::cubechildkeys(n, factor, parent))
            keys.push_back(child);
    const auto window = mesh::cubewindowgrid(n * factor, radius_m, keys);
    const std::vector<float> flat = landscape::upsample(constant, *levelzero, *window);
    const std::vector<float> curved = landscape::upsample(smooth, *levelzero, *window);
    double largest = 0.0;
    for (mesh::CellId c = 0; c < window->cellcount(); c++)
    {
        EXPECT_EQ(flat[c], 812.5f);
        largest = std::max(largest, std::fabs(static_cast<double>(curved[c]) - window->centre(c).z));
    }
    std::printf("[gather] cube N %d, level 1: largest error of z %.3g\n", n, largest);
    EXPECT_LT(largest, 2e-4);
}
