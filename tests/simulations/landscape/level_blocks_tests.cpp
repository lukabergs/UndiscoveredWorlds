// Tests of landscape/level_blocks.hpp and the G gather of upsampling.hpp (T8 W2-W2b, sphere-mesh.md
// section 8, "Cube rules"), suite Landscape: Landscape.LevelBlocks* and Landscape.UpsampleGather*.

#include "landscape_fixture.hpp"
#include "regional_test_world.hpp"

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/refinement_projection.hpp"
#include "simulations/landscape/regional_surface.hpp"
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

namespace
{
// The keys of the cube level-N cells within `rings` rings of a cell (on any faces), ascending.
std::vector<mesh::LevelKey> cubedisc(int n, mesh::LevelKey centre, int rings)
{
    const auto grid = mesh::makecubegrid(n, radius_m);
    const std::vector<int> ring = mesh::ringdistances(*grid, { static_cast<mesh::CellId>(centre) }, rings);
    std::vector<mesh::LevelKey> keys;
    for (mesh::CellId c = 0; c < grid->cellcount(); c++)
        if (ring[c] >= 0)
            keys.push_back(c);
    return keys;
}

// The children at factor of a set of parent keys at N, ascending.
std::vector<mesh::LevelKey> cubechildren(int n, int factor, const std::vector<mesh::LevelKey>& parents)
{
    std::vector<mesh::LevelKey> keys;
    for (const mesh::LevelKey parent : parents)
        for (const mesh::LevelKey child : mesh::cubechildkeys(n, factor, parent))
            keys.push_back(child);
    std::sort(keys.begin(), keys.end());
    return keys;
}

// The largest |block mean - level-0 value| over the blocks whose parent lies in `inner`.
double largestblockerror(const landscape::WindowBlocks& blocks, const std::vector<float>& fine, const std::vector<float>& field,
                         const std::vector<mesh::LevelKey>& inner)
{
    double largest = 0.0;
    for (std::size_t b = 0; b < blocks.count(); b++)
        if (std::binary_search(inner.begin(), inner.end(), blocks.parents[b]))
            largest = std::max(largest, std::fabs(landscape::blockmean(blocks, b, fine) -
                                                  static_cast<double>(field[blocks.parents[b]])));
    return largest;
}

// A deterministic value in [-1, 1) per key.
double jitter(mesh::LevelKey key)
{
    std::uint64_t h = key * 0x9e3779b97f4a7c15ull;
    h ^= h >> 29;
    h *= 0xbf58476d1ce4e5b9ull;
    h ^= h >> 32;
    return static_cast<double>(h % 2000u) / 1000.0 - 1.0;
}
}

TEST_F(Landscape, LevelBlocksNearestDomainKeyByRings)
{
    // R for nodes outside a domain: the domain key at the smallest ring distance (brute force by
    // ringdistances on the whole N 8 cube), ties to the smallest key; a domain key is itself.
    constexpr int n = 8;
    const auto grid = mesh::makecubegrid(n, radius_m);
    const std::vector<mesh::LevelKey> domain = cubedisc(n, mesh::cubekey(n, 1, 7, 0), 2);  // a cube corner
    const landscape::NearestDomainKey nearest(n, domain);
    int checked = 0;
    for (mesh::CellId key = 0; key < grid->cellcount(); key++)
    {
        const mesh::LevelKey got = domain[nearest(key)];
        if (std::binary_search(domain.begin(), domain.end(), mesh::LevelKey{ key }))
        {
            EXPECT_EQ(got, key);
            continue;
        }
        const std::vector<int> ring = mesh::ringdistances(*grid, { key }, 4 * n);
        int best = -1;
        mesh::LevelKey expected = mesh::nokey;
        for (const mesh::LevelKey cell : domain)
        {
            const int d = ring[static_cast<std::size_t>(cell)];
            if (d >= 0 && (best < 0 || d < best || (d == best && cell < expected)))
            {
                best = d;
                expected = cell;
            }
        }
        EXPECT_EQ(got, expected) << "key " << key;
        checked++;
        if (checked > 60)
            break;
    }
}

TEST_F(Landscape, LevelBlocksCubeBlockMeansConverge)
{
    // correctblockmeans on a cube window around a cube corner (level 0 N 8, k 4): a constant offset
    // of the gathered field is removed exactly (G of a constant residual is the constant), and
    // +-10 m of per-cell noise on the gathered field falls to a fraction after two passes over the
    // blocks two rings inside the domain (the gather reaches two level-0 rings).
    constexpr int n = 8, k = 4;
    const auto levelzero = mesh::makecubegrid(n, radius_m);
    std::vector<float> field(levelzero->cellcount());
    for (mesh::CellId c = 0; c < levelzero->cellcount(); c++)
        field[c] = static_cast<float>(1000.0 * levelzero->centre(c).x + 300.0 * levelzero->centre(c).z);
    const mesh::LevelKey corner = mesh::cubekey(n, 1, 7, 0);
    const auto window = mesh::cubewindowgrid(n * k, radius_m, cubechildren(n, k, cubedisc(n, corner, 4)));
    const std::vector<mesh::LevelKey> inner = cubedisc(n, corner, 2);
    const landscape::WindowBlocks blocks = landscape::windowblocks(*window, k);
    const std::vector<float> gathered = landscape::upsample(field, *levelzero, *window);

    std::vector<float> offset = gathered, plain = gathered;
    const double before = largestblockerror(blocks, gathered, field, inner);
    for (float& value : offset)
        value += 25.0f;
    landscape::correctblockmeans(field, *levelzero, *window, 1, offset);
    landscape::correctblockmeans(field, *levelzero, *window, 1, plain);
    double largest = 0.0;
    for (std::size_t cell = 0; cell < offset.size(); cell++)
        largest = std::max(largest, std::fabs(static_cast<double>(offset[cell]) - static_cast<double>(plain[cell])));
    EXPECT_LT(largest, 1e-3) << "the offset is removed";

    std::vector<float> noisy = gathered;
    for (mesh::CellId c = 0; c < window->cellcount(); c++)
        noisy[c] += static_cast<float>(10.0 * jitter(window->globalkey(c)));
    const double noisybefore = largestblockerror(blocks, noisy, field, inner);
    landscape::correctblockmeans(field, *levelzero, *window, 2, noisy);
    const double after = largestblockerror(blocks, noisy, field, inner);
    std::printf("[blocks] cube corner: gathered %.3g m, noisy %.3g -> %.3g m after two passes\n", before, noisybefore, after);
    EXPECT_LT(after, 0.2 * noisybefore);
}

TEST_F(Landscape, LevelBlocksCubeProlongAndParentMeans)
{
    // prolong on cube windows (N 8 -> 32, factor 4) across a face edge: a constant parent comes back
    // exactly; with the limit every child lies in the range of the central 2 x 2 nodes, and a step
    // parent field clamps some children; on a smooth parent field with +-20 m of per-cell noise
    // correctparentmeans brings the block means of the inner
    // parents near the parent values.
    constexpr int n = 8, factor = 4;
    const mesh::LevelKey edge = mesh::cubekey(n, 0, 7, 3);
    const std::vector<mesh::LevelKey> parentkeys = cubedisc(n, edge, 7);
    const auto parent = mesh::cubewindowgrid(n, radius_m, parentkeys);
    const auto level = mesh::cubewindowgrid(n * factor, radius_m, cubechildren(n, factor, cubedisc(n, edge, 3)));
    std::vector<float> flat(parent->cellcount(), 140.0f), step(parent->cellcount());
    for (mesh::CellId c = 0; c < parent->cellcount(); c++)
        step[c] = parent->centre(c).y > 0.2 ? 900.0f : 100.0f;
    const landscape::Prolongation constant = landscape::prolong(*parent, flat, *level, factor, true);
    for (const float value : constant.values)
        EXPECT_EQ(value, 140.0f);
    const landscape::Prolongation plain = landscape::prolong(*parent, step, *level, factor, false);
    const landscape::Prolongation limited = landscape::prolong(*parent, step, *level, factor, true);
    // The central nodes are cells or, beyond a face edge, ghosts (1-D Catmull-Rom values, which may
    // themselves overshoot a step: the rule of section 8).
    const auto nodevalue = [&](mesh::LevelKey key) { return step[parent->localofkey(key)]; };
    int overshoot = 0;
    for (mesh::CellId cell = 0; cell < level->cellcount(); cell++)
    {
        overshoot += plain.values[cell] > 900.0f || plain.values[cell] < 100.0f ? 1 : 0;
        double lowest = 1e30, highest = -1e30;
        for (int row = 1; row <= 2; row++)
            for (int column = 1; column <= 2; column++)
            {
                const double node =
                    mesh::applystencil(mesh::cubicstencilnode(n, level->centre(cell), column, row), nodevalue);
                lowest = std::min(lowest, node);
                highest = std::max(highest, node);
            }
        EXPECT_GE(limited.values[cell], static_cast<float>(lowest)) << "cell " << cell;
        EXPECT_LE(limited.values[cell], static_cast<float>(highest)) << "cell " << cell;
    }
    EXPECT_GT(overshoot, 0);
    EXPECT_GT(limited.clamped, 0);

    std::vector<float> smooth(parent->cellcount());
    for (mesh::CellId c = 0; c < parent->cellcount(); c++)
        smooth[c] = static_cast<float>(1000.0 * parent->centre(c).x + 300.0 * parent->centre(c).z);
    std::vector<float> fine = landscape::prolong(*parent, smooth, *level, factor, false).values;
    for (mesh::CellId c = 0; c < level->cellcount(); c++)
        fine[c] += static_cast<float>(20.0 * jitter(level->globalkey(c)));
    const landscape::WindowBlocks blocks = landscape::windowblocks(*level, factor);
    const std::vector<mesh::LevelKey> inner = cubedisc(n, edge, 1);
    std::vector<float> parentbykey(static_cast<std::size_t>(6 * n * n), 0.0f);
    for (mesh::CellId c = 0; c < parent->cellcount(); c++)
        parentbykey[parent->globalkey(c)] = smooth[c];
    const double before = largestblockerror(blocks, fine, parentbykey, inner);
    landscape::correctparentmeans(*parent, smooth, *level, factor, 2, fine);
    const double after = largestblockerror(blocks, fine, parentbykey, inner);
    std::printf("[blocks] cube edge: parent-mean error %.3g -> %.3g m after two passes\n", before, after);
    EXPECT_LT(after, 0.5 * before);
    EXPECT_THROW(landscape::prolong(*level, std::vector<float>(level->cellcount()), *level, factor, false),
                 std::invalid_argument);
}

TEST_F(Landscape, LevelBlocksCubeSmoothingAndSignedMeans)
{
    // smoothparentmeans on cube windows around a cube corner moves most of every kept parent's
    // lambda (rate 1, room 1 km); the signed block-mean correction keeps every cell on its side of the
    // sea level and lowers the block errors.
    constexpr int n = 8, factor = 4;
    const mesh::LevelKey corner = mesh::cubekey(n, 4, 0, 0);
    const auto parent = mesh::cubewindowgrid(n, radius_m, cubedisc(n, corner, 7));
    const auto level = mesh::cubewindowgrid(n * factor, radius_m, cubechildren(n, factor, cubedisc(n, corner, 3)));
    std::vector<float> parentfield(parent->cellcount());
    for (mesh::CellId c = 0; c < parent->cellcount(); c++)
        parentfield[c] = static_cast<float>(500.0 + 200.0 * jitter(parent->globalkey(c)));
    std::vector<float> field = landscape::prolong(*parent, parentfield, *level, factor, false).values;
    landscape::ProjectionBounds bounds;
    bounds.rate.assign(level->cellcount(), 1.0f);
    bounds.up_m.assign(level->cellcount(), 1000.0f);
    bounds.down_m.assign(level->cellcount(), 1000.0f);
    const std::vector<mesh::LevelKey> inner = cubedisc(n, corner, 1);
    std::vector<std::uint8_t> kept(level->cellcount(), 0);
    for (mesh::CellId c = 0; c < level->cellcount(); c++)
        kept[c] = std::binary_search(inner.begin(), inner.end(), mesh::cubeparentkey(n * factor, factor, level->globalkey(c)))
                      ? 1 : 0;
    const landscape::SmoothingReport report =
        landscape::smoothparentmeans(*level, field, bounds, *parent, parentfield, factor, 6, 1.0, 0, kept);
    std::printf("[blocks] cube corner smoothing: lambda %.3g -> %.3g m, step %.3g -> %.3g m\n", report.maxlambdabefore_m,
                report.maxlambdaafter_m, report.maxlambdastepbefore_m, report.maxlambdastepafter_m);
    EXPECT_LT(report.maxlambdaafter_m, 0.5 * report.maxlambdabefore_m);
    EXPECT_THROW(landscape::smoothparentmeans(*level, field, bounds, *parent, parentfield, factor, 1, 1.0, 1, kept),
                 std::invalid_argument);

    RegionalGlobalInputs inputs;
    inputs.grid = mesh::makecubegrid(n, radius_m);
    inputs.sealevel = 0.0f;
    inputs.elevation_m.resize(inputs.grid->cellcount());
    for (mesh::CellId c = 0; c < inputs.grid->cellcount(); c++)
        inputs.elevation_m[c] = static_cast<float>(400.0 * inputs.grid->centre(c).x - 100.0);
    std::vector<float> surface = landscape::upsample(inputs.elevation_m, *inputs.grid, *level);
    for (mesh::CellId c = 0; c < level->cellcount(); c++)
        surface[c] += static_cast<float>(30.0 * jitter(level->globalkey(c)));
    const std::vector<float> start = surface;
    const landscape::WindowBlocks blocks = landscape::windowblocks(*level, factor);
    const double before = largestblockerror(blocks, surface, inputs.elevation_m, inner);
    landscape::correctsignedblockmeans(inputs, *level, surface);
    const double after = largestblockerror(blocks, surface, inputs.elevation_m, inner);
    std::printf("[blocks] cube corner signed means: %.3g -> %.3g m\n", before, after);
    EXPECT_LT(after, before);
    for (std::size_t cell = 0; cell < surface.size(); cell++)
        EXPECT_EQ(surface[cell] > 0.0f, start[cell] > 0.0f) << "cell " << cell;
}
