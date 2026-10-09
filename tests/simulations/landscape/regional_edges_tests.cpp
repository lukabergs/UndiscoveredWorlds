// Tests of regional_edges.hpp (erosion plan Phase R2: the R6 contract; R8's smooth gate and
// straight-run breaker, and the RegionalOverlapEdges* variants with both on, on a world wider
// than the land-share rule's widened box), suite Landscape: TEST_F(Landscape,
// RegionalEdges<What>) on the synthetic world of regional_test_world.hpp.

#include <algorithm>

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_surface.hpp"

namespace
{
using landscape::FineWarp;

// The R2 coast (regional_surface_tests.cpp's coastparameters) with R8's features on (the
// smooth gate of reach 3, R8's default).
RegionalParameters edgeparameters(int straightmax = 4)
{
    RegionalParameters parameters;
    parameters.warp.amplitude_cells = 0.5;
    parameters.warp.wavelength_cells = 2.0;
    parameters.warp.fade_m = 400.0;
    parameters.coast.noiseamplitude_m = 150.0;
    parameters.coast.noisepersistence = 0.5;
    parameters.coast.noiserelieffraction = 0.25;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.edges.smoothgate = 3;
    parameters.edges.straightmax = straightmax;
    return parameters;
}

bool lakeglobal(const RegionalGlobalInputs& inputs, std::size_t global)
{
    return inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
}

// The longest straight shore run inside one block of the box (the breaker's definition:
// sea cells along a row or column of a block with land on the same side; a side beyond the
// box is not land), over the blocks that are not lake cells.
int longestblockrun(const RegionalGlobalInputs& inputs, const RegionalWindow& box, const std::vector<uint8_t>& mask)
{
    const int k = box.scale, finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    auto kind = [&](int x, int y) {
        if (x < 0 || y < 0 || x >= finewidth || y >= fineheight)
            return -1;
        return mask[at(x, y, finewidth)] != 0 ? 1 : 0;
    };
    const int sidex[4] = {0, 0, -1, 1}, sidey[4] = {-1, 1, 0, 0};
    int longest = 0;
    for (int by = 0; by < box.domainrows; by++)
        for (int bx = 0; bx < box.domaincolumns; bx++)
        {
            if (lakeglobal(inputs, at((box.domaincolumn0 + bx) % inputs.width, box.domainrow0 + by, inputs.width)))
                continue;
            for (int side = 0; side < 4; side++)
                for (int line = 0; line < k; line++)
                {
                    int run = 0;
                    for (int t = 0; t < k; t++)
                    {
                        const int x = side < 2 ? bx * k + t : bx * k + line, y = side < 2 ? by * k + line : by * k + t;
                        run = kind(x, y) == 1 && kind(x + sidex[side], y + sidey[side]) == 0 ? run + 1 : 0;
                        longest = std::max(longest, run);
                    }
                }
        }
    return longest;
}

// A global cell of the box's coast: land with a sea cell among its 8 neighbours, not a lake.
bool coastalland(const RegionalGlobalInputs& inputs, int i, int j)
{
    if (inputs.network.sea[at(i, j, inputs.width)] != 0 || j <= 0 || j + 1 >= inputs.height)
        return false;
    for (int dj = -1; dj <= 1; dj++)
        for (int di = -1; di <= 1; di++)
            if (inputs.network.sea[at((i + di + inputs.width) % inputs.width, j + dj, inputs.width)] != 0)
                return true;
    return false;
}
}

TEST_F(Landscape, RegionalEdgesGateIsTheThreeByThreeRuleByDefault)
{
    // With edges.smoothgate 0 the gate of every fine cell is nearseacells of its global cell,
    // and the breaker changes nothing at straightmax 0.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const RegionalParameters parameters;
    const std::vector<uint8_t> gate = landscape::seagate(inputs, window, FineWarp{}, parameters);
    const std::vector<uint8_t> nearsea = landscape::nearseacells(inputs);
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    ASSERT_EQ(gate.size(), static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight));
    int gated = 0, mismatches = 0;
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const uint8_t expected = nearsea[at(window.globalcolumn(x), window.globalrow(y), inputs.width)];
            if (gate[at(x, y, domainwidth)] != expected)
                mismatches++;
            gated += gate[at(x, y, domainwidth)] != 0 ? 1 : 0;
        }
    EXPECT_EQ(mismatches, 0);
    EXPECT_GT(gated, 0);
    EXPECT_LT(gated, domainwidth * domainheight);

    std::vector<uint8_t> mask(gate.size(), 0);
    for (std::size_t cell = 0; cell < mask.size(); cell++)
        mask[cell] = (cell % 7 == 0) ? 1 : 0;
    const std::vector<uint8_t> before = mask;
    landscape::breakstraightruns(inputs, window, parameters, gate, mask);
    EXPECT_EQ(mask, before);
}

TEST_F(Landscape, RegionalEdgesSmoothGateCoversItsReach)
{
    // edges.smoothgate n without a warp: 1 where a global sea cell's centre lies within n global
    // cells of the fine centre in each axis (half-open: -n <= d < n; n = 2 is bicubic's
    // stencil), 0 beyond and in lake cells; it contains the 3 x 3 rule and reaches into cells
    // that rule closes.
    const RegionalWindow window = windowc(8);  // the 64 x 32 world's window, short of the poles
    const RegionalGlobalInputs world = syntheticworld(64);
    const std::vector<uint8_t> nearsea = landscape::nearseacells(world);
    const int k = window.scale, domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    for (int reach : {2, 3})
    {
        SCOPED_TRACE("reach " + std::to_string(reach));
        RegionalParameters parameters;
        parameters.edges.smoothgate = reach;
        const std::vector<uint8_t> gate = landscape::seagate(world, window, FineWarp{}, parameters);
        int beyond3x3 = 0, missing3x3 = 0, wrong = 0;
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t global = at(window.globalcolumn(x), window.globalrow(y), world.width);
                const uint8_t value = gate[at(x, y, domainwidth)];
                if (nearsea[global] != 0)
                {
                    missing3x3 += value == 0 ? 1 : 0;
                    continue;
                }
                const double gx = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / k;
                const double gy = (static_cast<double>(window.absolutefinerow(y)) + 0.5) / k;
                bool within = false;
                for (int j = 0; j < world.height; j++)
                    for (int i = 0; i < world.width; i++)
                    {
                        if (world.network.sea[at(i, j, world.width)] == 0)
                            continue;
                        double dx = gx - (i + 0.5);
                        dx -= world.width * std::floor((dx + world.width / 2.0) / world.width);  // columns wrap
                        const double dy = gy - (j + 0.5);
                        within = within || (dx >= -reach && dx < reach && dy >= -reach && dy < reach);
                    }
                if (lakeglobal(world, global))
                    within = false;
                wrong += (value != 0) != within ? 1 : 0;
                beyond3x3 += value != 0 ? 1 : 0;
            }
        EXPECT_EQ(missing3x3, 0);
        EXPECT_EQ(wrong, 0);
        EXPECT_GT(beyond3x3, 0) << "the smooth gate must reach beyond the 3 x 3 rule somewhere";
    }

    // With the warp the 3 x 3 rule is still contained.
    const RegionalParameters warped = edgeparameters(0);
    const FineWarp warp = landscape::regionalwarp(window, world.seed, warped);
    ASSERT_FALSE(warp.empty());
    const std::vector<uint8_t> warpedgate = landscape::seagate(world, window, warp, warped);
    int warpedmissing = 0;
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
            if (nearsea[at(window.globalcolumn(x), window.globalrow(y), world.width)] != 0 && warpedgate[at(x, y, domainwidth)] == 0)
                warpedmissing++;
    EXPECT_EQ(warpedmissing, 0);
}

TEST_F(Landscape, RegionalEdgesInlandCellsFollowTheReach)
{
    // Reach 0 to 2: 1 - nearseacells (the R2 lift); reach n >= 3: 1 where no sea cell lies within
    // n - 1 cells (Chebyshev) or in a lake cell.
    RegionalGlobalInputs world = syntheticworld(64);
    for (int i = 0; i < world.width; i++)
        if (coastalland(world, i, 16))
        {
            world.lakesurface_m[at(i, 16, world.width)] = world.elevation_m[at(i, 16, world.width)] + 10.0f;
            break;
        }
    const std::vector<uint8_t> nearsea = landscape::nearseacells(world);
    for (int reach : {0, 2, 3, 4})
    {
        SCOPED_TRACE("reach " + std::to_string(reach));
        RegionalParameters parameters;
        parameters.edges.smoothgate = reach;
        const std::vector<float> inland = landscape::inlandcells(world, parameters);
        const int distance = std::max(1, reach - 1);
        int wrong = 0, ones = 0;
        for (int j = 0; j < world.height; j++)
            for (int i = 0; i < world.width; i++)
            {
                bool sea = false;
                for (int dj = -distance; dj <= distance; dj++)
                    for (int di = -distance; di <= distance; di++)
                        if (j + dj >= 0 && j + dj < world.height)
                            sea = sea || world.network.sea[at((i + di + world.width) % world.width, j + dj, world.width)] != 0;
                const std::size_t global = at(i, j, world.width);
                const float expected = sea && !lakeglobal(world, global) ? 0.0f : 1.0f;
                wrong += inland[global] != expected ? 1 : 0;
                ones += inland[global] != 0.0f ? 1 : 0;
                if (reach <= 2)
                    wrong += inland[global] != (nearsea[global] != 0 ? 0.0f : 1.0f) ? 1 : 0;
            }
        EXPECT_EQ(wrong, 0);
        EXPECT_GT(ones, 0);
    }
}

TEST_F(Landscape, RegionalEdgesSmoothGateExcludesLakeCells)
{
    // A lake cell beside the sea: its block is never gated, while a land cell beside it is.
    RegionalGlobalInputs world = syntheticworld(64);
    const RegionalWindow window = windowc(4);
    int lakei = -1, lakej = -1;
    for (int j = window.row0; j < window.row0 + window.rows && lakei < 0; j++)
        for (int i = window.column0; i < window.column0 + window.columns && lakei < 0; i++)
            if (coastalland(world, i % world.width, j))
            {
                lakei = i % world.width;
                lakej = j;
            }
    ASSERT_GE(lakei, 0);
    world.lakesurface_m[at(lakei, lakej, world.width)] = world.elevation_m[at(lakei, lakej, world.width)] + 10.0f;
    RegionalParameters parameters;
    parameters.edges.smoothgate = 3;
    for (const FineWarp& warp : {FineWarp{}, landscape::regionalwarp(window, world.seed, edgeparameters(0))})
    {
        const std::vector<uint8_t> gate = landscape::seagate(world, window, warp, warp.empty() ? parameters : edgeparameters(0));
        int inlake = 0, cells = 0;
        for (int y = 0; y < window.domainfineheight(); y++)
            for (int x = 0; x < window.domainfinewidth(); x++)
                if (window.globalcolumn(x) == lakei && window.globalrow(y) == lakej)
                {
                    cells++;
                    inlake += gate[at(x, y, window.domainfinewidth())] != 0 ? 1 : 0;
                }
        EXPECT_EQ(cells, window.scale * window.scale);
        EXPECT_EQ(inlake, 0);
    }
}

TEST_F(Landscape, RegionalEdgesBreakerBreaksLongRuns)
{
    // A straight coast along a fine row inside a block row (land north, sea south): every
    // block's run of k cells is longer than straightmax and gets chunks; afterwards no block
    // holds a straight run longer than straightmax, the mask changed, and a lake block is
    // untouched. The result does not depend on what lies beyond the block and its ring.
    RegionalGlobalInputs world = syntheticworld(64);
    const RegionalWindow box = windowc(16);
    const int k = box.scale, finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const int coastrow = 5 * k + 7;  // sea from here southwards, inside the sixth block row
    // A lake on one block of the coast's block row (a land cell of the world, made a lake).
    int lakebx = -1;
    for (int bx = 0; bx < box.domaincolumns && lakebx < 0; bx++)
        if (world.network.sea[at((box.domaincolumn0 + bx) % world.width, box.domainrow0 + 5, world.width)] == 0)
            lakebx = bx;
    ASSERT_GE(lakebx, 0);
    world.lakesurface_m[at((box.domaincolumn0 + lakebx) % world.width, box.domainrow0 + 5, world.width)] = 1000.0f;

    const RegionalParameters parameters = edgeparameters(4);
    std::vector<uint8_t> gate(static_cast<std::size_t>(finewidth) * fineheight, 1), mask(gate.size(), 0);
    for (int y = coastrow; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            mask[at(x, y, finewidth)] = 1;
    ASSERT_EQ(longestblockrun(world, box, mask), k);
    const std::vector<uint8_t> before = mask;
    landscape::breakstraightruns(world, box, parameters, gate, mask);
    EXPECT_NE(mask, before);
    EXPECT_LE(longestblockrun(world, box, mask), parameters.edges.straightmax);
    for (int y = 5 * k; y < 6 * k; y++)
        for (int x = lakebx * k; x < (lakebx + 1) * k; x++)
            EXPECT_EQ(mask[at(x, y, finewidth)], before[at(x, y, finewidth)]) << "lake block changed at " << x << ", " << y;
    // Every change lies within one fine cell row of the coast's side rows (chunks of two rows).
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            if (mask[at(x, y, finewidth)] != before[at(x, y, finewidth)])
                EXPECT_TRUE(y >= coastrow - 2 && y <= coastrow + 1) << x << ", " << y;

    // Locality: a different mask two blocks away from a block changes nothing in it.
    std::vector<uint8_t> other = before;
    for (int y = 0; y < 3 * k; y++)
        for (int x = 0; x < finewidth; x++)
            other[at(x, y, finewidth)] = ((x * 7 + y * 3) % 5 == 0) ? 1 : 0;
    landscape::breakstraightruns(world, box, parameters, gate, other);
    for (int y = 4 * k; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            ASSERT_EQ(other[at(x, y, finewidth)], mask[at(x, y, finewidth)]) << x << ", " << y;
}

TEST_F(Landscape, RegionalEdgesBreakerLeavesShortRunsAndClosedCells)
{
    const RegionalGlobalInputs world = syntheticworld(64);
    const RegionalWindow box = windowc(16);
    const int k = box.scale, finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    const int coastrow = 5 * k + 7;
    const RegionalParameters parameters = edgeparameters(4);
    std::vector<uint8_t> gate(static_cast<std::size_t>(finewidth) * fineheight, 1), mask(gate.size(), 0);

    // Runs of straightmax cells, split by a land cell every fifth column: untouched.
    for (int y = coastrow; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            mask[at(x, y, finewidth)] = (y == coastrow && x % 5 == 4) ? 0 : 1;
    ASSERT_LE(longestblockrun(world, box, mask), parameters.edges.straightmax);
    std::vector<uint8_t> before = mask;
    landscape::breakstraightruns(world, box, parameters, gate, mask);
    EXPECT_EQ(mask, before);

    // A long run with the gate closed on the land side: chunks may raise land, but no cell
    // whose gate is 0 becomes sea.
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
        {
            mask[at(x, y, finewidth)] = y >= coastrow ? 1 : 0;
            gate[at(x, y, finewidth)] = y >= coastrow ? 1 : 0;
        }
    before = mask;
    landscape::breakstraightruns(world, box, parameters, gate, mask);
    EXPECT_NE(mask, before);
    int opened = 0;
    for (std::size_t cell = 0; cell < mask.size(); cell++)
        opened += gate[cell] == 0 && before[cell] == 0 && mask[cell] != 0 ? 1 : 0;
    EXPECT_EQ(opened, 0);
}

TEST_F(Landscape, RegionalOverlapEdgesWideWorld)
{
    // RegionalOverlapCoastWideWorld (regional_surface_tests.cpp) with the smooth gate (reach 3)
    // and the breaker on (straightmax 4): on a 256 x 128 world the land-share box is neither a whole
    // turn nor at a pole; two windows across the continent's western coast, overlapping in
    // 130 .. 140 E, bit-identical at the smallest halo, at k = 8 (runs of up to 8 cells in a
    // block: the breaker acts), k = 4 and k = 2 (the share's reach of 1 and 2 cells).
    const RegionalParameters parameters = edgeparameters(4);
    const RegionalGlobalInputs inputs = syntheticworld(256);
    for (int scale : {8, 4, 2})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        // The windows straddle the equator because the land-share box with the breaker on is
        // wide: a share pass's reach is r = 1 + ceil(max(0, k + hookreach + 2 - k / 2) / k) blocks
        // (regional_surface.cpp, landsharesurface; hookreach 1 for the breaker), 2 at k 8 and 3
        // at k 4 and 2, so the margin is 42 or 50 global cells of 128 rows.
        const RegionalWindow wa = landscape::regionalwindow(256, scale, testhalo, 120.0, -8.0, 140.0, 8.0);
        const RegionalWindow wb = landscape::regionalwindow(256, scale, testhalo, 130.0, -4.0, 150.0, 8.0);
        const int hookreach = parameters.edges.straightmax > 0 ? 1 : 0;
        const int sharereach = 1 + (std::max(0, scale + hookreach + 2 - scale / 2) + scale - 1) / scale;
        const int margin = sharereach * parameters.coast.sharepasses + 2 + 2 * parameters.coast.signedpasses;
        for (const RegionalWindow& window : {wa, wb})
        {
            ASSERT_LT(window.domaincolumns + 2 * margin, 256);
            ASSERT_GT(window.domainrow0 - margin, 0);
            ASSERT_LT(window.domainrow0 + window.domainrows + margin, 128);
        }
        const OverlapCounts counts = compareoverlap(inputs, wa, wb, "edges wide A/B", parameters);
        EXPECT_GT(counts.sea, 0);
        EXPECT_LT(counts.sea, counts.shared);
        if (scale == 8)
        {
            // The breaker acts here: the sea without it differs.
            RegionalParameters without = parameters;
            without.edges.straightmax = 0;
            EXPECT_NE(landscape::runregionalpass(inputs, wa, parameters).sea, landscape::runregionalpass(inputs, wa, without).sea);
        }
    }
}

TEST_F(Landscape, RegionalOverlapEdgesDateLine)
{
    // Windows A and B of the 32 x 16 world (across the date line, domains at a pole: the box
    // is clamped there) with the smooth gate and the breaker, at k = 8 and k = 4.
    const RegionalParameters parameters = edgeparameters(4);
    const RegionalGlobalInputs inputs = syntheticworld(32);
    for (int scale : {8, 4})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts counts = compareoverlap(inputs, windowa(scale), windowb(scale), "edges A/B", parameters);
        EXPECT_GT(counts.sea, 0);
    }
}
