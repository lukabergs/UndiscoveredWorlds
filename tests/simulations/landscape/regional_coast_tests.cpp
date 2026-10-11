// Tests of regional_coast.hpp (erosion plan Phase R2: the R6 contract and R7's island sink,
// one-sided coast noise and coast factors), suite Landscape: TEST_F(Landscape,
// RegionalCoast<What>) and the RegionalOverlapIslands* variants with R7's features on, on the
// synthetic world of regional_test_world.hpp (one of them wider than the land-share rule's
// widened box). Multiscale plan S3.15: the shoreline displacement on a straight contour and the
// pocket fill (RegionalCoastShorelineDisplaced, RegionalCoastPocketsFilled; the two-window overlap
// with both on is refinement_levels_tests.cpp's RefinementLevelsOverlapCoastShoreline).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_surface.hpp"

namespace
{
using landscape::FineWarp;
using landscape::StartSurface;

// The R2 test coast (regional_surface_tests.cpp, coastparameters) with R7's features on: the
// one-sided noise, the island sink at half a global cell, the glacial and denudation terms.
RegionalParameters islandparameters()
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
    parameters.coast.onesided = 1;
    parameters.coast.islandmargin = 0.5;
    parameters.coast.glacialamplitude_m = 100.0;
    parameters.coast.glacialdropfraction = 0.9;
    parameters.coast.denudationamplitude_m = 50.0;
    parameters.coast.denudationreference_m = 1000.0;
    return parameters;
}

// Synthetic glaciation (0 within 10 degrees of the equator, 1 from 40 degrees) and
// denudation (0 to 1500 m, varying with longitude), so the coast factors vary along the coast.
RegionalGlobalInputs withcoastfields(RegionalGlobalInputs inputs)
{
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.glaciation.assign(size, 0.0f);
    inputs.denudation_m.assign(size, 0.0f);
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            const double latitude = std::fabs(latitude_rad(j, inputs.height)) * 180.0 / pi;
            inputs.glaciation[at(i, j, inputs.width)] = static_cast<float>(std::clamp((latitude - 10.0) / 30.0, 0.0, 1.0));
            inputs.denudation_m[at(i, j, inputs.width)] =
                static_cast<float>(750.0 * (1.0 + std::sin(3.0 * longitude_rad(i, inputs.width))));
        }
    return inputs;
}

// An all-sea world of 32 x 16 global cells for the sink's synthetic masks.
RegionalGlobalInputs seaworld()
{
    RegionalGlobalInputs inputs;
    inputs.width = 32;
    inputs.height = 16;
    inputs.sealevel = sealevel;
    inputs.network.width = 32;
    inputs.network.height = 16;
    inputs.network.sea.assign(32 * 16, 1);
    return inputs;
}

RegionalWindow sinkbox(int column0, int columns, int row0, int rows, int scale)
{
    RegionalWindow box;
    box.worldwidth = 32;
    box.worldheight = 16;
    box.scale = scale;
    box.domaincolumn0 = column0;
    box.domaincolumns = columns;
    box.domainrow0 = row0;
    box.domainrows = rows;
    return box;
}

// The kept window's fine land cells inside global sea cells, the land-share misses and the
// cells whose start surface disagrees with the sea mask (as checkshares counts them).
struct SeaLand
{
    int landinsea = 0, misses = 0, inconsistent = 0, sealevelsea = 0;
};

SeaLand sealand(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters)
{
    const FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const StartSurface start = landscape::regionalstartsurface(inputs, window, warp, parameters);
    std::vector<float> surface = start.start_m;
    const std::vector<uint8_t> domainsea = landscape::regionalsea(inputs, window, warp, parameters, surface);
    const std::vector<float> kept = keptpart(start.start_m, window);
    SeaLand result;
    result.misses = start.sharemisses;
    const int k = window.scale, width = inputs.width;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
        {
            const uint8_t sea = domainsea[at(window.keptfinecolumnoffset() + x, window.keptfinerowoffset() + y, window.domainfinewidth())];
            const int i = (window.column0 + x / k) % width, j = window.row0 + y / k;
            const std::size_t global = at(i, j, width);
            const float z = kept[at(x, y, window.finewidth())];
            if (sea == 0 && inputs.network.sea[global] != 0)
                result.landinsea++;
            if (sea != 0 && z > inputs.sealevel)
                result.inconsistent++;
            if (sea != 0 && z == inputs.sealevel)
                result.sealevelsea++;
            if (sea == 0 && inputs.network.sea[global] != 0 && z <= inputs.sealevel)
                result.inconsistent++;
        }
    return result;
}
}

TEST_F(Landscape, RegionalCoastFactorsFollowTheFields)
{
    // Inactive at the defaults and without the fields; with glaciation on one cell and the
    // coefficients on, that cell's amplitude and drop factor follow the formula and the
    // others stay 0 and 1.
    RegionalGlobalInputs inputs = syntheticworld(32);
    RegionalParameters parameters;
    EXPECT_FALSE(landscape::coastfactors(inputs, parameters).active);
    parameters.coast.glacialamplitude_m = 100.0;
    parameters.coast.glacialdropfraction = 0.9;
    parameters.coast.denudationamplitude_m = 40.0;
    parameters.coast.denudationreference_m = 1000.0;
    EXPECT_FALSE(landscape::coastfactors(inputs, parameters).active) << "empty fields mean no factors";

    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.glaciation.assign(size, 0.0f);
    inputs.denudation_m.assign(size, 0.0f);
    const std::size_t glacial = at(5, 3, inputs.width), eroded = at(9, 7, inputs.width);
    inputs.glaciation[glacial] = 0.5f;
    inputs.denudation_m[eroded] = 2500.0f;  // beyond the reference: capped at the full amplitude
    const landscape::CoastFactors factors = landscape::coastfactors(inputs, parameters);
    ASSERT_TRUE(factors.active);
    ASSERT_EQ(factors.amplitude_m.size(), size);
    ASSERT_EQ(factors.dropfactor.size(), size);
    EXPECT_FLOAT_EQ(factors.amplitude_m[glacial], 50.0f);
    EXPECT_FLOAT_EQ(factors.dropfactor[glacial], 0.55f);
    EXPECT_FLOAT_EQ(factors.amplitude_m[eroded], 40.0f);
    EXPECT_FLOAT_EQ(factors.dropfactor[eroded], 1.0f);
    int others = 0;
    for (std::size_t cell = 0; cell < size; cell++)
        if (cell != glacial && cell != eroded && (factors.amplitude_m[cell] != 0.0f || factors.dropfactor[cell] != 1.0f))
            others++;
    EXPECT_EQ(others, 0);
}

TEST_F(Landscape, RegionalCoastIslandSinkIsOffByDefault)
{
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const RegionalParameters parameters;
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    std::vector<uint8_t> gate(count, 1), mask(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        mask[cell] = (cell % 5 == 0) ? 1 : 0;
    const std::vector<uint8_t> before = mask;
    landscape::sinksmallislands(inputs, window, parameters, gate, mask);
    EXPECT_EQ(mask, before);
}

TEST_F(Landscape, RegionalCoastIslandSink)
{
    // A box of 6 x 6 global cells (columns 10 .. 15, rows 4 .. 9) at k = 8, margin 0.5 (a ring
    // at 4 fine cells), an all-sea world but for global cell (13, 7); one island per case:
    // a block whose ring is clear sinks (but for a cell whose gate is 0); one whose ring holds
    // land stays, as do one beside a world one-cell island, one beside a mountain island, one
    // whose ring reaches beyond the box and one in a global land cell.
    RegionalGlobalInputs inputs = seaworld();
    inputs.network.sea[at(13, 7, 32)] = 0;
    inputs.island.assign(32 * 16, 0);
    inputs.mountainisland.assign(32 * 16, 0);
    inputs.island[at(11, 9, 32)] = 1;          // beside global cell (12, 8)
    inputs.mountainisland[at(15, 9, 32)] = 1;  // beside global cell (14, 8)
    const RegionalWindow box = sinkbox(10, 6, 4, 6, 8);
    RegionalParameters parameters;
    parameters.coast.islandmargin = 0.5;
    const int finewidth = box.domainfinewidth();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(box.domainfineheight());
    std::vector<uint8_t> gate(count, 1), start(count, 1);
    auto land = [&](int x, int y) { start[at(x, y, finewidth)] = 0; };
    // Global (12, 5), block x 16 .. 23, y 8 .. 15: sinks; (21, 13) is gated off.
    land(19, 11), land(20, 11), land(19, 12), land(21, 13);
    gate[at(21, 13, finewidth)] = 0;
    // Global (14, 5), block x 32 .. 39: land on its ring at x = 43 (in a block whose own ring
    // leaves the box, so it stays too).
    land(35, 11), land(43, 11);
    land(19, 35);  // global (12, 8), beside the world island
    land(36, 36);  // global (14, 8), beside the mountain island
    land(3, 19);   // global (10, 6): the ring's west side lies beyond the box
    land(27, 27);  // global (13, 7), a land cell
    std::vector<uint8_t> expected = start;
    for (const auto& [x, y] : std::initializer_list<std::pair<int, int>>{{19, 11}, {20, 11}, {19, 12}})
        expected[at(x, y, finewidth)] = 1;
    std::vector<uint8_t> mask = start;
    landscape::sinksmallislands(inputs, box, parameters, gate, mask);
    EXPECT_EQ(mask, expected);
    for (const auto& [x, y] : std::initializer_list<std::pair<int, int>>{{21, 13}, {35, 11}, {43, 11}, {19, 35}, {36, 36}, {3, 19}, {27, 27}})
        EXPECT_EQ(mask[at(x, y, finewidth)], 0) << "kept land at fine (" << x << ", " << y << ")";

    // The same mask on the world without the island flags: the blocks beside them sink too.
    inputs.island.clear();
    inputs.mountainisland.clear();
    mask = start;
    landscape::sinksmallislands(inputs, box, parameters, gate, mask);
    EXPECT_EQ(mask[at(19, 35, finewidth)], 1);
    EXPECT_EQ(mask[at(36, 36, finewidth)], 1);
    EXPECT_EQ(mask[at(3, 19, finewidth)], 0);
    EXPECT_EQ(mask[at(35, 11, finewidth)], 0);
}

TEST_F(Landscape, RegionalCoastIslandSinkWholeTurn)
{
    // Around a whole turn the ring wraps: the island of block 0 (column 0) stays while the
    // ring's wrapped west side (the last fine columns) holds land, and sinks without it.
    const RegionalGlobalInputs inputs = seaworld();
    const RegionalWindow box = sinkbox(0, 32, 4, 6, 8);
    RegionalParameters parameters;
    parameters.coast.islandmargin = 0.5;
    const int finewidth = box.domainfinewidth();
    const std::size_t count = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(box.domainfineheight());
    const std::vector<uint8_t> gate(count, 1);
    std::vector<uint8_t> mask(count, 1);
    mask[at(3, 19, finewidth)] = 0;
    mask[at(finewidth - 4, 19, finewidth)] = 0;  // on block 0's ring, 4 cells west of it
    std::vector<uint8_t> wrapped = mask;
    landscape::sinksmallislands(inputs, box, parameters, gate, wrapped);
    EXPECT_EQ(wrapped[at(3, 19, finewidth)], 0);
    mask[at(finewidth - 4, 19, finewidth)] = 1;
    landscape::sinksmallislands(inputs, box, parameters, gate, mask);
    EXPECT_EQ(mask[at(3, 19, finewidth)], 1);
}

TEST_F(Landscape, RegionalCoastOneSidedAndSink)
{
    // On the 64 x 32 world (window C at k = 8, the R2 test coast): the one-sided noise and the
    // island sink each leave fewer fine land cells in global sea cells than the two-sided R2
    // coast (the noise may no longer lift a sea cell's decision above the sea level; the sink
    // removes blocks of isolated land); the land share holds, the surface agrees with the mask.
    const RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = windowc(8);
    RegionalParameters twosided = islandparameters();
    twosided.coast.onesided = 0;
    twosided.coast.islandmargin = 0.0;
    twosided.coast.glacialamplitude_m = twosided.coast.glacialdropfraction = twosided.coast.denudationamplitude_m = 0.0;
    RegionalParameters onesided = twosided;
    onesided.coast.onesided = 1;
    RegionalParameters sink = onesided;
    sink.coast.islandmargin = 0.5;
    const SeaLand a = sealand(inputs, window, twosided), b = sealand(inputs, window, onesided), c = sealand(inputs, window, sink);
    std::printf("[regional] fine land in global sea cells, k 8: two-sided %d, one-sided %d, with the sink %d; "
                "sea cells at the sea level %d, %d, %d\n",
                a.landinsea, b.landinsea, c.landinsea, a.sealevelsea, b.sealevelsea, c.sealevelsea);
    EXPECT_GT(a.landinsea, 0);
    EXPECT_LT(b.landinsea, a.landinsea);
    EXPECT_LE(c.landinsea, b.landinsea);
    for (const SeaLand& result : {a, b, c})
    {
        EXPECT_EQ(result.misses, 0);
        EXPECT_EQ(result.inconsistent, 0);
    }
    EXPECT_LE(c.sealevelsea, b.sealevelsea) << "sunk land takes the shore floor, not the sea level";
}

TEST_F(Landscape, RegionalOverlapIslands)
{
    // The pass with R7's features on (the one-sided noise, the sink, the glacial and
    // denudation terms on synthetic fields) is bit-identical on the overlaps of windows A/B
    // (across the date line, at a pole) and C/D (clear of the poles) at the smallest halo.
    const RegionalParameters parameters = islandparameters();
    const RegionalGlobalInputs columns = withcoastfields(syntheticworld(32));
    const RegionalGlobalInputs rows = withcoastfields(syntheticworld(64));
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts first = compareoverlap(columns, windowa(scale), windowb(scale), "islands A/B", parameters);
        EXPECT_EQ(first.shared, 4 * 5 * scale * scale);
        const OverlapCounts second = compareoverlap(rows, windowc(scale), windowd(scale), "islands C/D", parameters);
        EXPECT_EQ(second.shared, 12 * 4 * scale * scale);
        EXPECT_GT(second.sea, 0);
    }
}

TEST_F(Landscape, RegionalOverlapIslandsWideWorld)
{
    // On a 256 x 128 world, with R7's features on, the land-share box is neither a whole turn
    // nor at a pole (asserted with a bound on the margin that covers a share pass's reach with
    // the sink's ring: 3 global cells for k = 4, 4 for k = 2), so its edge is exercised: the
    // windows of RegionalOverlapCoastWideWorld, bit-identical at the smallest halo. Six share
    // and six sign-keeping passes keep the widened box clear of the poles.
    RegionalParameters parameters = islandparameters();
    parameters.coast.sharepasses = 6;
    parameters.coast.signedpasses = 6;
    const RegionalGlobalInputs inputs = withcoastfields(syntheticworld(256));
    for (int scale : {4, 2})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow wa = landscape::regionalwindow(256, scale, testhalo, 120.0, -5.0, 140.0, 15.0);
        const RegionalWindow wb = landscape::regionalwindow(256, scale, testhalo, 130.0, 0.0, 150.0, 20.0);
        const int margin = (scale <= 2 ? 4 : 3) * parameters.coast.sharepasses + 4 + 2 * parameters.coast.signedpasses;
        for (const RegionalWindow& window : {wa, wb})
        {
            ASSERT_LT(window.domaincolumns + 2 * margin, 256);
            ASSERT_GT(window.domainrow0 - margin, 0);
            ASSERT_LT(window.domainrow0 + window.domainrows + margin, 128);
        }
        const auto start = std::chrono::steady_clock::now();
        const OverlapCounts counts = compareoverlap(inputs, wa, wb, "islands wide A/B", parameters);
        std::printf("[regional] islands wide world, k %d: both windows in %.2f s\n", scale,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        EXPECT_GT(counts.sea, 0);
        EXPECT_LT(counts.sea, counts.shared);
    }
}

namespace
{
// Multiscale plan S3.15: a world of 32 x 16 global cells whose coast is a straight meridian, the
// west half sea at -65 m and the east half land at +14 m (window A's pair at 58 E: global cells
// (168, 58) and (169, 58) of erosion run 9), no lake.
RegionalGlobalInputs straightcoast()
{
    RegionalGlobalInputs inputs;
    inputs.width = 32;
    inputs.height = 16;
    inputs.radius_m = radius_m;
    inputs.sealevel = sealevel;
    inputs.seed = 1;
    const std::size_t size = 32 * 16;
    inputs.elevation_m.resize(size);
    inputs.lakesurface_m.assign(size, 0.0f);
    inputs.network.width = 32;
    inputs.network.height = 16;
    inputs.network.sea.resize(size);
    for (int j = 0; j < 16; j++)
        for (int i = 0; i < 32; i++)
        {
            inputs.elevation_m[at(i, j, 32)] = i < 16 ? -65.0f : 14.0f;
            inputs.network.sea[at(i, j, 32)] = i < 16 ? 1 : 0;
        }
    return inputs;
}

// The land-share coast without noise or the R2 warp, so the decision's contour is the straight
// interpolated one; `warpcells` the shoreline displacement.
RegionalParameters straightparameters(double warpcells)
{
    RegionalParameters parameters;
    parameters.noiseamplitude_m = 0.0;
    parameters.noiserelieffraction = 0.0;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.coast.onesided = 1;
    parameters.coast.warpcells = warpcells;
    parameters.coast.warppersistence = 0.7;
    return parameters;
}

// The kept window's sea mask (regionalsea on regionalstartsurface's surface), the share misses.
struct KeptSea
{
    std::vector<uint8_t> sea;
    int misses = 0;
};

KeptSea keptsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters)
{
    const FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const StartSurface start = landscape::regionalstartsurface(inputs, window, warp, parameters);
    std::vector<float> surface = start.start_m;
    const std::vector<uint8_t> domainsea = landscape::regionalsea(inputs, window, warp, parameters, surface);
    KeptSea result;
    result.misses = start.sharemisses;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
            result.sea.push_back(
                domainsea[at(window.keptfinecolumnoffset() + x, window.keptfinerowoffset() + y, window.domainfinewidth())]);
    return result;
}

// The longest run of coast cells (land with a 4-neighbour sea cell) down a column of the kept window.
int longestverticalcoast(const std::vector<uint8_t>& sea, int width, int height)
{
    int longest = 0;
    for (int x = 0; x < width; x++)
    {
        int run = 0;
        for (int y = 0; y < height; y++)
        {
            bool coast = false;
            if (sea[at(x, y, width)] == 0)
            {
                const int nx[4] = {x - 1, x + 1, x, x};
                const int ny[4] = {y, y, y - 1, y + 1};
                for (int n = 0; n < 4; n++)
                    coast = coast || (nx[n] >= 0 && nx[n] < width && ny[n] >= 0 && ny[n] < height && sea[at(nx[n], ny[n], width)] != 0);
            }
            run = coast ? run + 1 : 0;
            longest = std::max(longest, run);
        }
    }
    return longest;
}
}

TEST_F(Landscape, RegionalCoastShorelineDisplaced)
{
    // Multiscale plan S3.15 (a): on a straight interpolated contour (the coast between global columns 15
    // and 16, window columns 14 .. 17, rows 5 .. 10 at k = 16), the displacement breaks the straight
    // coast (its longest vertical run under one global cell, against the window's whole height without
    // it) while every global cell keeps its land share and the window's land share stays within 5 %.
    const RegionalGlobalInputs inputs = straightcoast();
    const RegionalWindow window = landscape::regionalwindow(32, 16, testhalo, -22.0, -30.0, 22.0, 30.0);
    ASSERT_EQ(window.column0, 14);
    ASSERT_EQ(window.columns, 4);
    const KeptSea straight = keptsea(inputs, window, straightparameters(0.0));
    const KeptSea displaced = keptsea(inputs, window, straightparameters(0.25));
    const int width = window.finewidth(), height = window.fineheight();
    const int before = longestverticalcoast(straight.sea, width, height);
    const int after = longestverticalcoast(displaced.sea, width, height);
    const auto landshare = [&](const std::vector<uint8_t>& sea) {
        return static_cast<double>(std::count(sea.begin(), sea.end(), uint8_t{0})) / static_cast<double>(sea.size());
    };
    std::printf("[regional] straight coast, k 16: longest vertical coast run %d cells without the displacement, %d with "
                "it; land share %.4f and %.4f\n",
                before, after, landshare(straight.sea), landshare(displaced.sea));
    EXPECT_EQ(before, height) << "the undisplaced contour is a straight meridian";
    EXPECT_LT(after, window.scale);
    EXPECT_NE(straight.sea, displaced.sea);
    EXPECT_EQ(straight.misses, 0);
    EXPECT_EQ(displaced.misses, 0);
    EXPECT_NEAR(landshare(displaced.sea), landshare(straight.sea), 0.05);
    // Every block keeps its share (land at least half in a land cell, at most half in a sea cell).
    for (int r = 0; r < window.rows; r++)
        for (int c = 0; c < window.columns; c++)
        {
            int land = 0;
            for (int y = r * 16; y < (r + 1) * 16; y++)
                for (int x = c * 16; x < (c + 1) * 16; x++)
                    land += displaced.sea[at(x, y, width)] == 0 ? 1 : 0;
            if (window.column0 + c < 16)
                EXPECT_LE(land, 128) << "global sea cell (" << window.column0 + c << ", " << window.row0 + r << ")";
            else
                EXPECT_GE(land, 128) << "global land cell (" << window.column0 + c << ", " << window.row0 + r << ")";
        }

    // The displacement's bound and band: |w| < bound() everywhere; a level whose band holds no octave
    // (k_l = parentscale) and warpcells 0 are inactive; the band of k 64 under k 16 is finer than k 16's.
    landscape::CoastParameters coast;
    coast.warpcells = 0.25;
    coast.warppersistence = 0.7;
    const landscape::ShorelineWarp first(32, 16, 16, 1, 1, coast), deeper(32, 16, 64, 16, 1, coast), empty(32, 16, 16, 16, 1, coast);
    ASSERT_TRUE(first.active());
    ASSERT_TRUE(deeper.active());
    EXPECT_FALSE(empty.active());
    EXPECT_LT(deeper.bound(), first.bound());
    double largest = 0.0;
    for (std::int64_t row = 0; row < 16 * 16; row += 3)
        for (std::int64_t column = 0; column < 32 * 16; column += 5)
        {
            double dx = 0.0, dy = 0.0;
            first.at(column, row, dx, dy);
            largest = std::max(largest, std::hypot(dx, dy));
        }
    EXPECT_LT(largest, first.bound());
    EXPECT_GT(largest, 0.25 * first.bound());
    coast.warpcells = 0.0;
    EXPECT_FALSE(landscape::ShorelineWarp(32, 16, 16, 1, 1, coast).active());
}

TEST_F(Landscape, RegionalCoastPocketsFilled)
{
    // Multiscale plan S3.15 (b): global land cells (10 .. 14, 4 .. 8) in an all-sea world, a box of
    // global columns 8 .. 16 and rows 2 .. 10 at k = 4, the global sea cells' blocks all sea. A pocket in
    // (12, 6), three global cells from any sea cell, and one in (11, 4) beside the sea cell (11, 3) but
    // cut from it by land are filled; a bay in (10, 6) open to the sea of (9, 6) and a sea cell of
    // (11, 4) touching (11, 3)'s sea are kept. Reach 0 fills nothing.
    RegionalGlobalInputs inputs = seaworld();
    for (int j = 4; j <= 8; j++)
        for (int i = 10; i <= 14; i++)
            inputs.network.sea[at(i, j, 32)] = 0;
    const RegionalWindow box = sinkbox(8, 9, 2, 9, 4);
    const int finewidth = box.domainfinewidth(), fineheight = box.domainfineheight();
    std::vector<uint8_t> sea(static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fineheight), 0);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            sea[at(x, y, finewidth)] = inputs.network.sea[at(8 + x / 4, 2 + y / 4, 32)];
    const std::initializer_list<std::pair<int, int>> pockets = {{17, 17}, {18, 17}, {17, 18}, {13, 9}, {14, 9}};
    const std::initializer_list<std::pair<int, int>> kept = {{8, 17}, {9, 17}, {12, 8}};
    for (const auto& [x, y] : pockets)
        sea[at(x, y, finewidth)] = 1;
    for (const auto& [x, y] : kept)
        sea[at(x, y, finewidth)] = 1;
    const std::vector<uint8_t> filled = landscape::seapockets(inputs, box, 2, sea);
    std::vector<uint8_t> expected(sea.size(), 0);
    for (const auto& [x, y] : pockets)
        expected[at(x, y, finewidth)] = 1;
    EXPECT_EQ(filled, expected);
    const landscape::PocketCount count = landscape::countpockets(filled, finewidth, 0, 0, finewidth, fineheight);
    EXPECT_EQ(count.pockets, 2);
    EXPECT_EQ(count.pocketcells, 5);
    // At reach 1 the bay is still open (its sea cell is next door); at reach 0 nothing is filled.
    EXPECT_EQ(landscape::seapockets(inputs, box, 1, sea), expected);
    EXPECT_EQ(landscape::seapockets(inputs, box, 0, sea), std::vector<uint8_t>(sea.size(), 0));
    // With a core that leaves out (9, 6)'s sea the bay is a pocket too.
    std::vector<uint8_t> core(sea.size(), 0);
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
            core[at(x, y, finewidth)] = inputs.network.sea[at(8 + x / 4, 2 + y / 4, 32)] != 0 && x / 4 != 1 ? 1 : 0;
    const std::vector<uint8_t> strict = landscape::seapockets(inputs, box, 1, sea, core);
    EXPECT_EQ(strict[at(8, 17, finewidth)], 1);
    EXPECT_EQ(strict[at(17, 17, finewidth)], 1);
    EXPECT_EQ(strict[at(12, 8, finewidth)], 0);

    // The basin of the pocket in (12, 6): the land below the level 4-connected to it within the reach
    // is raised with it; low land cut off by higher land is not.
    std::vector<uint8_t> after = sea;
    for (std::size_t cell = 0; cell < sea.size(); cell++)
        if (filled[cell] != 0)
            after[cell] = 0;
    std::vector<float> lowland(sea.size(), 50.0f);
    lowland[at(16, 17, finewidth)] = lowland[at(15, 17, finewidth)] = 5.0f;  // beside the pocket
    lowland[at(26, 26, finewidth)] = 5.0f;                                    // cut off
    const std::vector<uint8_t> basin = landscape::pocketbasins(inputs, box, 2, after, filled, lowland, 20.0f);
    for (const auto& [x, y] : pockets)
        EXPECT_EQ(basin[at(x, y, finewidth)], 1);
    EXPECT_EQ(basin[at(16, 17, finewidth)], 1);
    EXPECT_EQ(basin[at(15, 17, finewidth)], 1);
    EXPECT_EQ(basin[at(14, 17, finewidth)], 0);
    EXPECT_EQ(basin[at(26, 26, finewidth)], 0);
    EXPECT_EQ(std::count(basin.begin(), basin.end(), uint8_t{1}), 7);

    // A global lake cell at (11, 6), surface 10 m, beside the pocket's (12, 6): the lake's bed is never a
    // basin, neither in the lake cell (15, 17) nor below its surface in the neighbour (16, 17); the pocket
    // itself is still filled and raised; the pocket search does not cross the lake cell.
    RegionalGlobalInputs laked = inputs;
    laked.lakesurface_m.assign(32 * 16, 0.0f);
    laked.lakesurface_m[at(11, 6, 32)] = 10.0f;
    EXPECT_EQ(landscape::seapockets(laked, box, 2, sea), expected);
    const std::vector<uint8_t> lakebasin = landscape::pocketbasins(laked, box, 2, after, filled, lowland, 20.0f);
    for (const auto& [x, y] : pockets)
        EXPECT_EQ(lakebasin[at(x, y, finewidth)], 1);
    EXPECT_EQ(lakebasin[at(15, 17, finewidth)], 0);
    EXPECT_EQ(lakebasin[at(16, 17, finewidth)], 0);
    EXPECT_EQ(std::count(lakebasin.begin(), lakebasin.end(), uint8_t{1}), 5);

    // Through the land-share rule (the R2 test coast with R7's features on the synthetic world, window
    // C at k = 8): with the fill on, no kept fine sea cell of a global land cell is cut from the sea of
    // the global sea cells within two global cells, the shares hold and the surface agrees with the mask.
    const RegionalGlobalInputs world = syntheticworld(64);
    const RegionalWindow window = windowc(8);
    RegionalParameters parameters = islandparameters();
    parameters.coast.pocketreach = 2;
    const FineWarp warp = landscape::regionalwarp(window, world.seed, parameters);
    const StartSurface start = landscape::regionalstartsurface(world, window, warp, parameters);
    std::vector<float> surface = start.start_m;
    const std::vector<uint8_t> domainsea = landscape::regionalsea(world, window, warp, parameters, surface);
    const std::vector<uint8_t> left = landscape::seapockets(world, window, 2, domainsea);
    int keptleft = 0;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
            keptleft += left[at(window.keptfinecolumnoffset() + x, window.keptfinerowoffset() + y, window.domainfinewidth())];
    std::printf("[regional] pockets filled on window C, k 8: %d (%d cells), left %d (%d cells); kept cells still pockets "
                "against the global sea cells' sea: %d\n",
                start.pocketsfilled, start.pocketcellsfilled, start.pockets, start.pocketcells, keptleft);
    EXPECT_EQ(start.pockets, 0);
    EXPECT_EQ(keptleft, 0);
    EXPECT_EQ(start.sharemisses, 0);
    const SeaLand check = sealand(world, window, parameters);
    EXPECT_EQ(check.inconsistent, 0);
}
