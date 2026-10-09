// Tests of regional_coast.hpp (erosion plan Phase R2: the R6 contract and R7's island sink,
// one-sided coast noise and coast factors), suite Landscape: TEST_F(Landscape,
// RegionalCoast<What>) and the RegionalOverlapIslands* variants with R7's features on, on the
// synthetic world of regional_test_world.hpp (one of them wider than the land-share rule's
// widened box).

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
