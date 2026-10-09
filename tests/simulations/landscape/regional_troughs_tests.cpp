// Tests of landscape/regional_troughs.hpp (erosion plan Phase E5: the E5.0 contract, which
// step E5.1 implements: glacial troughs, the U-shaped valleys), suite Landscape:
// TEST_F(Landscape, RegionalTroughs<What>) on hand-built fine channels and the synthetic
// world of regional_test_world.hpp with withglacialfields. E5.1 adds the cross-section
// tests and the overlap variants with the troughs on (RegionalOverlapTroughs*), on a
// world wider than the widened box.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_troughs.hpp"

#include <limits>
#include <utility>

TEST_F(Landscape, RegionalTroughsAreOffByDefault)
{
    // At the defaults (width_km 0) the hook changes nothing, with or without the fields.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    RegionalParameters parameters;
    parameters.cirques.radius_km = 2.0;  // another feature on: the fields are active, the troughs still off
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    landscape::FineChannels channels;
    channels.channel.assign(count, 1);
    channels.discharge_m3s.assign(count, 5000.0f);
    channels.globalcell.assign(count, 0);
    const std::vector<uint8_t> sea(count, 0);
    std::vector<float> elevation(count, 1000.0f);
    landscape::carvetroughs(inputs, window, parameters, glacial, channels, sea, elevation);
    for (std::size_t cell = 0; cell < count; cell++)
        ASSERT_EQ(elevation[cell], 1000.0f);
}

namespace
{
// A hand-built trough bench on the domain of window A at k 8 (the half-width may reach 8
// fine cells): land at `land_m` everywhere, fully glaciated fine fields, no channel and no
// sea until a test places them; `global` is the land global cell with the largest ice
// catchment, `area_km2` that catchment, `cell_km` the fine cell's size on the ground.
struct TroughBench
{
    RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    RegionalWindow window = windowa(8);
    int width = window.domainfinewidth(), height = window.domainfineheight();
    std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    landscape::FineGlacial glacial;
    landscape::FineChannels channels;
    std::vector<uint8_t> sea;
    std::vector<float> elevation;
    int32_t global = -1;
    double area_km2 = 0.0;
    double cell_km = 0.0;
    RegionalParameters parameters;

    explicit TroughBench(float land_m = 2000.0f)
    {
        glacial.active = true;
        glacial.glaciation.assign(count, 1.0f);
        channels.channel.assign(count, 0);
        channels.globalcell.assign(count, -1);
        channels.discharge_m3s.assign(count, 0.0f);
        sea.assign(count, 0);
        elevation.assign(count, land_m);
        const std::vector<double> catchment = landscape::glacialcatchment_km2(inputs);
        for (std::size_t cell = 0; cell < catchment.size(); cell++)
            if (catchment[cell] > area_km2)
            {
                area_km2 = catchment[cell];
                global = static_cast<int32_t>(cell);
            }
        cell_km = landscape::finepatchmetric(window, inputs.radius_m).dy_m / 1000.0;
        // The reference at the axis's own catchment (ratio 1): the half-width is width_km.
        parameters.troughs.reference_km2 = area_km2;
        parameters.troughs.widthmax_km = 1.0e6;
    }

    std::size_t cell(int x, int y) const { return at(x, y, width); }
    void axis(int x, int y, float bed_m)
    {
        channels.channel[cell(x, y)] = 1;
        channels.globalcell[cell(x, y)] = global;
        elevation[cell(x, y)] = bed_m;
    }
    void carve() { landscape::carvetroughs(inputs, window, parameters, glacial, channels, sea, elevation); }
    // How many fine cells east of (x, y) were lowered, in a row.
    int reach(int x, int y, float land_m) const
    {
        int cells = 0;
        while (x + cells + 1 < width && elevation[cell(x + cells + 1, y)] < land_m)
            cells++;
        return cells;
    }
};

// The number of kept fine cells whose elevation differs between two passes.
int troughchangedcells(const RegionalFields& a, const RegionalFields& b)
{
    int changed = 0;
    for (std::size_t cell = 0; cell < a.elevation_m.size(); cell++)
        changed += a.elevation_m[cell] != b.elevation_m[cell] ? 1 : 0;
    return changed;
}

// The changed kept cells of `a`'s window that lie in the global cells `other`'s kept
// box shares with it (the overlap compareoverlap compares), so the overlap test shows
// the troughs act on the shared cells, not only somewhere in the window (E5.1's review).
int troughchangedonoverlap(const RegionalFields& a, const RegionalFields& b, const RegionalWindow& other)
{
    const RegionalWindow& window = a.window;
    const int k = window.scale, width = window.finewidth();
    int changed = 0;
    for (std::size_t cell = 0; cell < a.elevation_m.size(); cell++)
    {
        if (a.elevation_m[cell] == b.elevation_m[cell])
            continue;
        const int x = static_cast<int>(cell % static_cast<std::size_t>(width));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(width));
        const int i = (window.column0 + x / k) % window.worldwidth, j = window.row0 + y / k;
        const int column = ((i - other.column0) % window.worldwidth + window.worldwidth) % window.worldwidth;
        if (column < other.columns && j >= other.row0 && j < other.row0 + other.rows)
            changed++;
    }
    return changed;
}

// The troughs on the synthetic world, where a fine cell is 313 km at k 4 (32 columns) and
// a global land cell's ice catchment is of the order of 1e6 km2: a half-width of 600 km x
// (A_i / 1e6 km2)^0.5, so the troughs span one to k fine cells; walls of 1,000 m.
RegionalParameters troughparameters()
{
    RegionalParameters parameters;
    parameters.troughs.width_km = 600.0;
    parameters.troughs.reference_km2 = 1.0e6;
    parameters.troughs.widthmax_km = 2500.0;
    parameters.troughs.wall_m = 1000.0;
    return parameters;
}
}

TEST_F(Landscape, RegionalTroughsCrossSection)
{
    // An east-west channel row at bed 1,000 m in land of 2,000 m, half-width 4.5 fine cells,
    // walls of 810 m, exponent 2: the target at the perpendicular distance d is
    // 1,000 + 810 x (d / 4.5)^2 = 1,000 + 40 d^2: 1,040, 1,160, 1,360, 1,640 at d = 1 to 4;
    // d = 5 lies beyond the half-width (land kept). The nearest axis is the one straight
    // across (any other lies farther), so the row's middle shows the parabola on both sides.
    TroughBench bench;
    bench.parameters.troughs.width_km = 4.5 * bench.cell_km;
    bench.parameters.troughs.wall_m = 810.0;
    const int y0 = bench.height / 2;
    for (int x = 0; x < bench.width; x++)
        bench.axis(x, y0, 1000.0f);
    bench.carve();
    const int x0 = bench.width / 2;
    const double expected[] = {1000.0, 1040.0, 1160.0, 1360.0, 1640.0, 2000.0, 2000.0};
    for (int d = 0; d <= 6; d++)
    {
        SCOPED_TRACE("d = " + std::to_string(d));
        EXPECT_NEAR(bench.elevation[bench.cell(x0, y0 + d)], expected[d], 1.0e-3);
        EXPECT_NEAR(bench.elevation[bench.cell(x0, y0 - d)], expected[d], 1.0e-3);
    }
    // The exponent 4: a flatter floor and steeper walls, 1,000 + 810 x (d / 4.5)^4 =
    // 1,000 + 160 x (d / 3)^4: 1,001.98, 1,031.60, 1,160.00, 1,505.68 at d = 1 to 4.
    TroughBench flatter;
    flatter.parameters.troughs.width_km = 4.5 * flatter.cell_km;
    flatter.parameters.troughs.wall_m = 810.0;
    flatter.parameters.troughs.exponent = 4.0;
    for (int x = 0; x < flatter.width; x++)
        flatter.axis(x, y0, 1000.0f);
    flatter.carve();
    const double quartic[] = {1000.0, 1000.0 + 810.0 / 410.0625, 1000.0 + 810.0 * 16.0 / 410.0625,
                              1000.0 + 810.0 * 81.0 / 410.0625, 1000.0 + 810.0 * 256.0 / 410.0625};
    for (int d = 0; d <= 4; d++)
        EXPECT_NEAR(flatter.elevation[flatter.cell(x0, y0 + d)], quartic[d], 1.0e-3) << "d = " << d;
    EXPECT_NEAR(quartic[3], 1160.0, 1.0e-9);
    EXPECT_NEAR(quartic[4], 1505.68, 0.01);
}

TEST_F(Landscape, RegionalTroughsHalfWidth)
{
    // One axis cell; walls of 999 m above a bed of 1,000 m in land of 2,000 m, so every cell
    // within the half-width W (Euclidean, fine cells) is lowered and none beyond: the reach
    // along the row is floor(W).
    const float land = 2000.0f;
    const auto reachfor = [&](double width_cells, double reference_ratio, double widthmax_cells) {
        TroughBench bench(land);
        bench.parameters.troughs.width_km = width_cells * bench.cell_km;
        bench.parameters.troughs.reference_km2 = bench.area_km2 * reference_ratio;
        bench.parameters.troughs.widthmax_km = widthmax_cells * bench.cell_km;
        bench.parameters.troughs.wall_m = 999.0;
        const int x0 = bench.width / 2, y0 = bench.height / 2;
        bench.axis(x0, y0, 1000.0f);
        bench.carve();
        // The diagonal: (2, 2) lies at 2.83 cells, inside W = 2.6 x 2 = 5.2 and outside 2.6.
        // (Half-widths off whole cells, so the km round trip cannot move a cell across.)
        const bool diagonal = bench.elevation[bench.cell(x0 + 2, y0 + 2)] < land;
        return std::make_pair(bench.reach(x0, y0, land), diagonal);
    };
    // A_i at the reference: W = width_km = 2.6 cells.
    EXPECT_EQ(reachfor(2.6, 1.0, 100.0), std::make_pair(2, false));
    // A_i four times the reference: W = 2.6 x 4^0.5 = 5.2 cells.
    EXPECT_EQ(reachfor(2.6, 0.25, 100.0), std::make_pair(5, true));
    // A_i 16 times: W = 10.4 cells, capped at k = 8.
    EXPECT_EQ(reachfor(2.6, 1.0 / 16.0, 100.0).first, 8);
    // The cap widthmax_km at 3.5 cells: W = 3.5.
    EXPECT_EQ(reachfor(2.6, 0.25, 3.5).first, 3);
    // Under one fine cell: no trough.
    EXPECT_EQ(reachfor(0.9, 1.0, 100.0).first, 0);

    // An axis whose global cell has no ice catchment (a sea cell's 0) or no global cell
    // (-1) cuts nothing.
    TroughBench bench(land);
    bench.parameters.troughs.width_km = 3.0 * bench.cell_km;
    const std::vector<double> catchment = landscape::glacialcatchment_km2(bench.inputs);
    int32_t empty = -1;
    for (std::size_t cell = 0; cell < catchment.size() && empty < 0; cell++)
        if (catchment[cell] == 0.0)
            empty = static_cast<int32_t>(cell);
    ASSERT_GE(empty, 0);
    const int x0 = bench.width / 2, y0 = bench.height / 2;
    bench.axis(x0, y0, 1000.0f);
    bench.channels.globalcell[bench.cell(x0, y0)] = empty;
    bench.axis(x0, y0 + 20, 1000.0f);
    bench.channels.globalcell[bench.cell(x0, y0 + 20)] = -1;
    bench.carve();
    EXPECT_EQ(bench.reach(x0, y0, land), 0);
    EXPECT_EQ(bench.reach(x0, y0 + 20, land), 0);
}

TEST_F(Landscape, RegionalTroughsGlaciationThreshold)
{
    // The axis must be glaciated at troughs.glaciationmin (0.5) or more on the fine field;
    // inactive fields cut nothing.
    const float land = 2000.0f;
    for (const float glaciation : {0.49f, 0.5f, 1.0f})
    {
        TroughBench bench(land);
        bench.parameters.troughs.width_km = 3.5 * bench.cell_km;
        bench.glacial.glaciation.assign(bench.count, glaciation);
        const int x0 = bench.width / 2, y0 = bench.height / 2;
        bench.axis(x0, y0, 1000.0f);
        bench.carve();
        EXPECT_EQ(bench.reach(x0, y0, land), glaciation >= 0.5f ? 3 : 0) << "glaciation " << glaciation;
    }
    TroughBench inactive(land);
    inactive.parameters.troughs.width_km = 3.5 * inactive.cell_km;
    inactive.glacial.active = false;
    const int x0 = inactive.width / 2, y0 = inactive.height / 2;
    inactive.axis(x0, y0, 1000.0f);
    inactive.carve();
    EXPECT_EQ(inactive.reach(x0, y0, land), 0);
}

TEST_F(Landscape, RegionalTroughsKeepSeaChannelsAndTheSeaLevel)
{
    // Half-width 4 cells: a sea cell and another channel cell (a tributary at 1,800 m)
    // inside it are kept, an axis on a sea cell cuts nothing, the targets below the sea
    // level stop at the next float above it, and nothing is raised (a cell at 1,010 m
    // under a target of 1,040 m keeps its 1,010).
    const float land = 2000.0f;
    TroughBench bench(land);
    bench.parameters.troughs.width_km = 4.5 * bench.cell_km;
    bench.parameters.troughs.wall_m = 810.0;
    const int x0 = bench.width / 2, y0 = bench.height / 2;
    bench.axis(x0, y0, 1000.0f);
    bench.sea[bench.cell(x0 + 2, y0)] = 1;
    bench.axis(x0 - 2, y0, 1800.0f);
    bench.elevation[bench.cell(x0, y0 + 1)] = 1010.0f;
    const std::vector<float> before = bench.elevation;
    bench.carve();
    EXPECT_EQ(bench.elevation[bench.cell(x0 + 2, y0)], land);    // sea kept
    EXPECT_EQ(bench.elevation[bench.cell(x0 - 2, y0)], 1800.0f); // channel kept
    EXPECT_EQ(bench.elevation[bench.cell(x0, y0 + 1)], 1010.0f); // not raised to 1,040
    EXPECT_NEAR(bench.elevation[bench.cell(x0, y0 - 1)], 1040.0, 1.0e-3);
    EXPECT_NEAR(bench.elevation[bench.cell(x0 + 3, y0)], 1360.0, 1.0e-3);
    for (std::size_t cell = 0; cell < bench.count; cell++)
        ASSERT_LE(bench.elevation[cell], before[cell]);

    // An axis on a sea cell cuts nothing.
    TroughBench wet(land);
    wet.parameters.troughs.width_km = 4.5 * wet.cell_km;
    wet.axis(x0, y0, 1000.0f);
    wet.sea[wet.cell(x0, y0)] = 1;
    wet.carve();
    EXPECT_EQ(wet.reach(x0, y0, land), 0);

    // A bed at -50 m (the sea level 0): the targets -10 at d = 1 and 310 at d = 3; the first
    // stops at the next float above the sea level.
    TroughBench low(land);
    low.parameters.troughs.width_km = 4.5 * low.cell_km;
    low.parameters.troughs.wall_m = 810.0;
    low.axis(x0, y0, -50.0f);
    low.carve();
    const float floorlevel = std::nextafter(sealevel, std::numeric_limits<float>::infinity());
    EXPECT_EQ(low.elevation[low.cell(x0, y0 + 1)], floorlevel);
    EXPECT_GT(low.elevation[low.cell(x0, y0 + 1)], sealevel);
    EXPECT_NEAR(low.elevation[low.cell(x0, y0 + 3)], 310.0, 1.0e-3);
}

TEST_F(Landscape, RegionalOverlapTroughs)
{
    // The troughs on at the smallest halo the pass allows: the overlaps stay bit-identical
    // (compareoverlap asserts it), and the troughs act.
    const RegionalParameters parameters = troughparameters();
    const RegionalGlobalInputs columns = withglacialfields(syntheticworld(32));
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts counts = compareoverlap(columns, windowa(scale), windowb(scale), "troughs A/B", parameters);
        EXPECT_EQ(counts.shared, 4 * 5 * scale * scale);
        EXPECT_GT(counts.channel, 0);
        const RegionalFields on = landscape::runregionalpass(columns, windowb(scale), parameters);
        const RegionalFields off = landscape::runregionalpass(columns, windowb(scale), RegionalParameters{});
        const int changed = troughchangedcells(on, off);
        const int onoverlap = troughchangedonoverlap(on, off, windowa(scale));
        std::printf("[regional] troughs B, k %d: %d of %zu kept cells lowered, %d on the overlap with A\n", scale, changed,
                    on.elevation_m.size(), onoverlap);
        EXPECT_GT(changed, 0);
        EXPECT_GT(onoverlap, 0);
        for (std::size_t cell = 0; cell < on.elevation_m.size(); cell++)
            ASSERT_LE(on.elevation_m[cell], off.elevation_m[cell]);
    }
    const RegionalGlobalInputs rows = withglacialfields(syntheticworld(64));
    const OverlapCounts counts = compareoverlap(rows, windowc(4), windowd(4), "troughs C/D", parameters);
    EXPECT_GT(counts.sea, 0);
    EXPECT_GT(counts.lake, 0);
    EXPECT_GT(counts.channel, 0);
    const RegionalFields on = landscape::runregionalpass(rows, windowd(4), parameters);
    const RegionalFields off = landscape::runregionalpass(rows, windowd(4), RegionalParameters{});
    EXPECT_GT(troughchangedcells(on, off), 0);
    EXPECT_GT(troughchangedonoverlap(on, off, windowc(4)), 0);
}
