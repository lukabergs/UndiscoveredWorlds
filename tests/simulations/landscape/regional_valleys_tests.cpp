// Tests of regional_valleys.hpp (erosion plan Phase R2: the R6 contract and R9's incision
// and valley floors; the RegionalOverlapValleys* variants with R9's features on are in
// regional_pass_tests.cpp), suite Landscape: TEST_F(Landscape, RegionalValleys<What>) on
// hand-built global chains and fine channels and the synthetic world of regional_test_world.hpp.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_valleys.hpp"

#include <limits>
#include <utility>

TEST_F(Landscape, RegionalValleysAreOffByDefault)
{
    // No incision and no valleys at the defaults: the pass is the Phase R pass.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    EXPECT_FALSE(incision.active);
    EXPECT_TRUE(incision.incision_m.empty());
    EXPECT_TRUE(incision.filledelevation_m.empty());

    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    channels.direction.assign(count, 0);
    channels.bed_m.assign(count, 0.0f);
    channels.discharge_m3s.assign(count, 0.0f);
    channels.januaryflow.assign(count, 0);
    channels.julyflow.assign(count, 0);
    channels.globalcell.assign(count, -1);
    channels.channel[count / 2] = 1;
    channels.discharge_m3s[count / 2] = 5000.0f;
    const std::vector<uint8_t> sea(count, 0);
    std::vector<float> elevation(count);
    for (std::size_t cell = 0; cell < count; cell++)
        elevation[cell] = 100.0f + static_cast<float>(cell % 13);
    const std::vector<float> before = elevation;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, elevation);
    EXPECT_EQ(elevation, before);
}

namespace
{
// A global chain on an 8 x 4 world: row 1, columns 1 .. 5 land draining west (code 7)
// into the sea cell (0, 1); every other cell sea at -100 m. Filled elevations 5, 50, 100,
// 100, 100 m (cells 3 to 5 level: a donor's F equals its receiver's), sea level 0.
RegionalGlobalInputs chainworld()
{
    RegionalGlobalInputs inputs;
    inputs.width = 8;
    inputs.height = 4;
    inputs.radius_m = radius_m;
    inputs.sealevel = 0.0f;
    const std::size_t size = 32;
    landscape::GlobalNetwork& network = inputs.network;
    network.width = 8;
    network.height = 4;
    network.sea.assign(size, 1);
    network.direction.assign(size, 0);
    network.filledelevation_m.assign(size, -100.0f);
    network.drainagearea_km2.assign(size, 0.0f);
    network.discharge_m3s.assign(size, 0.0f);
    network.januaryflow.assign(size, 0);
    network.julyflow.assign(size, 0);
    const float filled[] = { 5.0f, 50.0f, 100.0f, 100.0f, 100.0f };
    for (int i = 1; i <= 5; i++)
    {
        const std::size_t cell = at(i, 1, 8);
        network.sea[cell] = 0;
        network.direction[cell] = 7;
        network.filledelevation_m[cell] = filled[i - 1];
    }
    inputs.elevation_m = network.filledelevation_m;
    return inputs;
}

std::vector<float> chainvalues(const std::vector<float>& field)
{
    std::vector<float> values;
    for (int i = 1; i <= 5; i++)
        values.push_back(field[at(i, 1, 8)]);
    return values;
}

// FineChannels on a domain of `count` cells, no channel set.
landscape::FineChannels nochannels(std::size_t count)
{
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    channels.direction.assign(count, 0);
    channels.bed_m.assign(count, 0.0f);
    channels.discharge_m3s.assign(count, 0.0f);
    channels.januaryflow.assign(count, 0);
    channels.julyflow.assign(count, 0);
    channels.globalcell.assign(count, -1);
    return channels;
}

// A straight channel along one domain row, its cells at the bed.
void addchannelrow(landscape::FineChannels& channels, std::vector<float>& elevation, int width, int row, float bed,
                   float discharge)
{
    for (int x = 0; x < width; x++)
    {
        const std::size_t cell = at(x, row, width);
        channels.channel[cell] = 1;
        channels.direction[cell] = 7;
        channels.bed_m[cell] = bed;
        channels.discharge_m3s[cell] = discharge;
        elevation[cell] = bed;
    }
}
}

TEST_F(Landscape, RegionalValleysIncisionIsConsistentDownstream)
{
    // Only the head cell (5) is uplifted, 30 m/Myr over 2 Myr: I = 60 m there. Cells 4 and
    // 3 stand level with it, so I' = 60 m propagates (F' = 40 m); cell 2 lies 50 m lower,
    // I' = 10 m (F' = 40 m); cell 1 lies 45 m lower still, I' = 0.
    RegionalGlobalInputs inputs = chainworld();
    inputs.upliftrate_m_myr.assign(32, 0.0f);
    inputs.upliftrate_m_myr[at(5, 1, 8)] = 30.0f;
    inputs.upliftrate_m_myr[at(4, 1, 8)] = -20.0f;  // subsidence counts as 0
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 2.0;
    parameters.valleys.upliftscale = 1.0;
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    ASSERT_TRUE(incision.active);
    EXPECT_EQ(chainvalues(incision.filledelevation_m), (std::vector<float>{ 5.0f, 40.0f, 40.0f, 40.0f, 40.0f }));
    EXPECT_EQ(chainvalues(incision.incision_m), (std::vector<float>{ 0.0f, 10.0f, 60.0f, 60.0f, 60.0f }));
    EXPECT_EQ(chainvalues(incision.rate_m_myr), (std::vector<float>{ 0.0f, 0.0f, 0.0f, 0.0f, 30.0f }));
    for (std::size_t cell = 0; cell < 32; cell++)
        if (inputs.network.sea[cell] != 0)
        {
            EXPECT_EQ(incision.incision_m[cell], 0.0f);
            EXPECT_EQ(incision.filledelevation_m[cell], inputs.network.filledelevation_m[cell]);
            EXPECT_EQ(incision.rate_m_myr[cell], 0.0f);
        }

    // The denudation rate (denudation / erosion duration) is the other source, the larger
    // one counts: 400 m over 20 Myr is 20 m/Myr at cell 5, under its uplift; 1000 m is
    // 50 m/Myr at cell 2, 100 m in 2 Myr, capped at F less the next float above the sea
    // level (F' = s+, I' = 50 m after rounding); cell 1 then follows down to s+ (I' = 5 m),
    // so F' never rises downstream.
    inputs.denudation_m.assign(32, 0.0f);
    inputs.denudation_m[at(5, 1, 8)] = 400.0f;
    inputs.denudation_m[at(2, 1, 8)] = 1000.0f;
    inputs.erosionduration_myr = 20.0;
    parameters.valleys.denudationratescale = 1.0;
    const landscape::GlobalIncision both = landscape::globalincision(inputs, parameters);
    ASSERT_TRUE(both.active);
    EXPECT_EQ(chainvalues(both.rate_m_myr), (std::vector<float>{ 0.0f, 50.0f, 0.0f, 0.0f, 30.0f }));
    const float above = std::nextafter(0.0f, 1.0f);
    EXPECT_EQ(chainvalues(both.filledelevation_m), (std::vector<float>{ above, above, 40.0f, 40.0f, 40.0f }));
    EXPECT_EQ(chainvalues(both.incision_m), (std::vector<float>{ 5.0f, 50.0f, 60.0f, 60.0f, 60.0f }));
}

TEST_F(Landscape, RegionalValleysIncisionStopsAtTheSeaLevel)
{
    // 1 km/Myr over 1 Myr everywhere (the local relief, 105 to 200 m against the -100 m
    // sea, does not bind): every channel cuts to the next float above the sea level, so no
    // land channel sits at the sea level itself.
    RegionalGlobalInputs inputs = chainworld();
    inputs.upliftrate_m_myr.assign(32, 1000.0f);
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 1.0;
    parameters.valleys.upliftscale = 1.0;
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    ASSERT_TRUE(incision.active);
    const float above = std::nextafter(0.0f, 1.0f);
    EXPECT_EQ(chainvalues(incision.filledelevation_m), (std::vector<float>{ above, above, above, above, above }));
    EXPECT_EQ(chainvalues(incision.incision_m), (std::vector<float>{ 5.0f, 50.0f, 100.0f, 100.0f, 100.0f }));

    // A sea level 5 m above cell 1's filled elevation (bad input for the channels, which
    // reject it): the cell keeps its F and nothing upstream drops to the sea level.
    inputs.sealevel = 10.0f;
    const landscape::GlobalIncision high = landscape::globalincision(inputs, parameters);
    const float tenplus = std::nextafter(10.0f, 11.0f);
    EXPECT_EQ(chainvalues(high.filledelevation_m), (std::vector<float>{ 5.0f, tenplus, tenplus, tenplus, tenplus }));
}

TEST_F(Landscape, RegionalValleysIncisionStopsAtTheLocalRelief)
{
    // The base elevation of every sea cell set to 90 m (the filled elevations unchanged):
    // the local relief is 10 m at cells 4 and 5, 50 m at cells 2 and 3 and 85 m at cell 1,
    // so 1 km of incision is capped there: F - I = 90, 90, 50 m at cells 5 to 3, and cells
    // 2 and 1 reach the next float above the sea level (their F less the relief cap is
    // lower, the cap at F - s+ binds first); F' = 50 at cell 3 under its donors' 90.
    RegionalGlobalInputs inputs = chainworld();
    for (std::size_t cell = 0; cell < 32; cell++)
        if (inputs.network.sea[cell] != 0)
            inputs.elevation_m[cell] = 90.0f;
    inputs.upliftrate_m_myr.assign(32, 1000.0f);
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 1.0;
    parameters.valleys.upliftscale = 1.0;
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    ASSERT_TRUE(incision.active);
    const float above = std::nextafter(0.0f, 1.0f);
    EXPECT_EQ(chainvalues(incision.filledelevation_m), (std::vector<float>{ above, above, 50.0f, 90.0f, 90.0f }));
    EXPECT_EQ(chainvalues(incision.incision_m), (std::vector<float>{ 5.0f, 50.0f, 50.0f, 10.0f, 10.0f }));
}

TEST_F(Landscape, RegionalValleysIncisionInactiveWithoutItsFields)
{
    RegionalGlobalInputs inputs = chainworld();
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 20.0;
    parameters.valleys.upliftscale = 1.0;
    parameters.valleys.denudationratescale = 1.0;
    // Empty fields: no rate.
    EXPECT_FALSE(landscape::globalincision(inputs, parameters).active);
    // Denudation without an erosion duration, uplift that is all subsidence: no rate.
    inputs.denudation_m.assign(32, 500.0f);
    inputs.upliftrate_m_myr.assign(32, -5.0f);
    EXPECT_FALSE(landscape::globalincision(inputs, parameters).active);
    // With the duration: 25 m/Myr. Then no elapsed time, and the denudation's scale 0.
    inputs.erosionduration_myr = 20.0;
    EXPECT_TRUE(landscape::globalincision(inputs, parameters).active);
    parameters.relief.elapsed_myr = 0.0;
    EXPECT_FALSE(landscape::globalincision(inputs, parameters).active);
    parameters.relief.elapsed_myr = 20.0;
    parameters.valleys.denudationratescale = 0.0;
    EXPECT_FALSE(landscape::globalincision(inputs, parameters).active);
}

TEST_F(Landscape, RegionalValleysFloorAndSides)
{
    // A straight channel along domain row 20 with bed 100 m and Q = 4000 m3/s: valley_floor_km
    // a hair over one fine cell (dy, 313 km at k = 4 on this coarse test world) gives the
    // half-width w = 1.0001 x sqrt(4000 / 1000) = 2.0002 fine cells, so rows 18 to 22 drop to
    // 100 m; with sides of slope s within 4 cells (k = 4), rows 17 and 23 drop to
    // 100 + (3 - w) x s x dy and rows 16 and 24 to 100 + (4 - w) s dy; rows beyond keep 1000 m.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    ASSERT_GT(height, 40);
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    landscape::FineChannels channels = nochannels(count);
    std::vector<float> elevation(count, 1000.0f);
    addchannelrow(channels, elevation, width, 20, 100.0f, 4000.0f);
    const std::vector<uint8_t> sea(count, 0);
    const double dy = landscape::finepatchmetric(window, inputs.radius_m).dy_m;
    RegionalParameters parameters;
    parameters.valleys.floor_km = 1.0001 * dy / 1000.0;
    parameters.valleys.reference_m3s = 1000.0;
    parameters.valleys.floormax_km = 4.0001 * dy / 1000.0;
    const double w = 2.0 * 1.0001;

    // No sides.
    std::vector<float> floor = elevation;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, floor);
    const int x = width / 2;
    for (int y = 10; y <= 30; y++)
        EXPECT_EQ(floor[at(x, y, width)], std::abs(y - 20) <= 2 ? 100.0f : 1000.0f) << "row " << y;

    // Sides: s = 0.0005, about 156 m per fine cell of this coarse test world.
    parameters.valleys.sideslope = 0.0005;
    parameters.valleys.reachcells = 4;
    std::vector<float> sides = elevation;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, sides);
    for (int y = 10; y <= 30; y++)
    {
        const int d = std::abs(y - 20);
        const double expected = d <= 2 ? 100.0 : (d <= 4 ? 100.0 + (d - w) * 0.0005 * dy : 1000.0);
        EXPECT_NEAR(sides[at(x, y, width)], expected, 1e-4) << "row " << y;
    }
    // The reach is clamped to k: 40 cells ask for more than one global cell.
    parameters.valleys.reachcells = 40;
    parameters.valleys.sideslope = 0.0001;
    std::vector<float> clamped = elevation;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, clamped);
    EXPECT_LT(clamped[at(x, 16, width)], 1000.0f);
    EXPECT_EQ(clamped[at(x, 15, width)], 1000.0f);

    // The floor's width follows sqrt(Q) and its cap: Q = 250 gives 0.5 (no floor cell),
    // Q = 1e6 gives 31.6, capped at 4 (the cap in km, a hair over 4 cells).
    parameters.valleys.sideslope = 0.0;
    parameters.valleys.reachcells = 0;
    for (float discharge : { 250.0f, 1.0e6f })
    {
        landscape::FineChannels other = nochannels(count);
        std::vector<float> field(count, 1000.0f);
        addchannelrow(other, field, width, 20, 100.0f, discharge);
        landscape::carvevalleys(inputs, window, parameters, other, sea, field);
        const int halfwidth = discharge < 1000.0f ? 0 : 4;
        for (int y = 10; y <= 30; y++)
            EXPECT_EQ(field[at(x, y, width)], std::abs(y - 20) <= halfwidth ? 100.0f : 1000.0f)
                << "Q " << discharge << ", row " << y;
    }

    // The floor is set in km on the ground: the same values at k = 8 (cells of half the size)
    // give a half-width of 4.0004 fine cells for Q = 4000 m3/s, rows 16 to 24.
    const RegionalWindow fine = windowa(8);
    const int finewidth = fine.domainfinewidth();
    const std::size_t finecount = static_cast<std::size_t>(finewidth) * static_cast<std::size_t>(fine.domainfineheight());
    landscape::FineChannels finechannels = nochannels(finecount);
    std::vector<float> finefield(finecount, 1000.0f);
    addchannelrow(finechannels, finefield, finewidth, 20, 100.0f, 4000.0f);
    landscape::carvevalleys(inputs, fine, parameters, finechannels, std::vector<uint8_t>(finecount, 0), finefield);
    for (int y = 10; y <= 30; y++)
        EXPECT_EQ(finefield[at(finewidth / 2, y, finewidth)], std::abs(y - 20) <= 4 ? 100.0f : 1000.0f) << "k 8, row " << y;
}

TEST_F(Landscape, RegionalValleysLowestTargetSeaAndChannelsKept)
{
    // Two channels 4 rows apart, beds 300 m (row 20) and 100 m (row 24), floors of 2
    // cells: rows 22 and 23 lie within both floors and take the lower bed. A channel cell
    // and a sea cell inside a floor keep their values. A third channel with its bed below
    // the sea level (row 32) floors its neighbours at the next float above the sea level.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const int width = window.domainfinewidth();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(window.domainfineheight());
    landscape::FineChannels channels = nochannels(count);
    std::vector<float> elevation(count, 1000.0f);
    addchannelrow(channels, elevation, width, 20, 300.0f, 4000.0f);
    addchannelrow(channels, elevation, width, 24, 100.0f, 4000.0f);
    addchannelrow(channels, elevation, width, 32, -5.0f, 4000.0f);
    std::vector<uint8_t> sea(count, 0);
    const int x = width / 2;
    sea[at(x, 21, width)] = 1;
    elevation[at(x, 21, width)] = -50.0f;
    RegionalParameters parameters;
    const double dy = landscape::finepatchmetric(window, inputs.radius_m).dy_m;
    parameters.valleys.floor_km = 1.0001 * dy / 1000.0;
    parameters.valleys.floormax_km = 4.0001 * dy / 1000.0;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, elevation);
    EXPECT_EQ(elevation[at(x, 17, width)], 1000.0f);
    EXPECT_EQ(elevation[at(x, 18, width)], 300.0f);
    EXPECT_EQ(elevation[at(x, 19, width)], 300.0f);
    EXPECT_EQ(elevation[at(x, 20, width)], 300.0f);  // a channel keeps its bed
    EXPECT_EQ(elevation[at(x, 21, width)], -50.0f);  // sea
    EXPECT_EQ(elevation[at(x + 1, 21, width)], 300.0f);
    EXPECT_EQ(elevation[at(x, 22, width)], 100.0f);  // both floors: the lower
    EXPECT_EQ(elevation[at(x, 23, width)], 100.0f);
    EXPECT_EQ(elevation[at(x, 26, width)], 100.0f);
    EXPECT_EQ(elevation[at(x, 27, width)], 1000.0f);
    EXPECT_EQ(elevation[at(x, 29, width)], 1000.0f);
    const float above = std::nextafter(0.0f, 1.0f);
    for (int y : { 30, 31, 33, 34 })
        EXPECT_EQ(elevation[at(x, y, width)], above) << "row " << y;
    EXPECT_EQ(elevation[at(x, 32, width)], -5.0f);
}

namespace
{
// One channel cell at the domain's centre (bed 100 m, Q = 1000 m3/s), every other cell at 1000 m.
struct PointChannel
{
    RegionalWindow window;
    int width = 0;
    int x = 0, y = 0;
    landscape::FineChannels channels;
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    double dy = 0.0;
};

PointChannel pointchannel(const RegionalGlobalInputs& inputs, int scale, int row = -1)
{
    PointChannel point;
    point.window = windowa(scale);
    point.width = point.window.domainfinewidth();
    const std::size_t count = static_cast<std::size_t>(point.width) * static_cast<std::size_t>(point.window.domainfineheight());
    point.channels = nochannels(count);
    point.elevation.assign(count, 1000.0f);
    point.sea.assign(count, 0);
    point.x = point.width / 2;
    point.y = row >= 0 ? row : point.window.domainfineheight() / 2;
    const std::size_t cell = at(point.x, point.y, point.width);
    point.channels.channel[cell] = 1;
    point.channels.direction[cell] = 7;
    point.channels.bed_m[cell] = 100.0f;
    point.channels.discharge_m3s[cell] = 1000.0f;
    point.elevation[cell] = 100.0f;
    point.dy = landscape::finepatchmetric(point.window, inputs.radius_m).dy_m;
    return point;
}
}

TEST_F(Landscape, RegionalValleysRoundFloorInTheLoop)
{
    // Multiscale plan S3.5: in the refinement loop (bandupper > 0) the floor is a disc of radius
    // w on the ground: a cell (dx, dy) away lies on it when (dx r)^2 + dy^2 <= w^2, r the axis
    // row's east-west over north-south cell size (the cosine of its latitude), so at high latitude
    // the floor spans more columns than rows (up to k); the single pass keeps R9's Chebyshev
    // square in cells. w = 3.5 cells (Q at the reference), no sides; at the domain's middle row and
    // at a row where r <= 0.8.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const landscape::PatchMetric metric = landscape::finepatchmetric(windowa(4), inputs.radius_m);
    int northern = -1;
    for (int y = 6; y < metric.height - 6 && northern < 0; y++)
        if (landscape::eastwestratio(metric, y) <= 0.8 && landscape::eastwestratio(metric, y) > 0.6)
            northern = y;
    ASSERT_GE(northern, 0);
    for (const int row : { -1, northern })
    {
        PointChannel point = pointchannel(inputs, 4, row);
        const double ratio = landscape::eastwestratio(metric, point.y);
        SCOPED_TRACE("row " + std::to_string(point.y) + ", ratio " + std::to_string(ratio));
        RegionalParameters parameters;
        parameters.valleys.floor_km = 3.5 * point.dy / 1000.0;
        parameters.valleys.reference_m3s = 1000.0;
        parameters.valleys.floormax_km = 3.9 * point.dy / 1000.0;
        std::vector<float> loop = point.elevation, single = point.elevation;
        landscape::carvevalleys(inputs, point.window, parameters, point.channels, point.sea, loop, 1.0,
                                std::numeric_limits<double>::infinity());
        landscape::carvevalleys(inputs, point.window, parameters, point.channels, point.sea, single);
        int floorloop = 0, floorsingle = 0, expected = 0;
        for (int dy = -6; dy <= 6; dy++)
            for (int dx = -6; dx <= 6; dx++)
            {
                const double east = dx * ratio;
                const bool disc = east * east + dy * dy <= 3.5 * 3.5 && std::abs(dx) <= 4;
                const bool square = std::abs(dx) <= 3 && std::abs(dy) <= 3;
                const float inloop = loop[at(point.x + dx, point.y + dy, point.width)];
                const float insingle = single[at(point.x + dx, point.y + dy, point.width)];
                EXPECT_EQ(inloop, disc ? 100.0f : 1000.0f) << "loop " << dx << ", " << dy;
                EXPECT_EQ(insingle, square ? 100.0f : 1000.0f) << "single " << dx << ", " << dy;
                expected += disc ? 1 : 0;
                floorloop += inloop == 100.0f ? 1 : 0;
                floorsingle += insingle == 100.0f ? 1 : 0;
            }
        EXPECT_EQ(floorloop, expected);
        EXPECT_EQ(floorsingle, 49);  // the 7 x 7 square
        // Round: (3, 3) at 4.2 north-south cells is off the floor in the loop.
        EXPECT_EQ(loop[at(point.x + 3, point.y + 3, point.width)], 1000.0f);
        if (row >= 0)
            EXPECT_EQ(loop[at(point.x + 4, point.y, point.width)], 100.0f);  // 4 r <= 3.2 cells east
    }
}

TEST_F(Landscape, RegionalValleysReachInKmAgainstTheBand)
{
    // Multiscale plan S3.5: in the loop valleys.reach_km replaces reachcells (40 here, ignored)
    // and the full width max(floor half-width, reach) is tested against both ends of the band.
    // A channel row with a floor of 0.2 cells and sides of slope s reaching 2.6 cells at k = 4:
    // carved in the band [1, 4) (rows 1 and 2 away at the bed + (d - 0.2) s dy, row 3 kept), not
    // in [4, 16) nor [0.25, 1). The reach is km on the ground: at k = 8 the same reach is 5.2 cells,
    // carved in [4, 16) only. The single pass never reads reach_km (reachcells clamped to k).
    const RegionalGlobalInputs inputs = syntheticworld(32);
    double reach = 2.6;  // in fine cells of k = 4
    const auto run = [&](int scale, double lower, double upper, double& dy, int& width, int& middle) {
        const RegionalWindow window = windowa(scale);
        width = window.domainfinewidth();
        const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(window.domainfineheight());
        landscape::FineChannels channels = nochannels(count);
        std::vector<float> elevation(count, 1000.0f);
        middle = window.domainfineheight() / 2;
        addchannelrow(channels, elevation, width, middle, 100.0f, 1000.0f);
        const double cell_m = landscape::finepatchmetric(windowa(4), inputs.radius_m).dy_m;
        dy = landscape::finepatchmetric(window, inputs.radius_m).dy_m;
        RegionalParameters parameters;
        parameters.valleys.floor_km = 0.2 * cell_m / 1000.0;
        parameters.valleys.reference_m3s = 1000.0;
        parameters.valleys.floormax_km = 4.0 * cell_m / 1000.0;
        parameters.valleys.sideslope = 0.0005;
        parameters.valleys.reachcells = 40;
        parameters.valleys.reach_km = reach * cell_m / 1000.0;
        landscape::carvevalleys(inputs, window, parameters, channels, std::vector<uint8_t>(count, 0), elevation, lower,
                                upper);
        return elevation;
    };
    double dy = 0.0;
    int width = 0, middle = 0;
    {
        const std::vector<float> carved = run(4, 1.0, 4.0, dy, width, middle);
        const int x = width / 2;
        const double w = 0.2;
        for (int d = 1; d <= 2; d++)
            for (int side : { -1, 1 })
                EXPECT_NEAR(carved[at(x, middle + side * d, width)], 100.0 + (d - w) * 0.0005 * dy, 1e-3) << "row " << d;
        EXPECT_EQ(carved[at(x, middle + 3, width)], 1000.0f);
        EXPECT_EQ(carved[at(x, middle - 3, width)], 1000.0f);
        EXPECT_EQ(carved[at(x, middle, width)], 100.0f);
    }
    for (const auto& [lower, upper] : std::vector<std::pair<double, double>>{{4.0, 16.0}, {0.25, 1.0}})
    {
        const std::vector<float> carved = run(4, lower, upper, dy, width, middle);
        int changed = 0;
        for (float value : carved)
            changed += value != 1000.0f && value != 100.0f ? 1 : 0;
        EXPECT_EQ(changed, 0) << "band [" << lower << ", " << upper << ")";
    }
    {
        const std::vector<float> carved = run(8, 1.0, 4.0, dy, width, middle);
        int changed = 0;
        for (float value : carved)
            changed += value != 1000.0f && value != 100.0f ? 1 : 0;
        EXPECT_EQ(changed, 0) << "k 8, band [1, 4)";
    }
    {
        const std::vector<float> carved = run(8, 4.0, 16.0, dy, width, middle);
        const int x = width / 2;
        const double w = 0.4;  // the same floor in km: twice the cells at k = 8
        for (int d = 1; d <= 5; d++)
            EXPECT_NEAR(carved[at(x, middle + d, width)], 100.0 + (d - w) * 0.0005 * dy, 1e-3) << "k 8, row " << d;
        EXPECT_EQ(carved[at(x, middle + 6, width)], 1000.0f);
    }
    {
        // The single pass: reachcells 40 clamped to k = 4, reach_km not read.
        const std::vector<float> carved = run(4, 0.0, 0.0, dy, width, middle);
        const int x = width / 2;
        EXPECT_LT(carved[at(x, middle + 4, width)], 1000.0f);
        EXPECT_EQ(carved[at(x, middle + 5, width)], 1000.0f);
    }
    {
        // Level 1 (upper end infinite): a reach of 6 cells, above any band, is carved and clamped to
        // k = 4 cells (one global cell, the halo rule); with the upper end 4 it is not carved.
        reach = 6.0;
        const std::vector<float> carved = run(4, 1.0, std::numeric_limits<double>::infinity(), dy, width, middle);
        const int x = width / 2;
        for (int d = 1; d <= 4; d++)
            EXPECT_NEAR(carved[at(x, middle + d, width)], 100.0 + (d - 0.2) * 0.0005 * dy, 1e-3) << "level 1, row " << d;
        EXPECT_EQ(carved[at(x, middle + 5, width)], 1000.0f);
        const std::vector<float> above = run(4, 1.0, 4.0, dy, width, middle);
        int changed = 0;
        for (float value : above)
            changed += value != 1000.0f && value != 100.0f ? 1 : 0;
        EXPECT_EQ(changed, 0) << "reach above the band [1, 4)";
    }
}

TEST_F(Landscape, RegionalValleysChannelCellKeepsItsMean)
{
    // Multiscale plan S3.5: a form never changes a channel cell; channelcellmeans then lowers a land
    // channel cell next to a carved cell to its cell mean on the carved field, b + min(h, cap)
    // (regional_solve.cpp's), never back to its bed under sloping sides and never up. A channel
    // row at bed 100 m whose cells sit at their mean 100 + min(900, 300) = 400 m before the form;
    // sides of slope s within 2.5 cells and a floor of 0.2 cells put rows 1 away at 100 + 0.8 s dy,
    // so the channel cells' mean is 100 + 0.8 s dy, above the bed.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    landscape::FineChannels channels = nochannels(count);
    std::vector<float> elevation(count, 1000.0f);
    const int middle = height / 2;
    addchannelrow(channels, elevation, width, middle, 100.0f, 1000.0f);
    const double cap = 300.0;
    for (int x = 0; x < width; x++)
        elevation[at(x, middle, width)] = 400.0f;
    elevation[at(3, middle, width)] = 150.0f;  // below its mean: never raised
    const std::vector<uint8_t> sea(count, 0);
    const double dy = landscape::finepatchmetric(window, inputs.radius_m).dy_m;
    RegionalParameters parameters;
    parameters.valleys.floor_km = 0.2 * dy / 1000.0;
    parameters.valleys.reference_m3s = 1000.0;
    parameters.valleys.floormax_km = 4.0 * dy / 1000.0;
    parameters.valleys.sideslope = 0.0005;
    parameters.valleys.reach_km = 2.5 * dy / 1000.0;
    const std::vector<float> before = elevation;
    std::vector<float> carved = elevation;
    landscape::carvevalleys(inputs, window, parameters, channels, sea, carved, 1.0, 4.0);
    for (int x = 0; x < width; x++)
        EXPECT_EQ(carved[at(x, middle, width)], before[at(x, middle, width)]) << "the form keeps channel column " << x;
    std::vector<float> nocap = carved;
    landscape::channelcellmeans(window, channels, sea, 0.0, before, nocap);
    EXPECT_EQ(nocap, carved);
    landscape::channelcellmeans(window, channels, sea, cap, before, carved);
    const double side = 0.8 * 0.0005 * dy;
    ASSERT_LT(side, cap);
    for (int x = 0; x < width; x++)
    {
        // The cell mean on the final field: the non-channel 8-neighbours' mean height above the bed.
        double sum = 0.0;
        int neighbours = 0;
        for (int dy8 = -1; dy8 <= 1; dy8 += 2)
            for (int dx = -1; dx <= 1; dx++)
            {
                if (x + dx < 0 || x + dx >= width)
                    continue;
                sum += static_cast<double>(carved[at(x + dx, middle + dy8, width)]) - 100.0;
                neighbours++;
            }
        const float mean = static_cast<float>(100.0 + std::min(std::max(0.0, sum / neighbours), cap));
        const float value = carved[at(x, middle, width)];
        if (x == 3)
        {
            EXPECT_EQ(value, 150.0f);
            continue;
        }
        EXPECT_EQ(value, mean) << "column " << x;
        EXPECT_GE(value, mean) << "column " << x;
        EXPECT_NEAR(value, 100.0 + side, 1e-3) << "column " << x;
        EXPECT_GT(value, 100.0f) << "column " << x;
    }
}

namespace
{
// The cells two level windows' kept boxes share, as index pairs into their kept fields.
std::vector<std::pair<std::size_t, std::size_t>> sharedlevelcells(const landscape::LevelWindow& a, const landscape::LevelWindow& b)
{
    std::vector<std::pair<std::size_t, std::size_t>> shared;
    const std::int64_t width = a.levelwidth();
    for (int y = 0; y < a.rows; y++)
        for (int x = 0; x < a.columns; x++)
        {
            const std::int64_t column = landscape::wraplevel(a.column0 + x, width), row = a.row0 + y;
            const std::int64_t xb = landscape::wraplevel(column - b.column0, width), yb = row - b.row0;
            if (xb >= b.columns || yb < 0 || yb >= b.rows)
                continue;
            shared.push_back({at(x, y, a.columns), at(static_cast<int>(xb), static_cast<int>(yb), b.columns)});
        }
    return shared;
}

template<typename T>
int differingcells(const std::vector<T>& a, const std::vector<T>& b,
                   const std::vector<std::pair<std::size_t, std::size_t>>& cells)
{
    int count = 0;
    for (const auto& [ca, cb] : cells)
        count += std::memcmp(&a[ca], &b[cb], sizeof(T)) != 0 ? 1 : 0;
    return count;
}
}

TEST_F(Landscape, RegionalValleysOverlapFormsPerLevel)
{
    // Multiscale plan S3.5: two levels at factor 4 (k = 16) on the 256-column synthetic world with
    // R9's incision, the valleys by reach in km (floors of 20 km x sqrt(Q / 1e3 m3/s) up to 40 km,
    // sides reaching 25 km: 2.56 level-2 cells), troughs (30 km), cirques (15 km: 1.54 level-2
    // cells), the channel cells' means (theta 0.3), the solve against its parent and the smoothing:
    // on windows P and Q (different column and row origins, Q across the date line) every level's
    // shared kept elevation and the final level's sea, lakes and rivers bit-identical; the forms act.
    RegionalGlobalInputs inputs = withglacialfields(syntheticworld(256));
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.upliftrate_m_myr.resize(size);
    inputs.denudation_m.resize(size);
    inputs.erosionduration_myr = 20.0;
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            const double c = std::cos(latitude_rad(j, inputs.height));
            inputs.upliftrate_m_myr[at(i, j, inputs.width)] = static_cast<float>(40.0 * c * c - 5.0);
            inputs.denudation_m[at(i, j, inputs.width)] = static_cast<float>(300.0 + 200.0 * std::sin(longitude_rad(i, inputs.width)));
        }
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 20.0;
    parameters.valleys.upliftscale = 1.0;
    parameters.valleys.denudationratescale = 1.0;
    parameters.valleys.floor_km = 20.0;
    parameters.valleys.reference_m3s = 1.0e3;
    parameters.valleys.floormax_km = 40.0;
    parameters.valleys.sideslope = 0.02;
    parameters.valleys.reachcells = 16;
    parameters.valleys.reach_km = 25.0;
    parameters.troughs.width_km = 30.0;
    parameters.troughs.widthmax_km = 38.0;  // 3.9 level-2 cells: the capped troughs at level 2
    parameters.troughs.glaciationmin = 0.1;
    parameters.cirques.radius_km = 15.0;
    parameters.cirques.candidates = 4;
    landscape::RefinementParameters refinement;
    refinement.smoothpasses = 6;
    refinement.parentlift = 2;
    refinement.channelcellslope = 0.3;
    refinement.prolongationlimit = 1;
    const int halo = landscape::refinementhalo(2, refinement.smoothpasses);
    const RegionalWindow windowp = landscape::regionalwindow(256, 16, halo, 160.0, -15.0, 180.0, 10.0);
    const RegionalWindow windowq = landscape::regionalwindow(256, 16, halo, 170.0, -10.0, -170.0, 15.0);
    const auto run = [&](const RegionalWindow& window, const RegionalParameters& settings) {
        landscape::LevelTarget target;
        target.window = window;
        target.level = 2;
        return landscape::runrefinement(inputs, target, settings, refinement);
    };
    const auto report = [](const landscape::RefinementResult& result, const char* label) {
        for (const landscape::LevelDiagnostics& level : result.levels)
            std::printf("[valleys] %s, level %d: valley %d, trough %d, cirque %d cells; lambda before %.1f m; forced %d\n",
                        label, level.level, level.valleycells, level.troughcells, level.cirquecells,
                        level.smoothing.maxlambdabefore_m, level.projection.forcedparents);
    };
    // Every level's shared kept elevation and the final level's fields bit-identical; returns the
    // final level's shared cells.
    const auto compare = [](const landscape::RefinementResult& p, const landscape::RefinementResult& q,
                            const std::string& label) {
        EXPECT_EQ(p.levels.size(), 2u) << label;
        EXPECT_EQ(q.levels.size(), 2u) << label;
        const auto first = sharedlevelcells(p.levels[0].window, q.levels[0].window);
        EXPECT_GT(first.size(), 0u) << label;
        EXPECT_EQ(differingcells(p.coarserelevation_m[0], q.coarserelevation_m[0], first), 0) << label << ": level 1";
        const auto cells = sharedlevelcells(p.levels[1].window, q.levels[1].window);
        EXPECT_GT(cells.size(), 0u) << label;
        int channels = 0;
        for (const auto& [ca, cb] : cells)
            channels += p.fields.riverdirection[ca] != 0 ? 1 : 0;
        EXPECT_GT(channels, 0) << label;
        EXPECT_EQ(differingcells(p.fields.elevation_m, q.fields.elevation_m, cells), 0) << label << ": elevation";
        EXPECT_EQ(differingcells(p.fields.sea, q.fields.sea, cells), 0) << label << ": sea";
        EXPECT_EQ(differingcells(p.fields.lakesurface_m, q.fields.lakesurface_m, cells), 0) << label << ": lakes";
        EXPECT_EQ(differingcells(p.fields.riverdirection, q.fields.riverdirection, cells), 0) << label << ": rivers";
        EXPECT_EQ(differingcells(p.fields.discharge_m3s, q.fields.discharge_m3s, cells), 0) << label << ": discharge";
        EXPECT_EQ(differingcells(p.fields.icethickness_m, q.fields.icethickness_m, cells), 0) << label << ": ice";
        return cells;
    };
    {
        const landscape::RefinementResult p = run(windowp, parameters);
        const landscape::RefinementResult q = run(windowq, parameters);
        report(p, "every form P");
        ASSERT_EQ(p.levels.size(), 2u);
        int valleys = 0, troughs = 0, cirques = 0;
        for (const landscape::LevelDiagnostics& level : p.levels)
        {
            valleys += level.valleycells;
            troughs += level.troughcells;
            cirques += level.cirquecells;
        }
        EXPECT_GT(p.levels[1].valleycells, 0);
        EXPECT_GT(p.levels[1].troughcells, 0);
        EXPECT_GT(troughs, 0);
        EXPECT_GT(cirques, 0);
        EXPECT_GT(valleys, 0);
        const auto cells = compare(p, q, "every form P/Q");
        // The troughs act inside the shared box: shared land cells more than 1 m below the run without them
        // (the projection puts most of a level-2 trough's volume back inside its parent cell, so few stay 50 m lower).
        RegionalParameters notroughs = parameters;
        notroughs.troughs.width_km = 0.0;
        const landscape::RefinementResult plain = run(windowp, notroughs);
        int lowered = 0, lowered1 = 0, raised1 = 0;
        for (const auto& [ca, cb] : cells)
        {
            const float with = p.fields.elevation_m[ca], without = plain.fields.elevation_m[ca];
            lowered += p.fields.sea[ca] == 0 && with < without - 50.0f ? 1 : 0;
            lowered1 += p.fields.sea[ca] == 0 && with < without - 1.0f ? 1 : 0;
            raised1 += p.fields.sea[ca] == 0 && with > without + 1.0f ? 1 : 0;
        }
        int keptlowered = 0;
        for (std::size_t cell = 0; cell < p.fields.elevation_m.size(); cell++)
            keptlowered += p.fields.sea[cell] == 0 && p.fields.elevation_m[cell] < plain.fields.elevation_m[cell] - 50.0f ? 1 : 0;
        std::printf("[valleys] every form P: %d of %zu shared cells more than 50 m below the run without troughs (1 m: %d "
                    "lowered, %d raised; the whole kept box: %d of %zu more than 50 m)\n",
                    lowered, cells.size(), lowered1, raised1, keptlowered, p.fields.elevation_m.size());
        EXPECT_GT(lowered1, 1000);
    }
    {
        // The level-1 clamp: sides reaching 200 km (5.1 level-1 cells, above every band) are carved at
        // level 1 only (its band has no upper end), clamped to one global cell, and the overlap holds.
        RegionalParameters wide = parameters;
        wide.valleys.reach_km = 200.0;
        wide.troughs.width_km = 0.0;
        wide.cirques.radius_km = 0.0;
        const landscape::RefinementResult p = run(windowp, wide);
        const landscape::RefinementResult q = run(windowq, wide);
        report(p, "valleys reaching 200 km P");
        ASSERT_EQ(p.levels.size(), 2u);
        EXPECT_GT(p.levels[0].valleycells, 0);
        EXPECT_EQ(p.levels[1].valleycells, 0);
        compare(p, q, "valleys reaching 200 km P/Q");
    }
}
