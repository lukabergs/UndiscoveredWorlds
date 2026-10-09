// Tests of regional_valleys.hpp (erosion plan Phase R2: the R6 contract and R9's incision
// and valley floors; the RegionalOverlapValleys* variants with R9's features on are in
// regional_pass_tests.cpp), suite Landscape: TEST_F(Landscape, RegionalValleys<What>) on
// hand-built global chains and fine channels and the synthetic world of regional_test_world.hpp.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_valleys.hpp"

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
