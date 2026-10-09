// Tests of regional_valleys.hpp (erosion plan Phase R2: the R6 contract; R9 adds the
// incision's and the valley floors' tests and the RegionalOverlapValleys* variants with its
// features on), suite Landscape: TEST_F(Landscape, RegionalValleys<What>) on the synthetic
// world of regional_test_world.hpp.

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
