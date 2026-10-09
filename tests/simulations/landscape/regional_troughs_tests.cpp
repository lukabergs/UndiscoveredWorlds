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
