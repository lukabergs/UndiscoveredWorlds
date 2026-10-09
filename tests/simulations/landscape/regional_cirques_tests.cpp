// Tests of landscape/regional_cirques.hpp (erosion plan Phase E5: the E5.0 contract, which
// step E5.2 implements: cirques, the parabolic bowls on glaciated high ground), suite
// Landscape: TEST_F(Landscape, RegionalCirques<What>) on the synthetic world of
// regional_test_world.hpp with withglacialfields. E5.2 adds the placement and bowl tests
// and the overlap variants with the cirques on (RegionalOverlapCirques*), on a world wider
// than the widened box.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_cirques.hpp"
#include "simulations/landscape/regional_glacial.hpp"

TEST_F(Landscape, RegionalCirquesAreOffByDefault)
{
    // At the defaults (radius_km 0) the hook changes nothing, with or without the fields.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    RegionalParameters parameters;
    parameters.troughs.width_km = 10.0;  // another feature on: the fields are active, the cirques still off
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    const std::vector<uint8_t> sea(count, 0);
    std::vector<float> elevation(count, 3000.0f);
    landscape::carvecirques(inputs, window, parameters, glacial, channels, sea, elevation);
    for (std::size_t cell = 0; cell < count; cell++)
        ASSERT_EQ(elevation[cell], 3000.0f);
}
