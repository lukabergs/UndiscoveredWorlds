// Tests of landscape/regional_ice.hpp (erosion plan Phase E5: the E5.0 contract, which step
// E5.3 implements: present-day ice on the fine cells above the present climate's ELA),
// suite Landscape: TEST_F(Landscape, RegionalIce<What>) on the synthetic world of
// regional_test_world.hpp with withglacialfields. E5.3 adds the rule's tests and the
// overlap variants with the ELA ice on (RegionalOverlapIce*), on a world wider than the
// widened box.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_ice.hpp"

TEST_F(Landscape, RegionalIceAboveTheElaIsOffByDefault)
{
    // At the defaults (thicknessfraction 0) the hook changes nothing, with or without the fields.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    RegionalParameters parameters;
    parameters.troughs.width_km = 10.0;  // another feature on: the fields are active, the ELA ice still off
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    const std::vector<uint8_t> sea(count, 0);
    const std::vector<float> elevation(count, 9000.0f);  // far above every ELA
    std::vector<float> ice(count, 0.0f);
    landscape::iceabovetheela(inputs, window, parameters, glacial, sea, elevation, ice);
    for (std::size_t cell = 0; cell < count; cell++)
        ASSERT_EQ(ice[cell], 0.0f);
}
