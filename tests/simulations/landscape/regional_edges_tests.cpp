// Tests of regional_edges.hpp (erosion plan Phase R2: the R6 contract; R8 adds the smooth
// gate's and the straight-run breaker's tests and the RegionalOverlapEdges* variants with its
// features on, on a world wider than the land-share rule's widened box), suite Landscape:
// TEST_F(Landscape, RegionalEdges<What>) on the synthetic world of regional_test_world.hpp.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_edges.hpp"

namespace
{
using landscape::FineWarp;
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
