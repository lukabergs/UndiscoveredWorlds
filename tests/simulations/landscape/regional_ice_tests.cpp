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

namespace
{
// The domain's cell count of a window.
std::size_t domaincount(const RegionalWindow& window)
{
    return static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
}

// The ELA ice on, edge noise as given.
RegionalParameters iceparameters(double fraction, double thickest, double edgenoise)
{
    RegionalParameters parameters;
    parameters.ice.thicknessfraction = fraction;
    parameters.ice.thicknessmax_m = thickest;
    parameters.ice.edgenoise_m = edgenoise;
    return parameters;
}

// A surface at the present ELA plus `above` metres on every cell.
std::vector<float> abovetheela(const landscape::FineGlacial& glacial, double above)
{
    std::vector<float> elevation(glacial.presentela_m.size());
    for (std::size_t cell = 0; cell < elevation.size(); cell++)
        elevation[cell] = static_cast<float>(static_cast<double>(glacial.presentela_m[cell]) + above);
    return elevation;
}

int icedcells(const std::vector<float>& ice)
{
    int iced = 0;
    for (float thickness : ice)
        iced += thickness > 0.0f ? 1 : 0;
    return iced;
}

// The number of kept fine cells whose ice thickness differs between two passes.
int changedice(const RegionalFields& a, const RegionalFields& b)
{
    int changed = 0;
    for (std::size_t cell = 0; cell < a.icethickness_m.size(); cell++)
        changed += a.icethickness_m[cell] != b.icethickness_m[cell] ? 1 : 0;
    return changed;
}
}

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

TEST_F(Landscape, RegionalIceAboveTheElaIsInertWithoutTheFields)
{
    // The ELA ice on but no ELA field (inactive FineGlacial): nothing is iced.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);
    const std::size_t count = domaincount(window);
    const RegionalParameters parameters = iceparameters(0.3, 1000.0, 0.0);
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_FALSE(glacial.active);
    const std::vector<uint8_t> sea(count, 0);
    const std::vector<float> elevation(count, 9000.0f);
    std::vector<float> ice(count, 0.0f);
    landscape::iceabovetheela(inputs, window, parameters, glacial, sea, elevation, ice);
    EXPECT_EQ(icedcells(ice), 0);
}

TEST_F(Landscape, RegionalIceThresholdAndThicknessLaw)
{
    // No edge noise, fraction 0.1, cap 150 m: on a hand-built surface at the present ELA plus
    // a height that cycles through -10, -0.5, 0, 5, 100, 1000 and 3000 m, the thickness is
    // 0 below the ELA, the 1 m floor at and just above it (0.1 x 5 = 0.5), 10 and 100 m in
    // the law's range and the cap 150 m above 1,500 m.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = domaincount(window);
    const RegionalParameters parameters = iceparameters(0.1, 150.0, 0.0);
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    const double heights[] = {-10.0, -0.5, 0.0, 5.0, 100.0, 1000.0, 3000.0};
    const float expected[] = {0.0f, 0.0f, 1.0f, 1.0f, 10.0f, 100.0f, 150.0f};
    std::vector<float> elevation(count);
    for (std::size_t cell = 0; cell < count; cell++)
        elevation[cell] = static_cast<float>(static_cast<double>(glacial.presentela_m[cell]) + heights[cell % 7]);
    const std::vector<uint8_t> sea(count, 0);
    std::vector<float> ice(count, 0.0f);
    landscape::iceabovetheela(inputs, window, parameters, glacial, sea, elevation, ice);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        SCOPED_TRACE("cell " + std::to_string(cell) + ", height " + std::to_string(heights[cell % 7]));
        // The float surface rounds the height by at most an ulp of the ELA (2.4e-4 m below 4,096 m).
        if (heights[cell % 7] == 0.0 || heights[cell % 7] < 0.0)
            ASSERT_EQ(ice[cell], expected[cell % 7]);
        else
            ASSERT_NEAR(ice[cell], expected[cell % 7], 1.0e-3);
    }
}

TEST_F(Landscape, RegionalIceKeepsExistingIceAndSkipsTheSea)
{
    // Existing ice is never removed or thinned: a cell below the ELA keeps its 40 m; a cell
    // 1,000 m above it (law: 300 m) keeps its 500 m and one with 50 m is raised to 300 m. A
    // sea cell above the ELA takes no ice and keeps the ice it has (the calving front of 5).
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = domaincount(window);
    const RegionalParameters parameters = iceparameters(0.3, 1000.0, 0.0);
    const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
    ASSERT_TRUE(glacial.active);
    std::vector<float> elevation = abovetheela(glacial, 1000.0);
    std::vector<uint8_t> sea(count, 0);
    std::vector<float> ice(count, 0.0f);
    elevation[0] = static_cast<float>(static_cast<double>(glacial.presentela_m[0]) - 200.0);
    ice[0] = 40.0f;
    ice[1] = 500.0f;
    ice[2] = 50.0f;
    sea[3] = 1;
    sea[4] = 1;
    ice[4] = 20.0f;
    landscape::iceabovetheela(inputs, window, parameters, glacial, sea, elevation, ice);
    EXPECT_EQ(ice[0], 40.0f);
    EXPECT_EQ(ice[1], 500.0f);
    EXPECT_NEAR(ice[2], 300.0f, 1.0e-3);
    EXPECT_EQ(ice[3], 0.0f);
    EXPECT_EQ(ice[4], 20.0f);
    EXPECT_NEAR(ice[5], 300.0f, 1.0e-3);
}

TEST_F(Landscape, RegionalIceEdgeNoise)
{
    // A surface exactly at the present ELA: without the edge noise every cell takes the 1 m
    // floor; with 100 m of noise only the cells where the noise is at most 0 do, a share
    // between a quarter and three quarters; 1,000 m above the ELA (beyond the noise's reach
    // at k 4) every cell is iced either way, 1,000 m below none.
    const RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = domaincount(window);
    const std::vector<uint8_t> sea(count, 0);
    for (double edgenoise : {0.0, 100.0})
    {
        SCOPED_TRACE("edge noise " + std::to_string(edgenoise));
        const RegionalParameters parameters = iceparameters(0.3, 1000.0, edgenoise);
        const landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, parameters);
        ASSERT_TRUE(glacial.active);
        std::vector<float> ice(count, 0.0f);
        landscape::iceabovetheela(inputs, window, parameters, glacial, sea, abovetheela(glacial, 0.0), ice);
        const int iced = icedcells(ice);
        if (edgenoise == 0.0)
            EXPECT_EQ(iced, static_cast<int>(count));
        else
        {
            EXPECT_GT(iced, static_cast<int>(count / 4));
            EXPECT_LT(iced, static_cast<int>(3 * count / 4));
        }
        for (float thickness : ice)
            ASSERT_TRUE(thickness == 0.0f || thickness == 1.0f);
        std::vector<float> high(count, 0.0f), low(count, 0.0f);
        landscape::iceabovetheela(inputs, window, parameters, glacial, sea, abovetheela(glacial, 1000.0), high);
        landscape::iceabovetheela(inputs, window, parameters, glacial, sea, abovetheela(glacial, -1000.0), low);
        EXPECT_EQ(icedcells(high), static_cast<int>(count));
        EXPECT_EQ(icedcells(low), 0);
    }
}

TEST_F(Landscape, RegionalOverlapIce)
{
    // The ELA ice on (fraction 0.3, edge noise 100 m) at the smallest halo the pass allows:
    // the overlaps stay bit-identical (compareoverlap asserts it), and the rule adds ice. On
    // the synthetic world the global ice (above 2,000 m) covers the fine cells above the
    // present ELA (1,800 to 2,400 m) and their neighbours, so the rule thickens that ice; on
    // the same world without global ice every iced cell is the rule's.
    const RegionalParameters parameters = iceparameters(0.3, 1000.0, 100.0);
    for (bool globalice : {true, false})
    {
        SCOPED_TRACE(globalice ? "with the global ice" : "without the global ice");
        RegionalGlobalInputs columns = withglacialfields(syntheticworld(32));
        RegionalGlobalInputs rows = withglacialfields(syntheticworld(64));
        if (!globalice)
        {
            std::fill(columns.icethickness_m.begin(), columns.icethickness_m.end(), 0.0f);
            std::fill(rows.icethickness_m.begin(), rows.icethickness_m.end(), 0.0f);
        }
        for (int scale : {4, 5})
        {
            SCOPED_TRACE("k = " + std::to_string(scale));
            const OverlapCounts counts = compareoverlap(columns, windowa(scale), windowb(scale), "ELA ice A/B", parameters);
            EXPECT_EQ(counts.shared, 4 * 5 * scale * scale);
            EXPECT_GT(counts.ice, 0);
            const RegionalFields on = landscape::runregionalpass(columns, windowb(scale), parameters);
            const RegionalFields off = landscape::runregionalpass(columns, windowb(scale), RegionalParameters{});
            EXPECT_GT(changedice(on, off), 0);
            if (!globalice)
            {
                EXPECT_EQ(icedcells(off.icethickness_m), 0);
                EXPECT_GT(icedcells(on.icethickness_m), 0);
            }
            EXPECT_EQ(on.elevation_m, off.elevation_m);  // the rule places ice, it shapes nothing
        }
        const OverlapCounts counts = compareoverlap(rows, windowc(4), windowd(4), "ELA ice C/D", parameters);
        EXPECT_GT(counts.sea, 0);
        EXPECT_GT(counts.ice, 0);
        const RegionalFields on = landscape::runregionalpass(rows, windowc(4), parameters);
        const RegionalFields off = landscape::runregionalpass(rows, windowc(4), RegionalParameters{});
        EXPECT_GT(changedice(on, off), 0);
    }
}
