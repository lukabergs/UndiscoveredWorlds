// Tests of regional_coast.hpp (erosion plan Phase R2: the R6 contract; R7 adds the island
// sink's tests and the RegionalOverlapIslands* variants with its features on, on a world wider
// than the land-share rule's widened box), suite Landscape: TEST_F(Landscape,
// RegionalCoast<What>) on the synthetic world of regional_test_world.hpp.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_coast.hpp"

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
