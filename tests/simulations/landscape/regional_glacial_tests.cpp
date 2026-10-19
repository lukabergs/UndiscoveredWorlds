// Tests of landscape/regional_glacial.hpp (erosion plan Phase E5, step E5.0: the fine ELA
// and glaciation fields shared by the glacial forms, and the ice catchment of the global
// network), suite Landscape: TEST_F(Landscape, RegionalGlacial<What>) on hand-built global
// chains and the synthetic world of regional_test_world.hpp with withglacialfields.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_glacial.hpp"

namespace
{
RegionalParameters troughson()
{
    RegionalParameters parameters;
    parameters.troughs.width_km = 10.0;
    return parameters;
}
}

TEST_F(Landscape, RegionalGlacialInactiveAtTheDefaultsAndWithoutTheField)
{
    // Nothing is on at the defaults, so the struct is inactive even with the fields; with a
    // feature on but no ELA field it is inactive too; the switch names the three features.
    const RegionalGlobalInputs plain = syntheticworld(32);
    const RegionalGlobalInputs glacial = withglacialfields(plain);
    const RegionalWindow window = windowa(4);
    EXPECT_FALSE(landscape::glacialformson(RegionalParameters{}));
    EXPECT_FALSE(landscape::regionalglacialfields(glacial, window, RegionalParameters{}).active);
    EXPECT_TRUE(landscape::regionalglacialfields(glacial, window, RegionalParameters{}).runela_m.empty());
    EXPECT_FALSE(landscape::regionalglacialfields(plain, window, troughson()).active);
    RegionalParameters cirques;
    cirques.cirques.radius_km = 2.0;
    RegionalParameters ice;
    ice.ice.thicknessfraction = 0.3;
    EXPECT_TRUE(landscape::glacialformson(troughson()));
    EXPECT_TRUE(landscape::glacialformson(cirques));
    EXPECT_TRUE(landscape::glacialformson(ice));
}

TEST_F(Landscape, RegionalGlacialFieldsAreTheUpsampledGlobalOnes)
{
    // Active with the ELA field and a feature on: domain-sized, the present ELA the run's
    // plus the offset on every cell, the glaciation within 0..1 and the ELA of a constant
    // field its constant; a block's mean ELA is near its global value (bicubic through the
    // centres) and a glaciation above 1 is clamped.
    RegionalGlobalInputs inputs = withglacialfields(syntheticworld(32));
    const RegionalWindow window = windowa(4);
    const std::size_t count = static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    landscape::FineGlacial glacial = landscape::regionalglacialfields(inputs, window, troughson());
    ASSERT_TRUE(glacial.active);
    ASSERT_EQ(glacial.runela_m.size(), count);
    ASSERT_EQ(glacial.presentela_m.size(), count);
    ASSERT_EQ(glacial.glaciation.size(), count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        EXPECT_NEAR(glacial.presentela_m[cell] - glacial.runela_m[cell], 500.0f, 0.01f);
        EXPECT_GE(glacial.glaciation[cell], 0.0f);
        EXPECT_LE(glacial.glaciation[cell], 1.0f);
        EXPECT_GT(glacial.runela_m[cell], 1200.0f);
        EXPECT_LT(glacial.runela_m[cell], 2000.0f);
    }
    // The kept window's first block: its fine mean against the global value.
    const int k = window.scale, domainwidth = window.domainfinewidth();
    const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
    double mean = 0.0;
    for (int y = y0; y < y0 + k; y++)
        for (int x = x0; x < x0 + k; x++)
            mean += glacial.runela_m[at(x, y, domainwidth)];
    mean /= static_cast<double>(k * k);
    const float global = inputs.equilibriumline_m[window.globalcell(x0, y0)];
    EXPECT_NEAR(mean, global, 30.0);

    std::fill(inputs.equilibriumline_m.begin(), inputs.equilibriumline_m.end(), 1500.0f);
    std::fill(inputs.glaciation.begin(), inputs.glaciation.end(), 3.0f);
    inputs.presentelaoffset_m = 0.0;
    glacial = landscape::regionalglacialfields(inputs, window, troughson());
    for (std::size_t cell = 0; cell < count; cell++)
    {
        EXPECT_NEAR(glacial.runela_m[cell], 1500.0f, 1e-3f);
        EXPECT_EQ(glacial.presentela_m[cell], glacial.runela_m[cell]);
        EXPECT_EQ(glacial.glaciation[cell], 1.0f);
    }
    // Without a glaciation field the fine glaciation is 0 everywhere, the ELA still there.
    inputs.glaciation.clear();
    glacial = landscape::regionalglacialfields(inputs, window, troughson());
    ASSERT_TRUE(glacial.active);
    ASSERT_EQ(glacial.glaciation.size(), count);
    for (std::size_t cell = 0; cell < count; cell++)
        EXPECT_EQ(glacial.glaciation[cell], 0.0f);
}

TEST_F(Landscape, RegionalGlacialCatchmentSumsUpstreamGlaciatedArea)
{
    // An 8 x 4 world: row 1 a chain 1 -> 2 -> 3 -> 4 -> 5 flowing east, cell (6, 1) sea
    // receiving (5, 1); row 2 cell (3, 2) flows north into (3, 1); glaciation 1 on (1, 1),
    // 0.5 on (3, 2), 0 elsewhere. The catchment of (3, 1) is area(1,1) + 0.5 area(3,2),
    // (5, 1) the same, (2, 1) area(1,1), (4, 2) and the sea 0; a world without the
    // glaciation field gives an empty vector.
    RegionalGlobalInputs inputs;
    inputs.width = 8;
    inputs.height = 4;
    inputs.grid = mesh::makelatlongrid(8, radius_m);
    inputs.radius_m = radius_m;
    inputs.sealevel = 0.0f;
    inputs.seed = 1;
    const std::size_t size = 32;
    inputs.elevation_m.assign(size, 100.0f);
    landscape::GlobalNetwork& network = inputs.network;
    network.width = 8;
    network.height = 4;
    network.sea.assign(size, 0);
    network.direction.assign(size, 0);
    network.filledelevation_m.assign(size, 100.0f);
    network.drainagearea_km2.assign(size, 1.0f);
    network.discharge_m3s.assign(size, 1.0f);
    network.januaryflow.assign(size, 0);
    network.julyflow.assign(size, 0);
    for (int x = 1; x <= 5; x++)
        network.direction[at(x, 1, 8)] = 3;  // east
    network.sea[at(6, 1, 8)] = 1;
    network.direction[at(3, 2, 8)] = 1;      // north
    EXPECT_TRUE(landscape::glacialcatchment_km2(inputs).empty());
    inputs.glaciation.assign(size, 0.0f);
    inputs.glaciation[at(1, 1, 8)] = 1.0f;
    inputs.glaciation[at(3, 2, 8)] = 0.5f;
    const std::vector<double> catchment = landscape::glacialcatchment_km2(inputs);
    ASSERT_EQ(catchment.size(), size);
    const landscape::SphereMetric metric(8, 4, radius_m);
    const double area11 = metric.area_m2(at(1, 1, 8)) / 1.0e6, area32 = metric.area_m2(at(3, 2, 8)) / 1.0e6;
    EXPECT_GT(area11, 1.0e5);  // a 45-degree cell of an Earth-sized sphere: millions of km2
    EXPECT_NEAR(catchment[at(1, 1, 8)], area11, 1e-6 * area11);
    EXPECT_NEAR(catchment[at(2, 1, 8)], area11, 1e-6 * area11);
    EXPECT_NEAR(catchment[at(3, 2, 8)], 0.5 * area32, 1e-6 * area32);
    EXPECT_NEAR(catchment[at(3, 1, 8)], area11 + 0.5 * area32, 1e-6 * area11);
    EXPECT_NEAR(catchment[at(5, 1, 8)], area11 + 0.5 * area32, 1e-6 * area11);
    EXPECT_EQ(catchment[at(4, 2, 8)], 0.0);
    EXPECT_EQ(catchment[at(6, 1, 8)], 0.0);
    EXPECT_EQ(catchment[at(0, 0, 8)], 0.0);
}

TEST_F(Landscape, RegionalGlacialFieldsLeaveTheDefaultPassUnchanged)
{
    // The Phase R2 parity, structural: the pass at the defaults on the world with the
    // glacial fields equals the pass on the world without them, bit for bit (nothing is
    // read while every glacial feature is off).
    const RegionalGlobalInputs plain = syntheticworld(32);
    const RegionalGlobalInputs glacial = withglacialfields(plain);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalFields a = landscape::runregionalpass(plain, windowb(scale), RegionalParameters{});
        const RegionalFields b = landscape::runregionalpass(glacial, windowb(scale), RegionalParameters{});
        ASSERT_EQ(a.elevation_m.size(), b.elevation_m.size());
        int changed = 0;
        for (std::size_t cell = 0; cell < a.elevation_m.size(); cell++)
            changed += (a.elevation_m[cell] != b.elevation_m[cell] || a.icethickness_m[cell] != b.icethickness_m[cell] ||
                        a.lakesurface_m[cell] != b.lakesurface_m[cell] || a.sea[cell] != b.sea[cell])
                           ? 1
                           : 0;
        EXPECT_EQ(changed, 0);
    }
}
