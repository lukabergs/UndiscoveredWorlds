#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

TEST(TectonicTransport, RigidTranslationRetainsResolvedReliefAndConservesMass) {
    constexpr uint32_t w = 64, h = 32;
    constexpr double pi = 3.141592653589793;
    for (const double wavelength : {4., 8., 16.}) {
        std::vector<double> initial(w * h);
        for (uint32_t y = 0; y < h; ++y) for (uint32_t x = 0; x < w; ++x)
            initial[y * w + x] = 1000 + 200 * std::sin(2 * pi * x / wavelength);
        auto field = initial;
        const std::vector<tectonics::transport::Displacement> v(w * h, {.16, -.12});
        for (int t = 0; t < 200; ++t) field = tectonics::transport::remap(w, h, field, v);
        double error = 0, amplitude = 0;
        for (uint32_t y = 0; y < h; ++y) for (uint32_t x = 0; x < w; ++x) {
            const double expected = initial[y * w + (x + w - 32) % w];
            error += std::pow(field[y * w + x] - expected, 2);
            amplitude += (field[y * w + x] - 1000) * std::sin(2 * pi * (static_cast<double>(x) - 32) / wavelength);
        }
        const double retained = amplitude / (w * h * 100);
        RecordProperty("retained_wavelength_" + std::to_string(static_cast<int>(wavelength)), retained);
        // Four-cell wavelengths quantify the unresolved end; they are not a
        // claim of preserved detail. No noise is injected to increase amplitude.
        if (wavelength >= 8) EXPECT_GT(retained, wavelength == 8 ? .60 : .90) << "wavelength=" << wavelength;
        EXPECT_LT(std::sqrt(error / (w * h)), wavelength == 4 ? 142 : wavelength == 8 ? 65 : 20);
        EXPECT_NEAR(std::accumulate(field.begin(), field.end(), 0.), w * h * 1000., 1e-6);
        EXPECT_GE(*std::min_element(field.begin(), field.end()), 800 - 1e-9);
        EXPECT_LE(*std::max_element(field.begin(), field.end()), 1200 + 1e-9);
    }
}

TEST(TectonicTransport, DiagonalVariableMotionPreservesPositiveMassAcrossBothSeams) {
    std::vector<double> field(128);
    std::vector<tectonics::transport::Displacement> v(128);
    for (size_t i = 0; i < field.size(); ++i) {
        field[i] = i % 5 == 0 ? 0 : static_cast<double>((i * 17) % 31);
        v[i] = {i % 2 ? .93 : -.89, i % 3 ? -.83 : .91};
    }
    const double mass = std::accumulate(field.begin(), field.end(), 0.);
    for (int i = 0; i < 100; ++i) field = tectonics::transport::remap(16, 8, field, v);
    EXPECT_GE(*std::min_element(field.begin(), field.end()), 0);
    EXPECT_NEAR(std::accumulate(field.begin(), field.end(), 0.), mass, mass * 1e-12);
    EXPECT_THROW(tectonics::transport::remap(16, 8, field, {}), std::runtime_error);
}

TEST(TectonicTransport, AgeUsesRockFluxesAndStaysWithinDonorAges) {
    std::vector<double> rock(128), moment(128);
    std::vector<tectonics::transport::Displacement> v(128);
    for (size_t i = 0; i < rock.size(); ++i) {
        rock[i] = i % 3 ? 100 + static_cast<double>((i * 13) % 27000) : 1e-20;
        moment[i] = rock[i] * (i % 5 ? 100 : 300);
        v[i] = {i % 2 ? .13 : -.11, i % 3 ? -.12 : .17};
    }
    const double initial_rock = std::accumulate(rock.begin(), rock.end(), 0.);
    const double initial_moment = std::accumulate(moment.begin(), moment.end(), 0.);
    for (int t = 0; t < 200; ++t) {
        std::vector<double> next_moment;
        rock = tectonics::transport::remap_with_moment(16, 8, rock, moment, v, next_moment);
        moment = std::move(next_moment);
        for (size_t i = 0; i < rock.size(); ++i) if (rock[i] > 0) {
            ASSERT_GE(moment[i] / rock[i], 100 - 1e-9) << t << ", cell " << i;
            ASSERT_LE(moment[i] / rock[i], 300 + 1e-9) << t << ", cell " << i;
        }
    }
    EXPECT_NEAR(std::accumulate(rock.begin(), rock.end(), 0.), initial_rock, initial_rock * 1e-12);
    EXPECT_NEAR(std::accumulate(moment.begin(), moment.end(), 0.), initial_moment, initial_moment * 1e-12);
}

TEST(TectonicTopology, FractionalRigidTranslationDoesNotLeaveTrailingFragments) {
    constexpr uint32_t w = 64, h = 32;
    std::vector<uint32_t> initial(w * h);
    for (uint32_t y = 0; y < h; ++y) for (uint32_t x = 0; x < w; ++x)
        initial[y * w + x] = (x >= 10 && x < 38 + (y / 8) % 2) ? 1U : 0U;
    tectonics::topology::Model model;
    model.initialize(w, h, 2, initial);
    std::vector<uint32_t> result;
    for (int i = 0; i < 80; ++i) result = model.advance({{.25, -.25}, {.25, -.25}});
    size_t different = 0;
    for (uint32_t y = 0; y < h; ++y) for (uint32_t x = 0; x < w; ++x)
        different += result[y * w + x] != initial[((y + 20) % h) * w + (x + w - 20) % w];
    EXPECT_LE(different, 16U); // Only unresolved one-cell corners may change.
    EXPECT_EQ(model.diagnostics().reassigned_cells, 0U);
}

TEST(TectonicTopology, TransformMotionLeavesStraightBoundaryConnected) {
    std::vector<uint32_t> initial(32 * 16);
    for (size_t i = 0; i < initial.size(); ++i) initial[i] = i % 32 < 16 ? 0U : 1U;
    tectonics::topology::Model model;
    model.initialize(32, 16, 2, initial);
    for (int i = 0; i < 200; ++i) EXPECT_EQ(model.advance({{0, .2}, {0, -.2}}), initial);
    EXPECT_EQ(model.diagnostics().reassigned_cells, 0U);
}

TEST(TectonicTopology, DetachedRasterIslandIsRepairedButPeriodicRegionIsConnected) {
    std::vector<uint32_t> initial(32 * 16);
    for (size_t i = 0; i < initial.size(); ++i) initial[i] = i % 32 >= 8 && i % 32 < 24 ? 1U : 0U;
    initial[8 * 32 + 16] = 0;
    tectonics::topology::Model model;
    model.initialize(32, 16, 2, initial);
    auto repaired = model.advance({{}, {}});
    initial[8 * 32 + 16] = 1;
    EXPECT_EQ(repaired, initial);
    EXPECT_EQ(model.diagnostics().max_detached_components, 1U);
    EXPECT_EQ(model.diagnostics().reassigned_cells, 1U);
    EXPECT_EQ(model.advance({{}, {}}), initial);
    EXPECT_EQ(model.diagnostics().reassigned_cells, 1U);
}
