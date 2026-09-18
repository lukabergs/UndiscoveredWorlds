#include "gtest/gtest.h"
#include "tectonic_boundary_graph.hpp"
#include "tectonic_provenance.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {
using namespace platec;
using contract::CrustClass;
using contract::GeologicRegime;

struct BoundaryFixture {
    static constexpr uint32_t width = 8, height = 8;
    std::vector<float> elevation = std::vector<float>(64, 0.5f);
    std::vector<uint32_t> owner = std::vector<uint32_t>(64), zero = owner;
    std::vector<uint8_t> convergence = std::vector<uint8_t>(64), divergence = convergence,
                         shear = convergence, regime = convergence;
    std::vector<uint8_t> crust = std::vector<uint8_t>(64, uint8_t(CrustClass::Oceanic));
    std::array<contract::PlateKinematics, 2> plates{};

    BoundaryFixture() {
        for (uint32_t y = 0; y < height; ++y)
            for (uint32_t x = 0; x < width; ++x) owner[y * width + x] = x < 4 ? 0 : 1;
        velocity(1.0f, 0.0f);
    }
    void velocity(float speed, float common_y) {
        for (size_t i = 0; i < plates.size(); ++i) {
            const float x = i == 0 ? speed : -speed;
            const float length = std::hypot(x, common_y);
            plates[i].linear_velocity_x = x;
            plates[i].linear_velocity_y = common_y;
            // Retain opposing headings even when stopped, as movement may do.
            plates[i].unit_velocity_x = length > 0 ? x / length : (i == 0 ? 1.0f : -1.0f);
            plates[i].unit_velocity_y = length > 0 ? common_y / length : 0.0f;
        }
    }
    void compute() {
        provenance::Inputs inputs{width, height, elevation.data(), owner.data(),
                                  plates.data(), 2, crust.data()};
        provenance::compute_maps(inputs, convergence.data(), divergence.data(), shear.data(), regime.data());
    }
    boundary_graph::Outputs graph() const {
        boundary_graph::Inputs inputs{width, height, owner.data(), convergence.data(),
            divergence.data(), shear.data(), regime.data(), crust.data(), zero.data(),
            zero.data(), nullptr, 0, 1, 1.0};
        return boundary_graph::build(inputs);
    }
};
}

TEST(BoundaryMotion, CommonTranslationDoesNotChangeRelativeSignals) {
    BoundaryFixture f;
    f.compute();
    const auto convergence = f.convergence, divergence = f.divergence, shear = f.shear;
    ASSERT_EQ(convergence[27], 100);
    f.velocity(1, 10);
    f.compute();
    EXPECT_EQ(f.convergence, convergence);
    EXPECT_EQ(f.divergence, divergence);
    EXPECT_EQ(f.shear, shear);
}

TEST(BoundaryMotion, SpeedAndStoppedPlatesAffectSignals) {
    BoundaryFixture f;
    f.velocity(0.5f, 0);
    f.compute();
    EXPECT_GT(f.convergence[27], 0);
    EXPECT_LT(f.convergence[27], 100);
    f.velocity(0, 0);
    f.compute();
    EXPECT_EQ(*std::max_element(f.convergence.begin(), f.convergence.end()), 0);
    EXPECT_EQ(*std::max_element(f.divergence.begin(), f.divergence.end()), 0);
    EXPECT_EQ(*std::max_element(f.shear.begin(), f.shear.end()), 0);
}

TEST(BoundaryMotion, RotationContributesAtTheSharedBoundaryPoint) {
    BoundaryFixture f;
    f.velocity(0, 0);
    f.plates[0].mass_center_x = 3.5f;
    f.plates[0].mass_center_y = 4.0f;
    f.plates[0].angular_velocity = 2.0f;
    f.compute();
    EXPECT_EQ(f.convergence[27], 100);
    EXPECT_EQ(f.convergence[28], 100);
}

TEST(BoundaryMotion, BothPeriodicSeamsCarryMotion) {
    BoundaryFixture f;
    f.compute();
    EXPECT_EQ(f.divergence[24], 100);
    EXPECT_EQ(f.divergence[31], 100);
    for (uint32_t y = 0; y < f.height; ++y)
        for (uint32_t x = 0; x < f.width; ++x) f.owner[y * f.width + x] = y < 4 ? 0 : 1;
    for (auto& plate : f.plates) {
        plate.linear_velocity_y = plate.linear_velocity_x;
        plate.linear_velocity_x = 0;
    }
    f.compute();
    EXPECT_EQ(f.divergence[3], 100);
    EXPECT_EQ(f.divergence[59], 100);
}

TEST(BoundaryMotion, OceanicConvergenceIsNeverContinentalCollision) {
    BoundaryFixture f;
    f.compute();
    EXPECT_EQ(f.regime[27], uint8_t(GeologicRegime::ConvergentArc));
    // Even stale raw hints must not override the actual crust pair in the graph.
    std::fill(f.regime.begin(), f.regime.end(), uint8_t(GeologicRegime::ContinentCollision));
    const auto graph = f.graph();
    ASSERT_FALSE(graph.segments.empty());
    for (const auto& segment : graph.segments)
        EXPECT_NE(segment.geologic_regime, GeologicRegime::ContinentCollision);
    for (auto regime : graph.boundary_regime)
        EXPECT_NE(regime, uint8_t(GeologicRegime::ContinentCollision));
}

TEST(BoundaryMotion, SubmergedContinentsCollideAndMixedCrustPreservesSides) {
    BoundaryFixture f;
    std::fill(f.crust.begin(), f.crust.end(), uint8_t(CrustClass::Continental));
    f.compute();
    EXPECT_EQ(f.regime[27], uint8_t(GeologicRegime::ContinentCollision));
    for (size_t i = 0; i < f.crust.size(); ++i)
        if (f.owner[i] == 1) f.crust[i] = uint8_t(CrustClass::Oceanic);
    f.compute();
    const auto graph = f.graph();
    EXPECT_EQ(graph.boundary_regime[27], uint8_t(GeologicRegime::ConvergentArc));
    EXPECT_EQ(graph.boundary_regime[28], uint8_t(GeologicRegime::TrenchAdjacent));
}

TEST(BoundaryMotion, NormalMotionIsIndependentOfBoundaryOrientation) {
    BoundaryFixture f;
    std::fill(f.convergence.begin(), f.convergence.end(), 100);
    for (int orientation = 0; orientation < 2; ++orientation) {
        const auto graph = f.graph();
        ASSERT_FALSE(graph.segments.empty());
        for (const auto& segment : graph.segments) EXPECT_FLOAT_EQ(segment.average_normal_motion, 1);
        for (uint32_t y = 0; y < f.height; ++y)
            for (uint32_t x = 0; x < f.width; ++x) f.owner[y * f.width + x] = y < 4 ? 0 : 1;
    }
}

TEST(BoundaryMotion, SubmergedContinentalRiftsAreNotMidOceanRidges) {
    BoundaryFixture f;
    std::fill(f.crust.begin(), f.crust.end(), uint8_t(CrustClass::Continental));
    f.velocity(-1, 0);
    f.compute();
    EXPECT_EQ(f.regime[27], uint8_t(GeologicRegime::DivergentRift));
    const auto graph = f.graph();
    EXPECT_EQ(graph.boundary_regime[27], uint8_t(GeologicRegime::DivergentRift));
}
