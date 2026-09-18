#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "tectonic_junction_coupling.hpp"
#include "tectonic_junctions.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

using platec::contract::BoundarySegment;
using platec::contract::BoundaryType;
using platec::contract::JunctionArmType;
using platec::contract::JunctionStability;
using platec::contract::PlateKinematics;
using platec::contract::TectonicJunction;

struct SyntheticJunctionInput {
    uint32_t width = 8;
    uint32_t height = 8;
    std::vector<uint32_t> plate_id;
    std::vector<uint32_t> boundary_segment_id;
    std::vector<BoundarySegment> segments;
    std::vector<PlateKinematics> plates;
};

struct TestVelocity {
    float x = 0.0f;
    float y = 0.0f;
};

size_t index_of(uint32_t x, uint32_t y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
}

float wrapped_delta_for_test(float value, float center, uint32_t period) {
    float delta = value - center;
    const float period_f = static_cast<float>(period);
    const float half_period = period_f * 0.5f;
    if (delta > half_period) {
        delta -= period_f;
    } else if (delta < -half_period) {
        delta += period_f;
    }
    return delta;
}

void set_plate_cell(SyntheticJunctionInput& input, uint32_t x, uint32_t y,
                    uint32_t plate_id) {
    input.plate_id[index_of(x, y, input.width)] = plate_id;
}

void mark_segment_cell(SyntheticJunctionInput& input, uint32_t x, uint32_t y,
                       uint32_t segment_id) {
    input.boundary_segment_id[index_of(x, y, input.width)] = segment_id;
}

void set_plate_velocity(SyntheticJunctionInput& input, uint32_t plate_id,
                        TestVelocity velocity) {
    if (input.plates.size() <= plate_id) {
        input.plates.resize(static_cast<size_t>(plate_id) + 1U);
    }
    input.plates[plate_id].linear_velocity_x = velocity.x;
    input.plates[plate_id].linear_velocity_y = velocity.y;
}

BoundarySegment segment(uint32_t id, uint32_t left_plate, uint32_t right_plate,
                        BoundaryType type) {
    BoundarySegment result;
    result.id = id;
    result.left_plate_id = left_plate;
    result.right_plate_id = right_plate;
    result.cell_count = 1;
    result.length_cells = 1.0f;
    result.boundary_type = type;
    return result;
}

SyntheticJunctionInput make_explicit_arm_input(
    BoundaryType ab, BoundaryType ac, BoundaryType bc,
    const std::array<TestVelocity, 3>& velocities) {
    SyntheticJunctionInput input;
    input.width = 12;
    input.height = 12;
    input.plate_id.assign(static_cast<size_t>(input.width) * input.height, 0U);
    input.boundary_segment_id.assign(input.plate_id.size(), platec::contract::kNoBoundaryId);

    for (uint32_t y = 0; y < input.height; ++y) {
        for (uint32_t x = 0; x < input.width; ++x) {
            const uint32_t plate_id = x <= 5U ? 0U : (y <= 5U ? 1U : 2U);
            set_plate_cell(input, x, y, plate_id);
        }
    }

    input.segments = {
        segment(10U, 0U, 1U, ab),
        segment(20U, 0U, 2U, ac),
        segment(30U, 1U, 2U, bc),
    };

    mark_segment_cell(input, 5U, 3U, 10U);
    mark_segment_cell(input, 5U, 4U, 10U);
    mark_segment_cell(input, 5U, 5U, 10U);
    mark_segment_cell(input, 3U, 6U, 20U);
    mark_segment_cell(input, 4U, 6U, 20U);
    mark_segment_cell(input, 5U, 6U, 20U);
    mark_segment_cell(input, 6U, 6U, 30U);
    mark_segment_cell(input, 7U, 7U, 30U);
    mark_segment_cell(input, 8U, 8U, 30U);

    for (size_t plate_id = 0; plate_id < velocities.size(); ++plate_id) {
        set_plate_velocity(input, static_cast<uint32_t>(plate_id), velocities[plate_id]);
    }

    return input;
}

SyntheticJunctionInput make_wrapped_candidate_input() {
    SyntheticJunctionInput input;
    input.width = 12;
    input.height = 12;
    input.plate_id.assign(static_cast<size_t>(input.width) * input.height, 0U);
    input.boundary_segment_id.assign(input.plate_id.size(), platec::contract::kNoBoundaryId);
    input.segments = {
        segment(10U, 0U, 1U, BoundaryType::Divergent),
        segment(20U, 0U, 2U, BoundaryType::Divergent),
        segment(30U, 1U, 2U, BoundaryType::Divergent),
    };

    set_plate_cell(input, 11U, 5U, 0U);
    set_plate_cell(input, 0U, 5U, 1U);
    set_plate_cell(input, 1U, 5U, 2U);
    set_plate_cell(input, 2U, 5U, 1U);
    set_plate_cell(input, 11U, 6U, 2U);
    set_plate_cell(input, 0U, 6U, 0U);
    set_plate_cell(input, 1U, 6U, 1U);
    set_plate_cell(input, 2U, 6U, 1U);

    mark_segment_cell(input, 11U, 3U, 10U);
    mark_segment_cell(input, 11U, 4U, 10U);
    mark_segment_cell(input, 11U, 5U, 10U);
    mark_segment_cell(input, 10U, 6U, 20U);
    mark_segment_cell(input, 11U, 6U, 20U);
    mark_segment_cell(input, 0U, 6U, 20U);
    mark_segment_cell(input, 0U, 5U, 30U);
    mark_segment_cell(input, 1U, 6U, 30U);
    mark_segment_cell(input, 2U, 7U, 30U);

    return input;
}

SyntheticJunctionInput make_clustered_candidate_input() {
    SyntheticJunctionInput input;
    input.width = 12;
    input.height = 12;
    input.plate_id.assign(static_cast<size_t>(input.width) * input.height, 0U);
    input.boundary_segment_id.assign(input.plate_id.size(), platec::contract::kNoBoundaryId);
    input.segments = {
        segment(10U, 0U, 1U, BoundaryType::Divergent),
        segment(20U, 0U, 2U, BoundaryType::Divergent),
        segment(30U, 1U, 2U, BoundaryType::Divergent),
    };

    set_plate_cell(input, 4U, 5U, 0U);
    set_plate_cell(input, 5U, 5U, 1U);
    set_plate_cell(input, 6U, 5U, 2U);
    set_plate_cell(input, 7U, 5U, 1U);
    set_plate_cell(input, 4U, 6U, 2U);
    set_plate_cell(input, 5U, 6U, 0U);
    set_plate_cell(input, 6U, 6U, 1U);
    set_plate_cell(input, 7U, 6U, 1U);

    mark_segment_cell(input, 4U, 5U, 10U);
    mark_segment_cell(input, 5U, 5U, 10U);
    mark_segment_cell(input, 4U, 6U, 20U);
    mark_segment_cell(input, 5U, 6U, 20U);
    mark_segment_cell(input, 6U, 5U, 30U);
    mark_segment_cell(input, 6U, 6U, 30U);

    return input;
}

SyntheticJunctionInput make_three_plate_input(BoundaryType ab, BoundaryType ac,
                                              BoundaryType bc, bool include_bc = true) {
    SyntheticJunctionInput input;
    input.plate_id.assign(static_cast<size_t>(input.width) * input.height, 0U);
    input.boundary_segment_id.assign(input.plate_id.size(), platec::contract::kNoBoundaryId);

    input.plate_id[index_of(4, 3, input.width)] = 1U;
    input.plate_id[index_of(3, 4, input.width)] = 2U;

    input.boundary_segment_id[index_of(4, 3, input.width)] = 10U;
    input.boundary_segment_id[index_of(3, 4, input.width)] = 20U;
    input.segments.push_back(segment(10U, 0U, 1U, ab));
    input.segments.push_back(segment(20U, 0U, 2U, ac));
    if (include_bc) {
        input.boundary_segment_id[index_of(4, 4, input.width)] = 30U;
        input.segments.push_back(segment(30U, 1U, 2U, bc));
    }

    return input;
}

SyntheticJunctionInput make_oriented_three_plate_input() {
    SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Convergent, BoundaryType::Transform);

    const auto mark = [&](uint32_t x, uint32_t y, uint32_t segment_id) {
        input.boundary_segment_id[index_of(x, y, input.width)] = segment_id;
    };
    mark(2, 3, 10U);
    mark(4, 3, 10U);
    mark(5, 3, 10U);
    mark(3, 2, 20U);
    mark(3, 4, 20U);
    mark(3, 5, 20U);
    mark(4, 4, 30U);
    mark(5, 5, 30U);
    mark(6, 6, 30U);

    input.plates.resize(3);
    input.plates[0].linear_velocity_x = 1.0f;
    input.plates[0].linear_velocity_y = 0.0f;
    input.plates[1].linear_velocity_x = 0.0f;
    input.plates[1].linear_velocity_y = 2.0f;
    input.plates[2].linear_velocity_x = -1.0f;
    input.plates[2].linear_velocity_y = 1.0f;
    return input;
}

SyntheticJunctionInput make_stable_ridge_junction_input() {
    SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent);

    const auto mark = [&](uint32_t x, uint32_t y, uint32_t segment_id) {
        input.boundary_segment_id[index_of(x, y, input.width)] = segment_id;
    };
    mark(2, 3, 10U);
    mark(4, 3, 10U);
    mark(5, 3, 10U);
    mark(1, 6, 20U);
    mark(2, 5, 20U);
    mark(3, 4, 20U);
    mark(4, 4, 30U);
    mark(5, 5, 30U);
    mark(6, 6, 30U);

    input.plates.resize(3);
    input.plates[0].linear_velocity_x = 1.0f;
    input.plates[0].linear_velocity_y = 0.0f;
    input.plates[1].linear_velocity_x = -1.0f;
    input.plates[1].linear_velocity_y = 0.0f;
    input.plates[2].linear_velocity_x = 0.0f;
    input.plates[2].linear_velocity_y = -1.0f;
    return input;
}

platec::tectonic_junctions::Outputs build(const SyntheticJunctionInput& input) {
    const platec::tectonic_junctions::Inputs detector_input{
        input.width,
        input.height,
        input.plate_id.data(),
        input.boundary_segment_id.data(),
        input.segments.data(),
        static_cast<uint32_t>(input.segments.size()),
        input.plates.data(),
        static_cast<uint32_t>(input.plates.size()),
    };
    return platec::tectonic_junctions::build(detector_input);
}

TectonicJunction synthetic_junction(JunctionStability stability,
                                    JunctionArmType arm0,
                                    JunctionArmType arm1,
                                    JunctionArmType arm2) {
    TectonicJunction junction;
    junction.id = 1U;
    junction.centroid_x = 8.0f;
    junction.centroid_y = 8.0f;
    junction.quality_score = 1.0f;
    junction.geometry_quality_score = 1.0f;
    junction.stability = stability;
    junction.stability_score = stability == JunctionStability::Unstable ? 0.0f : 1.0f;
    junction.arm_types[0] = arm0;
    junction.arm_types[1] = arm1;
    junction.arm_types[2] = arm2;
    junction.arm_strike_x[0] = 1.0f;
    junction.arm_strike_y[0] = 0.0f;
    junction.arm_strike_x[1] = 0.0f;
    junction.arm_strike_y[1] = 1.0f;
    junction.arm_strike_x[2] = 0.70710678f;
    junction.arm_strike_y[2] = 0.70710678f;
    junction.arm_orientation_quality[0] = 1.0f;
    junction.arm_orientation_quality[1] = 1.0f;
    junction.arm_orientation_quality[2] = 1.0f;
    return junction;
}

} // namespace

TEST(TectonicJunctions, DetectsRidgeRidgeRidgeContact) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent);

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_EQ(junction.id, 1U);
    EXPECT_EQ(junction.plate_ids[0], 0U);
    EXPECT_EQ(junction.plate_ids[1], 1U);
    EXPECT_EQ(junction.plate_ids[2], 2U);
    EXPECT_EQ(junction.boundary_segment_ids[0], 10U);
    EXPECT_EQ(junction.boundary_segment_ids[1], 20U);
    EXPECT_EQ(junction.boundary_segment_ids[2], 30U);
    EXPECT_EQ(junction.arm_types[0], JunctionArmType::Ridge);
    EXPECT_EQ(junction.arm_types[1], JunctionArmType::Ridge);
    EXPECT_EQ(junction.arm_types[2], JunctionArmType::Ridge);
    EXPECT_EQ(junction.candidate_cell_count, 1U);
    EXPECT_NEAR(junction.centroid_x, 3.5f, 1e-5f);
    EXPECT_NEAR(junction.centroid_y, 3.5f, 1e-5f);
    EXPECT_GT(junction.quality_score, 0.0f);
}

TEST(TectonicJunctions, ClassifiesFaultFaultFaultContact) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Transform);

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_EQ(output.junctions.front().arm_types[0], JunctionArmType::Fault);
    EXPECT_EQ(output.junctions.front().arm_types[1], JunctionArmType::Fault);
    EXPECT_EQ(output.junctions.front().arm_types[2], JunctionArmType::Fault);
    EXPECT_EQ(output.junctions.front().stability, JunctionStability::Unstable);
    EXPECT_NEAR(output.junctions.front().stability_score, 0.0f, 1e-6f);
}

TEST(TectonicJunctions, ClassifiesMixedRidgeTrenchFaultContact) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Convergent, BoundaryType::Transform);

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_EQ(output.junctions.front().arm_types[0], JunctionArmType::Ridge);
    EXPECT_EQ(output.junctions.front().arm_types[1], JunctionArmType::Trench);
    EXPECT_EQ(output.junctions.front().arm_types[2], JunctionArmType::Fault);
}

TEST(TectonicJunctions, EstimatesLocalArmOrientationAndVelocityMetrics) {
    const SyntheticJunctionInput input = make_oriented_three_plate_input();

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_GT(junction.arm_orientation_quality[0], 0.5f);
    EXPECT_GT(junction.arm_orientation_quality[1], 0.5f);
    EXPECT_GT(junction.arm_orientation_quality[2], 0.5f);
    EXPECT_NEAR(junction.arm_strike_x[0], 1.0f, 1e-5f);
    EXPECT_NEAR(junction.arm_strike_y[0], 0.0f, 1e-5f);
    EXPECT_NEAR(junction.arm_strike_x[1], 0.0f, 1e-5f);
    EXPECT_NEAR(junction.arm_strike_y[1], 1.0f, 1e-5f);
    EXPECT_NEAR(junction.arm_strike_x[2], 0.70710678f, 1e-5f);
    EXPECT_NEAR(junction.arm_strike_y[2], 0.70710678f, 1e-5f);
    EXPECT_NEAR(junction.junction_velocity_x, -0.25f, 1e-5f);
    EXPECT_NEAR(junction.junction_velocity_y, 1.25f, 1e-5f);
    EXPECT_GT(junction.velocity_closure_error, 0.0f);
    EXPECT_LT(junction.velocity_quality_score, 1.0f);
    EXPECT_GT(junction.velocity_quality_score, 0.5f);
}

TEST(TectonicJunctions, PenalizesVelocityQualityWhenArmModesConflict) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{1.0f, -1.0f},
         TestVelocity{-1.0f, 1.0f}});

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_GT(junction.velocity_closure_error, 0.01f);
    EXPECT_LT(junction.velocity_quality_score, 0.95f);
}

TEST(TectonicJunctions, KeepsCompatibleRidgeVelocityQualityHigh) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{0.0f, 0.0f}, TestVelocity{1.0f, 0.0f},
         TestVelocity{0.0f, 1.0f}});

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_NEAR(junction.velocity_closure_error, 0.0f, 1e-5f);
    EXPECT_NEAR(junction.velocity_quality_score, 1.0f, 1e-5f);
}

TEST(TectonicJunctions, DiagnosesStableRidgeJunctionFromVelocityConstraints) {
    const SyntheticJunctionInput input = make_stable_ridge_junction_input();

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_EQ(junction.stability, JunctionStability::Stable);
    EXPECT_NEAR(junction.velocity_constraint_error, 0.0f, 1e-5f);
    EXPECT_GT(junction.geometry_quality_score, 0.5f);
    EXPECT_GT(junction.stability_score, 0.5f);
}

TEST(TectonicJunctions, AllowsCompatibleTransformTransformTrenchJunction) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Convergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{1.0f, -1.0f},
         TestVelocity{-1.0f, 1.0f}});

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_EQ(junction.arm_types[0], JunctionArmType::Fault);
    EXPECT_EQ(junction.arm_types[1], JunctionArmType::Fault);
    EXPECT_EQ(junction.arm_types[2], JunctionArmType::Trench);
    EXPECT_TRUE(junction.stability == JunctionStability::Stable ||
                junction.stability == JunctionStability::ConditionallyStable);
    EXPECT_GT(junction.stability_score, 0.5f);
}

TEST(TectonicJunctions, RejectsInconsistentTransformTransformTrenchJunction) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Convergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{-1.0f, 1.0f},
         TestVelocity{1.0f, -1.0f}});

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_EQ(junction.arm_types[0], JunctionArmType::Fault);
    EXPECT_EQ(junction.arm_types[1], JunctionArmType::Fault);
    EXPECT_EQ(junction.arm_types[2], JunctionArmType::Trench);
    EXPECT_EQ(junction.stability, JunctionStability::Unstable);
    EXPECT_LT(junction.stability_score, 0.5f);
}

TEST(TectonicJunctions, TransformStabilityIsIndependentOfCommonTranslation) {
    SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Convergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{1.0f, -1.0f}, TestVelocity{-1.0f, 1.0f}});
    const auto original = build(input);
    ASSERT_EQ(original.junctions.size(), 1U);
    for (auto& plate : input.plates) {
        plate.linear_velocity_x += 17.0f;
        plate.linear_velocity_y -= 9.0f;
    }
    const auto translated = build(input);
    ASSERT_EQ(translated.junctions.size(), 1U);
    EXPECT_EQ(translated.junctions.front().stability, original.junctions.front().stability);
    EXPECT_GT(translated.junctions.front().stability_score, 0.5f);
}

TEST(TectonicJunctions, RequiresAllThreeIncidentBoundarySegments) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Convergent, BoundaryType::Transform, false);

    const platec::tectonic_junctions::Outputs output = build(input);

    EXPECT_TRUE(output.junctions.empty());
}

TEST(TectonicJunctions, MergesWrappedWorldCandidateCellsAtSeam) {
    const SyntheticJunctionInput input = make_wrapped_candidate_input();

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_GE(junction.candidate_cell_count, 2U);
    EXPECT_TRUE(junction.centroid_x < 1.0f ||
                junction.centroid_x > static_cast<float>(input.width) - 1.0f);
}

TEST(TectonicJunctions, UsesWrappedCentroidForAngularLocalVelocity) {
    SyntheticJunctionInput input = make_wrapped_candidate_input();
    const platec::tectonic_junctions::Outputs probe = build(input);
    ASSERT_EQ(probe.junctions.size(), 1U);
    const float centroid_x = probe.junctions.front().centroid_x;
    const float centroid_y = probe.junctions.front().centroid_y;
    const float mass_center_x =
        centroid_x < static_cast<float>(input.width) * 0.5f
            ? static_cast<float>(input.width) - 0.5f
            : 0.5f;

    input.plates.resize(3);
    for (PlateKinematics& plate : input.plates) {
        plate.angular_velocity = 1.0f;
        plate.mass_center_x = mass_center_x;
        plate.mass_center_y = centroid_y;
    }

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    const float dx = wrapped_delta_for_test(junction.centroid_x, mass_center_x, input.width);
    const float unwrapped_dx = junction.centroid_x - mass_center_x;
    ASSERT_GT(std::fabs(unwrapped_dx), static_cast<float>(input.width) * 0.5f);
    EXPECT_LE(std::fabs(dx), 1.0f);
    EXPECT_NEAR(junction.junction_velocity_x, 0.0f, 1e-4f);
    EXPECT_NEAR(junction.junction_velocity_y, dx, 1e-4f);
}

TEST(TectonicJunctions, ClustersNeighboringCandidatesWithSameSegments) {
    const SyntheticJunctionInput input = make_clustered_candidate_input();

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_GE(output.junctions.front().candidate_cell_count, 2U);
}

TEST(TectonicJunctions, DetectsTriplePointWhenBoundaryMapCollapsesOnePair) {
    SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{}, TestVelocity{}, TestVelocity{}});
    for (uint32_t& segment_id : input.boundary_segment_id) {
        if (segment_id == 30U) {
            segment_id = platec::contract::kNoBoundaryId;
        }
    }

    const platec::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_EQ(output.junctions.front().boundary_segment_ids[2], 30U);
    EXPECT_EQ(output.junctions.front().stability, JunctionStability::Indeterminate);
}

TEST(TectonicJunctions, DoesNotGuessBetweenUnrepresentedBoundarySegments) {
    SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{}, TestVelocity{}, TestVelocity{}});
    for (auto& id : input.boundary_segment_id) {
        if (id == 30U) {
            id = platec::contract::kNoBoundaryId;
        }
    }
    input.segments.push_back(segment(31U, 1U, 2U, BoundaryType::Divergent));
    EXPECT_TRUE(build(input).junctions.empty());
}

TEST(TectonicJunctions, DoesNotBorrowRemoteRasterSegmentForMissingArm) {
    SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{}, TestVelocity{}, TestVelocity{}});
    for (auto& id : input.boundary_segment_id) {
        if (id == 30U) {
            id = platec::contract::kNoBoundaryId;
        }
    }
    mark_segment_cell(input, 0U, 0U, 30U);
    EXPECT_TRUE(build(input).junctions.empty());
}

TEST(TectonicJunctions, RejectsFourPlateInstantaneousContact) {
    SyntheticJunctionInput input;
    input.width = 2;
    input.height = 2;
    input.plate_id = {0U, 1U, 2U, 3U};
    input.boundary_segment_id = {10U, 20U, 30U, 40U};
    input.segments = {
        segment(10U, 0U, 1U, BoundaryType::Divergent),
        segment(20U, 0U, 2U, BoundaryType::Convergent),
        segment(30U, 1U, 3U, BoundaryType::Transform),
        segment(40U, 2U, 3U, BoundaryType::PassiveMargin),
    };

    const platec::tectonic_junctions::Outputs output = build(input);

    EXPECT_TRUE(output.junctions.empty());
}

TEST(TectonicJunctionCoupling, StableRidgeJunctionCreatesRiftInfluenceAndTension) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Stable, JunctionArmType::Ridge,
                           JunctionArmType::Ridge, JunctionArmType::Ridge);
    const platec::tectonic_junction_coupling::Inputs input{16U, 16U, &junction, 1U};

    const platec::tectonic_junction_coupling::Outputs output =
        platec::tectonic_junction_coupling::build(input);

    const size_t east = index_of(10U, 8U, 16U);
    EXPECT_GT(output.rift_influence[east], 0.0f);
    EXPECT_EQ(output.transform_influence[east], 0.0f);
    EXPECT_EQ(output.subduction_influence[east], 0.0f);
    EXPECT_GT(output.rift_tension_x[east], 0.0f);
}

TEST(TectonicJunctionCoupling, StableTrenchJunctionCreatesSubductionInfluence) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Stable, JunctionArmType::Trench,
                           JunctionArmType::Trench, JunctionArmType::Ridge);
    const platec::tectonic_junction_coupling::Inputs input{16U, 16U, &junction, 1U};

    const platec::tectonic_junction_coupling::Outputs output =
        platec::tectonic_junction_coupling::build(input);

    const size_t center = index_of(8U, 8U, 16U);
    EXPECT_GT(output.subduction_influence[center], 0.0f);
    EXPECT_GT(output.rift_influence[center], 0.0f);
    EXPECT_EQ(output.transform_influence[center], 0.0f);
}

TEST(TectonicJunctionCoupling, UnstableFaultJunctionCreatesTransformInfluenceOnly) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Unstable, JunctionArmType::Fault,
                           JunctionArmType::Fault, JunctionArmType::Fault);
    const platec::tectonic_junction_coupling::Inputs input{16U, 16U, &junction, 1U};

    const platec::tectonic_junction_coupling::Outputs output =
        platec::tectonic_junction_coupling::build(input);

    const size_t center = index_of(8U, 8U, 16U);
    EXPECT_EQ(output.rift_influence[center], 0.0f);
    EXPECT_EQ(output.subduction_influence[center], 0.0f);
    EXPECT_GT(output.transform_influence[center], 0.0f);
}
