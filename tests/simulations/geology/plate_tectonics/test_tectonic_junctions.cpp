#ifndef NOMINMAX
#define NOMINMAX
#endif

// Triple junctions on the cube (W2-E, E3): the synthetic layouts of the earlier tests drawn on
// face 0 (column i east, row j south, the layout's centre at the face's centre (longitude 0,
// latitude 0)), the other faces plate 0; plate velocities (east, south in cells per update) as
// the rotations that give them at the face's centre; junction centroids, strikes and velocities
// in longitude, latitude and (east, north).
#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junction_coupling.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junctions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace {

using tectonics::contract::BoundarySegment;
using tectonics::contract::BoundaryType;
using tectonics::contract::JunctionArmType;
using tectonics::contract::JunctionStability;
using tectonics::contract::TectonicJunction;

struct TestVelocity {
    float x = 0.0f;  // east
    float y = 0.0f;  // south
};

struct SyntheticJunctionInput {
    uint32_t width = 8;
    uint32_t height = 8;
    std::vector<uint32_t> plate_id;
    std::vector<uint32_t> boundary_segment_id;
    std::vector<BoundarySegment> segments;
    std::vector<TestVelocity> plates;
    std::vector<mesh::Vec3> extra;  // per plate, a rotation (radians per update) added to the velocity's
};

size_t index_of(uint32_t x, uint32_t y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
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
    input.plates[plate_id] = velocity;
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
    input.boundary_segment_id.assign(input.plate_id.size(), tectonics::contract::kNoBoundaryId);

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
    input.boundary_segment_id.assign(input.plate_id.size(), tectonics::contract::kNoBoundaryId);
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
    input.boundary_segment_id.assign(input.plate_id.size(), tectonics::contract::kNoBoundaryId);
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
    input.boundary_segment_id.assign(input.plate_id.size(), tectonics::contract::kNoBoundaryId);

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

    input.plates = {TestVelocity{1.0f, 0.0f}, TestVelocity{0.0f, 2.0f}, TestVelocity{-1.0f, 1.0f}};
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

    input.plates = {TestVelocity{1.0f, 0.0f}, TestVelocity{-1.0f, 0.0f}, TestVelocity{0.0f, -1.0f}};
    return input;
}

// The layout on face 0 of a cube of N = max(width, height) cells per face edge.
struct OnCube {
    std::shared_ptr<const tectonics::grid::TectonicGrid> grid;
    std::vector<uint32_t> plate_id, segment_id;
    std::vector<mesh::Vec3> rotation;
};
OnCube oncube(const SyntheticJunctionInput& input) {
    const int n = static_cast<int>(std::max(input.width, input.height));
    OnCube cube;
    cube.grid = tectonics::grid::TectonicGrid::cube(n, 6371.0);
    cube.plate_id.assign(cube.grid->cellcount(), 0U);
    cube.segment_id.assign(cube.grid->cellcount(), tectonics::contract::kNoBoundaryId);
    for (uint32_t y = 0; y < input.height; ++y)
        for (uint32_t x = 0; x < input.width; ++x) {
            const auto c = static_cast<size_t>(y) * n + x;
            cube.plate_id[c] = input.plate_id[index_of(x, y, input.width)];
            cube.segment_id[c] = input.boundary_segment_id[index_of(x, y, input.width)];
        }
    // At (1, 0, 0): east (0, 1, 0), north (0, 0, 1); w = p x v / R gives v there.
    const double scale = cube.grid->cellspacing_m() / cube.grid->radius_m();
    uint32_t plates = 0;
    for (const uint32_t id : input.plate_id) plates = std::max(plates, id + 1);
    plates = std::max<uint32_t>(plates, static_cast<uint32_t>(input.plates.size()));
    cube.rotation.assign(plates, mesh::Vec3{});
    for (size_t i = 0; i < input.plates.size(); ++i)
        cube.rotation[i] = scale * mesh::Vec3{0.0, static_cast<double>(input.plates[i].y), static_cast<double>(input.plates[i].x)};
    for (size_t i = 0; i < input.extra.size() && i < cube.rotation.size(); ++i)
        cube.rotation[i] = cube.rotation[i] + input.extra[i];
    return cube;
}

tectonics::tectonic_junctions::Outputs build(const SyntheticJunctionInput& input) {
    const OnCube cube = oncube(input);
    const tectonics::tectonic_junctions::Inputs detector_input{
        cube.plate_id.data(),
        cube.segment_id.data(),
        input.segments.data(),
        static_cast<uint32_t>(input.segments.size()),
        static_cast<uint32_t>(cube.rotation.size()),
        cube.grid.get(),
        cube.rotation.data(),
    };
    return tectonics::tectonic_junctions::build(detector_input);
}

TectonicJunction synthetic_junction(JunctionStability stability,
                                    JunctionArmType arm0,
                                    JunctionArmType arm1,
                                    JunctionArmType arm2) {
    TectonicJunction junction;
    junction.id = 1U;
    junction.centroid_x = 0.0f;  // longitude, latitude: face 0's centre
    junction.centroid_y = 0.0f;
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

    const tectonics::tectonic_junctions::Outputs output = build(input);

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
    // The vertex between the layout's cells 3 and 4 is face 0's centre.
    EXPECT_NEAR(junction.centroid_x, 0.0f, 1e-4f);
    EXPECT_NEAR(junction.centroid_y, 0.0f, 1e-4f);
    EXPECT_GT(junction.quality_score, 0.0f);
}

TEST(TectonicJunctions, ClassifiesFaultFaultFaultContact) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Transform);

    const tectonics::tectonic_junctions::Outputs output = build(input);

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

    const tectonics::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_EQ(output.junctions.front().arm_types[0], JunctionArmType::Ridge);
    EXPECT_EQ(output.junctions.front().arm_types[1], JunctionArmType::Trench);
    EXPECT_EQ(output.junctions.front().arm_types[2], JunctionArmType::Fault);
}

TEST(TectonicJunctions, EstimatesLocalArmOrientationAndVelocityMetrics) {
    const SyntheticJunctionInput input = make_oriented_three_plate_input();

    const tectonics::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_GT(junction.arm_orientation_quality[0], 0.5f);
    EXPECT_GT(junction.arm_orientation_quality[1], 0.5f);
    EXPECT_GT(junction.arm_orientation_quality[2], 0.5f);
    // East-west, north-south and south-east arms (east, north; the strike's sign has east >= 0).
    EXPECT_NEAR(junction.arm_strike_x[0], 1.0f, 2e-2f);
    EXPECT_NEAR(junction.arm_strike_y[0], 0.0f, 2e-2f);
    EXPECT_NEAR(junction.arm_strike_x[1], 0.0f, 2e-2f);
    EXPECT_NEAR(std::fabs(junction.arm_strike_y[1]), 1.0f, 2e-2f);
    EXPECT_NEAR(junction.arm_strike_x[2], 0.70710678f, 3e-2f);
    EXPECT_NEAR(junction.arm_strike_y[2], -0.70710678f, 3e-2f);
    EXPECT_NEAR(junction.junction_velocity_x, -0.25f, 3e-2f);
    EXPECT_NEAR(junction.junction_velocity_y, -1.25f, 3e-2f);
    EXPECT_GT(junction.velocity_closure_error, 0.0f);
    EXPECT_LT(junction.velocity_quality_score, 1.0f);
    EXPECT_GT(junction.velocity_quality_score, 0.5f);
}

TEST(TectonicJunctions, PenalizesVelocityQualityWhenArmModesConflict) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{1.0f, -1.0f},
         TestVelocity{-1.0f, 1.0f}});

    const tectonics::tectonic_junctions::Outputs output = build(input);

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

    const tectonics::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_NEAR(junction.velocity_closure_error, 0.0f, 1e-5f);
    EXPECT_NEAR(junction.velocity_quality_score, 1.0f, 1e-5f);
}

TEST(TectonicJunctions, DiagnosesStableRidgeJunctionFromVelocityConstraints) {
    const SyntheticJunctionInput input = make_stable_ridge_junction_input();

    const tectonics::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_EQ(junction.stability, JunctionStability::Stable);
    EXPECT_NEAR(junction.velocity_constraint_error, 0.0f, 1e-5f);
    // About 0.5: arms at 45, 45 and 90 degrees against the ideal 60, read in (east, north) on
    // the cube's face, whose cells are not quite square.
    EXPECT_GT(junction.geometry_quality_score, 0.45f);
    EXPECT_GT(junction.stability_score, 0.45f);
}

TEST(TectonicJunctions, AllowsCompatibleTransformTransformTrenchJunction) {
    const SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Transform, BoundaryType::Transform, BoundaryType::Convergent,
        {TestVelocity{1.0f, 1.0f}, TestVelocity{1.0f, -1.0f},
         TestVelocity{-1.0f, 1.0f}});

    const tectonics::tectonic_junctions::Outputs output = build(input);

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

    const tectonics::tectonic_junctions::Outputs output = build(input);

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
    // A common rotation of every plate (17 cells per update east and 9 north at the junction).
    input.extra.assign(3, oncube(input).grid->cellspacing_m() / oncube(input).grid->radius_m() * mesh::Vec3{0.0, -9.0, 17.0});
    const auto translated = build(input);
    ASSERT_EQ(translated.junctions.size(), 1U);
    EXPECT_EQ(translated.junctions.front().stability, original.junctions.front().stability);
    EXPECT_GT(translated.junctions.front().stability_score, 0.5f);
}

TEST(TectonicJunctions, RequiresAllThreeIncidentBoundarySegments) {
    const SyntheticJunctionInput input = make_three_plate_input(
        BoundaryType::Divergent, BoundaryType::Convergent, BoundaryType::Transform, false);

    const tectonics::tectonic_junctions::Outputs output = build(input);

    EXPECT_TRUE(output.junctions.empty());
}

TEST(TectonicJunctions, MergesCubeCandidateCellsAcrossAFaceEdge) {
    // The wrapped pattern of make_wrapped_candidate_input on an N 12 cube: its columns 6..11
    // are face 0's, columns 0..5 face 1's (the same rows), so the seam is the face 0 / face 1
    // edge at 45 E.
    constexpr int N = 12;
    const auto table = tectonics::grid::TectonicGrid::cube(N, 6371.0);
    const SyntheticJunctionInput planar = make_wrapped_candidate_input();
    const auto cell = [&](uint32_t x, uint32_t y) {
        return static_cast<size_t>(x >= 6U ? 0 : N * N) + static_cast<size_t>(y) * N + x;
    };
    std::vector<uint32_t> plate_id(table->cellcount(), 0U);
    std::vector<uint32_t> segment_id(table->cellcount(), tectonics::contract::kNoBoundaryId);
    for (uint32_t y = 0; y < planar.height; ++y)
        for (uint32_t x = 0; x < planar.width; ++x) {
            plate_id[cell(x, y)] = planar.plate_id[index_of(x, y, planar.width)];
            segment_id[cell(x, y)] = planar.boundary_segment_id[index_of(x, y, planar.width)];
        }
    // Every plate turns eastward about +Z at one cell per update on the equator.
    const double spacing = std::sqrt(table->totalarea_m2() / static_cast<double>(table->cellcount()));
    const std::vector<mesh::Vec3> rotation(3, mesh::Vec3{0.0, 0.0, spacing / 6371000.0});
    const tectonics::tectonic_junctions::Inputs inputs{
        plate_id.data(), segment_id.data(), planar.segments.data(),
        static_cast<uint32_t>(planar.segments.size()), 3U, table.get(),
        rotation.data()};

    const tectonics::tectonic_junctions::Outputs output = tectonics::tectonic_junctions::build(inputs);

    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    // The layout's four candidate vertices: one on the face edge (cells of both faces), three
    // on face 1, merged into one junction between 45 E and the third column (60 E).
    EXPECT_EQ(junction.candidate_cell_count, 4U);
    EXPECT_EQ(junction.plate_ids[0], 0U);
    EXPECT_EQ(junction.plate_ids[1], 1U);
    EXPECT_EQ(junction.plate_ids[2], 2U);
    EXPECT_EQ(junction.boundary_segment_ids[0], 10U);
    EXPECT_EQ(junction.boundary_segment_ids[1], 20U);
    EXPECT_EQ(junction.boundary_segment_ids[2], 30U);
    EXPECT_GT(junction.centroid_x, 45.0f);
    EXPECT_LT(junction.centroid_x, 60.0f);
    EXPECT_GT(junction.centroid_y, 0.0f);
    EXPECT_LT(junction.centroid_y, 7.5f);
    // The common rotation's velocity at the junction: east at cos(latitude) cells per update.
    const float latitude = junction.centroid_y * 3.14159265f / 180.0f;
    EXPECT_NEAR(junction.junction_velocity_x, std::cos(latitude), 1e-4f);
    EXPECT_NEAR(junction.junction_velocity_y, 0.0f, 1e-4f);
}

TEST(TectonicJunctionCoupling, CubeRiftTensionPointsAwayAcrossAFaceEdge) {
    constexpr int N = 12;
    const auto table = tectonics::grid::TectonicGrid::cube(N, 6371.0);
    TectonicJunction junction =
        synthetic_junction(JunctionStability::Stable, JunctionArmType::Ridge,
                           JunctionArmType::Ridge, JunctionArmType::Ridge);
    // The vertex of face 0's column 11 and face 1's column 0 at rows 5 / 6 (the equator).
    junction.centroid_x = 45.0f;
    junction.centroid_y = 0.0f;
    const tectonics::tectonic_junction_coupling::Inputs input{&junction, 1U, table.get()};

    const tectonics::tectonic_junction_coupling::Outputs output =
        tectonics::tectonic_junction_coupling::build(input);

    ASSERT_EQ(output.rift_influence.size(), table->cellcount());
    const size_t face0 = 5U * N + 11U;              // north-west of the vertex
    const size_t face1 = static_cast<size_t>(N * N) + 5U * N;  // north-east of the vertex
    ASSERT_EQ(table->neighbour(static_cast<tectonics::grid::CellId>(face0), tectonics::grid::east), face1);
    EXPECT_GT(output.rift_influence[face0], 0.0f);
    EXPECT_GT(output.rift_influence[face1], 0.0f);
    EXPECT_NEAR(output.rift_influence[face0], output.rift_influence[face1], 1e-5f);
    EXPECT_LT(output.rift_tension_x[face0], 0.0f);
    EXPECT_GT(output.rift_tension_x[face1], 0.0f);
    EXPECT_GT(output.rift_tension_y[face0], 0.0f);
    EXPECT_GT(output.rift_tension_y[face1], 0.0f);
    EXPECT_EQ(output.rift_influence[3U * N * N], 0.0f);  // face 3, far away
}

// A junction on the date line (the layout drawn on face 2, centred on 180 degrees), its
// plates turning about a pole a cell from it: the junction's velocity is each plate's
// omega x p at the junction itself (their mean), whichever side of the date line the
// candidates lie on.
TEST(TectonicJunctions, UsesTheJunctionPointForRotationalVelocityAcrossTheDateLine) {
    constexpr int N = 12;
    const auto table = tectonics::grid::TectonicGrid::cube(N, 6371.0);
    const SyntheticJunctionInput planar = make_clustered_candidate_input();
    std::vector<uint32_t> plate_id(table->cellcount(), 0U);
    std::vector<uint32_t> segment_id(table->cellcount(), tectonics::contract::kNoBoundaryId);
    // Face 2 (-X): column i east, row j south as on face 0; its centre is longitude 180.
    for (uint32_t y = 0; y < planar.height; ++y)
        for (uint32_t x = 0; x < planar.width; ++x) {
            const size_t c = static_cast<size_t>(2 * N * N) + static_cast<size_t>(y) * N + x;
            plate_id[c] = planar.plate_id[index_of(x, y, planar.width)];
            segment_id[c] = planar.boundary_segment_id[index_of(x, y, planar.width)];
        }
    const double spacing = table->cellspacing_m(), radius = table->radius_m();
    // The pole one cell spacing east of the date line on the equator.
    const mesh::Vec3 pole = mesh::unitvector(0.0, -3.14159265358979323846 + spacing / radius);
    const std::vector<mesh::Vec3> rotation(3, (2 * spacing / radius) * pole);
    const tectonics::tectonic_junctions::Inputs inputs{
        plate_id.data(), segment_id.data(), planar.segments.data(),
        static_cast<uint32_t>(planar.segments.size()), 3U, table.get(), rotation.data()};
    const tectonics::tectonic_junctions::Outputs output = tectonics::tectonic_junctions::build(inputs);
    ASSERT_EQ(output.junctions.size(), 1U);
    const auto& junction = output.junctions.front();
    EXPECT_GT(std::fabs(junction.centroid_x), 170.0f);
    const mesh::Vec3 p = tectonics::boundary_cube::decodecentroid(junction.centroid_x, junction.centroid_y);
    const auto expected = tectonics::boundary_cube::velocity(*table, rotation.data(), 0, 3, p);
    EXPECT_GT(std::hypot(expected[0], expected[1]), .1f);
    EXPECT_NEAR(junction.junction_velocity_x, expected[0], 1e-4f);
    EXPECT_NEAR(junction.junction_velocity_y, expected[1], 1e-4f);
}

TEST(TectonicJunctions, ClustersNeighboringCandidatesWithSameSegments) {
    const SyntheticJunctionInput input = make_clustered_candidate_input();

    const tectonics::tectonic_junctions::Outputs output = build(input);

    ASSERT_EQ(output.junctions.size(), 1U);
    EXPECT_GE(output.junctions.front().candidate_cell_count, 2U);
}

TEST(TectonicJunctions, DetectsTriplePointWhenBoundaryMapCollapsesOnePair) {
    SyntheticJunctionInput input = make_explicit_arm_input(
        BoundaryType::Divergent, BoundaryType::Divergent, BoundaryType::Divergent,
        {TestVelocity{}, TestVelocity{}, TestVelocity{}});
    for (uint32_t& segment_id : input.boundary_segment_id) {
        if (segment_id == 30U) {
            segment_id = tectonics::contract::kNoBoundaryId;
        }
    }

    const tectonics::tectonic_junctions::Outputs output = build(input);

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
            id = tectonics::contract::kNoBoundaryId;
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
            id = tectonics::contract::kNoBoundaryId;
        }
    }
    mark_segment_cell(input, 0U, 0U, 30U);
    EXPECT_TRUE(build(input).junctions.empty());
}

TEST(TectonicJunctions, RejectsFourPlateInstantaneousContact) {
    // Four plates meet at face 0's centre vertex (the layout's cells 3 and 4 in both
    // directions); everything else is plate 4.
    SyntheticJunctionInput input;
    input.plate_id.assign(64, 4U);
    input.boundary_segment_id.assign(64, tectonics::contract::kNoBoundaryId);
    const uint32_t cells[4] = {static_cast<uint32_t>(index_of(3, 3, 8)), static_cast<uint32_t>(index_of(4, 3, 8)),
                               static_cast<uint32_t>(index_of(3, 4, 8)), static_cast<uint32_t>(index_of(4, 4, 8))};
    for (uint32_t k = 0; k < 4; ++k) {
        input.plate_id[cells[k]] = k;
        input.boundary_segment_id[cells[k]] = 10U * (k + 1);
    }
    input.segments = {
        segment(10U, 0U, 1U, BoundaryType::Divergent),
        segment(20U, 0U, 2U, BoundaryType::Convergent),
        segment(30U, 1U, 3U, BoundaryType::Transform),
        segment(40U, 2U, 3U, BoundaryType::PassiveMargin),
    };

    const tectonics::tectonic_junctions::Outputs output = build(input);

    EXPECT_TRUE(output.junctions.empty());
}

TEST(TectonicJunctionCoupling, StableRidgeJunctionCreatesRiftInfluenceAndTension) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Stable, JunctionArmType::Ridge,
                           JunctionArmType::Ridge, JunctionArmType::Ridge);
    const auto table = tectonics::grid::TectonicGrid::cube(16, 6371.0);
    const tectonics::tectonic_junction_coupling::Inputs input{&junction, 1U, table.get()};

    const tectonics::tectonic_junction_coupling::Outputs output =
        tectonics::tectonic_junction_coupling::build(input);

    const size_t east = index_of(10U, 8U, 16U);  // two cells east of the junction, face 0
    EXPECT_GT(output.rift_influence[east], 0.0f);
    EXPECT_EQ(output.transform_influence[east], 0.0f);
    EXPECT_EQ(output.subduction_influence[east], 0.0f);
    EXPECT_GT(output.rift_tension_x[east], 0.0f);
}

TEST(TectonicJunctionCoupling, StableTrenchJunctionCreatesSubductionInfluence) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Stable, JunctionArmType::Trench,
                           JunctionArmType::Trench, JunctionArmType::Ridge);
    const auto table = tectonics::grid::TectonicGrid::cube(16, 6371.0);
    const tectonics::tectonic_junction_coupling::Inputs input{&junction, 1U, table.get()};

    const tectonics::tectonic_junction_coupling::Outputs output =
        tectonics::tectonic_junction_coupling::build(input);

    const size_t center = index_of(8U, 8U, 16U);  // beside the junction, face 0
    EXPECT_GT(output.subduction_influence[center], 0.0f);
    EXPECT_GT(output.rift_influence[center], 0.0f);
    EXPECT_EQ(output.transform_influence[center], 0.0f);
}

TEST(TectonicJunctionCoupling, UnstableFaultJunctionCreatesTransformInfluenceOnly) {
    const TectonicJunction junction =
        synthetic_junction(JunctionStability::Unstable, JunctionArmType::Fault,
                           JunctionArmType::Fault, JunctionArmType::Fault);
    const auto table = tectonics::grid::TectonicGrid::cube(16, 6371.0);
    const tectonics::tectonic_junction_coupling::Inputs input{&junction, 1U, table.get()};

    const tectonics::tectonic_junction_coupling::Outputs output =
        tectonics::tectonic_junction_coupling::build(input);

    const size_t center = index_of(8U, 8U, 16U);  // beside the junction, face 0
    EXPECT_EQ(output.rift_influence[center], 0.0f);
    EXPECT_EQ(output.subduction_influence[center], 0.0f);
    EXPECT_GT(output.transform_influence[center], 0.0f);
}

// The candidate points are the cube's 6 N^2 + 2 mesh vertices: four cells each, three at the
// eight cube corners, every cell at four of them.
TEST(TectonicJunctions, CubeVerticesAreTheMeshVerticesEachOnce) {
    for (const int n : {2, 3, 8}) {
        const auto table = tectonics::grid::TectonicGrid::cube(n, 6371.0);
        const auto vertices = tectonics::tectonic_junctions::cube_vertices(*table);
        EXPECT_EQ(vertices.size(), 6U * n * n + 2U) << n;
        std::vector<int> corners(table->cellcount(), 0);
        size_t three = 0;
        for (const auto& vertex : vertices) {
            three += vertex.count == 3U;
            ASSERT_TRUE(vertex.count == 3U || vertex.count == 4U);
            for (uint32_t i = 0; i < vertex.count; ++i) ++corners[vertex.cells[i]];
        }
        EXPECT_EQ(three, 8U) << n;
        for (const int k : corners) EXPECT_EQ(k, 4);
    }
}

// Three plates in 120-degree sectors about the cube's diagonal (1, 1, 1), each holding one of
// the three cells at the corner (1, 1, 1), where a vertex has only three cells: a junction there
// with the three arms (the sectors meet again at (-1, -1, -1), where the one segment id per cell
// may leave an arm unrepresented).
TEST(TectonicJunctions, FindsAThreePlateJunctionAtACubeCorner) {
    constexpr int N = 8;
    const auto table = tectonics::grid::TectonicGrid::cube(N, 6371.0);
    const mesh::Vec3 d = mesh::normalise({1, 1, 1});
    const mesh::Vec3 e1 = mesh::normalise(mesh::Vec3{1, 0, 0} - mesh::dot({1, 0, 0}, d) * d), e2 = mesh::cross(d, e1);
    const double pi = 3.14159265358979323846;
    std::vector<uint32_t> plate_id(table->cellcount());
    for (tectonics::grid::CellId c = 0; c < table->cellcount(); ++c) {
        const mesh::Vec3 p = table->position(c);
        // Sectors centred at 30, 150 and 270 degrees from e1: the corner cells (at 0, 120 and
        // 240 about +d, at 60, 180 and 300 about -d) lie inside them.
        double angle = std::atan2(mesh::dot(p, e2), mesh::dot(p, e1)) + pi / 6;
        while (angle < 0) angle += 2 * pi;
        plate_id[c] = std::min(2U, static_cast<uint32_t>(angle / (2 * pi / 3)));
    }
    const std::vector<BoundarySegment> segments{segment(10U, 0U, 1U, BoundaryType::Divergent),
                                                segment(20U, 0U, 2U, BoundaryType::Divergent),
                                                segment(30U, 1U, 2U, BoundaryType::Divergent)};
    std::vector<uint32_t> segment_id(table->cellcount(), tectonics::contract::kNoBoundaryId);
    for (tectonics::grid::CellId c = 0; c < table->cellcount(); ++c)
        for (const auto m : table->neighbours(c))
            if (plate_id[m] != plate_id[c]) {
                const uint32_t lo = std::min(plate_id[c], plate_id[m]), hi = std::max(plate_id[c], plate_id[m]);
                segment_id[c] = lo == 0 ? (hi == 1 ? 10U : 20U) : 30U;
                break;
            }
    const std::vector<mesh::Vec3> rotation(3, mesh::Vec3{});
    const tectonics::tectonic_junctions::Inputs inputs{plate_id.data(), segment_id.data(), segments.data(), 3U, 3U,
                                                       table.get(), rotation.data()};
    const auto output = tectonics::tectonic_junctions::build(inputs);
    const double corner_lat = std::asin(1 / std::sqrt(3.0)) * 180 / pi;
    int found = 0;
    // Longitudes compared modulo 360 degrees.
    const auto near_longitude = [](double longitude, double target) {
        return std::abs(std::remainder(longitude - target, 360.0)) < 3;
    };
    for (const auto& junction : output.junctions) {
        if (!near_longitude(junction.centroid_x, 45) || std::abs(junction.centroid_y - corner_lat) >= 3) continue;
        ++found;
        EXPECT_EQ(junction.plate_ids[0], 0U);
        EXPECT_EQ(junction.plate_ids[1], 1U);
        EXPECT_EQ(junction.plate_ids[2], 2U);
        EXPECT_EQ(junction.boundary_segment_ids[0], 10U);
        EXPECT_EQ(junction.boundary_segment_ids[1], 20U);
        EXPECT_EQ(junction.boundary_segment_ids[2], 30U);
    }
    std::string where;
    for (const auto& junction : output.junctions)
        where += " (" + std::to_string(junction.centroid_x) + ", " + std::to_string(junction.centroid_y) + ")";
    EXPECT_EQ(found, 1) << output.junctions.size() << " junctions:" << where;
}
