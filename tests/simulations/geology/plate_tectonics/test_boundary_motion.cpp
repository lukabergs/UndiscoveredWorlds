#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_boundary_graph.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deforming_regions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace {
using namespace tectonics;
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

namespace {
// Face 0 of an N 4 cube is plate 0, spinning eastward about +Z at 4 cells per update on the
// equator; the other five faces are plate 1 at rest. Face 0's east edge converges with face 1,
// its west edge diverges from face 3, its top and bottom edges shear against faces 4 and 5,
// and its four corner cells touch two faces each.
struct CubeBoundaryFixture {
    static constexpr int N = 4;
    std::shared_ptr<const grid::TectonicGrid> table = grid::TectonicGrid::cube(N, 6371.0);
    size_t n = table->cellcount();
    std::vector<float> elevation = std::vector<float>(n, 0.5f);
    std::vector<uint32_t> owner = std::vector<uint32_t>(n, 1U), zero = std::vector<uint32_t>(n, 0U);
    std::vector<uint8_t> convergence = std::vector<uint8_t>(n), divergence = convergence,
                         shear = convergence, regime = convergence;
    std::vector<uint8_t> crust = std::vector<uint8_t>(n, uint8_t(CrustClass::Oceanic));
    std::array<contract::PlateKinematics, 2> plates{};
    std::array<mesh::Vec3, 2> rotation{};

    CubeBoundaryFixture() {
        for (grid::CellId c = 0; c < N * N; ++c) owner[c] = 0U;
        rotation[0] = {0.0, 0.0, 4.0 * spacing() / (6371.0 * 1000.0)};
    }
    double spacing() const { return std::sqrt(table->totalarea_m2() / static_cast<double>(n)); }
    static grid::CellId face0(int i, int j) { return static_cast<grid::CellId>(j * N + i); }
    void compute() {
        provenance::Inputs inputs{0, 0, elevation.data(), owner.data(), plates.data(), 2,
                                  crust.data(), 0.0, table.get(), rotation.data()};
        provenance::compute_maps(inputs, convergence.data(), divergence.data(), shear.data(), regime.data());
    }
    boundary_graph::Outputs graph() const {
        boundary_graph::Inputs inputs{0, 0, owner.data(), convergence.data(),
            divergence.data(), shear.data(), regime.data(), crust.data(), zero.data(),
            zero.data(), nullptr, 0, 1, 1.0, table.get()};
        return boundary_graph::build(inputs);
    }
};
}

TEST(BoundaryMotion, CubeContactsCarryMotionAcrossFaceEdgesAndCorners) {
    CubeBoundaryFixture f;
    f.compute();
    const int N = CubeBoundaryFixture::N;
    for (int j = 0; j < N; ++j) {
        const grid::CellId east_edge = f.face0(N - 1, j), west_edge = f.face0(0, j);
        const grid::CellId face1 = f.table->neighbour(east_edge, grid::east);
        const grid::CellId face3 = f.table->neighbour(west_edge, grid::west);
        ASSERT_EQ(face1 / (N * N), 1U);
        ASSERT_EQ(face3 / (N * N), 3U);
        EXPECT_EQ(f.convergence[east_edge], 100) << "j " << j;
        EXPECT_EQ(f.convergence[face1], 100) << "j " << j;
        EXPECT_EQ(f.divergence[west_edge], 100) << "j " << j;
        EXPECT_EQ(f.divergence[face3], 100) << "j " << j;
    }
    for (int i = 0; i < N; ++i) {
        const grid::CellId top = f.face0(i, 0), bottom = f.face0(i, N - 1);
        const grid::CellId face4 = f.table->neighbour(top, grid::north);
        const grid::CellId face5 = f.table->neighbour(bottom, grid::south);
        ASSERT_EQ(face4 / (N * N), 4U);
        ASSERT_EQ(face5 / (N * N), 5U);
        EXPECT_EQ(f.shear[top], 100) << "i " << i;
        EXPECT_EQ(f.shear[face4], 100) << "i " << i;
        EXPECT_EQ(f.shear[bottom], 100) << "i " << i;
        EXPECT_EQ(f.shear[face5], 100) << "i " << i;
    }
    // A cube corner cell: converging with face 1 and shearing against face 4.
    const grid::CellId corner = f.face0(N - 1, 0);
    EXPECT_EQ(f.convergence[corner], 100);
    EXPECT_EQ(f.shear[corner], 100);
    EXPECT_EQ(f.regime[corner], uint8_t(GeologicRegime::ConvergentArc));
    // Interior cells and the far faces stay quiet.
    EXPECT_EQ(f.convergence[f.face0(1, 1)] + f.divergence[f.face0(1, 1)] + f.shear[f.face0(1, 1)], 0);
    EXPECT_EQ(f.convergence[2 * N * N] + f.divergence[2 * N * N] + f.shear[2 * N * N], 0);

    // The convergent edge is one segment across the face edge, centred on (45 E, 0 N).
    const auto graph = f.graph();
    const uint32_t id = graph.boundary_segment_id[f.face0(N - 1, 1)];
    ASSERT_NE(id, contract::kNoBoundaryId);
    for (int j = 0; j < N; ++j) {
        EXPECT_EQ(graph.boundary_segment_id[f.face0(N - 1, j)], id) << "j " << j;
        EXPECT_EQ(graph.boundary_segment_id[f.table->neighbour(f.face0(N - 1, j), grid::east)], id) << "j " << j;
    }
    const auto segment = std::find_if(graph.segments.begin(), graph.segments.end(),
                                      [&](const contract::BoundarySegment& s) { return s.id == id; });
    ASSERT_NE(segment, graph.segments.end());
    EXPECT_EQ(segment->boundary_type, contract::BoundaryType::Convergent);
    std::string members;
    for (grid::CellId c = 0; c < f.n; ++c)
        if (graph.boundary_segment_id[c] == id)
            members += " " + std::to_string(c) + "(c" + std::to_string(f.convergence[c]) + " d" +
                       std::to_string(f.divergence[c]) + " s" + std::to_string(f.shear[c]) + " r" +
                       std::to_string(f.regime[c]) + ")";
    // Both columns of the face edge and, at the two cube corners, the face 4 and face 5 corner
    // cells, whose oblique contact with face 0's corner cells converges too.
    EXPECT_EQ(segment->cell_count, 2U * N + 2U) << members;
    EXPECT_EQ(graph.boundary_segment_id[f.table->neighbour(f.face0(N - 1, 0), grid::north)], id) << members;
    EXPECT_EQ(graph.boundary_segment_id[f.table->neighbour(f.face0(N - 1, N - 1), grid::south)], id) << members;
    EXPECT_NEAR(segment->centroid_x, 45.0f, 1e-3f);
    EXPECT_NEAR(segment->centroid_y, 0.0f, 1e-3f);
}

TEST(BoundaryMotion, CubeDeformingRegionSpansAFaceEdge) {
    const int N = 4;
    const auto table = grid::TectonicGrid::cube(N, 6371.0);
    const size_t n = table->cellcount();
    std::vector<uint32_t> owner(n, 1U), nearest(n, contract::kNoBoundaryId);
    std::vector<uint8_t> crust(n, uint8_t(CrustClass::Continental)), divergence(n, 100), zero(n, 0);
    std::vector<uint8_t> type(n, uint8_t(contract::BoundaryType::None)), regime(n, 0);
    std::vector<uint16_t> distance(n, 9);
    for (grid::CellId c = 0; c < N * N; ++c) owner[c] = 0U;
    for (int j = 0; j < N; ++j) {
        const grid::CellId west = static_cast<grid::CellId>(j * N + N - 1);
        for (const grid::CellId c : {west, table->neighbour(west, grid::east)}) {
            nearest[c] = 5U;
            distance[c] = 0;
        }
    }
    contract::BoundarySegment segment;
    segment.id = 5U;
    segment.left_plate_id = 0U;
    segment.right_plate_id = 1U;
    segment.boundary_type = contract::BoundaryType::Divergent;
    segment.geologic_regime = GeologicRegime::DivergentRift;
    segment.average_divergence_score = 100;
    std::array<contract::PlateKinematics, 2> plates{};
    const double spacing = std::sqrt(table->totalarea_m2() / static_cast<double>(n));
    // Plate 0 moves west at one cell per update on the equator; plate 1 is at rest.
    const std::array<mesh::Vec3, 2> rotation{mesh::Vec3{0.0, 0.0, -spacing / 6371000.0}, mesh::Vec3{}};
    const deforming_regions::Inputs inputs{0, 0, owner.data(), crust.data(), zero.data(),
        divergence.data(), zero.data(), type.data(), regime.data(), distance.data(), nearest.data(),
        &segment, 1U, plates.data(), 2U, table.get(), rotation.data()};
    const auto outputs = deforming_regions::build(inputs);
    ASSERT_EQ(outputs.regions.size(), 1U);
    EXPECT_EQ(outputs.regions.front().cell_count, 2U * N);
    EXPECT_NEAR(outputs.regions.front().centroid_x, 45.0f, 1e-3f);
    EXPECT_NEAR(outputs.regions.front().centroid_y, 0.0f, 1e-3f);
    // A face-0 cell (the owner is the left plate: weight 0.68) moves west, never north.
    const grid::CellId cell = static_cast<grid::CellId>(1 * N + N - 1);
    const double expected = -0.68 * std::cos(table->mesh()->latitude_rad(cell));
    EXPECT_EQ(outputs.deforming_region_id[cell], 5U);
    EXPECT_NEAR(outputs.deformation_velocity_x[cell], expected, 1e-4);
    EXPECT_NEAR(outputs.deformation_velocity_y[cell], 0.0f, 1e-5f);
    // The face-1 cell across the edge: weight 0.32 of plate 0.
    const grid::CellId across = table->neighbour(cell, grid::east);
    EXPECT_EQ(outputs.deforming_region_id[across], 5U);
    EXPECT_NEAR(outputs.deformation_velocity_x[across],
                -0.32 * std::cos(table->mesh()->latitude_rad(across)), 1e-4);
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
