// The engine's cell table on the cube (W2-E): steps across face edges, the edge table, rings,
// discs, Dijkstra, quantiles, centroids and the plate state the lithosphere keeps on it; the
// disabled run probe of the per-update measures and timing (E3 checklist item 9).
#include "gtest/gtest.h"

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_recipe.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <set>
#include <string>
#include <utility>

using tectonics::grid::CellId;
using tectonics::grid::TectonicGrid;
using namespace tectonics::grid;

namespace {

std::string environment(const char* name) {
#ifdef _MSC_VER
    char* value = nullptr;
    std::size_t size = 0;
    if (_dupenv_s(&value, &size, name) != 0 || value == nullptr) return {};
    std::string result(value);
    std::free(value);
    return result;
#else
    const char* value = std::getenv(name);
    return value ? value : "";
#endif
}

} // namespace

TEST(TectonicGrid, CubeStepsAreReciprocalAcrossFaceEdges) {
    const auto grid = TectonicGrid::cube(6, 6371.0);
    ASSERT_EQ(grid->cellcount(), 6u * 36u);
    for (CellId c = 0; c < grid->cellcount(); ++c)
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            const Step step = grid->ahead(c, k);
            ASSERT_NE(step.cell, nocell);
            // Turning round at the next cell leads back.
            EXPECT_EQ(grid->neighbour(step.cell, opposite(step.direction)), c) << c << " " << int(k);
        }
}

TEST(TectonicGrid, CubeStraightWalksCloseGreatCircles) {
    const int n = 5;
    const auto grid = TectonicGrid::cube(n, 6371.0);
    for (CellId start : {CellId(0), CellId(2 * n * n + 3), CellId(4 * n * n + n + 1)})
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            CellId c = start;
            std::uint8_t d = k;
            for (int i = 0; i < 4 * n; ++i) {
                const Step step = grid->ahead(c, d);
                c = step.cell;
                d = step.direction;
            }
            EXPECT_EQ(c, start) << start << " " << int(k);
            EXPECT_EQ(d, k);
        }
}

// Every adjacent pair once; the normal is the shared edge's perpendicular (the normalised cross
// product of its two vertices, tangent at the edge's midpoint, from a toward b), the centre
// line the tangent of the line between the centres; both agree within 30 degrees and exactly
// on the axes of a face.
TEST(TectonicGrid, CubeEdgeTableHasEveryPairOnceWithTheSharedEdgesPerpendicular) {
    const int n = 4;
    const auto grid = TectonicGrid::cube(n, 6371.0);
    const auto& edges = grid->edges();
    ASSERT_EQ(edges.size(), 2 * grid->cellcount());
    std::set<std::pair<CellId, CellId>> pairs;
    std::vector<int> degree(grid->cellcount(), 0);
    double worst = 1;
    for (const Edge& edge : edges) {
        EXPECT_LT(edge.a, edge.b);
        EXPECT_TRUE(pairs.insert({edge.a, edge.b}).second);
        EXPECT_EQ(grid->neighbour(edge.a, edge.direction), edge.b);
        EXPECT_EQ(grid->edges()[grid->edge(edge.a, edge.direction)].b, edge.b);
        ++degree[edge.a];
        ++degree[edge.b];
        const mesh::Vec3 pa = grid->position(edge.a), pb = grid->position(edge.b);
        EXPECT_NEAR(mesh::length(edge.normal), 1.0, 1e-12);
        EXPECT_NEAR(mesh::length(edge.midpoint), 1.0, 1e-12);
        EXPECT_NEAR(mesh::dot(edge.normal, edge.midpoint), 0.0, 1e-12);
        EXPECT_NEAR(mesh::dot(edge.centreline, mesh::normalise(pa + pb)), 0.0, 1e-12);
        EXPECT_GT(mesh::dot(edge.normal, pb - pa), 0.0);
        worst = std::min(worst, mesh::dot(edge.normal, edge.centreline));
        // The midpoint lies between the two centres.
        EXPECT_LT(mesh::greatcircle_rad(edge.midpoint, pa), mesh::greatcircle_rad(pa, pb));
        EXPECT_LT(mesh::greatcircle_rad(edge.midpoint, pb), mesh::greatcircle_rad(pa, pb));
        EXPECT_GT(edge.length_m, 0.0);
        EXPECT_GT(edge.distance_m, 0.0);
        EXPECT_GE(edge.length_m, grid->smallestedge_m(edge.a) - 1e-6);
    }
    for (int d : degree) EXPECT_EQ(d, 4);
    RecordProperty("worst_normal_centreline_cosine", std::to_string(worst));
    EXPECT_GT(worst, std::cos(30.0 * 3.14159265358979323846 / 180.0));
    // On a face's central row (odd N) normal and centre line agree.
    const auto odd = TectonicGrid::cube(5, 6371.0);
    const Edge& middle = odd->edges()[odd->edge(12, east)];
    EXPECT_NEAR(mesh::dot(middle.normal, middle.centreline), 1.0, 1e-12);
}

TEST(TectonicGrid, CubeQuantileCentroidOffsetsAndDistances) {
    const int n = 8;
    const auto grid = TectonicGrid::cube(n, 6371.0);
    std::vector<float> values(grid->cellcount());
    for (CellId c = 0; c < grid->cellcount(); ++c) values[c] = static_cast<float>(c);
    EXPECT_EQ(grid->quantile(values, 0.0f), 0.0f);
    EXPECT_EQ(grid->quantile(values, 1.0f), static_cast<float>(grid->cellcount() - 1));
    // Half the area lies below the median: faces 0-2 of six equal faces.
    const float median = grid->quantile(values, 0.5f);
    EXPECT_NEAR(median, 3.0f * n * n, 1.0f);

    std::vector<std::uint32_t> owners(grid->cellcount(), 1);
    for (CellId c = 0; c < static_cast<CellId>(n * n); ++c) owners[c] = 0;  // face 0, centred on +X
    const auto centroids = grid->centroids(owners, 2);
    EXPECT_NEAR(centroids[0].x, 1.0, 1e-12);
    EXPECT_NEAR(centroids[0].y, 0.0, 1e-12);
    EXPECT_NEAR(centroids[0].z, 0.0, 1e-12);

    const mesh::TangentFrame frame = mesh::tangentframe(centroids[0], 6371000.0);
    for (CellId c : {CellId(0), CellId(n * n / 2 + 3), CellId(n * n - 1)}) {
        const auto offset = grid->tangentoffset_m(centroids[0], c);
        const mesh::Vec3 back = mesh::fromtangent(frame, offset[0] / 1000.0, offset[1] / 1000.0);
        EXPECT_NEAR(mesh::length(back - grid->position(c)), 0.0, 1e-9) << c;
    }
    EXPECT_TRUE(std::isnan(grid->tangentoffset_m(centroids[0], 2 * n * n + 3)[0]));  // face 2, far side

    const CellId c = 3 * n + 4;
    EXPECT_EQ(grid->withinkm(c, 1500.0), mesh::cellswithin(*grid->mesh(), grid->position(c), 1500000.0));
    const Distances distances = grid->dijkstra({c});
    for (CellId other = 0; other < grid->cellcount(); ++other)
        EXPECT_GE(distances.distance_m[other] + 1e-6,
                  mesh::greatcircle_m(grid->position(c), grid->position(other), 6371000.0));
    // Rings of cardinal steps, ascending ids; a mask and a limit stop Dijkstra.
    const auto rings = grid->rings(c, 2);
    ASSERT_EQ(rings.size(), 3u);
    EXPECT_EQ(rings[1].size(), 4u);
    EXPECT_EQ(rings[2].size(), 8u);
    EXPECT_TRUE(std::is_sorted(rings[2].begin(), rings[2].end()));
    std::vector<std::uint8_t> passable(grid->cellcount(), 1);
    passable[grid->neighbour(c, east)] = 0;
    const Distances masked = grid->dijkstra({c}, passable);
    EXPECT_TRUE(std::isinf(masked.distance_m[grid->neighbour(c, east)]));
    const Distances limited = grid->dijkstra({c}, {}, 1.0);
    EXPECT_TRUE(std::isinf(limited.distance_m[grid->neighbour(c, east)]));
}

// A plate around a great circle (the equatorial band of cells) has a vanishing centre sum: its
// centroid falls back to its member cell nearest the sum, a unit vector of the plate.
TEST(TectonicGrid, CentroidOfAnEquatorialRingFallsBackToAMemberCell) {
    const int n = 8;
    const auto grid = TectonicGrid::cube(n, 6371.0);
    std::vector<std::uint32_t> owners(grid->cellcount(), 1);
    std::vector<mesh::Vec3> members;
    for (CellId c = 0; c < grid->cellcount(); ++c) {
        const mesh::Vec3 p = grid->position(c);
        if (std::abs(p.z) < .1) { owners[c] = 0; members.push_back(p); }
    }
    ASSERT_GT(members.size(), 0u);
    const auto centroids = grid->centroids(owners, 2);
    EXPECT_NEAR(mesh::length(centroids[0]), 1.0, 1e-12);
    EXPECT_TRUE(std::any_of(members.begin(), members.end(), [&](const mesh::Vec3& m) {
        return mesh::length(m - centroids[0]) < 1e-12; }));
    // The rest of the sphere, symmetric about the equator, keeps its sum's direction rule.
    EXPECT_NEAR(mesh::length(centroids[1]), 1.0, 1e-12);
    // The helper itself: an exactly zero sum takes the first member.
    EXPECT_EQ(tectonics::boundary_cube::meandirection({}, 4.0, members).x, members.front().x);
    EXPECT_EQ(mesh::length(tectonics::boundary_cube::meandirection({}, 0.0, {})), 0.0);
}

TEST(TectonicGrid, LithosphereKeepsPlateStateOnTheCube) {
    tectonics::scenario::Scenario scenario;
    scenario.seed = 1;
    scenario.sphere_face_cells = 8;
    scenario.plate_count = 5;
    scenario.delta_time_myr = 1.0;
    lithosphere world(scenario);
    for (int i = 0; i < 3 && !world.isFinished(); ++i) world.update();
    const TectonicGrid& grid = world.getGrid();
    EXPECT_EQ(grid.facecells(), 8);
    EXPECT_EQ(world.getWidth(), 8u);
    EXPECT_EQ(world.getHeight(), 48u);
    double area = 0;
    for (std::uint32_t p = 0; p < world.getPlateCount(); ++p) {
        const tectonics::PlateState& state = world.getPlate(p);
        const auto& published = world.getPlateKinematics()[p];
        EXPECT_NEAR(mesh::length(state.centroid), 1.0, 1e-12);
        const auto centre = tectonics::boundary_cube::encodecentroid(state.centroid);
        EXPECT_EQ(published.centroid_longitude_deg, centre[0]);
        EXPECT_EQ(published.centroid_latitude_deg, centre[1]);
        EXPECT_EQ(published.rotation_rad_s[0], state.rotation_rad_s.x);
        EXPECT_EQ(published.rotation_rad_s[1], state.rotation_rad_s.y);
        EXPECT_EQ(published.rotation_rad_s[2], state.rotation_rad_s.z);
        // The centroid's velocity omega x R c in cm/yr.
        const mesh::Vec3 v = mesh::cross(state.rotation_rad_s, grid.radius_m() * state.centroid);
        const double cm_yr = 100.0 * tectonics::dynamics::kSecondsPerMyr / 1.0e6;
        EXPECT_NEAR(std::hypot(published.velocity_east_cm_yr, published.velocity_north_cm_yr),
                    mesh::length(v) * cm_yr, 1e-5 * (1.0 + mesh::length(v) * cm_yr));
        EXPECT_EQ(state.stable_id, world.getTopologyDiagnostics().stable_ids[p]);
        EXPECT_GE(state.continental_share, 0.0f);
        EXPECT_LE(state.continental_share, 1.0f);
        EXPECT_TRUE(std::isfinite(state.rotation_rad_s.x));
        area += state.area_m2;
    }
    EXPECT_NEAR(area, grid.totalarea_m2(), 1e-6 * grid.totalarea_m2());
}

// The engine's per-update measures and wall time on one run (W2-E, E3 checklist item 9):
// UW_TECTONICS_MEASURES the JSON path; UW_TECTONICS_N, default 64; UW_TECTONICS_MYR, default 50.
// Seed 1, the baseline. Per update: the largest displacement over each cell's smallest edge,
// the transport's substeps, |1 - coverage| over rigid interiors by cell class, the force
// balance's condition, spin shares, net rotation, area-RMS and centroid speeds, the column
// mass residual and the update's seconds; then their extremes, a histogram and the phases.
TEST(TectonicMeasures, DISABLED_Run) {
    const std::string path = environment("UW_TECTONICS_MEASURES");
    ASSERT_FALSE(path.empty()) << "set UW_TECTONICS_MEASURES to the output path";
    const std::string n_text = environment("UW_TECTONICS_N");
    const std::string myr_text = environment("UW_TECTONICS_MYR");
    auto scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
    scenario.sphere_face_cells = !n_text.empty() ? static_cast<std::uint32_t>(std::stoul(n_text)) : 64u;
    scenario.evolution_duration_myr = !myr_text.empty() ? std::stod(myr_text) : 50.0;
    const auto start = std::chrono::steady_clock::now();
    lithosphere world(scenario);
    const auto created = std::chrono::steady_clock::now();
    const tectonics::PhaseTimings at_creation = world.getPhaseTimings();
    std::ofstream out(path);
    ASSERT_TRUE(out.good()) << path;
    out.precision(6);
    out << "{\n  \"grid\": \"cube N " << scenario.sphere_face_cells << "\",\n"
        << "  \"seed\": 1,\n  \"scenario\": \"tectonics baseline (configs/tectonics-baseline.json)\",\n"
        << "  \"duration_myr\": " << scenario.evolution_duration_myr << ",\n"
        << "  \"delta_time_myr\": " << scenario.delta_time_myr << ",\n"
        << "  \"cells\": " << world.getGrid().cellcount() << ",\n"
        << "  \"cell_spacing_km\": " << world.getGrid().cellspacing_m() / 1000 << ",\n"
        << "  \"updates\": [\n";
    double largest_ratio = 0, largest_residual = 0, largest_interior = 0, largest_edge = 0, largest_corner = 0;
    double largest_condition = 0, rms_sum = 0, largest_rms = 0, largest_net = 0, update_sum = 0, mean_share_sum = 0;
    std::array<int, 11> histogram{};
    int updates = 0;
    while (!world.isFinished()) {
        const auto before = std::chrono::steady_clock::now();
        world.update();
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - before).count();
        const auto& t = world.getTransportDiagnostics();
        const auto& m = world.getMotionDiagnostics();
        const double residual = std::abs(world.getColumnModel().budget().relative_mass_residual);
        largest_ratio = std::max(largest_ratio, t.max_displacement_edges);
        largest_residual = std::max(largest_residual, residual);
        largest_interior = std::max(largest_interior, t.coverage_interior_max);
        largest_edge = std::max(largest_edge, t.coverage_edge_max);
        largest_corner = std::max(largest_corner, t.coverage_corner_max);
        largest_condition = std::max(largest_condition, m.condition_number);
        largest_rms = std::max(largest_rms, m.rms_speed_cm_yr);
        largest_net = std::max(largest_net, m.net_rotation_cm_yr);
        rms_sum += m.rms_speed_cm_yr;
        mean_share_sum += m.spin_share_mean;
        update_sum += seconds;
        ++histogram[std::min<std::size_t>(10, static_cast<std::size_t>(t.max_displacement_edges * 10))];
        out << (updates ? ",\n" : "") << "    {\"time_myr\": " << world.getTimeMyr()
            << ", \"plates\": " << world.getPlateCount()
            << ", \"max_displacement_over_edge\": " << t.max_displacement_edges
            << ", \"substeps\": " << t.substeps
            << ", \"coverage_interior_max\": " << t.coverage_interior_max
            << ", \"coverage_interior_mean\": " << t.coverage_interior_mean
            << ", \"coverage_edge_max\": " << t.coverage_edge_max
            << ", \"coverage_edge_mean\": " << t.coverage_edge_mean
            << ", \"coverage_corner_max\": " << t.coverage_corner_max
            << ", \"coverage_corner_mean\": " << t.coverage_corner_mean
            << ", \"rigid_cells\": [" << t.interior_cells << ", " << t.edge_cells << ", " << t.corner_cells << "]"
            << ", \"condition\": " << m.condition_number
            << ", \"spin_share_mean\": " << m.spin_share_mean
            << ", \"spin_share_max\": " << m.spin_share_max
            << ", \"net_rotation_cm_yr\": " << m.net_rotation_cm_yr
            << ", \"rms_speed_cm_yr\": " << m.rms_speed_cm_yr
            << ", \"centroid_speed_mean_cm_yr\": " << m.mean_speed_cm_yr
            << ", \"centroid_speed_max_cm_yr\": " << m.max_speed_cm_yr
            << ", \"force_residual\": " << m.relative_force_residual
            << ", \"mass_residual\": " << residual
            << ", \"seconds\": " << seconds << "}";
        ++updates;
    }
    const auto end = std::chrono::steady_clock::now();
    const double creation = std::chrono::duration<double>(created - start).count();
    const tectonics::PhaseTimings& timings = world.getPhaseTimings();
    out << "\n  ],\n  \"summary\": {\n"
        << "    \"updates\": " << updates << ",\n"
        << "    \"plates_at_end\": " << world.getPlateCount() << ",\n"
        << "    \"max_displacement_over_edge\": " << largest_ratio << ",\n"
        << "    \"displacement_histogram_tenths\": [";
    for (std::size_t k = 0; k < histogram.size(); ++k) out << (k ? ", " : "") << histogram[k];
    out << "],\n"
        << "    \"coverage_interior_max\": " << largest_interior << ",\n"
        << "    \"coverage_edge_max\": " << largest_edge << ",\n"
        << "    \"coverage_corner_max\": " << largest_corner << ",\n"
        << "    \"condition_max\": " << largest_condition << ",\n"
        << "    \"spin_share_mean\": " << (updates ? mean_share_sum / updates : 0) << ",\n"
        << "    \"net_rotation_max_cm_yr\": " << largest_net << ",\n"
        << "    \"rms_speed_mean_cm_yr\": " << (updates ? rms_sum / updates : 0) << ",\n"
        << "    \"rms_speed_max_cm_yr\": " << largest_rms << ",\n"
        << "    \"rms_speed_earth_cm_yr\": 0.937,\n"
        << "    \"mass_residual_max\": " << largest_residual << ",\n"
        << "    \"creation_s\": " << creation << ",\n"
        << "    \"updates_s\": " << std::chrono::duration<double>(end - created).count() << ",\n"
        << "    \"seconds_per_update\": " << (updates ? update_sum / updates : 0) << "\n"
        << "  },\n  \"update_phases_s\": {\n";
    for (int phase = 0; phase < tectonics::PhaseTimings::count; ++phase)
        out << "    \"" << tectonics::PhaseTimings::name(phase) << "\": " << timings.seconds[phase] - at_creation.seconds[phase]
            << (phase + 1 < tectonics::PhaseTimings::count ? ",\n" : "\n");
    out << "  }\n}\n";
    EXPECT_LT(largest_residual, 1e-9);
}
