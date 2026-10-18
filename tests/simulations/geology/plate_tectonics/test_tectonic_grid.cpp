// The engine's cell table (W2-E, step E1): the planar scaffold that reproduces the torus,
// the cube over mesh::CellGrid, and the plate state the lithosphere keeps on it.
#include "gtest/gtest.h"

#include "core/mesh/mesh_geometry.hpp"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
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

CellId planarid(const TectonicGrid& grid, std::uint32_t x, std::uint32_t y) { return y * grid.width() + x; }

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

TEST(TectonicGrid, PlanarNeighboursWrapBothSeams) {
    const auto grid = TectonicGrid::planar(10, 6, 2000.0);
    EXPECT_EQ(grid->cellcount(), 60u);
    EXPECT_DOUBLE_EQ(grid->cell_m(), 2000.0 * 1000.0 / 10);
    // Across the x seam (the torus' wrapped rectangles of old).
    EXPECT_EQ(grid->neighbour(planarid(*grid, 0, 3), west), planarid(*grid, 9, 3));
    EXPECT_EQ(grid->neighbour(planarid(*grid, 9, 3), east), planarid(*grid, 0, 3));
    // Across the y seam.
    EXPECT_EQ(grid->neighbour(planarid(*grid, 4, 0), north), planarid(*grid, 4, 5));
    EXPECT_EQ(grid->neighbour(planarid(*grid, 4, 5), south), planarid(*grid, 4, 0));
    for (CellId c = 0; c < grid->cellcount(); ++c)
        for (std::uint8_t k = 0; k < directioncount; ++k) {
            const Step step = grid->ahead(c, k);
            EXPECT_EQ(step.cell, grid->neighbour(c, k));
            EXPECT_EQ(step.direction, k);
            EXPECT_EQ(grid->neighbour(step.cell, opposite(k)), c);
        }
}

TEST(TectonicGrid, PlanarEdgesAreTheTorusFaces) {
    const auto grid = TectonicGrid::planar(8, 4, 800.0);
    const auto& edges = grid->edges();
    ASSERT_EQ(edges.size(), 2 * grid->cellcount());
    for (CellId c = 0; c < grid->cellcount(); ++c) {
        const Edge& x = edges[2 * c];
        const Edge& y = edges[2 * c + 1];
        EXPECT_EQ(x.a, c);
        EXPECT_EQ(x.b, grid->neighbour(c, east));
        EXPECT_EQ(x.axis, 0);
        EXPECT_EQ(x.normal.x, 1.0);
        EXPECT_EQ(x.normal.y, 0.0);
        EXPECT_EQ(y.b, grid->neighbour(c, south));
        EXPECT_EQ(y.axis, 1);
        EXPECT_EQ(y.normal.y, 1.0);
        EXPECT_EQ(x.length_m, grid->cell_m());
        EXPECT_EQ(y.distance_m, grid->cell_m());
    }
    EXPECT_EQ(grid->area_m2(5), grid->cell_m() * grid->cell_m());
}

TEST(TectonicGrid, PlanarCentroidWrapsTheSeam) {
    // A plate across the x seam (columns 9, 0, 1) and across the y seam (rows 5, 0).
    const auto grid = TectonicGrid::planar(10, 6, 1000.0);
    std::vector<std::uint32_t> owners(grid->cellcount(), 1);
    for (std::uint32_t y : {5u, 0u})
        for (std::uint32_t x : {9u, 0u, 1u}) owners[planarid(*grid, x, y)] = 0;
    const auto centroids = grid->centroids(owners, 2);
    const double dx = std::remainder(centroids[0].x - 0.0, 10.0);
    const double dy = std::remainder(centroids[0].y - 5.5, 6.0);
    EXPECT_NEAR(dx, 0.0, 1e-9);
    EXPECT_NEAR(dy, 0.0, 1e-9);
    const auto offset = grid->tangentoffset_m(centroids[0], planarid(*grid, 9, 5));
    EXPECT_NEAR(offset[0], -grid->cell_m(), 1e-6);
    EXPECT_NEAR(offset[1], -0.5 * grid->cell_m(), 1e-6);
}

TEST(TectonicGrid, PlanarQuantileIsTheSortedIndex) {
    const auto grid = TectonicGrid::planar(5, 5, 100.0);
    std::vector<float> values(25);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = static_cast<float>((i * 7) % 25);
    for (const float q : {0.0f, 0.1f, 0.37f, 0.5f, 0.99f, 1.0f}) {
        std::vector<float> sorted = values;
        std::sort(sorted.begin(), sorted.end());
        const auto index = static_cast<std::size_t>(std::floor(q * static_cast<float>(sorted.size() - 1)));
        EXPECT_EQ(grid->quantile(values, q), sorted[index]) << q;
    }
}

TEST(TectonicGrid, PlanarRingsWithinKmAndDijkstra) {
    const auto grid = TectonicGrid::planar(12, 8, 1200.0);  // 100 km cells
    const CellId c = planarid(*grid, 0, 0);
    const auto rings = grid->rings(c, 2);
    ASSERT_EQ(rings.size(), 3u);
    EXPECT_EQ(rings[0].size(), 1u);
    EXPECT_EQ(rings[1].size(), 4u);
    EXPECT_EQ(rings[2].size(), 8u);
    EXPECT_TRUE(std::is_sorted(rings[2].begin(), rings[2].end()));
    // Centres within 150 km: the cell, its four cardinal and four diagonal neighbours (141 km).
    EXPECT_EQ(grid->withinkm(c, 150.0).size(), 9u);
    EXPECT_EQ(grid->withinkm(c, 99.0).size(), 1u);
    const Distances distances = grid->dijkstra({c, planarid(*grid, 6, 4)});
    EXPECT_EQ(distances.distance_m[planarid(*grid, 11, 7)], 2 * grid->cell_m());
    EXPECT_EQ(distances.source[planarid(*grid, 11, 7)], 0u);
    EXPECT_EQ(distances.source[planarid(*grid, 6, 5)], 1u);
    // A mask stops paths; a limit stops them too.
    std::vector<std::uint8_t> passable(grid->cellcount(), 1);
    passable[grid->neighbour(c, east)] = 0;
    const Distances masked = grid->dijkstra({c}, passable);
    EXPECT_TRUE(std::isinf(masked.distance_m[grid->neighbour(c, east)]));
    EXPECT_EQ(masked.distance_m[planarid(*grid, 2, 0)], 4 * grid->cell_m());  // round the masked cell
    const Distances limited = grid->dijkstra({c}, passable, 3 * grid->cell_m());
    EXPECT_TRUE(std::isinf(limited.distance_m[planarid(*grid, 2, 0)]));
    EXPECT_EQ(limited.distance_m[planarid(*grid, 1, 1)], 2 * grid->cell_m());
}

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

TEST(TectonicGrid, CubeEdgeTableHasEveryPairOnce) {
    const auto grid = TectonicGrid::cube(4, 6371.0);
    const auto& edges = grid->edges();
    ASSERT_EQ(edges.size(), 2 * grid->cellcount());
    std::set<std::pair<CellId, CellId>> pairs;
    std::vector<int> degree(grid->cellcount(), 0);
    for (const Edge& edge : edges) {
        EXPECT_LT(edge.a, edge.b);
        EXPECT_TRUE(pairs.insert({edge.a, edge.b}).second);
        EXPECT_EQ(grid->neighbour(edge.a, edge.direction), edge.b);
        ++degree[edge.a];
        ++degree[edge.b];
        const mesh::Vec3 pa = grid->position(edge.a), pb = grid->position(edge.b);
        const mesh::Vec3 middle = mesh::normalise(pa + pb);
        EXPECT_NEAR(mesh::length(edge.normal), 1.0, 1e-12);
        EXPECT_NEAR(mesh::dot(edge.normal, middle), 0.0, 1e-12);
        EXPECT_GT(mesh::dot(edge.normal, pb - pa), 0.0);
        EXPECT_GT(edge.length_m, 0.0);
        EXPECT_GT(edge.distance_m, 0.0);
    }
    for (int d : degree) EXPECT_EQ(d, 4);
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
}

TEST(TectonicGrid, LithosphereKeepsPlateStateOnTheScaffold) {
    tectonics::scenario::Scenario scenario;
    scenario.seed = 1;
    scenario.width = 48;
    scenario.height = 24;
    scenario.plate_count = 5;
    scenario.delta_time_myr = 1.0;
    lithosphere world(scenario);
    for (int i = 0; i < 3 && !world.isFinished(); ++i) world.update();
    const TectonicGrid& grid = world.getGrid();
    EXPECT_EQ(grid.kind(), TectonicGrid::Kind::planar);
    EXPECT_EQ(grid.width(), 48u);
    double area = 0;
    for (std::uint32_t p = 0; p < world.getPlateCount(); ++p) {
        const tectonics::PlateState& state = world.getPlate(p);
        const auto& published = world.getPlateKinematics()[p];
        EXPECT_EQ(published.linear_velocity_x, state.velocityX());
        EXPECT_EQ(published.unit_velocity_y, state.direction_y);
        EXPECT_EQ(published.angular_velocity, 0.0f);
        EXPECT_EQ(published.mass_center_x, static_cast<float>(state.centroid.x));
        EXPECT_EQ(state.stable_id, world.getTopologyDiagnostics().stable_ids[p]);
        EXPECT_GE(state.continental_share, 0.0f);
        EXPECT_LE(state.continental_share, 1.0f);
        area += state.area_m2;
    }
    EXPECT_NEAR(area, grid.totalarea_m2(), 1e-6 * grid.totalarea_m2());
}

// The engine's per-phase wall time on one run (UW_TECTONICS_TIMING: the JSON path;
// UW_TECTONICS_TIMING_WIDTH, default 256; UW_TECTONICS_TIMING_MYR, default 100). Seed 1, baseline.
TEST(TectonicTiming, DISABLED_Run) {
    const std::string path = environment("UW_TECTONICS_TIMING");
    ASSERT_FALSE(path.empty()) << "set UW_TECTONICS_TIMING to the output path";
    const std::string width_text = environment("UW_TECTONICS_TIMING_WIDTH");
    const std::string myr_text = environment("UW_TECTONICS_TIMING_MYR");
    auto scenario = tectonics::scenario::baseline();
    scenario.seed = 1;
    scenario.width = !width_text.empty() ? static_cast<std::uint32_t>(std::stoul(width_text)) : 256u;
    scenario.height = scenario.width / 2;
    scenario.evolution_duration_myr = !myr_text.empty() ? std::stod(myr_text) : 100.0;
    const auto start = std::chrono::steady_clock::now();
    lithosphere world(scenario);
    const auto created = std::chrono::steady_clock::now();
    const tectonics::PhaseTimings at_creation = world.getPhaseTimings();
    while (!world.isFinished()) world.update();
    const auto end = std::chrono::steady_clock::now();
    const double creation = std::chrono::duration<double>(created - start).count();
    const double updates = std::chrono::duration<double>(end - created).count();
    const tectonics::PhaseTimings& timings = world.getPhaseTimings();
    std::ofstream out(path);
    ASSERT_TRUE(out.good()) << path;
    out.precision(6);
    out << "{\n  \"grid\": \"planar scaffold " << scenario.width << " x " << scenario.height << "\",\n"
        << "  \"seed\": 1,\n  \"scenario\": \"tectonics baseline (configs/tectonics-baseline.json)\",\n"
        << "  \"duration_myr\": " << scenario.evolution_duration_myr << ",\n"
        << "  \"delta_time_myr\": " << scenario.delta_time_myr << ",\n"
        << "  \"updates\": " << timings.updates << ",\n"
        << "  \"plates_at_end\": " << world.getPlateCount() << ",\n"
        << "  \"total_s\": " << creation + updates << ",\n"
        << "  \"creation_s\": " << creation << ",\n"
        << "  \"updates_s\": " << updates << ",\n"
        << "  \"update_phases_s\": {\n";
    double attributed = 0;
    for (int phase = 0; phase < tectonics::PhaseTimings::count; ++phase) {
        const double seconds = timings.seconds[phase] - at_creation.seconds[phase];
        attributed += seconds;
        out << "    \"" << tectonics::PhaseTimings::name(phase) << "\": " << seconds
            << (phase + 1 < tectonics::PhaseTimings::count ? ",\n" : "\n");
    }
    out << "  },\n  \"update_unattributed_s\": " << updates - attributed << ",\n"
        << "  \"creation_phases_s\": {\n";
    for (int phase = 0; phase < tectonics::PhaseTimings::count; ++phase)
        out << "    \"" << tectonics::PhaseTimings::name(phase) << "\": " << at_creation.seconds[phase]
            << (phase + 1 < tectonics::PhaseTimings::count ? ",\n" : "\n");
    out << "  }\n}\n";
}
