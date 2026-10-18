// CubeKernels.*: the landscape kernels on a cube CellGrid through
// landscape::SphereMetric (docs/architecture/sphere-mesh.md, section 3; T8 step
// W2-L), at N 16 and 32: drainage of a bowl, Priority-Flood and flow accumulation
// on a seeded random field, the hillslope step visiting every face once, and the
// eroded volume of a radially symmetric valley centred in a face interior, on a
// face edge and on a cube corner (ValleyErosionIsRotationallyInvariant: what D8
// allows).

#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <utility>
#include <vector>
#include <gtest/gtest.h>

namespace
{
using landscape::SphereMetric;

constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
constexpr int facesizes[] = { 16, 32 };

SphereMetric cubemetric(int facecells) { return SphereMetric(mesh::makecubegrid(facecells, radius)); }

double angle(const mesh::Vec3& p, const mesh::Vec3& q) { return std::atan2(mesh::length(mesh::cross(p, q)), mesh::dot(p, q)); }

// Follows receivers to the cell that is its own receiver; -1 on a chain longer than the grid.
std::int32_t chainend(const landscape::FlowGraph& graph, std::size_t cell)
{
    std::size_t current = cell;
    for (std::size_t step = 0; step <= graph.size(); step++)
    {
        const std::size_t next = static_cast<std::size_t>(graph.receiver[current]);
        if (next == current)
            return static_cast<std::int32_t>(current);
        current = next;
    }
    return -1;
}

bool isneighbour(const SphereMetric& metric, std::size_t cell, std::size_t other)
{
    for (int k = 1; k <= landscape::directioncount; k++)
        if (metric.neighbourid(cell, k) == static_cast<std::int32_t>(other))
            return true;
    return false;
}

// Uniform terrain in [0, 1000) m from seed 1; cells below `seafraction` x 1000 m are sea.
struct Field
{
    std::vector<float> elevation;
    std::vector<std::uint8_t> sea;
};

Field randomfield(std::size_t count, double seafraction)
{
    std::mt19937 random(1);
    std::uniform_real_distribution<double> uniform(0.0, 1000.0);
    Field field;
    field.elevation.resize(count);
    field.sea.resize(count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const double value = uniform(random);
        field.elevation[cell] = static_cast<float>(value);
        field.sea[cell] = value < seafraction * 1000.0 ? 1 : 0;
    }
    return field;
}
}

TEST(CubeKernels, BowlDrainsToItsBottom)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const mesh::CellGrid& grid = *metric.grid;
        const std::size_t count = metric.size();
        // Bottoms in a face interior and in a cube corner cell.
        const mesh::CellId bottoms[] = { grid.locate(mesh::unitvector(0.1, 0.2)),
                                         grid.locate(mesh::normalise({ 1.0, 1.0, 1.0 })) };
        for (mesh::CellId bottom : bottoms)
        {
            const mesh::Vec3 centre = grid.centre(bottom);
            std::vector<float> elevation(count);
            for (std::size_t cell = 0; cell < count; cell++)
                elevation[cell] = static_cast<float>(1000.0 * angle(centre, grid.centre(static_cast<mesh::CellId>(cell))));
            const std::vector<std::uint8_t> sea(count, 0);
            const landscape::Depressions depressions = landscape::filldepressions(metric, elevation, sea, 0.0f);
            landscape::FlowGraph graph = landscape::buildflowgraph(metric, elevation, sea, 0.0f, depressions);
            landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);
            for (std::size_t cell = 0; cell < count; cell++)
            {
                ASSERT_EQ(chainend(graph, cell), static_cast<std::int32_t>(bottom)) << "N " << n << " cell " << cell;
                if (cell == bottom)
                    continue;
                const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
                ASSERT_NE(graph.direction[cell], 0) << "N " << n << " cell " << cell;
                ASSERT_EQ(metric.neighbourid(cell, graph.direction[cell]), graph.receiver[cell]);
                ASSERT_LT(elevation[receiver], elevation[cell]) << "N " << n << " cell " << cell;
            }
            EXPECT_NEAR(graph.area_m2[bottom], metric.totalarea_m2(), 1.0e-9 * metric.totalarea_m2()) << "N " << n;
            EXPECT_NEAR(graph.discharge_m3_yr[bottom], metric.totalarea_m2(), 1.0e-9 * metric.totalarea_m2());
        }
    }
}

TEST(CubeKernels, PriorityFloodLeavesNoPit)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.03);
        const landscape::Depressions depressions = landscape::filldepressions(metric, field.elevation, field.sea, 0.0f);
        const std::vector<float>& filled = depressions.filledelevation;
        std::size_t raised = 0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (field.sea[cell] != 0)
            {
                ASSERT_EQ(filled[cell], field.elevation[cell]);
                continue;
            }
            ASSERT_GE(filled[cell], field.elevation[cell]);
            raised += filled[cell] > field.elevation[cell] ? 1 : 0;
            // No pit: the flood receiver is a neighbour no higher on the filled surface,
            // and its chain reaches the sea.
            const std::int32_t receiver = depressions.floodreceiver[cell];
            ASSERT_GE(receiver, 0) << "N " << n << " cell " << cell;
            ASSERT_TRUE(isneighbour(metric, cell, static_cast<std::size_t>(receiver))) << "N " << n << " cell " << cell;
            ASSERT_LE(filled[static_cast<std::size_t>(receiver)], filled[cell]) << "N " << n << " cell " << cell;
            std::size_t current = cell;
            std::size_t steps = 0;
            while (field.sea[current] == 0 && steps++ <= count)
                current = static_cast<std::size_t>(depressions.floodreceiver[current]);
            ASSERT_NE(field.sea[current], 0) << "N " << n << " cell " << cell;
        }
        EXPECT_GT(raised, 0u) << "the random field has depressions to fill";
        // The routing on the filled surface has no cycle and every land cell reaches an outlet.
        const landscape::FlowGraph graph =
            landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions);
        for (std::size_t cell = 0; cell < count; cell++)
            if (field.sea[cell] == 0)
            {
                const std::int32_t end = chainend(graph, cell);
                ASSERT_GE(end, 0);
                ASSERT_GE(graph.outletsea[static_cast<std::size_t>(end)], 0) << "N " << n << " cell " << cell;
            }
    }
}

TEST(CubeKernels, FlowAccumulationSumsToTotalArea)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        for (double seafraction : { 0.03, 0.0 })
        {
            const Field field = randomfield(count, seafraction);
            const landscape::Depressions depressions =
                landscape::filldepressions(metric, field.elevation, field.sea, 0.0f);
            landscape::FlowGraph graph = landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions);
            landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);
            double land = 0.0, roots = 0.0;
            for (std::size_t cell = 0; cell < count; cell++)
            {
                if (graph.land[cell] == 0)
                    continue;
                land += metric.area_m2(cell);
                if (graph.receiver[cell] == static_cast<std::int32_t>(cell))
                    roots += graph.area_m2[cell];
                ASSERT_NEAR(graph.discharge_m3_yr[cell], graph.area_m2[cell], 1.0e-12 * graph.area_m2[cell]);
            }
            EXPECT_NEAR(roots, land, 1.0e-12 * land) << "N " << n << " sea fraction " << seafraction;
            if (seafraction == 0.0)
                EXPECT_NEAR(roots, metric.totalarea_m2(), 1.0e-12 * metric.totalarea_m2()) << "N " << n;
        }
    }
}

// Hillslope diffusion alone moves, across every face shared by two land cells, the
// limited flux computed once: the eroded volume is the sum over the distinct faces.
TEST(CubeKernels, HillslopeVisitsEveryFaceOnce)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.0);
        const landscape::Depressions depressions = landscape::filldepressions(metric, field.elevation, field.sea, 0.0f);
        landscape::FlowGraph graph = landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions);
        landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);

        landscape::StreamPowerParameters law;
        law.depositioncoefficient = 0.0;
        law.hillslopediffusivity_km2_myr = 1.0e3;
        std::vector<float> elevation = field.elevation, erosion(count), deposition(count);
        const landscape::StreamPowerBudget budget = landscape::stepstreampower(
            metric, graph, std::vector<float>(count, 0.0f), std::vector<float>(count, 0.0f), law, elevation, erosion,
            deposition);

        std::map<std::pair<std::size_t, std::size_t>, double> faces;
        const double diffusivity = law.hillslopediffusivity_km2_myr * 1.0e6;
        for (std::size_t i = 0; i < count; i++)
            for (int k : { 1, 3, 5, 7 })
            {
                const std::int32_t other = metric.neighbourid(i, k);
                ASSERT_GE(other, 0) << "every cardinal slot of a cube cell has a neighbour";
                const std::size_t j = static_cast<std::size_t>(other);
                ASSERT_EQ(metric.neighbourid(j, metric.reverse(i, k)), static_cast<std::int32_t>(i));
                const double zi = field.elevation[i], zj = field.elevation[j];
                if (zi == zj)
                    continue;
                const double drop = std::fabs(zi - zj);
                const double limit = 0.05 * drop * std::min(metric.area_m2(i), metric.area_m2(j));
                const double flux = diffusivity * law.step_myr * drop / metric.distance_m(i, k) * metric.facelength_m(i, k);
                faces.emplace(std::minmax(i, j), std::min(limit, flux));
            }
        EXPECT_EQ(faces.size(), 2 * count) << "a cube has two faces per cell";
        double expected = 0.0;
        for (const auto& face : faces)
            expected += face.second;
        EXPECT_NEAR(budget.eroded_km3 * 1.0e9, expected, 1.0e-9 * expected) << "N " << n;
        EXPECT_NEAR(budget.landdeposited_km3, budget.eroded_km3, 1.0e-12 * budget.eroded_km3);
    }
}

namespace
{
struct ValleyLaw
{
    double seaspacings = 1.5;  // the sea's radius in nominal cell spacings
    float erodibility = 100.0f;
    double areaexponent = 0.5;
    int steps = 3;
};

// The eroded volume (km^3) of a radially symmetric valley over the whole sphere
// centred at `centre`: elevation rising 2000 m per pi of great-circle distance from
// a sea of radius law.seaspacings nominal cell spacings (dy_m), eroded by
// detachment-limited steps without hillslope diffusion, rerouted each step.
double valleyvolume(const SphereMetric& metric, const mesh::Vec3& centre, const ValleyLaw& law)
{
    const mesh::CellGrid& grid = *metric.grid;
    const std::size_t count = metric.size();
    const double seaangle = law.seaspacings * metric.dy_m / metric.radius_m;
    std::vector<float> elevation(count);
    std::vector<std::uint8_t> sea(count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const double distance = angle(centre, grid.centre(static_cast<mesh::CellId>(cell)));
        sea[cell] = distance < seaangle ? 1 : 0;
        elevation[cell] = sea[cell] != 0 ? -100.0f : static_cast<float>(2000.0 * (distance - seaangle) / pi);
    }
    landscape::StreamPowerParameters parameters;
    parameters.depositioncoefficient = 0.0;
    parameters.hillslopediffusivity_km2_myr = 0.0;
    parameters.areaexponent = law.areaexponent;
    const std::vector<float> erodibility(count, law.erodibility), uplift(count, 0.0f), runoff(count, 1.0f);
    std::vector<float> erosion(count), deposition(count);
    double volume = 0.0;
    for (int step = 0; step < law.steps; step++)
    {
        const landscape::Depressions depressions = landscape::filldepressions(metric, elevation, sea, 0.0f);
        landscape::FlowGraph graph = landscape::buildflowgraph(metric, elevation, sea, 0.0f, depressions);
        landscape::accumulateflow(metric, runoff, graph);
        volume += landscape::stepstreampower(metric, graph, erodibility, uplift, parameters, elevation, erosion,
                                             deposition)
                      .eroded_km3;
    }
    return volume;
}
}

// Rotational invariance of the valley's eroded volume. D8 itself is anisotropic: the
// steepest of eight links follows a radial gradient only to within its angle to the
// nearest link, so a valley's volume depends on how its rays meet the grid (on the
// latitude-longitude grid, whose kernels are the reference, a valley centred on the
// equator erodes 3 % less than one centred on a pole at m = 0, 1 to 9 % less at
// m = 0.5; measured at W2-L). Two checks:
// - placements related by the cube's symmetries (face centres on an equatorial and a
//   polar face, midpoints of an equatorial and two polar face edges, three corners)
//   erode the same volume: within 1e-6 at m = 0 (the same links, through rotated
//   frames), within 1 % at m = 0.5 (D8 ties between symmetric links break by the
//   lowest code, so the networks differ a little);
// - a face interior, a face edge and a corner, the brief's comparison, agree within
//   2.5 % at m = 0, the anisotropy of D8 on the cube (1 to 2 % measured, less than
//   the latitude-longitude grid's).
TEST(CubeKernels, ValleyErosionIsRotationallyInvariant)
{
    const mesh::Vec3 placements[3][3] = {
        { { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, -1.0, 0.0 } },
        { mesh::normalise({ 1.0, 1.0, 0.0 }), mesh::normalise({ 1.0, 0.0, 1.0 }), mesh::normalise({ 0.0, 1.0, 1.0 }) },
        { mesh::normalise({ 1.0, 1.0, 1.0 }), mesh::normalise({ 1.0, -1.0, -1.0 }), mesh::normalise({ -1.0, -1.0, 1.0 }) },
    };
    const char* kinds[3] = { "face interior", "face edge", "cube corner" };
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        for (double m : { 0.5, 0.0 })
        {
            ValleyLaw law;
            law.areaexponent = m;
            law.erodibility = m == 0.0 ? 1000.0f : 100.0f;
            const double tolerance = m == 0.0 ? 1.0e-6 : 0.01;
            double first[3] = { 0.0, 0.0, 0.0 };
            for (int kind = 0; kind < 3; kind++)
                for (int placement = 0; placement < 3; placement++)
                {
                    const double volume = valleyvolume(metric, placements[kind][placement], law);
                    ASSERT_GT(volume, 0.0);
                    if (placement == 0)
                        first[kind] = volume;
                    else
                        EXPECT_NEAR(volume, first[kind], tolerance * first[kind])
                            << "N " << n << " m " << m << ": " << kinds[kind] << " placement " << placement;
                }
            std::cout << "N " << n << " m " << m << " eroded km3: face " << first[0] << ", edge " << first[1]
                      << ", corner " << first[2] << "\n";
            if (m == 0.0)
                for (int kind = 1; kind < 3; kind++)
                    EXPECT_NEAR(first[kind], first[0], 0.025 * first[0]) << "N " << n << ": " << kinds[kind];
        }
    }
}
