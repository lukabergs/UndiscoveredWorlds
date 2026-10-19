// CubeKernels.*: the landscape kernels on a cube CellGrid through
// landscape::SphereMetric (docs/architecture/sphere-mesh.md, section 3; T8 steps
// W2-L and its follow-ups), at N 16 and 32 (the hillslope bump at 32 and 64): drainage
// of a bowl, Priority-Flood, flow accumulation, randomised routing, relief roughness,
// the ocean mask, distributaries and sub-grid statistics on seeded fields; the
// hillslope step (rock conserved, a bump spreading alike at a face centre, edge and
// corner with the non-orthogonal correction); the eroded volume of a radially
// symmetric valley at those places (ValleyErosionIsRotationallyInvariant: what D8
// allows).

#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/distributaries.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/ocean_mask.hpp"
#include "simulations/landscape/relief_roughness.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"
#include "simulations/landscape/subgrid_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
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

// Hillslope diffusion alone on a random field: every cardinal slot has a neighbour
// that reaches the cell back through its reverse slot (two faces per cell), the step
// conserves the rock (sum of area x change 0) and its budget closes.
TEST(CubeKernels, HillslopeConservesRock)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.0);
        const landscape::Depressions depressions = landscape::filldepressions(metric, field.elevation, field.sea, 0.0f);
        landscape::FlowGraph graph = landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions);
        landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);

        std::size_t faces = 0;
        for (std::size_t i = 0; i < count; i++)
            for (int k : { 1, 3, 5, 7 })
            {
                const std::int32_t other = metric.neighbourid(i, k);
                ASSERT_GE(other, 0) << "every cardinal slot of a cube cell has a neighbour";
                const std::size_t j = static_cast<std::size_t>(other);
                ASSERT_EQ(metric.neighbourid(j, metric.reverse(i, k)), static_cast<std::int32_t>(i));
                faces += i < j ? 1 : 0;
            }
        EXPECT_EQ(faces, 2 * count) << "a cube has two faces per cell";

        landscape::StreamPowerParameters law;
        law.depositioncoefficient = 0.0;
        law.hillslopediffusivity_km2_myr = 1.0e3;
        std::vector<float> elevation = field.elevation, erosion(count), deposition(count);
        const landscape::StreamPowerBudget budget = landscape::stepstreampower(
            metric, graph, std::vector<float>(count, 0.0f), std::vector<float>(count, 0.0f), law, elevation, erosion,
            deposition);
        double change = 0.0, moved = 0.0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            change += (static_cast<double>(deposition[cell]) - static_cast<double>(erosion[cell])) * metric.area_m2(cell);
            moved += static_cast<double>(erosion[cell]) * metric.area_m2(cell);
        }
        ASSERT_GT(budget.eroded_km3, 0.0);
        EXPECT_NEAR(moved, budget.eroded_km3 * 1.0e9, 1.0e-6 * moved);
        EXPECT_NEAR(change, 0.0, 1.0e-6 * moved) << "N " << n;
        EXPECT_NEAR(budget.landdeposited_km3, budget.eroded_km3, 1.0e-12 * budget.eroded_km3);
    }
}

namespace
{
// The growth of a Gaussian bump's spread under hillslope diffusion alone: the bump
// (500 m, sigma 1200 km) on a 1000 m plain covering the sphere, centred at `centre`,
// diffused `steps` times with kappa dt = 2e8 m^2; returns the change of the
// excess-mass-weighted mean squared great-circle distance from the centre (m^2),
// 4 kappa t on a plane for any grid that diffuses isotropically.
double bumpspread(const SphereMetric& metric, const mesh::Vec3& centre, int steps)
{
    const mesh::CellGrid& grid = *metric.grid;
    const std::size_t count = metric.size();
    const double sigma = 1.2e6;
    std::vector<float> elevation(count);
    std::vector<double> distance2(count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const double s = radius * angle(centre, grid.centre(static_cast<mesh::CellId>(cell)));
        distance2[cell] = s * s;
        elevation[cell] = static_cast<float>(1000.0 + 500.0 * std::exp(-0.5 * s * s / (sigma * sigma)));
    }
    const auto moment = [&]() {
        double mass = 0.0, weighted = 0.0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            const double excess = (static_cast<double>(elevation[cell]) - 1000.0) * metric.area_m2(cell);
            mass += excess;
            weighted += excess * distance2[cell];
        }
        return weighted / mass;
    };
    const std::vector<std::uint8_t> sea(count, 0);
    const landscape::Depressions depressions = landscape::filldepressions(metric, elevation, sea, 0.0f);
    landscape::FlowGraph graph = landscape::buildflowgraph(metric, elevation, sea, 0.0f, depressions);
    landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);
    landscape::StreamPowerParameters law;
    law.depositioncoefficient = 0.0;
    law.hillslopediffusivity_km2_myr = 2.0e8 / 1.0e6;
    const std::vector<float> none(count, 0.0f);
    std::vector<float> erosion(count), deposition(count);
    const double before = moment();
    for (int step = 0; step < steps; step++)
        landscape::stepstreampower(metric, graph, none, none, law, elevation, erosion, deposition);
    return moment() - before;
}
}

// A radially symmetric bump diffused by the hillslope term spreads alike whether it
// sits at a face centre, on a face edge or on a cube corner (the non-orthogonal
// correction of stream_power.cpp; without it the two-point flux spreads it faster on
// an edge and slower at a corner): within 2 % at N 32 and 64.
TEST(CubeKernels, HillslopeBumpIsRotationallyInvariant)
{
    const mesh::Vec3 centres[3] = { { 1.0, 0.0, 0.0 }, mesh::normalise({ 1.0, 1.0, 0.0 }),
                                    mesh::normalise({ 1.0, 1.0, 1.0 }) };
    const char* kinds[3] = { "face centre", "face edge", "cube corner" };
    for (int n : { 32, 64 })
    {
        const SphereMetric metric = cubemetric(n);
        double spread[3];
        for (int c = 0; c < 3; c++)
            spread[c] = bumpspread(metric, centres[c], 100);
        std::cout << "N " << n << " bump spread m^2: face " << spread[0] << ", edge " << spread[1] << ", corner "
                  << spread[2] << " (4 kappa t " << 4.0 * 2.0e8 * 100 << ")\n";
        ASSERT_GT(spread[0], 0.0);
        for (int c = 1; c < 3; c++)
            EXPECT_NEAR(spread[c], spread[0], 0.02 * spread[0]) << "N " << n << ": " << kinds[c];
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
//   2.5 % at m = 0 and 8 % at m = 0.5, the anisotropy of D8 on the cube (1 to 2 % and
//   3 to 7 % measured, less than the latitude-longitude grid's).
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
            const double spread = m == 0.0 ? 0.025 : 0.08;
            for (int kind = 1; kind < 3; kind++)
                EXPECT_NEAR(first[kind], first[0], spread * first[0]) << "N " << n << " m " << m << ": " << kinds[kind];
        }
    }
}

// Relief roughness samples the noise at the cube cells' centres: sea cells keep
// their elevation, land changes by at most about the amplitude, and the change is as
// smooth across face edges as inside the faces (no seam).
TEST(CubeKernels, ReliefRoughnessAtCellCentres)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        Field field = randomfield(count, 0.03);
        for (std::size_t cell = 0; cell < count; cell++)
            if (field.sea[cell] == 0)
                field.elevation[cell] += 1000.0f;  // above the coastal taper
        landscape::ReliefRoughness roughness;
        roughness.amplitude_m = 50.0;
        roughness.wavelength_km = 4000.0;
        roughness.octaves = 4;
        roughness.seed = 1;
        std::vector<float> elevation = field.elevation;
        landscape::addreliefroughness(metric, roughness, field.sea, 0.0f, elevation);
        std::vector<double> delta(count);
        std::size_t changed = 0, land = 0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            delta[cell] = static_cast<double>(elevation[cell]) - static_cast<double>(field.elevation[cell]);
            if (field.sea[cell] != 0)
            {
                ASSERT_EQ(delta[cell], 0.0);
                continue;
            }
            land++;
            changed += delta[cell] != 0.0 ? 1 : 0;
            ASSERT_LE(std::fabs(delta[cell]), 2.0 * roughness.amplitude_m) << "N " << n << " cell " << cell;
        }
        EXPECT_GT(changed, land * 9 / 10);
        // Mean |change difference| per metre over cardinal links inside a face and across a face edge.
        double inside = 0.0, across = 0.0;
        std::size_t insidelinks = 0, acrosslinks = 0;
        for (std::size_t cell = 0; cell < count; cell++)
            for (int k : { 1, 3, 5, 7 })
            {
                const std::size_t other = static_cast<std::size_t>(metric.neighbourid(cell, k));
                if (field.sea[cell] != 0 || field.sea[other] != 0)
                    continue;
                const double rate = std::fabs(delta[other] - delta[cell]) / metric.distance_m(cell, k);
                const bool sameface = cell / static_cast<std::size_t>(n * n) == other / static_cast<std::size_t>(n * n);
                (sameface ? inside : across) += rate;
                (sameface ? insidelinks : acrosslinks)++;
            }
        ASSERT_GT(acrosslinks, 0u);
        const double ratio = (across / acrosslinks) / (inside / insidelinks);
        EXPECT_GT(ratio, 0.75) << "N " << n;
        EXPECT_LT(ratio, 1.33) << "N " << n;
    }
}

// Randomised receivers on the cube (the bearing offsets of tangentoffset): the graph
// has no cycle, every receiver is a neighbour lower on the filled surface (or the sea
// of an outlet), the choice differs from D8 somewhere, and accumulation still sums
// to the land area.
TEST(CubeKernels, RandomisedRoutingOnTheCube)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.03);
        const landscape::Depressions depressions = landscape::filldepressions(metric, field.elevation, field.sea, 0.0f);
        const landscape::FlowGraph steepest =
            landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions);
        landscape::RoutingOptions options;
        options.randomised = true;
        options.seed = 1;
        landscape::FlowGraph graph =
            landscape::buildflowgraph(metric, field.elevation, field.sea, 0.0f, depressions, options);
        landscape::accumulateflow(metric, std::vector<float>(count, 1.0f), graph);
        const std::vector<float>& filled = depressions.filledelevation;
        std::size_t differs = 0;
        double land = 0.0, roots = 0.0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (field.sea[cell] != 0)
                continue;
            land += metric.area_m2(cell);
            ASSERT_NE(graph.direction[cell], 0) << "N " << n << " cell " << cell;
            const std::int32_t target = metric.neighbourid(cell, graph.direction[cell]);
            ASSERT_GE(target, 0);
            if (graph.outletsea[cell] >= 0)
                ASSERT_EQ(target, graph.outletsea[cell]);
            else
            {
                ASSERT_EQ(target, graph.receiver[cell]);
                if (depressions.depressionid[cell] == 0)
                    ASSERT_LE(filled[static_cast<std::size_t>(target)], filled[cell]) << "N " << n << " cell " << cell;
            }
            differs += graph.direction[cell] != steepest.direction[cell] ? 1 : 0;
            if (graph.receiver[cell] == static_cast<std::int32_t>(cell))
                roots += graph.area_m2[cell];
        }
        EXPECT_GT(differs, 0u) << "N " << n;
        EXPECT_NEAR(roots, land, 1.0e-12 * land) << "N " << n;
    }
}

// The ocean mask joins a below-sea cap across the three faces of a cube corner into
// one ocean and leaves isolated below-sea cells out.
TEST(CubeKernels, OceanMaskAcrossACorner)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const mesh::CellGrid& grid = *metric.grid;
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.03);
        const mesh::Vec3 corner = mesh::normalise({ 1.0, 1.0, 1.0 });
        const double cap = 1.5e6 / radius;
        std::vector<std::uint8_t> belowsea(count, 0), incap(count, 0);
        for (std::size_t cell = 0; cell < count; cell++)
        {
            incap[cell] = angle(corner, grid.centre(static_cast<mesh::CellId>(cell))) < cap ? 1 : 0;
            belowsea[cell] = incap[cell] != 0 || field.sea[cell] != 0 ? 1 : 0;
        }
        const std::vector<std::uint8_t> ocean = landscape::oceanmask(metric, belowsea, 0.01);
        std::size_t faces[6] = { 0, 0, 0, 0, 0, 0 }, isolated = 0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (belowsea[cell] == 0)
            {
                ASSERT_EQ(ocean[cell], 0);
                continue;
            }
            if (incap[cell] != 0)
            {
                ASSERT_EQ(ocean[cell], 1) << "N " << n << " cell " << cell;
                faces[cell / static_cast<std::size_t>(n * n)]++;
                continue;
            }
            bool alone = true;
            for (int k = 1; k <= landscape::directioncount; k++)
            {
                const std::int32_t other = metric.neighbourid(cell, k);
                alone = alone && (other < 0 || belowsea[static_cast<std::size_t>(other)] == 0);
            }
            if (alone)
            {
                isolated++;
                ASSERT_EQ(ocean[cell], 0) << "N " << n << " cell " << cell;
            }
        }
        EXPECT_GT(faces[0], 0u);
        EXPECT_GT(faces[1], 0u);
        EXPECT_GT(faces[4], 0u);
        EXPECT_GT(isolated, 0u);
    }
}

// Distributaries on a delta ring around a bay centred on a face edge: the flow that
// enters the prograded land all reaches the sea (nothing stranded on a surface that
// falls towards the bay), and every delta cell's code points to a neighbour.
TEST(CubeKernels, DistributariesAroundAFaceEdgeBay)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const mesh::CellGrid& grid = *metric.grid;
        const std::size_t count = metric.size();
        const mesh::Vec3 centre = mesh::normalise({ 1.0, 1.0, 0.0 });
        const double spacing = metric.dy_m / radius;
        const double bay = 4.0 * spacing, delta = 7.0 * spacing;
        std::vector<std::uint8_t> sea(count), prograded(count);
        std::vector<float> surface(count);
        for (std::size_t cell = 0; cell < count; cell++)
        {
            const double distance = angle(centre, grid.centre(static_cast<mesh::CellId>(cell)));
            sea[cell] = distance < bay ? 1 : 0;
            prograded[cell] = sea[cell] == 0 && distance < delta ? 1 : 0;
            surface[cell] = sea[cell] != 0 ? -50.0f : static_cast<float>(5.0 + 100.0 * (distance - bay) / spacing);
        }
        const landscape::Depressions depressions = landscape::filldepressions(metric, surface, sea, 0.0f);
        const landscape::FlowGraph graph = landscape::buildflowgraph(metric, surface, sea, 0.0f, depressions);
        const std::vector<std::int32_t> january(count, 100), july(count, 50);
        const landscape::Distributaries result = landscape::splitdistributaries(
            metric, sea, prograded, depressions.filledelevation, 0.0, graph.direction, january, july, 1.0);
        EXPECT_GT(result.inflowjanuary, 0.0) << "N " << n;
        EXPECT_NEAR(result.seajanuary, result.inflowjanuary, 1.0e-9 * result.inflowjanuary);
        EXPECT_NEAR(result.seajuly, result.inflowjuly, 1.0e-9 * result.inflowjuly);
        EXPECT_EQ(result.strandedjanuary, 0.0);
        std::size_t cells = 0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (result.direction[cell] == 0)
                continue;
            ASSERT_NE(prograded[cell], 0);
            ASSERT_GE(metric.neighbourid(cell, result.direction[cell]), 0) << "N " << n << " cell " << cell;
            cells++;
        }
        EXPECT_GT(cells, 0u) << "N " << n;
    }
}

// Sub-grid statistics on the random field: zero at sea, the relief within the cell's
// local relief (its largest difference to the seven or eight neighbours, corner cells
// included) and the spacing at most the cell size.
TEST(CubeKernels, SubgridStatisticsOnTheCube)
{
    for (int n : facesizes)
    {
        const SphereMetric metric = cubemetric(n);
        const std::size_t count = metric.size();
        const Field field = randomfield(count, 0.03);
        landscape::SubgridParameters parameters;
        parameters.diffusivity_km2_myr = 200.0;
        parameters.referencedoverk = 1.0e4;
        parameters.referencerunoff_mm_yr = 500.0;
        const std::vector<float> erodibility(count, 1.0f), runoff(count, 500.0f), rate(count, 50.0f);
        std::vector<float> relief, spacing;
        landscape::subgridstatistics(metric, parameters, field.elevation, field.sea, erodibility, runoff, rate, relief,
                                     spacing);
        std::size_t positive = 0;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (field.sea[cell] != 0)
            {
                ASSERT_EQ(relief[cell], 0.0f);
                ASSERT_EQ(spacing[cell], 0.0f);
                continue;
            }
            float local = 0.0f;
            for (int k = 1; k <= landscape::directioncount; k++)
            {
                const std::int32_t other = metric.neighbourid(cell, k);
                if (other >= 0)
                    local = std::max(local, std::fabs(field.elevation[static_cast<std::size_t>(other)] - field.elevation[cell]));
            }
            ASSERT_GE(relief[cell], 0.0f);
            ASSERT_LE(relief[cell], local) << "N " << n << " cell " << cell;
            ASSERT_LE(spacing[cell], static_cast<float>(std::sqrt(metric.area_m2(cell)) / 1000.0) * (1.0f + 1.0e-6f));
            positive += relief[cell] > 0.0f ? 1 : 0;
        }
        EXPECT_GT(positive, 0u) << "N " << n;
    }
}
