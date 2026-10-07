// Landscape.TorusMetric*: the landscape kernels on a planar torus of square cells
// (torusmetric; docs/architecture/erosion-hydrology-plan.md, section 4, Phase D,
// step D2a): the metric's wrapping, a flow graph whose drainage crosses both
// wraps, and one stream power step on it.

#include "landscape_fixture.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/relief_roughness.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int width = 16;
constexpr int height = 8;
constexpr double cell_m = 156250.0;            // the tectonics engine's cell
constexpr double cellarea_m2 = 2.44140625e10;  // cell_m^2, exact in double
constexpr double runoff_m_yr = 0.5;
constexpr int landcells = width * height - 1;  // 127

// One sea cell at (0, 0), 100 m below the sea level 0; every land cell 100 m per
// cell of torus distance to it (the shortest over the wraps, centre to centre), so
// every land cell has a lower neighbour towards the sea cell and the eight
// neighbours of the sea cell, five of them across a wrap, are its outlets.
struct Torus
{
    landscape::SphereMetric metric = landscape::torusmetric(width, height, cell_m);
    std::size_t seacell = 0;
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    landscape::FlowGraph graph;
};

Torus maketorus(const landscape::RoutingOptions& options = {})
{
    Torus torus;
    const landscape::SphereMetric& metric = torus.metric;
    torus.seacell = metric.index(0, 0);
    torus.elevation.assign(metric.size(), 0.0f);
    torus.sea.assign(metric.size(), 0);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const int dx = std::min(x, width - x), dy = std::min(y, height - y);
            torus.elevation[metric.index(x, y)] = static_cast<float>(100.0 * std::hypot(dx, dy));
        }
    torus.sea[torus.seacell] = 1;
    torus.elevation[torus.seacell] = -100.0f;
    const landscape::Depressions depressions = landscape::filldepressions(metric, torus.elevation, torus.sea, 0.1f);
    torus.graph = landscape::buildflowgraph(metric, torus.elevation, torus.sea, 0.0f, depressions, options);
    landscape::accumulateflow(metric, std::vector<float>(metric.size(), static_cast<float>(runoff_m_yr)), torus.graph);
    return torus;
}

landscape::StreamPowerParameters makeparameters(double deposition, double diffusivity_km2_myr)
{
    landscape::StreamPowerParameters parameters;
    parameters.step_myr = 1.0;
    parameters.baselevel_m = 0.0;
    parameters.areaexponent = 0.5;
    parameters.slopeexponent = 1.0;
    parameters.depositioncoefficient = deposition;
    parameters.hillslopediffusivity_km2_myr = diffusivity_km2_myr;
    parameters.criticalslope = 0.0;
    parameters.sedimenttravel_km = 300.0;
    return parameters;
}

landscape::StreamPowerBudget step(Torus& torus, const landscape::StreamPowerParameters& parameters,
                                  float erodibility_k = 24.0f)
{
    const std::size_t cells = torus.metric.size();
    const std::vector<float> erodibility(cells, erodibility_k), uplift(cells, 0.0f);
    std::vector<float> erosion(cells, 0.0f), deposition(cells, 0.0f);
    return landscape::stepstreampower(torus.metric, torus.graph, erodibility, uplift, parameters, torus.elevation,
                                      erosion, deposition);
}
}

TEST_F(Landscape, TorusMetricWrapsRowsAndColumns)
{
    const landscape::SphereMetric metric = landscape::torusmetric(width, height, cell_m);
    EXPECT_TRUE(metric.wrapsrows);
    EXPECT_EQ(metric.radius_m, 0.0);
    EXPECT_EQ(metric.dy_m, cell_m);
    EXPECT_EQ(metric.area_m2(metric.index(5, 3)), cellarea_m2);
    EXPECT_NEAR(metric.totalarea_m2(), 3.125e12, 1e-9 * 3.125e12);  // 16 x 8 x cell^2

    const double diagonal = cell_m * std::sqrt(2.0);  // 220970.869... m
    for (int y : { 0, 3, height - 1 })
    {
        const std::size_t cell = metric.index(7, y);
        for (int direction : { 1, 3, 5, 7 })
        {
            EXPECT_EQ(metric.distance_m(cell, direction), cell_m) << "row " << y << " direction " << direction;
            EXPECT_EQ(metric.facelength_m(cell, direction), cell_m) << "row " << y << " direction " << direction;
        }
        for (int direction : { 2, 4, 6, 8 })
            EXPECT_NEAR(metric.distance_m(cell, direction), diagonal, 1e-9 * diagonal)
                << "row " << y << " direction " << direction;
    }

    EXPECT_EQ(metric.neighbour(metric.index(3, 0), 1), metric.index(3, height - 1));   // north of row 0
    EXPECT_EQ(metric.neighbour(metric.index(3, height - 1), 5), metric.index(3, 0));   // south of row H - 1
    EXPECT_EQ(metric.neighbour(metric.index(width - 1, 2), 3), metric.index(0, 2));    // east of column W - 1
    EXPECT_EQ(metric.neighbour(metric.index(0, 2), 7), metric.index(width - 1, 2));    // west of column 0
    EXPECT_EQ(metric.neighbour(metric.index(width - 1, 0), 2), metric.index(0, height - 1));  // north-east corner
    EXPECT_EQ(metric.neighbour(metric.index(0, height - 1), 6), metric.index(width - 1, 0));  // south-west corner

    // Relief roughness is sphere noise: a torus metric is refused.
    landscape::ReliefRoughness roughness;
    roughness.amplitude_m = 30.0;
    roughness.wavelength_km = 2000.0;
    std::vector<float> elevation(metric.size(), 100.0f);
    EXPECT_THROW(landscape::addreliefroughness(metric, roughness, std::vector<uint8_t>(metric.size(), 0), 0.0f, elevation),
                 std::invalid_argument);

    // The sphere is unchanged: a neighbour beyond a pole is the cell itself.
    const landscape::SphereMetric sphere(width, height, 6371000.0);
    EXPECT_FALSE(sphere.wrapsrows);
    EXPECT_EQ(sphere.neighbour(sphere.index(3, 0), 1), sphere.index(3, 0));
    EXPECT_EQ(sphere.neighbour(sphere.index(3, height - 1), 5), sphere.index(3, height - 1));
}

namespace
{
// The drainage invariants of a torus graph: the order lists every land cell once,
// each after its receiver; every receiver chain ends at an outlet into the sea
// cell; discharge = runoff x area everywhere; the outlets' areas sum to the land
// area and their discharges to runoff x land area. Returns the outlet count.
int expectdrainstosea(const Torus& torus)
{
    const landscape::SphereMetric& metric = torus.metric;
    const landscape::FlowGraph& graph = torus.graph;
    EXPECT_EQ(graph.size(), metric.size());
    EXPECT_EQ(graph.order.size(), static_cast<std::size_t>(landcells));
    if (graph.size() != metric.size())
        return -1;

    // The order lists every land cell once, each after its receiver.
    std::vector<int> position(metric.size(), -1);
    for (std::size_t k = 0; k < graph.order.size(); k++)
    {
        const std::size_t cell = static_cast<std::size_t>(graph.order[k]);
        EXPECT_EQ(position[cell], -1) << "cell " << cell;
        position[cell] = static_cast<int>(k);
    }

    double outletarea = 0.0, outletdischarge = 0.0;
    int outlets = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (cell == torus.seacell)
        {
            EXPECT_EQ(static_cast<int>(graph.land[cell]), 0);
            EXPECT_EQ(graph.area_m2[cell], 0.0);
            continue;
        }
        EXPECT_GE(position[cell], 0) << "cell " << cell;
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver != cell)
            EXPECT_LT(position[receiver], position[cell]) << "cell " << cell;
        // The receiver chain ends at an outlet into the sea cell.
        std::size_t walk = cell;
        for (int count = 0; count <= landcells && static_cast<std::size_t>(graph.receiver[walk]) != walk; count++)
            walk = static_cast<std::size_t>(graph.receiver[walk]);
        EXPECT_EQ(static_cast<std::size_t>(graph.receiver[walk]), walk) << "cell " << cell;
        EXPECT_EQ(graph.outletsea[walk], static_cast<int32_t>(torus.seacell)) << "cell " << cell;
        if (receiver == cell)
        {
            outlets++;
            outletarea += graph.area_m2[cell];
            outletdischarge += graph.discharge_m3_yr[cell];
        }
        EXPECT_NEAR(graph.discharge_m3_yr[cell], runoff_m_yr * graph.area_m2[cell], 1e-12 * graph.area_m2[cell]);
    }
    const double landarea = landcells * cellarea_m2;  // 3.1005859375e12 m^2
    EXPECT_NEAR(outletarea, landarea, 1e-12 * landarea);
    EXPECT_NEAR(outletdischarge, runoff_m_yr * landarea, 1e-12 * landarea);  // 1.55029296875e12 m^3/yr
    return outlets;
}
}

TEST_F(Landscape, TorusMetricFlowGraphDrainsAcrossTheWraps)
{
    const Torus torus = maketorus();
    const landscape::SphereMetric& metric = torus.metric;
    const landscape::FlowGraph& graph = torus.graph;
    EXPECT_EQ(expectdrainstosea(torus), 8);  // the eight neighbours of the sea cell

    // The randomised routing (its facet offsets use the neighbour's wrapped row):
    // the same invariants.
    landscape::RoutingOptions randomised;
    randomised.randomised = true;
    randomised.seed = 1;
    const Torus drawn = maketorus(randomised);
    EXPECT_EQ(expectdrainstosea(drawn), 8);
    int redirected = 0;  // the draw does pick other receivers than D8
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        redirected += drawn.graph.receiver[cell] != graph.receiver[cell] ? 1 : 0;
    EXPECT_GT(redirected, 0);

    // Outlets across the wraps: (0, 7) drains south and (15, 7) south-east into (0, 0).
    const std::size_t lastrow = metric.index(0, height - 1), corner = metric.index(width - 1, height - 1);
    EXPECT_EQ(graph.outletsea[lastrow], static_cast<int32_t>(torus.seacell));
    EXPECT_EQ(static_cast<int>(graph.direction[lastrow]), 5);
    EXPECT_FLOAT_EQ(graph.distance_m[lastrow], static_cast<float>(cell_m));
    EXPECT_EQ(graph.outletsea[corner], static_cast<int32_t>(torus.seacell));
    EXPECT_EQ(static_cast<int>(graph.direction[corner]), 4);
    EXPECT_FLOAT_EQ(graph.distance_m[corner], static_cast<float>(cell_m * std::sqrt(2.0)));
}

TEST_F(Landscape, TorusMetricStreamPowerStepClosesItsBudget)
{
    // Detachment-limited, no diffusion, no uplift: no land cell rises and all the
    // eroded rock is laid on the sea cell.
    {
        Torus torus = maketorus();
        const std::vector<float> before = torus.elevation;
        const landscape::StreamPowerBudget budget = step(torus, makeparameters(0.0, 0.0));
        for (std::size_t cell = 0; cell < torus.metric.size(); cell++)
        {
            ASSERT_TRUE(std::isfinite(torus.elevation[cell])) << "cell " << cell;
            if (cell != torus.seacell)
                EXPECT_LE(torus.elevation[cell], before[cell]) << "cell " << cell;
        }
        EXPECT_GT(torus.elevation[torus.seacell], before[torus.seacell]);
        EXPECT_GT(budget.eroded_km3, 0.0);
        EXPECT_EQ(budget.landdeposited_km3, 0.0);
        EXPECT_NEAR(budget.eroded_km3, budget.marinedeposited_km3, 1e-9 * budget.eroded_km3);
    }
    // Deposition and hillslope diffusion across the wrapped faces: still closed.
    {
        Torus torus = maketorus();
        const landscape::StreamPowerBudget budget = step(torus, makeparameters(1.0, 200.0));
        for (std::size_t cell = 0; cell < torus.metric.size(); cell++)
            ASSERT_TRUE(std::isfinite(torus.elevation[cell])) << "cell " << cell;
        EXPECT_GT(budget.eroded_km3, 0.0);
        EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3,
                    1e-9 * budget.eroded_km3);
    }
    // Diffusion only (K = 0, G = 0): the cone is symmetric under the point reflection
    // (x, y) -> ((W - x) % W, (H - y) % H) about the sea cell, and so is the
    // diffused surface, exactly, as long as every wrapped face is counted once (a
    // missing or doubled face between rows 0 and H - 1 breaks it at rows 1 and H - 1).
    {
        Torus torus = maketorus();
        const std::vector<float> before = torus.elevation;
        const landscape::StreamPowerBudget budget = step(torus, makeparameters(0.0, 2000.0), 0.0f);
        EXPECT_GT(budget.eroded_km3, 0.0);
        EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3, 1e-9 * budget.eroded_km3);
        EXPECT_EQ(budget.marinedeposited_km3, 0.0);
        const landscape::SphereMetric& metric = torus.metric;
        int changed = 0;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                const std::size_t cell = metric.index(x, y);
                const std::size_t mirror = metric.index((width - x) % width, (height - y) % height);
                EXPECT_EQ(torus.elevation[cell], torus.elevation[mirror]) << "x " << x << " y " << y;
                changed += torus.elevation[cell] != before[cell] ? 1 : 0;
            }
        EXPECT_GT(changed, 0);
    }
}
