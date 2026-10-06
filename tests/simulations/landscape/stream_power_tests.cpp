// Landscape.StreamPower*: the stream power step on hand-built flow graphs
// (docs/architecture/erosion-hydrology-plan.md, section 3, Tests; step A2).

#include "landscape_fixture.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
using landscape::FlowGraph;
using landscape::SphereMetric;
using landscape::StreamPowerBudget;
using landscape::StreamPowerParameters;

constexpr double earthradius_m = 6.371e6;
constexpr int gridwidth = 32;
constexpr int gridheight = 16;
constexpr int chainrow = 7;  // the row just north of the equator
constexpr int chainlength = gridwidth;
constexpr int westdirection = 7;
constexpr int southwestdirection = 6;
constexpr int southdirection = 5;
constexpr double erodibility_k = 24.0;
constexpr double areaexponent = 0.5;
// The law is scale-free (K Q^m / L does not depend on the cell size at fixed
// runoff), so the runoff alone sets the response time L / (K Q^m). 1e5 m/yr on
// cells of 1.6 million km^2 gives a few Myr at the head of the chain, so the
// chain converges well within the step cap, and a steady relief of about 4 km.
constexpr double runoff_m_yr = 1.0e5;
constexpr double seabed_m = -100.0;  // the first sea cell; each next one 100 m deeper

// A chain over row 7 (the whole row by default, else its first `length` columns):
// column x drains west to column x - 1; column 0
// is the outlet into the sea cell south-west of it (column 31 of row 8, across the
// date line). That sea cell heads a column of `seacells` sea cells down column 31,
// each draining south to the next, the last a seabed pit. Every other cell is sea
// and its own receiver. Area and discharge are accumulated by hand, discharge =
// runoff x drainage area.
struct Domain
{
    SphereMetric metric{ gridwidth, gridheight, earthradius_m };
    FlowGraph graph;
    std::vector<std::size_t> chain;      // chain[x]: the cell of column x, outlet first
    std::vector<std::size_t> seacolumn;  // the sea cells the outlet's sediment walks down
    std::vector<float> elevation, erodibility, uplift, erosion, deposition;
};

Domain makedomain(int seacells, float uplift_m_myr, int length = chainlength)
{
    Domain domain;
    const SphereMetric& metric = domain.metric;
    const std::size_t cells = metric.size();
    FlowGraph& graph = domain.graph;
    graph.receiver.resize(cells);
    for (std::size_t i = 0; i < cells; i++)
        graph.receiver[i] = static_cast<int32_t>(i);
    graph.direction.assign(cells, 0);
    graph.distance_m.assign(cells, 0.0f);
    graph.outletsea.assign(cells, -1);
    graph.land.assign(cells, 0);
    graph.area_m2.assign(cells, 0.0);
    graph.discharge_m3_yr.assign(cells, 0.0);
    domain.elevation.assign(cells, 0.0f);
    domain.erodibility.assign(cells, static_cast<float>(erodibility_k));
    domain.uplift.assign(cells, uplift_m_myr);  // sea cells too: the kernel must ignore them there
    domain.erosion.assign(cells, 0.0f);
    domain.deposition.assign(cells, 0.0f);

    for (int x = 0; x < length; x++)
    {
        const std::size_t cell = metric.index(x, chainrow);
        domain.chain.push_back(cell);
        graph.land[cell] = 1;
        graph.order.push_back(static_cast<int32_t>(cell));
        if (x > 0)
        {
            graph.receiver[cell] = static_cast<int32_t>(metric.index(x - 1, chainrow));
            graph.direction[cell] = westdirection;
            graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, westdirection));
        }
        else
        {
            const std::size_t sea = metric.neighbour(cell, southwestdirection);
            graph.direction[cell] = southwestdirection;
            graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, southwestdirection));
            graph.outletsea[cell] = static_cast<int32_t>(sea);
        }
    }
    for (int x = length - 1; x >= 0; x--)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double upstream = x + 1 < length ?graph.area_m2[domain.chain[static_cast<std::size_t>(x + 1)]] : 0.0;
        graph.area_m2[cell] = upstream + metric.area_m2(cell);
        graph.discharge_m3_yr[cell] = runoff_m_yr * graph.area_m2[cell];
    }

    std::size_t sea = static_cast<std::size_t>(graph.outletsea[domain.chain[0]]);
    for (int k = 0; k < seacells; k++)
    {
        domain.seacolumn.push_back(sea);
        domain.elevation[sea] = static_cast<float>(seabed_m * (k + 1));
        if (k + 1 < seacells)
        {
            const std::size_t next = metric.neighbour(sea, southdirection);
            graph.receiver[sea] = static_cast<int32_t>(next);
            graph.distance_m[sea] = static_cast<float>(metric.distance_m(sea, southdirection));
            sea = next;
        }
    }
    return domain;
}

// A hill on the chain: 100 m at the outlet rising by 100 m per cell upstream.
void sethill(Domain& domain)
{
    for (int x = 0; x < chainlength; x++)
        domain.elevation[domain.chain[static_cast<std::size_t>(x)]] = static_cast<float>(100.0 * (x + 1));
}

StreamPowerParameters makeparameters(double step_myr, double deposition, double diffusivity_km2_myr)
{
    StreamPowerParameters parameters;
    parameters.step_myr = step_myr;
    parameters.areaexponent = areaexponent;
    parameters.slopeexponent = 1.0;
    parameters.baselevel_m = 0.0;  // the sea surface, above the seabed of the sea cells
    parameters.depositioncoefficient = deposition;
    parameters.hillslopediffusivity_km2_myr = diffusivity_km2_myr;
    parameters.criticalslope = 0.0;
    parameters.sedimenttravel_km = 2000.0;  // about half the load stays in each sea cell of the column
    return parameters;
}

StreamPowerBudget step(Domain& domain, const StreamPowerParameters& parameters)
{
    return landscape::stepstreampower(domain.metric, domain.graph, domain.erodibility, domain.uplift, parameters,
                                      domain.elevation, domain.erosion, domain.deposition);
}

// Steps until the largest change of a chain cell is below 1e-6 m; returns the step
// count, or cap + 1 when the cap is reached first.
int runtosteadystate(Domain& domain, const StreamPowerParameters& parameters, int cap)
{
    for (int count = 1; count <= cap; count++)
    {
        const std::vector<float> before = domain.elevation;
        step(domain, parameters);
        double change = 0.0;
        for (std::size_t cell : domain.chain)
            change = std::max(change, std::abs(static_cast<double>(domain.elevation[cell]) - before[cell]));
        if (change < 1.0e-6)
            return count;
    }
    return cap + 1;
}

// The analytical steady drop of a link: z_i - z_r = U L_i / (K Q_i^m), Q in km^3/yr.
double steadydrop(const Domain& domain, std::size_t cell, double uplift_m_myr)
{
    const double discharge = domain.graph.discharge_m3_yr[cell] / 1.0e9;
    return uplift_m_myr * domain.graph.distance_m[cell] / (erodibility_k * std::pow(discharge, areaexponent));
}

double receiverelevation(const Domain& domain, int x, double baselevel)
{
    return x == 0 ? baselevel : static_cast<double>(domain.elevation[domain.chain[static_cast<std::size_t>(x - 1)]]);
}

// One float ulp at the magnitude of a value: the resolution of float storage.
double floatulp(double value)
{
    const float magnitude = static_cast<float>(std::abs(value));
    return static_cast<double>(std::nextafter(magnitude, std::numeric_limits<float>::infinity()) - magnitude);
}
}

TEST_F(Landscape, StreamPowerChainReachesAnalyticalSteadyState)
{
    const double uplift = 100.0;
    Domain domain = makedomain(1, static_cast<float>(uplift));
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    ASSERT_LE(runtosteadystate(domain, parameters, 5000), 5000);
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double drop = domain.elevation[cell] - receiverelevation(domain, x, parameters.baselevel_m);
        const double expected = steadydrop(domain, cell, uplift);
        EXPECT_NEAR(drop, expected, 1.0e-3 * expected) << "column " << x;
    }
}

TEST_F(Landscape, StreamPowerSteadyStateDoesNotDependOnStep)
{
    const double uplift = 100.0;
    Domain coarse = makedomain(1, static_cast<float>(uplift));
    Domain fine = makedomain(1, static_cast<float>(uplift));
    const StreamPowerParameters coarseparameters = makeparameters(1.0, 0.0, 0.0);
    const StreamPowerParameters fineparameters = makeparameters(0.1, 0.0, 0.0);
    ASSERT_LE(runtosteadystate(coarse, coarseparameters, 5000), 5000);
    ASSERT_LE(runtosteadystate(fine, fineparameters, 5000), 5000);
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = coarse.chain[static_cast<std::size_t>(x)];
        const double coarsedrop = coarse.elevation[cell] - receiverelevation(coarse, x, coarseparameters.baselevel_m);
        const double finedrop = fine.elevation[cell] - receiverelevation(fine, x, fineparameters.baselevel_m);
        EXPECT_NEAR(finedrop, coarsedrop, 1.0e-3 * coarsedrop) << "column " << x;
        EXPECT_NEAR(finedrop, steadydrop(fine, cell, uplift), 1.0e-3 * steadydrop(fine, cell, uplift)) << "column " << x;
        EXPECT_NEAR(fine.elevation[cell], coarse.elevation[cell], 1.0e-3 * coarse.elevation[cell]) << "column " << x;
    }
}

TEST_F(Landscape, StreamPowerHillLowersMonotonicallyAboveItsReceiver)
{
    Domain domain = makedomain(1, 0.0f);
    sethill(domain);
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    const std::vector<float> initial = domain.elevation;
    for (int count = 0; count < 20; count++)
    {
        const std::vector<float> before = domain.elevation;
        step(domain, parameters);
        for (int x = 0; x < chainlength; x++)
        {
            const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
            EXPECT_LE(domain.elevation[cell], before[cell]) << "step " << count << ", column " << x;
            EXPECT_GE(domain.elevation[cell], receiverelevation(domain, x, parameters.baselevel_m))
                << "step " << count << ", column " << x;
        }
    }
    for (std::size_t cell : domain.chain)
        EXPECT_LT(domain.elevation[cell], initial[cell]);
}

TEST_F(Landscape, StreamPowerDepositionBudgetCloses)
{
    const int seacells = 4;
    Domain domain = makedomain(seacells, 100.0f);
    sethill(domain);
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    const std::vector<float> before = domain.elevation;
    const StreamPowerBudget budget = step(domain, parameters);

    EXPECT_GT(budget.eroded_km3, 0.0);
    EXPECT_GT(budget.landdeposited_km3, 0.0);
    EXPECT_GT(budget.marinedeposited_km3, 0.0);
    EXPECT_GE(budget.depositionsweeps, 1);
    EXPECT_LT(budget.depositionsweeps, parameters.depositionmaxsweeps);  // converged, not capped
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);

    // The outlet's exported volume, from the per-cell outputs: everything that left
    // the land (hillslope transfers cancel in the sum). The outputs are float, so
    // this and the sea checks hold to float resolution, not to 1e-9; the 1e-9
    // closure is the one on the double budget above.
    double exported_m3 = 0.0;
    for (std::size_t cell : domain.chain)
    {
        EXPECT_GE(domain.erosion[cell], 0.0f);
        EXPECT_GE(domain.deposition[cell], 0.0f);
        exported_m3 += (static_cast<double>(domain.erosion[cell]) - domain.deposition[cell]) * domain.metric.area_m2(cell);
    }
    const double marine_m3 = budget.marinedeposited_km3 * 1.0e9;
    EXPECT_NEAR(exported_m3, marine_m3, 1.0e-5 * marine_m3);

    // Down the sea column each cell keeps 1 - exp(-d / travel) of what reaches it;
    // the pit keeps the rest. Sea cells change only by that deposit (no uplift).
    const double kept = 1.0 - std::exp(-domain.metric.dy_m / (parameters.sedimenttravel_km * 1000.0));
    double remaining = marine_m3;
    double seatotal_m3 = 0.0;
    for (int k = 0; k < seacells; k++)
    {
        const std::size_t cell = domain.seacolumn[static_cast<std::size_t>(k)];
        const double expected = k + 1 < seacells ? kept * remaining : remaining;
        remaining -= expected;
        const double volume = static_cast<double>(domain.deposition[cell]) * domain.metric.area_m2(cell);
        seatotal_m3 += volume;
        EXPECT_NEAR(volume, expected, 1.0e-6 * expected) << "sea cell " << k;
        EXPECT_EQ(domain.erosion[cell], 0.0f);
        EXPECT_NEAR(domain.elevation[cell], before[cell] + static_cast<double>(domain.deposition[cell]), 1.0e-2);
    }
    EXPECT_NEAR(seatotal_m3, marine_m3, 1.0e-6 * marine_m3);
}

TEST_F(Landscape, StreamPowerIsRepeatable)
{
    // The plan asks for identical results at 1 and 8 threads. The kernel is serial
    // (no parallel loop), so the thread count cannot change its result; this test
    // checks that two identical calls give bit-identical outputs.
    Domain first = makedomain(4, 100.0f);
    sethill(first);
    Domain second = first;
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    const StreamPowerBudget firstbudget = step(first, parameters);
    const StreamPowerBudget secondbudget = step(second, parameters);
    const std::size_t bytes = first.elevation.size() * sizeof(float);
    EXPECT_EQ(std::memcmp(first.elevation.data(), second.elevation.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(first.erosion.data(), second.erosion.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(first.deposition.data(), second.deposition.data(), bytes), 0);
    EXPECT_EQ(firstbudget.eroded_km3, secondbudget.eroded_km3);
    EXPECT_EQ(firstbudget.landdeposited_km3, secondbudget.landdeposited_km3);
    EXPECT_EQ(firstbudget.marinedeposited_km3, secondbudget.marinedeposited_km3);
    EXPECT_EQ(firstbudget.depositionsweeps, secondbudget.depositionsweeps);
}

TEST_F(Landscape, StreamPowerWithoutDepositionIsBraunWillett)
{
    // G = 0: z_new = (z_old + U dt + F z_r) / (1 + F) with F = K Q^m dt / L, z_r the
    // receiver's new elevation and the base level for the outlet. The expectation is
    // computed in double from the outlet up; the kernel stores float, so the match
    // is to one float ulp (2^-24 relative is coarser than 1e-9).
    const double uplift = 100.0;
    const double dt = 1.0;
    Domain domain = makedomain(1, static_cast<float>(uplift));
    sethill(domain);
    const StreamPowerParameters parameters = makeparameters(dt, 0.0, 0.0);
    const std::vector<float> before = domain.elevation;
    const StreamPowerBudget budget = step(domain, parameters);
    EXPECT_EQ(budget.depositionsweeps, 0);
    EXPECT_EQ(budget.landdeposited_km3, 0.0);

    double receiver = parameters.baselevel_m;
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double discharge = domain.graph.discharge_m3_yr[cell] / 1.0e9;
        const double rate = erodibility_k * std::pow(discharge, areaexponent) / domain.graph.distance_m[cell];
        const double factor = rate * dt;
        const double expected = (static_cast<double>(before[cell]) + uplift * dt + factor * receiver) / (1.0 + factor);
        EXPECT_NEAR(domain.elevation[cell], expected, floatulp(expected)) << "column " << x;
        receiver = expected;
    }
}

TEST_F(Landscape, StreamPowerDepositionClosedFormOnOneCell)
{
    // One land cell (an outlet), G = 1: S_up = 0, so the first sweep is exact. With
    // A = a, F' = k dt A / (A + G a) = k dt / (1 + G), b = z_old + U dt,
    // z_new = (b + F' z_r) / (1 + F'), E = k (z_new - z_r), D = G / (1 + G) E.
    // Checked to one float ulp, as the Braun-Willett test.
    const double uplift = 100.0;
    const double dt = 1.0;
    const double deposition = 1.0;
    Domain domain = makedomain(1, static_cast<float>(uplift), 1);
    const std::size_t cell = domain.chain[0];
    domain.elevation[cell] = 500.0f;
    const StreamPowerParameters parameters = makeparameters(dt, deposition, 0.0);
    const std::vector<float> before = domain.elevation;
    const StreamPowerBudget budget = step(domain, parameters);
    EXPECT_GE(budget.depositionsweeps, 1);

    const double discharge = domain.graph.discharge_m3_yr[cell] / 1.0e9;
    const double rate = erodibility_k * std::pow(discharge, areaexponent) / domain.graph.distance_m[cell];
    const double factor = rate * dt / (1.0 + deposition);
    const double base = static_cast<double>(before[cell]) + uplift * dt;
    const double receiver = parameters.baselevel_m;
    const double expected = (base + factor * receiver) / (1.0 + factor);
    const double incision = rate * (expected - receiver);
    const double expecteddeposition = deposition / (1.0 + deposition) * incision * dt;
    EXPECT_NEAR(domain.elevation[cell], expected, floatulp(expected));
    EXPECT_NEAR(domain.erosion[cell], incision * dt, floatulp(incision * dt));
    EXPECT_NEAR(domain.deposition[cell], expecteddeposition, floatulp(expecteddeposition));
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);
}

TEST_F(Landscape, StreamPowerLandRootKeepsItsSediment)
{
    // A land cell that is its own receiver without a sea cell (a region the flood
    // never reached) is accepted, is not incised (E = 0) and keeps what reaches it;
    // alone, it only rises by U dt and the budget is zero.
    const double uplift = 100.0;
    const double dt = 1.0;
    Domain domain = makedomain(1, static_cast<float>(uplift), 1);
    const std::size_t cell = domain.chain[0];
    domain.graph.outletsea[cell] = -1;
    domain.graph.direction[cell] = 0;
    domain.elevation[cell] = 100.0f;
    const StreamPowerParameters parameters = makeparameters(dt, 1.0, 0.0);
    StreamPowerBudget budget;
    ASSERT_NO_THROW(budget = step(domain, parameters));
    EXPECT_EQ(domain.elevation[cell], static_cast<float>(100.0 + uplift * dt));
    EXPECT_EQ(domain.erosion[cell], 0.0f);
    EXPECT_EQ(domain.deposition[cell], 0.0f);
    EXPECT_EQ(budget.eroded_km3, 0.0);
    EXPECT_EQ(budget.landdeposited_km3, 0.0);
    EXPECT_EQ(budget.marinedeposited_km3, 0.0);
    EXPECT_EQ(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3);
}

TEST_F(Landscape, StreamPowerSedimentFluxAtOutletIsTheExport)
{
    // Erosion plan C5a: the flux an outlet passes on is the load marine deposition
    // spreads, so flux x dt equals the step's marine volume; sea cells carry none.
    Domain domain = makedomain(4, 100.0f);
    sethill(domain);
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    std::vector<double> flux;
    const StreamPowerBudget budget = landscape::stepstreampower(domain.metric, domain.graph, domain.erodibility,
        domain.uplift, parameters, domain.elevation, domain.erosion, domain.deposition, &flux);
    ASSERT_EQ(flux.size(), domain.metric.size());
    const double marine_m3 = budget.marinedeposited_km3 * 1.0e9;
    ASSERT_GT(marine_m3, 0.0);
    EXPECT_NEAR(flux[domain.chain[0]] * parameters.step_myr, marine_m3, 1.0e-9 * marine_m3);
    for (std::size_t cell : domain.seacolumn)
        EXPECT_EQ(flux[cell], 0.0);
}

TEST_F(Landscape, StreamPowerSedimentFluxGrowsDownstreamWithoutDeposition)
{
    // G = 0 on a hill: every chain cell is incised and keeps nothing, so the flux each
    // passes on is its own incision plus everything upstream, growing towards the outlet.
    Domain domain = makedomain(1, 0.0f);
    sethill(domain);
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    std::vector<double> flux;
    landscape::stepstreampower(domain.metric, domain.graph, domain.erodibility, domain.uplift, parameters,
                               domain.elevation, domain.erosion, domain.deposition, &flux);
    EXPECT_GT(flux[domain.chain[static_cast<std::size_t>(chainlength - 1)]], 0.0);
    for (int x = 0; x + 1 < chainlength; x++)
        EXPECT_GT(flux[domain.chain[static_cast<std::size_t>(x)]], flux[domain.chain[static_cast<std::size_t>(x + 1)]])
            << "column " << x;
}

namespace
{
// The sea cells off the chain and off the sea column: their own receivers, at `depth_m`.
void setopensea(Domain& domain, double depth_m)
{
    for (std::size_t cell = 0; cell < domain.metric.size(); cell++)
    {
        if (domain.graph.land[cell] != 0
            || std::find(domain.seacolumn.begin(), domain.seacolumn.end(), cell) != domain.seacolumn.end())
            continue;
        domain.elevation[cell] = static_cast<float>(-depth_m);
    }
}

// Every sea cell at or below the delta top, and the sea's deposits sum to the
// step's marine volume.
void expectfilledandspilled(const Domain& domain, const StreamPowerBudget& budget, double top)
{
    double seatotal_m3 = 0.0;
    for (std::size_t cell = 0; cell < domain.metric.size(); cell++)
    {
        if (domain.graph.land[cell] != 0)
            continue;
        EXPECT_LE(domain.elevation[cell], top + 1.0e-4) << "sea cell " << cell;
        seatotal_m3 += static_cast<double>(domain.deposition[cell]) * domain.metric.area_m2(cell);
    }
    const double marine_m3 = budget.marinedeposited_km3 * 1.0e9;
    EXPECT_NEAR(seatotal_m3, marine_m3, 1.0e-6 * marine_m3);
    EXPECT_EQ(budget.strandedloads, 0);
    EXPECT_EQ(budget.stranded_km3, 0.0);
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);
}
}

TEST_F(Landscape, StreamPowerProgradationFillsShallowMouthCellsToTheDeltaTop)
{
    // Erosion plan C5b: a shallow sea column (1 to 4 m deep) in front of an eroding
    // hill, the open sea 1 km deep. Each sea cell on the path keeps no more than raises
    // it to the delta top, so every cell of the column, the pit included, stands at
    // baselevel + delta top (land under the whole-metre mask); the pit's remainder
    // spills into the open sea, and the budget still closes.
    const int seacells = 4;
    Domain domain = makedomain(seacells, 100.0f);
    sethill(domain);
    setopensea(domain, 1000.0);
    for (int k = 0; k < seacells; k++)
        domain.elevation[domain.seacolumn[static_cast<std::size_t>(k)]] = static_cast<float>(-1.0 * (k + 1));
    StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    parameters.deltatop_m = 1.5;
    const std::vector<float> before = domain.elevation;
    const StreamPowerBudget budget = step(domain, parameters);

    EXPECT_GT(budget.marinedeposited_km3, 0.0);
    const double top = parameters.baselevel_m + parameters.deltatop_m;
    for (int k = 0; k < seacells; k++)
    {
        const std::size_t cell = domain.seacolumn[static_cast<std::size_t>(k)];
        EXPECT_NEAR(domain.elevation[cell], top, 1.0e-4) << "sea cell " << k;
        EXPECT_NEAR(domain.deposition[cell], top - before[cell], 1.0e-4) << "sea cell " << k;
        EXPECT_GT(static_cast<double>(std::lround(domain.elevation[cell])), parameters.baselevel_m) << "sea cell " << k;
    }
    expectfilledandspilled(domain, budget, top);
}

TEST_F(Landscape, StreamPowerProgradationCountsTheEarlierOutletsDeposit)
{
    // Two outlets 2 m high on row 7 (columns 1 and 3) drain into one sea cell between
    // them (column 2, a pit 1 m deep); the open sea is 1 km deep. Their erodibility
    // is so high that each exports nearly its whole 2 m: the first fits into the pit's
    // 2.5 m of room, the second does not, so the pit stands at the delta top only if
    // the first outlet's deposit is counted, and the excess spills into the open sea.
    Domain domain = makedomain(1, 0.0f, 1);
    const SphereMetric& metric = domain.metric;
    FlowGraph& graph = domain.graph;
    for (std::size_t cell : domain.chain)
    {
        graph.land[cell] = 0;
        graph.outletsea[cell] = -1;
        graph.receiver[cell] = static_cast<int32_t>(cell);
        graph.direction[cell] = 0;
        graph.distance_m[cell] = 0.0f;
        graph.area_m2[cell] = 0.0;
        graph.discharge_m3_yr[cell] = 0.0;
    }
    graph.order.clear();
    domain.chain.clear();
    domain.seacolumn.clear();
    setopensea(domain, 1000.0);
    const std::size_t pit = metric.index(2, chainrow);
    domain.elevation[pit] = -1.0f;
    for (int x : { 1, 3 })
    {
        const std::size_t cell = metric.index(x, chainrow);
        const int direction = x < 2 ? 3 : westdirection;
        graph.land[cell] = 1;
        graph.order.push_back(static_cast<int32_t>(cell));
        graph.direction[cell] = static_cast<uint8_t>(direction);
        graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, direction));
        graph.outletsea[cell] = static_cast<int32_t>(pit);
        graph.area_m2[cell] = metric.area_m2(cell);
        graph.discharge_m3_yr[cell] = runoff_m_yr * graph.area_m2[cell];
        domain.elevation[cell] = 2.0f;
        domain.erodibility[cell] = 1.0e7f;
        domain.chain.push_back(cell);
    }
    StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    parameters.deltatop_m = 1.5;
    const StreamPowerBudget budget = step(domain, parameters);

    const double area = metric.area_m2(pit);
    const double room = 2.5 * area;
    const double first = static_cast<double>(domain.erosion[domain.chain[0]]) * area;
    const double second = static_cast<double>(domain.erosion[domain.chain[1]]) * area;
    ASSERT_LT(first, room);
    ASSERT_LT(second, room);
    ASSERT_GT(first + second, room);
    const double top = parameters.baselevel_m + parameters.deltatop_m;
    EXPECT_NEAR(domain.elevation[pit], top, 1.0e-4);
    const double marine_m3 = budget.marinedeposited_km3 * 1.0e9;
    EXPECT_NEAR(marine_m3, first + second, 1.0e-5 * marine_m3);
    EXPECT_NEAR(marine_m3 - static_cast<double>(domain.deposition[pit]) * area, first + second - room, 1.0e-3 * room);
    expectfilledandspilled(domain, budget, top);
}

TEST_F(Landscape, StreamPowerProgradationCappedWalkSpillsFromAShallowPit)
{
    // The sea column is 1 km deep down to a pit only 1 m deep; the open sea is 1 km
    // deep. The walk leaves its e-folding shares on the deep cells and ends at the
    // pit, which fills to the delta top; the rest spills into the open sea.
    const int seacells = 4;
    Domain domain = makedomain(seacells, 100.0f);
    sethill(domain);
    setopensea(domain, 1000.0);
    for (int k = 0; k + 1 < seacells; k++)
        domain.elevation[domain.seacolumn[static_cast<std::size_t>(k)]] = -1000.0f;
    const std::size_t pit = domain.seacolumn.back();
    domain.elevation[pit] = -1.0f;
    StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    parameters.deltatop_m = 1.5;
    const StreamPowerBudget budget = step(domain, parameters);

    const double top = parameters.baselevel_m + parameters.deltatop_m;
    EXPECT_NEAR(domain.elevation[pit], top, 1.0e-4);
    double spilled_m3 = 0.0;
    for (std::size_t cell = 0; cell < domain.metric.size(); cell++)
        if (domain.graph.land[cell] == 0
            && std::find(domain.seacolumn.begin(), domain.seacolumn.end(), cell) == domain.seacolumn.end())
            spilled_m3 += static_cast<double>(domain.deposition[cell]) * domain.metric.area_m2(cell);
    EXPECT_GT(spilled_m3, 0.0);
    expectfilledandspilled(domain, budget, top);
}

TEST_F(Landscape, StreamPowerProgradationBuildsNothingInADeepSea)
{
    // A deep sea (1 to 4 km) and a small load (one low outlet cell): no sea cell gets
    // near the delta top, so the deposits are those without the cap, bit for bit,
    // and every sea cell stays below the sea level.
    const int seacells = 4;
    Domain capped = makedomain(seacells, 0.0f, 1);
    capped.elevation[capped.chain[0]] = 10.0f;
    for (int k = 0; k < seacells; k++)
        capped.elevation[capped.seacolumn[static_cast<std::size_t>(k)]] = static_cast<float>(-1000.0 * (k + 1));
    Domain uncapped = capped;
    StreamPowerParameters parameters = makeparameters(1.0, 1.0, 0.0);
    const StreamPowerBudget uncappedbudget = step(uncapped, parameters);
    parameters.deltatop_m = 1.5;
    const StreamPowerBudget budget = step(capped, parameters);

    EXPECT_GT(budget.marinedeposited_km3, 0.0);
    EXPECT_EQ(budget.marinedeposited_km3, uncappedbudget.marinedeposited_km3);
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);
    for (std::size_t cell : capped.seacolumn)
    {
        EXPECT_EQ(capped.deposition[cell], uncapped.deposition[cell]);
        EXPECT_EQ(capped.elevation[cell], uncapped.elevation[cell]);
        EXPECT_LT(capped.elevation[cell], parameters.baselevel_m);
    }
}

TEST_F(Landscape, StreamPowerProgradationStrandsALoadInAnEnclosedSea)
{
    // A one-cell sea (column 2 of row 7, 1 m deep) enclosed by land, fed by an outlet
    // 10 m high at column 1 whose erodibility makes it export nearly all of its 10 m:
    // more than the sea cell's 2.5 m of room. The cell fills to the delta top, the
    // spill reaches no sea cell, and the remainder stays on the cell as one stranded load.
    Domain domain = makedomain(1, 0.0f, 1);
    const SphereMetric& metric = domain.metric;
    FlowGraph& graph = domain.graph;
    for (std::size_t cell : domain.chain)
    {
        graph.land[cell] = 0;
        graph.outletsea[cell] = -1;
        graph.direction[cell] = 0;
        graph.distance_m[cell] = 0.0f;
        graph.area_m2[cell] = 0.0;
        graph.discharge_m3_yr[cell] = 0.0;
    }
    graph.order.clear();
    domain.chain.clear();
    domain.seacolumn.clear();
    setopensea(domain, 1000.0);
    const std::size_t sea = metric.index(2, chainrow);
    domain.elevation[sea] = -1.0f;
    for (int direction = 1; direction <= landscape::directioncount; direction++)
    {
        const std::size_t cell = metric.neighbour(sea, direction);
        graph.land[cell] = 1;  // land around the sea cell, outside the order: unchanged by the step
        domain.elevation[cell] = 10.0f;
    }
    const std::size_t outlet = metric.index(1, chainrow);
    const int eastdirection = 3;
    graph.order.push_back(static_cast<int32_t>(outlet));
    graph.direction[outlet] = static_cast<uint8_t>(eastdirection);
    graph.distance_m[outlet] = static_cast<float>(metric.distance_m(outlet, eastdirection));
    graph.outletsea[outlet] = static_cast<int32_t>(sea);
    graph.area_m2[outlet] = metric.area_m2(outlet);
    graph.discharge_m3_yr[outlet] = runoff_m_yr * graph.area_m2[outlet];
    domain.erodibility[outlet] = 1.0e7f;
    StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    parameters.deltatop_m = 1.5;
    const StreamPowerBudget budget = step(domain, parameters);

    const double area = metric.area_m2(sea);
    const double room = 2.5 * area;
    const double load = static_cast<double>(domain.erosion[outlet]) * metric.area_m2(outlet);
    ASSERT_GT(load, room);
    EXPECT_EQ(budget.strandedloads, 1);
    EXPECT_NEAR(budget.stranded_km3, (load - room) / 1.0e9, 1.0e-6 * load / 1.0e9);
    EXPECT_NEAR(budget.marinedeposited_km3, load / 1.0e9, 1.0e-6 * load / 1.0e9);
    EXPECT_NEAR(domain.elevation[sea], -1.0 + load / area, 1.0e-4 * load / area);
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);
}

TEST_F(Landscape, StreamPowerRejectsInvalidInput)
{
    Domain domain = makedomain(1, 100.0f);
    StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    parameters.slopeexponent = 2.0;
    EXPECT_THROW(step(domain, parameters), std::invalid_argument);
    parameters.slopeexponent = 1.0;
    Domain landsea = domain;
    landsea.graph.outletsea[landsea.chain[0]] = static_cast<int32_t>(landsea.chain[1]);
    EXPECT_THROW(step(landsea, parameters), std::invalid_argument);
    domain.uplift.pop_back();
    EXPECT_THROW(step(domain, parameters), std::invalid_argument);
}
