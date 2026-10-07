// Suite Landscape, kernel Glacial (erosion plan Phase G, step G1): the ELA, the ice
// flux, the thickness, and the glacial term of stepstreampower (its tests are
// TEST_F(Landscape, StreamPowerGlacial<What>), so that both kernels' filters find them).

#include "landscape_fixture.hpp"

#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/glacial.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace
{
using landscape::FlowGraph;
using landscape::GlacialInput;
using landscape::GlacialParameters;
using landscape::SphereMetric;
using landscape::StreamPowerBudget;
using landscape::StreamPowerParameters;

// The chain domain of stream_power_tests.cpp, copied: row 7 of a 32 x 16 grid
// drains west to the outlet at column 0, whose sea cell (column 31 of row 8) heads
// a column of sea cells down column 31.
constexpr double earthradius_m = 6.371e6;
constexpr int gridwidth = 32;
constexpr int gridheight = 16;
constexpr int chainrow = 7;
constexpr int chainlength = gridwidth;
constexpr int westdirection = 7;
constexpr int southwestdirection = 6;
constexpr int southdirection = 5;
constexpr double erodibility_k = 24.0;
constexpr double areaexponent = 0.5;
constexpr double runoff_m_yr = 1.0e5;
constexpr double seabed_m = -100.0;
constexpr double flotation = 0.892;

struct Domain
{
    SphereMetric metric{ gridwidth, gridheight, earthradius_m };
    FlowGraph graph;
    std::vector<std::size_t> chain;
    std::vector<std::size_t> seacolumn;
    std::vector<float> elevation, erodibility, uplift, erosion, deposition;
};

Domain makedomain(int seacells, float uplift_m_myr)
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
    domain.uplift.assign(cells, uplift_m_myr);
    domain.erosion.assign(cells, 0.0f);
    domain.deposition.assign(cells, 0.0f);

    for (int x = 0; x < chainlength; x++)
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
    for (int x = chainlength - 1; x >= 0; x--)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double upstream =
            x + 1 < chainlength ? graph.area_m2[domain.chain[static_cast<std::size_t>(x + 1)]] : 0.0;
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

// 100 m at the outlet rising by 100 m per cell upstream.
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
    parameters.baselevel_m = 0.0;
    parameters.depositioncoefficient = deposition;
    parameters.hillslopediffusivity_km2_myr = diffusivity_km2_myr;
    parameters.criticalslope = 0.0;
    parameters.sedimenttravel_km = 2000.0;
    return parameters;
}

// The ice of a step: flux and thickness per cell, zero off the chain.
struct Ice
{
    std::vector<double> flux;
    std::vector<float> thickness;
    std::vector<float> glacialerosion;
    GlacialInput input;
};

// Ice flux = share(x) x discharge and thickness(x) on chain column x.
template <typename Share, typename Thickness>
void setice(const Domain& domain, Ice& ice, Share share, Thickness thickness)
{
    const std::size_t cells = domain.metric.size();
    ice.flux.assign(cells, 0.0);
    ice.thickness.assign(cells, 0.0f);
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        ice.flux[cell] = share(x) * domain.graph.discharge_m3_yr[cell];
        ice.thickness[cell] = static_cast<float>(thickness(x));
    }
    ice.input = GlacialInput{};
    ice.input.iceflux_m3_yr = &ice.flux;
    ice.input.thickness_m = &ice.thickness;
    ice.input.erodibilityratio = 2.0;
    ice.input.flotationratio = flotation;
    ice.input.glacialerosion_m = &ice.glacialerosion;
}

StreamPowerBudget step(Domain& domain, const StreamPowerParameters& parameters, const GlacialInput* glacial,
                       std::vector<double>* flux = nullptr)
{
    return landscape::stepstreampower(domain.metric, domain.graph, domain.erodibility, domain.uplift, parameters,
                                      domain.elevation, domain.erosion, domain.deposition, flux, glacial);
}

double bed(const Domain& domain, int x)
{
    return static_cast<double>(domain.elevation[domain.chain[static_cast<std::size_t>(x)]]);
}
}

TEST_F(Landscape, GlacialEquilibriumLineTemperature)
{
    EXPECT_NEAR(landscape::equilibriumlinetemperature(645.0), 0.0, 1.0e-12);
    EXPECT_NEAR(landscape::equilibriumlinetemperature(950.0), 1.0, 1.0e-12);
    EXPECT_NEAR(landscape::equilibriumlinetemperature(358.0), -1.0, 1.0e-12);
    // P = 645 + 296 T + 9 T^2 holds at the returned T.
    for (double precipitation : { 0.0, 150.0, 1500.0, 6000.0 })
    {
        const double t = landscape::equilibriumlinetemperature(precipitation);
        EXPECT_NEAR(645.0 + 296.0 * t + 9.0 * t * t, precipitation, 1.0e-9 * std::max(1.0, precipitation));
    }
    // Negative precipitation counts as none.
    EXPECT_EQ(landscape::equilibriumlinetemperature(-100.0), landscape::equilibriumlinetemperature(0.0));
    double previous = landscape::equilibriumlinetemperature(0.0);
    for (double precipitation = 50.0; precipitation <= 8000.0; precipitation += 50.0)
    {
        const double t = landscape::equilibriumlinetemperature(precipitation);
        EXPECT_GT(t, previous) << "P = " << precipitation;
        previous = t;
    }
}

TEST_F(Landscape, GlacialEquilibriumLine)
{
    const SphereMetric metric(gridwidth, gridheight, earthradius_m);
    const std::size_t cells = metric.size();
    const double lapse = 0.0065;
    const double cooling = 4.0;
    std::vector<float> elevation(cells, 0.0f), warmest(cells, 0.0f), precipitation(cells, 645.0f);
    // Cell 0: T_ela(950) = 1, so (5 - 4 - 1) / lapse = 0: the ELA is the elevation.
    elevation[0] = 1000.0f;
    warmest[0] = 5.0f;
    precipitation[0] = 950.0f;
    // Cell 1: T_ela(645) = 0, ELA = 0 + 6 / 0.0065 = 923.0769...
    elevation[1] = 0.0f;
    warmest[1] = 10.0f;
    precipitation[1] = 645.0f;
    // Cell 2: T_ela(358) = -1, ELA = 2000 + (-3 - 4 + 1) / 0.0065 = 1076.923...
    elevation[2] = 2000.0f;
    warmest[2] = -3.0f;
    precipitation[2] = 358.0f;
    std::vector<float> ela;
    landscape::equilibriumline(metric, elevation, warmest, precipitation, lapse, cooling, ela);
    ASSERT_EQ(ela.size(), cells);
    EXPECT_NEAR(ela[0], 1000.0, 1.0e-4);
    EXPECT_NEAR(ela[1], 923.0769231, 1.0e-4);
    EXPECT_NEAR(ela[2], 1076.9230769, 1.0e-4);

    // warmest = T0 - lapse z: the ELA is (T0 - cooling - T_ela(P)) / lapse on every cell.
    const double t0 = 12.0;
    for (std::size_t i = 0; i < cells; i++)
    {
        elevation[i] = static_cast<float>(10.0 * static_cast<double>(i % 400));
        warmest[i] = static_cast<float>(t0 - lapse * elevation[i]);
        precipitation[i] = 1200.0f;
    }
    landscape::equilibriumline(metric, elevation, warmest, precipitation, lapse, cooling, ela);
    const double expected = (t0 - cooling - landscape::equilibriumlinetemperature(1200.0)) / lapse;
    for (std::size_t i = 0; i < cells; i++)
        EXPECT_NEAR(ela[i], expected, 0.05) << "cell " << i;  // warmest is float: 1e-6 C over 0.0065 C/m

    EXPECT_THROW(landscape::equilibriumline(metric, elevation, warmest, precipitation, 0.0, cooling, ela),
                 std::invalid_argument);
    precipitation.pop_back();
    EXPECT_THROW(landscape::equilibriumline(metric, elevation, warmest, precipitation, lapse, cooling, ela),
                 std::invalid_argument);
}

TEST_F(Landscape, GlacialIceFluxAblatesToATerminus)
{
    // The hill (100 m at the outlet to 3200 m at the head) under an ELA of 1500 m with
    // a full-ice range of 500 m: in units of runoff x cell area (one row, so every
    // chain cell has the same area), the source of column x is min(0.2 x - 2.8, 1):
    // 1 from column 19 up, 0 at column 14, negative below. The flux is 13 at column
    // 19, peaks at 15 at columns 15 and 14, falls to 1.8 at column 3 (the terminus)
    // and is 0 from column 2 to the outlet.
    Domain domain = makedomain(4, 0.0f);
    sethill(domain);
    const std::size_t cells = domain.metric.size();
    const std::vector<float> ela(cells, 1500.0f);
    const std::vector<float> runoff(cells, static_cast<float>(runoff_m_yr));
    // The discharge as the routing accumulates it, cell by cell in the same order, so
    // that flux <= discharge holds without rounding slack.
    landscape::accumulateflow(domain.metric, runoff, domain.graph);
    std::vector<double> flux;
    landscape::accumulateice(domain.metric, domain.graph, domain.elevation, ela, runoff, 500.0, flux);
    ASSERT_EQ(flux.size(), cells);

    const double unit = runoff_m_yr * domain.metric.area_m2(domain.chain[0]);
    auto at = [&](int x) { return flux[domain.chain[static_cast<std::size_t>(x)]]; };
    EXPECT_NEAR(at(19), 13.0 * unit, 1.0e-9 * 13.0 * unit);
    EXPECT_NEAR(at(15), 15.0 * unit, 1.0e-9 * 15.0 * unit);
    EXPECT_NEAR(at(14), 15.0 * unit, 1.0e-9 * 15.0 * unit);
    EXPECT_NEAR(at(9), 12.0 * unit, 1.0e-9 * 12.0 * unit);
    EXPECT_NEAR(at(3), 1.8 * unit, 1.0e-9 * 1.8 * unit);
    for (int x = 30; x >= 15; x--)
        EXPECT_GT(at(x), at(x + 1)) << "column " << x;  // grows above the ELA
    for (int x = 13; x >= 3; x--)
        EXPECT_LT(at(x), at(x + 1)) << "column " << x;  // ablates below it
    for (int x = 2; x >= 0; x--)
        EXPECT_EQ(at(x), 0.0) << "column " << x;  // beyond the terminus
    for (std::size_t cell = 0; cell < cells; cell++)
    {
        if (domain.graph.land[cell] == 0)
            EXPECT_EQ(flux[cell], 0.0) << "sea cell " << cell;
        else
            EXPECT_LE(flux[cell], domain.graph.discharge_m3_yr[cell]) << "cell " << cell;
    }

    std::vector<float> shortrunoff = runoff;
    shortrunoff.pop_back();
    EXPECT_THROW(landscape::accumulateice(domain.metric, domain.graph, domain.elevation, ela, shortrunoff, 500.0, flux),
                 std::invalid_argument);
}

TEST_F(Landscape, GlacialThickness)
{
    // Defaults: h = 0.25 x 400 m x (A / 1 km^2)^0.3, capped at 3000 m. A flux of
    // 1e6 m^3/yr is 1 km^2 (100 m); 1e8 is 100 km^2 (100 x 10^0.6 = 398.107 m);
    // 1e12 is 1e6 km^2 (6310 m, capped).
    const GlacialParameters parameters;
    const std::vector<double> flux = { 1.0e6, 1.0e8, 1.0e12, 0.0, -5.0 };
    std::vector<float> thickness;
    landscape::icethickness(flux, parameters, thickness);
    ASSERT_EQ(thickness.size(), flux.size());
    EXPECT_NEAR(thickness[0], 100.0, 1.0e-4);
    EXPECT_NEAR(thickness[1], 398.1071706, 1.0e-4);
    EXPECT_EQ(thickness[2], 3000.0f);
    EXPECT_EQ(thickness[3], 0.0f);
    EXPECT_EQ(thickness[4], 0.0f);
}

TEST_F(Landscape, StreamPowerGlacialWithoutIceIsUnchanged)
{
    Domain plain = makedomain(4, 100.0f);
    sethill(plain);
    Domain iced = plain;
    Ice ice;
    setice(iced, ice, [](int) { return 0.0; }, [](int) { return 0.0; });
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    std::vector<double> plainflux, icedflux;
    const StreamPowerBudget plainbudget = step(plain, parameters, nullptr, &plainflux);
    const StreamPowerBudget icedbudget = step(iced, parameters, &ice.input, &icedflux);

    const std::size_t bytes = plain.elevation.size() * sizeof(float);
    EXPECT_EQ(std::memcmp(plain.elevation.data(), iced.elevation.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(plain.erosion.data(), iced.erosion.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(plain.deposition.data(), iced.deposition.data(), bytes), 0);
    ASSERT_EQ(plainflux.size(), icedflux.size());
    EXPECT_EQ(std::memcmp(plainflux.data(), icedflux.data(), plainflux.size() * sizeof(double)), 0);
    EXPECT_EQ(plainbudget.eroded_km3, icedbudget.eroded_km3);
    EXPECT_EQ(plainbudget.landdeposited_km3, icedbudget.landdeposited_km3);
    EXPECT_EQ(plainbudget.marinedeposited_km3, icedbudget.marinedeposited_km3);
    EXPECT_EQ(plainbudget.depositionsweeps, icedbudget.depositionsweeps);
    EXPECT_GT(plainbudget.eroded_km3, 0.0);
    ASSERT_EQ(ice.glacialerosion.size(), plain.metric.size());
    for (float value : ice.glacialerosion)
        EXPECT_EQ(value, 0.0f);
}

TEST_F(Landscape, StreamPowerGlacialSteadyState)
{
    // gamma = 1 everywhere, no thickness: K_eff = erodibilityratio x K, so each link's
    // steady drop is U L / (K_eff Q^m).
    const double uplift = 100.0;
    Domain domain = makedomain(1, static_cast<float>(uplift));
    Ice ice;
    setice(domain, ice, [](int) { return 1.0; }, [](int) { return 0.0; });
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    int count = 1;
    for (; count <= 5000; count++)
    {
        const std::vector<float> before = domain.elevation;
        step(domain, parameters, &ice.input);
        double change = 0.0;
        for (std::size_t cell : domain.chain)
            change = std::max(change, std::abs(static_cast<double>(domain.elevation[cell]) - before[cell]));
        if (change < 1.0e-6)
            break;
    }
    ASSERT_LE(count, 5000);
    const double keff = ice.input.erodibilityratio * erodibility_k;
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double receiver = x == 0 ? parameters.baselevel_m : bed(domain, x - 1);
        const double discharge = domain.graph.discharge_m3_yr[cell] / 1.0e9;
        const double expected = uplift * domain.graph.distance_m[cell] / (keff * std::pow(discharge, areaexponent));
        EXPECT_NEAR(bed(domain, x) - receiver, expected, 1.0e-3 * expected) << "column " << x;
        // All of the incision is glacial.
        EXPECT_NEAR(ice.glacialerosion[cell], domain.erosion[cell], 1.0e-6 * domain.erosion[cell]) << "column " << x;
    }
}

TEST_F(Landscape, StreamPowerGlacialOverdeepensAndFloats)
{
    // A coastal glacier (gamma = 1) whose thickness grows upstream, h = 100 + 20 x m,
    // on the hill, without uplift: the ice surface relaxes towards the sea level until
    // each bed stands on its flotation floor -0.892 h, so the bed of the outlet ends
    // below the sea (a fjord, -89.2 m) and every bed below its receiver's (overdeepened).
    Domain domain = makedomain(4, 0.0f);
    sethill(domain);
    Ice ice;
    setice(domain, ice, [](int) { return 1.0; }, [](int x) { return 100.0 + 20.0 * x; });
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    StreamPowerBudget budget;
    for (int count = 0; count < 500; count++)
        budget = step(domain, parameters, &ice.input);
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3,
                1.0e-9 * std::max(budget.eroded_km3, 1.0e-30));

    const double outletfloor = parameters.baselevel_m - flotation * ice.thickness[domain.chain[0]];
    EXPECT_LT(bed(domain, 0), parameters.baselevel_m);
    EXPECT_GE(bed(domain, 0), outletfloor - 1.0e-3);
    EXPECT_NEAR(bed(domain, 0), outletfloor, 1.0e-3);
    int overdeepened = 0;
    for (int x = 1; x < chainlength; x++)
    {
        const double floor = parameters.baselevel_m - flotation * ice.thickness[domain.chain[static_cast<std::size_t>(x)]];
        EXPECT_GE(bed(domain, x), floor - 1.0e-3) << "column " << x;
        if (bed(domain, x) < bed(domain, x - 1))
            overdeepened++;
    }
    EXPECT_EQ(overdeepened, chainlength - 1);
}

TEST_F(Landscape, StreamPowerGlacialBudgetCloses)
{
    // Ice over the whole chain (gamma 1 from column 16 up, 0.5 below; h = 200 + 10 x m),
    // deposition and hillslope diffusion on, a slow uplift: over 30 steps the outlet
    // reaches its flotation floor, every step's budget closes and the glacial part of
    // the erosion stays within [0, erosion].
    Domain domain = makedomain(4, 20.0f);
    sethill(domain);
    Ice ice;
    setice(domain, ice, [](int x) { return x >= 16 ? 1.0 : 0.5; }, [](int x) { return 200.0 + 10.0 * x; });
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    for (int count = 0; count < 30; count++)
    {
        const StreamPowerBudget budget = step(domain, parameters, &ice.input);
        ASSERT_GT(budget.eroded_km3, 0.0) << "step " << count;
        EXPECT_GT(budget.landdeposited_km3, 0.0) << "step " << count;
        EXPECT_GT(budget.marinedeposited_km3, 0.0) << "step " << count;
        EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3,
                    1.0e-9 * budget.eroded_km3)
            << "step " << count;
        int glaciated = 0;
        for (std::size_t cell = 0; cell < domain.metric.size(); cell++)
        {
            EXPECT_GE(ice.glacialerosion[cell], 0.0f) << "step " << count << ", cell " << cell;
            EXPECT_LE(ice.glacialerosion[cell], domain.erosion[cell]) << "step " << count << ", cell " << cell;
            if (ice.glacialerosion[cell] > 0.0f)
                glaciated++;
        }
        EXPECT_GT(glaciated, 0) << "step " << count;
    }
    EXPECT_LT(bed(domain, 0), parameters.baselevel_m);
}

TEST_F(Landscape, StreamPowerGlacialIsRepeatable)
{
    // The kernel is serial (no parallel loop), so the thread count cannot change its
    // result; as StreamPowerIsRepeatable, two identical calls give bit-identical outputs.
    Domain first = makedomain(4, 100.0f);
    sethill(first);
    Domain second = first;
    Ice firstice, secondice;
    auto share = [](int x) { return x >= 16 ? 1.0 : 0.5; };
    auto thickness = [](int x) { return 200.0 + 10.0 * x; };
    setice(first, firstice, share, thickness);
    setice(second, secondice, share, thickness);
    const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 200.0);
    std::vector<double> firstflux, secondflux;
    const StreamPowerBudget firstbudget = step(first, parameters, &firstice.input, &firstflux);
    const StreamPowerBudget secondbudget = step(second, parameters, &secondice.input, &secondflux);
    const std::size_t bytes = first.elevation.size() * sizeof(float);
    EXPECT_EQ(std::memcmp(first.elevation.data(), second.elevation.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(first.erosion.data(), second.erosion.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(first.deposition.data(), second.deposition.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(firstice.glacialerosion.data(), secondice.glacialerosion.data(), bytes), 0);
    EXPECT_EQ(std::memcmp(firstflux.data(), secondflux.data(), firstflux.size() * sizeof(double)), 0);
    EXPECT_EQ(firstbudget.eroded_km3, secondbudget.eroded_km3);
    EXPECT_EQ(firstbudget.landdeposited_km3, secondbudget.landdeposited_km3);
    EXPECT_EQ(firstbudget.marinedeposited_km3, secondbudget.marinedeposited_km3);
    EXPECT_EQ(firstbudget.depositionsweeps, secondbudget.depositionsweeps);
}

TEST_F(Landscape, StreamPowerGlacialRejectsInvalidInput)
{
    Domain domain = makedomain(1, 100.0f);
    Ice ice;
    setice(domain, ice, [](int) { return 1.0; }, [](int) { return 0.0; });
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    ice.flux.pop_back();
    EXPECT_THROW(step(domain, parameters, &ice.input), std::invalid_argument);
    setice(domain, ice, [](int) { return 1.0; }, [](int) { return 0.0; });
    ice.thickness.pop_back();
    EXPECT_THROW(step(domain, parameters, &ice.input), std::invalid_argument);
    setice(domain, ice, [](int) { return 1.0; }, [](int) { return 0.0; });
    ice.input.thickness_m = nullptr;
    EXPECT_THROW(step(domain, parameters, &ice.input), std::invalid_argument);
}

TEST_F(Landscape, StreamPowerGlacialPartialIceSteadyState)
{
    // gamma = 0.5 everywhere, no thickness: K_eff = 0.5^m x 2 K + 0.5^m x K = 3 / sqrt(2) K
    // (2.1213 K; a linear mix would give 1.5 K), and the glacial share of the incision
    // is 0.5^m x 2 K / K_eff = 2/3.
    const double uplift = 100.0;
    Domain domain = makedomain(1, static_cast<float>(uplift));
    Ice ice;
    setice(domain, ice, [](int) { return 0.5; }, [](int) { return 0.0; });
    const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
    int count = 1;
    for (; count <= 5000; count++)
    {
        const std::vector<float> before = domain.elevation;
        step(domain, parameters, &ice.input);
        double change = 0.0;
        for (std::size_t cell : domain.chain)
            change = std::max(change, std::abs(static_cast<double>(domain.elevation[cell]) - before[cell]));
        if (change < 1.0e-6)
            break;
    }
    ASSERT_LE(count, 5000);
    const double keff = 3.0 / std::sqrt(2.0) * erodibility_k;
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double receiver = x == 0 ? parameters.baselevel_m : bed(domain, x - 1);
        const double discharge = domain.graph.discharge_m3_yr[cell] / 1.0e9;
        const double expected = uplift * domain.graph.distance_m[cell] / (keff * std::pow(discharge, areaexponent));
        EXPECT_NEAR(bed(domain, x) - receiver, expected, 1.0e-3 * expected) << "column " << x;
        ASSERT_GT(domain.erosion[cell], 0.0f) << "column " << x;
        EXPECT_NEAR(static_cast<double>(ice.glacialerosion[cell]) / domain.erosion[cell], 2.0 / 3.0, 1.0e-6)
            << "column " << x;
    }
}

TEST_F(Landscape, StreamPowerGlacialFloorKeepsColumnBalance)
{
    // One 20 Myr step of a glacier 100 m thick on the hill (gamma = 1, G = 0, no
    // diffusion): the outlet's surface would fall below its flotation floor, so the
    // floor binds there (bed -89.2 m), and every chain cell still changes by exactly
    // U dt - erosion + deposition (the floor reduces the incision, it adds no rock).
    const double uplift = 1.0;
    const double dt = 20.0;
    Domain domain = makedomain(4, static_cast<float>(uplift));
    sethill(domain);
    Ice ice;
    setice(domain, ice, [](int) { return 1.0; }, [](int) { return 100.0; });
    const StreamPowerParameters parameters = makeparameters(dt, 0.0, 0.0);
    const std::vector<float> before = domain.elevation;
    const StreamPowerBudget budget = step(domain, parameters, &ice.input);
    for (int x = 0; x < chainlength; x++)
    {
        const std::size_t cell = domain.chain[static_cast<std::size_t>(x)];
        const double change = static_cast<double>(domain.elevation[cell]) - before[cell];
        const double balance = uplift * dt - domain.erosion[cell] + domain.deposition[cell];
        EXPECT_NEAR(change, balance, 2.0e-3) << "column " << x;
    }
    const double floor = parameters.baselevel_m - flotation * 100.0;
    EXPECT_NEAR(bed(domain, 0), floor, 1.0e-3);
    EXPECT_GT(domain.erosion[domain.chain[0]], 0.0f);
    EXPECT_NEAR(budget.eroded_km3, budget.landdeposited_km3 + budget.marinedeposited_km3, 1.0e-9 * budget.eroded_km3);
}

TEST_F(Landscape, StreamPowerGlacialBedBelowFloorStays)
{
    // The outlet's bed starts below its flotation floor (h = 100 m, floor -89.2 m); the
    // hill upstream is ice-free. It is not incised: its erosion is 0 and it changes by
    // its deposition only.
    // (a) Bed -500 m, G = 1: its surface (-400 m) is below the sea, and the sediment
    //     settling on it is its whole change.
    {
        Domain domain = makedomain(4, 0.0f);
        sethill(domain);
        const std::size_t outlet = domain.chain[0];
        domain.elevation[outlet] = -500.0f;
        Ice ice;
        setice(domain, ice, [](int) { return 0.0; }, [](int x) { return x == 0 ? 100.0 : 0.0; });
        const StreamPowerParameters parameters = makeparameters(1.0, 1.0, 0.0);
        const std::vector<float> before = domain.elevation;
        step(domain, parameters, &ice.input);
        EXPECT_EQ(domain.erosion[outlet], 0.0f);
        EXPECT_GT(domain.deposition[outlet], 0.0f);
        EXPECT_LE(domain.elevation[outlet], before[outlet] + domain.deposition[outlet] + 1.0e-3);
        EXPECT_EQ(ice.glacialerosion[outlet], 0.0f);
    }
    // (b) Bed -95 m, G = 0: its surface (+5 m) is above the sea but below the floor's
    //     surface (+10.8 m), so the solve would incise it; it stays where it was.
    {
        Domain domain = makedomain(4, 0.0f);
        sethill(domain);
        const std::size_t outlet = domain.chain[0];
        domain.elevation[outlet] = -95.0f;
        Ice ice;
        setice(domain, ice, [](int) { return 1.0; }, [](int x) { return x == 0 ? 100.0 : 0.0; });
        const StreamPowerParameters parameters = makeparameters(1.0, 0.0, 0.0);
        step(domain, parameters, &ice.input);
        EXPECT_EQ(domain.erosion[outlet], 0.0f);
        EXPECT_EQ(domain.deposition[outlet], 0.0f);
        EXPECT_EQ(domain.elevation[outlet], -95.0f);
        EXPECT_EQ(ice.glacialerosion[outlet], 0.0f);
    }
}

TEST_F(Landscape, GlacialIceFluxCalves)
{
    // The ELA (0 m) below the whole hill: every source is positive, min(z / 500, 1) of
    // runoff x area, so nothing ablates and the outlet carries the sum of the chain's
    // sources, which calves: no sea cell carries ice.
    Domain domain = makedomain(4, 0.0f);
    sethill(domain);
    const std::size_t cells = domain.metric.size();
    const std::vector<float> ela(cells, 0.0f);
    const std::vector<float> runoff(cells, static_cast<float>(runoff_m_yr));
    const double range = 500.0;
    std::vector<double> flux;
    landscape::accumulateice(domain.metric, domain.graph, domain.elevation, ela, runoff, range, flux);
    double sources = 0.0;
    for (std::size_t cell : domain.chain)
        sources += runoff_m_yr * domain.metric.area_m2(cell)
            * std::min((static_cast<double>(domain.elevation[cell]) - ela[cell]) / range, 1.0);
    EXPECT_NEAR(flux[domain.chain[0]], sources, 1.0e-12 * sources);
    // 0.2 + 0.4 + 0.6 + 0.8 + 28 x 1 cell sources.
    EXPECT_NEAR(sources, 30.0 * runoff_m_yr * domain.metric.area_m2(domain.chain[0]), 1.0e-12 * sources);
    for (std::size_t cell = 0; cell < cells; cell++)
        if (domain.graph.land[cell] == 0)
            EXPECT_EQ(flux[cell], 0.0) << "sea cell " << cell;

    FlowGraph broken = domain.graph;
    broken.receiver[domain.chain[5]] = static_cast<int32_t>(cells);
    EXPECT_THROW(landscape::accumulateice(domain.metric, broken, domain.elevation, ela, runoff, range, flux),
                 std::invalid_argument);
    broken = domain.graph;
    broken.order.push_back(-1);
    EXPECT_THROW(landscape::accumulateice(domain.metric, broken, domain.elevation, ela, runoff, range, flux),
                 std::invalid_argument);
}
