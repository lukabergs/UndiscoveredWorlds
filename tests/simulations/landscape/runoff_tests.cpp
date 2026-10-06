// Landscape.Runoff*: the Budyko runoff and the transmission losses on synthetic grids
// (docs/architecture/erosion-hydrology-plan.md, Phase C, step C2).

#include "landscape_fixture.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/runoff.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr double fuw = 2.6;

// A hand-built chain in row 1 of an 8 x 4 grid on a 50 km sphere (cells of about
// 1.4e9 m2, so the losses are a sizeable share of the flow): column 0 is sea,
// columns 3 -> 2 -> 1 drain west and column 1 is the outlet into column 0. Every
// link is 100 km long; the elevations 100, 10 and 1 m over a sea level of 0 give
// slopes of 9e-4 (column 3), 9e-5 (column 2) and 1e-5 (column 1).
struct Chain
{
    landscape::SphereMetric metric{8, 4, 50000.0};
    landscape::FlowGraph graph;
    std::vector<float> elevation;
    std::size_t sea = 0, outlet = 0, middle = 0, head = 0;
    double length = 100000.0;
    double area = 0.0;

    Chain()
    {
        const std::size_t size = metric.size();
        sea = metric.index(0, 1);
        outlet = metric.index(1, 1);
        middle = metric.index(2, 1);
        head = metric.index(3, 1);
        area = metric.area_m2(outlet);
        elevation.assign(size, -100.0f);
        elevation[outlet] = 1.0f;
        elevation[middle] = 10.0f;
        elevation[head] = 100.0f;
        graph.receiver.resize(size);
        for (std::size_t cell = 0; cell < size; cell++)
            graph.receiver[cell] = static_cast<int32_t>(cell);
        graph.receiver[middle] = static_cast<int32_t>(outlet);
        graph.receiver[head] = static_cast<int32_t>(middle);
        graph.direction.assign(size, 0);
        graph.direction[outlet] = graph.direction[middle] = graph.direction[head] = 7;
        graph.distance_m.assign(size, 0.0f);
        graph.distance_m[outlet] = graph.distance_m[middle] = graph.distance_m[head] = static_cast<float>(length);
        graph.outletsea.assign(size, -1);
        graph.outletsea[outlet] = static_cast<int32_t>(sea);
        graph.land.assign(size, 0);
        graph.land[outlet] = graph.land[middle] = graph.land[head] = 1;
        graph.order = { static_cast<int32_t>(outlet), static_cast<int32_t>(middle), static_cast<int32_t>(head) };
    }

    std::vector<float> uniform(float value) const
    {
        std::vector<float> field(metric.size(), 0.0f);
        field[outlet] = field[middle] = field[head] = value;
        return field;
    }
};

landscape::TransmissionLoss chainloss(const Chain& chain, float deficit)
{
    landscape::TransmissionLoss loss;
    loss.deficit_m_yr = chain.uniform(deficit);
    loss.openwater.assign(chain.metric.size(), 0);
    loss.widthcoefficient = 5.0;
    loss.widthexponent = 0.5;
    loss.floodplainfactor = 10.0;
    loss.floodplainslope = 2.0e-4;
    loss.sealevel = 0.0f;
    return loss;
}

// Column 0 is sea at -100 m; land rises 100 m per column, with a 3 x 3 bowl at
// columns 6..8, rows 3..5 lowered by 500 m (flood receivers), on a 16 x 8 grid.
struct Terrain
{
    landscape::SphereMetric metric{16, 8, 6371000.0};
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
    landscape::FlowGraph graph;
    std::vector<float> runoff;

    Terrain()
    {
        elevation.resize(metric.size());
        sea.resize(metric.size());
        runoff.resize(metric.size());
        for (int y = 0; y < metric.height; y++)
            for (int x = 0; x < metric.width; x++)
            {
                const std::size_t cell = metric.index(x, y);
                sea[cell] = x == 0 ? 1 : 0;
                elevation[cell] = x == 0 ? -100.0f : 100.0f * static_cast<float>(x);
                if (x >= 6 && x <= 8 && y >= 3 && y <= 5)
                    elevation[cell] -= 500.0f;
                runoff[cell] = 0.1f + 0.013f * static_cast<float>((x * 7 + y * 3) % 11);
            }
        const landscape::Depressions depressions = landscape::filldepressions(metric, elevation, sea, 0.1f);
        graph = landscape::buildflowgraph(metric, elevation, sea, 0.0f, depressions);
    }
};
}

TEST_F(Landscape, RunoffBudykoFuValues)
{
    // runoff = (P^w + PET^w)^(1/w) - PET with w = 2.6; e.g. P = PET = 1: 2^(1/2.6) - 1 =
    // exp(0.2665950) - 1 = 0.3055117.
    EXPECT_NEAR(landscape::budykorunoff(1.0, 1.0, fuw), 0.30551169770980957, 1e-12);
    EXPECT_NEAR(landscape::budykorunoff(1.0, 0.5, fuw), 0.5604767505428636, 1e-12);
    EXPECT_NEAR(landscape::budykorunoff(1.0, 2.0, fuw), 0.12095350108572678, 1e-12);
    EXPECT_NEAR(landscape::budykorunoff(0.5, 2.0, fuw), 0.020754046842914597, 1e-12);
    EXPECT_NEAR(landscape::budykorunoff(0.2, 2.0, fuw), 0.0019307289583454512, 1e-12);
    EXPECT_NEAR(landscape::budykorunoff(2.0, 0.1, fuw), 1.9003186568834893, 1e-12);
}

TEST_F(Landscape, RunoffBudykoLimits)
{
    EXPECT_EQ(landscape::budykorunoff(0.8, 0.0, fuw), 0.8);  // no evaporation: all the precipitation
    EXPECT_EQ(landscape::budykorunoff(0.0, 1.5, fuw), 0.0);  // no precipitation: no runoff
    EXPECT_EQ(landscape::budykorunoff(-1.0, 1.5, fuw), 0.0);
    // Far into the arid limit the runoff tends to 0 and stays non-negative; at PET << P it tends to P - PET.
    const double arid = landscape::budykorunoff(0.01, 10.0, fuw);
    EXPECT_GT(arid, 0.0);
    EXPECT_LT(arid, 1.0e-6);
    EXPECT_NEAR(landscape::budykorunoff(10.0, 0.01, fuw), 10.0 - 0.01, 1.0e-6);
    EXPECT_THROW(landscape::budykorunoff(1.0, 1.0, 1.0), std::invalid_argument);
}

TEST_F(Landscape, RunoffLossesChainExactDischarge)
{
    Chain chain;
    const std::vector<float> runoff = chain.uniform(0.5f);
    const landscape::TransmissionLoss loss = chainloss(chain, 0.2f);
    const double total = landscape::accumulateflowwithlosses(chain.metric, runoff, chain.elevation, loss, chain.graph);

    // Head: slope 9e-4, a channel; middle and outlet: slopes below 2e-4, floodplains.
    const double r = static_cast<double>(0.5f);
    const double d = static_cast<double>(0.2f);
    const double qhead = r * chain.area;
    const double losshead = std::min(qhead, 5.0 * std::sqrt(qhead / secondsperyear) * chain.length * d);
    const double qmiddle = (qhead - losshead) + r * chain.area;
    const double lossmiddle = std::min(qmiddle, 10.0 * 5.0 * std::sqrt(qmiddle / secondsperyear) * chain.length * d);
    const double qoutlet = (qmiddle - lossmiddle) + r * chain.area;
    const double lossoutlet = std::min(qoutlet, 10.0 * 5.0 * std::sqrt(qoutlet / secondsperyear) * chain.length * d);

    EXPECT_GT(losshead, 0.0);
    EXPECT_LT(lossoutlet, qoutlet);  // the losses are partial here, not capped
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.head], qhead - losshead);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.middle], qmiddle - lossmiddle);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.outlet], qoutlet - lossoutlet);
    EXPECT_DOUBLE_EQ(total, losshead + lossmiddle + lossoutlet);
    // Water is conserved: the runoff either leaves at the outlet or is lost.
    EXPECT_NEAR(chain.graph.discharge_m3_yr[chain.outlet] + total, 3.0 * r * chain.area, 1e-6 * chain.area);
    // Drainage area is unaffected by the losses.
    EXPECT_DOUBLE_EQ(chain.graph.area_m2[chain.outlet], 3.0 * chain.area);
    EXPECT_EQ(chain.graph.discharge_m3_yr[chain.sea], 0.0);
}

TEST_F(Landscape, RunoffLossesCappedAtFlow)
{
    Chain chain;
    const std::vector<float> runoff = chain.uniform(0.5f);
    const landscape::TransmissionLoss loss = chainloss(chain, 1.0e6f);
    const double total = landscape::accumulateflowwithlosses(chain.metric, runoff, chain.elevation, loss, chain.graph);
    for (std::size_t cell : { chain.head, chain.middle, chain.outlet })
        EXPECT_EQ(chain.graph.discharge_m3_yr[cell], 0.0) << "cell " << cell;
    EXPECT_DOUBLE_EQ(total, 3.0 * static_cast<double>(0.5f) * chain.area);
}

TEST_F(Landscape, RunoffLossesOpenWaterLosesDeficitTimesArea)
{
    Chain chain;
    const std::vector<float> runoff = chain.uniform(0.5f);
    landscape::TransmissionLoss loss = chainloss(chain, 0.0f);
    loss.deficit_m_yr[chain.middle] = 0.2f;
    loss.openwater[chain.middle] = 1;
    const double total = landscape::accumulateflowwithlosses(chain.metric, runoff, chain.elevation, loss, chain.graph);
    const double r = static_cast<double>(0.5f);
    const double lake = static_cast<double>(0.2f) * chain.area;
    EXPECT_DOUBLE_EQ(total, lake);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.head], r * chain.area);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.middle], 2.0 * r * chain.area - lake);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.outlet], 3.0 * r * chain.area - lake);
}

TEST_F(Landscape, RunoffLossesZeroDeficitReproducesAccumulateFlow)
{
    Terrain terrain;
    landscape::FlowGraph plain = terrain.graph;
    landscape::accumulateflow(terrain.metric, terrain.runoff, plain);

    landscape::TransmissionLoss loss;
    loss.deficit_m_yr.assign(terrain.metric.size(), 0.0f);
    loss.widthcoefficient = 5.0;
    loss.floodplainfactor = 10.0;
    loss.floodplainslope = 2.0e-4;
    landscape::FlowGraph lossy = terrain.graph;
    const double total = landscape::accumulateflowwithlosses(terrain.metric, terrain.runoff, terrain.elevation, loss, lossy);
    EXPECT_EQ(total, 0.0);
    ASSERT_EQ(lossy.discharge_m3_yr.size(), plain.discharge_m3_yr.size());
    for (std::size_t cell = 0; cell < plain.size(); cell++)
    {
        EXPECT_EQ(lossy.discharge_m3_yr[cell], plain.discharge_m3_yr[cell]) << "cell " << cell;
        EXPECT_EQ(lossy.area_m2[cell], plain.area_m2[cell]) << "cell " << cell;
    }
}

TEST_F(Landscape, RunoffLossesDeterministic)
{
    Terrain terrain;
    landscape::TransmissionLoss loss;
    loss.deficit_m_yr.resize(terrain.metric.size());
    loss.openwater.assign(terrain.metric.size(), 0);
    for (std::size_t cell = 0; cell < loss.deficit_m_yr.size(); cell++)
        loss.deficit_m_yr[cell] = 0.05f * static_cast<float>(cell % 5);
    loss.openwater[terrain.metric.index(7, 4)] = 1;  // the bowl's centre
    loss.widthcoefficient = 5.0;
    loss.floodplainfactor = 10.0;
    loss.floodplainslope = 2.0e-4;

    landscape::FlowGraph first = terrain.graph;
    landscape::FlowGraph second = terrain.graph;
    const double firsttotal = landscape::accumulateflowwithlosses(terrain.metric, terrain.runoff, terrain.elevation, loss, first);
    const double secondtotal = landscape::accumulateflowwithlosses(terrain.metric, terrain.runoff, terrain.elevation, loss, second);
    EXPECT_GT(firsttotal, 0.0);
    EXPECT_EQ(firsttotal, secondtotal);
    EXPECT_EQ(first.discharge_m3_yr, second.discharge_m3_yr);
    for (double value : first.discharge_m3_yr)
        EXPECT_GE(value, 0.0);
}

TEST_F(Landscape, RunoffOpenWaterEvaporationLaw)
{
    // E = 0.075 x max(0, T + 2) m/yr (landscape.lakes): 2.025 at 25 C, 0.9 at 10 C, 0.15 at 0 C, 0 at -2 C and below.
    EXPECT_DOUBLE_EQ(landscape::openwaterevaporation(25.0, 0.075, -2.0), 0.075 * 27.0);
    EXPECT_DOUBLE_EQ(landscape::openwaterevaporation(10.0, 0.075, -2.0), 0.075 * 12.0);
    EXPECT_DOUBLE_EQ(landscape::openwaterevaporation(0.0, 0.075, -2.0), 0.075 * 2.0);
    EXPECT_EQ(landscape::openwaterevaporation(-2.0, 0.075, -2.0), 0.0);
    EXPECT_EQ(landscape::openwaterevaporation(-30.0, 0.075, -2.0), 0.0);
    // Applied per season and averaged, a cold winter does not cancel a warm summer:
    // seasons -20, 0, 20, 0 C give (0 + 0.15 + 1.65 + 0.15) / 4 = 0.4875 m/yr, the annual
    // mean temperature (0 C) only 0.15.
    double mean = 0.0;
    for (double temperature : { -20.0, 0.0, 20.0, 0.0 })
        mean += landscape::openwaterevaporation(temperature, 0.075, -2.0);
    EXPECT_NEAR(mean / 4.0, 0.4875, 1e-12);
}

TEST_F(Landscape, RunoffLossesSkipDepressionCells)
{
    Chain chain;
    const std::vector<float> runoff = chain.uniform(0.5f);
    landscape::TransmissionLoss loss = chainloss(chain, 0.2f);
    loss.depressionid.assign(chain.metric.size(), 0);
    loss.depressionid[chain.middle] = 3;  // flooded: the lake balance owns its evaporation
    loss.openwater[chain.middle] = 1;     // even as open water it loses nothing here
    const double total = landscape::accumulateflowwithlosses(chain.metric, runoff, chain.elevation, loss, chain.graph);

    const double r = static_cast<double>(0.5f);
    const double d = static_cast<double>(0.2f);
    const double qhead = r * chain.area;
    const double losshead = std::min(qhead, 5.0 * std::sqrt(qhead / secondsperyear) * chain.length * d);
    const double qmiddle = (qhead - losshead) + r * chain.area;
    const double qoutlet = qmiddle + r * chain.area;
    const double lossoutlet = std::min(qoutlet, 10.0 * 5.0 * std::sqrt(qoutlet / secondsperyear) * chain.length * d);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.middle], qmiddle);
    EXPECT_DOUBLE_EQ(chain.graph.discharge_m3_yr[chain.outlet], qoutlet - lossoutlet);
    EXPECT_DOUBLE_EQ(total, losshead + lossoutlet);
}
