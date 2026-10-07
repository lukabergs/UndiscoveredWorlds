// Landscape.LakeBalance*: lake levels by water balance on synthetic depressions
// (docs/architecture/erosion-hydrology-plan.md, Phase B, step B2).

#include "landscape_fixture.hpp"
#include "simulations/landscape/lake_balance.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

using landscape::LakeBalance;
using landscape::LakeBalanceInputs;
using landscape::LakeLevel;
using landscape::SphereMetric;

namespace
{
constexpr double earthradius_m = 6371000.0;
constexpr int bowlx = 10;
constexpr int bowly = 6;
constexpr float spill_m = 200.0f;
constexpr float plain_m = 300.0f;
constexpr double ignoreddischarge_m3_yr = 1.0e12;  // on every cell that is not an external donor
// The bowl's terrain, row by row: distinct, all below the spill; the lowest is (0, 1).
constexpr float bowlfloor_m[9] = { 150.0f, 120.0f, 180.0f, 110.0f, 140.0f, 170.0f, 130.0f, 160.0f, 190.0f };

SphereMetric makemetric() { return SphereMetric(32, 16, earthradius_m); }

// A plain at 300 m with no depression, no receivers, no runoff and P = E.
LakeBalanceInputs makeplain(const SphereMetric& metric)
{
    const std::size_t size = metric.size();
    LakeBalanceInputs inputs;
    inputs.metric = &metric;
    inputs.elevation.assign(size, plain_m);
    inputs.filledelevation.assign(size, plain_m);
    inputs.depressionid.assign(size, 0);
    inputs.receiver.assign(size, -1);
    inputs.discharge_m3_yr.assign(size, ignoreddischarge_m3_yr);
    inputs.runoff_m_yr.assign(size, 0.0f);
    inputs.precipitation_m_yr.assign(size, 0.5f);
    inputs.evaporation_m_yr.assign(size, 0.5f);
    return inputs;
}

std::size_t lowestbowlcell(const SphereMetric& metric, int x0, int y0) { return metric.index(x0, y0 + 1); }

// A 3 x 3 bowl at (x0, y0) with id `id`: every cell drains to the lowest one (internal
// flow, which brings no inflow), and the cell west of the lowest drains into it with
// `donor_m3_yr` (the one external donor).
void addbowl(LakeBalanceInputs& inputs, const SphereMetric& metric, int x0, int y0, int32_t id, double donor_m3_yr)
{
    const std::size_t lowest = lowestbowlcell(metric, x0, y0);
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 3; i++)
        {
            const std::size_t cell = metric.index(x0 + i, y0 + j);
            inputs.elevation[cell] = bowlfloor_m[j * 3 + i];
            inputs.filledelevation[cell] = spill_m;
            inputs.depressionid[cell] = id;
            inputs.receiver[cell] = static_cast<int32_t>(lowest);
        }
    }
    const std::size_t donor = metric.index(x0 - 1, y0 + 1);
    inputs.receiver[donor] = static_cast<int32_t>(lowest);
    inputs.discharge_m3_yr[donor] = donor_m3_yr;
}

// The bowl's cells by terrain ascending.
std::vector<std::size_t> bowlcells(const LakeBalanceInputs& inputs, int32_t id)
{
    std::vector<std::size_t> cells;
    for (std::size_t cell = 0; cell < inputs.depressionid.size(); cell++)
        if (inputs.depressionid[cell] == id)
            cells.push_back(cell);
    std::sort(cells.begin(), cells.end(), [&](std::size_t a, std::size_t b) {
        return inputs.elevation[a] < inputs.elevation[b];
    });
    return cells;
}

double areaof(const SphereMetric& metric, const std::vector<std::size_t>& cells, std::size_t count)
{
    double area = 0.0;
    for (std::size_t index = 0; index < count; index++)
        area += metric.area_m2(cells[index]);
    return area;
}

// The area-weighted mean of E - P over the cells.
double netloss(const LakeBalanceInputs& inputs, const std::vector<std::size_t>& cells)
{
    double area = 0.0;
    double precipitation = 0.0;
    double evaporation = 0.0;
    for (std::size_t cell : cells)
    {
        const double a = inputs.metric->area_m2(cell);
        area += a;
        precipitation += static_cast<double>(inputs.precipitation_m_yr[cell]) * a;
        evaporation += static_cast<double>(inputs.evaporation_m_yr[cell]) * a;
    }
    return std::max(0.0, evaporation / area - precipitation / area);
}

void setuniformclimate(LakeBalanceInputs& inputs, float runoff, float precipitation, float evaporation)
{
    std::fill(inputs.runoff_m_yr.begin(), inputs.runoff_m_yr.end(), runoff);
    std::fill(inputs.precipitation_m_yr.begin(), inputs.precipitation_m_yr.end(), precipitation);
    std::fill(inputs.evaporation_m_yr.begin(), inputs.evaporation_m_yr.end(), evaporation);
}
}

TEST_F(Landscape, LakeBalanceClosedBowlReachesItsLevel)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    constexpr float runoff = 0.01f;
    setuniformclimate(inputs, runoff, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    const double loss = netloss(inputs, cells);
    const double depressionarea = areaof(metric, cells, cells.size());
    // A* halfway into the fourth lowest cell: the lake covers the four lowest cells.
    const double equilibrium = areaof(metric, cells, 3) + 0.5 * metric.area_m2(cells[3]);
    const double ownrunoff = static_cast<double>(runoff) * depressionarea;
    const double donor = equilibrium * loss - ownrunoff;
    ASSERT_GT(donor, 0.0);
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = donor;

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    const LakeLevel& lake = balance.lakes[0];
    const float level = 0.5f * (inputs.elevation[cells[3]] + inputs.elevation[cells[4]]);
    EXPECT_EQ(lake.id, 7);
    EXPECT_EQ(lake.spill_m, spill_m);
    EXPECT_FALSE(lake.overflowing);
    EXPECT_EQ(lake.lakecells, 4);
    EXPECT_EQ(lake.level_m, level);
    EXPECT_EQ(lake.lowestcell, static_cast<int32_t>(lowestbowlcell(metric, bowlx, bowly)));
    EXPECT_NEAR(lake.inflow_m3_yr, donor + ownrunoff, 1.0e-9 * (donor + ownrunoff));
    EXPECT_NEAR(lake.netloss_m_yr, loss, 1.0e-12);
    EXPECT_NEAR(lake.area_m2, areaof(metric, cells, 4), 1.0e-9 * depressionarea);
    EXPECT_NEAR(lake.depressionarea_m2, depressionarea, 1.0e-9 * depressionarea);
    ASSERT_EQ(balance.surface.size(), metric.size());
    for (std::size_t index = 0; index < cells.size(); index++)
        EXPECT_EQ(balance.surface[cells[index]], index < 4 ? level : 0.0f) << "bowl cell " << index;
    std::size_t lakecells = 0;
    for (float surface : balance.surface)
        lakecells += surface != 0.0f ? 1 : 0;
    EXPECT_EQ(lakecells, 4u);
}

TEST_F(Landscape, LakeBalanceLargeInflowOverflows)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.01f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    const double depressionarea = areaof(metric, cells, cells.size());
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = 2.0 * depressionarea * netloss(inputs, cells);

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    const LakeLevel& lake = balance.lakes[0];
    EXPECT_TRUE(lake.overflowing);
    EXPECT_EQ(lake.level_m, spill_m);
    EXPECT_EQ(lake.lakecells, 9);
    EXPECT_NEAR(lake.area_m2, depressionarea, 1.0e-9 * depressionarea);
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        EXPECT_EQ(balance.surface[cell], inputs.depressionid[cell] == 7 ? spill_m : 0.0f) << "cell " << cell;
}

TEST_F(Landscape, LakeBalanceDryDepressionHoldsNoLake)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    // A* is half the lowest cell's area.
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = 0.5 * metric.area_m2(cells[0]) * netloss(inputs, cells);

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    const LakeLevel& lake = balance.lakes[0];
    EXPECT_FALSE(lake.overflowing);
    EXPECT_EQ(lake.lakecells, 0);
    EXPECT_EQ(lake.level_m, 0.0f);
    EXPECT_EQ(lake.area_m2, 0.0);
    EXPECT_EQ(lake.lowestcell, static_cast<int32_t>(cells[0]));
    EXPECT_EQ(lake.lowestcell, static_cast<int32_t>(lowestbowlcell(metric, bowlx, bowly)));
    for (float surface : balance.surface)
        EXPECT_EQ(surface, 0.0f);
}

TEST_F(Landscape, LakeBalanceNoNetLossOverflows)
{
    const SphereMetric metric = makemetric();
    for (const float evaporation : { 0.3f, 0.1f })  // E = P, then E < P; no inflow at all
    {
        LakeBalanceInputs inputs = makeplain(metric);
        setuniformclimate(inputs, 0.0f, 0.3f, evaporation);
        addbowl(inputs, metric, bowlx, bowly, 7, 0.0);

        const LakeBalance balance = landscape::balancelakes(inputs);

        ASSERT_EQ(balance.lakes.size(), 1u);
        const LakeLevel& lake = balance.lakes[0];
        EXPECT_EQ(lake.inflow_m3_yr, 0.0);
        EXPECT_EQ(lake.netloss_m_yr, 0.0);
        EXPECT_TRUE(lake.overflowing) << "E = " << evaporation;
        EXPECT_EQ(lake.level_m, spill_m);
        EXPECT_EQ(lake.lakecells, 9);
        for (std::size_t cell = 0; cell < metric.size(); cell++)
            EXPECT_EQ(balance.surface[cell], inputs.depressionid[cell] == 7 ? spill_m : 0.0f) << "cell " << cell;
    }
}

TEST_F(Landscape, LakeBalanceOneLevelPerDepressionInAscendingId)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    // Id 9 is found first in cell order, id 3 later.
    addbowl(inputs, metric, 20, 3, 9, 0.0);
    addbowl(inputs, metric, 5, 9, 3, 0.0);
    const std::vector<std::size_t> cells3 = bowlcells(inputs, 3);
    const std::vector<std::size_t> cells9 = bowlcells(inputs, 9);
    // Id 3 gets an overflowing inflow, id 9 none (dry).
    const double area3 = areaof(metric, cells3, cells3.size());
    inputs.discharge_m3_yr[metric.index(4, 10)] = 2.0 * area3 * netloss(inputs, cells3);
    inputs.discharge_m3_yr[metric.index(19, 4)] = 0.0;

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 2u);
    EXPECT_EQ(balance.lakes[0].id, 3);
    EXPECT_EQ(balance.lakes[1].id, 9);
    EXPECT_TRUE(balance.lakes[0].overflowing);
    EXPECT_FALSE(balance.lakes[1].overflowing);
    EXPECT_EQ(balance.lakes[1].lakecells, 0);
    EXPECT_EQ(balance.lakes[1].lowestcell, static_cast<int32_t>(cells9[0]));
    EXPECT_NEAR(balance.lakes[0].depressionarea_m2, area3, 1.0e-9 * area3);
    const double area9 = areaof(metric, cells9, cells9.size());
    EXPECT_NEAR(balance.lakes[1].depressionarea_m2, area9, 1.0e-9 * area9);
}

TEST_F(Landscape, LakeBalanceIsDeterministic)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.01f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    addbowl(inputs, metric, 20, 3, 2, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = areaof(metric, cells, 5) * netloss(inputs, cells);

    const LakeBalance first = landscape::balancelakes(inputs);
    const LakeBalance second = landscape::balancelakes(inputs);

    ASSERT_EQ(first.surface.size(), second.surface.size());
    EXPECT_EQ(std::memcmp(first.surface.data(), second.surface.data(), first.surface.size() * sizeof(float)), 0);
    ASSERT_EQ(first.lakes.size(), second.lakes.size());
    for (std::size_t index = 0; index < first.lakes.size(); index++)
    {
        EXPECT_EQ(first.lakes[index].id, second.lakes[index].id);
        EXPECT_EQ(first.lakes[index].level_m, second.lakes[index].level_m);
        EXPECT_EQ(first.lakes[index].lakecells, second.lakes[index].lakecells);
        EXPECT_EQ(first.lakes[index].inflow_m3_yr, second.lakes[index].inflow_m3_yr);
    }
}

TEST_F(Landscape, LakeBalanceRejectsSizeMismatch)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    inputs.evaporation_m_yr.pop_back();
    EXPECT_THROW(landscape::balancelakes(inputs), std::invalid_argument);
}

TEST_F(Landscape, LakeBalanceClosedLakeCoveringEveryCell)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 4, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 4);
    const double loss = netloss(inputs, cells);
    // A* halfway into the highest cell: every cell is lake, the lake stays closed.
    const double equilibrium = areaof(metric, cells, 8) + 0.5 * metric.area_m2(cells[8]);
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = equilibrium * loss;

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    const LakeLevel& lake = balance.lakes[0];
    EXPECT_FALSE(lake.overflowing);
    EXPECT_EQ(lake.lakecells, 9);
    EXPECT_EQ(lake.level_m, 0.5f * (inputs.elevation[cells[8]] + spill_m));
}

TEST_F(Landscape, LakeBalanceTiedCellsFloodTogether)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 5, 0.0);
    std::vector<std::size_t> cells = bowlcells(inputs, 5);
    // The fifth lowest cell takes the fourth lowest's terrain.
    inputs.elevation[cells[4]] = inputs.elevation[cells[3]];
    cells = bowlcells(inputs, 5);
    const double loss = netloss(inputs, cells);
    const double equilibrium = areaof(metric, cells, 3) + 0.5 * metric.area_m2(cells[3]);
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = equilibrium * loss;

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    const LakeLevel& lake = balance.lakes[0];
    EXPECT_EQ(lake.lakecells, 5);
    EXPECT_EQ(lake.level_m, 0.5f * (inputs.elevation[cells[4]] + inputs.elevation[cells[5]]));
}

TEST_F(Landscape, LakeBalanceFlatDownstreamIsNotADonor)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 6, 0.0);
    const std::size_t lowest = lowestbowlcell(metric, bowlx, bowly);
    // A shallow flooded cell east of the bowl at the spill level, id 0: the bowl's lowest
    // cell drains to it and it drains back into the bowl with a huge discharge.
    const std::size_t shallow = metric.index(bowlx + 3, bowly + 1);
    inputs.elevation[shallow] = spill_m - 0.05f;
    inputs.filledelevation[shallow] = spill_m;
    inputs.receiver[lowest] = static_cast<int32_t>(shallow);
    inputs.receiver[shallow] = static_cast<int32_t>(metric.index(bowlx + 2, bowly + 1));
    inputs.discharge_m3_yr[shallow] = 1.0e15;

    const LakeBalance balance = landscape::balancelakes(inputs);

    ASSERT_EQ(balance.lakes.size(), 1u);
    EXPECT_EQ(balance.lakes[0].inflow_m3_yr, 0.0);
    EXPECT_EQ(balance.lakes[0].lakecells, 0);
}

namespace
{
// The cells at or below `sealevel` in whole metres (the rule of planet::basebelowsea).
std::vector<uint8_t> belowsea(const LakeBalanceInputs& inputs, float sealevel)
{
    std::vector<uint8_t> below(inputs.elevation.size(), 0);
    for (std::size_t cell = 0; cell < below.size(); cell++)
        below[cell] = static_cast<float>(std::lround(inputs.elevation[cell])) <= sealevel ? 1 : 0;
    return below;
}

// The bowl of LakeBalanceClosedBowlReachesItsLevel: a closed lake on its four lowest cells.
LakeBalanceInputs closedbowl(const SphereMetric& metric)
{
    LakeBalanceInputs inputs = makeplain(metric);
    constexpr float runoff = 0.01f;
    setuniformclimate(inputs, runoff, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    const double equilibrium = areaof(metric, cells, 3) + 0.5 * metric.area_m2(cells[3]);
    const double ownrunoff = static_cast<double>(runoff) * areaof(metric, cells, cells.size());
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] = equilibrium * netloss(inputs, cells) - ownrunoff;
    return inputs;
}
}

// An enclosed below-sea bowl (erosion plan H3; sea level 175 m, spill 200 m, floor 110 to
// 190 m) with a closed lake on its four lowest cells (110 to 140 m): its other cells at or
// below the sea level (150, 160, 170 m) are salt pans, the two above it (180, 190 m) and
// the lake are not. A second pass on the same inputs gives the same lake and pans.
TEST_F(Landscape, LakeBalanceSaltFlatsAroundAClosedLakeBelowTheSea)
{
    const SphereMetric metric = makemetric();
    const LakeBalanceInputs inputs = closedbowl(metric);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    const std::vector<uint8_t> below = belowsea(inputs, 175.0f);

    const LakeBalance balance = landscape::balancelakes(inputs);
    const std::vector<uint8_t> pans = landscape::saltpans(inputs, balance, below);

    ASSERT_EQ(balance.lakes.size(), 1u);
    EXPECT_FALSE(balance.lakes[0].overflowing);
    EXPECT_EQ(balance.lakes[0].lakecells, 4);
    const float level = 0.5f * (inputs.elevation[cells[3]] + inputs.elevation[cells[4]]);
    EXPECT_LT(level, 175.0f);  // the lake surface lies below the sea level
    for (std::size_t index = 0; index < cells.size(); index++)
    {
        EXPECT_EQ(balance.surface[cells[index]], index < 4 ? level : 0.0f) << "bowl cell " << index;
        EXPECT_EQ(pans[cells[index]], index >= 4 && index <= 6 ? 1 : 0) << "bowl cell " << index;
    }
    std::size_t pancount = 0;
    for (uint8_t pan : pans)
        pancount += pan;
    EXPECT_EQ(pancount, 3u);

    const LakeBalance again = landscape::balancelakes(inputs);
    EXPECT_EQ(again.surface, balance.surface);
    EXPECT_EQ(landscape::saltpans(inputs, again, below), pans);
}

// A dry bowl below the sea: its lowest cell and every cell at or below the sea level are pans.
TEST_F(Landscape, LakeBalanceDryBowlBelowTheSeaIsAPan)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.0f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    const LakeBalance balance = landscape::balancelakes(inputs);
    ASSERT_EQ(balance.lakes[0].lakecells, 0);
    const std::vector<uint8_t> pans = landscape::saltpans(inputs, balance, belowsea(inputs, 175.0f));
    for (std::size_t index = 0; index < cells.size(); index++)
        EXPECT_EQ(pans[cells[index]], index <= 6 ? 1 : 0) << "bowl cell " << index;
    // Above the sea only the lowest cell.
    const std::vector<uint8_t> high = landscape::saltpans(inputs, balance, belowsea(inputs, 100.0f));
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        EXPECT_EQ(high[cell], cell == cells[0] ? 1 : 0) << "cell " << cell;
}

// An overflowing (open) depression gets no pans, above the sea or not.
TEST_F(Landscape, LakeBalanceOpenDepressionHasNoPans)
{
    const SphereMetric metric = makemetric();
    LakeBalanceInputs inputs = makeplain(metric);
    setuniformclimate(inputs, 0.01f, 0.3f, 1.2f);
    addbowl(inputs, metric, bowlx, bowly, 7, 0.0);
    const std::vector<std::size_t> cells = bowlcells(inputs, 7);
    inputs.discharge_m3_yr[metric.index(bowlx - 1, bowly + 1)] =
        2.0 * areaof(metric, cells, cells.size()) * netloss(inputs, cells);
    const LakeBalance balance = landscape::balancelakes(inputs);
    ASSERT_TRUE(balance.lakes[0].overflowing);
    for (float sealevel : { 100.0f, 175.0f })
        for (uint8_t pan : landscape::saltpans(inputs, balance, belowsea(inputs, sealevel)))
            EXPECT_EQ(pan, 0);
}
