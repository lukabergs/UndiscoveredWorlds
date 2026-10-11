// Tests of landscape/analytical_relief.hpp (erosion plan Phase E, step E2c), suite Landscape:
// Landscape.AnalyticalRelief* on hand-built patch metrics (a chain of N x 1 cells, which
// PatchMetric allows: a neighbour beyond the edge is the cell itself, and square grids).

#include "landscape_fixture.hpp"
#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/regional_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
using landscape::AnalyticalReliefInput;
using landscape::AnalyticalReliefParameters;
using landscape::AnalyticalReliefResult;
using landscape::PatchMetric;

// A W x H patch of square cells of side cell_m (dx = dy, area cell^2).
PatchMetric squaremetric(int width, int height, double cell_m)
{
    PatchMetric metric;
    metric.width = width;
    metric.height = height;
    metric.radius_m = 6.371e6;
    metric.dy_m = cell_m;
    metric.rowdx_m.assign(static_cast<std::size_t>(height), cell_m);
    metric.rowarea_m2.assign(static_cast<std::size_t>(height), cell_m * cell_m);
    metric.rowlatitude_rad.assign(static_cast<std::size_t>(height), 0.0);
    return metric;
}

// A deterministic value in [0, 1) per cell.
double wobble(int x, int y)
{
    std::uint32_t h = static_cast<std::uint32_t>(x) * 73856093u ^ static_cast<std::uint32_t>(y) * 19349663u;
    h ^= h >> 13;
    h *= 0x5bd1e995u;
    h ^= h >> 15;
    return static_cast<double>(h % 1000u) / 1000.0;
}

// The outer ring of a patch.
std::vector<uint8_t> ring(const PatchMetric& metric)
{
    std::vector<uint8_t> pinned(metric.size(), 0);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
            if (x == 0 || y == 0 || x == metric.width - 1 || y == metric.height - 1)
                pinned[metric.index(x, y)] = 1;
    return pinned;
}

// Follows the directions from every unpinned cell; true when each reaches a pinned
// cell within size() links (so no loop) through valid neighbours.
bool everypathreachesapinnedcell(const PatchMetric& metric, const std::vector<uint8_t>& pinned,
                                 const std::vector<uint8_t>& direction)
{
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        std::size_t current = cell;
        std::size_t steps = 0;
        while (pinned[current] == 0)
        {
            const int code = direction[current];
            if (code < 1 || code > 8)
                return false;
            const std::size_t next = metric.neighbour(current, code);
            if (next == current || ++steps > metric.size())
                return false;
            current = next;
        }
        if (pinned[cell] != 0 && direction[cell] != 0)
            return false;
    }
    return true;
}

// The pinned outlet at x = 0 of a chain of N cells 1000 m apart, uniform runoff
// 1 m/yr on cells of 1e6 m^2, start surface 100 m + 1 m per cell (slope 0.001).
// With `branch`, a second row 50 m above the chain without runoff, draining north
// into it.
constexpr int chainlength = 6;
constexpr double chainspacing = 1000.0;
constexpr double chainbase = 100.0;

AnalyticalReliefInput chaininput(bool branch = false)
{
    AnalyticalReliefInput input;
    for (int row = 0; row < (branch ? 2 : 1); row++)
        for (int i = 0; i < chainlength; i++)
        {
            input.start_m.push_back(static_cast<float>(chainbase + i + 50.0 * row));
            input.pinned.push_back(i == 0 && row == 0 ? 1 : 0);
            input.runoff_m_yr.push_back(row == 0 ? 1.0f : 0.0f);
        }
    return input;
}

// A dome without pits on a 40 x 24 patch of 500 m cells: 50 m per cell of distance
// to the edge plus up to 3 m of wobble, the ring pinned at its start, runoff 0.5 m/yr.
PatchMetric domemetric() { return squaremetric(40, 24, 500.0); }
AnalyticalReliefInput domeinput(const PatchMetric& metric)
{
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const int edge = std::min(std::min(x, y), std::min(metric.width - 1 - x, metric.height - 1 - y));
            input.start_m.push_back(static_cast<float>(50.0 * edge + 3.0 * wobble(x, y)));
            input.runoff_m_yr.push_back(0.5f);
        }
    return input;
}

// A bowl on a 30 x 20 patch of 1000 m cells: the interior is one depression below
// the pinned ring, runoff 1 m/yr.
PatchMetric bowlmetric() { return squaremetric(30, 20, 1000.0); }
AnalyticalReliefInput bowlinput(const PatchMetric& metric)
{
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const double dx = x - 14.5, dy = y - 9.5;
            input.start_m.push_back(static_cast<float>(2.0 * (dx * dx + dy * dy) + 20.0 * wobble(x, y)));
            input.runoff_m_yr.push_back(1.0f);
        }
    return input;
}

int cellsabovestart(const AnalyticalReliefInput& input, const AnalyticalReliefResult& result)
{
    int count = 0;
    for (std::size_t cell = 0; cell < input.start_m.size(); cell++)
        if (result.elevation_m[cell] > input.start_m[cell])
            count++;
    return count;
}
}

TEST_F(Landscape, AnalyticalReliefChainFollowsTheCharacteristicSolution)
{
    // Row 0 is the chain; row 1, without runoff, drains north into it: its cells carry
    // no discharge, never respond (Q^m = 0) and keep their start exactly. The floats
    // near 100 m resolve 7.6e-6 m, so elevations match within 4e-6 m.
    const PatchMetric metric = squaremetric(chainlength, 2, chainspacing);
    const AnalyticalReliefInput input = chaininput(true);
    AnalyticalReliefParameters parameters;
    parameters.erodibility = 750.0;
    parameters.areaexponent = 0.5;
    parameters.thresholdslope = 1.0;  // steeper than the start's 0.001

    // Independently: cell i (i >= 1) drains (N - i) cells of 1e6 m^2 at 1 m/yr, so
    // Q_i = (N - i) 1e-3 km^3/yr, the knickpoint speed on the link i -> i - 1 is
    // K Q_i^m and the response time at the cell's position is the sum of the link
    // times. A cell at position x takes z0 at the position whose response time is
    // tau(x) - t, piecewise linear in position, and z0(x) = 100 + 0.001 x.
    std::vector<double> tau(chainlength, 0.0);
    for (int i = 1; i < chainlength; i++)
        tau[static_cast<std::size_t>(i)] = tau[static_cast<std::size_t>(i - 1)]
                                           + chainspacing / (750.0 * std::sqrt((chainlength - i) * 1.0e-3));
    // The times in the test's own terms: tau = 18.856, 39.938, 64.281, 94.095, 136.259 Myr.
    EXPECT_NEAR(tau[1], 18.856181, 1e-6);
    EXPECT_NEAR(tau[5], 136.259199, 1e-6);

    for (double elapsed : { 30.0, 100.0 })
    {
        parameters.elapsed_myr = elapsed;
        const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
        EXPECT_TRUE(result.converged);
        EXPECT_EQ(result.changedreceivers, 0);
        EXPECT_EQ(result.iterations, parameters.maxiterations);
        EXPECT_EQ(result.elevation_m[0], 100.0f);
        EXPECT_EQ(result.direction[0], 0);
        EXPECT_EQ(result.discharge_m3_yr[0], 0.0);
        for (int i = 1; i < chainlength; i++)
        {
            const std::size_t cell = static_cast<std::size_t>(i);
            EXPECT_EQ(result.direction[cell], 7) << "cell " << i;  // west
            EXPECT_DOUBLE_EQ(result.discharge_m3_yr[cell], (chainlength - i) * 1.0e6) << "cell " << i;
            double expected = chainbase;
            const double target = tau[cell] - elapsed;
            if (target > 0.0)
            {
                int link = 1;
                while (tau[static_cast<std::size_t>(link)] < target)
                    link++;
                const double t0 = tau[static_cast<std::size_t>(link - 1)], t1 = tau[static_cast<std::size_t>(link)];
                const double position = chainspacing * ((link - 1) + (target - t0) / (t1 - t0));
                expected = chainbase + 0.001 * position;
            }
            EXPECT_NEAR(result.elevation_m[cell], expected, 4e-6) << "cell " << i << ", t = " << elapsed;
        }
        for (int i = 0; i < chainlength; i++)
        {
            const std::size_t cell = metric.index(i, 1);
            EXPECT_EQ(result.elevation_m[cell], input.start_m[cell]) << "branch cell " << i;
            EXPECT_EQ(result.discharge_m3_yr[cell], 0.0) << "branch cell " << i;
            EXPECT_EQ(result.direction[cell], 1) << "branch cell " << i;  // north
        }
    }
    // t = 30: cell 1 has responded (the outlet's 100 m, not its own 101 m), cell 2
    // lies 9.938 Myr up the first link.
    parameters.elapsed_myr = 30.0;
    const AnalyticalReliefResult early = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(early.elevation_m[1], 100.0f);
    EXPECT_NEAR(early.elevation_m[2], 100.0 + (39.938032 - 30.0) / 18.856181, 4e-6);
    // t = 100: cells 1 to 4 at the outlet's 100 m (not a receiver's start: cell 4's
    // receiver started at 103 m), cell 5 lies 36.259 Myr above cell 1.
    parameters.elapsed_myr = 100.0;
    const AnalyticalReliefResult late = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(late.elevation_m[4], 100.0f);
    EXPECT_NEAR(late.elevation_m[5], 101.0 + (136.259199 - 100.0 - 18.856181) / (39.938032 - 18.856181), 4e-6);
}

TEST_F(Landscape, AnalyticalReliefUpliftLiftsTheRespondedProfile)
{
    // Erosion plan R9: with a constant uplift U and t at least every response time, cell i
    // takes min(z0(i), z0(outlet) + U tau(i)), the steady-state profile of slope
    // U / (K Q^m) up to where it meets the start surface. U = 0.04 m/Myr: U tau = 0.754,
    // 1.598, 2.571, 3.764 and 5.450 m, so cells 1 to 4 lie on the profile and cell 5 is
    // clamped at its start (105 m). Within the reach only: at t = 30 cells 2 to 5 (tau > t)
    // advect the start as without uplift.
    const PatchMetric metric = squaremetric(chainlength, 1, chainspacing);
    const AnalyticalReliefInput plain = chaininput();
    AnalyticalReliefInput input = plain;
    constexpr double uplift = 0.04;
    input.uplift_m_myr.assign(input.start_m.size(), static_cast<float>(uplift));
    AnalyticalReliefParameters parameters;
    parameters.thresholdslope = 1.0;
    std::vector<double> tau(chainlength, 0.0);
    for (int i = 1; i < chainlength; i++)
        tau[static_cast<std::size_t>(i)] = tau[static_cast<std::size_t>(i - 1)]
                                           + chainspacing / (750.0 * std::sqrt((chainlength - i) * 1.0e-3));

    parameters.elapsed_myr = 200.0;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.elevation_m[0], 100.0f);
    for (int i = 1; i < chainlength; i++)
    {
        const std::size_t cell = static_cast<std::size_t>(i);
        EXPECT_EQ(result.direction[cell], 7) << "cell " << i;
        const double expected = std::min(chainbase + i, chainbase + uplift * tau[cell]);
        EXPECT_NEAR(result.elevation_m[cell], expected, 8e-6) << "cell " << i;
    }
    EXPECT_LT(result.elevation_m[4], 104.0f);
    EXPECT_EQ(result.elevation_m[5], 105.0f);

    parameters.elapsed_myr = 30.0;
    const AnalyticalReliefResult early = landscape::solveanalyticalrelief(metric, input, parameters);
    const AnalyticalReliefResult earlyplain = landscape::solveanalyticalrelief(metric, plain, parameters);
    EXPECT_NEAR(early.elevation_m[1], chainbase + uplift * tau[1], 8e-6);
    for (int i = 2; i < chainlength; i++)
        EXPECT_EQ(early.elevation_m[static_cast<std::size_t>(i)], earlyplain.elevation_m[static_cast<std::size_t>(i)])
            << "cell " << i;
}

TEST_F(Landscape, AnalyticalReliefZeroUpliftChangesNothing)
{
    // U = 0 everywhere gives the solution without uplift bit for bit, on the chain at
    // three times and on the dome (whose ring outlets differ in height) and the bowl.
    const PatchMetric chain = squaremetric(chainlength, 2, chainspacing);
    const AnalyticalReliefInput plain = chaininput(true);
    AnalyticalReliefInput zero = plain;
    zero.uplift_m_myr.assign(plain.start_m.size(), 0.0f);
    AnalyticalReliefParameters parameters;
    parameters.thresholdslope = 1.0;
    for (double elapsed : { 0.0, 30.0, 100.0 })
    {
        parameters.elapsed_myr = elapsed;
        EXPECT_EQ(landscape::solveanalyticalrelief(chain, zero, parameters).elevation_m,
                  landscape::solveanalyticalrelief(chain, plain, parameters).elevation_m)
            << "t = " << elapsed;
    }
    for (const PatchMetric& metric : { domemetric(), bowlmetric() })
    {
        const AnalyticalReliefInput base = metric.width == 40 ? domeinput(metric) : bowlinput(metric);
        AnalyticalReliefInput lifted = base;
        lifted.uplift_m_myr.assign(base.start_m.size(), 0.0f);
        AnalyticalReliefParameters defaults;
        defaults.elapsed_myr = 5.0;
        const AnalyticalReliefResult a = landscape::solveanalyticalrelief(metric, base, defaults);
        const AnalyticalReliefResult b = landscape::solveanalyticalrelief(metric, lifted, defaults);
        EXPECT_EQ(a.elevation_m, b.elevation_m);
        EXPECT_EQ(a.direction, b.direction);
    }
}

TEST_F(Landscape, AnalyticalReliefZeroTimeKeepsTheStartBelowTheThreshold)
{
    // A cliff on a chain: with t = 0 only the threshold slope (0.5) acts, from the
    // outlet upwards: 1500 -> 100 + 500, 1600 -> 600 + 500.
    const PatchMetric chain = squaremetric(4, 1, 1000.0);
    AnalyticalReliefInput cliff;
    cliff.start_m = { 0.0f, 100.0f, 1500.0f, 1600.0f };
    cliff.pinned = { 1, 0, 0, 0 };
    cliff.runoff_m_yr.assign(4, 1.0f);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 0.0;
    parameters.thresholdslope = 0.5;
    const AnalyticalReliefResult lowered = landscape::solveanalyticalrelief(chain, cliff, parameters);
    EXPECT_EQ(lowered.elevation_m, (std::vector<float>{ 0.0f, 100.0f, 600.0f, 1100.0f }));

    // A threshold steeper than every start slope: the start exactly, pits included.
    const PatchMetric metric = squaremetric(24, 16, 500.0);
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            input.start_m.push_back(static_cast<float>(40.0 * std::min(std::min(x, y), 6) + 100.0 * wobble(x, y)));
            input.runoff_m_yr.push_back(static_cast<float>(0.2 + wobble(y, x)));
        }
    parameters.thresholdslope = 10.0;  // every start slope is below (40 + 100) m / 500 m = 0.28
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(result.elevation_m, input.start_m);
    EXPECT_TRUE(result.converged);
}

TEST_F(Landscape, AnalyticalReliefKeepsPinnedCellsLowersAndRespectsTheThreshold)
{
    // The dome: slopes up to 0.106 against a threshold of 0.08, a partial response at
    // t = 5 Myr.
    const PatchMetric metric = domemetric();
    const AnalyticalReliefInput input = domeinput(metric);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 5.0;
    parameters.thresholdslope = 0.08;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
    int deep = 0, shallow = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (input.pinned[cell] != 0)
        {
            EXPECT_EQ(result.elevation_m[cell], input.start_m[cell]) << "cell " << cell;
            continue;
        }
        EXPECT_LE(result.elevation_m[cell], input.start_m[cell]) << "cell " << cell;
        const std::size_t receiver = metric.neighbour(cell, result.direction[cell]);
        const double slope
            = (static_cast<double>(result.elevation_m[cell]) - static_cast<double>(result.elevation_m[receiver]))
              / metric.distance_m(cell, result.direction[cell]);
        EXPECT_LE(slope, parameters.thresholdslope + 1e-9) << "cell " << cell;
        if (result.elevation_m[cell] < input.start_m[cell] - 25.0f)
            deep++;
        else
            shallow++;
    }
    // Cells lowered by more than half a step (the knickpoint has passed) and cells
    // lowered less (it has not): the test is not vacuous.
    EXPECT_GT(deep, 0);
    EXPECT_GT(shallow, 0);
}

TEST_F(Landscape, AnalyticalReliefBowlDrainsToTheRing)
{
    // The bowl's cells drain through the flood's directions towards the ring's lowest
    // cell, at 2 Myr with the default damping.
    const PatchMetric metric = bowlmetric();
    const AnalyticalReliefInput input = bowlinput(metric);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 2.0;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(result.iterations, parameters.maxiterations);
    EXPECT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
    // Every unpinned cell's runoff reaches the ring: the discharge into the ring is the total.
    double intoring = 0.0, total = 0.0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (input.pinned[cell] != 0)
        {
            EXPECT_EQ(result.discharge_m3_yr[cell], 0.0);
            continue;
        }
        total += 1.0e6;
        if (input.pinned[metric.neighbour(cell, result.direction[cell])] != 0)
            intoring += result.discharge_m3_yr[cell];
    }
    EXPECT_DOUBLE_EQ(intoring, total);
}

TEST_F(Landscape, AnalyticalReliefRunsEveryPassAndReportsTheLastTwo)
{
    // Every pass runs (the passes are a prefix of a longer run's), and converged says
    // whether the last two passes' receivers agree: the receivers of a run with n
    // passes against those of a run with n - 1. The bowl at 2 Myr with the default
    // damping agrees at pass 11 and not at passes 10 and 12 (it does not settle: 1 or
    // 2 cells keep flipping, agreeing again at 14, 23, 26, 35, ...).
    const PatchMetric metric = bowlmetric();
    const AnalyticalReliefInput input = bowlinput(metric);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 2.0;
    std::vector<uint8_t> previous;
    int agreeing = 0, disagreeing = 0;
    for (int passes = 9; passes <= 12; passes++)
    {
        parameters.maxiterations = passes;
        const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
        EXPECT_EQ(result.iterations, passes);
        if (!previous.empty())
        {
            EXPECT_EQ(result.converged, result.direction == previous) << passes << " passes";
            // The count is the cells whose receiver changed: 0 exactly when converged.
            EXPECT_EQ(result.converged, result.changedreceivers == 0) << passes << " passes";
            EXPECT_GE(result.changedreceivers, 0);
            (result.converged ? agreeing : disagreeing)++;
        }
        previous = result.direction;
    }
    // Which passes agree is a property of the defaults (damping, threshold), not of the
    // rule: reported, not asserted (E2c's re-review).
    std::cout << "[ info ] passes 10 to 12 agreeing with the previous: " << agreeing << ", disagreeing: " << disagreeing << "\n";
}

TEST_F(Landscape, AnalyticalReliefNoCatchmentCrossesADiagonalChannel)
{
    // Pinned cells on the diagonal x = y of a 12 x 12 patch at 0 m, the surface
    // 100 m x (x - y) on both sides. Cell (y + 1, y) would descend most steeply to
    // (y, y + 1) across the channel (200 m over 1414 m against 100 m over 1000 m to a
    // pinned cell), passing between the pinned corners (y, y) and (y + 1, y + 1); the
    // side below the line is one depression whose flood could leak the same way.
    const int size = 12;
    const PatchMetric metric = squaremetric(size, size, 1000.0);
    AnalyticalReliefInput input;
    for (int y = 0; y < size; y++)
        for (int x = 0; x < size; x++)
        {
            input.start_m.push_back(static_cast<float>(100.0 * (x - y)));
            input.pinned.push_back(x == y ? 1 : 0);
            input.runoff_m_yr.push_back(1.0f);
        }
    for (double elapsed : { 0.0, 2.0 })
    {
        AnalyticalReliefParameters parameters;
        parameters.elapsed_myr = elapsed;
        const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
        EXPECT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
        double above = 0.0, below = 0.0;  // discharge into the channel from each side
        for (int y = 0; y < size; y++)
            for (int x = 0; x < size; x++)
            {
                const std::size_t cell = metric.index(x, y);
                if (input.pinned[cell] != 0)
                    continue;
                const std::size_t receiver = metric.neighbour(cell, result.direction[cell]);
                const int rx = static_cast<int>(receiver % size), ry = static_cast<int>(receiver / size);
                if (input.pinned[receiver] != 0)
                    (x > y ? above : below) += result.discharge_m3_yr[cell];
                else
                    EXPECT_EQ(x > y, rx > ry) << "cell (" << x << ", " << y << "), t = " << elapsed;
            }
        // 66 cells on each side, 1e6 m^3/yr each.
        EXPECT_DOUBLE_EQ(above, 66.0e6) << "t = " << elapsed;
        EXPECT_DOUBLE_EQ(below, 66.0e6) << "t = " << elapsed;
    }
}

TEST_F(Landscape, AnalyticalReliefSubBasinsMatchAcrossOverlappingPatches)
{
    // A 60 x 40 patch of 500 m cells with pits everywhere, its ring pinned and a
    // channel network: full columns x = 8, 20, 38, 50, the row y = 14 between x = 20
    // and 38, and the 8-connected diagonal (21 + i, 15 + i), i = 0..16, below it.
    // Solved whole and as two sub-patches, columns 0..41 and 18..59, each with its own
    // ring pinned at the start, keys the whole patch's indices. S3.12: a flat-floored pit,
    // the square x = 25..31, y = 4..10 at 100 m (its surroundings 112 m and more), lies in
    // the basin above y = 14, so the flat resolution's cells are compared too.
    const int width = 60, height = 40;
    const PatchMetric whole = squaremetric(width, height, 500.0);
    AnalyticalReliefInput input;
    input.pinned = ring(whole);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const bool square = x >= 25 && x <= 31 && y >= 4 && y <= 10;
            input.start_m.push_back(square ? 100.0f
                                           : static_cast<float>(4.0 * (height - 1 - y) + 60.0 * wobble(x, y)));
            input.runoff_m_yr.push_back(static_cast<float>(0.3 + wobble(y + 3, x)));
            input.key.push_back(whole.index(x, y));
            if (x == 8 || x == 20 || x == 38 || x == 50 || (y == 14 && x >= 20 && x <= 38))
                input.pinned[whole.index(x, y)] = 1;
        }
    for (int i = 0; i <= 16; i++)
        input.pinned[whole.index(21 + i, 15 + i)] = 1;
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 3.0;
    const AnalyticalReliefResult full = landscape::solveanalyticalrelief(whole, input, parameters);

    const int subwidth = 42, leftcolumn = 0, rightcolumn = 18;
    const PatchMetric sub = squaremetric(subwidth, height, 500.0);
    const auto cut = [&](int column0) {
        AnalyticalReliefInput part;
        part.pinned = ring(sub);
        for (int y = 0; y < height; y++)
            for (int x = 0; x < subwidth; x++)
            {
                const std::size_t cell = whole.index(column0 + x, y);
                part.start_m.push_back(input.start_m[cell]);
                part.runoff_m_yr.push_back(input.runoff_m_yr[cell]);
                part.key.push_back(input.key[cell]);
                if (input.pinned[cell] != 0)
                    part.pinned[sub.index(x, y)] = 1;
            }
        return landscape::solveanalyticalrelief(sub, part, parameters);
    };
    const AnalyticalReliefResult left = cut(leftcolumn);
    const AnalyticalReliefResult right = cut(rightcolumn);

    // Sub-basins of the whole patch: unpinned cells connected by orthogonal steps and
    // by diagonal steps not between two pinned corners. Compared: the cells of every
    // sub-basin lying wholly inside both sub-patches' interiors (columns 19..40).
    std::vector<int> basin(whole.size(), -1);
    std::vector<uint8_t> inside;
    for (std::size_t seed = 0; seed < whole.size(); seed++)
    {
        if (input.pinned[seed] != 0 || basin[seed] >= 0)
            continue;
        const int id = static_cast<int>(inside.size());
        bool contained = true;
        std::vector<std::size_t> stack{ seed };
        basin[seed] = id;
        while (!stack.empty())
        {
            const std::size_t cell = stack.back();
            stack.pop_back();
            const int x = static_cast<int>(cell % width), y = static_cast<int>(cell / width);
            contained = contained && x >= rightcolumn + 1 && x <= leftcolumn + subwidth - 2;
            for (int direction = 1; direction <= 8; direction++)
            {
                const std::size_t next = whole.neighbour(cell, direction);
                if (next == cell || input.pinned[next] != 0 || basin[next] >= 0)
                    continue;
                const int nx = static_cast<int>(next % width), ny = static_cast<int>(next / width);
                if (nx != x && ny != y && input.pinned[whole.index(nx, y)] != 0 && input.pinned[whole.index(x, ny)] != 0)
                    continue;
                basin[next] = id;
                stack.push_back(next);
            }
        }
        inside.push_back(contained ? 1 : 0);
    }
    int compared = 0, differing = 0, squarecompared = 0;
    for (int y = 0; y < height; y++)
        for (int x = rightcolumn; x < leftcolumn + subwidth; x++)
        {
            const std::size_t cell = whole.index(x, y);
            const std::size_t l = sub.index(x - leftcolumn, y), r = sub.index(x - rightcolumn, y);
            if (input.pinned[cell] == 0 && inside[static_cast<std::size_t>(basin[cell])] != 0)
            {
                compared++;
                squarecompared += x >= 25 && x <= 31 && y >= 4 && y <= 10 ? 1 : 0;
                EXPECT_EQ(full.elevation_m[cell], left.elevation_m[l]) << "(" << x << ", " << y << ")";
                EXPECT_EQ(full.elevation_m[cell], right.elevation_m[r]) << "(" << x << ", " << y << ")";
                EXPECT_EQ(full.direction[cell], left.direction[l]) << "(" << x << ", " << y << ")";
                EXPECT_EQ(full.direction[cell], right.direction[r]) << "(" << x << ", " << y << ")";
                EXPECT_EQ(full.discharge_m3_yr[cell], left.discharge_m3_yr[l]) << "(" << x << ", " << y << ")";
                EXPECT_EQ(full.discharge_m3_yr[cell], right.discharge_m3_yr[r]) << "(" << x << ", " << y << ")";
            }
            else if (left.elevation_m[l] != right.elevation_m[r] || left.direction[l] != right.direction[r])
                differing++;
        }
    // The basin above y = 14 (17 x 13 cells) and the two below it (17 x 24 less the
    // 17 diagonal cells). Outside them, the two sub-patches differ on 102 of the 190
    // unpinned overlap cells (columns 18, 19 and 39..41), whose sub-basins each cut.
    EXPECT_EQ(compared, 221 + 391);
    EXPECT_EQ(squarecompared, 49);
    EXPECT_GT(differing, 0);
}

TEST_F(Landscape, AnalyticalReliefNoCellRises)
{
    // The cases that rose before the clamp (min of the start and the advected value):
    // the bowl's depression (120 cells at 0.5 and 2 Myr, about 110 at 20 Myr) and the
    // dome at 20 Myr with thresholds 0.104 and 1 (22 and 31 cells), damped and undamped.
    const PatchMetric bowl = bowlmetric();
    const AnalyticalReliefInput bowlstart = bowlinput(bowl);
    const PatchMetric dome = domemetric();
    const AnalyticalReliefInput domestart = domeinput(dome);
    for (double damping : { 1.0, 0.5 })
        for (double elapsed : { 0.5, 2.0, 20.0 })
        {
            AnalyticalReliefParameters parameters;
            parameters.damping = damping;
            parameters.elapsed_myr = elapsed;
            EXPECT_EQ(cellsabovestart(bowlstart, landscape::solveanalyticalrelief(bowl, bowlstart, parameters)), 0)
                << "bowl, t = " << elapsed << ", damping " << damping;
            for (double threshold : { 0.104, 1.0 })
            {
                parameters.thresholdslope = threshold;
                EXPECT_EQ(cellsabovestart(domestart, landscape::solveanalyticalrelief(dome, domestart, parameters)), 0)
                    << "dome, t = " << elapsed << ", threshold " << threshold << ", damping " << damping;
            }
        }
}

TEST_F(Landscape, AnalyticalReliefPitDrainsThroughItsRim)
{
    // Chain 0 (pinned), 18, 20, 5, 30: cell 2's steepest drop is into the pit (cell 3),
    // whose flood receiver is cell 2; the pit lies in the filled surface's depression,
    // so cell 2 drains west and the pit through its rim: every direction is west.
    const PatchMetric metric = squaremetric(5, 1, 1000.0);
    AnalyticalReliefInput input;
    input.start_m = { 0.0f, 18.0f, 20.0f, 5.0f, 30.0f };
    input.pinned = { 1, 0, 0, 0, 0 };
    input.runoff_m_yr.assign(5, 1.0f);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 0.0;
    parameters.thresholdslope = 1.0;
    parameters.maxiterations = 1;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(result.direction, (std::vector<uint8_t>{ 0, 7, 7, 7, 7 }));
    EXPECT_EQ(result.elevation_m, input.start_m);
    EXPECT_EQ(result.iterations, 1);
    EXPECT_FALSE(result.converged);
    EXPECT_EQ(result.changedreceivers, 0);  // one pass: nothing to compare
    EXPECT_EQ(result.discharge_m3_yr, (std::vector<double>{ 0.0, 4.0e6, 3.0e6, 2.0e6, 1.0e6 }));
}

TEST_F(Landscape, AnalyticalReliefTieTakesTheLowestCode)
{
    // The centre of a 3 x 3 patch above a flat pinned ring: four equal cardinal drops,
    // the lowest code (1, north) wins.
    const PatchMetric metric = squaremetric(3, 3, 1000.0);
    AnalyticalReliefInput input;
    input.start_m.assign(9, 0.0f);
    input.start_m[4] = 10.0f;
    input.pinned = ring(metric);
    input.runoff_m_yr.assign(9, 1.0f);
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, AnalyticalReliefParameters{});
    EXPECT_EQ(result.direction[4], 1);
}

TEST_F(Landscape, AnalyticalReliefIsRepeatable)
{
    // Pits everywhere: a rough surface, the ring and an interior column pinned.
    const PatchMetric metric = squaremetric(32, 16, 750.0);
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
    {
        input.pinned[metric.index(20, y)] = 1;
        for (int x = 0; x < metric.width; x++)
        {
            input.start_m.push_back(static_cast<float>(30.0 * std::min(x, 31 - x) + 200.0 * wobble(x, y)));
            input.runoff_m_yr.push_back(static_cast<float>(0.1 + wobble(y + 7, x)));
        }
    }
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 3.0;
    const AnalyticalReliefResult first = landscape::solveanalyticalrelief(metric, input, parameters);
    const AnalyticalReliefResult second = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(first.elevation_m, second.elevation_m);
    EXPECT_EQ(first.direction, second.direction);
    EXPECT_EQ(first.discharge_m3_yr, second.discharge_m3_yr);
    EXPECT_EQ(first.iterations, second.iterations);
    EXPECT_EQ(first.converged, second.converged);
    EXPECT_TRUE(everypathreachesapinnedcell(metric, input.pinned, first.direction));
}

TEST_F(Landscape, AnalyticalReliefKeysWithoutTiesChangeNothing)
{
    // The patch's 32 columns as a window crossing the date line on a fine grid of
    // kW = 64 columns: local column x is absolute column (x + 48) mod 64, so the keys
    // of columns 0..15 (absolute 48..63) exceed those of columns 16..31 (absolute
    // 0..15) in each row. With distinct elevations no tie arises: bit-identical.
    const PatchMetric metric = squaremetric(32, 16, 750.0);
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
    {
        input.pinned[metric.index(20, y)] = 1;
        for (int x = 0; x < metric.width; x++)
        {
            input.start_m.push_back(static_cast<float>(30.0 * std::min(x, 31 - x) + 200.0 * wobble(x, y)
                                                       + 0.01 * (y * metric.width + x)));
            input.runoff_m_yr.push_back(static_cast<float>(0.1 + wobble(y + 7, x)));
        }
    }
    std::vector<float> sorted = input.start_m;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_EQ(std::adjacent_find(sorted.begin(), sorted.end()), sorted.end());  // no two equal levels

    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 3.0;
    AnalyticalReliefInput local = input;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        local.key.push_back(cell);
    AnalyticalReliefInput wrapped = input;
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
            wrapped.key.push_back(static_cast<std::uint64_t>(y) * 64u + static_cast<std::uint64_t>((x + 48) % 64));
    const AnalyticalReliefResult byindex = landscape::solveanalyticalrelief(metric, local, parameters);
    const AnalyticalReliefResult bykey = landscape::solveanalyticalrelief(metric, wrapped, parameters);
    const AnalyticalReliefResult unkeyed = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(byindex.elevation_m, bykey.elevation_m);
    EXPECT_EQ(byindex.direction, bykey.direction);
    EXPECT_EQ(byindex.discharge_m3_yr, bykey.discharge_m3_yr);
    EXPECT_EQ(byindex.iterations, bykey.iterations);
    EXPECT_EQ(byindex.elevation_m, unkeyed.elevation_m);
    EXPECT_EQ(byindex.direction, unkeyed.direction);
}

TEST_F(Landscape, AnalyticalReliefFlatTiesFollowTheKey)
{
    // A flat chain of five cells at 0 m with both ends pinned: cells 1 to 3 are one flat whose low
    // edges are the pinned ends (S3.12). Cells 1 and 3 drain to their end (increment 2,000 m to 0);
    // cell 2 (4,000 m) ties between them and takes the neighbour with the smaller hashed key (the
    // splitmix64 finaliser of hashedkey). Keys equal to the local index: h(3) < h(1), so cell 2
    // drains east (3). As a window crossing the date line (absolute columns 60, 61, 62, 63, 0 of
    // kW = 64): h(61) < h(63), so it drains west (7).
    const PatchMetric metric = squaremetric(5, 1, 1000.0);
    AnalyticalReliefInput input;
    input.start_m.assign(5, 0.0f);
    input.pinned = { 1, 0, 0, 0, 1 };
    input.runoff_m_yr.assign(5, 1.0f);
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 1.0;
    input.key = { 0, 1, 2, 3, 4 };
    const AnalyticalReliefResult byindex = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(byindex.direction, (std::vector<uint8_t>{ 0, 7, 3, 3, 0 }));
    EXPECT_EQ(byindex.discharge_m3_yr, (std::vector<double>{ 0.0, 1.0e6, 1.0e6, 2.0e6, 0.0 }));
    input.key = { 60, 61, 62, 63, 0 };
    const AnalyticalReliefResult bykey = landscape::solveanalyticalrelief(metric, input, parameters);
    EXPECT_EQ(bykey.direction, (std::vector<uint8_t>{ 0, 7, 7, 3, 0 }));
    EXPECT_EQ(bykey.discharge_m3_yr, (std::vector<double>{ 0.0, 2.0e6, 1.0e6, 1.0e6, 0.0 }));
    EXPECT_EQ(bykey.elevation_m, input.start_m);
}

TEST_F(Landscape, AnalyticalReliefPitDrainsToItsSpillAlongShortChains)
{
    // S3.12. A plane falling west (20 m + 0.5 m per column, 100 m cells, the ring pinned) with a
    // round pit of radius 9 cells around (30, 20), 5 m + 0.05 m per squared cell of radius: the
    // flood fills it to the spill level 30 m, whose cells (20, 19..21) drain west; the pit's only
    // cell next to them is its west tip (21, 20). Every pit cell's chain converges on the tip and
    // leaves the pit to the spill (20, 20), along a length at most 1.5 x the straight distance to it
    // (the resolution's gradient towards the spill and away from the rim).
    const PatchMetric metric = squaremetric(48, 40, 100.0);
    const int cx = 30, cy = 20, radius = 9;
    const auto inpit = [&](int x, int y) { return (x - cx) * (x - cx) + (y - cy) * (y - cy) <= radius * radius; };
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const int r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
            input.start_m.push_back(static_cast<float>(inpit(x, y) ? 5.0 + 0.05 * r2 : 20.0 + 0.5 * x));
            input.runoff_m_yr.push_back(1.0f);
        }
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 0.0;
    parameters.maxiterations = 1;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    ASSERT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
    int cells = 0;
    double worst = 0.0;
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            if (!inpit(x, y))
                continue;
            cells++;
            std::size_t current = metric.index(x, y);
            double length = 0.0;
            int cxl = x, cyl = y;
            while (inpit(cxl, cyl))
            {
                length += metric.distance_m(current, result.direction[current]);
                current = metric.neighbour(current, result.direction[current]);
                cxl = static_cast<int>(current % static_cast<std::size_t>(metric.width));
                cyl = static_cast<int>(current / static_cast<std::size_t>(metric.width));
            }
            EXPECT_EQ(cxl, cx - radius - 1) << "(" << x << ", " << y << ") leaves the pit at (" << cxl << ", " << cyl << ")";
            EXPECT_EQ(cyl, cy) << "(" << x << ", " << y << ") leaves the pit at (" << cxl << ", " << cyl << ")";
            const double straight = 100.0 * std::hypot(x - (cx - radius - 1), y - cy);
            EXPECT_LE(length, 1.5 * straight) << "(" << x << ", " << y << ")";
            worst = std::max(worst, length / straight);
        }
    EXPECT_GT(cells, 200);
    std::printf("[relief] pit: %d cells, longest chain %.3f x the straight distance to the spill\n", cells, worst);
}

TEST_F(Landscape, AnalyticalReliefFlatBetweenSlopesDrainsInsideItself)
{
    // S3.12. A valley floor at 10 m (rows 5..9, columns 1..23, 100 m cells) between hillsides
    // rising 10 m per row away from it (20 m and more), draining to column 0 (pinned; 5 m on the
    // floor). Columns 2..23 of the floor are one flat; its low edges are the floor cells of
    // column 1, which drain west. Every flat cell's chain stays on the floor, never steps up and
    // reaches column 0; away from the low edges the gradient away from the hillsides (Barnes et
    // al. 2014) turns rows 5 and 6 south-west, rows 8 and 9 north-west, and row 7 runs west.
    const PatchMetric metric = squaremetric(24, 15, 100.0);
    const auto onfloor = [](int x, int y) { return x >= 1 && y >= 5 && y <= 9; };
    AnalyticalReliefInput input;
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            double z = 20.0 + 10.0 * (std::abs(y - 7) - 3) + 0.1 * x;
            if (y >= 5 && y <= 9)
                z = x == 0 ? 5.0 : 10.0;
            input.start_m.push_back(static_cast<float>(z));
            input.pinned.push_back(x == 0 ? 1 : 0);
            input.runoff_m_yr.push_back(1.0f);
        }
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 0.0;
    parameters.maxiterations = 1;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    ASSERT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
    for (int y = 5; y <= 9; y++)
        for (int x = 2; x < metric.width; x++)
        {
            std::size_t current = metric.index(x, y);
            while (input.pinned[current] == 0)
            {
                const std::size_t next = metric.neighbour(current, result.direction[current]);
                const int nx = static_cast<int>(next % static_cast<std::size_t>(metric.width));
                const int ny = static_cast<int>(next / static_cast<std::size_t>(metric.width));
                EXPECT_TRUE(onfloor(nx, ny) || nx == 0) << "(" << x << ", " << y << ") leaves the floor at (" << nx << ", " << ny << ")";
                EXPECT_LE(input.start_m[next], input.start_m[current]) << "(" << x << ", " << y << ") steps up";
                if (!(onfloor(nx, ny) || nx == 0))
                    break;
                current = next;
            }
            if (x >= 4)
            {
                const int expected = y < 7 ? 6 : (y > 7 ? 8 : 7);
                EXPECT_EQ(result.direction[metric.index(x, y)], expected) << "(" << x << ", " << y << ")";
            }
        }
}

TEST_F(Landscape, AnalyticalReliefSquareFlatTiesFollowTheHashedKey)
{
    // S3.12. A 15 x 15 patch of 100 m cells at 0 m, its ring pinned: the interior is one flat whose
    // low edges are the ring. The increment is twice the distance to the ring, so each cell drains
    // to a nearest side by a cardinal step; on the square's diagonals two sides (four at the centre)
    // tie, and the hashed key decides: not the lowest code everywhere (a raster order), not one code
    // everywhere. The same cells (same absolute keys) at an offset inside a larger patch whose other
    // cells are pinned give the same result: the resolution reads only the flat and its edges.
    constexpr int size = 15, offsetx = 4, offsety = 3, bigwidth = 21, bigheight = 19;
    const auto absolutekey = [](int x, int y) {
        return static_cast<std::uint64_t>(y + 200) * 1000u + static_cast<std::uint64_t>(x + 100);
    };
    const PatchMetric metric = squaremetric(size, size, 100.0);
    AnalyticalReliefInput input;
    input.pinned = ring(metric);
    input.start_m.assign(metric.size(), 0.0f);
    input.runoff_m_yr.assign(metric.size(), 1.0f);
    for (int y = 0; y < size; y++)
        for (int x = 0; x < size; x++)
            input.key.push_back(absolutekey(x, y));
    const PatchMetric big = squaremetric(bigwidth, bigheight, 100.0);
    AnalyticalReliefInput shifted;
    shifted.start_m.assign(big.size(), 0.0f);
    shifted.runoff_m_yr.assign(big.size(), 1.0f);
    for (int y = 0; y < bigheight; y++)
        for (int x = 0; x < bigwidth; x++)
        {
            const int sx = x - offsetx, sy = y - offsety;
            const bool interior = sx >= 1 && sx <= size - 2 && sy >= 1 && sy <= size - 2;
            shifted.pinned.push_back(interior ? 0 : 1);
            shifted.key.push_back(absolutekey(sx + 50, sy + 50) + (interior ? 0u : 1000000u));
        }
    for (int y = 0; y < size; y++)
        for (int x = 0; x < size; x++)
            shifted.key[big.index(x + offsetx, y + offsety)] = absolutekey(x, y);
    AnalyticalReliefParameters parameters;
    parameters.maxiterations = 1;
    const AnalyticalReliefResult result = landscape::solveanalyticalrelief(metric, input, parameters);
    const AnalyticalReliefResult moved = landscape::solveanalyticalrelief(big, shifted, parameters);
    ASSERT_TRUE(everypathreachesapinnedcell(metric, input.pinned, result.direction));
    int lowest = 0, other = 0, tied = 0;
    std::vector<int> codes(9, 0);
    for (int y = 1; y < size - 1; y++)
        for (int x = 1; x < size - 1; x++)
        {
            const std::size_t cell = metric.index(x, y), there = big.index(x + offsetx, y + offsety);
            EXPECT_EQ(result.direction[cell], moved.direction[there]) << "(" << x << ", " << y << ")";
            EXPECT_EQ(result.discharge_m3_yr[cell], moved.discharge_m3_yr[there]) << "(" << x << ", " << y << ")";
            const int nearest = std::min(std::min(x, y), std::min(size - 1 - x, size - 1 - y));
            std::vector<int> options;  // cardinal codes towards a nearest side, increasing
            if (y == nearest)
                options.push_back(1);
            if (size - 1 - x == nearest)
                options.push_back(3);
            if (size - 1 - y == nearest)
                options.push_back(5);
            if (x == nearest)
                options.push_back(7);
            const int code = result.direction[cell];
            codes[static_cast<std::size_t>(code)]++;
            EXPECT_NE(std::find(options.begin(), options.end(), code), options.end()) << "(" << x << ", " << y << ")";
            if (options.size() < 2)
                continue;
            tied++;
            (code == options.front() ? lowest : other)++;
        }
    EXPECT_EQ(tied, 25);
    EXPECT_GT(lowest, 0);
    EXPECT_GT(other, 0);
    for (int code = 1; code <= 8; code++)
        EXPECT_LT(codes[static_cast<std::size_t>(code)], (size - 2) * (size - 2)) << "code " << code;
    std::printf("[relief] square flat: %d tied cells, %d take the lowest code, %d another\n", tied, lowest, other);
}

TEST_F(Landscape, AnalyticalReliefRejectsBadInput)
{
    const PatchMetric metric = squaremetric(chainlength, 1, chainspacing);
    const AnalyticalReliefInput good = chaininput();
    const AnalyticalReliefParameters defaults;
    EXPECT_NO_THROW(landscape::solveanalyticalrelief(metric, good, defaults));

    AnalyticalReliefInput input = good;
    input.start_m.pop_back();
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input = good;
    input.pinned.push_back(0);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input = good;
    input.runoff_m_yr.clear();
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input = good;
    input.pinned.assign(input.pinned.size(), 0);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input = good;
    input.key.assign(good.start_m.size() - 1, 0);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input = good;
    input.key = { 5, 9, 2, 7, 9, 0 };  // 9 twice
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input.key = { 5, 9, 2, 7, 8, 0 };
    EXPECT_NO_THROW(landscape::solveanalyticalrelief(metric, input, defaults));
    input = good;
    input.uplift_m_myr.assign(good.start_m.size() - 1, 1.0f);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input.uplift_m_myr.assign(good.start_m.size() + 1, 1.0f);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, defaults), std::invalid_argument);
    input.uplift_m_myr.assign(good.start_m.size(), 1.0f);
    EXPECT_NO_THROW(landscape::solveanalyticalrelief(metric, input, defaults));

    AnalyticalReliefParameters parameters = defaults;
    parameters.erodibility = 0.0;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters = defaults;
    parameters.areaexponent = -0.1;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters = defaults;
    parameters.elapsed_myr = -1.0;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters.elapsed_myr = std::numeric_limits<double>::infinity();
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters.elapsed_myr = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters = defaults;
    parameters.thresholdslope = 0.0;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters = defaults;
    parameters.maxiterations = 0;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters = defaults;
    parameters.damping = 0.0;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters.damping = 1.5;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters.damping = 1.0;
    EXPECT_NO_THROW(landscape::solveanalyticalrelief(metric, good, parameters));

    // S3.9: a negative catchment reach, or a positive one without the global cells.
    parameters = defaults;
    parameters.catchmentreach = -1;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    parameters.catchmentreach = 1;
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, good, parameters), std::invalid_argument);
    input = good;
    input.globalcolumn.assign(good.start_m.size(), 0);
    input.globalrow.assign(good.start_m.size(), 0);
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, parameters), std::invalid_argument);  // globalwidth 0
    input.globalwidth = 4;
    EXPECT_NO_THROW(landscape::solveanalyticalrelief(metric, input, parameters));
    input.globalrow.pop_back();
    EXPECT_THROW(landscape::solveanalyticalrelief(metric, input, parameters), std::invalid_argument);
}

namespace
{
// S3.9: a patch solved whole and as its sub-patch of the columns from `subcolumn` on (its own ring
// pinned at the start, keys and global cells the whole patch's), compared on the unpinned cells at
// least `blocks` global cells east of the cut (columns from subcolumn + blocks x block to the last
// but one).
struct CutComparison
{
    int compared = 0, elevations = 0, directions = 0, discharges = 0;
    int unreachedwhole = 0, unreachedcut = 0;
};

CutComparison comparecut(const PatchMetric& whole, const AnalyticalReliefInput& input, int subcolumn, int block,
                         const AnalyticalReliefParameters& parameters, const std::vector<int>& distances,
                         std::vector<CutComparison>& results, AnalyticalReliefResult* wholeresult = nullptr)
{
    const int width = whole.width, height = whole.height, subwidth = width - subcolumn;
    const PatchMetric sub = squaremetric(subwidth, height, whole.dy_m);
    AnalyticalReliefInput part;
    part.pinned = ring(sub);
    part.globalwidth = input.globalwidth;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < subwidth; x++)
        {
            const std::size_t cell = whole.index(subcolumn + x, y);
            part.start_m.push_back(input.start_m[cell]);
            part.runoff_m_yr.push_back(input.runoff_m_yr[cell]);
            part.key.push_back(input.key[cell]);
            if (!input.globalcolumn.empty())
            {
                part.globalcolumn.push_back(input.globalcolumn[cell]);
                part.globalrow.push_back(input.globalrow[cell]);
            }
            if (input.pinned[cell] != 0)
                part.pinned[sub.index(x, y)] = 1;
        }
    const AnalyticalReliefResult full = landscape::solveanalyticalrelief(whole, input, parameters);
    const AnalyticalReliefResult cut = landscape::solveanalyticalrelief(sub, part, parameters);
    results.clear();
    for (int blocks : distances)
    {
        CutComparison counts;
        counts.unreachedwhole = full.unreachedcells;
        counts.unreachedcut = cut.unreachedcells;
        for (int y = 1; y < height - 1; y++)
            for (int x = subcolumn + blocks * block; x < width - 1; x++)
            {
                const std::size_t cell = whole.index(x, y), other = sub.index(x - subcolumn, y);
                if (input.pinned[cell] != 0)
                    continue;
                counts.compared++;
                counts.elevations += full.elevation_m[cell] != cut.elevation_m[other] ? 1 : 0;
                counts.directions += full.direction[cell] != cut.direction[other] ? 1 : 0;
                counts.discharges += full.discharge_m3_yr[cell] != cut.discharge_m3_yr[other] ? 1 : 0;
            }
        results.push_back(counts);
    }
    if (wholeresult != nullptr)
        *wholeresult = full;
    return results.back();
}
}

TEST_F(Landscape, AnalyticalReliefBoundedCatchmentsCutAWindowIndependentTrough)
{
    // S3.9 (multiscale plan §3, "Seamless solve: a bounded catchment", verify (a)). A 96 x 40 patch
    // of 1,000 m cells in blocks of 8 x 8 cells (the global cells): a trough along row 20 sloping
    // down westwards across eleven blocks to one pinned cell at (4, 20), walled to its west; pinned
    // lines at 600 m along rows 8 and 31 above its rims; slopes of 650 m and more beyond them
    // falling towards the lines; the ring pinned. Solved whole and as the sub-patch of columns
    // 32 .. 95, whose ring (column 32) cuts the trough. With catchment reach 0 the trough drains to the
    // far cell in the whole patch and to the cut in the sub-patch, so the cells two blocks inside the cut
    // differ; with reach 1 the cells at least 4 blocks inside (the bound, 4r) are identical, and
    // every outlet lies within one block of its cell's block.
    constexpr int block = 8, width = 96, height = 40, subcolumn = 32;
    const PatchMetric whole = squaremetric(width, height, 1000.0);
    AnalyticalReliefInput input;
    input.pinned = ring(whole);
    input.globalwidth = 1000;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            double z = 0.0;
            if (y == 8 || y == 31)
                z = 600.0;
            else if (y < 8)
                z = 650.0 + 5.0 * (8 - y) + 3.0 * wobble(x, y);
            else if (y > 31)
                z = 650.0 + 5.0 * (y - 31) + 3.0 * wobble(x, y);
            else if (x < 4)
                z = 700.0 + 3.0 * wobble(x, y);
            else
                z = 10.0 + 2.0 * (x - 4) + 20.0 * std::abs(y - 20) + 3.0 * wobble(x, y);
            input.start_m.push_back(static_cast<float>(z));
            input.runoff_m_yr.push_back(1.0f);
            input.key.push_back(whole.index(x, y));
            input.globalcolumn.push_back(x / block);
            input.globalrow.push_back(y / block);
            if (y == 8 || y == 31 || (x == 4 && y == 20))
                input.pinned[whole.index(x, y)] = 1;
        }

    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 200.0;
    const std::vector<int> distances{ 2, 3, 4 };
    std::vector<CutComparison> results;
    for (int reach : { 0, 1 })
    {
        parameters.catchmentreach = reach;
        AnalyticalReliefResult full;
        comparecut(whole, input, subcolumn, block, parameters, distances, results, &full);
        for (std::size_t d = 0; d < distances.size(); d++)
            std::printf("[relief] trough, catchment reach %d, %d blocks inside the cut: %d cells compared, %d elevations, "
                        "%d directions and %d discharges differ; unreached %d whole, %d cut\n",
                        reach, distances[d], results[d].compared, results[d].elevations, results[d].directions,
                        results[d].discharges, results[d].unreachedwhole, results[d].unreachedcut);
        if (reach == 0)
        {
            EXPECT_GT(results[0].elevations, 0) << "the unbounded solve: the cut changes the trough two blocks inside it";
            EXPECT_TRUE(full.unreached.empty());
            continue;
        }
        const CutComparison& bound = results.back();
        EXPECT_GT(bound.compared, 0);
        EXPECT_EQ(bound.elevations, 0) << "catchment reach " << reach;
        EXPECT_EQ(bound.directions, 0) << "catchment reach " << reach;
        EXPECT_EQ(bound.discharges, 0) << "catchment reach " << reach;
        EXPECT_EQ(full.unreached.size(), whole.size());
        for (std::size_t cell = 0; cell < whole.size(); cell++)
        {
            const std::size_t outlet = static_cast<std::size_t>(full.outlet[cell]);
            EXPECT_LE(std::abs(input.globalcolumn[cell] - input.globalcolumn[outlet]), reach);
            EXPECT_LE(std::abs(input.globalrow[cell] - input.globalrow[outlet]), reach);
        }
    }
}

TEST_F(Landscape, AnalyticalReliefBoundedCatchmentsDominoCorridor)
{
    // S3.9's review: the adversarial case of a labelling race. A 96 x 24 patch of 1,000 m cells in
    // blocks of 8 x 8: walls of 1,000 m, the ring pinned, a flat corridor at 100 m along row 12 (columns
    // 1 .. 94) and one pinned cell at 200 m above it at the west end of every block ((8b + 1, 11)). A race
    // seeded at all pinned cells in the order of the row-major key lets the west-most seed fill the
    // corridor up to its reach, the next block's seed find its entry taken and the one after start
    // again, so the labels alternate block by block from the patch's west end and a cut flips them
    // all, at any distance. With one flood per global cell and the hashed tie the labels are local:
    // solved whole and as the sub-patch of the columns 16 .. 95 (its ring cutting the corridor), the
    // cells at least 4 blocks east of the cut (the bound, 4r) are identical.
    constexpr int block = 8, width = 96, height = 24, subcolumn = 16;
    const PatchMetric whole = squaremetric(width, height, 1000.0);
    AnalyticalReliefInput input;
    input.pinned = ring(whole);
    input.globalwidth = 1000;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const bool corridor = y == 12 && x >= 1 && x <= width - 2;
            const bool seed = y == 11 && x % block == 1;
            input.start_m.push_back(corridor ? 100.0f : (seed ? 200.0f : 1000.0f));
            input.runoff_m_yr.push_back(1.0f);
            input.key.push_back(whole.index(x, y));
            input.globalcolumn.push_back(x / block);
            input.globalrow.push_back(y / block);
            if (seed)
                input.pinned[whole.index(x, y)] = 1;
        }
    AnalyticalReliefParameters parameters;
    parameters.elapsed_myr = 5.0;
    parameters.catchmentreach = 1;
    const std::vector<int> distances{ 1, 2, 3, 4 };
    std::vector<CutComparison> results;
    comparecut(whole, input, subcolumn, block, parameters, distances, results);
    for (std::size_t d = 0; d < distances.size(); d++)
        std::printf("[relief] domino corridor, %d blocks east of the cut: %d cells compared, %d elevations, %d directions "
                    "and %d discharges differ; unreached %d whole, %d cut\n",
                    distances[d], results[d].compared, results[d].elevations, results[d].directions, results[d].discharges,
                    results[d].unreachedwhole, results[d].unreachedcut);
    const CutComparison& bound = results.back();
    EXPECT_GT(bound.compared, 0);
    EXPECT_EQ(bound.elevations, 0);
    EXPECT_EQ(bound.directions, 0);
    EXPECT_EQ(bound.discharges, 0);
}
