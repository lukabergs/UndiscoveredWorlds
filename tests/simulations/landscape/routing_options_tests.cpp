// Landscape.RoutingOptions*: the randomised receivers and the flat resolution of
// buildflowgraph on synthetic grids (docs/architecture/erosion-hydrology-plan.md,
// Phase C, step C3). Landscape.PathRouting*: the path-corrected receivers (D8-LTD;
// master plan 5.2, terrain-roadmap.md "Item 3, stage 1 design").

#include "landscape_fixture.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;

struct Terrain
{
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
};

// A plane descending at azimuth `azimuth` (radians anticlockwise from east) with
// slope `slope`, the east term scaled by the dx of row `referencerow`, so that in
// row y the descent is s cos(azimuth) dx_ref / dx_y east and s sin(azimuth) north;
// the last column is sea.
Terrain azimuthplane(const landscape::SphereMetric& metric, double azimuth, double slope, int referencerow)
{
    Terrain terrain;
    terrain.elevation.resize(metric.size());
    terrain.sea.resize(metric.size());
    const double dxref = metric.rowdx_m[static_cast<std::size_t>(referencerow)];
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            const bool sea = x == metric.width - 1;
            terrain.sea[cell] = sea ? 1 : 0;
            const double height = 10000.0 - slope * std::cos(azimuth) * x * dxref + slope * std::sin(azimuth) * y * metric.dy_m;
            terrain.elevation[cell] = sea ? -100.0f : static_cast<float>(height);
        }
    return terrain;
}

// Column 0 is sea at -100 m; land elsewhere rises 100 m per column (as in landscape_tests.cpp).
Terrain tiltedplane(const landscape::SphereMetric& metric)
{
    Terrain terrain;
    terrain.elevation.resize(metric.size());
    terrain.sea.resize(metric.size());
    for (int y = 0; y < metric.height; y++)
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            terrain.sea[cell] = x == 0 ? 1 : 0;
            terrain.elevation[cell] = x == 0 ? -100.0f : 100.0f * static_cast<float>(x);
        }
    return terrain;
}

bool inbowl(int x, int y) { return x >= 6 && x <= 8 && y >= 3 && y <= 5; }

// The plane with a 3 x 3 bowl at columns 6..8, rows 3..5, lowered by 500 m.
Terrain bowlplane(const landscape::SphereMetric& metric)
{
    Terrain terrain = tiltedplane(metric);
    for (int y = 3; y <= 5; y++)
        for (int x = 6; x <= 8; x++)
            terrain.elevation[metric.index(x, y)] -= 500.0f;
    return terrain;
}

// The plane with a shelf: columns 4..9 all at 400 m, so columns 5..9 have no lower neighbour.
Terrain shelfplane(const landscape::SphereMetric& metric)
{
    Terrain terrain = tiltedplane(metric);
    for (int y = 0; y < metric.height; y++)
        for (int x = 4; x <= 9; x++)
            terrain.elevation[metric.index(x, y)] = 400.0f;
    return terrain;
}

landscape::FlowGraph route(const landscape::SphereMetric& metric, const Terrain& terrain,
                           const landscape::RoutingOptions& options)
{
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    return landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions, options);
}

landscape::RoutingOptions options(bool randomised, std::uint64_t seed, bool resolveflats)
{
    landscape::RoutingOptions result;
    result.randomised = randomised;
    result.seed = seed;
    result.resolveflats = resolveflats;
    return result;
}

// The cells on the receiver chain from a cell to its end (the cell that is its own
// receiver); empty when the chain is longer than the grid (a cycle).
std::vector<std::size_t> chain(const landscape::FlowGraph& graph, std::size_t cell)
{
    std::vector<std::size_t> cells(1, cell);
    for (std::size_t step = 0; step <= graph.size(); step++)
    {
        const std::size_t next = static_cast<std::size_t>(graph.receiver[cells.back()]);
        if (next == cells.back())
            return cells;
        cells.push_back(next);
    }
    return {};
}

// Every land cell appears once in `order`, after its receiver.
void expectorder(const landscape::FlowGraph& graph)
{
    std::vector<int> position(graph.size(), -1);
    for (std::size_t k = 0; k < graph.order.size(); k++)
    {
        const std::size_t cell = static_cast<std::size_t>(graph.order[k]);
        ASSERT_LT(cell, graph.size());
        EXPECT_EQ(static_cast<int>(graph.land[cell]), 1) << "cell " << cell;
        EXPECT_EQ(position[cell], -1) << "cell " << cell << " listed twice";
        position[cell] = static_cast<int>(k);
    }
    for (std::size_t cell = 0; cell < graph.size(); cell++)
    {
        if (graph.land[cell] == 0)
            continue;
        ASSERT_GE(position[cell], 0) << "land cell " << cell << " missing from order";
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver != cell)
            EXPECT_LT(position[receiver], position[cell]) << "cell " << cell;
    }
}

// Every land cell has a direction, its chain ends at an outlet into the sea, and
// each receiver is strictly lower on the filled surface or (flat cells) at its level.
void expectdrained(const landscape::SphereMetric& metric, const Terrain& terrain, const landscape::FlowGraph& graph,
                   bool strictlylower)
{
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    const std::vector<float>& filled = depressions.filledelevation;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (terrain.sea[cell] != 0)
            continue;
        EXPECT_GE(static_cast<int>(graph.direction[cell]), 1) << "cell " << cell;
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver == cell)
        {
            ASSERT_GE(graph.outletsea[cell], 0) << "cell " << cell;
            EXPECT_EQ(static_cast<int>(terrain.sea[static_cast<std::size_t>(graph.outletsea[cell])]), 1) << "cell " << cell;
            EXPECT_EQ(metric.neighbour(cell, graph.direction[cell]), static_cast<std::size_t>(graph.outletsea[cell]));
        }
        else
        {
            EXPECT_EQ(metric.neighbour(cell, graph.direction[cell]), receiver) << "cell " << cell;
            if (strictlylower)
                EXPECT_LT(filled[receiver], filled[cell]) << "cell " << cell;
            else
                EXPECT_LE(filled[receiver], filled[cell]) << "cell " << cell;
        }
        const std::vector<std::size_t> cells = chain(graph, cell);
        ASSERT_FALSE(cells.empty()) << "cell " << cell << " is on a cycle";
        EXPECT_GE(graph.outletsea[cells.back()], 0) << "cell " << cell;
    }
    expectorder(graph);
}

// A neighbour's offset on the tangent plane in metres (east, north): the row's dx
// east-west, dy north-south, the mean dx of the two rows for a diagonal.
struct Displacement
{
    double east = 0.0;
    double north = 0.0;
};

Displacement neighbouroffset(const landscape::SphereMetric& metric, std::size_t cell, int direction)
{
    const int y = metric.row(cell);
    const int ox = landscape::directionoffsetx(direction), oy = landscape::directionoffsety(direction);
    double dx = metric.rowdx_m[static_cast<std::size_t>(y)];
    if (ox != 0 && oy != 0)
        dx = 0.5 * (dx + metric.rowdx_m[static_cast<std::size_t>(y + oy)]);
    return { ox * dx, -oy * metric.dy_m };
}

// Per row of a band, the deviation (degrees) of the summed receiver displacement
// from the gradient azimuth, for deterministic D8, for the randomised receivers
// (four seeds summed), and for the expectation of the angle rule (the second edge
// drawn with probability angle(p1, g) / angle(p1, p2)), on a 512 x 256 plane
// descending at 30 degrees north of east at the band's middle row.
struct BandDeviations
{
    std::vector<double> d8, randomised, anglerule;
};

BandDeviations banddeviations(int firstrow, int lastrow)
{
    const landscape::SphereMetric metric(512, 256, radius);
    const double azimuth = pi / 6.0;
    const int referencerow = (firstrow + lastrow) / 2;
    const Terrain terrain = azimuthplane(metric, azimuth, 1.0e-4, referencerow);
    const landscape::FlowGraph d8 = route(metric, terrain, {});
    std::vector<landscape::FlowGraph> randomised;
    // Seeds that differ above the cell index's bits: the draw hashes seed ^ cell, so
    // seeds 1..4 would only permute the draws among neighbouring cells.
    for (std::uint64_t seed = 1; seed <= 4; seed++)
        randomised.push_back(route(metric, terrain, options(true, seed << 32, false)));

    const double degrees = 180.0 / pi;
    const double dxref = metric.rowdx_m[static_cast<std::size_t>(referencerow)];
    BandDeviations result;
    for (int y = firstrow; y <= lastrow; y++)
    {
        const double ge = std::cos(azimuth) * dxref / metric.rowdx_m[static_cast<std::size_t>(y)];
        const double gn = std::sin(azimuth);
        const double expected = std::atan2(gn, ge);
        Displacement d8sum, randomisedsum, anglesum;
        for (int x = 2; x <= metric.width - 4; x++)
        {
            const std::size_t cell = metric.index(x, y);
            const Displacement chosen = neighbouroffset(metric, cell, d8.direction[cell]);
            d8sum.east += chosen.east;
            d8sum.north += chosen.north;
            for (const landscape::FlowGraph& graph : randomised)
            {
                const Displacement drawn = neighbouroffset(metric, cell, graph.direction[cell]);
                randomisedsum.east += drawn.east;
                randomisedsum.north += drawn.north;
            }
            // The gradient lies in the facet of the east (p1) and north-east (p2) neighbours.
            const Displacement p1 = neighbouroffset(metric, cell, 3), p2 = neighbouroffset(metric, cell, 2);
            const double span = std::atan2(p2.north, p2.east);
            EXPECT_LT(expected, span) << "row " << y;
            const double share = expected / span;
            anglesum.east += (1.0 - share) * p1.east + share * p2.east;
            anglesum.north += (1.0 - share) * p1.north + share * p2.north;
        }
        result.d8.push_back((std::atan2(d8sum.north, d8sum.east) - expected) * degrees);
        result.randomised.push_back((std::atan2(randomisedsum.north, randomisedsum.east) - expected) * degrees);
        result.anglerule.push_back((std::atan2(anglesum.north, anglesum.east) - expected) * degrees);
    }
    return result;
}

double mean(const std::vector<double>& values)
{
    double sum = 0.0;
    for (double value : values)
        sum += value;
    return values.empty() ? 0.0 : sum / static_cast<double>(values.size());
}

double meanmagnitude(const std::vector<double>& values)
{
    double sum = 0.0;
    for (double value : values)
        sum += std::fabs(value);
    return values.empty() ? 0.0 : sum / static_cast<double>(values.size());
}

double largestmagnitude(const std::vector<double>& values)
{
    double largest = 0.0;
    for (double value : values)
        largest = std::max(largest, std::fabs(value));
    return largest;
}

void expectsamegraph(const landscape::FlowGraph& a, const landscape::FlowGraph& b)
{
    EXPECT_EQ(a.receiver, b.receiver);
    EXPECT_EQ(a.direction, b.direction);
    EXPECT_EQ(a.distance_m, b.distance_m);
    EXPECT_EQ(a.outletsea, b.outletsea);
    EXPECT_EQ(a.land, b.land);
    EXPECT_EQ(a.order, b.order);
}
}

// The mean displacement of the randomised receivers follows the gradient. Bands on
// 512 x 256: rows 50..77 (54.5 to 35.5 degrees north) and rows 22..42 (74.2 to 60.1
// degrees north). Per row, the summed displacement of the receivers (four seeds,
// about 2000 draws; standard error about 0.7 degrees at 45 N, 1.2 at 70 N) points
// within 3 degrees of the gradient azimuth (measured: at most 1.25 and 2.29), and
// the band mean of the signed deviations is within 1 degree (0.22 and 0.48; part of
// it is the analytic azimuth against the facet's, which places the diagonal
// neighbour at the mean dx of the two rows). The angle rule's expectation leans to
// the diagonal by 7.9 and 18.6 degrees on the band means (10.8 and 25.3 at most per
// row); D8 by 23.7 and 28.3 degrees on average per row.
TEST_F(Landscape, RoutingOptionsRandomisedDisplacementFollowsGradient)
{
    struct Band
    {
        const char* name;
        int firstrow, lastrow;
    };
    constexpr double rowtolerance = 3.0, bandtolerance = 1.0;
    for (const Band& band : { Band{ "35-55N", 50, 77 }, Band{ "60-75N", 22, 42 } })
    {
        const BandDeviations deviations = banddeviations(band.firstrow, band.lastrow);
        for (std::size_t row = 0; row < deviations.randomised.size(); row++)
            EXPECT_LT(std::fabs(deviations.randomised[row]), rowtolerance)
                << band.name << " row " << band.firstrow + static_cast<int>(row);
        EXPECT_LT(std::fabs(mean(deviations.randomised)), bandtolerance) << band.name;
        // The angle rule and D8 would fail both tests.
        EXPECT_GT(mean(deviations.anglerule), 5.0 * bandtolerance) << band.name;
        EXPECT_GT(largestmagnitude(deviations.anglerule), 3.0 * rowtolerance) << band.name;
        EXPECT_GT(meanmagnitude(deviations.d8), 5.0 * rowtolerance) << band.name;
    }
}

TEST_F(Landscape, RoutingOptionsRandomisedReceiversLowerOrSea)
{
    const landscape::SphereMetric metric(128, 64, radius);
    const Terrain terrain = azimuthplane(metric, pi / 6.0, 1.0e-4, 16);
    expectdrained(metric, terrain, route(metric, terrain, options(true, 1, false)), true);
    expectdrained(metric, terrain, route(metric, terrain, options(true, 1, true)), true);

    const landscape::SphereMetric small(16, 8, radius);
    const Terrain bowl = bowlplane(small);
    expectdrained(small, bowl, route(small, bowl, options(true, 7, false)), false);
    const Terrain plane = tiltedplane(small);
    expectdrained(small, plane, route(small, plane, options(true, 7, false)), true);
}

TEST_F(Landscape, RoutingOptionsOffMatchesDefault)
{
    const landscape::SphereMetric small(16, 8, radius);
    for (const Terrain& terrain : { bowlplane(small), shelfplane(small) })
    {
        const landscape::Depressions depressions = landscape::filldepressions(small, terrain.elevation, terrain.sea, 0.1f);
        expectsamegraph(landscape::buildflowgraph(small, terrain.elevation, terrain.sea, 0.0f, depressions),
                        landscape::buildflowgraph(small, terrain.elevation, terrain.sea, 0.0f, depressions,
                                                  options(false, 12345, false)));
    }
    const landscape::SphereMetric metric(128, 64, radius);
    const Terrain terrain = azimuthplane(metric, pi / 6.0, 1.0e-4, 16);
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    expectsamegraph(landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions),
                    landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions,
                                              options(false, 12345, false)));
}

TEST_F(Landscape, RoutingOptionsRandomisedSeedDetermines)
{
    const landscape::SphereMetric metric(128, 64, radius);
    const Terrain terrain = azimuthplane(metric, pi / 6.0, 1.0e-4, 16);
    const landscape::FlowGraph first = route(metric, terrain, options(true, 11, true));
    expectsamegraph(first, route(metric, terrain, options(true, 11, true)));
    const landscape::FlowGraph other = route(metric, terrain, options(true, 12, true));
    std::size_t differing = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        differing += first.direction[cell] != other.direction[cell] ? 1 : 0;
    EXPECT_GT(differing, metric.size() / 10) << "cells whose direction changes with the seed";
}

// The bowl's nine cells, all at the spill elevation, drain through the spill cell
// to the sea, in more than one direction (the flood order alone gives one chain).
TEST_F(Landscape, RoutingOptionsFlatBowlDrainsThroughSpill)
{
    const landscape::SphereMetric metric(16, 8, radius);
    const Terrain terrain = bowlplane(metric);
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    ASSERT_EQ(depressions.depressions.size(), 1u);
    const std::size_t spill = static_cast<std::size_t>(depressions.depressions[0].spillcell);
    for (bool randomised : { false, true })
    {
        const landscape::FlowGraph graph = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f,
                                                                     depressions, options(randomised, 3, true));
        expectdrained(metric, terrain, graph, false);
        std::set<int> directions;
        for (std::size_t cell = 0; cell < metric.size(); cell++)
        {
            if (!inbowl(metric.column(cell), metric.row(cell)))
                continue;
            directions.insert(graph.direction[cell]);
            const std::vector<std::size_t> cells = chain(graph, cell);
            ASSERT_FALSE(cells.empty()) << "cell " << cell;
            bool spilled = false;
            for (std::size_t k = 1; k < cells.size(); k++)
            {
                spilled = spilled || cells[k] == spill;
                if (!spilled)
                    EXPECT_TRUE(inbowl(metric.column(cells[k]), metric.row(cells[k])))
                        << "cell " << cell << " leaves the bowl before the spill";
            }
            EXPECT_TRUE(spilled) << "cell " << cell;
        }
        EXPECT_GE(directions.size(), 3u) << (randomised ? "randomised" : "deterministic");
    }
    // The depression itself is the kernel's: buildflowgraph takes it as const.
    EXPECT_EQ(depressions.depressions[0].spillcell, static_cast<int32_t>(spill));
}

// A shelf at its own terrain (no fill): its cells drain across it, never east, to
// the cells that step down at its lower edge (column 4).
TEST_F(Landscape, RoutingOptionsShelfDrainsToLowerEdge)
{
    const landscape::SphereMetric metric(16, 8, radius);
    const Terrain terrain = shelfplane(metric);
    const landscape::FlowGraph graph = route(metric, terrain, options(false, 0, true));
    expectdrained(metric, terrain, graph, false);
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const int x = metric.column(cell);
        if (x < 5 || x > 9)
            continue;
        // The increment depends on the column only: west is the steepest descent on it.
        EXPECT_EQ(static_cast<int>(graph.direction[cell]), 7) << "cell " << cell;
        const std::vector<std::size_t> cells = chain(graph, cell);
        ASSERT_FALSE(cells.empty()) << "cell " << cell;
        // Along the shelf the chain only moves west until it leaves it at column 4.
        bool left = false;
        for (std::size_t k = 1; k < cells.size() && !left; k++)
        {
            const int column = metric.column(cells[k]);
            EXPECT_LE(column, metric.column(cells[k - 1])) << "cell " << cell;
            left = column <= 4;
        }
        EXPECT_TRUE(left) << "cell " << cell;
    }
}

namespace
{
landscape::RoutingOptions pathoptions(double damping, bool resolveflats)
{
    landscape::RoutingOptions result;
    result.pathcorrected = true;
    result.pathdamping = damping;
    result.resolveflats = resolveflats;
    return result;
}

// A gridless 72 x 32 window of rows (H != W / 2, so no grid) with square cells:
// every row's dx set to dy (pi R / 32, 9.8 km at R = 100 km), area dx dy.
landscape::SphereMetric squaremetric()
{
    landscape::SphereMetric metric(72, 32, 1.0e5);
    for (double& dx : metric.rowdx_m)
        dx = metric.dy_m;
    for (double& area : metric.rowarea_m2)
        area = metric.dy_m * metric.dy_m;
    return metric;
}

// Seeded noise of 0..800 m on a ramp rising 30 m per cell away from a sea of 8
// cells per row (by index, so on any metric); about a tenth of the land cells
// lowered by 400 m into pits, so the fill has depressions and flats.
Terrain roughterrain(const landscape::SphereMetric& metric, std::uint32_t seed)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> noise(0.0, 800.0), pick(0.0, 1.0);
    Terrain terrain;
    terrain.elevation.resize(metric.size());
    terrain.sea.resize(metric.size());
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const std::size_t along = cell % 32;
        const bool sea = along < 3;
        terrain.sea[cell] = sea ? 1 : 0;
        double height = 30.0 * static_cast<double>(along) + noise(generator);
        if (pick(generator) < 0.1)
            height -= 400.0;
        terrain.elevation[cell] = sea ? -100.0f : static_cast<float>(std::max(height, 1.0));
    }
    return terrain;
}

// The heading (degrees anticlockwise from east) of the summed tangent offsets of
// the first `count` steps of the receiver chain from `cell` (fewer, returned in
// `steps`, when the chain reaches its outlet or row 0 first). A fixed count, not
// "until row 0", so the chord does not end on a north step by construction.
double chainheading(const landscape::SphereMetric& metric, const landscape::FlowGraph& graph, std::size_t cell,
                    int count, int& steps)
{
    double east = 0.0, north = 0.0;
    steps = 0;
    for (; steps < count && metric.row(cell) > 0;)
    {
        const landscape::TangentOffset offset = metric.tangentoffset(cell, graph.direction[cell]);
        east += offset.east;
        north += offset.north;
        steps++;
        const std::size_t next = static_cast<std::size_t>(graph.receiver[cell]);
        if (next == cell)
            break;
        cell = next;
    }
    return std::atan2(north, east) * 180.0 / pi;
}

// The path-corrected graph `on` against the D8 graph `off` on the same inputs:
// every land cell drains along a chain without a cycle (each chain ends at a root),
// `order` lists every land cell after its receiver, a changed receiver is a
// neighbour strictly lower on the filled surface or (an outlet) a sea neighbour,
// and the cells of a depression without a sea neighbour keep their receiver.
// Returns the number of changed cells.
std::size_t expectpathtree(const landscape::SphereMetric& metric, const Terrain& terrain,
                           const landscape::Depressions& depressions, const landscape::FlowGraph& on,
                           const landscape::FlowGraph& off)
{
    const std::vector<float>& filled = depressions.filledelevation;
    std::size_t changed = 0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (terrain.sea[cell] != 0)
        {
            EXPECT_EQ(on.receiver[cell], off.receiver[cell]) << "sea cell " << cell;
            continue;
        }
        EXPECT_FALSE(chain(on, cell).empty()) << "cell " << cell << " is on a cycle";
        const bool same = on.receiver[cell] == off.receiver[cell] && on.direction[cell] == off.direction[cell]
                          && on.outletsea[cell] == off.outletsea[cell];
        if (same)
            continue;
        changed++;
        bool seaneighbour = false;
        for (int direction = 1; direction <= landscape::directioncount; direction++)
        {
            const std::int32_t next = metric.neighbourid(cell, direction);
            seaneighbour = seaneighbour || (next >= 0 && terrain.sea[static_cast<std::size_t>(next)] != 0);
        }
        EXPECT_FALSE(depressions.depressionid[cell] != 0 && !seaneighbour) << "flooded cell " << cell << " changed";
        const std::int32_t target = on.direction[cell] == 0 ? -1 : metric.neighbourid(cell, on.direction[cell]);
        EXPECT_GE(target, 0) << "cell " << cell << " has no receiver direction";
        if (target < 0)
            continue;
        const std::size_t receiver = static_cast<std::size_t>(on.receiver[cell]);
        if (receiver == cell)
        {
            EXPECT_EQ(on.outletsea[cell], target) << "cell " << cell;
            EXPECT_EQ(static_cast<int>(terrain.sea[static_cast<std::size_t>(target)]), 1) << "cell " << cell;
        }
        else
        {
            EXPECT_EQ(static_cast<std::size_t>(target), receiver) << "cell " << cell;
            EXPECT_EQ(static_cast<int>(terrain.sea[receiver]), 0) << "cell " << cell;
            EXPECT_LT(filled[receiver], filled[cell]) << "cell " << cell;
        }
    }
    expectorder(on);
    return changed;
}
}

// (a) Off is the previous routing: pathcorrected false gives the default graph for
// any pathdamping, with and without the randomised receivers and the flat
// resolution, on the tilted plane, the bowl, the shelf and a rough terrain.
TEST_F(Landscape, PathRoutingOffMatchesDefault)
{
    const landscape::SphereMetric small(16, 8, radius);
    const landscape::SphereMetric metric(32, 16, radius);
    struct Case
    {
        const landscape::SphereMetric* metric;
        Terrain terrain;
    };
    const Case cases[] = { { &small, bowlplane(small) }, { &small, shelfplane(small) },
                           { &metric, azimuthplane(metric, pi / 8.0, 1.0e-4, 8) }, { &metric, roughterrain(metric, 1) } };
    for (const Case& item : cases)
    {
        const landscape::Depressions depressions
            = landscape::filldepressions(*item.metric, item.terrain.elevation, item.terrain.sea, 0.0f);
        for (bool randomised : { false, true })
            for (bool resolveflats : { false, true })
            {
                landscape::RoutingOptions plain = options(randomised, 5, resolveflats);
                landscape::RoutingOptions damped = plain;
                damped.pathdamping = 0.3;
                expectsamegraph(landscape::buildflowgraph(*item.metric, item.terrain.elevation, item.terrain.sea, 0.0f,
                                                          depressions, plain),
                                landscape::buildflowgraph(*item.metric, item.terrain.elevation, item.terrain.sea, 0.0f,
                                                          depressions, damped));
            }
        expectsamegraph(
            landscape::buildflowgraph(*item.metric, item.terrain.elevation, item.terrain.sea, 0.0f, depressions),
            landscape::buildflowgraph(*item.metric, item.terrain.elevation, item.terrain.sea, 0.0f, depressions,
                                      options(false, 0, false)));
    }
}

// Exclusive with the randomised receivers; the damping lies in [0, 1].
TEST_F(Landscape, PathRoutingRejectsInvalidOptions)
{
    const landscape::SphereMetric small(16, 8, radius);
    const Terrain terrain = tiltedplane(small);
    const landscape::Depressions depressions = landscape::filldepressions(small, terrain.elevation, terrain.sea, 0.0f);
    landscape::RoutingOptions both = pathoptions(1.0, true);
    both.randomised = true;
    EXPECT_THROW(landscape::buildflowgraph(small, terrain.elevation, terrain.sea, 0.0f, depressions, both),
                 std::invalid_argument);
    for (double damping : { -0.1, 1.5, std::nan("") })
        EXPECT_THROW(landscape::buildflowgraph(small, terrain.elevation, terrain.sea, 0.0f, depressions,
                                               pathoptions(damping, true)),
                     std::invalid_argument)
            << "damping " << damping;
    EXPECT_NO_THROW(landscape::buildflowgraph(small, terrain.elevation, terrain.sea, 0.0f, depressions,
                                              pathoptions(0.0, true)));
}

// (b) A uniform plane on square cells (squaremetric, 2304 cells) descending at 22.5
// degrees north of east, and at 22 degrees: the 64-step chords of the chains from
// the nine cells at columns 2..4 of the three southern rows head at the plane's
// azimuth with path routing (damping 1) within 1.5 degrees per chain (one lattice
// step in 64 is 0.8 degrees) and 0.75 on their mean (measured: mean 22.79 and
// 22.28, largest chain error 1.13 and 0.87); the local choice (damping 0, always the
// nearer facet neighbour E) heads at 0. D8 at 22 degrees (E steeper than NE by
// 0.7 %) locks to the axis, heading at 0; at 22.5 degrees E and NE tie and D8's
// choice is the rounding of the float heights, so it is not asserted there.
TEST_F(Landscape, PathRoutingPlaneFollowsHeading)
{
    const landscape::SphereMetric metric = squaremetric();
    constexpr int chord = 64;
    constexpr double chaintolerance = 1.5, meantolerance = 0.75;
    for (double degrees : { 22.5, 22.0 })
    {
        const Terrain terrain = azimuthplane(metric, degrees * pi / 180.0, 0.01, 0);
        const landscape::Depressions depressions
            = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.0f);
        const auto graph = [&](const landscape::RoutingOptions& routing) {
            return landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions, routing);
        };
        const landscape::FlowGraph d8 = graph(options(false, 0, true));
        const landscape::FlowGraph path = graph(pathoptions(1.0, true));
        const landscape::FlowGraph local = graph(pathoptions(0.0, true));
        EXPECT_GT(expectpathtree(metric, terrain, depressions, path, d8), 0u) << degrees;
        double largest = 0.0, sum = 0.0;
        int chains = 0;
        for (int y = metric.height - 3; y < metric.height; y++)
            for (int x = 2; x <= 4; x++)
            {
                const std::size_t cell = metric.index(x, y);
                int steps = 0;
                const double heading = chainheading(metric, path, cell, chord, steps);
                EXPECT_EQ(steps, chord) << degrees << " cell (" << x << ", " << y << ")";
                EXPECT_NEAR(heading, degrees, chaintolerance) << degrees << " cell (" << x << ", " << y << ")";
                largest = std::max(largest, std::fabs(heading - degrees));
                sum += heading;
                chains++;
                EXPECT_NEAR(chainheading(metric, local, cell, chord, steps), 0.0, 1.0e-9)
                    << degrees << " cell (" << x << ", " << y << ")";
                if (degrees < 22.5)
                    EXPECT_NEAR(chainheading(metric, d8, cell, chord, steps), 0.0, 1.0e-9)
                        << degrees << " cell (" << x << ", " << y << ")";
            }
        const double meanheading = sum / chains;
        EXPECT_NEAR(meanheading, degrees, meantolerance) << degrees;
        std::cout << "PathRoutingPlaneFollowsHeading: azimuth " << degrees << ", mean heading " << meanheading
                  << ", largest chain error " << largest << " degrees\n";
    }
}

// (c) On a rough terrain with pits (seeded, fills with minimum depth 0 and 0.1 m, with and
// without the flat resolution, damping 1 and 0.5) on the 32 x 16 latitude-longitude
// grid: the path-corrected graph is a tree whose changed receivers are strictly
// lower neighbours or sea, flooded cells keep theirs, and some receivers change.
TEST_F(Landscape, PathRoutingReceiversLowerOrSea)
{
    const landscape::SphereMetric metric(32, 16, radius);
    for (std::uint32_t seed : { 1u, 2u })
    {
        const Terrain terrain = roughterrain(metric, seed);
        for (float minimumdepth : { 0.0f, 0.1f })
        {
            const landscape::Depressions depressions
                = landscape::filldepressions(metric, terrain.elevation, terrain.sea, minimumdepth);
            for (bool resolveflats : { false, true })
                for (double damping : { 1.0, 0.5 })
                {
                    const landscape::FlowGraph off = landscape::buildflowgraph(
                        metric, terrain.elevation, terrain.sea, 0.0f, depressions, options(false, 0, resolveflats));
                    const landscape::FlowGraph on = landscape::buildflowgraph(
                        metric, terrain.elevation, terrain.sea, 0.0f, depressions, pathoptions(damping, resolveflats));
                    EXPECT_GT(expectpathtree(metric, terrain, depressions, on, off), 0u)
                        << "seed " << seed << " minimum depth " << minimumdepth;
                    expectsamegraph(on, landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f,
                                                                  depressions, pathoptions(damping, resolveflats)));
                }
        }
    }
}

// (d) The same on the cube's tables (N 8, 384 cells; the metric's tangent offsets
// along the slots' bearings, frames rotating across face edges): a terrain falling
// away from one axis with seeded noise and pits, the sea round the opposite pole.
TEST_F(Landscape, PathRoutingCubeReceiversLowerOrSea)
{
    const landscape::SphereMetric metric(mesh::makecubegrid(8, radius));
    const mesh::Vec3 axis = mesh::normalise(mesh::Vec3{ 0.3, 0.5, 0.8 });
    std::mt19937 generator(1);
    std::uniform_real_distribution<double> noise(0.0, 60.0), pick(0.0, 1.0);
    Terrain terrain;
    terrain.elevation.resize(metric.size());
    terrain.sea.resize(metric.size());
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const double along = mesh::dot(mesh::normalise(metric.grid->centre(static_cast<mesh::CellId>(cell))), axis);
        const bool sea = along < -0.6;
        terrain.sea[cell] = sea ? 1 : 0;
        double height = 2000.0 * (1.0 + along) + noise(generator);
        if (pick(generator) < 0.08)
            height -= 150.0;
        terrain.elevation[cell] = sea ? -100.0f : static_cast<float>(std::max(height, 1.0));
    }
    std::size_t changed = 0;
    for (float minimumdepth : { 0.0f, 0.1f })
    {
        const landscape::Depressions depressions
            = landscape::filldepressions(metric, terrain.elevation, terrain.sea, minimumdepth);
        const landscape::FlowGraph off
            = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions, options(false, 0, true));
        expectsamegraph(off, landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions,
                                                       [] {
                                                           landscape::RoutingOptions routing = options(false, 0, true);
                                                           routing.pathdamping = 0.25;
                                                           return routing;
                                                       }()));
        const landscape::FlowGraph on = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f,
                                                                  depressions, pathoptions(1.0, true));
        changed += expectpathtree(metric, terrain, depressions, on, off);
    }
    EXPECT_GT(changed, 10u);
}

// (e) A flat lake whose spill has a lower index than the lake's cells: on the
// 22.5-degree plane of (b), a 3 x 3 block at columns bx..bx + 2, rows 20..22 lowered
// by 3000 m (bx 10..15, six placements); the fill (minimum depth 0) raises it to
// exactly the spill's elevation, the spill being the block's north-east diagonal
// (bx + 3, row 19), and the lake drains into it through the block's corner
// (bx + 2, row 20), a donor at the spill's filled elevation and a higher index; the
// spill's west neighbour, on the plane, is a tributary in some placements.
// pathcorrectionorder lists every donor before its receiver, so the area summed
// along it equals accumulateflow's at every cell (relative 1e-12); the lake is the
// spill's largest donor and passes 0, so the spill takes the local choice, E (the
// nearer facet neighbour), in every placement. Visiting the spill before the
// corner (increasing index at one elevation) let it inherit the tributary's
// deviation and dropped the lake's area below it.
TEST_F(Landscape, PathRoutingLakeDonorsBeforeSpill)
{
    const landscape::SphereMetric metric = squaremetric();
    const int by = 20;
    std::size_t tributaries = 0;
    for (int bx = 10; bx <= 15; bx++)
    {
        Terrain terrain = azimuthplane(metric, pi / 8.0, 0.01, 0);
        for (int y = by; y <= by + 2; y++)
            for (int x = bx; x <= bx + 2; x++)
                terrain.elevation[metric.index(x, y)] -= 3000.0f;
        const landscape::Depressions depressions
            = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.0f);
        const std::vector<float>& filled = depressions.filledelevation;
        const std::size_t spill = metric.index(bx + 3, by - 1), corner = metric.index(bx + 2, by);
        ASSERT_EQ(depressions.depressions.size(), 1u) << "bx " << bx;
        EXPECT_EQ(depressions.depressions[0].spillcell, static_cast<int32_t>(spill)) << "bx " << bx;
        const landscape::FlowGraph off = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f,
                                                                   depressions, options(false, 0, true));
        const landscape::FlowGraph on = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f,
                                                                  depressions, pathoptions(1.0, true));
        expectpathtree(metric, terrain, depressions, on, off);
        ASSERT_EQ(on.receiver[corner], static_cast<int32_t>(spill)) << "bx " << bx;
        EXPECT_EQ(filled[corner], filled[spill]) << "bx " << bx;
        EXPECT_GT(corner, spill);

        // Every donor before its receiver, and the area along the order is accumulateflow's.
        const std::vector<int32_t> order = landscape::pathcorrectionorder(on, filled);
        std::vector<int> position(metric.size(), -1);
        for (std::size_t k = 0; k < order.size(); k++)
            position[static_cast<std::size_t>(order[k])] = static_cast<int>(k);
        std::vector<double> area(metric.size(), 0.0);
        for (int32_t member : order)
            area[static_cast<std::size_t>(member)] = metric.area_m2(static_cast<std::size_t>(member));
        for (int32_t member : order)
        {
            const std::size_t cell = static_cast<std::size_t>(member);
            const std::size_t receiver = static_cast<std::size_t>(on.receiver[cell]);
            if (receiver == cell)
                continue;
            EXPECT_LT(position[cell], position[receiver]) << "bx " << bx << " cell " << cell;
            area[receiver] += area[cell];
        }
        landscape::FlowGraph accumulated = on;
        landscape::accumulateflow(metric, std::vector<float>(metric.size(), 1.0f), accumulated);
        for (std::size_t cell = 0; cell < metric.size(); cell++)
            EXPECT_NEAR(area[cell], accumulated.area_m2[cell], 1.0e-12 * accumulated.area_m2[cell])
                << "bx " << bx << " cell " << cell;

        // The lake (9 cells and the plane above it) is the spill's largest donor; the spill goes east.
        std::size_t largestdonor = spill;
        int donors = 0;
        for (std::size_t cell = 0; cell < metric.size(); cell++)
        {
            if (terrain.sea[cell] != 0 || static_cast<std::size_t>(on.receiver[cell]) != spill || cell == spill)
                continue;
            donors++;
            if (largestdonor == spill || accumulated.area_m2[cell] > accumulated.area_m2[largestdonor])
                largestdonor = cell;
        }
        EXPECT_EQ(largestdonor, corner) << "bx " << bx;
        EXPECT_GE(accumulated.area_m2[corner], 9.0 * metric.area_m2(corner)) << "bx " << bx;
        EXPECT_EQ(static_cast<int>(on.direction[spill]), 3) << "bx " << bx;
        tributaries += donors > 1 ? 1 : 0;
    }
    EXPECT_GT(tributaries, 0u);
}
