// Landscape.RoutingOptions*: the randomised receivers and the flat resolution of
// buildflowgraph on synthetic grids (docs/architecture/erosion-hydrology-plan.md,
// Phase C, step C3).

#include "landscape_fixture.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
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
