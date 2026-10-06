// Landscape.SphereMetric*, Landscape.FlowGraph*, Landscape.Depressions*: the
// metric, the drainage tree and the Priority-Flood kernels on synthetic grids
// (docs/architecture/erosion-hydrology-plan.md, section 3, Tests; step A1).

#include "landscape_fixture.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
constexpr int width = 16;
constexpr int height = 8;

// Column 0 is sea at -100 m; land elsewhere rises 100 m per column. Longitude
// wraps, so column 15 also borders the sea column (to the east).
struct Terrain
{
    std::vector<float> elevation;
    std::vector<uint8_t> sea;
};

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

// Follows receivers from a cell to the cell that is its own receiver; -1 when
// the chain is longer than the grid (a cycle).
int32_t chainend(const landscape::FlowGraph& graph, std::size_t cell)
{
    std::size_t current = cell;
    for (std::size_t step = 0; step <= graph.size(); step++)
    {
        const std::size_t next = static_cast<std::size_t>(graph.receiver[current]);
        if (next == current)
            return static_cast<int32_t>(current);
        current = next;
    }
    return -1;
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
}

TEST_F(Landscape, SphereMetricAreasSumToSphere)
{
    const landscape::SphereMetric metric(width, height, radius);
    const double sphere = 4.0 * pi * radius * radius;
    double total = 0.0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
        total += metric.area_m2(cell);
    EXPECT_NEAR(total, sphere, 1e-9 * sphere);
    EXPECT_NEAR(metric.totalarea_m2(), sphere, 1e-9 * sphere);
}

TEST_F(Landscape, SphereMetricEquatorSpacing)
{
    const landscape::SphereMetric metric(width, height, radius);
    const double dx = 2.0 * pi * radius / width * std::cos(pi / (2.0 * height));
    const double dy = pi * radius / height;
    EXPECT_NEAR(metric.rowdx_m[height / 2 - 1], dx, 1e-12 * dx);
    EXPECT_NEAR(metric.rowdx_m[height / 2], dx, 1e-12 * dx);
    EXPECT_NEAR(metric.dy_m, dy, 1e-12 * dy);
}

TEST_F(Landscape, SphereMetricNoNeighbourBeyondPoles)
{
    const landscape::SphereMetric metric(width, height, radius);
    for (int x = 0; x < width; x++)
    {
        const std::size_t north = metric.index(x, 0);
        const std::size_t south = metric.index(x, height - 1);
        for (int direction : { 1, 2, 8 })
            EXPECT_EQ(metric.neighbour(north, direction), north) << "column " << x << " direction " << direction;
        for (int direction : { 4, 5, 6 })
            EXPECT_EQ(metric.neighbour(south, direction), south) << "column " << x << " direction " << direction;
    }
}

TEST_F(Landscape, SphereMetricLongitudeWraps)
{
    const landscape::SphereMetric metric(width, height, radius);
    for (int y = 0; y < height; y++)
    {
        EXPECT_EQ(metric.neighbour(metric.index(0, y), 7), metric.index(width - 1, y)) << "row " << y;
        EXPECT_EQ(metric.neighbour(metric.index(width - 1, y), 3), metric.index(0, y)) << "row " << y;
    }
}

TEST_F(Landscape, SphereMetricOppositeDirectionRoundTrips)
{
    const landscape::SphereMetric metric(width, height, radius);
    EXPECT_EQ(landscape::oppositedirection(0), 0);
    const std::size_t cell = metric.index(5, 3);
    for (int direction = 1; direction <= landscape::directioncount; direction++)
    {
        const int opposite = landscape::oppositedirection(direction);
        EXPECT_NE(opposite, direction);
        EXPECT_EQ(landscape::oppositedirection(opposite), direction);
        EXPECT_EQ(landscape::directionoffsetx(opposite), -landscape::directionoffsetx(direction));
        EXPECT_EQ(landscape::directionoffsety(opposite), -landscape::directionoffsety(direction));
        EXPECT_EQ(metric.neighbour(metric.neighbour(cell, direction), opposite), cell) << "direction " << direction;
    }
}

TEST_F(Landscape, SphereMetricDiagonalDistance)
{
    const landscape::SphereMetric metric(width, height, radius);
    for (int y = 1; y < height - 1; y++)
    {
        const std::size_t cell = metric.index(4, y);
        for (int direction : { 2, 4, 6, 8 })
        {
            const int other = y + landscape::directionoffsety(direction);
            const double meandx = 0.5 * (metric.rowdx_m[static_cast<std::size_t>(y)]
                                         + metric.rowdx_m[static_cast<std::size_t>(other)]);
            const double expected = std::hypot(meandx, metric.dy_m);
            EXPECT_NEAR(metric.distance_m(cell, direction), expected, 1e-12 * expected)
                << "row " << y << " direction " << direction;
        }
    }
}

TEST_F(Landscape, FlowGraphTiltedPlaneDrainsToSea)
{
    const landscape::SphereMetric metric(width, height, radius);
    const Terrain terrain = tiltedplane(metric);
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    EXPECT_TRUE(depressions.depressions.empty());
    landscape::FlowGraph graph = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions);
    ASSERT_EQ(graph.size(), metric.size());

    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (terrain.sea[cell] != 0)
        {
            EXPECT_EQ(static_cast<int>(graph.land[cell]), 0);
            EXPECT_EQ(static_cast<int>(graph.direction[cell]), 0);
            continue;
        }
        EXPECT_EQ(static_cast<int>(graph.land[cell]), 1);
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver == cell)
        {
            // Outlets: column 1 drains west; column 15 borders the sea column
            // across the date line and drains east.
            const int x = metric.column(cell);
            ASSERT_GE(graph.outletsea[cell], 0) << "cell " << cell;
            EXPECT_EQ(metric.column(static_cast<std::size_t>(graph.outletsea[cell])), 0) << "cell " << cell;
            EXPECT_TRUE(x == 1 || x == width - 1) << "cell " << cell;
            EXPECT_EQ(static_cast<int>(graph.direction[cell]), x == 1 ? 7 : 3) << "cell " << cell;
        }
        else
        {
            EXPECT_LT(terrain.elevation[receiver], terrain.elevation[cell]) << "cell " << cell;
            EXPECT_EQ(graph.outletsea[cell], -1) << "cell " << cell;
        }
    }
    expectorder(graph);

    const std::vector<float> runoff(metric.size(), 0.5f);
    landscape::accumulateflow(metric, runoff, graph);
    for (int y = 0; y < height; y++)
    {
        // Columns 1..14 drain to the column-1 outlet, column 15 is its own outlet:
        // the land of the band (15 cells) reaches the sea.
        const double cellarea = metric.rowarea_m2[static_cast<std::size_t>(y)];
        const std::size_t west = metric.index(1, y);
        const std::size_t east = metric.index(width - 1, y);
        EXPECT_NEAR(graph.area_m2[west], 14.0 * cellarea, 1e-12 * 14.0 * cellarea) << "row " << y;
        EXPECT_NEAR(graph.area_m2[east], cellarea, 1e-12 * cellarea) << "row " << y;
        EXPECT_NEAR(graph.area_m2[west] + graph.area_m2[east], 15.0 * cellarea, 1e-12 * 15.0 * cellarea)
            << "row " << y;
        EXPECT_NEAR(graph.discharge_m3_yr[west], 0.5 * graph.area_m2[west], 1e-12 * graph.area_m2[west]);
        EXPECT_NEAR(graph.discharge_m3_yr[east], 0.5 * graph.area_m2[east], 1e-12 * graph.area_m2[east]);
        const std::size_t seacell = metric.index(0, y);
        EXPECT_EQ(graph.area_m2[seacell], 0.0);
        EXPECT_EQ(graph.discharge_m3_yr[seacell], 0.0);
    }
}

TEST_F(Landscape, DepressionsBowlFillsToSpill)
{
    const landscape::SphereMetric metric(width, height, radius);
    const Terrain terrain = bowlplane(metric);
    const landscape::Depressions result = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);

    ASSERT_EQ(result.depressions.size(), 1u);
    const landscape::Depression& depression = result.depressions[0];
    EXPECT_EQ(depression.id, 1);
    EXPECT_EQ(depression.spillelevation_m, 500.0f);

    // The rim at column 5 (500 m) is the lowest; its cells tie, so the flood
    // enters the bowl from the first one popped (lowest index) that borders it.
    int32_t expectedspill = -1;
    for (int y = 0; y < height && expectedspill < 0; y++)
    {
        const std::size_t cell = metric.index(5, y);
        for (int direction = 1; direction <= landscape::directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next != cell && inbowl(metric.column(next), metric.row(next)))
            {
                expectedspill = static_cast<int32_t>(cell);
                break;
            }
        }
    }
    ASSERT_EQ(expectedspill, static_cast<int32_t>(metric.index(5, 2)));
    EXPECT_EQ(depression.spillcell, expectedspill);
    EXPECT_EQ(result.filledelevation[static_cast<std::size_t>(expectedspill)],
              terrain.elevation[static_cast<std::size_t>(expectedspill)]);

    double area = 0.0, volume = 0.0;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const int x = metric.column(cell), y = metric.row(cell);
        if (terrain.sea[cell] != 0)
        {
            EXPECT_EQ(result.floodreceiver[cell], -1) << "cell " << cell;
            EXPECT_EQ(result.depressionid[cell], 0) << "cell " << cell;
            continue;
        }
        EXPECT_NE(result.floodreceiver[cell], -1) << "land cell " << cell << " undrained";
        if (!inbowl(x, y))
        {
            EXPECT_EQ(result.depressionid[cell], 0) << "cell " << cell;
            EXPECT_EQ(result.filledelevation[cell], terrain.elevation[cell]) << "cell " << cell;
            continue;
        }
        EXPECT_EQ(result.depressionid[cell], 1) << "cell " << cell;
        EXPECT_EQ(result.filledelevation[cell], depression.spillelevation_m) << "cell " << cell;
        const double cellarea = metric.area_m2(cell);
        area += cellarea;
        volume += (500.0 - static_cast<double>(terrain.elevation[cell])) * cellarea;

        // The flood receivers lead through the spill cell to the sea.
        bool spilled = false;
        std::size_t current = cell;
        for (std::size_t step = 0; step <= metric.size() && result.floodreceiver[current] >= 0; step++)
        {
            current = static_cast<std::size_t>(result.floodreceiver[current]);
            spilled = spilled || current == static_cast<std::size_t>(expectedspill);
        }
        EXPECT_TRUE(spilled) << "cell " << cell;
        EXPECT_EQ(static_cast<int>(terrain.sea[current]), 1) << "cell " << cell;
    }
    EXPECT_NEAR(depression.area_m2, area, 1e-12 * area);
    EXPECT_NEAR(depression.volume_m3, volume, 1e-12 * volume);
}

TEST_F(Landscape, DepressionsBowlDrainsThroughFlowGraph)
{
    const landscape::SphereMetric metric(width, height, radius);
    const Terrain terrain = bowlplane(metric);
    const landscape::Depressions depressions = landscape::filldepressions(metric, terrain.elevation, terrain.sea, 0.1f);
    const landscape::FlowGraph graph = landscape::buildflowgraph(metric, terrain.elevation, terrain.sea, 0.0f, depressions);
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        if (terrain.sea[cell] != 0)
            continue;
        EXPECT_GE(static_cast<int>(graph.direction[cell]), 1) << "cell " << cell;
        EXPECT_LE(static_cast<int>(graph.direction[cell]), 8) << "cell " << cell;
        const int32_t end = chainend(graph, cell);
        ASSERT_GE(end, 0) << "cell " << cell << " is on a cycle";
        EXPECT_GE(graph.outletsea[static_cast<std::size_t>(end)], 0) << "cell " << cell;
    }
    expectorder(graph);
}
