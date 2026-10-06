// Landscape.Distributaries*: the split of river flow over prograded land
// (docs/architecture/erosion-hydrology-plan.md, step C5b).

#include "landscape_fixture.hpp"
#include "simulations/landscape/distributaries.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace
{
using landscape::Distributaries;
using landscape::SphereMetric;

constexpr double earthradius_m = 6.371e6;
constexpr int gridwidth = 16;
constexpr int gridheight = 8;
constexpr double sealevel_m = 0.0;
constexpr int northdirection = 1;
constexpr int northeastdirection = 2;
constexpr int eastdirection = 3;
constexpr int southdirection = 5;

// Rows 0 to 2 are land off the prograded land (10 m, no flow direction), the rest is
// sea; tests mark prograded land below row 2 and one river cell on row 2.
struct Plan
{
    SphereMetric metric{ gridwidth, gridheight, earthradius_m };
    std::vector<uint8_t> sea, prograded, direction;
    std::vector<float> surface;
    std::vector<int32_t> january, july;

    Plan()
    {
        const std::size_t size = metric.size();
        sea.assign(size, 1);
        prograded.assign(size, 0);
        direction.assign(size, 0);
        surface.assign(size, -100.0f);
        january.assign(size, 0);
        july.assign(size, 0);
        for (int y = 0; y <= 2; y++)
            for (int x = 0; x < gridwidth; x++)
                setland(x, y, 10.0f, false);
    }
    std::size_t at(int x, int y) const { return metric.index(x, y); }
    void setland(int x, int y, float height, bool isprograded)
    {
        const std::size_t cell = at(x, y);
        sea[cell] = 0;
        prograded[cell] = isprograded ? 1 : 0;
        surface[cell] = height;
    }
    void setriver(int x, int y, int code, int jan, int jul)
    {
        const std::size_t cell = at(x, y);
        direction[cell] = static_cast<uint8_t>(code);
        january[cell] = jan;
        july[cell] = jul;
    }
    Distributaries split(double minimumflow) const
    {
        return landscape::splitdistributaries(metric, sea, prograded, surface, sealevel_m, direction, january, july,
                                              minimumflow);
    }
};

// A river at (5, 2) flows south onto an apex at (5, 3), 2 m high, between two
// prograded cells of the same height; the row below (columns 4 to 6) is prograded
// land 1 m high; the rest of rows 3 to 7 is sea.
Plan makefan()
{
    Plan plan;
    for (int x = 4; x <= 6; x++)
    {
        plan.setland(x, 3, 2.0f, true);
        plan.setland(x, 4, 1.0f, true);
    }
    plan.setriver(5, 2, southdirection, 100, 60);
    return plan;
}
}

TEST_F(Landscape, DistributariesConserveTheFlowThatEnters)
{
    const Plan plan = makefan();
    const Distributaries result = plan.split(1.0);
    EXPECT_DOUBLE_EQ(result.inflowjanuary, 100.0);
    EXPECT_DOUBLE_EQ(result.inflowjuly, 60.0);
    EXPECT_NEAR(result.seajanuary, 100.0, 1.0e-9);
    EXPECT_NEAR(result.seajuly, 60.0, 1.0e-9);
    EXPECT_EQ(result.strandedjanuary, 0.0);
    EXPECT_EQ(result.strandedjuly, 0.0);
    // Nothing is written off the prograded land.
    for (std::size_t cell = 0; cell < plan.metric.size(); cell++)
        if (plan.prograded[cell] == 0)
        {
            EXPECT_EQ(result.direction[cell], 0) << "cell " << cell;
            EXPECT_EQ(result.january[cell], 0) << "cell " << cell;
        }
}

TEST_F(Landscape, DistributariesPointToTheLargestSupplier)
{
    // The apex takes the river from the north; it splits evenly over the three lower
    // cells below it (the cells beside it are as high), so (5, 4) is fed from the north
    // and (4, 4) from the north-east, each a third of the flow.
    const Plan plan = makefan();
    const Distributaries result = plan.split(1.0);
    const std::size_t apex = plan.at(5, 3);
    EXPECT_EQ(result.direction[apex], northdirection);
    EXPECT_EQ(result.january[apex], 100);
    EXPECT_EQ(result.july[apex], 60);
    EXPECT_EQ(result.direction[plan.at(5, 4)], northdirection);
    EXPECT_EQ(result.direction[plan.at(4, 4)], northeastdirection);
    EXPECT_EQ(result.january[plan.at(4, 4)], 33);
    EXPECT_EQ(result.july[plan.at(4, 4)], 20);
    // The cells beside the apex receive nothing: no delta cell.
    EXPECT_EQ(result.direction[plan.at(4, 3)], 0);
    EXPECT_EQ(result.direction[plan.at(6, 3)], 0);
}

TEST_F(Landscape, DistributariesKeepTheThreshold)
{
    // The apex carries 100 and 60, a mean of 80: a delta cell at a minimum of 80, none at 81.
    const Plan plan = makefan();
    const std::size_t apex = plan.at(5, 3);
    EXPECT_EQ(plan.split(80.0).direction[apex], northdirection);
    const Distributaries above = plan.split(81.0);
    EXPECT_EQ(above.direction[apex], 0);
    EXPECT_EQ(above.january[apex], 0);
    EXPECT_EQ(above.july[apex], 0);
    EXPECT_NEAR(above.seajanuary, 100.0, 1.0e-9);  // the flow still reaches the sea
}

TEST_F(Landscape, DistributariesPassAFlatCellToItsReceiver)
{
    // A prograded cell with no lower neighbour hands everything to its routing
    // receiver, which comes after it: the apex (5, 3), walled in by land off the
    // prograded land and as high as (5, 4) below it, passes its whole flow south to
    // (5, 4), which splits it between (5, 5) and the sea.
    Plan plan;
    plan.setland(5, 3, 2.0f, true);
    plan.setland(5, 4, 2.0f, true);
    plan.setland(5, 5, 1.0f, true);
    for (int x : { 4, 6 })
        for (int y = 3; y <= 4; y++)
            plan.setland(x, y, 10.0f, false);
    plan.setriver(5, 2, southdirection, 50, 50);
    plan.direction[plan.at(5, 3)] = southdirection;
    plan.direction[plan.at(5, 4)] = southdirection;
    const Distributaries result = plan.split(1.0);
    EXPECT_NEAR(result.seajanuary, 50.0, 1.0e-9);
    EXPECT_EQ(result.strandedjanuary, 0.0);
    EXPECT_EQ(result.january[plan.at(5, 4)], 50);
    EXPECT_EQ(result.direction[plan.at(5, 4)], northdirection);
    EXPECT_EQ(result.direction[plan.at(5, 5)], northdirection);
}

TEST_F(Landscape, DistributariesCrossTheDateLine)
{
    // A river at (0, 2) flows south onto an apex at (0, 3), 2 m high; (15, 3), across
    // the date line, is prograded land 1 m high; (1, 3) is sea; rows 4 and below near
    // the apex are land off the prograded land. The apex splits 1 : 2 between (15, 3)
    // and the sea, and (15, 3) points east, back across the date line.
    Plan plan;
    plan.setland(0, 3, 2.0f, true);
    plan.setland(gridwidth - 1, 3, 1.0f, true);
    for (int x : { gridwidth - 2, gridwidth - 1, 0, 1 })
        plan.setland(x, 4, 50.0f, false);
    plan.setriver(0, 2, southdirection, 90, 30);
    const Distributaries result = plan.split(1.0);
    const std::size_t west = plan.at(gridwidth - 1, 3);
    EXPECT_EQ(result.direction[plan.at(0, 3)], northdirection);
    EXPECT_EQ(result.direction[west], eastdirection);
    EXPECT_EQ(result.january[west], 30);
    EXPECT_EQ(result.july[west], 10);
    EXPECT_NEAR(result.seajanuary, 90.0, 1.0e-9);
    EXPECT_NEAR(result.seajuly, 30.0, 1.0e-9);
}

TEST_F(Landscape, DistributariesStrandACellWithOnlyOrdinaryLandBelow)
{
    // A prograded cell (5, 3), 2 m high, fed by a river from the north; its lower
    // neighbours are land off the prograded land (1 m) and its routing receiver is
    // one of them, so it has no link: its whole inflow is stranded.
    Plan plan;
    plan.setland(5, 3, 2.0f, true);
    for (int x = 4; x <= 6; x++)
        for (int y = 3; y <= 4; y++)
            if (x != 5 || y != 3)
                plan.setland(x, y, 1.0f, false);
    plan.setriver(5, 2, southdirection, 40, 20);
    plan.direction[plan.at(5, 3)] = southdirection;
    const Distributaries result = plan.split(1.0);
    EXPECT_DOUBLE_EQ(result.inflowjanuary, 40.0);
    EXPECT_DOUBLE_EQ(result.inflowjuly, 20.0);
    EXPECT_DOUBLE_EQ(result.strandedjanuary, result.inflowjanuary);
    EXPECT_DOUBLE_EQ(result.strandedjuly, result.inflowjuly);
    EXPECT_EQ(result.seajanuary, 0.0);
    EXPECT_EQ(result.seajuly, 0.0);
    const std::size_t cell = plan.at(5, 3);
    EXPECT_EQ(result.direction[cell], 0);
    EXPECT_EQ(result.january[cell], 0);
    EXPECT_EQ(result.july[cell], 0);
}
