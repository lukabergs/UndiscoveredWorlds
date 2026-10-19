// Tests of landscape/regional_grid.hpp (erosion plan Phase E, step E2a), suite Landscape:
// Landscape.RegionalGrid* (the window) and Landscape.PatchMetric* (E0's patch metric).

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
}

// A 16 x 8 world: columns of 22.5 degrees from -180, rows of 22.5 degrees from 90.
TEST_F(Landscape, RegionalGridSnapsOutwardsToWholeCells)
{
    // Longitudes 10 .. 30 lie in columns 8 (0 .. 22.5) and 9 (22.5 .. 45); latitudes
    // 5 .. 20 in row 3 (22.5 .. 0).
    const landscape::RegionalWindow window = landscape::regionalwindow(16, 4, 2, 10.0, 5.0, 30.0, 20.0);
    EXPECT_EQ(window.worldwidth, 16);
    EXPECT_EQ(window.worldheight, 8);
    EXPECT_EQ(window.column0, 8);
    EXPECT_EQ(window.columns, 2);
    EXPECT_EQ(window.row0, 3);
    EXPECT_EQ(window.rows, 1);
    EXPECT_EQ(window.domaincolumn0, 6);
    EXPECT_EQ(window.domaincolumns, 6);
    EXPECT_EQ(window.domainrow0, 1);
    EXPECT_EQ(window.domainrows, 5);
    EXPECT_EQ(window.finewidth(), 8);
    EXPECT_EQ(window.fineheight(), 4);
    EXPECT_EQ(window.domainfinewidth(), 24);
    EXPECT_EQ(window.domainfineheight(), 20);
    EXPECT_EQ(window.keptfinecolumnoffset(), 8);
    EXPECT_EQ(window.keptfinerowoffset(), 8);

    // A box on cell edges is not widened: -180 .. -135 is columns 0 and 1, 90 .. 45 rows 0 and 1.
    const landscape::RegionalWindow edges = landscape::regionalwindow(16, 1, 0, -180.0, 45.0, -135.0, 90.0);
    EXPECT_EQ(edges.column0, 0);
    EXPECT_EQ(edges.columns, 2);
    EXPECT_EQ(edges.row0, 0);
    EXPECT_EQ(edges.rows, 2);
}

TEST_F(Landscape, RegionalGridCrossesTheDateLine)
{
    // West 170 lies in column 15, east -170 in column 0.
    const landscape::RegionalWindow window = landscape::regionalwindow(16, 2, 0, 170.0, -10.0, -170.0, 10.0);
    EXPECT_EQ(window.column0, 15);
    EXPECT_EQ(window.columns, 2);
    EXPECT_EQ(window.domaincolumn0, 15);
    EXPECT_EQ(window.row0, 3);
    EXPECT_EQ(window.rows, 2);
    EXPECT_EQ(window.absolutefinecolumn(0), 30);
    EXPECT_EQ(window.absolutefinecolumn(1), 31);
    EXPECT_EQ(window.absolutefinecolumn(2), 0);
    EXPECT_EQ(window.absolutefinecolumn(3), 1);
    EXPECT_EQ(ancestorcolumn(window, 1), 15);
    EXPECT_EQ(ancestorcolumn(window, 2), 0);

    // With a halo of 3 the domain starts at column 12 and the kept box 3 columns in.
    const landscape::RegionalWindow haloed = landscape::regionalwindow(16, 2, 3, 170.0, -10.0, -170.0, 10.0);
    EXPECT_EQ(haloed.domaincolumn0, 12);
    EXPECT_EQ(haloed.domaincolumns, 8);
    EXPECT_EQ(haloed.keptfinecolumnoffset(), 6);
    EXPECT_EQ(haloed.absolutefinecolumn(0), 24);
    EXPECT_EQ(haloed.absolutefinecolumn(haloed.domainfinewidth() - 1), 7);
    EXPECT_EQ(ancestorcolumn(haloed, haloed.domainfinewidth() - 1), 3);
}

TEST_F(Landscape, RegionalGridLongitudesWrapWithinATurn)
{
    // -350 .. -330 is 10 .. 30: columns 8 and 9.
    const landscape::RegionalWindow shifted = landscape::regionalwindow(16, 2, 1, -350.0, 5.0, -330.0, 20.0);
    EXPECT_EQ(shifted.column0, 8);
    EXPECT_EQ(shifted.columns, 2);
    // -180 .. 180 and 0 .. 360 are the whole circle.
    const landscape::RegionalWindow whole = landscape::regionalwindow(16, 2, 0, -180.0, 5.0, 180.0, 20.0);
    EXPECT_EQ(whole.column0, 0);
    EXPECT_EQ(whole.columns, 16);
    const landscape::RegionalWindow turn = landscape::regionalwindow(16, 2, 0, 0.0, 5.0, 360.0, 20.0);
    EXPECT_EQ(turn.column0, 8);
    EXPECT_EQ(turn.columns, 16);
    // 180 .. -179 crosses the date line into column 0 only.
    const landscape::RegionalWindow sliver = landscape::regionalwindow(16, 2, 0, 180.0, 5.0, -179.0, 20.0);
    EXPECT_EQ(sliver.column0, 0);
    EXPECT_EQ(sliver.columns, 1);
}

TEST_F(Landscape, RegionalGridClipsTheHaloAtThePoles)
{
    // Latitudes 70 .. 85 lie in row 0: the halo of 2 has no rows to the north.
    const landscape::RegionalWindow north = landscape::regionalwindow(16, 3, 2, 0.0, 70.0, 20.0, 85.0);
    EXPECT_EQ(north.row0, 0);
    EXPECT_EQ(north.rows, 1);
    EXPECT_EQ(north.domainrow0, 0);
    EXPECT_EQ(north.domainrows, 3);
    EXPECT_EQ(north.keptfinerowoffset(), 0);
    EXPECT_EQ(north.absolutefinerow(0), 0);

    // Latitudes -89 .. -80 lie in row 7: the halo has no rows to the south.
    const landscape::RegionalWindow south = landscape::regionalwindow(16, 3, 2, 0.0, -89.0, 20.0, -80.0);
    EXPECT_EQ(south.row0, 7);
    EXPECT_EQ(south.rows, 1);
    EXPECT_EQ(south.domainrow0, 5);
    EXPECT_EQ(south.domainrows, 3);
    EXPECT_EQ(south.keptfinerowoffset(), 6);
    EXPECT_EQ(south.absolutefinerow(south.domainfineheight() - 1), 23);
    EXPECT_EQ(ancestorrow(south, south.domainfineheight() - 1), 7);

    // Columns are never clipped.
    EXPECT_EQ(north.domaincolumns, north.columns + 4);
}

TEST_F(Landscape, RegionalGridFineCentresAverageToTheGlobalCentre)
{
    // A window over global columns 14, 15, 0 (crossing the date line) and rows 2 .. 4.
    constexpr int k = 4;
    const landscape::RegionalWindow window = landscape::regionalwindow(16, k, 0, 140.0, -20.0, -170.0, 30.0);
    ASSERT_EQ(window.domaincolumn0, 14);
    ASSERT_EQ(window.domaincolumns, 3);
    ASSERT_EQ(window.domainrow0, 2);
    ASSERT_EQ(window.domainrows, 3);
    for (int dr = 0; dr < window.domainrows; dr++)
    {
        const int j = window.domainrow0 + dr;
        double latitude = 0.0;
        for (int y = dr * k; y < (dr + 1) * k; y++)
        {
            EXPECT_EQ(ancestorrow(window, y), j);
            latitude += window.finelatitude_rad(y);
        }
        const double north = pi / 2.0 - pi * j / 8.0, south = pi / 2.0 - pi * (j + 1) / 8.0;
        EXPECT_NEAR(latitude / k, 0.5 * (north + south), 1e-12) << "row " << j;
    }
    for (int dc = 0; dc < window.domaincolumns; dc++)
    {
        const int i = (window.domaincolumn0 + dc) % 16;
        double longitude = 0.0;
        for (int x = dc * k; x < (dc + 1) * k; x++)
        {
            EXPECT_EQ(ancestorcolumn(window, x), i);
            longitude += window.finelongitude_rad(x);
        }
        EXPECT_NEAR(longitude / k, 2.0 * pi * (i + 0.5) / 16.0, 1e-12) << "column " << i;
    }
}

TEST_F(Landscape, RegionalGridThrowsOnBadInput)
{
    EXPECT_THROW(landscape::regionalwindow(16, 2, 1, 10.0, 0.0, 10.0, 20.0), std::invalid_argument);    // west == east
    EXPECT_THROW(landscape::regionalwindow(16, 2, 1, 0.0, 20.0, 10.0, 20.0), std::invalid_argument);    // south == north
    EXPECT_THROW(landscape::regionalwindow(16, 2, 1, 0.0, 30.0, 10.0, 20.0), std::invalid_argument);    // south > north
    EXPECT_THROW(landscape::regionalwindow(16, 2, 1, 0.0, -91.0, 10.0, 20.0), std::invalid_argument);   // south < -90
    EXPECT_THROW(landscape::regionalwindow(16, 2, 1, 0.0, 0.0, 10.0, 91.0), std::invalid_argument);     // north > 90
    EXPECT_THROW(landscape::regionalwindow(16, 0, 1, 0.0, 0.0, 10.0, 20.0), std::invalid_argument);     // scale < 1
    EXPECT_THROW(landscape::regionalwindow(16, 2, -1, 0.0, 0.0, 10.0, 20.0), std::invalid_argument);    // halo < 0
    EXPECT_THROW(landscape::regionalwindow(2, 2, 0, 0.0, 0.0, 10.0, 20.0), std::invalid_argument);      // width < 4
    EXPECT_THROW(landscape::regionalwindow(15, 2, 0, 0.0, 0.0, 10.0, 20.0), std::invalid_argument);     // odd width
    EXPECT_THROW(landscape::regionalwindow(16, 2, 0, 180.0, 0.0, -180.0, 20.0), std::invalid_argument); // one meridian
    EXPECT_THROW(landscape::regionalwindow(16, 2, 0, -360.5, 0.0, 10.0, 20.0), std::invalid_argument);  // |west| > 360
    EXPECT_THROW(landscape::regionalwindow(16, 2, 0, 0.0, 0.0, 1e300, 20.0), std::invalid_argument);    // |east| > 360
    // Columns 8 .. 17 of 16 (10 columns) with a halo of 4: 18 > 16.
    EXPECT_THROW(landscape::regionalwindow(16, 2, 4, 1.0, 0.0, -136.0, 20.0), std::invalid_argument);
    EXPECT_NO_THROW(landscape::regionalwindow(16, 2, 3, 1.0, 0.0, -136.0, 20.0));  // 16 columns: the whole turn
}

TEST_F(Landscape, PatchMetricOfTheWholeSphereSumsToItsArea)
{
    const landscape::RegionalWindow window = landscape::regionalwindow(16, 2, 0, -180.0, -90.0, 180.0, 90.0);
    ASSERT_EQ(window.domaincolumns, 16);
    ASSERT_EQ(window.domainrows, 8);
    const landscape::PatchMetric metric = landscape::finepatchmetric(window, radius);
    ASSERT_EQ(metric.width, 32);
    ASSERT_EQ(metric.height, 16);
    double total = 0.0;
    for (double area : metric.rowarea_m2)
        total += area * metric.width;
    const double sphere = 4.0 * pi * radius * radius;
    EXPECT_NEAR(total / sphere, 1.0, 1e-9);
}

TEST_F(Landscape, PatchMetricNeighbourBeyondTheEdgeIsTheCell)
{
    const landscape::RegionalWindow window = landscape::regionalwindow(16, 2, 1, 10.0, 5.0, 30.0, 20.0);
    const landscape::PatchMetric metric = landscape::finepatchmetric(window, radius);
    const std::size_t northwest = metric.index(0, 0);
    for (int direction : { 1, 2, 6, 7, 8 })
        EXPECT_EQ(metric.neighbour(northwest, direction), northwest) << "direction " << direction;
    EXPECT_EQ(metric.neighbour(northwest, 3), metric.index(1, 0));
    EXPECT_EQ(metric.neighbour(northwest, 5), metric.index(0, 1));
    const std::size_t southeast = metric.index(metric.width - 1, metric.height - 1);
    for (int direction : { 2, 3, 4, 5, 6 })
        EXPECT_EQ(metric.neighbour(southeast, direction), southeast) << "direction " << direction;
    EXPECT_EQ(metric.neighbour(southeast, 7), metric.index(metric.width - 2, metric.height - 1));
}

TEST_F(Landscape, PatchMetricRowsEqualThoseOfThePlanetWideFineGrid)
{
    constexpr int k = 3;
    const landscape::RegionalWindow window = landscape::regionalwindow(16, k, 1, 10.0, -30.0, 30.0, 20.0);
    ASSERT_GT(window.domainrow0, 0);
    const landscape::PatchMetric metric = landscape::finepatchmetric(window, radius);
    const landscape::SphereMetric planet(16 * k, 8 * k, radius);
    EXPECT_EQ(metric.dy_m, planet.dy_m);
    for (int y = 0; y < metric.height; y++)
    {
        const std::size_t row = static_cast<std::size_t>(y);
        const std::size_t absolute = static_cast<std::size_t>(window.absolutefinerow(y));
        EXPECT_EQ(metric.rowdx_m[row], planet.rowdx_m[absolute]) << "row " << y;
        EXPECT_EQ(metric.rowarea_m2[row], planet.rowarea_m2[absolute]) << "row " << y;
        EXPECT_EQ(metric.rowlatitude_rad[row], planet.rowlatitude_rad[absolute]) << "row " << y;
        EXPECT_EQ(window.finelatitude_rad(y), planet.rowlatitude_rad[absolute]) << "row " << y;
    }
}
