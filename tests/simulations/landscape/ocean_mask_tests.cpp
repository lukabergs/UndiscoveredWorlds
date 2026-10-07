// Landscape.OceanMask*: the sea of a below-sea mask
// (docs/architecture/erosion-hydrology-plan.md, Phase H, step H3).
//
// On the 16 x 8 grid a cell's area is proportional to sin(north edge) - sin(south
// edge) of its row: 0.076 in the polar rows (0 and 7), 0.217 (rows 1 and 6), 0.324
// (rows 2 and 5) and 0.383 (rows 3 and 4); the planet sums to 16 x 2 = 32.

#include "landscape_fixture.hpp"
#include "simulations/landscape/ocean_mask.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int width = 16;
constexpr int height = 8;
constexpr double radius = 6371000.0;
constexpr double largestonly = 1.0;  // a fraction no region reaches: only the largest is sea

landscape::SphereMetric grid() { return landscape::SphereMetric(width, height, radius); }

// Marks the cells of columns x0..x1 and rows y0..y1 (inclusive) below the sea.
void mark(const landscape::SphereMetric& metric, std::vector<uint8_t>& below, int x0, int x1, int y0, int y1)
{
    for (int y = y0; y <= y1; y++)
        for (int x = x0; x <= x1; x++)
            below[metric.index(x, y)] = 1;
}
}

// An ocean over the three northern rows and a two-cell basin enclosed by land (1.4
// percent of the planet): under a 5 percent threshold the basin is not sea, every
// ocean cell is; with a threshold of 0 the basin is sea too.
TEST_F(Landscape, OceanMaskEnclosedBasinIsNotSea)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, width - 1, 0, 2);
    mark(metric, below, 8, 9, 6, 6);
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, 0.05);
    for (int x = 0; x < width; x++)
        for (int y = 0; y <= 2; y++)
            EXPECT_EQ(ocean[metric.index(x, y)], 1) << x << "," << y;
    EXPECT_EQ(ocean[metric.index(8, 6)], 0);
    EXPECT_EQ(ocean[metric.index(9, 6)], 0);
    for (int x = 0; x < width; x++)
        for (int y = 3; y < height; y++)
            if (!(y == 6 && (x == 8 || x == 9)))
                EXPECT_EQ(ocean[metric.index(x, y)], 0) << x << "," << y;
    EXPECT_EQ(landscape::oceanmask(metric, below, 0.0), below);
}

// A basin cell that touches the ocean only at a corner is part of it (8-connected).
TEST_F(Landscape, OceanMaskDiagonalTouchIsSea)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, width - 1, 0, 2);
    mark(metric, below, 4, 4, 3, 3);  // a bay below the ocean
    mark(metric, below, 5, 5, 4, 4);  // its only below-sea neighbour is (4, 3), a corner
    mark(metric, below, 6, 6, 5, 5);  // and so on, one more corner
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, largestonly);
    EXPECT_EQ(ocean[metric.index(4, 3)], 1);
    EXPECT_EQ(ocean[metric.index(5, 4)], 1);
    EXPECT_EQ(ocean[metric.index(6, 5)], 1);
    EXPECT_EQ(ocean, below);
}

// Two blocks on either side of the date line (areas 4.24 and 2.83) are one region of
// 7.07 through the longitude wrap, larger than a block of 5.45 in the middle: the
// wrapped region is the ocean, the middle block is not.
TEST_F(Landscape, OceanMaskWrapsLongitude)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, 2, 2, 5);
    mark(metric, below, width - 2, width - 1, 2, 5);
    mark(metric, below, 6, 10, 2, 4);
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, largestonly);
    for (int y = 2; y <= 5; y++)
    {
        for (int x : { 0, 1, 2, width - 2, width - 1 })
            EXPECT_EQ(ocean[metric.index(x, y)], 1) << x << "," << y;
        for (int x = 6; x <= 10; x++)
            EXPECT_EQ(ocean[metric.index(x, y)], 0) << x << "," << y;
    }
}

// The largest region is the largest by area: a polar sea of 16 cells (1.22) loses
// to a low-latitude sea of 4 cells (1.53).
TEST_F(Landscape, OceanMaskLargestByArea)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, width - 1, 0, 0);
    mark(metric, below, 6, 9, 3, 3);
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, largestonly);
    for (int x = 0; x < width; x++)
        EXPECT_EQ(ocean[metric.index(x, 0)], 0) << x;
    for (int x = 6; x <= 9; x++)
        EXPECT_EQ(ocean[metric.index(x, 3)], 1) << x;
}

// Under a 10 percent threshold two separate oceans (26.5 and 22.1 percent) both stay
// sea; a one-cell enclosed sea (0.7 percent) does not.
TEST_F(Landscape, OceanMaskKeepsEveryRegionAboveTheThreshold)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, 5, 2, 5);
    mark(metric, below, 9, 13, 2, 5);
    mark(metric, below, 7, 7, 6, 6);
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, 0.1);
    for (int y = 2; y <= 5; y++)
    {
        for (int x = 0; x <= 5; x++)
            EXPECT_EQ(ocean[metric.index(x, y)], 1) << x << "," << y;
        for (int x = 9; x <= 13; x++)
            EXPECT_EQ(ocean[metric.index(x, y)], 1) << x << "," << y;
    }
    EXPECT_EQ(ocean[metric.index(7, 6)], 0);
}

// No neighbour beyond a pole: two cells of the polar row on opposite meridians are
// separate regions of equal area; the first in cell order (the lower index) wins the tie.
TEST_F(Landscape, OceanMaskDoesNotCrossThePole)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    below[metric.index(0, 0)] = 1;
    below[metric.index(width / 2, 0)] = 1;
    const std::vector<uint8_t> ocean = landscape::oceanmask(metric, below, largestonly);
    EXPECT_EQ(ocean[metric.index(0, 0)], 1);
    EXPECT_EQ(ocean[metric.index(width / 2, 0)], 0);
}

// One below-sea region comes back unchanged; no below-sea cell gives no sea; a mask
// of the wrong size is rejected.
TEST_F(Landscape, OceanMaskOneRegionUnchangedAndEdgeCases)
{
    const landscape::SphereMetric metric = grid();
    std::vector<uint8_t> below(metric.size(), 0);
    mark(metric, below, 0, width - 1, 0, 3);
    mark(metric, below, 3, 7, 4, 6);
    EXPECT_EQ(landscape::oceanmask(metric, below, 0.002), below);

    const std::vector<uint8_t> land(metric.size(), 0);
    EXPECT_EQ(landscape::oceanmask(metric, land, 0.002), land);

    EXPECT_THROW(landscape::oceanmask(metric, std::vector<uint8_t>(metric.size() + 1, 0), 0.002), std::invalid_argument);
}
