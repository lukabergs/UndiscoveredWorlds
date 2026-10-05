#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_world_grid.hpp"

TEST(TectonicWorldGrid, IdentityPreservesFieldsAndRegistration) {
    tectonics::WorldGrid grid(4, 2, 4, 2);
    const float heights[]{0, 10, 20, 30, 100, 110, 120, 130};
    for (unsigned y = 0; y < 2; ++y) for (unsigned x = 0; x < 4; ++x) {
        EXPECT_EQ(grid.nativeIndex(x, y), y * 4 + x);
        EXPECT_FLOAT_EQ(grid.height(heights, x, y), heights[y * 4 + x]);
        EXPECT_FLOAT_EQ(grid.worldCoordinate(static_cast<float>(x)), static_cast<float>(x));
    }
}

TEST(TectonicWorldGrid, ResamplingWrapsLongitudeAndClampsPolarRows) {
    tectonics::WorldGrid grid(4, 2, 8, 4);
    const float heights[]{0, 10, 20, 30, 100, 110, 120, 130};
    EXPECT_DOUBLE_EQ(grid.scale(), 2);
    EXPECT_EQ(grid.nativeIndex(0, 0), 0u);
    EXPECT_EQ(grid.nativeIndex(7, 3), 7u);
    EXPECT_FLOAT_EQ(grid.height(heights, 0, 0), 7.5f);
    EXPECT_FLOAT_EQ(grid.height(heights, 7, 3), 122.5f);
    EXPECT_FLOAT_EQ(grid.height(heights, 2, 1), 32.5f);
    EXPECT_FLOAT_EQ(grid.worldCoordinate(0), .5f);
    EXPECT_THROW((tectonics::WorldGrid(0, 2, 8, 4)), std::invalid_argument);
    EXPECT_THROW((tectonics::WorldGrid(4, 2, 8, 3)), std::invalid_argument);
}
