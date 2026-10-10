// Tests of landscape/refinement_levels.hpp (docs/architecture/multiscale-plan.md, section
// 3, "Tests"), suite Landscape: Landscape.RefinementLevels*. S0.4 tests the cell sizes and
// the level assignment; S2.3 adds the overlap, outlet-discharge and level-consistency tests
// of the loop.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "simulations/landscape/refinement_levels.hpp"

#include <stdexcept>
#include <string>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
}

// A 256-column world (128 rows): level 0 cells of pi R / 128 = 156.4 km; at factor 4 the
// levels 1 to 4 are 39.1, 9.77, 2.44 and 0.611 km.
TEST_F(Landscape, RefinementLevelsCellSizesAndAssignment)
{
    const double global_km = pi * radius / 128.0 / 1000.0;
    EXPECT_NEAR(landscape::levelcellsize_km(radius, 256, 4, 0), global_km, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius, 256, 4, 1), global_km / 4.0, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius, 256, 4, 4), global_km / 256.0, 1e-12);
    EXPECT_NEAR(landscape::levelcellsize_km(radius, 256, 2, 3), global_km / 8.0, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius, 256, 4, 2), 9.7727, 1e-3);

    // The coarsest level whose cells are at most size / detail_cells, within 1 .. levels;
    // detail_cells 2 (the default): a detail is resolved by two cells of its level.
    EXPECT_EQ(landscape::detaillevel(10.0, 2.0, radius, 256, 4, 4), 3);    // 2.44 <= 5 < 9.77
    EXPECT_EQ(landscape::detaillevel(80.0, 2.0, radius, 256, 4, 4), 1);    // 39.1 <= 40
    EXPECT_EQ(landscape::detaillevel(78.0, 2.0, radius, 256, 4, 4), 2);    // 39.1 > 39
    EXPECT_EQ(landscape::detaillevel(5.0, 1.0, radius, 256, 4, 4), 3);     // one cell per detail: 2.44 <= 5
    EXPECT_EQ(landscape::detaillevel(1000.0, 2.0, radius, 256, 4, 4), 1);  // coarser than a global cell: level 1
    EXPECT_EQ(landscape::detaillevel(0.1, 2.0, radius, 256, 4, 4), 4);     // finer than the finest level: clamped
    EXPECT_EQ(landscape::detaillevel(0.1, 2.0, radius, 256, 4, 2), 2);
    EXPECT_EQ(landscape::detaillevel(10.0, 2.0, radius, 256, 2, 8), 5);    // factor 2: 78.2, 39.1, 19.5, 9.77, 4.89
    EXPECT_THROW(landscape::detaillevel(0.0, 2.0, radius, 256, 4, 4), std::invalid_argument);
    EXPECT_THROW(landscape::detaillevel(5.0, 0.0, radius, 256, 4, 4), std::invalid_argument);
    EXPECT_THROW(landscape::detaillevel(5.0, 2.0, radius, 256, 4, 0), std::invalid_argument);
    EXPECT_THROW(landscape::levelcellsize_km(radius, 255, 4, 1), std::invalid_argument);
    EXPECT_THROW(landscape::levelcellsize_km(radius, 256, 1, 1), std::invalid_argument);

    // The steps' fixed order within a level (section 3, "Kernels").
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::noiseoctave), 0);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::presmoothing), 1);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::solve), 5);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::projection), 9);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::water), landscape::levelstepcount - 1);
    EXPECT_EQ(std::string(landscape::levelstepname(landscape::LevelStep::glacialforms)), "glacial_forms");
}
