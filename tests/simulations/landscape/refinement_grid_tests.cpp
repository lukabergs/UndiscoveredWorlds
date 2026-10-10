// Tests of landscape/refinement_grid.hpp (docs/architecture/multiscale-plan.md, section 3,
// "Tests"), suite Landscape: Landscape.RefinementGrid*. S0.4 tests the inline accessors;
// S2.1 adds snapping, the halo chain and its conversion between levels.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "refinement_test_windows.hpp"
#include "simulations/landscape/refinement_grid.hpp"

#include <cstdint>

// A 256-column world at level 10 of factor 4: 4^10 = 1048576 cells per global cell edge,
// a level grid of 268435456 x 134217728 cells, so the absolute keys exceed 2^32. The
// window crosses the date line and its domain reaches one parent cell (4 level cells)
// beyond the kept box.
TEST_F(Landscape, RefinementGridAccessorsWrapAndKeysExceed32Bits)
{
    const std::int64_t width = 268435456;
    const landscape::LevelWindow window = handlevelwindow(256, 4, 10, 1, width - 8, 16, 1000, 8);
    EXPECT_EQ(window.scale(), 1048576);
    EXPECT_EQ(window.levelwidth(), width);
    EXPECT_EQ(window.levelheight(), 134217728);
    EXPECT_EQ(window.halocells(), 4);
    EXPECT_EQ(window.domaincolumn0, width - 12);
    EXPECT_EQ(window.domaincolumns, 24);
    EXPECT_EQ(window.domainrow0, 996);
    EXPECT_EQ(window.domainrows, 16);
    EXPECT_EQ(window.keptsize(), 128u);
    EXPECT_EQ(window.domainsize(), 384u);
    EXPECT_EQ(window.keptcolumnoffset(), 4);
    EXPECT_EQ(window.keptrowoffset(), 4);
    EXPECT_EQ(window.index(3, 2), 51u);

    // The column wraps at the date line; the row does not.
    EXPECT_EQ(window.absolutecolumn(0), width - 12);
    EXPECT_EQ(window.absolutecolumn(11), width - 1);
    EXPECT_EQ(window.absolutecolumn(12), 0);
    EXPECT_EQ(window.absolutecolumn(23), 11);
    EXPECT_EQ(window.absoluterow(0), 996);
    EXPECT_EQ(window.absoluterow(15), 1011);

    // Keys: row x width + column in 64 bits, above 2^32 already at row 996.
    const std::uint64_t key = window.absolutekey(12, 0);
    EXPECT_EQ(key, 996ull * 268435456ull);
    EXPECT_GT(key, 1ull << 32);
    EXPECT_EQ(window.absolutekey(11, 15), 1011ull * 268435456ull + 268435455ull);
    EXPECT_EQ(landscape::levelkey(width - 1, 1011, width), window.absolutekey(11, 15));

    // Domain positions of absolute cells, across the date line.
    EXPECT_EQ(window.localcolumn(width - 12), 0);
    EXPECT_EQ(window.localcolumn(0), 12);
    EXPECT_TRUE(window.containscolumn(11));
    EXPECT_FALSE(window.containscolumn(12));
    EXPECT_FALSE(window.containscolumn(width - 13));
    EXPECT_EQ(window.localrow(1000), 4);
    EXPECT_TRUE(window.containsrow(996));
    EXPECT_TRUE(window.containsrow(1011));
    EXPECT_FALSE(window.containsrow(995));
    EXPECT_FALSE(window.containsrow(1012));

    // Parents at level 9 and the global cells: four level cells share a parent.
    EXPECT_EQ(window.parentcolumn(11), (width - 1) / 4);
    EXPECT_EQ(window.parentcolumn(12), 0);
    EXPECT_EQ(window.parentcolumn(15), 0);
    EXPECT_EQ(window.parentcolumn(16), 1);
    EXPECT_EQ(window.parentrow(0), 249);
    EXPECT_EQ(window.parentrow(3), 249);
    EXPECT_EQ(window.parentrow(4), 250);
    EXPECT_EQ(window.domainparentcolumns(), 6);
    EXPECT_EQ(window.domainparentrows(), 4);
    EXPECT_EQ(window.globalcolumn(11), 255);
    EXPECT_EQ(window.globalcolumn(12), 0);
    EXPECT_EQ(window.globalrow(15), 0);

    // Level and scale.
    EXPECT_EQ(landscape::levelscale(4, 0), 1);
    EXPECT_EQ(landscape::levelscale(2, 5), 32);
    EXPECT_EQ(landscape::levelofscale(1048576, 4), 10);
    EXPECT_EQ(landscape::levelofscale(1048576, 2), 20);
    EXPECT_EQ(landscape::levelofscale(1, 4), 0);
    EXPECT_EQ(landscape::levelofscale(12, 4), -1);
    EXPECT_EQ(landscape::levelofscale(8, 4), -1);

    // The halo per level in parent cells: the target's halo at level 1, plus the channels'
    // reach of 2 global cells (2 factor^(L - 1) parent cells) less the 2 already in it.
    EXPECT_EQ(landscape::levelhalo(16, 4, 0), 0);
    EXPECT_EQ(landscape::levelhalo(16, 4, 1), 16);
    EXPECT_EQ(landscape::levelhalo(16, 4, 2), 22);
    EXPECT_EQ(landscape::levelhalo(16, 4, 3), 46);
    EXPECT_EQ(landscape::levelhalo(16, 2, 3), 22);

    // At the north pole the domain's rows clip.
    const landscape::LevelWindow polar = handlevelwindow(256, 4, 2, 2, 64, 16, 0, 16);
    EXPECT_EQ(polar.domainrow0, 0);
    EXPECT_EQ(polar.domainrows, 24);
    EXPECT_EQ(polar.keptrowoffset(), 0);
}
