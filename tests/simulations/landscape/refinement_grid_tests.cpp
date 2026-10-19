// Tests of landscape/refinement_grid.hpp (docs/architecture/multiscale-plan.md, section 3,
// "Tests"), suite Landscape: Landscape.RefinementGrid*. S0.4 tests the inline accessors;
// S2.1 adds snapping, the halo chain and its conversion between levels.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "refinement_test_windows.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/refinement_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
// A target window of regionalwindow's rules built by hand (global.scale is not read).
landscape::RegionalWindow handtarget(int worldwidth, int halo, int column0, int columns, int row0, int rows)
{
    landscape::RegionalWindow window;
    window.worldwidth = worldwidth;
    window.worldheight = worldwidth / 2;
    window.scale = 1;
    window.halo = halo;
    window.column0 = column0;
    window.columns = columns;
    window.row0 = row0;
    window.rows = rows;
    window.domaincolumn0 = static_cast<int>(landscape::wraplevel(column0 - halo, worldwidth));
    window.domaincolumns = columns + 2 * halo;
    window.domainrow0 = std::max(0, row0 - halo);
    window.domainrows = std::min(window.worldheight, row0 + rows + halo) - window.domainrow0;
    return window;
}

// The domain cells of `child` whose parent cell, widened by `reach` parent cells on every
// side (rows beyond a pole skipped), is not inside `parent`'s kept box, tested with
// containscolumn and containsrow of a window whose domain is that kept box.
std::size_t outsideparentkeptbox(const landscape::LevelWindow& child, const landscape::LevelWindow& parent, int reach)
{
    const landscape::LevelWindow kept = handlevelwindow(parent.worldwidth, parent.factor, parent.level, 0, parent.column0,
                                                        parent.columns, parent.row0, parent.rows);
    std::vector<bool> columnok(static_cast<std::size_t>(child.domaincolumns), true);
    std::vector<bool> rowok(static_cast<std::size_t>(child.domainrows), true);
    for (int x = 0; x < child.domaincolumns; x++)
        for (int d = -reach; d <= reach; d++)
            if (!kept.containscolumn(landscape::wraplevel(child.parentcolumn(x) + d, kept.levelwidth())))
                columnok[static_cast<std::size_t>(x)] = false;
    for (int y = 0; y < child.domainrows; y++)
        for (int d = -reach; d <= reach; d++)
        {
            const std::int64_t row = child.parentrow(y) + d;
            if (row >= 0 && row < kept.levelheight() && !kept.containsrow(row))
                rowok[static_cast<std::size_t>(y)] = false;
        }
    std::size_t outside = 0;
    for (int y = 0; y < child.domainrows; y++)
        for (int x = 0; x < child.domaincolumns; x++)
            if (!columnok[static_cast<std::size_t>(x)] || !rowok[static_cast<std::size_t>(y)])
                outside++;
    return outside;
}

void expectsamewindow(const landscape::LevelWindow& a, const landscape::LevelWindow& b)
{
    EXPECT_EQ(a.level, b.level);
    EXPECT_EQ(a.halo, b.halo);
    EXPECT_EQ(a.column0, b.column0);
    EXPECT_EQ(a.columns, b.columns);
    EXPECT_EQ(a.row0, b.row0);
    EXPECT_EQ(a.rows, b.rows);
    EXPECT_EQ(a.domaincolumn0, b.domaincolumn0);
    EXPECT_EQ(a.domaincolumns, b.domaincolumns);
    EXPECT_EQ(a.domainrow0, b.domainrow0);
    EXPECT_EQ(a.domainrows, b.domainrows);
}

void expectlevel(const landscape::HaloChainLevel& level, int keptcolumns, int keptrows, int domaincolumns, int domainrows,
                 int halo, double reach)
{
    EXPECT_EQ(level.keptcolumns, keptcolumns);
    EXPECT_EQ(level.keptrows, keptrows);
    EXPECT_EQ(level.domaincolumns, domaincolumns);
    EXPECT_EQ(level.domainrows, domainrows);
    EXPECT_EQ(level.halo, halo);
    EXPECT_DOUBLE_EQ(level.west_cells, reach);
    EXPECT_DOUBLE_EQ(level.east_cells, reach);
    EXPECT_DOUBLE_EQ(level.north_cells, reach);
    EXPECT_DOUBLE_EQ(level.south_cells, reach);
}
}

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
    EXPECT_EQ(ancestorcolumn(window, 11), 255);
    EXPECT_EQ(ancestorcolumn(window, 12), 0);
    EXPECT_EQ(ancestorrow(window, 15), 0);

    // Level and scale.
    EXPECT_EQ(landscape::levelscale(4, 0), 1);
    EXPECT_EQ(landscape::levelscale(2, 5), 32);
    EXPECT_EQ(landscape::levelofscale(1048576, 4), 10);
    EXPECT_EQ(landscape::levelofscale(1048576, 2), 20);
    EXPECT_EQ(landscape::levelofscale(1, 4), 0);
    EXPECT_EQ(landscape::levelofscale(12, 4), -1);
    EXPECT_EQ(landscape::levelofscale(8, 4), -1);

    // The halo per level in parent cells, whole global cells: the target's halo at level 1;
    // deeper, the parent-cell reaches (16 - 5 = 11) rounded up to whole global cells plus the
    // channels' 3 global cells (2 until the S3.11 confluences) and the solve's margin of 2 (S3.8 by
    // trial; S3.9's bound 4r - 2 at catchment reach r = 1, beyond the link's 2 cells; factor 4: 8
    // and 6 global cells; factor 2, level 3: 3 + 5 = 8; a target halo under 5 keeps 5 global cells
    // deeper, 4 at level 2 of halo 0).
    EXPECT_EQ(landscape::channelreachglobalcells, 3);
    EXPECT_EQ(landscape::globalreachcells, 5);
    EXPECT_EQ(landscape::solvemargin(1), 2);
    EXPECT_EQ(landscape::solvemargin(2), 6);
    EXPECT_EQ(landscape::solvemargin(0), landscape::unboundedsolvemarginglobalcells);
    EXPECT_EQ(landscape::levelhalo(16, 4, 0), 0);
    EXPECT_EQ(landscape::levelhalo(16, 4, 1), 16);
    EXPECT_EQ(landscape::levelhalo(16, 4, 2), 32);
    EXPECT_EQ(landscape::levelhalo(16, 4, 3), 96);
    EXPECT_EQ(landscape::levelhalo(16, 2, 3), 32);
    EXPECT_EQ(landscape::levelhalo(2, 4, 3), 80);
    EXPECT_EQ(landscape::levelhalo(0, 4, 2), 16);
    EXPECT_EQ(landscape::levelhalo(3, 4, 1), 3);
    // The loop's defaults (refinementhalo 31): 31, 48 and 112 parent cells (31, 12, 7 global cells).
    EXPECT_EQ(landscape::levelhalo(31, 4, 2), 48);
    EXPECT_EQ(landscape::levelhalo(31, 4, 3), 112);
    // A margin of 0: 16 - 3 = 13 rounded up plus 3 global cells (factor 4: 7 and 4 global cells).
    EXPECT_EQ(landscape::levelhalo(16, 4, 2, 0), 28);
    EXPECT_EQ(landscape::levelhalo(16, 4, 3, 0), 64);
    EXPECT_EQ(landscape::levelhalo(31, 4, 2, 0), 40);
    EXPECT_EQ(landscape::levelhalo(31, 4, 3, 0), 80);
    constexpr int unbounded = landscape::unboundedsolvemarginglobalcells;
    EXPECT_EQ(unbounded, 2);

    // At the north pole the domain's rows clip.
    const landscape::LevelWindow polar = handlevelwindow(256, 4, 2, 2, 64, 16, 0, 16);
    EXPECT_EQ(polar.domainrow0, 0);
    EXPECT_EQ(polar.domainrows, 24);
    EXPECT_EQ(polar.keptrowoffset(), 0);
}

// Snapping outwards to whole parent cells with floor division, also for negative firsts;
// a zero count stays zero at the snapped first.
TEST_F(Landscape, RefinementGridSnapsOutwardsWithNegativeFirsts)
{
    const auto snap = [](std::int64_t first, std::int64_t count, int factor)
    { return landscape::snapoutwards(landscape::LevelSpan{first, count}, factor); };
    EXPECT_EQ(snap(5, 6, 4).first, 4);
    EXPECT_EQ(snap(5, 6, 4).count, 8);
    EXPECT_EQ(snap(8, 4, 4).first, 8);
    EXPECT_EQ(snap(8, 4, 4).count, 4);
    EXPECT_EQ(snap(-3, 2, 4).first, -4);
    EXPECT_EQ(snap(-3, 2, 4).count, 4);
    EXPECT_EQ(snap(-4, 4, 4).first, -4);
    EXPECT_EQ(snap(-4, 4, 4).count, 4);
    EXPECT_EQ(snap(-5, 1, 4).first, -8);
    EXPECT_EQ(snap(-5, 1, 4).count, 4);
    EXPECT_EQ(snap(-1, 3, 2).first, -2);
    EXPECT_EQ(snap(-1, 3, 2).count, 4);
    EXPECT_EQ(snap(-7, 9, 4).first, -8);
    EXPECT_EQ(snap(-7, 9, 4).count, 12);

    EXPECT_EQ(snap(-9, 0, 4).first, -12);
    EXPECT_EQ(snap(-9, 0, 4).count, 0);
    EXPECT_EQ(snap(7, 0, 2).first, 6);
    EXPECT_EQ(snap(7, 0, 2).count, 0);
    EXPECT_EQ(snap(8, 0, 4).first, 8);
    EXPECT_EQ(snap(8, 0, 4).count, 0);

    EXPECT_THROW(snap(0, 4, 1), std::invalid_argument);
    EXPECT_THROW(snap(0, -1, 4), std::invalid_argument);
}

// The 16 x 8 world of RegionalGridSnapsOutwardsToWholeCells: kept box columns 8 .. 9,
// row 3, halo 2 global cells.
TEST_F(Landscape, RefinementGridLevelWindowScalesTheKeptBox)
{
    const landscape::RegionalWindow global = landscape::regionalwindow(16, 1, 2, 10.0, 5.0, 30.0, 20.0);
    // The hand-computed deeper levels take the unbounded solve's margin (S3.8, catchment reach 0);
    // a margin of 0 keeps the channels' 3 global cells at level 2 (S3.11; 2 before): halo 12 parent cells.
    constexpr int unbounded = landscape::unboundedsolvemarginglobalcells;
    EXPECT_EQ(landscape::levelwindow(global, 4, 2, 0).halo, 12);

    // Level 0: the kept box itself, halo 0.
    const landscape::LevelWindow zero = landscape::levelwindow(global, 4, 0);
    expectsamewindow(zero, handlevelwindow(16, 4, 0, 0, 8, 2, 3, 1));

    // Level 1 at factor 4: s 4, halo 2 parent cells (8 level cells); the global domain.
    const landscape::LevelWindow one = landscape::levelwindow(global, 4, 1);
    expectsamewindow(one, handlevelwindow(16, 4, 1, 2, 32, 8, 12, 4));
    EXPECT_EQ(one.domaincolumn0, 24);
    EXPECT_EQ(one.domaincolumns, 24);
    EXPECT_EQ(one.domainrow0, 4);
    EXPECT_EQ(one.domainrows, 20);

    // Level 2 at factor 4: s 16, halo 4 (5 + ceil((2 - 5) / 4)) = 20 parent cells (80 level
    // cells, 5 global cells: the channels' 3 and the margin's 2; the rows clip at both poles).
    const landscape::LevelWindow two = landscape::levelwindow(global, 4, 2, unbounded);
    expectsamewindow(two, handlevelwindow(16, 4, 2, 20, 128, 32, 48, 16));
    EXPECT_EQ(two.domaincolumn0, 48);
    EXPECT_EQ(two.domaincolumns, 192);
    EXPECT_EQ(two.domainrow0, 0);
    EXPECT_EQ(two.domainrows, 128);

    // Level 3 at factor 4: s 64, halo 16 (5 + ceil((2 - 5) / 16)) = 80 parent cells (320 level
    // cells, 5 global cells).
    const landscape::LevelWindow three = landscape::levelwindow(global, 4, 3, unbounded);
    expectsamewindow(three, handlevelwindow(16, 4, 3, 80, 512, 128, 192, 64));
    EXPECT_EQ(three.levelwidth(), 1024);
    EXPECT_EQ(three.domaincolumn0, 192);
    EXPECT_EQ(three.domaincolumns, 768);
    EXPECT_EQ(three.domainrow0, 0);
    EXPECT_EQ(three.domainrows, 512);

    // Level 2 at factor 2: s 4, halo 2 (5 + ceil((2 - 5) / 2)) = 8 parent cells (16 level cells).
    const landscape::LevelWindow twobytwo = landscape::levelwindow(global, 2, 2, unbounded);
    expectsamewindow(twobytwo, handlevelwindow(16, 2, 2, 8, 32, 8, 12, 4));

    // At level 1 the domain is regionalwindow's domain, also across the date line and at
    // the poles.
    const std::vector<landscape::RegionalWindow> windows = {
        global, landscape::regionalwindow(64, 1, 3, 170.0, -10.0, -170.0, 10.0),
        landscape::regionalwindow(64, 1, 3, -20.0, 70.0, 20.0, 90.0),
        landscape::regionalwindow(64, 1, 3, -20.0, -90.0, 20.0, -75.0)};
    for (const landscape::RegionalWindow& window : windows)
        for (int factor : {2, 4})
        {
            const landscape::LevelWindow level = landscape::levelwindow(window, factor, 1);
            EXPECT_EQ(level.halo, window.halo);
            EXPECT_EQ(level.domaincolumn0, static_cast<std::int64_t>(window.domaincolumn0) * factor);
            EXPECT_EQ(level.domaincolumns, window.domaincolumns * factor);
            EXPECT_EQ(level.domainrow0, static_cast<std::int64_t>(window.domainrow0) * factor);
            EXPECT_EQ(level.domainrows, window.domainrows * factor);
        }
}

// The chain's nesting: every level's domain, each cell's parent widened by the reach, lies
// inside its parent's kept box; level 0's domain is its kept box; level >= 1 kept boxes
// are whole parent cells. Windows across the date line and at both poles, factors 2 and 4.
TEST_F(Landscape, RefinementGridHaloChainNestsAcrossTheDateLineAndPoles)
{
    const landscape::RegionalWindow dateline = landscape::regionalwindow(64, 1, 3, 170.0, -10.0, -170.0, 10.0);
    EXPECT_EQ(dateline.column0, 62);
    EXPECT_EQ(dateline.columns, 4);
    EXPECT_EQ(dateline.row0, 14);
    EXPECT_EQ(dateline.rows, 4);
    const landscape::RegionalWindow north = landscape::regionalwindow(64, 1, 3, -20.0, 70.0, 20.0, 90.0);
    EXPECT_EQ(north.row0, 0);
    EXPECT_EQ(north.rows, 4);
    const landscape::RegionalWindow south = landscape::regionalwindow(64, 1, 3, -20.0, -90.0, 20.0, -75.0);
    EXPECT_EQ(south.row0, 29);
    EXPECT_EQ(south.rows, 3);
    // An odd halo, so that the parents' kept boxes need snapping.
    const landscape::RegionalWindow odd = landscape::regionalwindow(64, 1, 5, 170.0, -10.0, -170.0, 10.0);

    const int reach = 2, top = 3;
    for (const landscape::RegionalWindow* target : {&dateline, &north, &south, &odd})
        for (int factor : {2, 4})
        {
            SCOPED_TRACE(testing::Message() << "factor " << factor << ", target column0 " << target->column0
                                            << ", row0 " << target->row0 << ", halo " << target->halo);
            const std::vector<landscape::LevelWindow> chain = landscape::halochain(*target, factor, top, reach);
            ASSERT_EQ(chain.size(), static_cast<std::size_t>(top) + 1);
            expectsamewindow(chain[static_cast<std::size_t>(top)], landscape::levelwindow(*target, factor, top));
            for (int l = 0; l <= top; l++)
            {
                const landscape::LevelWindow& window = chain[static_cast<std::size_t>(l)];
                EXPECT_EQ(window.level, l);
                EXPECT_EQ(window.factor, factor);
                EXPECT_EQ(window.halo, landscape::levelhalo(target->halo, factor, l));
                if (l >= 1)
                {
                    // The kept box and the domain are whole global cells.
                    const std::int64_t scale = window.scale();
                    EXPECT_EQ(window.column0 % scale, 0);
                    EXPECT_EQ(window.columns % scale, 0);
                    EXPECT_EQ(window.row0 % scale, 0);
                    EXPECT_EQ(window.rows % scale, 0);
                    EXPECT_EQ(window.domaincolumn0 % scale, 0);
                    EXPECT_EQ(window.domaincolumns % scale, 0);
                    EXPECT_EQ(window.domainrow0 % scale, 0);
                    EXPECT_EQ(window.domainrows % scale, 0);
                    const landscape::RegionalWindow global = landscape::regionalwindowof(window);
                    EXPECT_EQ(global.domainfinewidth(), window.domaincolumns);
                    EXPECT_EQ(global.domainfineheight(), window.domainrows);
                    EXPECT_EQ(global.absolutefinecolumn(0), window.absolutecolumn(0));
                    EXPECT_EQ(global.absolutefinerow(0), window.absoluterow(0));
                    EXPECT_EQ(outsideparentkeptbox(window, chain[static_cast<std::size_t>(l) - 1], reach), 0u);
                }
                if (target == &dateline || target == &odd)
                    EXPECT_GT(window.domaincolumn0 + window.domaincolumns, window.levelwidth());
                if (target == &north)
                    EXPECT_EQ(window.domainrow0, 0);
                if (target == &south)
                    EXPECT_EQ(window.domainrow0 + window.domainrows, window.levelheight());
            }
            expectsamewindow(chain[0], handlevelwindow(64, factor, 0, 0, chain[0].column0, chain[0].columns,
                                                       chain[0].row0, chain[0].rows));
        }

    // The odd halo at factor 2, with the unbounded solve's margin (S3.8): level 3 (s 8, halo 4 (4 + ceil(1 / 4)) = 20 parent cells,
    // 5 global cells) has the domain columns 456 .. 567 and rows 72 .. 183; their parents
    // 228 .. 283 and 36 .. 91, widened by 2, are 226 .. 285 and 34 .. 93, snapped to whole
    // global cells (4 level-2 cells) 224 .. 287 and 32 .. 95; level 2's halo is
    // 2 (4 + ceil(1 / 2)) = 10 parent cells.
    const std::vector<landscape::LevelWindow> chain = landscape::halochain(odd, 2, 3, reach, landscape::unboundedsolvemarginglobalcells);
    EXPECT_EQ(chain[3].domaincolumn0, 456);
    EXPECT_EQ(chain[3].domaincolumns, 112);
    EXPECT_EQ(chain[3].domainrow0, 72);
    EXPECT_EQ(chain[3].domainrows, 112);
    EXPECT_EQ(chain[2].column0, 224);
    EXPECT_EQ(chain[2].columns, 64);
    EXPECT_EQ(chain[2].row0, 32);
    EXPECT_EQ(chain[2].rows, 64);
    EXPECT_EQ(chain[2].halo, 10);
    EXPECT_EQ(chain[2].domaincolumn0, 204);
    EXPECT_EQ(chain[2].domaincolumns, 104);

    // Pole clipping of a level window at factor 4: rows 0 .. 3 and 29 .. 31 at level 1.
    const landscape::LevelWindow northone = landscape::levelwindow(north, 4, 1);
    EXPECT_EQ(northone.row0, 0);
    EXPECT_EQ(northone.rows, 16);
    EXPECT_EQ(northone.domainrow0, 0);
    EXPECT_EQ(northone.domainrows, 28);
    const landscape::LevelWindow southone = landscape::levelwindow(south, 4, 1);
    EXPECT_EQ(southone.row0, 116);
    EXPECT_EQ(southone.rows, 12);
    EXPECT_EQ(southone.domainrow0, 104);
    EXPECT_EQ(southone.domainrows, 24);
}

TEST_F(Landscape, RefinementGridRefusesInvalidWindowsAndChains)
{
    const landscape::RegionalWindow good = handtarget(16, 2, 8, 2, 3, 1);
    EXPECT_NO_THROW(landscape::levelwindow(good, 4, 0));
    EXPECT_THROW(landscape::levelwindow(good, 1, 1), std::invalid_argument);
    EXPECT_THROW(landscape::levelwindow(good, 4, -1), std::invalid_argument);

    landscape::RegionalWindow bad = good;
    bad.worldwidth = 15;
    bad.worldheight = 7;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.worldwidth = 2;
    bad.worldheight = 1;
    bad.column0 = 0;
    bad.columns = 1;
    bad.row0 = 0;
    bad.halo = 0;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.worldheight = 9;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.column0 = 16;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.column0 = -1;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.columns = 0;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.row0 = -1;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.row0 = 7;
    bad.rows = 2;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);
    bad = good;
    bad.halo = -1;
    EXPECT_THROW(landscape::levelwindow(bad, 4, 1), std::invalid_argument);

    // A domain wider than the level grid: 6 columns and a halo of 6 on 16 columns are
    // 18 global cells at level 1; at level 2 with the unbounded solve's margin (S3.8) the halo is 4 (4 + ceil(2 / 4)) = 20 parent
    // cells, 5 global cells, so the domain just fits (256 of 256 columns), but the chain's
    // level 1 does not.
    // S3.11 (the channels keep 3 global cells): 4 columns and a halo of 7, 18 global cells at
    // level 1; at level 2 4 (5 + ceil(2 / 4)) = 24 parent cells, 6 global cells, 16 of 16.
    const landscape::RegionalWindow wide = handtarget(16, 7, 0, 4, 2, 2);
    EXPECT_THROW(landscape::levelwindow(wide, 4, 1), std::invalid_argument);
    EXPECT_EQ(landscape::levelwindow(wide, 4, 2, landscape::unboundedsolvemarginglobalcells).domaincolumns, 256);
    EXPECT_THROW(landscape::halochain(wide, 4, 2, 2, landscape::unboundedsolvemarginglobalcells), std::invalid_argument);

    // The level grid below 2^63 cells: 256 x 128 x 4^(2 L) is 2^59 at level 11 and 2^63 at 12.
    const landscape::RegionalWindow large = handtarget(256, 16, 164, 20, 48, 16);
    EXPECT_EQ(landscape::levelwindow(large, 4, 11).levelwidth(), 256ll << 22);
    EXPECT_THROW(landscape::levelwindow(large, 4, 12), std::invalid_argument);
    EXPECT_THROW(landscape::levelwindow(large, 4, 40), std::invalid_argument);

    // The chain: level >= 1, reach >= 0, and levelwindow's refusals.
    EXPECT_THROW(landscape::halochain(good, 4, 0, 2), std::invalid_argument);
    EXPECT_THROW(landscape::halochain(good, 4, 2, -1), std::invalid_argument);
    EXPECT_THROW(landscape::halochain(good, 1, 2, 2), std::invalid_argument);
    EXPECT_THROW(landscape::halochain(large, 4, 12, 2), std::invalid_argument);

    // The report: a non-empty chain whose element l is at level l (on 64 columns: since S3.8 the
    // two-level chain of `good` is wider than its 16-column world at level 0).
    const landscape::RegionalWindow roomy = handtarget(64, 2, 8, 2, 3, 1);
    EXPECT_THROW(landscape::halochainreport(roomy, {}), std::invalid_argument);
    std::vector<landscape::LevelWindow> chain = landscape::halochain(roomy, 4, 2, 2);
    EXPECT_NO_THROW(landscape::halochainreport(roomy, chain));
    std::swap(chain[1], chain[2]);
    EXPECT_THROW(landscape::halochainreport(roomy, chain), std::invalid_argument);
}

// Window A of the plan (--pixels 164,48,183,63 on a 256-column world: columns 164 .. 183,
// rows 48 .. 63, halo 16) at factor 4 with the Catmull-Rom reach of 2 and the unbounded solve's
// margin of 2 (S3.8; catchment reach 0). Hand-computed:
// at k = 16 (level 2) the halo is 4 (4 + ceil(12 / 4)) = 28 parent cells (112 level cells,
// 7 global cells; S3.8); the parents of level 2's domain (columns 2512 .. 3055, rows
// 656 .. 1135) are 628 .. 763 and 164 .. 283, widened by 2 to 626 .. 765 and 162 .. 285,
// snapped to whole global cells 624 .. 767 and 160 .. 287, so level 1 keeps 144 x 128 and
// computes 272 x 256 (halo 16 parent cells, 24 global cells beyond the target), and level 0
// keeps and computes the 72 x 68 global cells 138 .. 209 by 22 .. 89 (26 beyond the target).
TEST_F(Landscape, RefinementGridHaloChainReportOfWindowA)
{
    const landscape::RegionalWindow target = handtarget(256, 16, 164, 20, 48, 16);

    const std::vector<landscape::LevelWindow> sixteen =
        landscape::halochain(target, 4, 2, 2, landscape::unboundedsolvemarginglobalcells);
    for (std::size_t l = 1; l < sixteen.size(); l++)
        EXPECT_EQ(outsideparentkeptbox(sixteen[l], sixteen[l - 1], 2), 0u);
    const landscape::HaloChainReport sixteenreport = landscape::halochainreport(target, sixteen);
    ASSERT_EQ(sixteenreport.levels.size(), 3u);
    // S3.11: the channels keep 3 global cells (2 before): level 2's halo 4 (5 + ceil(11 / 4)) = 32
    // parent cells (8 global cells); the coarser levels reach 25 and 27 global cells.
    expectlevel(sixteenreport.levels[0], 74, 70, 74, 70, 0, 27.0);
    expectlevel(sixteenreport.levels[1], 152, 136, 280, 264, 16, 25.0);
    expectlevel(sixteenreport.levels[2], 320, 256, 576, 512, 32, 8.0);
    EXPECT_EQ(sixteenreport.domaincells, 5180 + 73920 + 294912);
    EXPECT_EQ(sixteenreport.domaincells, 374012);

    // At k = 64 (level 3) the halo is 16 (4 + ceil(12 / 16)) = 80 parent cells (320 level
    // cells, 5 global cells; 48 and 3 before S3.8); the coarser levels reach 13, 30 and 32
    // global cells beyond the target (level 2's kept box 2542 .. 3025 by 686 .. 1105 snaps to
    // 2528 .. 3039 by 672 .. 1119).
    const std::vector<landscape::LevelWindow> sixtyfour =
        landscape::halochain(target, 4, 3, 2, landscape::unboundedsolvemarginglobalcells);
    for (std::size_t l = 1; l < sixtyfour.size(); l++)
        EXPECT_EQ(outsideparentkeptbox(sixtyfour[l], sixtyfour[l - 1], 2), 0u);
    const landscape::HaloChainReport sixtyfourreport = landscape::halochainreport(target, sixtyfour);
    ASSERT_EQ(sixtyfourreport.levels.size(), 4u);
    // S3.11: level 3's halo 16 (5 + ceil(11 / 16)) = 96 parent cells (6 global cells); the
    // coarser levels reach 15, 32 and 34 global cells.
    expectlevel(sixtyfourreport.levels[0], 88, 84, 88, 84, 0, 34.0);
    expectlevel(sixtyfourreport.levels[1], 208, 192, 336, 320, 16, 32.0);
    expectlevel(sixtyfourreport.levels[2], 544, 480, 800, 736, 32, 15.0);
    expectlevel(sixtyfourreport.levels[3], 1280, 1024, 2048, 1792, 96, 6.0);
    EXPECT_EQ(sixtyfourreport.domaincells, 7392 + 107520 + 588800 + 3670016);
    EXPECT_EQ(sixtyfourreport.domaincells, 4373728);
    for (int l = 0; l <= 3; l++)
        EXPECT_EQ(sixtyfourreport.levels[static_cast<std::size_t>(l)].level, l);
}

// A level window as a RegionalWindow: at level 1 it is regionalwindow's window of the same
// box with k = factor; at levels 2 and 3 its fine grid is the level window's domain.
TEST_F(Landscape, RefinementGridRegionalWindowOfALevelWindow)
{
    const struct
    {
        int width, halo;
        double west, south, east, north;
    } boxes[] = {{16, 2, 10.0, 5.0, 30.0, 20.0},
                 {64, 3, 170.0, -10.0, -170.0, 10.0},
                 {64, 3, -20.0, 70.0, 20.0, 90.0},
                 {64, 3, -20.0, -90.0, 20.0, -75.0}};
    for (const auto& box : boxes)
        for (int factor : {2, 4})
        {
            SCOPED_TRACE(testing::Message() << "width " << box.width << ", west " << box.west << ", south " << box.south
                                            << ", factor " << factor);
            const landscape::RegionalWindow expected =
                landscape::regionalwindow(box.width, factor, box.halo, box.west, box.south, box.east, box.north);
            const landscape::RegionalWindow one = landscape::regionalwindowof(landscape::levelwindow(expected, factor, 1));
            EXPECT_EQ(one.worldwidth, expected.worldwidth);
            EXPECT_EQ(one.worldheight, expected.worldheight);
            EXPECT_EQ(one.scale, expected.scale);
            EXPECT_EQ(one.halo, expected.halo);
            EXPECT_EQ(one.column0, expected.column0);
            EXPECT_EQ(one.columns, expected.columns);
            EXPECT_EQ(one.row0, expected.row0);
            EXPECT_EQ(one.rows, expected.rows);
            EXPECT_EQ(one.domaincolumn0, expected.domaincolumn0);
            EXPECT_EQ(one.domaincolumns, expected.domaincolumns);
            EXPECT_EQ(one.domainrow0, expected.domainrow0);
            EXPECT_EQ(one.domainrows, expected.domainrows);

            for (int level : {2, 3})
            {
                const landscape::LevelWindow window = landscape::levelwindow(expected, factor, level);
                const landscape::RegionalWindow global = landscape::regionalwindowof(window);
                EXPECT_EQ(global.scale, window.scale());
                EXPECT_EQ(global.finewidth(), window.columns);
                EXPECT_EQ(global.fineheight(), window.rows);
                EXPECT_EQ(global.domainfinewidth(), window.domaincolumns);
                EXPECT_EQ(global.domainfineheight(), window.domainrows);
                EXPECT_EQ(global.keptfinecolumnoffset(), window.keptcolumnoffset());
                EXPECT_EQ(global.keptfinerowoffset(), window.keptrowoffset());
                EXPECT_EQ(global.halo * window.scale(), window.halocells());
            }
        }

    // Level 0, and a level-2 window of whole parent cells that are not whole global cells.
    EXPECT_THROW(landscape::regionalwindowof(handlevelwindow(16, 4, 0, 0, 8, 2, 3, 1)), std::invalid_argument);
    EXPECT_THROW(landscape::regionalwindowof(handlevelwindow(16, 4, 2, 4, 4, 16, 16, 16)), std::invalid_argument);
    EXPECT_THROW(landscape::regionalwindowof(handlevelwindow(16, 4, 2, 1, 16, 16, 16, 16)), std::invalid_argument);
    EXPECT_NO_THROW(landscape::regionalwindowof(handlevelwindow(16, 4, 2, 4, 16, 16, 16, 16)));
}

// --- The chain on key sets (sphere-mesh.md section 8, "Halo by rings"; T8 W2-W1b). ---

namespace
{
constexpr double keyradius = 6371000.0;

// The level-0 cells of a level window's kept box or domain (whole global cells), ascending.
std::vector<mesh::LevelKey> globalcellsof(const landscape::LevelWindow& window, bool domain)
{
    const std::int64_t s = window.scale();
    const std::int64_t column0 = (domain ? window.domaincolumn0 : window.column0) / s;
    const std::int64_t row0 = (domain ? window.domainrow0 : window.row0) / s;
    const std::int64_t columns = (domain ? window.domaincolumns : window.columns) / s;
    const std::int64_t rows = (domain ? window.domainrows : window.rows) / s;
    std::vector<mesh::LevelKey> keys;
    for (std::int64_t y = 0; y < rows; y++)
        for (std::int64_t x = 0; x < columns; x++)
            keys.push_back(landscape::levelkey(landscape::wraplevel(column0 + x, window.worldwidth), row0 + y, window.worldwidth));
    std::sort(keys.begin(), keys.end());
    return keys;
}

std::vector<mesh::LevelKey> targetcells(const landscape::RegionalWindow& target)
{
    std::vector<mesh::LevelKey> keys;
    for (int y = 0; y < target.rows; y++)
        for (int x = 0; x < target.columns; x++)
            keys.push_back(landscape::levelkey((target.column0 + x) % target.worldwidth, target.row0 + y, target.worldwidth));
    std::sort(keys.begin(), keys.end());
    return keys;
}
}

// On lat-lon the chain on keys is halochain's, window for window and cell for cell, for
// every test window (date line, both poles, an odd halo, window A, the region-base box).
TEST_F(Landscape, RefinementGridKeyChainEqualsTheHaloChainOnLatLon)
{
    const std::vector<landscape::RegionalWindow> targets = {
        landscape::regionalwindow(16, 1, 2, 10.0, 5.0, 30.0, 20.0),
        landscape::regionalwindow(64, 1, 3, 170.0, -10.0, -170.0, 10.0),
        landscape::regionalwindow(64, 1, 3, -20.0, 70.0, 20.0, 90.0),
        landscape::regionalwindow(64, 1, 3, -20.0, -90.0, 20.0, -75.0),
        landscape::regionalwindow(64, 1, 5, 170.0, -10.0, -170.0, 10.0),
        handtarget(256, 16, 164, 20, 48, 16),
        handtarget(128, 31, 102, 8, 4, 6),  // --pixels 102,4,109,9 with the default halo
        handtarget(16, 4, 3, 4, 2, 3),      // level 0 over every column (reach 2)
        handtarget(16, 6, 3, 4, 2, 3)};     // level 1 over every column (reach 0)
    int fullwidth = 0;  // windows whose domain spans every column
    for (const landscape::RegionalWindow& target : targets)
    {
        const auto base = mesh::makelatlongrid(target.worldwidth, keyradius);
        const std::vector<mesh::LevelKey> kept0 = targetcells(target);
        EXPECT_EQ(landscape::cellsinpixelbox(*base, target.worldwidth, target.column0, target.row0,
                                             (target.column0 + target.columns - 1) % target.worldwidth, target.row0 + target.rows - 1),
                  kept0);
        for (int factor : {2, 4})
            for (int top = 1; top <= 3; top++)
                for (int margin : {landscape::solvemarginglobalcells, landscape::unboundedsolvemarginglobalcells})
                for (int reach : {2, 0})
                {
                    SCOPED_TRACE(testing::Message() << "width " << target.worldwidth << ", column0 " << target.column0 << ", row0 "
                                                    << target.row0 << ", halo " << target.halo << ", factor " << factor
                                                    << ", level " << top << ", margin " << margin << ", reach " << reach);
                    std::vector<landscape::LevelWindow> chain;
                    try
                    {
                        chain = landscape::halochain(target, factor, top, reach, margin);
                    }
                    catch (const std::invalid_argument&)
                    {
                        continue;  // a domain wider than the planet: no chain to compare
                    }
                    const std::vector<landscape::KeyLevelWindow> keys =
                        landscape::keyhalochain(*base, factor, top, kept0, target.halo, reach, margin);
                    ASSERT_EQ(keys.size(), chain.size());
                    for (std::size_t l = 0; l < chain.size(); l++)
                    {
                        const landscape::LevelWindow& window = chain[l];
                        const landscape::KeyLevelWindow& keyed = keys[l];
                        EXPECT_EQ(keyed.level, window.level);
                        EXPECT_EQ(keyed.halo, window.halo);
                        if (l >= 1)
                            EXPECT_EQ(static_cast<std::int64_t>(keyed.halorings) * landscape::levelscale(factor, window.level - 1),
                                      window.halo);
                        EXPECT_EQ(keyed.keptglobal, globalcellsof(window, false)) << "level " << l;
                        EXPECT_EQ(keyed.domainglobal, globalcellsof(window, true)) << "level " << l;
                        EXPECT_EQ(keyed.keptsize(), window.keptsize());
                        EXPECT_EQ(keyed.domainsize(), window.domainsize());
                        EXPECT_EQ(keyed.keptwest, window.column0 / window.scale()) << "level " << l;
                        EXPECT_EQ(keyed.domainwest, window.domaincolumn0 / window.scale()) << "level " << l;
                        fullwidth += window.domaincolumns == window.levelwidth() ? 1 : 0;
                        if (window.domainsize() > 300000)
                            continue;
                        // The domain grid in the window's own row-major order.
                        const auto grid = landscape::keywindowgrid(keyed, keyradius);
                        ASSERT_EQ(grid->cellcount(), window.domainsize());
                        for (int y = 0; y < window.domainrows; y++)
                            for (int x = 0; x < window.domaincolumns; x++)
                                ASSERT_EQ(grid->globalkey(static_cast<mesh::CellId>(window.index(x, y))), window.absolutekey(x, y));
                    }
                }
    }
    EXPECT_GT(fullwidth, 0);  // the start of a domain over every column is halochain's too
}

// On the cube the chain nests: each level's kept box holds its domain's parents widened by
// the reach (checked on the level keys themselves, not on the level-0 rule), every domain is
// its kept box widened by its rings, for boxes on a face edge, at a cube corner and at a pole.
TEST_F(Landscape, RefinementGridKeyChainNestsOnTheCube)
{
    const int n = 8, exportcolumns = 4 * n, factor = 4, top = 2, halo = 6, reach = 2;
    const auto base = mesh::makecubegrid(n, keyradius);
    const struct
    {
        const char* name;
        int x0, y0, x1, y1;
    } boxes[] = {{"face edge 0|1", 19, 7, 20, 8}, {"corner 0|1|4", 19, 4, 20, 5}, {"north pole", 0, 0, 31, 1}, {"date line", 31, 9, 0, 10}};
    for (const auto& box : boxes)
    {
        SCOPED_TRACE(box.name);
        const std::vector<mesh::LevelKey> kept0 = landscape::cellsinpixelbox(*base, exportcolumns, box.x0, box.y0, box.x1, box.y1);
        ASSERT_FALSE(kept0.empty());
        const std::vector<landscape::KeyLevelWindow> chain = landscape::keyhalochain(*base, factor, top, kept0, halo, reach);
        ASSERT_EQ(chain.size(), static_cast<std::size_t>(top) + 1);
        EXPECT_EQ(chain[top].keptglobal, kept0);
        EXPECT_EQ(chain[0].domainglobal, chain[0].keptglobal);
        for (int l = 0; l <= top; l++)
        {
            const landscape::KeyLevelWindow& window = chain[static_cast<std::size_t>(l)];
            EXPECT_TRUE(std::includes(window.domainglobal.begin(), window.domainglobal.end(), window.keptglobal.begin(),
                                      window.keptglobal.end()));
            EXPECT_EQ(window.domainglobal, landscape::dilatecells(*base, window.keptglobal, window.halorings));
            if (l == top)
                continue;
            // The child's domain at this level, widened by the reach in this level's cells.
            landscape::KeyLevelWindow parents = window;
            parents.domainglobal = chain[static_cast<std::size_t>(l) + 1].domainglobal;
            const std::vector<mesh::LevelKey> needed = mesh::dilatekeys(parents.levelcells(), parents.domainkeys(), reach);
            const std::vector<mesh::LevelKey> kept = window.keptkeys();
            EXPECT_TRUE(std::includes(kept.begin(), kept.end(), needed.begin(), needed.end())) << "level " << l;
        }
        const auto grid = landscape::keywindowgrid(chain[top], keyradius);
        EXPECT_EQ(grid->cellcount(), chain[top].domainsize());
        EXPECT_EQ(grid->facecells(), n * 16);
    }
}

// The region-base box (--pixels 102,4,109,9, scale 4, halo 31) on the N 32 cube: cells and
// bytes per level, the rings' size in km, and the first W2-W worker's count (275 level-1
// cells whose centres lie in the box; 52,854 within 31 cells' km of them).
TEST_F(Landscape, RefinementGridKeyChainOfTheRegionBaseBoxOnTheCube)
{
    constexpr double pi = 3.14159265358979323846;
    const int n = 32, exportcolumns = 4 * n, factor = 4, halo = 31, reach = 2;
    const auto base = mesh::makecubegrid(n, keyradius);
    const std::vector<mesh::LevelKey> kept0 = landscape::cellsinpixelbox(*base, exportcolumns, 102, 4, 109, 9);
    const std::vector<landscape::KeyLevelWindow> chain = landscape::keyhalochain(*base, factor, 1, kept0, halo, reach);
    const double ringkm = 2.0 * pi * keyradius / exportcolumns / 1000.0;
    for (const landscape::KeyLevelWindow& window : chain)
    {
        const double sphere = 6.0 * static_cast<double>(window.levelcells()) * static_cast<double>(window.levelcells());
        std::cout << "level " << window.level << ": kept " << window.keptsize() << " cells, domain " << window.domainsize()
                  << " cells (" << 100.0 * static_cast<double>(window.domainsize()) / sphere << " % of the level), halo "
                  << window.halorings << " rings (" << window.halorings * ringkm << " km at the equator); bytes compact "
                  << window.domainsize() * 80 << ", stored " << window.domainsize() * 288 << "\n";
    }
    EXPECT_EQ(chain[1].keptsize(), kept0.size() * 16);

    // The first worker's numbers: level-1 centres in the box and the km halo around them.
    const auto level1 = mesh::makecubegrid(n * factor, keyradius);
    const std::vector<mesh::LevelKey> direct = landscape::cellsinpixelbox(*level1, exportcolumns, 102, 4, 109, 9);
    EXPECT_EQ(direct.size(), 275u);
    const double limit = halo * 2.0 * pi * keyradius / exportcolumns;
    std::size_t within = 0;
    for (mesh::CellId c = 0; c < level1->cellcount(); c++)
    {
        const mesh::Vec3 p = level1->centre(c);
        for (mesh::LevelKey key : direct)
        {
            const mesh::Vec3 q = level1->centre(static_cast<mesh::CellId>(key));
            if (keyradius * std::atan2(mesh::length(mesh::cross(p, q)), mesh::dot(p, q)) <= limit)
            {
                within++;
                break;
            }
        }
    }
    std::cout << "first worker's rule: " << direct.size() << " level-1 centres in the box, " << within << " level-1 cells within "
              << limit / 1000.0 << " km of them; kept_0 " << kept0.size() << " level-0 cells (" << kept0.size() * 16
              << " level-1 descendants)\n";
    EXPECT_EQ(within, 52854u);
}
