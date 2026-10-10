// Tests of landscape/refinement_grid.hpp (docs/architecture/multiscale-plan.md, section 3,
// "Tests"), suite Landscape: Landscape.RefinementGrid*. S0.4 tests the inline accessors;
// S2.1 adds snapping, the halo chain and its conversion between levels.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "refinement_test_windows.hpp"
#include "simulations/landscape/refinement_grid.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

    // The halo per level in parent cells, whole global cells: the target's halo at level 1;
    // deeper, the pre-smoothing's 14 rounded up to whole global cells plus the channels'
    // 2 global cells (factor 4: 6 and 3 global cells; factor 2, level 3: 4 + 2 = 6).
    EXPECT_EQ(landscape::levelhalo(16, 4, 0), 0);
    EXPECT_EQ(landscape::levelhalo(16, 4, 1), 16);
    EXPECT_EQ(landscape::levelhalo(16, 4, 2), 24);
    EXPECT_EQ(landscape::levelhalo(16, 4, 3), 48);
    EXPECT_EQ(landscape::levelhalo(16, 2, 3), 24);
    EXPECT_EQ(landscape::levelhalo(2, 4, 3), 32);
    EXPECT_EQ(landscape::levelhalo(0, 4, 2), 8);
    EXPECT_EQ(landscape::levelhalo(3, 4, 1), 3);

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

    // Level 2 at factor 4: s 16, halo 2 + 2 (4 - 1) = 8 parent cells (32 level cells,
    // 2 global cells).
    const landscape::LevelWindow two = landscape::levelwindow(global, 4, 2);
    expectsamewindow(two, handlevelwindow(16, 4, 2, 8, 128, 32, 48, 16));
    EXPECT_EQ(two.domaincolumn0, 96);
    EXPECT_EQ(two.domaincolumns, 96);
    EXPECT_EQ(two.domainrow0, 16);
    EXPECT_EQ(two.domainrows, 80);

    // Level 3 at factor 4: s 64, halo 2 + 2 (16 - 1) = 32 parent cells (128 level cells).
    const landscape::LevelWindow three = landscape::levelwindow(global, 4, 3);
    expectsamewindow(three, handlevelwindow(16, 4, 3, 32, 512, 128, 192, 64));
    EXPECT_EQ(three.levelwidth(), 1024);
    EXPECT_EQ(three.domaincolumn0, 384);
    EXPECT_EQ(three.domaincolumns, 384);
    EXPECT_EQ(three.domainrow0, 64);
    EXPECT_EQ(three.domainrows, 320);

    // Level 2 at factor 2: s 4, halo 2 + 2 (2 - 1) = 4 parent cells (8 level cells).
    const landscape::LevelWindow twobytwo = landscape::levelwindow(global, 2, 2);
    expectsamewindow(twobytwo, handlevelwindow(16, 2, 2, 4, 32, 8, 12, 4));

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

    // The odd halo at factor 2: level 3 (s 8, halo 4 (2 + ceil(3 / 4)) = 12 parent cells,
    // 3 global cells) has the domain columns 472 .. 551 and rows 88 .. 167; their parents
    // 236 .. 275 and 44 .. 83, widened by 2, are 234 .. 277 and 42 .. 85, snapped to whole
    // global cells (4 level-2 cells) 232 .. 279 and 40 .. 87; level 2's halo is
    // 2 (2 + ceil(3 / 2)) = 8 parent cells.
    const std::vector<landscape::LevelWindow> chain = landscape::halochain(odd, 2, 3, reach);
    EXPECT_EQ(chain[3].domaincolumn0, 472);
    EXPECT_EQ(chain[3].domaincolumns, 80);
    EXPECT_EQ(chain[3].domainrow0, 88);
    EXPECT_EQ(chain[3].domainrows, 80);
    EXPECT_EQ(chain[2].column0, 232);
    EXPECT_EQ(chain[2].columns, 48);
    EXPECT_EQ(chain[2].row0, 40);
    EXPECT_EQ(chain[2].rows, 48);
    EXPECT_EQ(chain[2].halo, 8);
    EXPECT_EQ(chain[2].domaincolumn0, 216);
    EXPECT_EQ(chain[2].domaincolumns, 80);

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

    // A domain wider than the level grid: 10 columns and a halo of 4 on 16 columns are
    // 18 global cells at level 1; at level 2 the halo is 4 (2 + ceil(2 / 4)) = 12 parent
    // cells, 3 global cells, so the domain just fits (256 of 256 columns), but the chain's
    // level 1 does not.
    const landscape::RegionalWindow wide = handtarget(16, 4, 0, 10, 2, 2);
    EXPECT_THROW(landscape::levelwindow(wide, 4, 1), std::invalid_argument);
    EXPECT_EQ(landscape::levelwindow(wide, 4, 2).domaincolumns, 256);
    EXPECT_THROW(landscape::halochain(wide, 4, 2, 2), std::invalid_argument);

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

    // The report: a non-empty chain whose element l is at level l.
    EXPECT_THROW(landscape::halochainreport(good, {}), std::invalid_argument);
    std::vector<landscape::LevelWindow> chain = landscape::halochain(good, 4, 2, 2);
    EXPECT_NO_THROW(landscape::halochainreport(good, chain));
    std::swap(chain[1], chain[2]);
    EXPECT_THROW(landscape::halochainreport(good, chain), std::invalid_argument);
}

// Window A of the plan (--pixels 164,48,183,63 on a 256-column world: columns 164 .. 183,
// rows 48 .. 63, halo 16) at factor 4 with the Catmull-Rom reach of 2. Hand-computed:
// at k = 16 (level 2) the halo is 24 parent cells (96 level cells, 6 global cells); the
// parents of level 2's domain (columns 2528 .. 3039, rows 672 .. 1119) are 632 .. 759 and
// 168 .. 279, widened by 2 to 630 .. 761 and 166 .. 281, snapped to whole global cells
// 628 .. 763 and 164 .. 283, so level 1 keeps 136 x 120 and computes 264 x 248 (halo 16
// parent cells, 23 global cells beyond the target), and level 0 keeps and computes the
// 70 x 66 global cells 139 .. 208 by 23 .. 88 (25 beyond the target).
TEST_F(Landscape, RefinementGridHaloChainReportOfWindowA)
{
    const landscape::RegionalWindow target = handtarget(256, 16, 164, 20, 48, 16);

    const std::vector<landscape::LevelWindow> sixteen = landscape::halochain(target, 4, 2, 2);
    for (std::size_t l = 1; l < sixteen.size(); l++)
        EXPECT_EQ(outsideparentkeptbox(sixteen[l], sixteen[l - 1], 2), 0u);
    const landscape::HaloChainReport sixteenreport = landscape::halochainreport(target, sixteen);
    ASSERT_EQ(sixteenreport.levels.size(), 3u);
    expectlevel(sixteenreport.levels[0], 70, 66, 70, 66, 0, 25.0);
    expectlevel(sixteenreport.levels[1], 136, 120, 264, 248, 16, 23.0);
    expectlevel(sixteenreport.levels[2], 320, 256, 512, 448, 24, 6.0);
    EXPECT_EQ(sixteenreport.domaincells, 4620 + 65472 + 229376);
    EXPECT_EQ(sixteenreport.domaincells, 299468);

    // At k = 64 (level 3) the halo is 48 parent cells (192 level cells, 3 global cells);
    // the coarser levels reach 10, 27 and 29 global cells beyond the target (level 2's kept
    // box 2574 .. 2993 by 718 .. 1073 snaps to 2560 .. 3007 by 704 .. 1087).
    const std::vector<landscape::LevelWindow> sixtyfour = landscape::halochain(target, 4, 3, 2);
    for (std::size_t l = 1; l < sixtyfour.size(); l++)
        EXPECT_EQ(outsideparentkeptbox(sixtyfour[l], sixtyfour[l - 1], 2), 0u);
    const landscape::HaloChainReport sixtyfourreport = landscape::halochainreport(target, sixtyfour);
    ASSERT_EQ(sixtyfourreport.levels.size(), 4u);
    expectlevel(sixtyfourreport.levels[0], 78, 74, 78, 74, 0, 29.0);
    expectlevel(sixtyfourreport.levels[1], 168, 152, 296, 280, 16, 27.0);
    expectlevel(sixtyfourreport.levels[2], 448, 384, 640, 576, 24, 10.0);
    expectlevel(sixtyfourreport.levels[3], 1280, 1024, 1664, 1408, 48, 3.0);
    EXPECT_EQ(sixtyfourreport.domaincells, 5772 + 82880 + 368640 + 2342912);
    EXPECT_EQ(sixtyfourreport.domaincells, 2800204);
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
