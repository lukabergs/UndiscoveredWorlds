#pragma once

// Level grids of the refinement hierarchy (docs/architecture/multiscale-plan.md, section 3,
// "Fields and units" and "Kernels"; contracts written in step S0.4, the declared
// functions implemented in S2.1). The world's W x H cell-centred grid (H = W / 2) is
// level 0; level L has s_L = factor^L cells per global cell edge, so its planet-wide grid
// is width_L = W s_L columns by height_L = H s_L rows. Level-L cell (i, j) has the parent
// (i div factor, j div factor) at level L - 1 and the absolute key j width_L + i, which
// exceeds 32 bits at deep levels (at 1 m on a 256-column world width_L is about 67 M):
// absolute positions and keys are 64-bit throughout. Columns wrap modulo width_L, rows
// stop at the poles. A level window is a kept box of level cells, aligned to whole parent
// cells, widened by a halo counted in parent cells (the domain, computed and then
// discarded), as RegionalWindow (regional_grid.hpp) is for a single scale; regional_grid
// stays for the k-based callers until S2.3 retires them. Degrees and metres follow
// regional_grid.hpp's conventions with k = s_L. Planet-free.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
// value modulo period, in 0 .. period - 1 for any sign of value (period > 0).
inline std::int64_t wraplevel(std::int64_t value, std::int64_t period) { return ((value % period) + period) % period; }

// s_L = factor^level, the level's cells per global cell edge (1 at level 0).
inline std::int64_t levelscale(int factor, int level)
{
    std::int64_t scale = 1;
    for (int step = 0; step < level; step++)
        scale *= factor;
    return scale;
}

// The level L with factor^L == scale, or -1 when scale is not a power of factor
// (factor >= 2, scale >= 1).
inline int levelofscale(std::int64_t scale, int factor)
{
    if (factor < 2 || scale < 1)
        return -1;
    int level = 0;
    while (scale % factor == 0)
    {
        scale /= factor;
        level++;
    }
    return scale == 1 ? level : -1;
}

// The global cells around a cell whose pinned channels the channels and the solve of a
// level depend on, at every level: the "+ 2" of the single-scale rule halo >= 2 x
// meanpasses + 2 (regional_pass.hpp). At level L that is 2 factor^(L - 1) parent cells.
constexpr int channelreachglobalcells = 2;

// The halo of level L >= 1 in parent cells for a target whose halo is targethalo global
// cells: factor^(L - 1) (channelreachglobalcells + ceil((targethalo -
// channelreachglobalcells) / factor^(L - 1))), a whole number of global cells at every
// level. The pre-smoothing's reach (targethalo - 2 parent cells at level 1, the
// 2 x meanpasses + 2 rule), rounded up to whole global cells, plus the channels' and the
// solve's 2 global cells. Level 1 takes targethalo (today's single-scale domain); factor 4,
// targethalo 16: level 2 24 parent cells (6 global cells), level 3 48 (3 global cells).
// Level 0: 0 (the global fields are given, not computed).
inline int levelhalo(int targethalo, int factor, int level)
{
    if (level < 1)
        return 0;
    const std::int64_t scale = levelscale(factor, level - 1);
    const std::int64_t excess = targethalo - channelreachglobalcells;
    const std::int64_t globalcells = excess >= 0 ? (excess + scale - 1) / scale : -(-excess / scale);  // ceil
    return static_cast<int>(scale * (channelreachglobalcells + globalcells));
}

// The absolute key of level cell (column, row) on a level grid of levelwidth columns:
// row levelwidth + column, unsigned 64-bit (the analogue of the solve key I + J kW of
// regional_pass.cpp).
inline std::uint64_t levelkey(std::int64_t column, std::int64_t row, std::int64_t levelwidth)
{
    return static_cast<std::uint64_t>(row) * static_cast<std::uint64_t>(levelwidth) + static_cast<std::uint64_t>(column);
}

// A window on level `level`'s planet-wide grid. Absolute positions are int64; the
// extents are int (a window is held in memory, so its extents fit). A valid window has
// worldwidth >= 4 and even, worldheight = worldwidth / 2, factor >= 2, level >= 0,
// halo >= 0; at level >= 1 column0, columns, row0 and rows are multiples of factor (the
// kept box is whole parent cells); 0 <= column0 < levelwidth(); 0 <= row0,
// row0 + rows <= levelheight(); and the domain fields follow from the others by the rule
// below. The accessors assume a valid window and do not check it.
struct LevelWindow
{
    int worldwidth = 0;   // W global columns
    int worldheight = 0;  // H = W / 2 global rows
    int factor = 0;       // level cells per parent cell edge, at least 2 (refinement.factor)
    int level = 0;        // L; 0 is the global grid
    // The halo in parent cells (level L - 1 cells), at level 0 in level-0 cells: the
    // domain extends halocells() level cells beyond the kept box on every side. Windows
    // built by levelwindow and halochain take levelhalo(target halo, factor, level).
    int halo = 0;

    // The kept box in level cells: columns column0 .. column0 + columns - 1, each taken
    // modulo levelwidth() (the box may cross the date line), rows row0 .. row0 + rows - 1.
    std::int64_t column0 = 0;
    int columns = 0;
    std::int64_t row0 = 0;
    int rows = 0;

    // The domain in level cells: the kept box widened by halocells() on every side, rows
    // clipped to 0 .. levelheight() - 1 (domainrow0 = max(0, row0 - halocells()),
    // domainrows = min(levelheight(), row0 + rows + halocells()) - domainrow0); columns
    // never clipped (domaincolumn0 = column0 - halocells() modulo levelwidth(),
    // domaincolumns = columns + 2 halocells(), at most levelwidth()). At level >= 1 the
    // domain is whole parent cells, since levelheight() is a multiple of factor.
    std::int64_t domaincolumn0 = 0;
    int domaincolumns = 0;
    std::int64_t domainrow0 = 0;
    int domainrows = 0;

    std::int64_t scale() const { return levelscale(factor, level); }
    std::int64_t levelwidth() const { return static_cast<std::int64_t>(worldwidth) * scale(); }
    std::int64_t levelheight() const { return static_cast<std::int64_t>(worldheight) * scale(); }
    int halocells() const { return level >= 1 ? halo * factor : halo; }

    std::size_t keptsize() const { return static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows); }
    std::size_t domainsize() const { return static_cast<std::size_t>(domaincolumns) * static_cast<std::size_t>(domainrows); }
    // Row-major index of domain cell (x, y), north row first:
    // 0 <= x < domaincolumns, 0 <= y < domainrows.
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(domaincolumns) + static_cast<std::size_t>(x);
    }

    // The kept box starts at domain cell (keptcolumnoffset(), keptrowoffset()).
    int keptcolumnoffset() const { return static_cast<int>(wraplevel(column0 - domaincolumn0, levelwidth())); }
    int keptrowoffset() const { return static_cast<int>(row0 - domainrow0); }

    // Absolute level indices of domain cell (x, y): the column wraps modulo levelwidth().
    std::int64_t absolutecolumn(int x) const { return wraplevel(domaincolumn0 + x, levelwidth()); }
    std::int64_t absoluterow(int y) const { return domainrow0 + y; }
    std::uint64_t absolutekey(int x, int y) const { return levelkey(absolutecolumn(x), absoluterow(y), levelwidth()); }

    // The domain position of an absolute level cell: the column's offset eastwards from
    // domaincolumn0 modulo levelwidth(), the row's offset from domainrow0. The cell lies in
    // the domain when 0 <= localcolumn < domaincolumns and 0 <= localrow < domainrows.
    int localcolumn(std::int64_t column) const { return static_cast<int>(wraplevel(column - domaincolumn0, levelwidth())); }
    std::int64_t localrow(std::int64_t row) const { return row - domainrow0; }
    bool containscolumn(std::int64_t column) const { return localcolumn(column) < domaincolumns; }
    bool containsrow(std::int64_t row) const { return localrow(row) >= 0 && localrow(row) < domainrows; }

    // The parent (level L - 1) of domain cell (x, y), level >= 1: absolute column
    // absolutecolumn(x) div factor (already in 0 .. levelwidth() / factor - 1) and row
    // absoluterow(y) div factor. Since the domain is whole parent cells, the domain's
    // cells x = factor p .. factor p + factor - 1 share one parent, the domain's p-th
    // parent column (p = x div factor; likewise rows).
    std::int64_t parentcolumn(int x) const { return absolutecolumn(x) / factor; }
    std::int64_t parentrow(int y) const { return absoluterow(y) / factor; }
    // The parent cells the domain spans: domaincolumns / factor by domainrows / factor
    // (level >= 1).
    int domainparentcolumns() const { return domaincolumns / factor; }
    int domainparentrows() const { return domainrows / factor; }

    // The global (level 0) cell holding domain cell (x, y).
    int globalcolumn(int x) const { return static_cast<int>(absolutecolumn(x) / scale()); }
    int globalrow(int y) const { return static_cast<int>(absoluterow(y) / scale()); }
};

// --- Declared in S0.4, implemented in S2.1 (they throw std::logic_error until then). ---

// A half-open interval of level cells first .. first + count - 1 (first may be negative
// or beyond the level width before wrapping).
struct LevelSpan
{
    std::int64_t first = 0;
    std::int64_t count = 0;
};

// The span snapped outwards to whole parent cells: first' = floor(first / factor) factor
// (floor division, also for negative first), end' = ceil((first + count) / factor) factor,
// count' = end' - first'. first' is not wrapped (the caller wraps a column modulo the level
// width, which keeps it aligned since the level width is a multiple of factor at level
// >= 1). count == 0 gives count' == 0 at first'. Throws std::invalid_argument when
// factor < 2 or count < 0.
LevelSpan snapoutwards(LevelSpan span, int factor);

// The level window of the global window's kept box at `level`: column0 =
// global.column0 s_L, columns = global.columns s_L, row0 = global.row0 s_L, rows =
// global.rows s_L (whole global cells are whole parent cells at every level), halo =
// levelhalo(global.halo, factor, level) parent cells (at level 1 global.halo: the domain is
// the global window's domain; at level L >= 2 it reaches channelreachglobalcells +
// ceil((global.halo - channelreachglobalcells) / factor^(L - 1)) global cells; at level 0
// halo 0), the domain by
// the rule of LevelWindow. global.scale is not read. Throws std::invalid_argument when
// factor < 2, level < 0, global is not a window of regionalwindow's rules (worldwidth
// >= 4 and even, worldheight = worldwidth / 2, the kept box inside the grid, halo >= 0),
// the domain is wider than the level grid (domaincolumns > levelwidth()), or
// levelwidth() x levelheight() does not fit below 2^63.
LevelWindow levelwindow(const RegionalWindow& global, int factor, int level);

// The halo chain of a target: the windows of levels 0 .. level a pass at that level
// computes, index = level. Element `level` is levelwindow(target, factor, level). For
// l = level - 1 down to 1, window l is built from window l + 1 (call it c): its kept box
// spans the parent cells of c's domain widened by `reach` level-l cells on every side
// (columns c.domaincolumn0 / factor - reach, count c.domaincolumns / factor + 2 reach;
// rows c.domainrow0 / factor - reach, count c.domainrows / factor + 2 reach, then
// clipped to 0 .. levelheight(l)), snapped outwards to whole global cells (multiples of
// factor^l, as snapoutwards does for one factor), the column wrapped modulo levelwidth(l),
// rows clipped again to the grid; its halo is levelhalo(target.halo, factor, l) parent
// cells and its domain follows by the rule of LevelWindow. So every window of the chain
// has its kept box and its domain in whole global cells (regionalwindowof).
// Window 0 is built from window 1 the same way without the snap (level-0 cells are global
// cells) and with halo 0: the global fields are given, not computed, so level 0's
// domain is its kept box. `reach` is the parent cells the prolongation reads beyond a
// child's parent cell (2 for upsampling.hpp's Catmull-Rom). So every level's domain lies
// inside its parent's kept box, and the parent's halo protects that kept box from the
// parent's own edge. level >= 1; throws std::invalid_argument as levelwindow, when
// level < 1 or reach < 0, and when a level's domain is wider than its grid (whole-planet
// chains are not supported, as regionalwindow refuses a domain wider than the planet).
std::vector<LevelWindow> halochain(const RegionalWindow& target, int factor, int level, int reach);

// The RegionalWindow of a level window whose kept box and domain are whole global cells
// (every window of a halo chain at level >= 1): scale = factor^level, worldwidth and
// worldheight copied, the kept box and the domain in global cells (column0 =
// window.column0 / scale, columns = window.columns / scale, likewise the rows and the
// domain fields) and halo = halocells() / scale global cells, so its fine grid is the level
// window's domain. Throws std::invalid_argument when level < 1, factor^level does not fit
// an int, or a position, an extent or halocells() is not a multiple of factor^level.
RegionalWindow regionalwindowof(const LevelWindow& window);

// The extent of one level of a halo chain, for the run's report.
struct HaloChainLevel
{
    int level = 0;
    int keptcolumns = 0, keptrows = 0;      // level cells
    int domaincolumns = 0, domainrows = 0;  // level cells
    // How far the domain reaches beyond the target's kept box, in global cells (level
    // cells / s_L): west = (target.column0 s_L - domaincolumn0) modulo levelwidth(), east =
    // domaincolumns - west - target.columns s_L, north = target.row0 s_L - domainrow0,
    // south = domainrows - north - target.rows s_L (smaller at a pole, where rows clip).
    // Each level's halo is levelhalo(target.halo, factor, level) parent cells, so the
    // target level's own domain reaches channelreachglobalcells + ceil((target.halo -
    // channelreachglobalcells) / factor^(L - 1)) global cells, and the coarser levels reach
    // further by the prolongation's reach and their own halos.
    int halo = 0;                           // the level's halo in parent cells
    double west_cells = 0.0, east_cells = 0.0, north_cells = 0.0, south_cells = 0.0;
};

struct HaloChainReport
{
    std::vector<HaloChainLevel> levels;  // index = level, 0 .. the target's level
    std::int64_t domaincells = 0;        // sum of domainsize() over the levels
};

// The report of a chain built by halochain for `target`. Throws std::invalid_argument
// when chain is empty or its element l is not at level l.
HaloChainReport halochainreport(const RegionalWindow& target, const std::vector<LevelWindow>& chain);
}
