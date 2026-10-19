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
#include <memory>
#include <vector>

#include "core/mesh/cell_grid.hpp"
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
// level depend on, at every level: the "+ 3" of the single-scale rule halo >= 2 x
// meanpasses + 3 (regional_pass.hpp). At level L that is 3 factor^(L - 1) parent cells.
// linkreachglobalcells: a link's two cells (2 until S3.10, which the bounded solve's 4r covers);
// multiscale plan S3.11's confluences add one: a piece's prefix end reads the owners along its
// path, whose beds read the start surface one global cell beyond (channel_network.hpp, Seams).
constexpr int linkreachglobalcells = 2;
constexpr int channelreachglobalcells = linkreachglobalcells + 1;

// The solve's margin in global cells beyond channelreachglobalcells at every level.
// - With bounded catchments (S3.9, regional.catchment_reach_global_cells r >= 1; analytical_relief.hpp)
//   the solve's reach is bounded by construction, not by trial: a cell's label is a function of the
//   cells within 2r global cells (one flood per global cell within r of it, each confined to r around
//   its own); its catchment lies within r of its outlet's global cell, so within 2r of it, and the
//   catchment's extent reads the labels of those cells: the solve at a cell reads the cells within 4r
//   global cells (the outlet's bed, its link's two cells, lies within that too). Beyond the link's
//   reach that is max(0, 4r - linkreachglobalcells): 2 at r = 1 (the design's "within 2" counted
//   the catchment, not the labels deciding it); S3.11's confluence cell comes on top of it.
// - Without them (r = 0, the solve up to S3.8) the margin of S3.8, found by trial, not derived
//   (unboundedsolvemarginglobalcells): the unbounded solve's sub-basins, the unpinned cells
//   connected between the pinned channels, reach the domain's outer ring at every level. On the
//   synthetic world's P'/Q' pair at k 64 with the valleys on (regional_tributaries_tests.cpp, halo 13)
//   a level-3 domain of 3 global cells beyond the kept box failed (263 kept elevations off, the
//   solve's) and 4 passed; that margin gives 5.
constexpr int unboundedsolvemarginglobalcells = 2;
constexpr int solvemargin(int catchmentreach)
{
    return catchmentreach < 1 ? unboundedsolvemarginglobalcells
                              : (4 * catchmentreach > linkreachglobalcells ? 4 * catchmentreach - linkreachglobalcells : 0);
}
constexpr int solvemarginglobalcells = solvemargin(1);  // the default key's

// The reaches kept in global cells at every level: the channels' and the solve's with its
// margin (here the default key's; levelhalo takes the margin). The others are parent cells at
// every level: the presmoothing's and the smoothing's 2 per pass, the tributaries' tributaryreach,
// the lake flood's (2 global cells at level 1, 1 parent cell + 3 level cells at level >= 2;
// refinement_levels.hpp, lakefloodreach), the forms' and the channel-cell mean's
// (refinement_levels.hpp, formreach and channelmeanreach).
constexpr int globalreachcells = channelreachglobalcells + solvemarginglobalcells;

// The halo of level L >= 1 in parent cells for a target whose halo is targethalo global
// cells, with g = channelreachglobalcells + solvemargin: factor^(L - 1) (g + ceil((targethalo - g) /
// factor^(L - 1))), a whole number of global cells at every level: the reaches in global cells
// (g) kept in global cells, the rest of the target's halo (targethalo - g: the reaches in parent
// cells) as parent cells, rounded up to whole global cells. Level 1 takes targethalo (today's
// single-scale domain); factor 4, targethalo 16, margin 0 (S3.9): level 2 28 parent cells (7
// global cells), level 3 64 (4 global cells); margin 2 (S3.8, the unbounded solve): 32 (8) and 96 (6)
// (S3.11: the channels keep 3 global cells; until S3.10 2: 24, 48, 28, 80).
// Up to S3.8 only the channels' 2 global cells were kept with no margin.
// Level 0: 0 (the global fields are given, not computed).
inline int levelhalo(int targethalo, int factor, int level, int solvemargin = solvemarginglobalcells)
{
    if (level < 1)
        return 0;
    const std::int64_t reach = channelreachglobalcells + solvemargin;
    const std::int64_t scale = levelscale(factor, level - 1);
    const std::int64_t excess = targethalo - reach;
    const std::int64_t globalcells = excess >= 0 ? (excess + scale - 1) / scale : -(-excess / scale);  // ceil
    return static_cast<int>(scale * (reach + globalcells));
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

    // The level-0 ancestor of domain cell (x, y): the id on the world's grid (on lat-lon
    // J W + I) of the global cell holding it.
    std::size_t globalcell(int x, int y) const
    {
        return static_cast<std::size_t>(absoluterow(y) / scale()) * static_cast<std::size_t>(worldwidth) +
               static_cast<std::size_t>(absolutecolumn(x) / scale());
    }
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
// levelhalo(global.halo, factor, level, solvemargin) parent cells (at level 1 global.halo: the
// domain is the global window's domain; at level L >= 2 it reaches g + ceil((global.halo - g) /
// factor^(L - 1)) global cells, g = channelreachglobalcells + solvemargin; at level 0
// halo 0), the domain by
// the rule of LevelWindow. global.scale is not read. Throws std::invalid_argument when
// factor < 2, level < 0, solvemargin < 0, global is not a window of regionalwindow's rules (worldwidth
// >= 4 and even, worldheight = worldwidth / 2, the kept box inside the grid, halo >= 0),
// the domain is wider than the level grid (domaincolumns > levelwidth()), or
// levelwidth() x levelheight() does not fit below 2^63.
LevelWindow levelwindow(const RegionalWindow& global, int factor, int level, int solvemargin = solvemarginglobalcells);

// The halo chain of a target: the windows of levels 0 .. level a pass at that level
// computes, index = level. Element `level` is levelwindow(target, factor, level, solvemargin). For
// l = level - 1 down to 1, window l is built from window l + 1 (call it c): its kept box
// spans the parent cells of c's domain widened by `reach` level-l cells on every side
// (columns c.domaincolumn0 / factor - reach, count c.domaincolumns / factor + 2 reach;
// rows c.domainrow0 / factor - reach, count c.domainrows / factor + 2 reach, then
// clipped to 0 .. levelheight(l)), snapped outwards to whole global cells (multiples of
// factor^l, as snapoutwards does for one factor), the column wrapped modulo levelwidth(l),
// rows clipped again to the grid; its halo is levelhalo(target.halo, factor, l, solvemargin) parent
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
std::vector<LevelWindow> halochain(const RegionalWindow& target, int factor, int level, int reach,
                                   int solvemargin = solvemarginglobalcells);

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
    // target level's own domain reaches globalreachcells + ceil((target.halo -
    // globalreachcells) / factor^(L - 1)) global cells, and the coarser levels reach
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

// --- The chain on key sets (docs/architecture/sphere-mesh.md section 8, "Halo by rings";
// T8 W2-W1b): the same windows on any grid kind, built on the level-0 grid's table. ---

// levelhalo in rings of level-0 cells: levelhalo / factor^(L - 1), a whole number by
// construction (targethalo at level 1); 0 at level 0.
int levelhalorings(int targethalo, int factor, int level, int solvemargin = solvemarginglobalcells);

// A level window whose kept box and domain are unions of whole level-0 cells (every
// window of a chain), held as the ascending keys of those level-0 cells; the level's
// own keys are their descendants (on the cube the per-face quadtree, on lat-lon the
// s_L x s_L blocks of the level grid, key row W s_L + column).
struct KeyLevelWindow
{
    mesh::GridKind kind = mesh::GridKind::cube;
    int basecells = 0;  // level 0: the cube's N or the lat-lon W
    int factor = 0;
    int level = 0;
    int halo = 0;       // levelhalo, parent cells
    int halorings = 0;  // the same in level-0 rings (8-connected hops on the level-0 table)
    std::vector<mesh::LevelKey> keptglobal;    // ascending level-0 keys
    std::vector<mesh::LevelKey> domainglobal;  // ascending level-0 keys, a superset of keptglobal
    // Lat-lon: the level-0 column where the kept box and the domain start eastwards
    // (halochain's column0 and domaincolumn0 / s_L), also for a set over every column,
    // whose start only the chain knows (keyhalochain; a full-width target box starts at 0).
    std::int64_t keptwest = 0;
    std::int64_t domainwest = 0;

    std::int64_t scale() const { return levelscale(factor, level); }
    // The level's cells per face edge (cube) or columns (lat-lon).
    std::int64_t levelcells() const { return static_cast<std::int64_t>(basecells) * scale(); }
    std::size_t keptsize() const { return keptglobal.size() * static_cast<std::size_t>(scale() * scale()); }
    std::size_t domainsize() const { return domainglobal.size() * static_cast<std::size_t>(scale() * scale()); }
    // The level keys of the kept box and of the domain, ascending.
    std::vector<mesh::LevelKey> keptkeys() const;
    std::vector<mesh::LevelKey> domainkeys() const;
};

// The level-0 cells of a full grid whose centres lie in a box of the lat-lon export's
// pixels (exportcolumns x exportcolumns / 2; columns x0 .. x1 modulo the export width,
// x0 > x1 crossing the date line; rows y0 .. y1 from the north): longitudes from
// -180 + 360 x0 / W eastwards over the box's columns, latitudes from 90 - 180 (y1 + 1) / H
// (included) to 90 - 180 y0 / H (excluded). On a lat-lon grid of the export's width these
// are the box's own cells. Ascending keys (cell ids). Throws std::invalid_argument for a
// subgrid, an odd or too small export width or rows outside 0 .. H - 1 or y0 > y1.
std::vector<mesh::LevelKey> cellsinpixelbox(const mesh::CellGrid& base, int exportcolumns, int x0, int y0, int x1,
                                            int y1);

// The keys within `rings` rings of a key set on a full grid's table (ascending).
std::vector<mesh::LevelKey> dilatecells(const mesh::CellGrid& base, const std::vector<mesh::LevelKey>& keys, int rings);

// The halo chain on keys: windows of levels 0 .. level (index = level) for a target whose
// kept box is the level-0 cells kept0 (on the cube the cells whose centres lie in the
// --pixels box, cellsinpixelbox; the target's kept level cells are their descendants)
// and whose halo is targethalo global cells, on base, the world's full level-0 grid:
// - the target level's domain: kept0 dilated by levelhalorings(targethalo, ..., level) rings;
// - level l = level - 1 .. 1: the kept box is the child's domain dilated by
//   ceil(reach / factor^l) rings (reach level-l cells snapped outwards to whole level-0
//   cells: 1 ring while reach < factor^l), the domain that dilated by levelhalorings(l);
// - level 0: the child's domain dilated by reach rings; its domain is its kept box.
// On lat-lon this is halochain's chain cell for cell (rings of the 8-connected table are
// the rectangles widened by whole cells, columns wrapping, rows stopping at the poles).
// level >= 1; throws std::invalid_argument for a subgrid base, level < 1, factor < 2,
// reach < 0, targethalo < 0, solvemargin < 0, a key out of range, or kept0 empty.
std::vector<KeyLevelWindow> keyhalochain(const mesh::CellGrid& base, int factor, int level,
                                         const std::vector<mesh::LevelKey>& kept0, int targethalo, int reach,
                                         int solvemargin = solvemarginglobalcells);

// The cell grid of a key window's domain (regionalwindowof's analogue): on the cube
// mesh::cubewindowgrid over the domain's level keys (compact or stored); on lat-lon the
// domain must be a box of level-0 cells (as every lat-lon chain window is; columns as
// the wrapped interval outside the widest gap, rows contiguous) and the grid is
// mesh::latlonwindowgrid of that box at the level (compact not read). Throws
// std::invalid_argument for a lat-lon domain that is not a box.
std::shared_ptr<const mesh::CellGrid> keywindowgrid(const KeyLevelWindow& window, double radius_m, bool compact = true);
}
