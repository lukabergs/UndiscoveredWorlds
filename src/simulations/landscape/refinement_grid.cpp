#include "simulations/landscape/refinement_grid.hpp"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <stdexcept>
#include <string>

// The accessors are inline in the header (step S0.4); the window construction, the
// snapping and the halo chain are step S2.1 (docs/architecture/multiscale-plan.md).

namespace landscape
{
namespace
{
// Floor and ceiling of value / divisor for any sign of value (divisor > 0).
std::int64_t floordivide(std::int64_t value, std::int64_t divisor)
{
    const std::int64_t quotient = value / divisor;
    return (value % divisor != 0 && value < 0) ? quotient - 1 : quotient;
}

std::int64_t ceildivide(std::int64_t value, std::int64_t divisor)
{
    const std::int64_t quotient = value / divisor;
    return (value % divisor != 0 && value > 0) ? quotient + 1 : quotient;
}

// The level window with the given kept box (column0 already wrapped, the box already
// inside the rows) and halo in parent cells, its domain by the rule of LevelWindow.
// Throws std::invalid_argument when the domain is wider than the level grid or an
// extent does not fit an int. The caller has checked that the level grid fits below 2^63.
LevelWindow buildwindow(int worldwidth, int factor, int level, std::int64_t halo, std::int64_t column0,
                        std::int64_t columns, std::int64_t row0, std::int64_t rows, const char* caller)
{
    LevelWindow window;
    window.worldwidth = worldwidth;
    window.worldheight = worldwidth / 2;
    window.factor = factor;
    window.level = level;
    const std::int64_t levelwidth = window.levelwidth();
    const std::int64_t levelheight = window.levelheight();
    const std::int64_t cells = level >= 1 ? halo * factor : halo;
    const std::int64_t domaincolumns = columns + 2 * cells;
    if (domaincolumns > levelwidth)
        throw std::invalid_argument(std::string(caller) + ": the domain (box and halo) is wider than the level grid");
    const std::int64_t domainrow0 = std::max<std::int64_t>(0, row0 - cells);
    const std::int64_t domainrows = std::min(levelheight, row0 + rows + cells) - domainrow0;
    if (domaincolumns > INT_MAX || domainrows > INT_MAX)
        throw std::invalid_argument(std::string(caller) + ": the domain's extents do not fit an int");

    window.halo = static_cast<int>(halo);
    window.column0 = column0;
    window.columns = static_cast<int>(columns);
    window.row0 = row0;
    window.rows = static_cast<int>(rows);
    window.domaincolumn0 = wraplevel(column0 - cells, levelwidth);
    window.domaincolumns = static_cast<int>(domaincolumns);
    window.domainrow0 = domainrow0;
    window.domainrows = static_cast<int>(domainrows);
    return window;
}

// levelhalo in 64 bits (levelhalo's int overflows at deep levels before the domain check
// can refuse the window): factor^(L - 1) (4 + ceil((targethalo - 4) / factor^(L - 1))), 4 =
// globalreachcells.
std::int64_t levelhalo64(int targethalo, int factor, int level)
{
    if (level < 1)
        return 0;
    const std::int64_t scale = levelscale(factor, level - 1);
    return scale * (globalreachcells + ceildivide(targethalo - globalreachcells, scale));
}

// The span snapped outwards to multiples of unit (unit >= 1), the rule of snapoutwards.
LevelSpan snapto(LevelSpan span, std::int64_t unit)
{
    LevelSpan snapped;
    snapped.first = floordivide(span.first, unit) * unit;
    if (span.count == 0)
        return snapped;
    const std::int64_t end = ceildivide(span.first + span.count, unit) * unit;
    snapped.count = end - snapped.first;
    return snapped;
}
}

LevelSpan snapoutwards(LevelSpan span, int factor)
{
    if (factor < 2)
        throw std::invalid_argument("snapoutwards: the factor must be at least 2");
    if (span.count < 0)
        throw std::invalid_argument("snapoutwards: the count must not be negative");
    return snapto(span, factor);
}

LevelWindow levelwindow(const RegionalWindow& global, int factor, int level)
{
    if (factor < 2)
        throw std::invalid_argument("levelwindow: the factor must be at least 2");
    if (level < 0)
        throw std::invalid_argument("levelwindow: the level must not be negative");
    const int width = global.worldwidth, height = global.worldheight;
    if (width < 4 || width % 2 != 0 || height != width / 2)
        throw std::invalid_argument("levelwindow: the world must be W x W / 2 with W even and at least 4");
    if (global.column0 < 0 || global.column0 >= width || global.columns < 1 || global.columns > width ||
        global.row0 < 0 || global.rows < 1 || global.row0 + global.rows > height)
        throw std::invalid_argument("levelwindow: the kept box must lie inside the grid");
    if (global.halo < 0)
        throw std::invalid_argument("levelwindow: the halo must not be negative");

    // s_L, then width_L x height_L below 2^63, checked step by step against overflow.
    constexpr std::int64_t largest = INT64_MAX;
    std::int64_t scale = 1;
    for (int step = 0; step < level; step++)
    {
        if (scale > largest / factor)
            throw std::invalid_argument("levelwindow: the level grid does not fit below 2^63 cells");
        scale *= factor;
    }
    if (scale > largest / width)
        throw std::invalid_argument("levelwindow: the level grid does not fit below 2^63 cells");
    const std::int64_t levelwidth = width * scale, levelheight = height * scale;
    if (levelwidth > largest / levelheight)
        throw std::invalid_argument("levelwindow: the level grid does not fit below 2^63 cells");

    return buildwindow(width, factor, level, levelhalo64(global.halo, factor, level), global.column0 * scale,
                       global.columns * scale, global.row0 * scale, global.rows * scale, "levelwindow");
}

std::vector<LevelWindow> halochain(const RegionalWindow& target, int factor, int level, int reach)
{
    if (level < 1)
        throw std::invalid_argument("halochain: the level must be at least 1");
    if (reach < 0)
        throw std::invalid_argument("halochain: the reach must not be negative");
    std::vector<LevelWindow> chain(static_cast<std::size_t>(level) + 1);
    chain[static_cast<std::size_t>(level)] = levelwindow(target, factor, level);
    for (int l = level - 1; l >= 0; l--)
    {
        const LevelWindow& child = chain[static_cast<std::size_t>(l) + 1];
        const std::int64_t scale = levelscale(factor, l);
        const std::int64_t levelwidth = static_cast<std::int64_t>(target.worldwidth) * scale;
        const std::int64_t levelheight = static_cast<std::int64_t>(target.worldheight) * scale;

        // The parent cells of the child's domain, widened by the prolongation's reach.
        LevelSpan columns{child.domaincolumn0 / factor - reach, child.domaincolumns / factor + 2 * static_cast<std::int64_t>(reach)};
        const std::int64_t rowfirst = child.domainrow0 / factor - reach;
        const std::int64_t rowend = rowfirst + child.domainrows / factor + 2 * static_cast<std::int64_t>(reach);
        LevelSpan rows{std::max<std::int64_t>(0, rowfirst), 0};
        rows.count = std::min(levelheight, rowend) - rows.first;

        // Whole global cells (multiples of s_l) at level >= 1; level-0 cells are global cells.
        if (l >= 1)
        {
            columns = snapto(columns, scale);
            rows = snapto(rows, scale);
            const std::int64_t end = std::min(levelheight, rows.first + rows.count);
            rows.first = std::max<std::int64_t>(0, rows.first);
            rows.count = end - rows.first;
        }
        if (columns.count > levelwidth)
            throw std::invalid_argument("halochain: a level's kept box is wider than its grid");
        chain[static_cast<std::size_t>(l)] =
            buildwindow(target.worldwidth, factor, l, levelhalo64(target.halo, factor, l),
                        wraplevel(columns.first, levelwidth), columns.count, rows.first, rows.count, "halochain");
    }
    return chain;
}

RegionalWindow regionalwindowof(const LevelWindow& window)
{
    if (window.level < 1)
        throw std::invalid_argument("regionalwindowof: the level must be at least 1");
    const std::int64_t scale = window.scale();
    if (scale > INT_MAX)
        throw std::invalid_argument("regionalwindowof: factor^level does not fit an int");
    const std::int64_t halocells = window.halocells();
    for (const std::int64_t value : {window.column0, static_cast<std::int64_t>(window.columns), window.row0,
                                     static_cast<std::int64_t>(window.rows), window.domaincolumn0,
                                     static_cast<std::int64_t>(window.domaincolumns), window.domainrow0,
                                     static_cast<std::int64_t>(window.domainrows), halocells})
        if (value % scale != 0)
            throw std::invalid_argument("regionalwindowof: the window is not whole global cells");

    RegionalWindow global;
    global.worldwidth = window.worldwidth;
    global.worldheight = window.worldheight;
    global.scale = static_cast<int>(scale);
    global.halo = static_cast<int>(halocells / scale);
    global.column0 = static_cast<int>(window.column0 / scale);
    global.columns = static_cast<int>(window.columns / scale);
    global.row0 = static_cast<int>(window.row0 / scale);
    global.rows = static_cast<int>(window.rows / scale);
    global.domaincolumn0 = static_cast<int>(window.domaincolumn0 / scale);
    global.domaincolumns = static_cast<int>(window.domaincolumns / scale);
    global.domainrow0 = static_cast<int>(window.domainrow0 / scale);
    global.domainrows = static_cast<int>(window.domainrows / scale);
    return global;
}

HaloChainReport halochainreport(const RegionalWindow& target, const std::vector<LevelWindow>& chain)
{
    if (chain.empty())
        throw std::invalid_argument("halochainreport: the chain is empty");
    HaloChainReport report;
    report.levels.resize(chain.size());
    for (std::size_t l = 0; l < chain.size(); l++)
    {
        const LevelWindow& window = chain[l];
        if (window.level != static_cast<int>(l))
            throw std::invalid_argument("halochainreport: chain element l must be at level l");
        const std::int64_t scale = window.scale();
        const double cells = static_cast<double>(scale);
        // The target's kept box in level cells and the domain's reach beyond it.
        const std::int64_t west = wraplevel(target.column0 * scale - window.domaincolumn0, window.levelwidth());
        const std::int64_t east = window.domaincolumns - west - target.columns * scale;
        const std::int64_t north = target.row0 * scale - window.domainrow0;
        const std::int64_t south = window.domainrows - north - target.rows * scale;

        HaloChainLevel& entry = report.levels[l];
        entry.level = window.level;
        entry.keptcolumns = window.columns;
        entry.keptrows = window.rows;
        entry.domaincolumns = window.domaincolumns;
        entry.domainrows = window.domainrows;
        entry.halo = window.halo;
        entry.west_cells = static_cast<double>(west) / cells;
        entry.east_cells = static_cast<double>(east) / cells;
        entry.north_cells = static_cast<double>(north) / cells;
        entry.south_cells = static_cast<double>(south) / cells;
        report.domaincells += static_cast<std::int64_t>(window.domainsize());
    }
    return report;
}
}
