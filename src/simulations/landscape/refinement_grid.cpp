#include "simulations/landscape/refinement_grid.hpp"

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
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
// can refuse the window): factor^(L - 1) (g + ceil((targethalo - g) / factor^(L - 1))), g =
// channelreachglobalcells + solvemargin.
std::int64_t levelhalo64(int targethalo, int factor, int level, int solvemargin)
{
    if (level < 1)
        return 0;
    const std::int64_t reach = channelreachglobalcells + solvemargin;
    const std::int64_t scale = levelscale(factor, level - 1);
    return scale * (reach + ceildivide(targethalo - reach, scale));
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

LevelWindow levelwindow(const RegionalWindow& global, int factor, int level, int solvemargin)
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
    if (solvemargin < 0)
        throw std::invalid_argument("levelwindow: the solve's margin must not be negative");

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

    return buildwindow(width, factor, level, levelhalo64(global.halo, factor, level, solvemargin), global.column0 * scale,
                       global.columns * scale, global.row0 * scale, global.rows * scale, "levelwindow");
}

std::vector<LevelWindow> halochain(const RegionalWindow& target, int factor, int level, int reach, int solvemargin)
{
    if (level < 1)
        throw std::invalid_argument("halochain: the level must be at least 1");
    if (reach < 0)
        throw std::invalid_argument("halochain: the reach must not be negative");
    std::vector<LevelWindow> chain(static_cast<std::size_t>(level) + 1);
    chain[static_cast<std::size_t>(level)] = levelwindow(target, factor, level, solvemargin);
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
            buildwindow(target.worldwidth, factor, l, levelhalo64(target.halo, factor, l, solvemargin),
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

int levelhalorings(int targethalo, int factor, int level, int solvemargin)
{
    if (level < 1)
        return 0;
    return static_cast<int>(levelhalo(targethalo, factor, level, solvemargin) / levelscale(factor, level - 1));
}

namespace
{
std::vector<mesh::LevelKey> descendants(const KeyLevelWindow& window, const std::vector<mesh::LevelKey>& global)
{
    const std::int64_t s = window.scale(), n = window.basecells, fine = window.levelcells();
    std::vector<mesh::LevelKey> keys;
    keys.reserve(global.size() * static_cast<std::size_t>(s * s));
    for (mesh::LevelKey key : global)
    {
        if (window.kind == mesh::GridKind::cube)
        {
            const mesh::CubeCell cell = mesh::cubecell(n, key);
            for (std::int64_t dj = 0; dj < s; dj++)
                for (std::int64_t di = 0; di < s; di++)
                    keys.push_back(mesh::cubekey(fine, cell.face, cell.column * s + di, cell.row * s + dj));
        }
        else
        {
            const std::int64_t column = static_cast<std::int64_t>(key % static_cast<mesh::LevelKey>(n));
            const std::int64_t row = static_cast<std::int64_t>(key / static_cast<mesh::LevelKey>(n));
            for (std::int64_t dj = 0; dj < s; dj++)
                for (std::int64_t di = 0; di < s; di++)
                    keys.push_back(levelkey(column * s + di, row * s + dj, fine));
        }
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}
}

std::vector<mesh::LevelKey> KeyLevelWindow::keptkeys() const { return descendants(*this, keptglobal); }
std::vector<mesh::LevelKey> KeyLevelWindow::domainkeys() const { return descendants(*this, domainglobal); }

std::vector<mesh::LevelKey> cellsinpixelbox(const mesh::CellGrid& base, int exportcolumns, int x0, int y0, int x1, int y1)
{
    const int height = exportcolumns / 2;
    if (base.issubgrid())
        throw std::invalid_argument("cellsinpixelbox needs a full grid");
    if (exportcolumns < 4 || exportcolumns % 2 != 0 || y0 < 0 || y1 >= height || y0 > y1)
        throw std::invalid_argument("cellsinpixelbox: an even export width of at least 4 and rows 0 <= y0 <= y1 < W / 2");
    constexpr double pi = 3.14159265358979323846;
    const double west = -pi + 2.0 * pi * static_cast<double>(wraplevel(x0, exportcolumns)) / exportcolumns;
    const double width = 2.0 * pi * static_cast<double>(wraplevel(static_cast<std::int64_t>(x1) - x0, exportcolumns) + 1) / exportcolumns;
    const double north = 0.5 * pi - pi * y0 / height, south = 0.5 * pi - pi * (y1 + 1) / height;
    std::vector<mesh::LevelKey> keys;
    for (mesh::CellId c = 0; c < base.cellcount(); c++)
    {
        const double latitude = base.latitude_rad(c);
        if (latitude < south || latitude >= north)
            continue;
        double offset = std::fmod(base.longitude_rad(c) - west, 2.0 * pi);
        if (offset < 0.0)
            offset += 2.0 * pi;
        if (offset < width)
            keys.push_back(c);
    }
    return keys;
}

std::vector<mesh::LevelKey> dilatecells(const mesh::CellGrid& base, const std::vector<mesh::LevelKey>& keys, int rings)
{
    std::vector<mesh::CellId> seeds;
    seeds.reserve(keys.size());
    for (mesh::LevelKey key : keys)
        seeds.push_back(static_cast<mesh::CellId>(key));
    const std::vector<int> ring = mesh::ringdistances(base, seeds, rings);
    std::vector<mesh::LevelKey> result;
    for (mesh::CellId c = 0; c < base.cellcount(); c++)
        if (ring[c] >= 0)
            result.push_back(c);
    return result;
}

std::vector<KeyLevelWindow> keyhalochain(const mesh::CellGrid& base, int factor, int level, const std::vector<mesh::LevelKey>& kept0,
                                         int targethalo, int reach, int solvemargin)
{
    if (base.issubgrid() || level < 1 || factor < 2 || reach < 0 || targethalo < 0 || solvemargin < 0 || kept0.empty())
        throw std::invalid_argument("keyhalochain needs a full base grid, level >= 1, factor >= 2, non-negative reach, halo and "
                                    "margin, and a kept box");
    std::vector<mesh::LevelKey> kept = kept0;
    std::sort(kept.begin(), kept.end());
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    if (kept.back() >= base.cellcount())
        throw std::invalid_argument("keyhalochain: a kept key is out of range");
    std::vector<KeyLevelWindow> chain(static_cast<std::size_t>(level) + 1);
    for (int l = level; l >= 0; l--)
    {
        KeyLevelWindow& window = chain[static_cast<std::size_t>(l)];
        window.kind = base.kind();
        window.basecells = base.kind() == mesh::GridKind::cube ? base.facecells() : base.columns();
        window.factor = factor;
        window.level = l;
        window.halo = levelhalo(targethalo, factor, l, solvemargin);
        window.halorings = levelhalorings(targethalo, factor, l, solvemargin);
        if (l == level)
            window.keptglobal = kept;
        else
        {
            const std::vector<mesh::LevelKey>& child = chain[static_cast<std::size_t>(l) + 1].domainglobal;
            const std::int64_t scale = levelscale(factor, l);
            const int rings = l == 0 ? reach : static_cast<int>((reach + scale - 1) / scale);
            window.keptglobal = dilatecells(base, child, rings);
        }
        window.domainglobal = l == 0 ? window.keptglobal : dilatecells(base, window.keptglobal, window.halorings);
    }
    return chain;
}

std::shared_ptr<const mesh::CellGrid> keywindowgrid(const KeyLevelWindow& window, double radius_m, bool compact)
{
    if (window.kind == mesh::GridKind::cube)
        return mesh::cubewindowgrid(static_cast<int>(window.levelcells()), radius_m, window.domainkeys(), compact);
    const int width = window.basecells;
    std::vector<int> columns, rows;
    for (mesh::LevelKey key : window.domainglobal)
    {
        columns.push_back(static_cast<int>(key % static_cast<mesh::LevelKey>(width)));
        rows.push_back(static_cast<int>(key / static_cast<mesh::LevelKey>(width)));
    }
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    if (columns.empty())
        throw std::invalid_argument("keywindowgrid: an empty domain");
    // The columns as the wrapped interval outside the widest gap.
    std::size_t start = 0;
    int widest = columns.front() + width - columns.back() - 1;
    for (std::size_t k = 1; k < columns.size(); k++)
        if (columns[k] - columns[k - 1] - 1 > widest)
        {
            widest = columns[k] - columns[k - 1] - 1;
            start = k;
        }
    const std::size_t count = columns.size();
    if (static_cast<int>(count) + widest != width || rows.back() - rows.front() + 1 != static_cast<int>(rows.size()) ||
        window.domainglobal.size() != count * rows.size())
        throw std::invalid_argument("keywindowgrid: a lat-lon domain that is not a box");
    const std::int64_t s = window.scale();
    mesh::LatLonWindowBox box;
    box.column0 = static_cast<std::int64_t>(columns[start]) * s;
    box.columns = static_cast<int>(static_cast<std::int64_t>(count) * s);
    box.row0 = static_cast<int>(rows.front() * s);
    box.rows = static_cast<int>(static_cast<std::int64_t>(rows.size()) * s);
    return mesh::latlonwindowgrid(static_cast<int>(window.levelcells()), radius_m, box);
}
}
