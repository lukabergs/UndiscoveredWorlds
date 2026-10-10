#include "simulations/landscape/regional_tributaries.hpp"

// Tributaries per level (multiscale plan S3.3; the rules in the header). Every position below is
// an absolute level cell (column unwrapped, wrapped modulo the level's width when looked up or
// hashed) or an offset from one, so a tributary is drawn alike in every window that holds its
// neighbourhood.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double secondsperyear = 365.25 * 86400.0;
// Hash salts.
constexpr std::int64_t candidatesalt = 0x54524942;  // a candidate: angle, length, priority
constexpr std::int64_t nodesalt = 0x4e4f4445;       // an inherited parent cell's node
constexpr std::int64_t slotsalt = 0x534c4f54;       // a slot's opening
constexpr std::size_t nocell = static_cast<std::size_t>(-1);

enum CellKind : uint8_t
{
    freecell = 0,
    linkcell = 1,
    inheritedcell = 2,
    growncell = 3,
};

double unitinterval(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }
double signedunit(std::uint64_t hash) { return 2.0 * unitinterval(hash) - 1.0; }

// An absolute level cell; the column may lie outside 0 .. planetwidth - 1 (unwrapped).
struct Cell
{
    std::int64_t x = 0;
    std::int64_t y = 0;
};

bool samecell(const Cell& a, const Cell& b) { return a.x == b.x && a.y == b.y; }
bool neighbours8(const Cell& a, const Cell& b)
{
    return !samecell(a, b) && std::llabs(a.x - b.x) <= 1 && std::llabs(a.y - b.y) <= 1;
}

// The level's domain grid in absolute cells.
struct Grid
{
    std::int64_t x0 = 0;           // absolute column of domain cell (0, 0), 0 .. planetwidth - 1
    std::int64_t y0 = 0;
    int width = 0;
    int height = 0;
    std::int64_t planetwidth = 0;  // the level's width in cells

    bool index(std::int64_t x, std::int64_t y, std::size_t& out) const
    {
        const std::int64_t row = y - y0;
        if (row < 0 || row >= height)
            return false;
        const std::int64_t column = wraplevel(x - x0, planetwidth);
        if (column >= width)
            return false;
        out = static_cast<std::size_t>(row) * static_cast<std::size_t>(width) + static_cast<std::size_t>(column);
        return true;
    }
    Cell cell(std::size_t index) const
    {
        return { x0 + static_cast<std::int64_t>(index % static_cast<std::size_t>(width)),
                 y0 + static_cast<std::int64_t>(index / static_cast<std::size_t>(width)) };
    }
    std::uint64_t key(const Cell& cell) const { return levelkey(wraplevel(cell.x, planetwidth), cell.y, planetwidth); }
    // The offset b - a in columns, wrapped to the shorter way round the planet.
    std::int64_t dx(std::int64_t a, std::int64_t b) const
    {
        std::int64_t d = wraplevel(b - a, planetwidth);
        if (d > planetwidth / 2)
            d -= planetwidth;
        return d;
    }
};

Grid gridof(const RegionalWindow& window)
{
    Grid grid;
    grid.planetwidth = static_cast<std::int64_t>(window.worldwidth) * window.scale;
    grid.x0 = static_cast<std::int64_t>(window.domaincolumn0) * window.scale;
    grid.y0 = static_cast<std::int64_t>(window.domainrow0) * window.scale;
    grid.width = window.domainfinewidth();
    grid.height = window.domainfineheight();
    return grid;
}

int directioncode(std::int64_t ox, std::int64_t oy)
{
    for (int direction = 1; direction <= directioncount; direction++)
        if (directionoffsetx(direction) == ox && directionoffsety(direction) == oy)
            return direction;
    return 0;
}

// The domain cell a channel cell drains into (nocell without a direction or outside the domain).
std::size_t downstream(const Grid& grid, const FineChannels& channels, std::size_t index)
{
    const int direction = channels.direction[index];
    if (direction < 1 || direction > directioncount)
        return nocell;
    const Cell cell = grid.cell(index);
    std::size_t next = nocell;
    if (!grid.index(cell.x + directionoffsetx(direction), cell.y + directionoffsety(direction), next))
        return nocell;
    return next;
}

// The 8-connected chain of cells from `from` along the offsets of `points` (cell-unit offsets
// from from's centre), as legpath builds its chains: a cell already on the chain cuts it back to
// that cell, and a cell that neighbours the cell before the last replaces the last (an L corner
// becomes a diagonal step).
struct Chain
{
    std::vector<Cell> cells;
    void add(const Cell& cell)
    {
        if (samecell(cell, cells.back()))
            return;
        const auto found = std::find_if(cells.begin(), cells.end(), [&](const Cell& other) { return samecell(other, cell); });
        if (found != cells.end())
        {
            cells.erase(found + 1, cells.end());
            return;
        }
        while (cells.size() >= 2 && neighbours8(cells[cells.size() - 2], cell))
            cells.pop_back();
        cells.push_back(cell);
    }
};

// The straight line from cell a to the cell at offset (ox, oy) from it, both ends included.
std::vector<Cell> linecells(const Cell& a, std::int64_t ox, std::int64_t oy)
{
    Chain chain;
    chain.cells.push_back(a);
    const std::int64_t span = std::max(std::llabs(ox), std::llabs(oy));
    const int samples = static_cast<int>(4 * span) + 1;
    for (int s = 1; s <= samples; s++)
    {
        if (s == samples)
        {
            chain.add({ a.x + ox, a.y + oy });
            continue;
        }
        const double t = static_cast<double>(s) / samples;
        chain.add({ a.x + static_cast<std::int64_t>(std::floor(0.5 + t * static_cast<double>(ox))),
                    a.y + static_cast<std::int64_t>(std::floor(0.5 + t * static_cast<double>(oy))) });
    }
    return chain.cells;
}

// The ray from a's centre along the unit vector (ux, uy), `length` cells long.
std::vector<Cell> raycells(const Cell& a, double ux, double uy, double length)
{
    Chain chain;
    chain.cells.push_back(a);
    const int samples = static_cast<int>(std::ceil(4.0 * length));
    for (int s = 1; s <= samples; s++)
    {
        const double t = length * static_cast<double>(s) / samples;
        chain.add({ a.x + static_cast<std::int64_t>(std::floor(0.5 + t * ux)),
                    a.y + static_cast<std::int64_t>(std::floor(0.5 + t * uy)) });
    }
    return chain.cells;
}

bool istributary(const FineChannels& channels, const std::vector<uint8_t>& sea, std::size_t cell)
{
    return channels.channel[cell] != 0 && channels.globalcell[cell] < 0 && sea[cell] == 0;
}

std::vector<uint8_t> orders(const Grid& grid, const FineChannels& channels, const std::vector<uint8_t>& sea)
{
    const std::size_t count = channels.channel.size();
    std::vector<uint8_t> order(count, 0);
    std::vector<int32_t> pending(count, 0);
    std::vector<uint8_t> highest(count, 0), times(count, 0);
    std::vector<std::size_t> next(count, nocell);
    std::vector<std::size_t> queue;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (!istributary(channels, sea, cell))
            continue;
        const std::size_t down = downstream(grid, channels, cell);
        if (down != nocell && istributary(channels, sea, down))
        {
            next[cell] = down;
            pending[down]++;
        }
    }
    for (std::size_t cell = 0; cell < count; cell++)
        if (istributary(channels, sea, cell) && pending[cell] == 0)
            queue.push_back(cell);
    for (std::size_t head = 0; head < queue.size(); head++)
    {
        const std::size_t cell = queue[head];
        const int value = highest[cell] == 0 ? 1 : (times[cell] >= 2 ? highest[cell] + 1 : highest[cell]);
        order[cell] = static_cast<uint8_t>(std::min(value, 255));
        const std::size_t down = next[cell];
        if (down == nocell)
            continue;
        if (order[cell] > highest[down])
        {
            highest[down] = order[cell];
            times[down] = 1;
        }
        else if (order[cell] == highest[down] && times[down] < 2)
            times[down]++;
        if (--pending[down] == 0)
            queue.push_back(down);
    }
    return order;
}

// A grown candidate.
struct Candidate
{
    std::size_t junction = nocell;
    std::vector<std::size_t> path;  // domain cells, junction excluded, from the junction outwards
    std::vector<double> along;      // their path length from the junction, cells
    std::uint64_t priority = 0;
    std::uint64_t slot = 0;
};

// The two highest priorities covering a cell, with their candidates.
struct Cover
{
    std::uint64_t first = 0, second = 0;
    int firstid = -1, secondid = -1;
    void add(std::uint64_t priority, int id)
    {
        if (id == firstid || id == secondid)
            return;
        if (firstid < 0 || priority > first)
        {
            second = first;
            secondid = firstid;
            first = priority;
            firstid = id;
        }
        else if (secondid < 0 || priority > second)
        {
            second = priority;
            secondid = id;
        }
    }
};
}

std::vector<uint8_t> tributaryorders(const RegionalWindow& window, const FineChannels& channels, const std::vector<uint8_t>& sea)
{
    const Grid grid = gridof(window);
    const std::size_t count = static_cast<std::size_t>(grid.width) * static_cast<std::size_t>(grid.height);
    if (channels.channel.size() != count || channels.direction.size() != count || channels.globalcell.size() != count ||
        sea.size() != count)
        throw std::invalid_argument("tributaryorders: the channels and the sea must hold the window's domain cells");
    return orders(grid, channels, sea);
}

void tributarieslevel(const LevelContext& context, LevelState& state) { growtributaries(context, state); }

TributaryReport growtributaries(const LevelContext& context, LevelState& state)
{
    TributaryReport report;
    state.tributaries = 0;
    state.tributarycells = 0;
    state.tributarybedsabove = 0;
    const TributaryParameters& parameters = context.parameters.tributaries;
    if (parameters.enabled == 0)
        return report;
    if (!(parameters.spacing_km > 0.0) || !std::isfinite(parameters.spacing_km) || !(parameters.hackcoefficient > 0.0) ||
        !std::isfinite(parameters.hackcoefficient) || !(parameters.hackexponent > 0.0) || !(parameters.hackexponent < 1.0) ||
        !(parameters.junctionangle_deg > 0.0) || !(parameters.junctionangle_deg < 180.0) ||
        !(parameters.anglespread >= 0.0) || !(parameters.anglespread < 1.0) || !(parameters.nodemargin >= 0.0) ||
        !(parameters.nodemargin < 0.5) || !(parameters.lengthminfraction >= 0.0) || !(parameters.lengthminfraction <= 1.0))
        throw std::invalid_argument("tributarieslevel: spacing_km and hackcoefficient must be positive and finite, "
                                    "hackexponent in (0, 1), junctionangle_deg in (0, 180), anglespread in [0, 1), "
                                    "nodemargin in [0, 0.5) and lengthminfraction in [0, 1]");

    const RegionalWindow& window = context.window;
    const RegionalGlobalInputs& inputs = context.inputs;
    const int k = window.scale;
    const int edge = context.factor >= 2 ? context.factor : k;
    if (edge < 3)
        return report;  // no length lies in [2, edge)
    const Grid grid = gridof(window);
    const std::size_t count = static_cast<std::size_t>(grid.width) * static_cast<std::size_t>(grid.height);
    FineChannels& channels = state.channels;
    const std::vector<uint8_t>& sea = state.sea;
    const std::vector<float>& surface = context.incision.active ? state.bedsurface_m : state.surface_m;
    if (channels.channel.size() != count || sea.size() != count || surface.size() != count ||
        state.surface_m.size() != count)
        throw std::invalid_argument("tributarieslevel: the channels, the sea and the surface must hold the domain's cells");
    const std::uint64_t seed = context.parameters.channels.seed;
    const int level = context.level;
    const double cellsize_km = pi * inputs.radius_m / (static_cast<double>(window.worldheight) * k) / 1000.0;
    const std::int64_t parentwidth = grid.planetwidth / edge;  // the level's parent cells around the planet

    std::vector<uint8_t> kind(count, freecell);
    for (std::size_t cell = 0; cell < count; cell++)
        if (channels.channel[cell] != 0)
            kind[cell] = linkcell;

    const auto island = [&](std::size_t cell) { return sea[cell] == 0; };
    // The land link 8-neighbour of `cell` with the lowest bed, then the smallest key; nocell if none.
    const auto linkneighbour = [&](const Cell& cell) {
        std::size_t best = nocell;
        std::uint64_t bestkey = 0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            std::size_t other = nocell;
            const Cell neighbour{ cell.x + directionoffsetx(direction), cell.y + directionoffsety(direction) };
            if (!grid.index(neighbour.x, neighbour.y, other) || kind[other] != linkcell || !island(other))
                continue;
            const std::uint64_t key = grid.key(neighbour);
            if (best == nocell || channels.bed_m[other] < channels.bed_m[best] ||
                (channels.bed_m[other] == channels.bed_m[best] && key < bestkey))
            {
                best = other;
                bestkey = key;
            }
        }
        return best;
    };
    const auto stepdirection = [&](std::size_t from, std::size_t to) {
        const Cell a = grid.cell(from), b = grid.cell(to);
        return static_cast<uint8_t>(directioncode(grid.dx(a.x, b.x), b.y - a.y));
    };
    // The parent cell (column wrapped) of a level cell, and whether two parent cells are within
    // one of each other (the parent cell and its 8 neighbours).
    const auto parentof = [&](const Cell& cell) {
        return Cell{ wraplevel(cell.x, grid.planetwidth) / edge, cell.y / edge };
    };
    const auto nearparent = [&](const Cell& a, const Cell& b) {
        std::int64_t dx = wraplevel(a.x - b.x, parentwidth);
        if (dx > parentwidth / 2)
            dx -= parentwidth;
        return std::llabs(dx) <= 1 && std::llabs(a.y - b.y) <= 1;
    };

    // 1. Inherited (header, item 1).
    if (context.parent != nullptr)
    {
        const ParentFields& parent = *context.parent;
        const LevelWindow& pw = parent.window;
        const FineChannels& pch = parent.channels;
        const std::int64_t plw = pw.levelwidth();
        const int f = context.factor;
        if (pch.channel.size() != pw.domainsize() || parent.sea.size() != pw.domainsize() ||
            pch.bed_m.size() != pw.domainsize())
            throw std::invalid_argument("tributarieslevel: the parent's channels and sea must hold its domain's cells");
        // A parent cell in the parent's kept box: its domain index.
        const auto parentindex = [&](std::int64_t column, std::int64_t row, std::size_t& out) {
            column = wraplevel(column, plw);
            if (wraplevel(column - pw.column0, plw) >= pw.columns || row < pw.row0 || row >= pw.row0 + pw.rows)
                return false;
            out = pw.index(pw.localcolumn(column), static_cast<int>(pw.localrow(row)));
            return true;
        };
        const auto parenttributary = [&](std::size_t p) {
            return pch.channel[p] != 0 && pch.globalcell[p] < 0 && parent.sea[p] == 0;
        };
        struct ParentCell
        {
            std::int64_t column = 0, row = 0;
            std::size_t index = 0;
            std::size_t down = nocell;  // the downstream parent cell (domain index, kept box, a channel), nocell if none
            int depth = 0;              // tributary cells down to the link: 0 unknown, -1 none (not drawn)
            float linkbed = 0.0f;       // the parent's bed of that link cell
        };
        std::vector<ParentCell> cells;
        std::unordered_map<std::size_t, int> lookup;
        for (int y = 0; y < pw.rows; y++)
            for (int x = 0; x < pw.columns; x++)
            {
                const std::int64_t column = wraplevel(pw.column0 + x, plw), row = pw.row0 + y;
                std::size_t p = 0;
                if (!parentindex(column, row, p) || !parenttributary(p))
                    continue;
                ParentCell cell;
                cell.column = column;
                cell.row = row;
                cell.index = p;
                const int direction = pch.direction[p];
                std::size_t q = 0;
                if (direction >= 1 && direction <= directioncount &&
                    parentindex(column + directionoffsetx(direction), row + directionoffsety(direction), q) &&
                    pch.channel[q] != 0)
                    cell.down = q;
                lookup.emplace(p, static_cast<int>(cells.size()));
                cells.push_back(cell);
            }
        report.inheritedparents = static_cast<int>(cells.size());
        // The parent's chains (in its kept box, which every window shares where it matters): the
        // depth and the bed of the link each tributary cell reaches.
        for (std::size_t c = 0; c < cells.size(); c++)
        {
            std::vector<std::size_t> pending;
            std::size_t current = c;
            int base = -1;
            float linkbed = 0.0f;
            while (true)
            {
                if (cells[current].depth != 0)
                {
                    base = cells[current].depth;
                    linkbed = cells[current].linkbed;
                    break;
                }
                pending.push_back(current);
                if (pending.size() > cells.size())
                {
                    base = -1;
                    break;
                }
                const std::size_t q = cells[current].down;
                if (q == nocell)
                {
                    base = -1;
                    break;
                }
                const auto found = lookup.find(q);
                if (found == lookup.end())
                {
                    // A link of the parent (land: a tributary never drains into a sea channel cell).
                    base = parent.sea[q] == 0 && std::isfinite(pch.bed_m[q]) ? 0 : -1;
                    linkbed = pch.bed_m[q];
                    break;
                }
                current = static_cast<std::size_t>(found->second);
            }
            for (auto it = pending.rbegin(); it != pending.rend(); ++it)
            {
                base = base < 0 ? -1 : base + 1;
                cells[*it].depth = base;
                cells[*it].linkbed = linkbed;
            }
        }
        const int margin = static_cast<int>(std::floor(parameters.nodemargin * f));
        const int span = std::max(1, f - 2 * margin);
        const auto nodeof = [&](const ParentCell& cell) {
            const std::uint64_t hash =
                deterministiccombine(seed, nodesalt, level, static_cast<std::int64_t>(levelkey(cell.column, cell.row, plw)));
            return Cell{ cell.column * f + margin + static_cast<std::int64_t>(hash % static_cast<std::uint64_t>(span)),
                         cell.row * f + margin + static_cast<std::int64_t>((hash >> 32) % static_cast<std::uint64_t>(span)) };
        };
        // Every segment on its own (two phases, so no segment's cells depend on another's).
        struct Segment
        {
            std::vector<std::size_t> cells;
            std::vector<float> beds;
            std::size_t receiver = nocell;
            float discharge = 0.0f;
            int32_t januaryflow = 0, julyflow = 0;
            int depth = 0;
            std::uint64_t key = 0;
        };
        std::vector<Segment> segments;
        for (const ParentCell& cell : cells)
        {
            if (cell.depth < 1)
                continue;
            const Cell node = nodeof(cell);
            const Cell parentcell{ cell.column, cell.row };
            Cell target;
            bool hastarget = false;
            bool tonode = false;
            const auto found = lookup.find(cell.down);
            if (found != lookup.end())
            {
                target = nodeof(cells[static_cast<std::size_t>(found->second)]);
                hastarget = true;
                tonode = true;
            }
            else
            {
                // The nearest land link cell among the children of the parent cell and its 8 neighbours.
                std::int64_t bestdistance = 0;
                std::uint64_t bestkey = 0;
                for (std::int64_t y = (cell.row - 1) * f; y < (cell.row + 2) * f; y++)
                    for (std::int64_t x = (cell.column - 1) * f; x < (cell.column + 2) * f; x++)
                    {
                        std::size_t other = 0;
                        if (!grid.index(x, y, other) || kind[other] != linkcell || !island(other))
                            continue;
                        const std::int64_t ox = grid.dx(node.x, x), oy = y - node.y;
                        const std::int64_t distance = ox * ox + oy * oy;
                        const std::uint64_t key = grid.key({ x, y });
                        if (!hastarget || distance < bestdistance || (distance == bestdistance && key < bestkey))
                        {
                            target = { node.x + ox, y };
                            bestdistance = distance;
                            bestkey = key;
                            hastarget = true;
                        }
                    }
            }
            if (!hastarget)
                continue;
            const std::vector<Cell> line = linecells(node, grid.dx(node.x, target.x), target.y - node.y);
            Segment segment;
            bool drawn = true;
            for (std::size_t m = 0; m < line.size(); m++)
            {
                const Cell& step = line[m];
                std::size_t index = 0;
                if (!grid.index(step.x, step.y, index) || !nearparent(parentof(step), parentcell))
                {
                    drawn = false;
                    break;
                }
                if (sea[index] != 0 || kind[index] == linkcell || (tonode && m + 1 == line.size()))
                {
                    // The sea, a link or the downstream node receives it (a node in the sea: not drawn).
                    segment.receiver = index;
                    drawn = m > 0 || sea[index] == 0;
                    break;
                }
                segment.cells.push_back(index);
                const std::size_t link = linkneighbour(step);
                if (link != nocell)
                {
                    segment.receiver = link;
                    break;
                }
            }
            if (!drawn || segment.receiver == nocell)
                continue;
            // Beds: the parent's link bed, held at or below the running minimum of the surface.
            float running = std::numeric_limits<float>::infinity();
            for (const std::size_t index : segment.cells)
            {
                if (std::isfinite(surface[index]))
                    running = std::min(running, surface[index]);
                segment.beds.push_back(std::min(cell.linkbed, running));
            }
            segment.discharge = pch.discharge_m3s[cell.index];
            segment.januaryflow = pch.januaryflow[cell.index];
            segment.julyflow = pch.julyflow[cell.index];
            if (tonode)
            {
                segment.discharge = std::min(segment.discharge, pch.discharge_m3s[cell.down]);
                segment.januaryflow = std::min(segment.januaryflow, pch.januaryflow[cell.down]);
                segment.julyflow = std::min(segment.julyflow, pch.julyflow[cell.down]);
            }
            else if (kind[segment.receiver] == linkcell)
            {
                segment.discharge = std::min(segment.discharge, channels.discharge_m3s[segment.receiver]);
                segment.januaryflow = std::min(segment.januaryflow, channels.januaryflow[segment.receiver]);
                segment.julyflow = std::min(segment.julyflow, channels.julyflow[segment.receiver]);
            }
            segment.depth = cell.depth;
            segment.key = levelkey(cell.column, cell.row, plw);
            segments.push_back(std::move(segment));
        }
        report.deadparents = report.inheritedparents - static_cast<int>(segments.size());
        // A cell on several segments belongs to the most downstream one (the smaller depth), then
        // to the smaller parent key.
        std::vector<int32_t> owner(count, -1);
        for (std::size_t s = 0; s < segments.size(); s++)
            for (const std::size_t index : segments[s].cells)
            {
                const int32_t current = owner[index];
                const Segment& a = segments[s];
                if (current < 0 || a.depth < segments[static_cast<std::size_t>(current)].depth ||
                    (a.depth == segments[static_cast<std::size_t>(current)].depth &&
                     a.key < segments[static_cast<std::size_t>(current)].key))
                    owner[index] = static_cast<int32_t>(s);
            }
        for (std::size_t s = 0; s < segments.size(); s++)
        {
            const Segment& segment = segments[s];
            for (std::size_t m = 0; m < segment.cells.size(); m++)
            {
                const std::size_t index = segment.cells[m];
                if (owner[index] != static_cast<int32_t>(s))
                    continue;
                const std::size_t next = m + 1 < segment.cells.size() ? segment.cells[m + 1] : segment.receiver;
                kind[index] = inheritedcell;
                channels.channel[index] = 1;
                channels.direction[index] = stepdirection(index, next);
                channels.bed_m[index] = segment.beds[m];
                channels.discharge_m3s[index] = segment.discharge;
                channels.januaryflow[index] = segment.januaryflow;
                channels.julyflow[index] = segment.julyflow;
                channels.globalcell[index] = -1;
                report.inheritedcells++;
            }
        }
        // An inherited cell draining onto neither a channel nor the sea (its downstream node was
        // not drawn) ends there (direction 0); a loop of inherited cells (segments crossing) is cut
        // at its cell of the smallest key. Both are pure functions of the cells involved.
        std::vector<uint8_t> visit(count, 0);  // 0 unseen, 1 on the current walk, 2 done
        for (std::size_t start = 0; start < count; start++)
        {
            if (kind[start] != inheritedcell || visit[start] != 0)
                continue;
            std::vector<std::size_t> walk;
            std::size_t current = start;
            while (current != nocell && kind[current] == inheritedcell && visit[current] == 0)
            {
                visit[current] = 1;
                walk.push_back(current);
                const std::size_t next = downstream(grid, channels, current);
                if (next == nocell || (channels.channel[next] == 0 && sea[next] == 0))
                {
                    channels.direction[current] = 0;
                    report.inheritedends++;
                    current = nocell;
                    break;
                }
                current = next;
            }
            if (current != nocell && kind[current] == inheritedcell && visit[current] == 1)
            {
                // A loop: from current around to the walk's end.
                std::size_t cut = current;
                std::uint64_t cutkey = grid.key(grid.cell(current));
                for (auto it = std::find(walk.begin(), walk.end(), current); it != walk.end(); ++it)
                {
                    const std::uint64_t key = grid.key(grid.cell(*it));
                    if (key < cutkey)
                    {
                        cut = *it;
                        cutkey = key;
                    }
                }
                channels.direction[cut] = 0;
                report.inheritedends++;
            }
            for (const std::size_t cell : walk)
                visit[cell] = 2;
        }
    }

    // 2. Grown (header, item 2). A junction's rank: links first, then its Strahler order among the
    // inherited cells inside its parent cell and the 8 neighbours.
    const auto localorder = [&](std::size_t junction) {
        const Cell home = parentof(grid.cell(junction));
        std::unordered_map<std::size_t, int> memo;
        // Post-order over the upstream inherited cells inside the block.
        std::vector<std::pair<std::size_t, bool>> stack{ { junction, false } };
        while (!stack.empty())
        {
            const auto [cell, expanded] = stack.back();
            stack.pop_back();
            const Cell c = grid.cell(cell);
            std::vector<std::size_t> upstream;
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const Cell neighbour{ c.x + directionoffsetx(direction), c.y + directionoffsety(direction) };
                std::size_t other = 0;
                if (!grid.index(neighbour.x, neighbour.y, other) || kind[other] != inheritedcell ||
                    !nearparent(parentof(neighbour), home) || downstream(grid, channels, other) != cell)
                    continue;
                upstream.push_back(other);
            }
            if (!expanded)
            {
                stack.push_back({ cell, true });
                for (const std::size_t other : upstream)
                    if (memo.find(other) == memo.end())
                        stack.push_back({ other, false });
                continue;
            }
            int highest = 0, times = 0;
            for (const std::size_t other : upstream)
            {
                const auto found = memo.find(other);
                const int value = found == memo.end() ? 1 : found->second;
                if (value > highest)
                {
                    highest = value;
                    times = 1;
                }
                else if (value == highest)
                    times++;
            }
            memo[cell] = highest == 0 ? 1 : (times >= 2 ? highest + 1 : highest);
        }
        return memo[junction];
    };
    const double hackexponent = parameters.hackexponent;
    const double hackcoefficient = parameters.hackcoefficient;
    const double angle = parameters.junctionangle_deg * pi / 180.0;
    const bool spacingfield = !inputs.valleyspacing_km.empty() &&
                              inputs.valleyspacing_km.size() == static_cast<std::size_t>(inputs.width) * inputs.height;
    std::vector<Candidate> candidates;
    std::unordered_map<std::uint64_t, int> slots;  // slot -> the best valid candidate
    for (std::size_t junction = 0; junction < count; junction++)
    {
        if (kind[junction] != linkcell && kind[junction] != inheritedcell)
            continue;
        if (!island(junction) || !std::isfinite(channels.bed_m[junction]))
            continue;
        const std::size_t next = downstream(grid, channels, junction);
        if (next == nocell)
            continue;
        const Cell j = grid.cell(junction);
        const float junctionbed = channels.bed_m[junction];
        // The main upstream channel neighbour: the largest discharge, then the smallest key.
        std::size_t up = nocell;
        std::uint64_t upkey = 0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const Cell neighbour{ j.x + directionoffsetx(direction), j.y + directionoffsety(direction) };
            std::size_t other = 0;
            if (!grid.index(neighbour.x, neighbour.y, other) || channels.channel[other] == 0 ||
                downstream(grid, channels, other) != junction)
                continue;
            const std::uint64_t key = grid.key(neighbour);
            if (up == nocell || channels.discharge_m3s[other] > channels.discharge_m3s[up] ||
                (channels.discharge_m3s[other] == channels.discharge_m3s[up] && key < upkey))
            {
                up = other;
                upkey = key;
            }
        }
        const Cell n = grid.cell(next);
        const Cell from = up != nocell ? grid.cell(up) : j;
        const double fx = static_cast<double>(grid.dx(from.x, n.x)), fy = static_cast<double>(n.y - from.y);
        const double flength = std::hypot(fx, fy);
        if (!(flength > 0.0))
            continue;
        // Upstream unit vector.
        const double upx = -fx / flength, upy = -fy / flength;
        const std::size_t globalindex =
            static_cast<std::size_t>(j.y / k) * static_cast<std::size_t>(inputs.width) +
            static_cast<std::size_t>(wraplevel(j.x, grid.planetwidth) / k);
        double spacing = parameters.spacing_km;
        if (spacingfield && inputs.valleyspacing_km[globalindex] > 0.0f && std::isfinite(inputs.valleyspacing_km[globalindex]))
            spacing = inputs.valleyspacing_km[globalindex];
        const double strip_km = std::max(spacing, edge * cellsize_km);
        const double hack_km = std::pow(hackcoefficient * std::pow(strip_km, hackexponent), 1.0 / (1.0 - hackexponent));
        const double lengthcells = std::min(static_cast<double>(edge), hack_km / cellsize_km);
        const std::uint64_t junctionkey = grid.key(j);
        int rank = -1;  // computed once a candidate is valid
        for (int side = -1; side <= 1; side += 2)
        {
            const std::uint64_t hash = deterministiccombine(seed, candidatesalt, level, static_cast<std::int64_t>(junctionkey), side);
            const double theta = side * angle * (1.0 + parameters.anglespread * signedunit(hash));
            const double ux = upx * std::cos(theta) - upy * std::sin(theta);
            const double uy = upx * std::sin(theta) + upy * std::cos(theta);
            const std::uint64_t lengthhash = deterministicmix64(hash);
            const double lengthfactor =
                parameters.lengthminfraction + (1.0 - parameters.lengthminfraction) * unitinterval(lengthhash);
            const double target = std::max(2.0, lengthcells * lengthfactor);
            const std::vector<Cell> ray = raycells(j, ux, uy, target);
            Candidate candidate;
            candidate.junction = junction;
            double along = 0.0;
            for (std::size_t m = 1; m < ray.size(); m++)
            {
                const Cell& cell = ray[m];
                const bool diagonal = cell.x != ray[m - 1].x && cell.y != ray[m - 1].y;
                const double length = along + (diagonal ? std::sqrt(2.0) : 1.0);
                if (!(length < static_cast<double>(edge)))
                    break;
                std::size_t index = 0;
                if (!grid.index(cell.x, cell.y, index) || sea[index] != 0 || channels.channel[index] != 0)
                    break;
                // The ground under the junction's bed: the bed would stand above it (item 3).
                if (!(surface[index] >= junctionbed))
                    break;
                bool touches = false;
                for (int direction = 1; direction <= directioncount && !touches; direction++)
                {
                    const Cell neighbour{ cell.x + directionoffsetx(direction), cell.y + directionoffsety(direction) };
                    std::size_t other = 0;
                    if (!grid.index(neighbour.x, neighbour.y, other) || channels.channel[other] == 0)
                        continue;
                    touches = m >= 2 || !(std::llabs(grid.dx(j.x, neighbour.x)) <= 1 && std::llabs(neighbour.y - j.y) <= 1);
                }
                if (touches)
                    break;
                along = length;
                candidate.path.push_back(index);
                candidate.along.push_back(along);
            }
            if (candidate.path.empty() || candidate.along.back() < 2.0)
                continue;
            report.candidates++;
            if (rank < 0)
                rank = kind[junction] == linkcell ? 255 : std::min(254, localorder(junction));
            const Cell slotcell = parentof(j);
            candidate.slot = levelkey(slotcell.x, slotcell.y, parentwidth) * 2u + (side > 0 ? 1u : 0u);
            candidate.priority = (static_cast<std::uint64_t>(rank) << 56) | (deterministicmix64(lengthhash) >> 8);
            // The slot's opening (the density target).
            const double open = std::min(1.0, edge * cellsize_km / spacing);
            if (open < 1.0 && !(unitinterval(deterministiccombine(seed, slotsalt, level, static_cast<std::int64_t>(candidate.slot))) < open))
                continue;
            const int id = static_cast<int>(candidates.size());
            const auto slot = slots.find(candidate.slot);
            if (slot == slots.end())
                slots.emplace(candidate.slot, id);
            else if (candidate.priority > candidates[static_cast<std::size_t>(slot->second)].priority)
                slot->second = id;
            candidates.push_back(std::move(candidate));
        }
    }
    // The winners, one round of the conflict rule.
    std::vector<int> winners;
    for (const auto& slot : slots)
        winners.push_back(slot.second);
    std::sort(winners.begin(), winners.end());
    report.slotwinners = static_cast<int>(winners.size());
    std::unordered_map<std::size_t, Cover> covers;
    for (const int id : winners)
    {
        const Candidate& candidate = candidates[static_cast<std::size_t>(id)];
        for (const std::size_t index : candidate.path)
        {
            const Cell cell = grid.cell(index);
            for (std::int64_t oy = -1; oy <= 1; oy++)
                for (std::int64_t ox = -1; ox <= 1; ox++)
                {
                    std::size_t other = 0;
                    if (grid.index(cell.x + ox, cell.y + oy, other))
                        covers[other].add(candidate.priority, id);
                }
        }
    }
    std::vector<int> grown;
    for (const int id : winners)
    {
        const Candidate& candidate = candidates[static_cast<std::size_t>(id)];
        bool beaten = false;
        for (const std::size_t index : candidate.path)
        {
            const Cover& cover = covers[index];
            const std::uint64_t other = cover.firstid == id ? (cover.secondid >= 0 ? cover.second : 0) : cover.first;
            const int otherid = cover.firstid == id ? cover.secondid : cover.firstid;
            if (otherid >= 0 && other > candidate.priority)
            {
                beaten = true;
                break;
            }
        }
        if (!beaten)
            grown.push_back(id);
    }

    // 3. The grown tributaries' cells: beds b + w (max(r, b) - b), r the running minimum of the
    // surface from the tip (at least b along the path: item 2's cut), w the path length from J
    // over the tip's.
    const bool withrunoff = inputs.runoff_mm_yr.size() == static_cast<std::size_t>(inputs.width) * inputs.height;
    for (const int id : grown)
    {
        const Candidate& candidate = candidates[static_cast<std::size_t>(id)];
        const std::size_t junction = candidate.junction;
        const Cell j = grid.cell(junction);
        const std::size_t globalindex =
            static_cast<std::size_t>(j.y / k) * static_cast<std::size_t>(inputs.width) +
            static_cast<std::size_t>(wraplevel(j.x, grid.planetwidth) / k);
        const double runoff_m = withrunoff ? std::max(0.0, static_cast<double>(inputs.runoff_mm_yr[globalindex])) / 1000.0 : 0.0;
        const double junctiondischarge = std::max(0.0, static_cast<double>(channels.discharge_m3s[junction]));
        const double floor = channels.bed_m[junction];
        const std::size_t size = candidate.path.size();
        const double total = candidate.along.back();
        float running = std::numeric_limits<float>::infinity();
        for (std::size_t t = 0; t < size; t++)  // from the tip down
        {
            const std::size_t m = size - 1 - t;
            const std::size_t index = candidate.path[m];
            const std::size_t next = m > 0 ? candidate.path[m - 1] : junction;
            running = std::min(running, surface[index]);
            const double weight = candidate.along[m] / total;
            const double bed = floor + weight * (std::max(static_cast<double>(running), floor) - floor);
            const double fromtip = total - candidate.along[m] + 1.0;  // cells, the tip's own included
            const double area_km2 = std::pow(fromtip * cellsize_km / hackcoefficient, 1.0 / hackexponent);
            const double discharge = std::min(junctiondischarge, runoff_m * area_km2 * 1.0e6 / secondsperyear);
            const double ratio = junctiondischarge > 0.0 ? discharge / junctiondischarge : 0.0;
            kind[index] = growncell;
            channels.channel[index] = 1;
            channels.direction[index] = stepdirection(index, next);
            channels.bed_m[index] = static_cast<float>(bed);
            channels.discharge_m3s[index] = static_cast<float>(discharge);
            channels.januaryflow[index] = static_cast<int32_t>(channels.januaryflow[junction] * ratio);
            channels.julyflow[index] = static_cast<int32_t>(channels.julyflow[junction] * ratio);
            channels.globalcell[index] = -1;
            report.growncells++;
        }
    }
    report.grown = static_cast<int>(grown.size());
    for (std::size_t cell = 0; cell < count; cell++)
        if ((kind[cell] == inheritedcell || kind[cell] == growncell) && channels.bed_m[cell] > state.surface_m[cell])
            report.bedsabove++;
    const std::vector<uint8_t> finalorders = orders(grid, channels, sea);
    for (const uint8_t value : finalorders)
        report.maxorder = std::max<int>(report.maxorder, value);
    state.tributaries = report.grown;
    state.tributarycells = report.inheritedcells + report.growncells;
    state.tributarybedsabove = report.bedsabove;
    return report;
}
}
