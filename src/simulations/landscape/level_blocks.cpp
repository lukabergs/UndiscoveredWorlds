#include "simulations/landscape/level_blocks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/refinement_grid.hpp"

namespace landscape
{
namespace
{
// The level's size on lat-lon: its columns (rows are columns / 2).
std::uint64_t latloncolumns(const mesh::CellGrid& grid, int factor)
{
    const std::uint64_t columns = static_cast<std::uint64_t>(grid.columns());
    if (factor < 1 || columns == 0 || columns % static_cast<std::uint64_t>(factor) != 0)
        throw std::invalid_argument("level blocks: the factor must divide the level's columns");
    return columns;
}

void checkcube(const mesh::CellGrid& grid, int factor)
{
    if (factor < 1 || grid.facecells() % factor != 0)
        throw std::invalid_argument("level blocks: the factor must divide the level's face cells");
}
}

mesh::LevelKey parentkey(const mesh::CellGrid& grid, mesh::LevelKey key, int factor)
{
    if (grid.kind() == mesh::GridKind::cube)
    {
        checkcube(grid, factor);
        return mesh::cubeparentkey(grid.facecells(), factor, key);
    }
    const std::uint64_t columns = latloncolumns(grid, factor), f = static_cast<std::uint64_t>(factor);
    return (key / columns / f) * (columns / f) + (key % columns) / f;
}

std::vector<mesh::LevelKey> childkeys(const mesh::CellGrid& grid, mesh::LevelKey parent, int factor)
{
    if (grid.kind() == mesh::GridKind::cube)
    {
        checkcube(grid, factor);
        return mesh::cubechildkeys(grid.facecells() / factor, factor, parent);
    }
    const std::uint64_t columns = latloncolumns(grid, factor), f = static_cast<std::uint64_t>(factor);
    const std::uint64_t parentcolumns = columns / f;
    const std::uint64_t row0 = (parent / parentcolumns) * f, column0 = (parent % parentcolumns) * f;
    std::vector<mesh::LevelKey> children;
    children.reserve(static_cast<std::size_t>(f * f));
    for (std::uint64_t dj = 0; dj < f; dj++)
        for (std::uint64_t di = 0; di < f; di++)
            children.push_back((row0 + dj) * columns + column0 + di);
    return children;
}

bool firstchild(const mesh::CellGrid& grid, mesh::LevelKey key, int factor)
{
    const std::int64_t f = factor;
    if (grid.kind() == mesh::GridKind::cube)
    {
        checkcube(grid, factor);
        const mesh::CubeCell cell = mesh::cubecell(grid.facecells(), key);
        return cell.column % f == 0 && cell.row % f == 0;
    }
    const std::uint64_t columns = latloncolumns(grid, factor);
    return static_cast<std::int64_t>(key % columns) % f == 0 && static_cast<std::int64_t>(key / columns) % f == 0;
}

WindowBlocks windowblocks(const mesh::CellGrid& window, int factor)
{
    WindowBlocks blocks;
    blocks.factor = factor;
    const std::size_t count = window.cellcount();
    blocks.children.reserve(count);
    for (mesh::CellId local = 0; local < count; local++)
    {
        const mesh::LevelKey key = window.globalkey(local);
        if (!firstchild(window, key, factor))
            continue;
        const mesh::LevelKey parent = parentkey(window, key, factor);
        blocks.parents.push_back(parent);
        for (const mesh::LevelKey child : childkeys(window, parent, factor))
        {
            const mesh::CellId id = window.localofkey(child);
            if (id == mesh::nocell)
                throw std::invalid_argument("windowblocks: a parent's child lies outside the window");
            blocks.children.push_back(id);
        }
    }
    if (blocks.children.size() != count)
        throw std::invalid_argument("windowblocks: the window is not whole parents");
    return blocks;
}

double blockmean(const WindowBlocks& blocks, std::size_t block, const std::vector<float>& field)
{
    const mesh::CellId* children = blocks.childrenof(block);
    double sum = 0.0;
    for (std::size_t n = 0; n < blocks.perparent(); n++)
        sum += static_cast<double>(field[children[n]]);
    return sum / (static_cast<double>(blocks.factor) * blocks.factor);
}

NearestDomainKey::NearestDomainKey(std::int64_t facecells, std::vector<mesh::LevelKey> domain)
    : facecells_(facecells), domain_(std::move(domain))
{
    std::sort(domain_.begin(), domain_.end());
    domain_.erase(std::unique(domain_.begin(), domain_.end()), domain_.end());
    if (facecells_ < 1 || domain_.empty())
        throw std::invalid_argument("NearestDomainKey needs a cube level and a non-empty domain");
}

bool NearestDomainKey::contains(mesh::LevelKey key) const
{
    return std::binary_search(domain_.begin(), domain_.end(), key);
}

std::size_t NearestDomainKey::operator()(mesh::LevelKey key) const
{
    const auto found = std::lower_bound(domain_.begin(), domain_.end(), key);
    if (found != domain_.end() && *found == key)
        return static_cast<std::size_t>(found - domain_.begin());
    const auto cached = std::lower_bound(cache_.begin(), cache_.end(), std::make_pair(key, std::size_t{ 0 }));
    if (cached != cache_.end() && cached->first == key)
        return cached->second;
    // Ring by ring from the key; the first ring holding domain keys gives its smallest.
    std::vector<mesh::LevelKey> ring{ key }, visited{ key }, next;
    std::size_t result = domain_.size();
    for (std::int64_t step = 0; step < 4 * facecells_ && result == domain_.size(); step++)
    {
        next.clear();
        for (const mesh::LevelKey cell : ring)
            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const mesh::LevelKey other = mesh::cubekeyneighbour(facecells_, cell, slot);
                if (other == mesh::nokey || std::find(visited.begin(), visited.end(), other) != visited.end() ||
                    std::find(next.begin(), next.end(), other) != next.end())
                    continue;
                next.push_back(other);
            }
        mesh::LevelKey best = mesh::nokey;
        for (const mesh::LevelKey cell : next)
            if (contains(cell) && cell < best)
                best = cell;
        if (best != mesh::nokey)
            result = static_cast<std::size_t>(std::lower_bound(domain_.begin(), domain_.end(), best) - domain_.begin());
        visited.insert(visited.end(), next.begin(), next.end());
        ring.swap(next);
    }
    if (result == domain_.size())
        throw std::logic_error("NearestDomainKey: no domain key reached");
    cache_.insert(cached, std::make_pair(key, result));
    return result;
}

std::vector<double> gatherblocks(const WindowBlocks& blocks, const std::vector<double>& values, const mesh::CellGrid& window)
{
    if (values.size() != blocks.count())
        throw std::invalid_argument("gatherblocks: a cube window, one value per block and a whole factor");
    return gatherblocks(gatherplan(blocks, window), values);
}

GatherPlan gatherplan(const WindowBlocks& blocks, const mesh::CellGrid& window)
{
    if (window.kind() != mesh::GridKind::cube || blocks.factor < 1 || window.facecells() % blocks.factor != 0)
        throw std::invalid_argument("gatherblocks: a cube window, one value per block and a whole factor");
    GatherPlan plan;
    plan.blocks = blocks.count();
    plan.cells = window.cellcount();
    plan.parentcells = window.facecells() / blocks.factor;
    plan.nearest = std::make_shared<const NearestDomainKey>(plan.parentcells, blocks.parents);
    // The blocks by the parents' ascending order (NearestDomainKey's indices).
    plan.order.resize(blocks.count());
    for (std::size_t b = 0; b < blocks.count(); b++)
        plan.order[b] = static_cast<std::uint32_t>(b);
    std::sort(plan.order.begin(), plan.order.end(),
              [&](std::uint32_t x, std::uint32_t y) { return blocks.parents[x] < blocks.parents[y]; });
    // Per face the bounding box of the blocks' parents, each parent cell holding its block.
    plan.faces.resize(mesh::cubefacecount);
    std::int64_t low[mesh::cubefacecount][2], high[mesh::cubefacecount][2];
    for (int f = 0; f < mesh::cubefacecount; f++)
    {
        low[f][0] = low[f][1] = plan.parentcells;
        high[f][0] = high[f][1] = -1;
    }
    for (const mesh::LevelKey parent : blocks.parents)
    {
        const mesh::CubeCell cell = mesh::cubecell(plan.parentcells, parent);
        low[cell.face][0] = std::min(low[cell.face][0], cell.column);
        low[cell.face][1] = std::min(low[cell.face][1], cell.row);
        high[cell.face][0] = std::max(high[cell.face][0], cell.column);
        high[cell.face][1] = std::max(high[cell.face][1], cell.row);
    }
    for (int f = 0; f < mesh::cubefacecount; f++)
    {
        GatherPlan::FaceTable& table = plan.faces[static_cast<std::size_t>(f)];
        if (high[f][0] < 0)
            continue;
        table.column0 = low[f][0];
        table.row0 = low[f][1];
        table.columns = high[f][0] - low[f][0] + 1;
        table.rows = high[f][1] - low[f][1] + 1;
        table.block.assign(static_cast<std::size_t>(table.columns * table.rows), -1);
    }
    for (std::size_t b = 0; b < blocks.count(); b++)
    {
        const mesh::CubeCell cell = mesh::cubecell(plan.parentcells, blocks.parents[b]);
        GatherPlan::FaceTable& table = plan.faces[static_cast<std::size_t>(cell.face)];
        table.block[static_cast<std::size_t>((cell.row - table.row0) * table.columns + (cell.column - table.column0))] =
            static_cast<std::int32_t>(b);
    }
    // The stencils at the window's centres, kept with the window grid (T8 W2-W3b).
    plan.stencils = &mesh::centrestencils(window, plan.parentcells);
    return plan;
}

std::size_t GatherPlan::blockof(mesh::LevelKey key) const
{
    const mesh::CubeCell cell = mesh::cubecell(parentcells, key);
    const FaceTable& table = faces[static_cast<std::size_t>(cell.face)];
    const std::int64_t x = cell.column - table.column0, y = cell.row - table.row0;
    if (x >= 0 && y >= 0 && x < table.columns && y < table.rows)
    {
        const std::int32_t block = table.block[static_cast<std::size_t>(y * table.columns + x)];
        if (block >= 0)
            return static_cast<std::size_t>(block);
    }
    return order[(*nearest)(key)];
}

std::vector<double> gatherblocks(const GatherPlan& plan, const std::vector<double>& values)
{
    if (values.size() != plan.blocks)
        throw std::invalid_argument("gatherblocks: one value per block of the plan");
    std::vector<double> gathered(plan.cells, 0.0);
    if (plan.stencils == nullptr)
        return gathered;
    // mesh::applystencil's arithmetic: relative to the first node.
    mesh::CentreStencils::Interior nodes;
    const mesh::LevelKey n = static_cast<mesh::LevelKey>(plan.parentcells);
    for (mesh::CellId c = 0; c < plan.cells; c++)
    {
        bool started = false;
        double origin = 0.0, sum = 0.0;
        if (plan.stencils->interior(c, nodes))
        {
            // CentreStencils::visit's on-face loop, the block from the face table when the node is a parent.
            const GatherPlan::FaceTable& table = plan.faces[static_cast<std::size_t>(nodes.face)];
            const mesh::LevelKey base = static_cast<mesh::LevelKey>(nodes.face) * n * n;
            for (int m = 0; m < 4; m++)
                for (int l = 0; l < 4; l++)
                {
                    const double w = nodes.wx[l] * nodes.wy[m];
                    if (w == 0.0)
                        continue;
                    const std::int64_t column = nodes.i0 + l, row = nodes.j0 + m;
                    const std::int64_t x = column - table.column0, y = row - table.row0;
                    std::int32_t block = -1;
                    if (x >= 0 && y >= 0 && x < table.columns && y < table.rows)
                        block = table.block[static_cast<std::size_t>(y * table.columns + x)];
                    const double value =
                        values[block >= 0 ? static_cast<std::size_t>(block)
                                          : plan.order[(*plan.nearest)(base + static_cast<mesh::LevelKey>(row) * n +
                                                                       static_cast<mesh::LevelKey>(column))]];
                    if (!started)
                    {
                        origin = value;
                        started = true;
                        continue;
                    }
                    sum += w * (value - origin);
                }
            gathered[c] = started ? origin + sum : 0.0;
            continue;
        }
        plan.stencils->visit(c, [&](mesh::LevelKey key, double weight) {
            const double value = values[plan.blockof(key)];
            if (!started)
            {
                origin = value;
                started = true;
                return;
            }
            sum += weight * (value - origin);
        });
        gathered[c] = started ? origin + sum : 0.0;
    }
    return gathered;
}

ExportPixel exportpixel(const mesh::Vec3& p, std::int64_t columns)
{
    // CellGrid::locatekey's lat-lon expressions in 64 bits.
    const double pi = 3.14159265358979323846;
    const std::int64_t rows = columns / 2;
    const double latitude = std::asin(std::clamp(p.z, -1.0, 1.0)), longitude = mesh::longitudeof(p);
    ExportPixel pixel;
    pixel.row = std::clamp(static_cast<std::int64_t>(std::floor((0.5 * pi - latitude) / (pi / static_cast<double>(rows)))),
                           std::int64_t{ 0 }, rows - 1);
    std::int64_t x = static_cast<std::int64_t>(std::floor((longitude + pi) / (2.0 * pi / static_cast<double>(columns))));
    x = std::clamp(x, std::int64_t{ 0 }, columns);
    pixel.column = x == columns ? 0 : x;
    return pixel;
}

std::shared_ptr<const mesh::CellGrid> levelkeygrid(const LevelWindow& window)
{
    mesh::LatLonWindowBox box;
    box.column0 = window.domaincolumn0;
    box.columns = window.domaincolumns;
    box.row0 = static_cast<int>(window.domainrow0);
    box.rows = window.domainrows;
    return mesh::latlonwindowgrid(static_cast<int>(window.levelwidth()), 1.0, box);
}
}
