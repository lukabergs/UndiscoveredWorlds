#include "simulations/landscape/level_blocks.hpp"

#include <algorithm>
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
    if (window.kind() != mesh::GridKind::cube || values.size() != blocks.count() || blocks.factor < 1 ||
        window.facecells() % blocks.factor != 0)
        throw std::invalid_argument("gatherblocks: a cube window, one value per block and a whole factor");
    const std::int64_t parentcells = window.facecells() / blocks.factor;
    const NearestDomainKey nearest(parentcells, blocks.parents);
    // values by the parents' ascending order (NearestDomainKey's indices).
    std::vector<std::size_t> order(blocks.count());
    for (std::size_t b = 0; b < blocks.count(); b++)
        order[b] = b;
    std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return blocks.parents[x] < blocks.parents[y]; });
    std::vector<double> sorted(blocks.count());
    for (std::size_t m = 0; m < order.size(); m++)
        sorted[m] = values[order[m]];
    std::vector<double> gathered(window.cellcount());
    for (mesh::CellId c = 0; c < window.cellcount(); c++)
        gathered[c] = mesh::applystencil(mesh::cubicstencil(parentcells, window.centre(c)),
                                         [&](mesh::LevelKey key) { return sorted[nearest(key)]; });
    return gathered;
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
