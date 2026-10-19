#include "simulations/landscape/level_blocks.hpp"

#include <cstdint>
#include <stdexcept>

#include "core/mesh/cube_geometry.hpp"
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
