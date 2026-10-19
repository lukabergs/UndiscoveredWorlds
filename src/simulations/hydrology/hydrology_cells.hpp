#pragma once

// Cell queries of the hydrology files on the world's cell grid (docs/architecture/
// sphere-mesh.md, section 3): neighbourhoods from the neighbour table, never from
// x and y offsets, so the same code runs on a cube and on the lat-lon grid. A slot
// without a neighbour (beyond a lat-lon pole, a cube corner's missing diagonal) is
// skipped, as the (x, y) loops skipped the rows beyond a pole.

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "world/cell_queries.hpp"
#include "world/planet.hpp"

namespace hydrologycells
{
// The 3 x 3 neighbourhood in the order of the legacy loops (x outer, y inner):
// north-west, west, south-west, north, the cell itself (0), south, north-east,
// east, south-east. Where the first match decides, the order keeps the lat-lon
// results of those loops.
constexpr std::array<int, 9> legacysquare = { 8, 7, 6, 1, 0, 5, 2, 3, 4 };

// The neighbourhood's cells in that order, the missing neighbours left out.
inline std::vector<mesh::CellId> square(const mesh::CellGrid& grid, mesh::CellId c)
{
    std::vector<mesh::CellId> cells;
    cells.reserve(legacysquare.size());
    for (const int slot : legacysquare)
    {
        if (slot == 0)
        {
            cells.push_back(c);
            continue;
        }
        const std::int32_t neighbour = grid.neighbour(c, slot);
        if (neighbour >= 0)
            cells.push_back(static_cast<mesh::CellId>(neighbour));
    }
    return cells;
}

// The neighbour in a slot, or the cell itself where none exists.
inline mesh::CellId neighbourorself(const mesh::CellGrid& grid, mesh::CellId c, int slot)
{
    const std::int32_t neighbour = grid.neighbour(c, slot);
    return neighbour < 0 ? c : static_cast<mesh::CellId>(neighbour);
}

// The cells in the visiting order of the legacy (x, y) loops: column by column on
// a lat-lon grid (x outer, y inner), so that loops which draw from the shared
// random sequence or read what they wrote keep their lat-lon results; id order on
// a cube.
inline std::vector<mesh::CellId> visitorder(const mesh::CellGrid& grid)
{
    std::vector<mesh::CellId> cells;
    cells.reserve(grid.cellcount());
    if (grid.kind() == mesh::GridKind::latlon && !grid.issubgrid())
    {
        for (int x = 0; x < grid.columns(); x++)
            for (int y = 0; y < grid.rows(); y++)
                cells.push_back(grid.cellof(x, y));
        return cells;
    }
    for (std::size_t c = 0; c < grid.cellcount(); c++)
        cells.push_back(static_cast<mesh::CellId>(c));
    return cells;
}

// The shared planet-field queries (world/cell_queries.hpp) under the hydrology files' names.
inline bool islake(const planet& world, mesh::CellId c)
{
    return worldcells::islakecell(world, c);
}

inline bool isdelta(const planet& world, mesh::CellId c)
{
    return worldcells::isdeltacell(world, c);
}

// Open water as planet::seawrap reads it: the base at or below the sea level and no lake surface.
inline bool openwater(const planet& world, mesh::CellId c)
{
    return world.basebelowsea(c) && world.lakesurfacefloat(c) == 0.0f;
}

// planet::outline on the cell grid: a land cell (planet::sea) with open water
// (openwater) across one of its four edges. Unlike planet::outline, the lat-lon
// polar rows are not left out (their missing slot is skipped).
inline bool outline(const planet& world, mesh::CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    if (world.sea(c))
        return false;
    for (const int slot : { 7, 1, 3, 5 })
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);
        if (neighbour >= 0 && openwater(world, static_cast<mesh::CellId>(neighbour)))
            return true;
    }
    return false;
}
}
