#pragma once

// Index helpers shared by the regional pass's translation units (erosion plan Phase R2,
// step R6): not part of the pass's interface. regional_surface.cpp keeps its own copies
// (its expressions are what the overlap tests established); the newer units use these.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "core/mesh/cell_grid.hpp"

namespace landscape::regionaldetail
{
// Whether `test` holds for a cell of the 3 x 3 neighbourhood of a level-0 cell on the world's grid,
// the cell included (T8 W2-W2c, sphere-mesh.md section 8: gates by table neighbours): the cell, then
// its slot neighbours 1..8 (none skipped). On the lat-lon grid these are the cells of the raster
// aroundcell below (columns wrap, rows beyond the poles do not exist), so an `any` or a maximum over
// them is the same.
template<typename Test>
bool aroundcell(const mesh::CellGrid& grid, std::size_t cell, Test test)
{
    if (test(cell))
        return true;
    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t other = grid.neighbour(static_cast<mesh::CellId>(cell), slot);
        if (other >= 0 && test(static_cast<std::size_t>(other)))
            return true;
    }
    return false;
}

inline std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

// The cells of a grid within `rings` rings (8-connected hops on its table) of a cell, the cell
// first, each once with its ring (T8 W2-W2d: the cube forms' neighbourhoods). `stamp` is scratch of
// one value per grid cell, `id` a value not yet stamped (the caller counts it up per call).
inline void cellswithin(const mesh::CellGrid& grid, mesh::CellId start, int rings, std::vector<std::uint32_t>& stamp,
                        std::uint32_t id, std::vector<std::pair<mesh::CellId, int>>& out)
{
    out.clear();
    if (stamp.size() != grid.cellcount())
        stamp.assign(grid.cellcount(), 0);
    out.emplace_back(start, 0);
    stamp[start] = id;
    for (std::size_t head = 0; head < out.size(); head++)
    {
        const auto [cell, ring] = out[head];
        if (ring >= rings)
            continue;
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t other = grid.neighbour(cell, slot);
            if (other < 0 || stamp[static_cast<std::size_t>(other)] == id)
                continue;
            stamp[static_cast<std::size_t>(other)] = id;
            out.emplace_back(static_cast<mesh::CellId>(other), ring + 1);
        }
    }
}

// The great-circle distance between two cells' centres in units of `cell_m` (atan2 of the cross and
// dot products, accurate for neighbours).
inline double centredistance(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b, double cell_m)
{
    const mesh::Vec3 p = grid.centre(a), q = grid.centre(b);
    return grid.radius_m() * std::atan2(mesh::length(mesh::cross(p, q)), mesh::dot(p, q)) / cell_m;
}

// Columns wrap around the planet; rows beyond the poles do not exist.
inline int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

// Whether `test` holds for a cell of the 3 x 3 neighbourhood of global cell (i, j), the
// cell included (columns wrap, rows beyond the poles are skipped).
template<typename Test>
bool aroundcell(int width, int height, int i, int j, Test test)
{
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
            if (test(at(wrapcolumn(i + di, width), row, width)))
                return true;
    }
    return false;
}
}
