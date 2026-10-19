#pragma once

// Per-cell queries of the planet's fields shared by the terrain, resource and hydrology
// layers (docs/architecture/sphere-mesh.md sections 3 and 4): cells are ids, neighbours come
// from the neighbour table (slots 1..8 in the cell's own frame, -1 where none exists). On the
// lat-lon grid every query reads the cells the legacy (x, y) loops read (x wrapped, y clipped
// at the poles). One copy for every layer (T8 W2-T3); a layer's own queries stay in its own
// header (physical_layers::cellqueries in resources/physical_layers_internal.hpp).

#include <algorithm>
#include <cstdint>
#include <cstdlib>

#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "world/planet.hpp"

namespace worldcells
{
using mesh::CellId;

inline bool iswetlandcell(const planet& world, CellId c)
{
    const int special = world.special(c);
    return special >= 130 && special < 140;
}

inline bool islakecell(const planet& world, CellId c)
{
    return world.truelake(c) != 0 || world.riftlakesurfacefloat(c) != 0;
}

inline bool issaltfeature(const planet& world, CellId c)
{
    const int special = world.special(c);
    return special == 100 || special == 110;
}

inline bool isdeltacell(const planet& world, CellId c)
{
    return world.deltadir(c) != 0 || world.deltajan(c) != 0 || world.deltajul(c) != 0;
}

inline float annualrainfall(const planet& world, CellId c)
{
    return (world.summerrainfloat(c) + world.winterrainfloat(c)) / 2.0f;
}

inline int averageriverflow(const planet& world, CellId c)
{
    const int riverflow = world.riveraveflow(c);
    const int deltaflow = (std::abs(world.deltajan(c)) + std::abs(world.deltajul(c))) / 2;
    return std::max(riverflow, deltaflow);
}

// f(cell) for every cell within one (radius 1) or two (radius 2) neighbour steps of c,
// c included: the 3 x 3 and 5 x 5 boxes of the legacy grid (rows beyond a pole left out).
// A cell may be visited more than once, so f must be idempotent (a maximum, an any).
// A callback that returns bool stops the walk when it returns true.
template <typename F>
inline bool anyinbox(const planet& world, CellId c, int radius, F&& f)
{
    const mesh::CellGrid& grid = world.grid();
    for (const int first : boxslots)
    {
        const std::int32_t inner = first == 0 ? static_cast<std::int32_t>(c) : grid.neighbour(c, first);
        if (inner < 0)
            continue;
        if (radius < 2)
        {
            if (f(static_cast<CellId>(inner)))
                return true;
            continue;
        }
        for (const int second : boxslots)
        {
            const std::int32_t outer = second == 0 ? inner : grid.neighbour(static_cast<CellId>(inner), second);
            if (outer >= 0 && f(static_cast<CellId>(outer)))
                return true;
        }
    }
    return false;
}

// Whether a cell within one or two neighbour steps satisfies the predicate (anyinbox).
template <typename Predicate>
inline bool hasfeatureinradius(const planet& world, CellId c, int radius, Predicate predicate)
{
    return anyinbox(world, c, radius, [&](CellId n) { return static_cast<bool>(predicate(n)); });
}
}
