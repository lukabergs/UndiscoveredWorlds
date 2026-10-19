#ifndef physical_layers_internal_hpp
#define physical_layers_internal_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>
#include "core/grid_coordinates.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "world/planet.hpp"

// Shared planet-grid queries for derived physical layers.
namespace physical_layers::detail
{
inline float clampunit(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

inline int clampscore(float value)
{
    return static_cast<int>(std::round(std::clamp(value, 0.0f, 100.0f)));
}

// The (x, y) forms hydrology's drainage_basins.cpp reads until its W2-H port (which no
// longer includes this header): coordinates use the planet's inclusive bounds; longitude
// wraps and latitude never wraps. The resource layers use physical_layers::cellqueries.
constexpr std::array<std::array<int, 2>, 8> drainageoffsets =
{ {
    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
} };

inline twointegers wrappeddestination(const planet& world, int x, int y, int dir)
{
    twointegers destination = getdestination(x, y, dir);
    destination.x = wrap(destination.x, world.width());
    destination.y = std::clamp(destination.y, 0, world.height());
    return destination;
}

inline bool islakecell(const planet& world, int x, int y)
{
    return world.truelake(x, y) != 0 || world.riftlakesurfacefloat(x, y) != 0;
}

inline bool isdeltacell(const planet& world, int x, int y)
{
    return world.deltadir(x, y) != 0 || world.deltajan(x, y) != 0 || world.deltajul(x, y) != 0;
}
}

// The same queries on the world's cell grid (docs/architecture/sphere-mesh.md sections 3
// and 4): cells are ids, neighbours come from the neighbour table (slots 1..8 in the cell's
// own frame, -1 where none exists), a direction code is a slot of its cell. On the lat-lon
// grid every query reads the cells the (x, y) loops read (x wrapped, y clipped at the poles).
// A namespace of its own: terrain's surface_potentials.cpp declares same-named cell queries
// and reads physical_layers::detail.
namespace physical_layers::cellqueries
{
using mesh::CellId;

// Cells per worker of the parallel loops over cells.
constexpr int cellsperworker = 4096;

// fn(first, last) over disjoint inclusive ranges of cell ids covering the grid.
template <typename Fn>
inline void parallelforcells(const planet& world, Fn&& fn)
{
    parallelforrows(0, static_cast<int>(world.cellcount()) - 1, [&](int first, int last)
    {
        fn(static_cast<CellId>(first), static_cast<CellId>(last));
    }, cellsperworker);
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

// planet::outline and planet::coast on any grid (terrain's cell-grid forms).
inline bool isoutline(const planet& world, CellId c)
{
    if (world.grid().kind() == mesh::GridKind::latlon)
        return world.outline(world.grid().column(c), world.grid().row(c));
    return terrain::coastgrid::outline(world, c);
}

inline bool iscoast(const planet& world, CellId c)
{
    if (world.grid().kind() == mesh::GridKind::latlon)
        return world.coast(world.grid().column(c), world.grid().row(c));
    return terrain::coastgrid::coast(world, c);
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

template <typename Predicate>
inline bool hasfeatureinradius(const planet& world, CellId c, int radius, Predicate predicate)
{
    return anyinbox(world, c, radius, [&](CellId n) { return static_cast<bool>(predicate(n)); });
}

// Downstream cell of each land or delta cell, by cell id; -1 where no neighbour exists
// in the direction or no lower neighbour exists. Rivers and deltas give the direction;
// elsewhere the steepest descent does. In an eroded world (erosion_duration > 0) a land
// cell drains along the routing's flow_direction (an outlet's code points into its sea
// cell); a land cell without a code there (sea when the routing ran, or never reached)
// and a sea delta cell keep the rule above. slots, when given, receives the slot (1..8)
// of each receiver (0 without one).
inline std::vector<int32_t> drainagereceivers(const planet& world, std::vector<std::int8_t>* slots = nullptr)
{
    const mesh::CellGrid& grid = world.grid();
    const std::size_t count = grid.cellcount();
    std::vector<int32_t> receiver(count, -1);
    if (slots != nullptr)
        slots->assign(count, 0);
    const worldgrid<std::uint8_t>* flow = world.erosionduration() > 0.0f ? &world.flowdirection() : nullptr;
    for (CellId c = 0; c < count; ++c)
    {
        if (world.sea(c) != 0 && !isdeltacell(world, c)) continue;
        int direction = 0;
        if (flow != nullptr && world.sea(c) == 0)
            direction = (*flow)[c];
        if (direction < 1 || direction > 8)
        {
            direction = world.riverdir(c);
            if (world.deltadir(c) > 0)
                direction = world.deltadir(c);
        }
        if (direction < 1 || direction > 8)
        {
            float lowest = world.nomfloat(c);
            for (int candidate = 1; candidate <= 8; ++candidate)
            {
                const std::int32_t next = grid.neighbour(c, candidate);
                if (next < 0) continue;
                const float elevation = world.nomfloat(static_cast<CellId>(next));
                if (elevation < lowest) { lowest = elevation; direction = candidate; }
            }
        }
        if (direction < 1 || direction > 8) continue;
        const std::int32_t next = grid.neighbour(c, direction);
        if (next < 0 || static_cast<CellId>(next) == c) continue;
        receiver[c] = next;
        if (slots != nullptr)
            (*slots)[c] = static_cast<std::int8_t>(direction);
    }
    return receiver;
}
}

#endif
