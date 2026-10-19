#ifndef physical_layers_internal_hpp
#define physical_layers_internal_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <vector>
#include "core/grid_coordinates.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "world/cell_queries.hpp"
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

}

// The same queries on the world's cell grid (docs/architecture/sphere-mesh.md sections 3
// and 4): cells are ids, neighbours come from the neighbour table (slots 1..8 in the cell's
// own frame, -1 where none exists), a direction code is a slot of its cell. On the lat-lon
// grid every query reads the cells the (x, y) loops read (x wrapped, y clipped at the poles).
// The planet-field queries and the box walk are the shared ones of world/cell_queries.hpp;
// a namespace of its own, read by the resource files and terrain's surface_potentials.cpp.
namespace physical_layers::cellqueries
{
using mesh::CellId;
using worldcells::annualrainfall;
using worldcells::anyinbox;
using worldcells::averageriverflow;
using worldcells::hasfeatureinradius;
using worldcells::isdeltacell;
using worldcells::islakecell;
using worldcells::issaltfeature;

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

// planet::outline and planet::coast on any grid (terrain's cell-grid forms). On lat-lon
// they differ from the planet's only in the first and last columns, whose west and east
// neighbours planet::wrapx takes one column short (it wraps by W - 1).
inline bool isoutline(const planet& world, CellId c)
{
    return terrain::coastgrid::outline(world, c);
}

inline bool iscoast(const planet& world, CellId c)
{
    return terrain::coastgrid::coast(world, c);
}

// f(cell) for the cells within steps legacy pixels (terrain.sphere_noise.legacy_columns) of c's
// centre, c included: the legacy box of half-width steps as the disc of steps + 0.5 pixels
// (terrain::cellswithinkm), the same in km at every resolution (T8 W2-T3; at 128 columns,
// 312.7 km cells, the discs of one or two pixels hold only c). A callback that returns bool
// stops the walk when it returns true; the result says whether it did.
template <typename F>
inline bool anywithinpixels(const planet& world, CellId c, int steps, F&& f)
{
    for (const CellId n : terrain::cellswithinkm(world, c, (steps + 0.5) * terrain::legacycellkm(world)))
    {
        if constexpr (std::is_same_v<decltype(f(CellId{})), bool>)
        {
            if (f(n))
                return true;
        }
        else
            f(n);
    }
    return false;
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
