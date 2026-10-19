// Questions about drainage on the world grid: flatness, flow destinations, inflows, lakes nearby.
// Moved from misc.cpp (Jonathan Hill, 2022); on the cell grid since the sphere mesh (W2-H).

#include "simulations/hydrology/flow_queries.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "simulations/hydrology/hydrology_cells.hpp"
#include "world/planet.hpp"

namespace
{
// The lat-lon cell of legacy (x, y) coordinates: x taken round the globe, nocell
// beyond a pole (the legacy loops skipped those rows).
mesh::CellId latloncell(const planet& world, int x, int y)
{
    const mesh::CellGrid& grid = world.grid();
    if (y < 0 || y >= grid.rows())
        return mesh::nocell;
    const int columns = grid.columns();
    return grid.cellof(((x % columns) + columns) % columns, y);
}

twointegers coordinates(const planet& world, mesh::CellId c)
{
    twointegers point;
    point.x = world.grid().column(c);
    point.y = world.grid().row(c);
    return point;
}
}

// This function checks how flat the land is at a given point: the mean difference of the
// water-level elevation to the west, north-west, north, north-east, south-west, south and
// south-east neighbours, divided by their number plus one (the legacy count, which had a
// place for the east neighbour but never added it).

int getflatness(planet& world, mesh::CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    const float thiselev = getflatelevation(world, c);

    float flatness = 0.0f;
    float crount = 1.0f;

    for (const int slot : { 7, 8, 1, 2, 6, 5, 4 })
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        const float thatelev = getflatelevation(world, static_cast<mesh::CellId>(neighbour));
        flatness = flatness + std::fabs(thiselev - thatelev);
        crount = crount + 1.0f;
    }

    flatness = flatness / crount;

    return (int)flatness;
}

// This function gets the elevation of a point, incorporating water level.

float getflatelevation(planet& world, mesh::CellId c)
{
    float elev;

    float surface = world.lakesurfacefloat(c);

    if (surface > 0)
        elev = surface;
    else
    {
        if (world.sea(c) == 1)
            elev = static_cast<float>(world.sealevel());
        else
            elev = world.nomfloat(c);
    }

    return (elev);
}

// This returns the cell that the river in the given cell is flowing into.

mesh::CellId getflowdestination(planet& world, mesh::CellId c, int dir)
{
    if (dir == 0)
        dir = world.riverdir(c);

    return hydrologycells::neighbourorself(world.grid(), c, dir);
}

// This function finds the cell with the largest flow into the current one.

std::int32_t getupstreamcell(planet& world, mesh::CellId c)
{
    std::int32_t upcell = -1;
    int largest = 0;

    for (const mesh::CellId candidate : hydrologycells::square(world.grid(), c))
    {
        const int riverjan = world.riverjan(candidate);
        const int riverjul = world.riverjul(candidate);

        if (riverjan + riverjul > largest) // If this is larger than the current largest inflow
        {
            if (getflowdestination(world, candidate, 0) == c) // If this is actually flowing into our cell
            {
                upcell = static_cast<std::int32_t>(candidate);
                largest = riverjan + riverjul;
            }
        }
    }
    return (upcell);
}

// This function checks to see whether a cell has any water flowing into it: a neighbour
// whose river direction is the slot by which it reaches this cell.

int checkwaterinflow(planet& world, mesh::CellId c)
{
    const mesh::CellGrid& grid = world.grid();

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour >= 0 && world.riverdir(static_cast<mesh::CellId>(neighbour)) == grid.reverse(c, slot))
            return 1;
    }

    return 0;
}

// This function gets the total water inflow (January, July) into a given cell.

twointegers gettotalinflow(planet& world, mesh::CellId c)
{
    twointegers total;

    total.x = 0;
    total.y = 0;

    for (const mesh::CellId candidate : hydrologycells::square(world.grid(), c))
    {
        if (getflowdestination(world, candidate, 0) == c)
        {
            total.x = total.x + world.riverjan(candidate);
            total.y = total.y + world.riverjul(candidate);
        }
    }
    return total;
}

// This function tells whether a given cell is in or next to a lake: the first lake surface
// found in the 3 x 3 neighbourhood (in the legacy order), then ring by ring outwards for a
// larger dist (each ring in ascending id).

float nearlake(planet& world, mesh::CellId c, int dist, bool rift)
{
    const mesh::CellGrid& grid = world.grid();

    const auto lakeat = [&](mesh::CellId cell)
    {
        if (world.lakesurfacefloat(cell) != 0)
            return world.lakesurfacefloat(cell);

        if (rift == 1 && world.riftlakesurfacefloat(cell) != 0)
            return world.riftlakesurfacefloat(cell);

        return 0.0f;
    };

    if (dist < 1)
        return lakeat(c);

    std::vector<mesh::CellId> seen = hydrologycells::square(grid, c);

    for (const mesh::CellId cell : seen)
    {
        const float surface = lakeat(cell);

        if (surface != 0)
            return surface;
    }

    std::vector<mesh::CellId> ring = seen;
    std::sort(seen.begin(), seen.end());

    for (int step = 2; step <= dist; step++)
    {
        std::vector<mesh::CellId> next;

        for (const mesh::CellId cell : ring)
        {
            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const std::int32_t neighbour = grid.neighbour(cell, slot);

                if (neighbour < 0)
                    continue;

                const mesh::CellId candidate = static_cast<mesh::CellId>(neighbour);

                if (!std::binary_search(seen.begin(), seen.end(), candidate) &&
                    std::find(next.begin(), next.end(), candidate) == next.end())
                    next.push_back(candidate);
            }
        }

        std::sort(next.begin(), next.end());

        for (const mesh::CellId cell : next)
        {
            const float surface = lakeat(cell);

            if (surface != 0)
                return surface;
        }

        seen.insert(seen.end(), next.begin(), next.end());
        std::sort(seen.begin(), seen.end());
        ring = next;
    }
    return 0.0f;
}

// This function tells us whether a cell is a lake cell on the edge.

int getlakeedge(planet& world, mesh::CellId c)
{
    if (world.lakesurfacefloat(c) == 0)
        return 0;

    for (const mesh::CellId candidate : hydrologycells::square(world.grid(), c))
    {
        if (world.lakesurfacefloat(candidate) == 0)
            return 1;
    }
    return 0;
}

/// The (x, y) forms: lat-lon worlds only; a row beyond a pole has no cell (nothing flows,
// no lake, no flatness).

int getflatness(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    return c == mesh::nocell ? 0 : getflatness(world, c);
}

float getflatelevation(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    return c == mesh::nocell ? 0.0f : getflatelevation(world, c);
}

twointegers getflowdestination(planet& world, int x, int y, int dir)
{
    const mesh::CellId c = latloncell(world, x, y);
    if (c == mesh::nocell)
    {
        twointegers point;
        point.x = x;
        point.y = std::clamp(y, 0, world.grid().rows() - 1);
        return point;
    }
    return coordinates(world, getflowdestination(world, c, dir));
}

twointegers getupstreamcell(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    const std::int32_t upcell = c == mesh::nocell ? -1 : getupstreamcell(world, c);

    if (upcell < 0)
    {
        twointegers none;
        none.x = -1;
        none.y = -1;
        return none;
    }
    return coordinates(world, static_cast<mesh::CellId>(upcell));
}

int checkwaterinflow(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    return c == mesh::nocell ? 0 : checkwaterinflow(world, c);
}

twointegers gettotalinflow(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    if (c == mesh::nocell)
    {
        twointegers none;
        none.x = 0;
        none.y = 0;
        return none;
    }
    return gettotalinflow(world, c);
}

float nearlake(planet& world, int x, int y, int dist, bool rift)
{
    const mesh::CellId c = latloncell(world, x, y);
    return c == mesh::nocell ? 0.0f : nearlake(world, c, dist, rift);
}

int getlakeedge(planet& world, int x, int y)
{
    const mesh::CellId c = latloncell(world, x, y);
    return c == mesh::nocell ? 0 : getlakeedge(world, c);
}
