#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "simulations/hydrology/hydrology_cells.hpp"
#include "simulations/hydrology/tides.hpp"
#include "simulations/terrain/terrain_elevation.hpp"

using namespace std;

// On the world's cell grid: casts follow the neighbour table and measure their fetch in km.

namespace
{
// The legacy cast was 1000 cells at most, and each open-water cell of it added 0.01 to the
// range: lengths of legacy pixels (terrain::legacycellkm, a cell of the 2048-column world).
constexpr double castlimitpixels = 1000.0;
constexpr double rangeperpixel = 0.01;

// The open-water fetch in km of one cast from c: along a slot, straight on through each cell
// (the slot opposite the one by which the cast entered, in that cell's own frame), adding the
// centre distances of the cells below the sea level until land, a missing neighbour or the
// cast's limit.
double castfetch_km(const planet& world, mesh::CellId c, int slot, double limit_km)
{
    const mesh::CellGrid& grid = world.grid();
    double fetch = 0.0;
    mesh::CellId current = c;
    int direction = slot;

    while (fetch < limit_km)
    {
        const std::int32_t neighbour = grid.neighbour(current, direction);

        if (neighbour < 0)
            break;

        const mesh::CellId next = static_cast<mesh::CellId>(neighbour);

        if (!world.basebelowsea(next))
            break;

        fetch += grid.distance_m(current, direction) / 1000.0;
        direction = (grid.reverse(current, direction) + 3) % 8 + 1;
        current = next;
    }

    return std::min(fetch, limit_km);
}

// This function finds the tidal range of a given point: casts in each of the eight directions
// until they hit land; the longer each cast, the more tidal range it generates.

int gettidalrange(planet& world, mesh::CellId c)
{
    const double pixel_km = terrain::legacycellkm(world);
    double total = 0.0;

    for (int slot = 1; slot <= mesh::slotcount; slot++)
        total += castfetch_km(world, c, slot, castlimitpixels * pixel_km) / pixel_km;

    total = total * rangeperpixel;

    total = total * world.lunar(); // Lunar pull multiplies tides.

    int tide = (int)total; // In theory this could be as high as 80, but never in practice assuming lunar = 1.0. The highest in the world is 16 metres, so divide this value by 4 to get it in metres.

    return (tide);
}
}

// This function calculates the tides: each sea cell next to land gets its range, and its
// neighbours without one take it too.

void createtidalmap(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int lunar = (int)world.lunar();

    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        world.settide(c, -1);

    for (const mesh::CellId c : hydrologycells::visitorder(grid))
    {
        if (world.sea(c) == 1)
        {
            bool landfound = 0;

            for (const mesh::CellId cell : hydrologycells::square(grid, c))
            {
                if (world.sea(cell) == 0)
                {
                    landfound = 1;
                    break;
                }
            }

            if (landfound == 1)
            {
                int tide = gettidalrange(world, c);

                for (const mesh::CellId cell : hydrologycells::square(grid, c))
                {
                    if (world.tide(cell) == -1)
                        world.settide(cell, tide);
                }

                if (tide < lunar * 2)
                    tide = lunar * 2;

                world.settide(c, tide);
            }
        }
    }
}
