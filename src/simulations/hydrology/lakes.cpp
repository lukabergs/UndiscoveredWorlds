#include <vector>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"
#include "simulations/hydrology/hydrology_cells.hpp"
#include "simulations/hydrology/lakes.hpp"

using namespace std;

namespace
{
// This checks whether a given tile is on the edge of a given lake: a neighbour (any of the
// eight) outside it. thislake is the caller's lat-lon array, read at the neighbours' columns
// and rows.

int lakeoutline(planet& world, vector<vector<int>>& thislake, int lakeno, int x, int y)
{
    const mesh::CellGrid& grid = world.grid();
    const mesh::CellId c = grid.cellof(x, y);

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        const mesh::CellId n = static_cast<mesh::CellId>(neighbour);

        if (thislake[grid.column(n)][grid.row(n)] != lakeno)
            return(1);
    }

    return (0);
}
}

// This function tells us whether a cell is next to a lake: the first lake surface of its
// 3 x 3 neighbourhood (in the legacy order).

float nexttolake(planet& world, mesh::CellId c)
{
    for (const mesh::CellId candidate : hydrologycells::square(world.grid(), c))
    {
        if (world.lakesurfacefloat(candidate) != 0)
            return world.lakesurfacefloat(candidate);
    }
    return 0.0f;
}

float nexttolake(planet& world, int x, int y)
{
    const mesh::CellGrid& grid = world.grid();
    if (y < 0 || y >= grid.rows())
        return 0.0f;
    const int columns = grid.columns();
    return nexttolake(world, grid.cellof(((x % columns) + columns) % columns, y));
}

// Select the river source cell for a lake or terrain water feature. Lat-lon worlds only: the
// lake is the caller's (x, y) array and box.
void makelakestartpoint(planet& world, vector<vector<int>>& thislake, int lakeno, int leftx, int lefty, int rightx, int righty)
{
    twointegers startpoint;
    startpoint.x = -1;
    startpoint.y = -1;

    int crount = 0;

    do
    {
        for (int i = leftx; i <= rightx; i++)
        {
            for (int j = lefty; j <= righty; j++)
            {
                if (thislake[i][j] == lakeno)
                {
                    if (lakeoutline(world, thislake, lakeno, i, j) == 1) // If this is on the edge of the lake
                    {
                        if (random(1, 3) == 1)
                        {
                            startpoint.x = i;
                            startpoint.y = j;
                        }
                    }
                }
            }
        }

        crount++;

        if (crount > 1000) // Something's gone horribly wrong
            return;

    } while (startpoint.x == -1);

    world.setlakestart(startpoint.x, startpoint.y, 1);
}
