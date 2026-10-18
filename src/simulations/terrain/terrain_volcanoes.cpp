// Isolated volcanoes and mountainous islands.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the cell grid
// (docs/architecture/sphere-mesh.md sections 3 and 12) since T8 W2-T1.

#include "simulations/terrain/terrain_internal.hpp"
#include "simulations/terrain/terrain_volcanoes.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_ridges.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;

namespace
{
// The neighbour slot of a step of dx columns east and dy rows south (each -1, 0 or 1; not both 0).
int stepslot(int dx, int dy)
{
    if (dy < 0)
        return dx < 0 ? 8 : (dx == 0 ? 1 : 2);
    if (dy > 0)
        return dx < 0 ? 6 : (dx == 0 ? 5 : 4);
    return dx < 0 ? 7 : 3;
}

// Marks a cell and its eight neighbours.
template <typename Fn>
void forbox(const planet& world, mesh::CellId c, Fn&& fn)
{
    for (const int slot : boxslots)
    {
        const std::int32_t n = slot == 0 ? static_cast<std::int32_t>(c) : world.grid().neighbour(c, slot);
        if (n >= 0)
            fn(static_cast<mesh::CellId>(n));
    }
}

// Whether a continental-shelf cell lies within the shelf check of five legacy pixels.
bool nearshelf(const planet& world, mesh::CellId c, const worldgrid<bool>& shelves)
{
    const int shelfcheck = 5;
    for (const mesh::CellId n : terrain::cellswithinkm(world, c, (shelfcheck + 0.5) * terrain::legacycellkm(world)))
    {
        if (shelves[n] == 1)
            return true;
    }
    return false;
}
}

// This makes a small, mountainous island.

void terrain::detail::makemountainisland(planet& world, int x, int y, float peakheight)
{
    const mesh::CellId c = world.grid().cellof(x, y);
    if (c != mesh::nocell)
        terrain::makemountainisland(world, c, peakheight);
}

void terrain::makemountainisland(planet& world, mesh::CellId c, float peakheight)
{
    int sealevel = world.sealevel();

    int iterations = random(1, 3);

    for (int n = 1; n <= iterations; n++)
    {
        int landheight = sealevel + random(5, 50);

        int dx = 0;
        int dy = 0;

        do
        {
            dx = -2 + random(1, 3);
            dy = -2 + random(1, 3);

        } while (dx == 0 && dy == 0);

        const int dir = stepslot(dx, dy);
        const std::int32_t neighbour = world.grid().neighbour(c, dir);

        if (neighbour < 0)
            return;

        const mesh::CellId next = static_cast<mesh::CellId>(neighbour);

        if (world.mountainheightfloat(c) < peakheight)
            world.setmountainheight(c, peakheight);

        if (world.mountainheightfloat(next) < peakheight)
            world.setmountainheight(next, peakheight);

        // Mark it all out as mountain islands!
        forbox(world, c, [&](mesh::CellId m) { world.setmountainisland(m, 1); });
        forbox(world, next, [&](mesh::CellId m) { world.setmountainisland(m, 1); });

        if (world.nomfloat(c) < landheight)
            world.setnomfloat(c, static_cast<float>(landheight));

        if (world.nomfloat(next) < landheight)
            world.setnomfloat(next, static_cast<float>(landheight));

        // If there isn't a ridge between them, make one. The legacy code gives the neighbour the
        // same direction code (dir + 8, less 8), not the reverse; kept as it was.

        if (getridge(world, c, dir) == 0)
        {
            const int code = terrain::detail::getcode(dir);

            world.setmountainridge(c, world.mountainridge(c) + code);
            world.setmountainridge(next, world.mountainridge(next) + code);
        }

        c = next;

        peakheight = peakheight + randomsign(random((int)peakheight / 4, (int)peakheight / 2));
    }
}

// This makes an isolated volcano, together with a line of extinct neighbours.

void createisolatedvolcano(planet& world, int x, int y, worldgrid<bool>& shelves, vector<vector<int>>& volcanodirection, int peakheight, bool strato)
{
    // The legacy direction grid: column x for the east shift, column x + width / 2 for the north shift.
    const mesh::CellGrid& grid = world.grid();
    const mesh::CellId c = grid.cellof(x, y);
    if (c == mesh::nocell)
        return;
    const int width = world.width();
    vector<float> eastdirection(world.cellcount()), northdirection(world.cellcount());
    for (mesh::CellId cell = 0; cell < world.cellcount(); cell++)
    {
        const int column = grid.column(cell);
        const int row = grid.row(cell);
        int other = column + width / 2;
        if (other > width)
            other = other - width;
        eastdirection[cell] = static_cast<float>(volcanodirection[column][row]);
        northdirection[cell] = static_cast<float>(volcanodirection[other][row]);
    }
    terrain::createisolatedvolcano(world, c, shelves, eastdirection, northdirection, peakheight, strato);
}

void terrain::createisolatedvolcano(planet& world, mesh::CellId c, const worldgrid<bool>& shelves,
                                    const vector<float>& eastdirection, const vector<float>& northdirection,
                                    int peakheight, bool strato)
{
    int maxelev = world.maxelevation();

    float maxshift = 6;
    float shiftdiv = maxshift / (maxelev / 2);

    const double pixel_km = terrain::legacycellkm(world);

    if (world.sea(c) == 1 && nearshelf(world, c, shelves)) // Don't do any submarine volcanoes near continental shelves.
        return;

    bool active = 1;

    int total = random(1, 6);

    for (int n = 1; n != total; n++)
    {
        if (world.sea(c) == 1) // No shading around undersea volcanoes.
            forbox(world, c, [&](mesh::CellId m) { world.setnoshade(m, 1); });

        float thispeakheight = (float)peakheight;

        thispeakheight = thispeakheight / 100.0f;

        thispeakheight = thispeakheight * (float)random(80, 120);

        if (active == 0)
            thispeakheight = (float)(0 - peakheight);

        // Remove any volcanoes immediately around.
        forbox(world, c, [&](mesh::CellId m) { world.setvolcano(m, 0); });

        world.setvolcano(c, thispeakheight);
        world.setstrato(c, strato);

        active = 0;

        float xshift = eastdirection[c] - (float)(maxelev / 2);
        float yshift = northdirection[c] - (float)(maxelev / 2);

        xshift = xshift * shiftdiv + (float)randomsign(random(0, 2));
        yshift = yshift * shiftdiv + (float)randomsign(random(0, 2));

        // Whole legacy pixels east and south, as the legacy (x, y) step.
        c = terrain::movedcell(world, c, static_cast<double>((int)xshift) * pixel_km, -static_cast<double>((int)yshift) * pixel_km);

        if (world.sea(c) == 1 && nearshelf(world, c, shelves)) // Don't do any submarine volcanoes near continental shelves.
            return;
    }
}
