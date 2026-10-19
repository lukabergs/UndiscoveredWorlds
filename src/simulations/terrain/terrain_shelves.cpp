// Continental shelf masks and their connected gaps.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the sphere since T8 W2-T2
// (docs/architecture/sphere-mesh.md sections 3 and 11).

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_shapes.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;
using namespace terrain::coastgrid;
using mesh::CellId;

namespace
{
void removeshelfgaps(planet& world, worldgrid<bool>& shelves);

// The shelves on the sphere: the legacy construction in legacy pixels of km
// (tuning::terrain::sphere_noise_coast): the warp and the circle radii from sphere noise, the
// joining curve (a straight line in the legacy, its control point the midpoint) as a great-circle
// arc, the panel test within 50 pixels as a disc, the circles' chance per pixel area.
void spherecontinentalshelves(planet& world, worldgrid<bool>& shelves, int pointdist)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = terrain::legacycellkm(world);
    int maxelev = world.maxelevation();
    int midelev = maxelev / 2;

    vector<char> outline(count, 0);

    const float maximum = static_cast<float>(maxelev);
    const vector<float> fractal = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x5e1f0001), terrain::terrainfractal(), maximum);
    const vector<float> fractalx = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x5e1f0002), terrain::terrainfractal(), maximum);
    const vector<float> fractaly = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x5e1f0003), terrain::terrainfractal(), maximum);

    int maxwarp = 60;
    int warpdiv = max(1, maxelev / maxwarp);
    int maxradius = 20;
    int raddiv = max(1, maxelev / maxradius);

    const vector<CellId> order = columnmajorcells(grid);

    for (const CellId c : order)
    {
        if (!terrain::coastgrid::outline(world, c))
            continue;

        // The centre of this circle, offset to make the shelves more interesting (east and south in pixels).
        const int warpdistx = (static_cast<int>(fractalx[c]) - midelev) / warpdiv;
        const int warpdisty = (static_cast<int>(fractaly[c]) - midelev) / warpdiv;
        const CellId target = terrain::movedcell(world, c, warpdistx * pixelkm, -warpdisty * pixelkm);

        // A line from the outline to that centre, so that there is no gap in the shelf.
        for (const CellId n : mesh::cellsalongarc(grid, grid.centre(c), grid.centre(target)))
            outline[n] = 1;

        // Now we can draw a circle.
        const int radius = static_cast<int>(fractal[c]) / raddiv;
        if (radius > 0)
            for (const CellId n : terrain::cellswithinkm(world, target, (radius + 0.5) * pixelkm))
                outline[n] = 1;
    }

    // Now, every panel of the voronoi map that has any outline within 50 pixels is continental shelf.

    vector<int> voronoi;

    makeshelvesvoronoi(world, voronoi, pointdist);

    int panels = 0;
    for (const int id : voronoi)
        panels = max(panels, id);
    vector<int> start(static_cast<size_t>(panels) + 2, 0);
    for (const int id : voronoi)
        start[static_cast<size_t>(id) + 1]++;
    for (size_t k = 1; k < start.size(); k++)
        start[k] += start[k - 1];
    vector<CellId> members(count);
    {
        vector<int> fillat(start.begin(), start.end() - 1);
        for (CellId c = 0; c < count; c++)
            members[static_cast<size_t>(fillat[static_cast<size_t>(voronoi[c])]++)] = c;
    }

    const double panelreach_m = 50.5 * pixelkm * 1000.0;

    for (const CellId c : order)
    {
        if (outline[c] == 1 && shelves[c] == 0)
        {
            const int id = voronoi[c];
            const mesh::Vec3 p = grid.centre(c);
            for (int k = start[static_cast<size_t>(id)]; k < start[static_cast<size_t>(id) + 1]; k++)
            {
                const CellId n = members[static_cast<size_t>(k)];
                if (mesh::greatcircle_m(p, grid.centre(n), grid.radius_m()) <= panelreach_m)
                    shelves[n] = 1;
            }
        }
    }

    // And anywhere that's land is continental shelf too, technically.

    for (CellId c = 0; c < count; c++)
    {
        if (world.sea(c) == 0)
            shelves[c] = 1;
    }

    // Now add a circles of shelf to the shelves, to make the edges a bit more varied.

    vector<char> circles(count, 0);

    const double circlechance = 1.0 / 40.0; // per legacy pixel
    const double mindist = 5.0; // Minimum distance to non-shelf sea (pixels).
    const double radius = 8.0;

    for (const CellId c : order)
    {
        if (shelves[c] == 1 && uniformdraw() < cellchance(world, c, circlechance))
        {
            bool nearsea = 0;

            for (const CellId n : terrain::cellswithinkm(world, c, (mindist + 0.5) * pixelkm))
            {
                if (shelves[n] == 0)
                {
                    nearsea = 1;
                    break;
                }
            }

            if (nearsea == 0)
                for (const CellId n : terrain::cellswithinkm(world, c, (radius + 0.5) * pixelkm))
                    circles[n] = 1;
        }
    }

    for (CellId c = 0; c < count; c++)
    {
        if (circles[c] == 1)
            shelves[c] = 1;
    }
}
}

// This function makes continental shelves.

void makecontinentalshelves(planet& world, worldgrid<bool>& shelves, int pointdist)
{
    spherecontinentalshelves(world, shelves, pointdist);

    // Now we need to remove any holes in the continental shelves!

    removeshelfgaps(world, shelves);
}

namespace
{
// This function removes any gaps in the continental shelves: every non-shelf cell outside the
// largest 4-connected area of non-shelf (by area) becomes shelf.

void removeshelfgaps(planet& world, worldgrid<bool>& shelves)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();

    vector<double> area(count, 0.0);
    vector<char> checked(count, 0);
    vector<CellId> component;
    const vector<CellId> order = columnmajorcells(grid);

    auto enterable = [&](CellId c) { return shelves[c] == 0; };

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const CellId first = *it;
        if (shelves[first] != 0 || checked[first] || area[first] != 0.0)
            continue;

        component.clear();
        queue<CellId> q;
        q.push(first);
        checked[first] = 1;
        double size = 0.0;
        while (!q.empty())
        {
            const CellId c = q.front();
            q.pop();
            component.push_back(c);
            size += grid.area_m2(c);
            for (int slot = 1; slot <= mesh::slotcount; slot += 2)
            {
                const std::int32_t n = grid.neighbour(c, slot);
                if (n >= 0 && !checked[n] && enterable(static_cast<CellId>(n)))
                {
                    checked[n] = 1;
                    q.push(static_cast<CellId>(n));
                }
            }
        }
        for (const CellId c : component)
            area[c] = size;
    }

    // Now we just find the largest area of non-shelf.

    double largest = 0.0;

    for (const CellId c : order)
    {
        if (area[c] > largest)
            largest = area[c];
    }

    // Now we simply remove any sea tiles that aren't part of the largest area of sea.

    for (CellId c = 0; c < count; c++)
    {
        if (area[c] != largest)
            shelves[c] = 1;
    }
}
}
