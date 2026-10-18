// Continental shelf masks and their connected gaps.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the cell grid since T8 W2-T2
// (docs/architecture/sphere-mesh.md sections 3 and 11).

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;
using namespace terrain::coastgrid;
using mesh::CellId;

namespace
{
void removeshelfgaps(planet& world, worldgrid<bool>& shelves);

// The legacy shelves (lat-lon): outline circles at fractal-warped centres joined to the coast
// by a curve, the Voronoi panels that hold outline, random circles, all in pixels.
void latloncontinentalshelves(planet& world, worldgrid<bool>& shelves, int pointdist)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int midelev = maxelev / 2;
    const mesh::CellGrid& grid = world.grid();

    // First, sort out the outline.

    vector<vector<bool>> outline(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f; //0.4f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    int maxwarp = 60;
    int warpdiv = maxelev / maxwarp;
    int maxradius = 20;
    int raddiv = maxelev / maxradius;

    twofloats pt, mm1, mm2, mm3;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.outline(i, j) == 1)
            {
                // First, work out where the centre of this circle will be. It's getting offset to make the shelves more interesting.

                int warpdistx = 0;
                int warpdisty = 0;

                int warpx = i + width / 4;

                if (warpx > width)
                    warpx = warpx - width;

                int warpy = j;

                if (fractal[warpx][warpy] > maxelev / 2)
                    warpdistx = (fractal[warpx][warpy] - midelev) / warpdiv;
                else
                    warpdistx = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

                int ii = i + warpdistx;

                warpx = warpx + width / 2;

                if (warpx > width)
                    warpx = warpx - width;

                if (fractal[warpx][warpy] > maxelev / 2)
                    warpdisty = (fractal[warpx][warpy] - midelev) / warpdiv;
                else
                    warpdisty = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

                int jj = j + warpdisty;

                // ii and jj are the centre of the circle we're drawing. But we want to draw a line from the actual outline to that centre, to ensure that there isn't a gap in the shelf.

                mm1.x = (float) i;
                mm1.y = (float) j;

                mm2.x = (float) (i + ii) / 2;
                mm2.y = (float) (j + jj) / 2;

                mm3.x = (float) ii;
                mm3.y = (float) jj;

                for (int n = 1; n <= 2; n++) // Two halves of the curve.
                {
                    for (float t = 0.0f; t <= 1.0f; t = t + 0.01f)
                    {
                        if (n == 1)
                            pt = curvepos(mm1, mm1, mm2, mm3, t);
                        else
                            pt = curvepos(mm1, mm2, mm3, mm3, t);

                        int x = (int)pt.x;
                        int y = (int)pt.y;

                        if (x<0 || x>width)
                            x = wrap(x, width);

                        if (y < 0)
                            y = 0;

                        if (y > height)
                            y = height;

                        outline[x][y] = 1;
                    }
                }

                // Now we can draw a circle.

                if (jj < 0)
                    jj = 0;

                if (jj > height)
                    jj = height;

                int radius = fractal[i][j] / raddiv;

                for (int x = -radius; x <= radius; x++)
                {
                    int xx = ii + x;

                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    for (int y = -radius; y <= radius; y++)
                    {
                        int yy = jj + y;

                        if (yy >= 0 && yy <= height)
                        {
                            if (x * x + y * y < radius * radius + radius)
                                outline[xx][yy] = 1;
                        }
                    }
                }
            }
        }
    }

    // Now, we need a voronoi map.

    vector<int> voronoi;

    makeshelvesvoronoi(world, voronoi, pointdist);

    // Now, every panel of the voronoi map that has any outline in it is continental shelf.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (outline[i][j] == 1 && shelves[i][j] == 0)
            {
                int thiscell = voronoi[grid.cellof(i, j)];

                for (int k = i - 50; k <= i + 50; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 50; l <= j + 50; l++)
                    {
                        if (l >= 0 && l < height)
                        {
                            if (voronoi[grid.cellof(kk, l)] == thiscell)
                                shelves[kk][l] = 1;
                        }
                    }
                }
            }
        }
    }

    // And anywhere that's land is continental shelf too, technically.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (world.sea(i, j) == 0)
                shelves[i][j] = 1;
        }
    }

    // Now add a circles of shelf to the shelves, to make the edges a bit more varied.

    vector<vector<bool>> circles(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int circlechance = 40;
    int mindist = 5; // Minimum distance to non-shelf sea.
    int radius = 8;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (shelves[i][j] == 1 && random(1, circlechance) == 1)
            {
                bool nearsea = 0;

                for (int k = -mindist; k <= mindist; k++)
                {
                    for (int l = -mindist; l <= mindist; l++)
                    {
                        int ii = i + k;
                        int jj = j + l;

                        if (ii<0 || ii>width)
                            ii = wrap(ii, width);

                        if (jj >= 0 && jj <= height)
                        {
                            if (shelves[ii][jj] == 0)
                            {
                                nearsea = 1;
                                k = mindist;
                                l = mindist;
                            }
                        }
                    }
                }

                if (nearsea == 0)
                {

                    for (int x = -radius; x <= radius; x++)
                    {
                        int xx = i + x;

                        if (xx<0 || xx>width)
                            xx = wrap(xx, width);

                        for (int y = -radius; y <= radius; y++)
                        {
                            int yy = j + y;

                            if (yy >= 0 && yy <= height)
                            {
                                if (x * x + y * y < radius * radius + radius)
                                    circles[xx][yy] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (circles[i][j] == 1)
                shelves[i][j] = 1;
        }
    }
}

// The shelves on the sphere: the same construction in legacy pixels of km
// (tuning::terrain::sphere_noise_coast): the warp and the circle radii from sphere noise, the
// joining curve (a straight line in the legacy, its control point the midpoint) as a great-circle
// arc, the panel test within 50 pixels as a disc, the circles' chance per pixel area.
void spherecontinentalshelves(planet& world, worldgrid<bool>& shelves, int pointdist)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = legacypixelkm(world);
    int maxelev = world.maxelevation();
    int midelev = maxelev / 2;

    vector<char> outline(count, 0);

    const vector<int> fractal = spherefractal(world, 0x5e1f0001, grain8noise(), maxelev);
    const vector<int> fractalx = spherefractal(world, 0x5e1f0002, grain8noise(), maxelev);
    const vector<int> fractaly = spherefractal(world, 0x5e1f0003, grain8noise(), maxelev);

    int maxwarp = 60;
    int warpdiv = max(1, maxelev / maxwarp);
    int maxradius = 20;
    int raddiv = max(1, maxelev / maxradius);

    const vector<CellId> order = columnorder(grid);

    for (const CellId c : order)
    {
        if (!terrain::coastgrid::outline(world, c))
            continue;

        // The centre of this circle, offset to make the shelves more interesting (east and south in pixels).
        const int warpdistx = (fractalx[c] - midelev) / warpdiv;
        const int warpdisty = (fractaly[c] - midelev) / warpdiv;
        const mesh::Vec3 target = movedpoint(world, c, warpdistx * pixelkm, -warpdisty * pixelkm);

        // A line from the outline to that centre, so that there is no gap in the shelf.
        for (const CellId n : mesh::cellsalongarc(grid, grid.centre(c), target))
            outline[n] = 1;

        // Now we can draw a circle.
        const int radius = fractal[c] / raddiv;
        if (radius > 0)
            for (const CellId n : cellswithinkm(world, target, (radius + 0.5) * pixelkm))
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

            for (const CellId n : cellswithinkm(world, grid.centre(c), (mindist + 0.5) * pixelkm))
            {
                if (shelves[n] == 0)
                {
                    nearsea = 1;
                    break;
                }
            }

            if (nearsea == 0)
                for (const CellId n : cellswithinkm(world, grid.centre(c), (radius + 0.5) * pixelkm))
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
    if (world.grid().kind() == mesh::GridKind::latlon)
        latloncontinentalshelves(world, shelves, pointdist);
    else
        spherecontinentalshelves(world, shelves, pointdist);

    // Now we need to remove any holes in the continental shelves!

    removeshelfgaps(world, shelves);
}

namespace
{
// This function removes any gaps in the continental shelves: every non-shelf cell outside the
// largest 4-connected area of non-shelf becomes shelf. On lat-lon the legacy flood: started from
// interior cells, never entering the last row, the largest by cell count.

void removeshelfgaps(planet& world, worldgrid<bool>& shelves)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const bool legacy = grid.kind() == mesh::GridKind::latlon;

    vector<double> area(count, 0.0);
    vector<char> checked(count, 0);
    vector<CellId> component;
    const vector<CellId> order = columnorder(grid);

    auto enterable = [&](CellId c)
    {
        return shelves[c] == 0 && (!legacy || grid.neighbour(c, 5) >= 0);
    };

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const CellId first = *it;
        if (shelves[first] != 0 || checked[first] || area[first] != 0.0)
            continue;
        if (legacy)
        {
            const int x = grid.column(first);
            const int y = grid.row(first);
            if (x < 1 || x > grid.columns() - 2 || y < 1 || y > grid.rows() - 2)
                continue;
        }

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
            size += legacy ? 1.0 : grid.area_m2(c);
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
