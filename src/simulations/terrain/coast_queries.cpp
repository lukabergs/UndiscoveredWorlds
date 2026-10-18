// Questions about coasts on the world grid: distance to land and sea, coastal cells, land in one quadrant.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022); on the cell grid since T8 W2-T2.

#include "simulations/terrain/coast_queries.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <utility>

#include "core/classes.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"

using namespace std;
using namespace terrain::coastgrid;
using mesh::CellId;

namespace
{
// planet::sea at a neighbour: out of the grid (-1) reads as land, as sea(x, y) does.
bool seaat(const planet& world, std::int32_t n)
{
    return n >= 0 && world.sea(static_cast<CellId>(n));
}

// The lat-lon cell of (x, y) with x wrapped, -1 outside the rows (or on another grid kind).
std::int32_t latloncell(const planet& world, int x, int y)
{
    const mesh::CellGrid& grid = world.grid();
    if (grid.kind() != mesh::GridKind::latlon || y < 0 || y >= grid.rows())
        return -1;
    x = wrap(x, grid.columns() - 1);
    const CellId c = grid.cellof(x, y);
    return c == mesh::nocell ? -1 : static_cast<std::int32_t>(c);
}

// The cardinal slot whose bearing is nearest the target bearing (0 north, pi south).
int slottoward(const mesh::CellGrid& grid, CellId c, double target)
{
    int best = 0;
    double bestdifference = 10.0;
    for (int slot = 1; slot <= mesh::slotcount; slot += 2)
    {
        if (grid.neighbour(c, slot) < 0)
            continue;
        double difference = fabs(grid.bearing_rad(c, slot) - target);
        difference = min(difference, 2.0 * 3.14159265358979323846 - difference);
        if (difference < bestdifference)
        {
            bestdifference = difference;
            best = slot;
        }
    }
    return best;
}
}

// This finds the distance to the nearest land in the direction of the equator.

int terrain::coastgrid::landdistance(const planet& world, CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    const bool north = grid.latitude_rad(c) > 0.0;

    // Reverse the direction of search if the equator is colder than the poles!
    bool towardequator = true;
    if (north && world.northpolartemp() > world.eqtemp())
        towardequator = false;
    if (!north && world.southpolartemp() > world.eqtemp())
        towardequator = false;
    const double target = (north == towardequator) ? 3.14159265358979323846 : 0.0;

    const int limit = static_cast<int>(3.14159265358979323846 * grid.radius_m() / max(1.0, grid.minimumspacing_m())) + 2;
    std::int32_t at = static_cast<std::int32_t>(c);

    for (int dist = 0; dist < limit; dist++)
    {
        const CellId here = static_cast<CellId>(at);
        if (towardequator && (grid.latitude_rad(here) > 0.0) != north)
            return -1;

        const int slot = slottoward(grid, here, target);
        if (slot == 0)
            return -1;
        int turn = 0;
        const std::int32_t next = step(grid, here, slot, turn);

        if (world.sea(here) == 0)
        {
            int turn2 = turn;
            const std::int32_t after = next >= 0 ? step(grid, static_cast<CellId>(next), slot, turn2) : -1;
            if (!seaat(world, next) && !seaat(world, after)) // To avoid little islands showing up
            {
                const std::int32_t left = next >= 0 ? grid.neighbour(static_cast<CellId>(next), turned(slot, turn - 2)) : -1;
                const std::int32_t right = next >= 0 ? grid.neighbour(static_cast<CellId>(next), turned(slot, turn + 2)) : -1;
                if (!seaat(world, left) && !seaat(world, right))
                    return dist;
            }
        }

        if (next < 0)
            return -1;
        at = next;
    }

    return -1; // No land found at all.
}

int landdistance(planet& world, int x, int y)
{
    const std::int32_t c = latloncell(world, x, y);
    if (c < 0 || x < 0 || x > world.width())
        return -1;
    return terrain::coastgrid::landdistance(world, static_cast<CellId>(c));
}

// This finds the nearest sea (or land, if land=1) to the cell specified.

std::int32_t terrain::coastgrid::nearestsea(const planet& world, CellId c, bool land, double limit_km)
{
    const mesh::CellGrid& grid = world.grid();
    const mesh::Vec3 p = grid.centre(c);
    std::int32_t best = -1;
    double bestdistance = numeric_limits<double>::max();
    for (const CellId n : cellswithinkm(world, p, limit_km))
    {
        if (n == c || (world.sea(n) == 1) == land)
            continue;
        const double distance = mesh::greatcircle_rad(p, grid.centre(n));
        if (distance < bestdistance)
        {
            bestdistance = distance;
            best = static_cast<std::int32_t>(n);
        }
    }
    return best;
}

twointegers nearestsea(planet& world, int i, int j, bool land, int limit, int grain)
{
    (void)grain;
    twointegers seapoint;
    seapoint.x = -1;
    seapoint.y = -1;

    const std::int32_t c = latloncell(world, i, j);
    if (c < 0)
        return seapoint;
    const mesh::CellGrid& grid = world.grid();
    const double rowkm = 3.14159265358979323846 * grid.radius_m() / 1000.0 / static_cast<double>(grid.rows());
    const std::int32_t found = terrain::coastgrid::nearestsea(world, static_cast<CellId>(c), land, static_cast<double>(limit) * rowkm);
    if (found >= 0)
    {
        seapoint.x = grid.column(static_cast<CellId>(found));
        seapoint.y = grid.row(static_cast<CellId>(found));
    }
    return seapoint;
}

// This function checks whether a global tile is either sea next to land or land next to sea.

bool terrain::coastgrid::vaguelycoastal(const planet& world, CellId c)
{
    const bool sea = world.sea(c) == 1;

    bool found = 0;

    box(world.grid(), c, 1, [&](CellId n, int, int)
    {
        if ((world.sea(n) == 1) != sea)
            found = 1;
        return found == 0;
    });

    return found;
}

bool vaguelycoastal(planet& world, int x, int y)
{
    const std::int32_t c = latloncell(world, x, y);
    if (c < 0 || x < 0 || x > world.width())
        return false;
    return terrain::coastgrid::vaguelycoastal(world, static_cast<CellId>(c));
}

// These check whether a tile is land bordered by sea on three sides of one corner, with land
// at the opposite corner. On lat-lon the legacy forms read the opposite corner of the
// northwest and southwest cases without wrapping x (land beyond the column edge).

bool terrain::coastgrid::northwestlandonly(const planet& world, std::int32_t c)
{
    if (c < 0)
        return 0;
    const mesh::CellGrid& grid = world.grid();
    const CellId cell = static_cast<CellId>(c);

    if (world.sea(cell) == 1)
        return 0;

    if (grid.neighbour(cell, 5) < 0 || grid.neighbour(cell, 1) < 0)
        return 0;

    if (!seaat(world, grid.neighbour(cell, 3)) || !seaat(world, grid.neighbour(cell, 5)) || !seaat(world, grid.neighbour(cell, 4)))
        return 0;

    if (seaat(world, rawneighbour(grid, cell, 8)))
        return 0;

    return 1;
}

bool terrain::coastgrid::northeastlandonly(const planet& world, std::int32_t c)
{
    if (c < 0)
        return 0;
    const mesh::CellGrid& grid = world.grid();
    const CellId cell = static_cast<CellId>(c);

    if (world.sea(cell) == 1)
        return 0;

    if (grid.neighbour(cell, 5) < 0 || grid.neighbour(cell, 1) < 0)
        return 0;

    if (!seaat(world, grid.neighbour(cell, 7)) || !seaat(world, grid.neighbour(cell, 5)) || !seaat(world, grid.neighbour(cell, 6)))
        return 0;

    if (seaat(world, grid.neighbour(cell, 2)))
        return 0;

    return 1;
}

bool terrain::coastgrid::southwestlandonly(const planet& world, std::int32_t c)
{
    if (c < 0)
        return 0;
    const mesh::CellGrid& grid = world.grid();
    const CellId cell = static_cast<CellId>(c);

    if (world.sea(cell) == 1)
        return 0;

    if (grid.neighbour(cell, 1) < 0 || grid.neighbour(cell, 5) < 0)
        return 0;

    if (!seaat(world, grid.neighbour(cell, 3)) || !seaat(world, grid.neighbour(cell, 1)) || !seaat(world, grid.neighbour(cell, 2)))
        return 0;

    if (seaat(world, rawneighbour(grid, cell, 4)))
        return 0;

    return 1;
}

bool terrain::coastgrid::southeastlandonly(const planet& world, std::int32_t c)
{
    if (c < 0)
        return 0;
    const mesh::CellGrid& grid = world.grid();
    const CellId cell = static_cast<CellId>(c);

    if (world.sea(cell) == 1)
        return 0;

    if (grid.neighbour(cell, 1) < 0 || grid.neighbour(cell, 5) < 0)
        return 0;

    if (!seaat(world, grid.neighbour(cell, 7)) || !seaat(world, grid.neighbour(cell, 1)) || !seaat(world, grid.neighbour(cell, 8)))
        return 0;

    if (seaat(world, grid.neighbour(cell, 4)))
        return 0;

    return 1;
}

bool northwestlandonly(planet& world, int x, int y)
{
    return terrain::coastgrid::northwestlandonly(world, latloncell(world, x, y));
}

bool northeastlandonly(planet& world, int x, int y)
{
    return terrain::coastgrid::northeastlandonly(world, latloncell(world, x, y));
}

bool southwestlandonly(planet& world, int x, int y)
{
    return terrain::coastgrid::southwestlandonly(world, latloncell(world, x, y));
}

bool southeastlandonly(planet& world, int x, int y)
{
    return terrain::coastgrid::southeastlandonly(world, latloncell(world, x, y));
}

namespace
{
// The legacy lattice Voronoi (lat-lon): points every pointdist pixels, moved by up to three
// quarters of that, and each cell's nearest point within growing square searches.
void latlonshelvesvoronoi(planet& world, vector<int>& voronoi, int pointdist)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    int grain = 16; // Level of detail on this fractal map.
    float valuemod = 0.4f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0); // Its warp is not applied; the call keeps the draws.

    vector<vector<int>> pointsmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int pointnumber = 1;
    int movedist = (pointdist / 4) * 3;

    for (int i = 0; i <= width; i = i + pointdist)
    {
        for (int j = 0; j <= height; j = j + pointdist)
        {
            int ii = i + randomsign(random(1, movedist));
            int jj = j + randomsign(random(1, movedist));

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            if (jj < 0)
                jj = 0;

            if (jj > height)
                jj = height;

            pointsmap[ii][jj] = pointnumber;
            pointnumber++;
        }
    }

    const mesh::CellGrid& grid = world.grid();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int closest = 0;
            int dist = 10000000;

            int checkx = i;
            int checky = j;

            int distcheck = 10;

            do
            {
                for (int x = checkx - distcheck; x <= checkx + distcheck; x++)
                {
                    int xx = x;

                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    for (int y = checky - distcheck; y <= checky + distcheck; y++)
                    {
                        if (y >= 0 && y <= height)
                        {
                            if (pointsmap[xx][y] != 0)
                            {
                                int xdist = checkx - xx;
                                int ydist = checky - y;
                                int pointdistance = xdist * xdist + ydist * ydist;

                                if (pointdistance < dist)
                                {
                                    closest = pointsmap[xx][y];
                                    dist = pointdistance;
                                }
                            }
                        }
                    }
                }

                distcheck = distcheck + 10;

            } while (closest == 0);

            voronoi[grid.cellof(i, j)] = static_cast<short>(closest); // the legacy short panels
        }
    }
}

// The Voronoi panels on the sphere: seeds drawn by area, one per (pointdist legacy pixels)^2,
// and each cell the seed nearest its centre (a Dijkstra-ordered propagation over the neighbour
// table with the great-circle distance to the candidate seed; ties to the earlier seed).
void sphereshelvesvoronoi(planet& world, vector<int>& voronoi, int pointdist)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double spacing_m = static_cast<double>(max(1, pointdist)) * legacypixelkm(world) * 1000.0;
    const int seeds = max(1, static_cast<int>(lround(grid.totalarea_m2() / (spacing_m * spacing_m))));

    vector<mesh::Vec3> points;
    points.reserve(static_cast<size_t>(seeds) + 1);
    points.push_back(mesh::Vec3{});
    for (int k = 1; k <= seeds; k++)
        points.push_back(randompoint());

    vector<double> best(count, numeric_limits<double>::max());
    fill(voronoi.begin(), voronoi.end(), 0);
    for (int k = 1; k <= seeds; k++)
    {
        const CellId c = grid.locate(points[static_cast<size_t>(k)]);
        if (c == mesh::nocell)
            continue;
        const double d = mesh::greatcircle_rad(grid.centre(c), points[static_cast<size_t>(k)]);
        if (d < best[c])
        {
            best[c] = d;
            voronoi[c] = k;
        }
    }

    using Entry = pair<double, CellId>;
    priority_queue<Entry, vector<Entry>, greater<Entry>> q;
    for (CellId c = 0; c < count; c++)
        if (voronoi[c] != 0)
            q.push({ best[c], c });

    while (!q.empty())
    {
        const Entry top = q.top();
        q.pop();
        if (top.first > best[top.second])
            continue;
        const int seed = voronoi[top.second];
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t n = grid.neighbour(top.second, slot);
            if (n < 0)
                continue;
            const double d = mesh::greatcircle_rad(grid.centre(static_cast<CellId>(n)), points[static_cast<size_t>(seed)]);
            if (d < best[n] || (d == best[n] && seed < voronoi[n]))
            {
                best[n] = d;
                voronoi[n] = seed;
                q.push({ d, static_cast<CellId>(n) });
            }
        }
    }
}
}

// This does the same thing, for the continental shelves.

void makeshelvesvoronoi(planet& world, vector<int>& voronoi, int pointdist)
{
    voronoi.assign(world.cellcount(), 0);
    if (world.grid().kind() == mesh::GridKind::latlon)
        latlonshelvesvoronoi(world, voronoi, pointdist);
    else
        sphereshelvesvoronoi(world, voronoi, pointdist);
}
