// Procedural ocean ridges, faults, and trenches.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the sphere since T8 W2-T2.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <queue>

#include "core/classes.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_shapes.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_ridges.hpp"
#include "simulations/terrain/terrain_parameters.hpp"

using namespace std;
using mesh::CellId;

// The ocean features on the cell grid (T8 W2-T2; docs/architecture/sphere-mesh.md sections 3 and
// 11): the legacy construction in legacy pixels of km (tuning::terrain::sphere_noise_coast). The
// lattice of ridge points is a coarse cube grid of grid_size pixels per cell (its neighbour
// table gives each point's north, northeast, east and southeast neighbours); lines are the
// legacy Catmull-Rom curves in the tangent plane at their middle, joined by great-circle arcs;
// circles are the cells within the radius; fault seeds are drawn by area; the fault's area is
// the quadrilateral in the tangent plane at its middle. Ridge angles are relative to each cell's
// local east and north (sphere-mesh.md section 5); every chance is per legacy pixel area.

namespace
{
namespace cg = terrain::coastgrid;

bool sphereshelfedge(const planet& world, const worldgrid<bool>& shelves, CellId c)
{
    if (shelves[c] == 0)
        return false;
    bool found = false;
    cg::box(world.grid(), c, 1, [&](CellId n, int, int)
    {
        if (shelves[n] == 0)
            found = true;
        return !found;
    });
    return found;
}

// The legacy position (x east, y south, in pixels) of p in a tangent frame.
twofloats legacyplane(const mesh::TangentFrame& frame, const mesh::Vec3& p, double pixelkm)
{
    double east = 0.0, north = 0.0;
    cg::totangent(frame, p, east, north);
    twofloats point;
    point.x = static_cast<float>(east / pixelkm);
    point.y = static_cast<float>(-north / pixelkm);
    return point;
}

mesh::Vec3 fromlegacyplane(const mesh::TangentFrame& frame, twofloats point, double pixelkm)
{
    return mesh::fromtangent(frame, static_cast<double>(point.x) * pixelkm, -static_cast<double>(point.y) * pixelkm);
}

// drawoceanridgeline on the sphere: the curve through a randomly moved middle point, in the
// tangent plane at the middle of the two points, its samples joined by great-circle arcs.
void spheredrawoceanridgeline(const planet& world, const mesh::Vec3& from, const mesh::Vec3& to, vector<int>& array, int value, int keyfrom, int keyto)
{
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);
    const std::uint64_t lineseed = deterministiccontextseed(world.seed(), 0x51a9e245);

    const mesh::Vec3 sum = from + to;
    if (mesh::length(sum) < 1e-9)
        return;
    const mesh::TangentFrame frame = mesh::tangentframe(mesh::normalise(sum), grid.radius_m());

    const twofloats mm1 = legacyplane(frame, from, pixelkm);
    const twofloats mm3 = legacyplane(frame, to, pixelkm);

    int distance = static_cast<int>(fabs(mm1.x - mm3.x) + fabs(mm1.y - mm3.y));

    int minvar = 2;
    int maxvar = distance / 8;

    if (maxvar < 3)
        maxvar = 3;

    twofloats mm2;
    mm2.x = (mm1.x + mm3.x) / 2.0f + deterministicsignedrandom(lineseed ^ 0x1008ull, minvar, maxvar, keyfrom, keyto, value, 0);
    mm2.y = (mm1.y + mm3.y) / 2.0f + deterministicsignedrandom(lineseed ^ 0x1009ull, minvar, maxvar, keyfrom, keyto, value, 1);

    bool started = false;
    mesh::Vec3 previous;

    for (int n = 1; n <= 2; n++) // Two halves of the curve.
    {
        for (float t = 0.0f; t <= 1.0f; t = t + 0.01f)
        {
            const twofloats pt = n == 1 ? curvepos(mm1, mm1, mm2, mm3, t) : curvepos(mm1, mm2, mm3, mm3, t);
            const mesh::Vec3 p = fromlegacyplane(frame, pt, pixelkm);
            if (!started)
            {
                const CellId c = grid.locate(p);
                if (c != mesh::nocell)
                    array[c] = value;
                started = true;
            }
            else
            {
                for (const CellId c : mesh::cellsalongarc(grid, previous, p))
                    array[c] = value;
            }
            previous = p;
        }
    }
}

// The smallest legacy circle radius r >= 1 (offsets with k^2 + l^2 < r^2 + r) holding a point d pixels away.
int smallestradius(double d)
{
    const double d2 = d * d;
    int r = max(1, static_cast<int>(ceil((-1.0 + sqrt(1.0 + 4.0 * d2)) / 2.0)));
    while (static_cast<double>(r * r + r) <= d2)
        r++;
    while (r > 1 && static_cast<double>((r - 1) * (r - 1) + (r - 1)) > d2)
        r--;
    return r;
}

// A straight walk along the rift: from start, dist times to the first unvisited rift cell among
// the eight neighbours (in the legacy box order); false when the rift ends first (when required).
bool walkrift(const planet& world, CellId start, int dist, vector<CellId>& visited, CellId& end, bool required)
{
    const mesh::CellGrid& grid = world.grid();
    CellId at = start;
    for (int n = 1; n <= dist; n++)
    {
        if (find(visited.begin(), visited.end(), at) == visited.end())
            visited.push_back(at);
        std::int32_t next = -1;
        cg::box(grid, at, 1, [&](CellId m, int, int)
        {
            if (world.oceanrifts(m) != 0 && find(visited.begin(), visited.end(), m) == visited.end())
            {
                next = static_cast<std::int32_t>(m);
                return false;
            }
            return true;
        });
        if (next < 0)
        {
            if (required)
                return false;
            continue;
        }
        at = static_cast<CellId>(next);
    }
    end = at;
    return true;
}

bool insidequad(const twofloats quad[4], twofloats p)
{
    bool inside = false;
    for (int a = 0, b = 3; a < 4; b = a++)
    {
        if ((quad[a].y > p.y) != (quad[b].y > p.y))
        {
            const float x = quad[b].x + (p.y - quad[b].y) * (quad[a].x - quad[b].x) / (quad[a].y - quad[b].y);
            if (p.x < x)
                inside = !inside;
        }
    }
    return inside;
}

// createoceanfault on the sphere.
void spherecreateoceanfault(planet& world, CellId mid, int mindist, int maxdist, vector<int>& ridgesmap, int masksize)
{
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);
    const std::uint64_t faultseed = deterministiccontextseed(world.seed(), 0x7f4a7c15);

    int minshiftdist = 2;
    int maxshiftdist = 6;

    int dist = deterministicrandom(faultseed ^ 0x100aull, mindist, maxdist, mid, masksize); // Distance from the midpoint that the edge points should be.

    // Now go along the line twice - once in each direction.

    vector<CellId> visited;
    CellId ends[2] = { mid, mid };
    for (int dir = 0; dir < 2; dir++)
        if (!walkrift(world, mid, dist, visited, ends[dir], true)) // This means we reached the end of the rift line, so we can't do it.
            return;

    const mesh::TangentFrame frame = mesh::tangentframe(grid.centre(mid), grid.radius_m());
    const twofloats p1 = legacyplane(frame, grid.centre(ends[0]), pixelkm);
    const twofloats p2 = legacyplane(frame, grid.centre(ends[1]), pixelkm);

    // Now we need to find the coordinates of the four corners of the area to be moved.

    const float length = 70.0f; // Distance to go out from the central rift (pixels).

    twofloats corners[4]; // 1, 2 from the first end, 4, 3 from the second

    for (int n = 0; n < 2; n++) // Do it for each of the two end points
    {
        const twofloats p = n == 0 ? p1 : p2;

        float angle = (float)(0 - world.oceanridgeangle(ends[n]));

        while (angle < 0.0f)
            angle = angle + 360.0f;

        for (int nn = 0; nn < 2; nn++) // Two halves - one on each side of the rift.
        {
            float fangle = angle;

            if (nn == 1)
            {
                fangle = angle + 180.0f;

                if (fangle > 360.0f)
                    fangle = fangle - 360.0f;
            }

            fangle = fangle * 0.01745329f;

            twofloats end;
            end.x = p.x + (float)(int)(length * sin(fangle));
            end.y = p.y + (float)(int)(length * cos(fangle));

            corners[n == 0 ? nn : 3 - nn] = end;
        }
    }

    // Check that they don't cross over each other:

    if (p1.x < p2.x)
    {
        if (corners[0].x > corners[3].x || corners[1].x > corners[2].x)
            return;
    }
    else
    {
        if (corners[0].x < corners[3].x || corners[1].x < corners[2].x)
            return;
    }

    if (p1.y < p2.y)
    {
        if (corners[0].y > corners[3].y || corners[1].y > corners[2].y)
            return;
    }
    else
    {
        if (corners[0].y < corners[3].y || corners[1].y < corners[2].y)
            return;
    }

    // Now we need to create a mask of this area.

    vector<CellId> mask;
    vector<char> inmask(grid.cellcount(), 0);
    for (const CellId c : terrain::cellswithinkm(world, mid, (masksize / 2 + 1) * pixelkm))
    {
        if (mesh::dot(grid.centre(c), frame.origin) <= 0.0)
            continue;
        if (insidequad(corners, legacyplane(frame, grid.centre(c), pixelkm)))
        {
            mask.push_back(c);
            inmask[c] = 1;
        }
    }

    if (inmask[mid] == 0) // The shape does not hold its own middle - too weird to deal with!
        return;

    // Now we need to remove any links to ridges outside this area.

    for (const CellId c : mask)
    {
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t n = grid.neighbour(c, slot);
            if (n < 0 || inmask[n] == 0)
                terrain::deleteoceanridge(world, c, slot);
        }
    }

    // Now we need to record all of the ridge information for this area

    const size_t count = mask.size();
    vector<int> oldridgesmap(count), oldoceanrifts(count), oldoceanridges(count), oldoceanridgeangle(count), oldoceanridgeoffset(count);
    vector<float> oldoceanridgeheights(count);

    for (size_t k = 0; k < count; k++)
    {
        const CellId c = mask[k];
        oldridgesmap[k] = ridgesmap[c];
        oldoceanrifts[k] = world.oceanrifts(c);
        oldoceanridges[k] = world.oceanridges(c);
        oldoceanridgeheights[k] = world.oceanridgeheightsfloat(c);
        oldoceanridgeangle[k] = world.oceanridgeangle(c);
        oldoceanridgeoffset[k] = world.oceanridgeoffset(c);
    }

    // Now we need to work out the direction and distance that all of this is going to be shifted.

    float angle = (float)(0 - (world.oceanridgeangle(ends[0]) + world.oceanridgeangle(ends[1])) / 2);

    while (angle < 0.0f)
        angle = angle + 360.0f;

    if (deterministicrandom(faultseed ^ 0x100bull, 1, 2, mid, ends[0], ends[1]) == 1)
    {
        angle = angle + 180.0f;

        while (angle > 360.0f)
            angle = angle - 360.0f;
    }

    float fangle = angle * 0.01745329f;

    float shiftdist = (float)deterministicrandom(faultseed ^ 0x100cull, minshiftdist, maxshiftdist, mid, ends[0], ends[1]);

    int xshift = (int)(shiftdist * (float)sin(fangle));
    int yshift = (int)(shiftdist * (float)cos(fangle));

    // Now we rewrite all that information, offset by the new shift values.

    for (size_t k = 0; k < count; k++)
    {
        const CellId target = terrain::movedcell(world, mask[k], xshift * pixelkm, -yshift * pixelkm);
        ridgesmap[target] = oldridgesmap[k];
        world.setoceanrifts(target, oldoceanrifts[k]);
        world.setoceanridges(target, oldoceanridges[k]);
        world.setoceanridgeheights(target, oldoceanridgeheights[k]);
        world.setoceanridgeangle(target, oldoceanridgeangle[k]);
        world.setoceanridgeoffset(target, oldoceanridgeoffset[k]);
    }
}

void sphereoceanridges(planet& world, worldgrid<bool>& shelves)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = terrain::legacycellkm(world);
    int maxelev = world.maxelevation();
    const std::uint64_t ridgeseed = deterministiccontextseed(world.seed(), 0x0cea6e51);
    fast_srand(deterministicfastseed(ridgeseed));

    const vector<CellId> order = columnmajorcells(grid);

    vector<int> nearestshelfdist(count, 0);
    vector<std::int32_t> nearestshelf(count, -1);
    vector<char> boundaries(count, 0);
    vector<int> ridges(count, 0);
    vector<int> ridgesmap(count, 0);
    vector<int> ridgedistances(count, 0);

    // We need the edge points - points where the coastal shelves meet ocean - and every ocean
    // point's closest edge point, by a frontier expansion (diagonals only sometimes).

    vector<CellId> frontier;
    for (const CellId c : order)
    {
        if (sphereshelfedge(world, shelves, c))
        {
            nearestshelfdist[c] = 1;
            nearestshelf[c] = static_cast<std::int32_t>(c);
            frontier.push_back(c);
        }
    }

    int sweep = 2; // Because the first sweep was just setting up the edge points themselves

    while (!frontier.empty())
    {
        vector<CellId> nextfrontier;
        nextfrontier.reserve(frontier.size() * 2);

        for (const CellId c : frontier)
        {
            cg::box(grid, c, 1, [&](CellId n, int dx, int dy)
            {
                if (dx == 0 || dy == 0 || deterministicrandom(ridgeseed ^ 0x1001ull, 1, 2, c, n, sweep) == 1) // Only sometimes do diagonals, otherwise the end result looks too angular.
                {
                    if (shelves[n] == 0 && nearestshelfdist[n] == 0)
                    {
                        nearestshelfdist[n] = sweep;
                        nearestshelf[n] = nearestshelf[c];
                        nextfrontier.push_back(n);
                    }
                }
            });
        }

        frontier.swap(nextfrontier);
        sweep++;
    }

    // Now we need to find where the zones meet: neighbours whose edge points lie far apart.

    const double maxdiff_m = tuning::terrain::oceanridges::boundaryMaxSourceDifference * pixelkm * 1000.0;

    for (const CellId c : order)
    {
        if (nearestshelfdist[c] <= 0 || shelves[c] != 0 || nearestshelf[c] < 0)
            continue;
        const mesh::Vec3 source = grid.centre(static_cast<CellId>(nearestshelf[c]));
        cg::box(grid, c, 1, [&](CellId n, int, int)
        {
            if (shelves[n] == 0 && boundaries[n] == 0 && nearestshelf[n] >= 0
                && mesh::greatcircle_m(source, grid.centre(static_cast<CellId>(nearestshelf[n])), grid.radius_m()) > maxdiff_m)
                boundaries[c] = 1;
        });
    }

    // We've got the boundaries marked out. Now we need a grid of points: a coarse cube of
    // grid_size pixels per cell, a point where a boundary lies within half a cell of its centre.

    const int gridsize = max(1, tuning::terrain::oceanridges::gridSize);
    const int halfgrid = gridsize / 2 - 1;
    const int coarsecells = max(1, static_cast<int>(lround(tuning::terrain::sphere_noise::legacyColumns / (4.0 * gridsize))));
    const shared_ptr<const mesh::CellGrid> coarse = mesh::makecubegrid(coarsecells, grid.radius_m());
    const size_t points = coarse->cellcount();
    const vector<CellId> coarseorder = columnmajorcells(*coarse);

    vector<int> gridvalue(points, 0);

    for (const CellId k : coarseorder)
    {
        const mesh::Vec3 p = coarse->centre(k);
        double best = 1e30;
        for (const CellId c : terrain::cellswithinkm(world, grid.locate(p), (halfgrid + 0.5) * pixelkm))
        {
            if (boundaries[c] == 0)
                continue;
            const double d = mesh::greatcircle_rad(p, grid.centre(c));
            if (d < best)
            {
                best = d;
                gridvalue[k] = nearestshelfdist[c];
            }
        }
    }

    auto coarsepoint = [&](CellId k, int slot) -> std::int32_t
    {
        const std::int32_t n = coarse->neighbour(k, slot);
        return n >= 0 && gridvalue[n] != 0 ? n : -1;
    };

    for (const CellId k : coarseorder) // This one is to get rid of extra points around the diagonals.
    {
        if (gridvalue[k] == 0)
            continue;
        if (coarsepoint(k, 3) >= 0 && coarsepoint(k, 5) >= 0)
            gridvalue[k] = 0;
        else if (coarsepoint(k, 3) >= 0 && coarsepoint(k, 1) >= 0)
            gridvalue[k] = 0;
    }

    // Now we have a grid of points. We need to assign each one a unique number.

    vector<int> gridnumber(points, 0);
    vector<CellId> pointcell(1, 0);
    int totalpoints = 0;

    for (const CellId k : coarseorder)
    {
        if (gridvalue[k] != 0)
        {
            totalpoints++;
            gridnumber[k] = totalpoints;
            pointcell.push_back(k);
        }
    }

    // Each point's location, its distance to the closest shelf, and the points north, northeast,
    // east and southeast of it; then each point shifted by a sphere noise and a random amount.

    // The shift fields at the cells under the points (grain 8: the terrain fractal's set).
    const vector<float> shiftx = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x0cea0001), terrain::terrainfractal(), (float)maxelev);
    const vector<float> shifty = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x0cea0002), terrain::terrainfractal(), (float)maxelev);

    int maxshift = tuning::terrain::oceanridges::maxShift;
    int maxadditionalshift = tuning::terrain::oceanridges::maxAdditionalShift;
    int minadditionalshift = tuning::terrain::oceanridges::minAdditionalShift;

    float div = (float)maxelev / (float)(maxshift * 2);

    vector<mesh::Vec3> position(static_cast<size_t>(totalpoints) + 1);

    for (int thispoint = 1; thispoint <= totalpoints; thispoint++)
    {
        const CellId k = pointcell[static_cast<size_t>(thispoint)];
        const mesh::Vec3 p = coarse->centre(k);

        const CellId under = grid.locate(p);
        float xshift = (float)(int)shiftx[under] / div - (float)maxshift;
        float yshift = (float)(int)shifty[under] / div - (float)maxshift;

        xshift = max(-(float)maxshift, min((float)maxshift, xshift));
        yshift = max(-(float)maxshift, min((float)maxshift, yshift));

        int xextrashift = deterministicsignedrandom(ridgeseed ^ 0x1003ull, minadditionalshift, maxadditionalshift, thispoint, k, 0);
        int yextrashift = deterministicsignedrandom(ridgeseed ^ 0x1004ull, minadditionalshift, maxadditionalshift, thispoint, k, 1);

        const double east = static_cast<double>((int)xshift + xextrashift) * pixelkm;
        const double south = static_cast<double>((int)yshift + yextrashift) * pixelkm;
        position[static_cast<size_t>(thispoint)] = mesh::fromtangent(mesh::tangentframe(p, grid.radius_m()), east, -south);
    }

    // Now join up the points.

    for (int thispoint = 1; thispoint <= totalpoints; thispoint++)
    {
        const CellId k = pointcell[static_cast<size_t>(thispoint)];
        for (int slot = 1; slot <= 4; slot++)
        {
            const std::int32_t other = coarsepoint(k, slot);
            if (other < 0)
                continue;
            const int otherpoint = gridnumber[other];
            spheredrawoceanridgeline(world, position[static_cast<size_t>(thispoint)], position[static_cast<size_t>(otherpoint)], ridges, gridvalue[k], thispoint, otherpoint);
        }
    }

    // Now get rid of excess points around the diagonals.

    for (const CellId c : order)
    {
        if (ridges[c] == 0)
            continue;
        const std::int32_t east = grid.neighbour(c, 3);
        const std::int32_t south = grid.neighbour(c, 5);
        const std::int32_t north = grid.neighbour(c, 1);
        const std::int32_t southeast = grid.neighbour(c, 4);
        const std::int32_t northeast = grid.neighbour(c, 2);

        if (south >= 0 && southeast >= 0 && ridges[southeast] != 0)
        {
            ridges[south] = 0;
            if (east >= 0)
                ridges[east] = 0;
        }

        if (north >= 0 && northeast >= 0 && ridges[northeast] != 0)
        {
            ridges[north] = 0;
            if (east >= 0)
                ridges[east] = 0;
        }
    }

    // Now make sure each ridge point is adjacent to no more than one other ridge point.

    bool found = 0;

    do
    {
        found = 0;

        for (const CellId c : order)
        {
            if (ridges[c] == 0)
                continue;

            int neighbours = -1; // Because we'll add one for itself

            cg::box(grid, c, 1, [&](CellId n, int, int)
            {
                if (ridges[n] != 0)
                    neighbours++;
            });

            if (neighbours > 2)
            {
                found = 1;

                cg::box(grid, c, 1, [&](CellId n, int, int)
                {
                    if (ridges[n] != 0)
                    {
                        ridges[n] = 0;
                        return false;
                    }
                    return true;
                });
            }
        }

    } while (found == 1);

    vector<CellId> ridgecells;

    for (const CellId c : order)
        if (ridges[c] != 0)
            ridgecells.push_back(c);

    // Now raise the land around the ridges.

    const int firstmaxradius = tuning::terrain::oceanridges::firstPassMaxRadius;
    int heightmult = tuning::terrain::oceanridges::firstPassHeightMultiplier;
    int maxvolcanoradius = tuning::terrain::oceanridges::maxVolcanoRadius;

    vector<float> thisheightperthousand(static_cast<size_t>(firstmaxradius) + 1); // Create a lookup table for these values, for speed.

    thisheightperthousand[1] = 1.1037f;

    for (int n = 2; n <= firstmaxradius; n++)
        thisheightperthousand[static_cast<size_t>(n)] = thisheightperthousand[static_cast<size_t>(n) - 1] * thisheightperthousand[1];

    for (const CellId r : ridgecells)
    {
        const mesh::Vec3 p = grid.centre(r);
        const float riftheightdiv = (float)(ridges[r] * heightmult) / 1000.0f;

        for (const CellId n : terrain::cellswithinkm(world, r, (firstmaxradius + 0.5) * pixelkm))
        {
            const int radius = smallestradius(terrain::distancekm(world, r, n) / pixelkm);
            if (radius > firstmaxradius)
                continue;
            const int mult = firstmaxradius + 1 - radius;
            const int thisheight = (int)(riftheightdiv * thisheightperthousand[static_cast<size_t>(mult)]);

            if (ridgesmap[n] < thisheight)
            {
                ridgesmap[n] = thisheight;
                ridgedistances[n] = mult;
            }
        }
    }

    for (CellId c = 0; c < count; c++)
        if (ridgesmap[c] > 0)
            world.setoceanridgeheights(c, static_cast<float>(ridgesmap[c]));

    for (CellId c = 0; c < count; c++)
    {
        if (ridgedistances[c] != 0 && ridgedistances[c] <= maxvolcanoradius)
        {
            const int thisheight = ridgesmap[c];

            if (thisheight > 0)
            {
                unsigned int hash = static_cast<unsigned int>(world.seed());
                hash ^= static_cast<unsigned int>(c) + 0x9e3779b9u + (hash << 6) + (hash >> 2);
                hash ^= static_cast<unsigned int>(thisheight) + 0x9e3779b9u + (hash << 6) + (hash >> 2);

                const double perpixel = static_cast<double>(std::max(1, thisheight / 10)) / 4000000.0;

                if (static_cast<double>(hash % 4000000u) / 4000000.0 < cg::cellchance(world, c, perpixel))
                {
                    const int minvolcano = std::max(1, thisheight / 4);
                    const int maxvolcano = std::max(minvolcano, (thisheight / 4) * 3);
                    const unsigned int span = static_cast<unsigned int>(maxvolcano - minvolcano + 1);

                    hash ^= 0x85ebca6bu + (hash << 6) + (hash >> 2);
                    world.setvolcano(c, static_cast<float>(minvolcano + static_cast<int>(hash % span)));
                }
            }
        }
    }

    const int maxradius = tuning::terrain::oceanridges::secondPassMaxRadius;
    heightmult = tuning::terrain::oceanridges::secondPassHeightMultiplier;
    vector<int> secondpass(count, 0);

    for (const CellId r : ridgecells)
    {
        const mesh::Vec3 p = grid.centre(r);
        const int ridgeheight = ridges[r] * heightmult;
        const int heightdivs = ridgeheight / maxradius;

        for (const CellId n : terrain::cellswithinkm(world, r, (maxradius + 0.5) * pixelkm))
        {
            const int radius = smallestradius(terrain::distancekm(world, r, n) / pixelkm);
            if (radius > maxradius)
                continue;
            const int mult = maxradius + 1 - radius;
            int thisheight = heightdivs * mult;
            thisheight = (thisheight * 3 + ridgeheight / radius) / 4; // 3, 4

            if (secondpass[n] < thisheight)
                secondpass[n] = thisheight;
        }
    }

    for (CellId c = 0; c < count; c++)
        ridgesmap[c] = max(ridgesmap[c], secondpass[c]);

    // Add the rift in the middle

    vector<CellId> riftcells;

    for (const CellId r : ridgecells)
    {
        ridgesmap[r] = ridgesmap[r] / 2;
        world.setoceanrifts(r, ridgesmap[r]);

        if (ridgesmap[r] != 0)
            riftcells.push_back(r);
    }

    // Now we have to work out the angles of lines crossing the rifts, relative to each cell's
    // east and north, blended with an angle field (the legacy fractal's values wrapped at 360).

    const vector<float> anglefield = terrain::spherefractalwrapped(world, deterministiccontextseed(world.seed(), 0x0cea0003), cg::ridgeanglefractal(), 360.0f);

    short dist = tuning::terrain::oceanridges::ridgeAngleSearchDistance;

    for (const CellId c : riftcells)
    {
        // Now go along the line twice - once in each direction.

        vector<CellId> visited;
        CellId ends[2] = { c, c };
        for (int dir = 0; dir < 2; dir++)
            walkrift(world, c, dist, visited, ends[dir], false);

        // Now we have the two nearby points, we just find the angle between them.

        const mesh::Vec3 along = grid.centre(ends[1]) - grid.centre(ends[0]);
        const double east = mesh::dot(along, grid.east(c));
        const double south = -mesh::dot(along, grid.north(c));

        float angle = (float)atan2(south, east) * 180.0f / 3.14159265358979323846f; // This gives us the angle of the ridge

        angle = angle + 180.0f; // Because we want the line that crosses the ridge

        while (angle > 360.0f)
            angle = angle - 360.0f;

        int angle2 = wrappedaverage((int)angle, 1 + (int)anglefield[c], 360);

        while (angle2 > 180)
            angle2 = angle2 - 180;

        world.setoceanridgeangle(c, angle2);
    }

    // Mark out the ridge mountains, following the contours of the raised land

    for (const CellId c : order)
    {
        if (ridgedistances[c] == 0)
            continue;
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t n = grid.neighbour(c, slot);
            if (n >= 0 && ridgedistances[n] == ridgedistances[c] && !terrain::getoceanridge(world, c, slot))
                world.setoceanridges(c, world.oceanridges(c) + (1 << (slot - 1)));
        }
    }

    // Now remove extraneous ridges

    for (const CellId c : order)
    {
        if (world.oceanridges(c) == 0)
            continue;

        short total = 0;

        for (int dir = 1; dir <= 8; dir++)
            if (terrain::getoceanridge(world, c, dir))
                total++;

        if (total > 2)
        {
            int extra = total - 2;

            for (int n = 1; n <= extra; n++)
            {
                int a = 1;
                int b = 8;
                int step = 1;

                if (deterministicrandom(ridgeseed ^ 0x1005ull, 1, 2, c, n, extra) == 1)
                {
                    a = 8;
                    b = 1;
                    step = -1;
                }

                for (int dir = a; dir != b; dir = dir + step)
                {
                    if (terrain::getoceanridge(world, c, dir))
                        terrain::deleteoceanridge(world, c, dir);
                }
            }
        }
    }

    // Now a field that will be used to displace the ridges in the regional map.

    int maxdisplace = tuning::terrain::oceanridges::regionalDisplacement;

    const vector<float> ridgefractal = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x0cea0004), cg::finefractal(), (float)maxelev);

    div = (float)maxelev / ((float)maxdisplace * 2.0f);

    for (CellId c = 0; c < count; c++)
    {
        float amount = (float)(int)ridgefractal[c];

        amount = amount / div;
        amount = amount - (float)maxdisplace;

        if (amount < 0.0f - (float)maxdisplace)
            amount = 0.0f - (float)maxdisplace;

        if (amount > (float)maxdisplace)
            amount = (float)maxdisplace;

        world.setoceanridgeoffset(c, (int)amount);
    }

    // Now we need to disrupt the rifts by adding some faults, at points drawn by area, one per
    // fault_step pixels squared.

    const double faultspacing_m = tuning::terrain::oceanridges::faultStep * pixelkm * 1000.0;
    const int faultseeds = max(1, static_cast<int>(lround(grid.totalarea_m2() / (faultspacing_m * faultspacing_m))));

    for (int n = 1; n <= tuning::terrain::oceanridges::faultPasses; n++)
    {
        int lookdist = tuning::terrain::oceanridges::faultLookDistance;

        for (int k = 0; k < faultseeds; k++)
        {
            const double u = mesh::uniformfrombits(deterministiccombine(ridgeseed ^ 0x1006ull, n, k, 0));
            const double v = mesh::uniformfrombits(deterministiccombine(ridgeseed ^ 0x1007ull, n, k, 1));
            const CellId c = grid.locate(mesh::uniformpoint(u, v));

            if (c == mesh::nocell || world.sea(c) != 1)
                continue;

            for (const CellId r : terrain::cellswithinkm(world, c, (lookdist + 0.5) * pixelkm))
            {
                if (world.oceanrifts(r) != 0)
                {
                    int maxwidth = 30;
                    int minwidth = 6;

                    int masksize = 200 + maxwidth;

                    spherecreateoceanfault(world, r, minwidth, maxwidth, ridgesmap, masksize);
                    break;
                }
            }
        }
    }

    // Remove any ridge-related stuff that's on continental plates

    for (CellId c = 0; c < count; c++)
    {
        if (shelves[c] == 1 || world.sea(c) == 0)
        {
            ridgesmap[c] = 0;
            world.setoceanrifts(c, 0);
            world.setoceanridges(c, 0);
            world.setoceanridgeheights(c, 0.0f);
            world.setoceanridgeangle(c, 0);
            world.setoceanridgeoffset(c, 0);
        }
    }

    // Draw the raised land onto the actual sea bed

    for (CellId c = 0; c < count; c++)
        world.setnomfloat(c, world.nomfloat(c) + static_cast<float>(ridgesmap[c]));
}

// createoceantrenches on the sphere: circles of trench at shelf edges where the noise is high.
void sphereoceantrenches(planet& world, worldgrid<bool>& shelves)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = terrain::legacycellkm(world);
    int maxelev = world.maxelevation();

    int trenchmin = maxelev / 2; // Values in the fractal higher than this will spawn trenches.

    int maxradius = 20;

    int div = max(1, (maxelev - trenchmin) / maxradius);

    vector<char> trenchmap(count, 0);
    const vector<float> fractal = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x7e4c0001), cg::trenchfractal(), (float)maxelev);
    const vector<CellId> order = columnmajorcells(grid);

    for (const CellId c : order)
    {
        if (shelves[c] == 1 && (int)fractal[c] > trenchmin && sphereshelfedge(world, shelves, c))
        {
            int radius = ((int)fractal[c] - trenchmin) / div;

            if (radius > 0)
                for (const CellId n : terrain::cellswithinkm(world, c, (radius + 0.5) * pixelkm))
                    trenchmap[n] = 1;
        }
    }

    // Now just add the trenches to the world map.

    for (const CellId c : order)
    {
        if (trenchmap[c] == 1 && world.sea(c) == 1 && shelves[c] == 0 && world.oceanridges(c) == 0)
        {
            float newval = world.nomfloat(c) - static_cast<float>(random(4900, 5100));

            if (newval < 1.0f)
                newval = 1.0f;

            world.setnomfloat(c, newval);
        }
    }
}
}

// This creates mid-ocean ridges.

void createoceanridges(planet& world, worldgrid<bool>& shelves)
{
    sphereoceanridges(world, shelves);
}

// This creates the ocean trenches.

void createoceantrenches(planet& world, worldgrid<bool>& shelves)
{
    sphereoceantrenches(world, shelves);
}
