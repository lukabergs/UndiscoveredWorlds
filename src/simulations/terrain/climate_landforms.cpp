// Climate-driven dunes, salt pans and surface roughness.
// On the cell grid since T8 W2-T2 (docs/architecture/sphere-mesh.md sections 3, 4 and 11): the
// ergs and salt pans are drawn on the sphere (positions by area, templates through the tangent
// frame, distances in legacy pixels of km, chances per pixel area).

#include <algorithm>
#include <cmath>
#include <mutex>
#include <queue>
#include <vector>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"
#include "simulations/terrain/climate_landforms.hpp"
#include "core/grid_shapes.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_ridges.hpp"

using namespace std;
using mesh::CellId;
namespace cg = terrain::coastgrid;

namespace
{
// This ensures that the elevation of each erg remains constant: each 4-connected group of erg
// cells (special 120) takes its lowest base elevation (at least 2 and above sea level) as its
// surface.

void checkergelevation(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    int sealevel = world.sealevel();

    vector<char> checked(count, 0);
    vector<char> thiserg(count, 0);
    vector<CellId> erg;

    for (const CellId start : columnmajorcells(grid))
    {
        if (world.special(start) != 120 || checked[start] != 0)
            continue;

        // First, find the lowest elevation within this erg.

        float lowestelev = 1000000.0f;
        erg.clear();
        queue<CellId> q;
        q.push(start);
        thiserg[start] = 1;

        while (!q.empty())
        {
            const CellId c = q.front();
            q.pop();
            erg.push_back(c);
            lowestelev = min(lowestelev, world.nomfloat(c));

            for (int slot = 1; slot <= mesh::slotcount; slot += 2)
            {
                const std::int32_t n = grid.neighbour(c, slot);
                if (n < 0)
                    continue;
                if (world.special(static_cast<CellId>(n)) == 120 && thiserg[n] == 0) // If this node is erg
                {
                    thiserg[n] = 1;
                    q.push(static_cast<CellId>(n));
                }
            }
        }

        if (lowestelev < 2)
            lowestelev = 2;

        if (lowestelev <= sealevel)
            lowestelev = sealevel + 1;

        // Now apply that to the whole erg.

        for (const CellId c : erg)
        {
            thiserg[c] = 0;
            checked[c] = 1;
            world.setlakesurface(c, static_cast<float>(wholemetres(lowestelev)));
        }
    }
}

// The start point of a special lake on the regional map (makelakestartpoint on cells): one of
// its edge cells (a neighbour of another lake number, or none), drawn as the legacy loop drew it.
void makelakestartcell(planet& world, const vector<int>& thislake, int lakeno, const vector<CellId>& cells)
{
    const mesh::CellGrid& grid = world.grid();
    vector<CellId> edge;
    for (const CellId c : cells)
    {
        bool onedge = false;
        for (int slot = 1; slot <= mesh::slotcount && !onedge; slot++)
        {
            const std::int32_t n = grid.neighbour(c, slot);
            if (n >= 0 && thislake[n] != lakeno)
                onedge = true;
        }
        if (onedge)
            edge.push_back(c);
    }
    if (edge.empty())
        return;

    std::int32_t startpoint = -1;
    for (int crount = 0; crount <= 1000 && startpoint < 0; crount++)
        for (const CellId c : edge)
            if (random(1, 3) == 1)
                startpoint = static_cast<std::int32_t>(c);

    if (startpoint >= 0)
        world.setlakestart(static_cast<CellId>(startpoint), 1);
}

// drawspeciallake on the sphere: the template stamped at the cell's centre one legacy pixel per
// template pixel; the climate, outline and spacing tests within 1 and 4 pixels (in km).
void spheredrawspeciallake(planet& world, int shapenumber, CellId centre, int lakeno, vector<int>& thislake, boolshapetemplate laketemplate[], bool canoverlap, int special)
{
    const mesh::CellGrid& grid = world.grid();

    if (grid.neighbour(centre, 1) < 0 || grid.neighbour(centre, 5) < 0)
    {
        world.setlakesurface(centre, 0.0f);
        return;
    }

    const double pixelkm = terrain::legacycellkm(world);
    const double minlakedistance = 4.0; // Minimum distance between lakes (pixels).
    int maxmountainlakeheight = 1000; // Ignore mountains smaller than this.

    const bool leftr = random(0, 1) == 1; // If it's 1 then we reverse it left-right
    const bool downr = random(0, 1) == 1; // If it's 1 then we reverse it top-bottom

    const vector<CellId> cells = cg::stamptemplate(world, grid.centre(centre), laketemplate[shapenumber], leftr, downr);

    // First, check that the location is clear. Also, get the height of the lowest point in it, and make that the surface height.

    float surfaceheight = 1000000000.0f;

    for (const CellId c : cells)
    {
        surfaceheight = min(surfaceheight, world.nomfloat(c));

        if (world.sea(c) != 0 || world.mountainheightfloat(c) >= maxmountainlakeheight) // Don't try to put them on top of sea or mountains...
            return;

        if (special == 110 || special == 120)
        {
            for (const CellId n : terrain::cellswithinkm(world, c, 1.5 * pixelkm))
                if (world.climate(n) != 5 || world.riverjan(n) > 0 || world.riverjul(n) > 0)
                    return; // Don't put lakes onto inappropriate climates.
        }

        if (cg::outline(world, c))
            return;

        for (const CellId n : terrain::cellswithinkm(world, c, (minlakedistance + 0.5) * pixelkm))
        {
            if (thislake[n] != 0 && thislake[n] != lakeno && canoverlap == 0) // There's another lake here already!
                return;

            if (world.lakesurfacefloat(n) != 0 && thislake[n] != lakeno && (canoverlap == 0 || (special != 120 || world.special(n) != 120))) // There's another kind of lake here already! (Doesn't apply when both are ergs.)
                return;

            if (world.sea(n) == 1) // Sea here.
                return;
        }
    }

    if (cells.empty())
        return;

    // Now, do the actual thing.

    for (const CellId c : cells)
    {
        world.setlakesurface(c, static_cast<float>(wholemetres(surfaceheight))); // Put the lake on this bit of the lake map.
        world.setnomfloat(c, surfaceheight);
        world.setspecial(c, special);

        thislake[c] = lakeno; // Mark it on the keeping-track array too.

        if (world.mountainheightfloat(c) != 0) // Remove any mountains that might be here.
        {
            for (int dir = 1; dir <= 8; dir++)
                terrain::deleteridge(world, c, dir);
        }
    }

    // Now mark a "start point" somewhere on the edge of the "lake", which we will use when it comes to the regional map.

    makelakestartcell(world, thislake, lakeno, cells);
}

// createergs on the sphere: cluster templates at positions drawn by area, an erg's chance per
// legacy pixel area, the spacing to lakes within 2 pixels.
void spherecreateergs(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate shape[])
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = terrain::legacycellkm(world);
    int worldsize = world.size();

    vector<short> ergprobability(count, 0); // This is the probability of making an erg here.
    vector<short> overlapprobability(count, 0); // This is the probability of this erg being able to overlap its neighbours.

    int clusternumber = random(2, 50); // 50; Number of possible clusters of erg on the map.

    if (worldsize == 1)
        clusternumber = clusternumber * 4;

    if (worldsize == 2)
        clusternumber = clusternumber * 16;

    for (int n = 0; n < clusternumber; n++) // Draw some shapes on this array, to allow for erg clusters.
    {
        int shapenumber = random(1, 11);

        int probability = random(20, 200);

        int overlap = random(1, 100);

        const mesh::Vec3 centre = grid.centre(terrain::randomcellbyarea(world));

        const bool leftr = random(0, 1) == 1; // If it's 1 then we reverse it left-right
        const bool downr = random(0, 1) == 1; // If it's 1 then we reverse it top-bottom

        for (const CellId c : cg::stamptemplate(world, centre, shape[shapenumber], leftr, downr))
        {
            ergprobability[c] = static_cast<short>(probability);
            overlapprobability[c] = static_cast<short>(overlap);
        }
    }

    const double mindist = 2.0; // Minimum distance to existing lakes etc. (pixels)

    int lakeno = 0;

    vector<int> thislake(count, 0);

    for (const CellId c : columnmajorcells(grid))
    {
        if (world.climate(c) == 5 && ergprobability[c] != 0) // If this is a hot desert
        {
            if (cg::uniformdraw() < cg::cellchance(world, c, 1.0 / static_cast<double>(ergprobability[c])))
            {
                bool goahead = 1;

                for (const CellId n : terrain::cellswithinkm(world, c, (mindist + 0.5) * pixelkm)) // Check that there aren't any lakes too close.
                {
                    if (world.truelake(n) != 0)
                        goahead = 0;
                }

                if (goahead == 1)
                {
                    bool canoverlap = 0;

                    if (random(1, 100) < overlapprobability[c])
                        canoverlap = 1;

                    int shapenumber;

                    if (random(1, 4) == 1)
                    {
                        shapenumber = random(0, 11);
                        spheredrawspeciallake(world, shapenumber, c, lakeno, thislake, smalllake, canoverlap, 120); // 120 is the code for ergs.
                    }
                    else
                    {
                        if (random(1, 2) == 1)
                            shapenumber = random(0, 2);
                        else
                            shapenumber = random(0, 9);
                        spheredrawspeciallake(world, shapenumber, c, lakeno, thislake, largelake, canoverlap, 120); // 120 is the code for ergs.
                    }

                    lakeno++;
                }
            }
        }
    }
}

// createsaltpans on the sphere: a pan's chance per legacy pixel area, the spacing within 2 pixels.
void spherecreatesaltpans(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[])
{
    (void)largelake;
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);

    const double saltchance = 1.0 / 1000.0; // Probability of salt pans in any given desert pixel.
    const double mindist = 2.0; // Minimum distance to existing lakes etc. (pixels)

    int lakeno = 0;

    vector<int> thislake(grid.cellcount(), 0);

    for (const CellId c : columnmajorcells(grid))
    {
        if (world.climate(c) == 5) // If this is a hot desert
        {
            if (cg::uniformdraw() < cg::cellchance(world, c, saltchance))
            {
                bool goahead = 1;

                for (const CellId n : terrain::cellswithinkm(world, c, (mindist + 0.5) * pixelkm)) // Check that there aren't any other lakes too close.
                {
                    if (world.lakesurfacefloat(n) != 0)
                        goahead = 0;
                }

                if (goahead == 1)
                {
                    int shapenumber;

                    if (random(1, 4) != 1)
                        shapenumber = random(0, 3);
                    else
                        shapenumber = random(3, 11);

                    spheredrawspeciallake(world, shapenumber, c, lakeno, thislake, smalllake, 0, 110); // 110 is the code for salt pans.

                    lakeno++;
                }
            }
        }
    }
}
}

// This function puts ergs in deserts.

void createergs(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate shape[])
{
    spherecreateergs(world, smalllake, largelake, shape);

    // Some ergs may have overlapped, so now we need to make sure their elevation remains constant.

    checkergelevation(world);
}

// This function puts salt pans in deserts.

void createsaltpans(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[])
{
    spherecreatesaltpans(world, smalllake, largelake);
}

// This refines the roughness map so that it takes into account more factors, producing a roughness factor that will be used for the fractal terrain in the regional map.

void refineroughnessmap(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int maxelev = world.maxelevation();

    // This array will hold the valuemod for every tile on the map. The valuemod is used in the regional map to determine how rough the diamond-square routine makes the terrain.

    vector<float> valuemod(grid.cellcount(), 0.0f);

    // A neighbour, or the cell itself where there is none (the legacy clip of y at the poles).
    auto orself = [&](CellId c, int slot)
    {
        const std::int32_t n = grid.neighbour(c, slot);
        return n < 0 ? c : static_cast<CellId>(n);
    };

    for (const CellId c : columnmajorcells(grid))
    {
        if (world.sea(c) == 0)
        {
            if (world.basebelowsea(c) && world.basebelowsea(orself(c, 1)) && world.basebelowsea(orself(c, 5)) && world.basebelowsea(orself(c, 3)) && world.basebelowsea(orself(c, 7))) // Sea beds will be very smooth.
            {
                float m = 0.08f / (float)maxelev;
                float n = world.roughness(c);
                float extrabit = 0.01f * (float)random(1, 7);

                valuemod[c] = (m * n) + extrabit;
            }
            else // On land, the higher it is, the rougher it is.
                valuemod[c] = world.roughness(c) / (float)maxelev;
        }
        else // Sea beds - very smooth.
        {
            float m = 0.08f / (float)maxelev;
            float n = world.roughness(c);
            float extrabit = 0.01f * (float)random(1, 7);

            valuemod[c] = (m * n) + extrabit;
        }
    }

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (valuemod[c] < 0.0f)
            world.setroughness(c, 0.0f);
        else
            world.setroughness(c, valuemod[c]);
    }
}
