// Crater placement and impact relief used by editing tools.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp. On the cell grid since T8
// W2-T3 (docs/architecture/sphere-mesh.md sections 3, 11 and 12): centres drawn by area, radii
// in legacy pixels (terrain.sphere_noise.legacy_columns) as km, the density a sphere noise
// (terrain.sphere_noise.crater_density_*); the neighbour-scale tests (slope, the blur's ring)
// read the neighbour table.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "world/cell_queries.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;
using mesh::CellId;

namespace
{
// Salt of the craters' density noise (it replaces a fractal drawn from the ambient stream;
// combined with the number of earlier craters, so each press of the button draws another).
constexpr int cratersalt = 0x43a7;

void makecrater(planet& world, vector<int>& squareroot, vector<uint8_t>& placed, int thiscraterno, CellId centre, int size);

// The squared distance between two cell centres in legacy pixels, rounded to a whole number
// (the legacy loops' i * i + j * j).
int pixeldistancesquared(const planet& world, CellId a, CellId b, double pixelkm)
{
    const double d = terrain::distancekm(world, a, b) / pixelkm;
    return static_cast<int>(lround(d * d));
}

// The export pixel of a cell: its column and row on lat-lon; elsewhere the pixel of its centre
// on the 4 N x 2 N lat-lon export.
void exportpixel(const mesh::CellGrid& grid, CellId c, int& x, int& y)
{
    if (grid.kind() == mesh::GridKind::latlon)
    {
        x = grid.column(c);
        y = grid.row(c);
        return;
    }

    const double pi = 3.14159265358979323846;
    const int columns = 4 * grid.facecells();
    const int rows = columns / 2;
    x = clamp(static_cast<int>(floor((grid.longitude_rad(c) + pi) / (2.0 * pi) * columns)), 0, columns - 1);
    y = clamp(static_cast<int>(floor((0.5 * pi - grid.latitude_rad(c)) / pi * rows)), 0, rows - 1);
}
}

// This function creates craters.

void createcratermap(planet& world, int cratertotal, vector<int>& squareroot, bool custom)
{
    const CellId count = static_cast<CellId>(world.cellcount());
    int maxelev = world.maxelevation();

    if (world.size() == 0)
        cratertotal = cratertotal / 4;

    if (world.size() == 2)
        cratertotal = cratertotal * 4;

    int oldcraterno = world.craterno(); // Because this may not be the first time we've done this for this world (e.g. the "craters" button is clicked multiple times in the custom world screen).

    if (oldcraterno == MAXCRATERS)
        return;

    cratertotal = cratertotal + oldcraterno;

    if (cratertotal > MAXCRATERS)
        cratertotal = MAXCRATERS;

    world.setcraterno(cratertotal);

    int checktype = 0; // What additional kinds of checks there might be on crater placement.

    if (custom == 0) // If we haven't come here from the custom world page, we might restrict where craters can appear.
    {
        if (random(1, 10) == 1) // Based on slopes.
            checktype = 1;

        if (random(1, 4) != 1) // Based on elevation (high).
            checktype = 2;

        if (random(1, 8) == 1) // Based on elevation (low).
            checktype = 3;
    }

    int totalsea = 0;

    for (CellId c = 0; c < count; c++)
    {
        if (world.sea(c))
            totalsea++;
    }

    float lowestelev = static_cast<float>(maxelev);
    float highestelev = 0.0f;

    for (CellId c = 0; c < count; c++)
    {
        float thiselev = world.nomfloat(c);

        if (thiselev < lowestelev)
            lowestelev = thiselev;

        if (thiselev > highestelev)
            highestelev = thiselev;
    }

    int boundary = random(wholemetres(lowestelev), wholemetres(highestelev));

    // A broad fractal, to vary the density of craters (a sphere noise for the legacy grain-2 fractal).

    const vector<float> broad = terrain::spherefractal(world, deterministiccombine(deterministiccontextseed(world.seed(), cratersalt), oldcraterno),
        terrain::craterdensityfractal(), static_cast<float>(maxelev));

    // The centres of the earlier craters: two craters never share one.

    vector<uint8_t> centred(count, 0);

    for (int n = 0; n < oldcraterno; n++)
    {
        const CellId c = world.cratercell(n);

        if (c < count)
            centred[c] = 1;
    }

    vector<uint8_t> placed(count, 0);

    for (int craterno = oldcraterno; craterno < cratertotal; craterno++)
    {
        const CellId c = terrain::randomcellbyarea(world);

        if (random(1, maxelev) < broad[c])
        {
            bool makethiscrater = 0;

            float thiselev = world.nomfloat(c);

            if (checktype == 0)
                makethiscrater = 1;

            if (checktype == 1)
            {
                float slope = 0.0f;

                worldcells::anyinbox(world, c, 1, [&](CellId n)
                {
                    slope = max(slope, thiselev - world.nomfloat(n));
                    return false;
                });

                if (random(1, 60) < slope) // The steeper it is here, the more likely there are to be craters.
                    makethiscrater = 1;
            }

            if (checktype == 2)
            {
                if (thiselev > boundary)
                    makethiscrater = 1;
            }

            if (checktype == 3)
            {
                if (thiselev < boundary)
                    makethiscrater = 1;
            }

            if (centred[c]) // Don't allow two craters with precisely the same centre.
                makethiscrater = 0;

            if (makethiscrater)
            {
                int range = 10; // To ensure that most craters are small.

                if (random(1, 6) == 1)
                {
                    range = 25;

                    if (random(1, 8) == 1)
                    {
                        range = 50;

                        if (random(1, 20) == 1)
                        {
                            range = 70;

                            if (random(1, 50) == 1)
                                range = 100;
                        }
                    }
                }

                float thissize = (float)random(1, range);
                thissize = thissize / 100.0f;
                thissize = thissize * 15.0f; // In practice, can't draw craters of radius>15 on the regional map.

                if (thissize < 1.0f)
                    thissize = 1.0f;

                makecrater(world, squareroot, placed, craterno, c, (int)thissize);
                centred[c] = 1;
            }
        }
    }

    // If doing this has created sea when before there was none, probably lower the sea level.

    if (totalsea < 100)
    {
        int newtotalsea = 0;

        for (CellId c = 0; c < count; c++)
        {
            if (world.sea(c))
                newtotalsea++;
        }

        if (newtotalsea > 100)
        {
            if (random(1, 10) != 1)
            {
                float lowest = static_cast<float>(maxelev);

                for (CellId c = 0; c < count; c++)
                {
                    if (world.nomfloat(c) < lowest)
                        lowest = world.nomfloat(c);
                }

                int sealevel = wholemetres(lowest) - random(10, 500);

                if (sealevel < 1)
                    sealevel = 1;

                world.setsealevel(sealevel);
            }
        }
    }
}

namespace
{
// This function creates a crater on the world map: a floor over the cells within size legacy
// pixels of the centre cell, a rim on its outer pixel, a central peak and a blurred edge.
// placed (one flag per cell, all clear) marks the floor while the crater is drawn.

void makecrater(planet& world, vector<int>& squareroot, vector<uint8_t>& placed, int thiscraterno, CellId centre, int size)
{
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);

    int centreelev = size * 20;
    float rimelev = (float)(size * 250);

    float var = ((float)random(80, 120))/100.0f;

    centreelev = (int)((float)centreelev * var);

    var = ((float)random(80, 120)) / 100.0f;

    rimelev = rimelev * var;

    int centrex = 0;
    int centrey = 0;
    exportpixel(grid, centre, centrex, centrey);

    world.setcraterx(thiscraterno, centrex);
    world.setcratery(thiscraterno, centrey);
    world.setcraterelev(thiscraterno, centreelev);
    world.setcraterradius(thiscraterno, size);
    world.setcratercell(thiscraterno, centre);

    int sizecheck = size * size + size;
    int sizechecksmall = (size - 1) * (size - 1) + size - 1;

    vector<CellId> placedcells; // The cells of the floor, to clear placed afterwards.

    // First, find the base elevation for this crater.

    float baseelev = world.nomfloat(centre);

    // Now make a circular depression, with a rim.

    int depthmult = random(30, 70); // The higher this is, the deeper and more sloping the crater will be.

    float deepestelev = baseelev - static_cast<float>(size * depthmult);

    if (deepestelev < 1.0f)
        deepestelev = 1.0f;

    for (const CellId c : terrain::cellswithinkm(world, centre, sqrt(sizecheck + 0.5) * pixelkm))
    {
        int thischeck = pixeldistancesquared(world, centre, c, pixelkm);

        if (thischeck <= sizecheck)
        {
            // Work out how deep the crater should be at this point.

            int depth = squareroot[sizecheck - thischeck];

            depth = depth * depthmult;

            float thiselev = baseelev - depth;

            if (thiselev < 1.0f)
                thiselev = 1.0f;

            float oldelev = world.nomfloat(c);

            if (oldelev >= deepestelev)
            {
                placed[c] = 1;
                placedcells.push_back(c);

                world.setnomfloat(c, thiselev); // This is the flat crater floor.

                world.setmountainheight(c, 0); // Problem here.
                world.setmountainridge(c, 0);
                world.setvolcano(c, 0);
                world.setcraterrim(c, 0);
                world.setcratercentre(c, 0);

                if (thischeck <= sizecheck && thischeck >= sizechecksmall)
                {
                    world.setcraterrim(c, rimelev); // This is the ring of mountains around the rim.
                }
            }
        }
    }

    world.setcratercentre(centre, centreelev); // Do this only now, as we wiped any that the crater covers when making the crater floor just now.

    // Now raise the middle a bit around the central peak.

    int peakrad = (int)((float)size / 4.0f);

    if (peakrad < 1)
        peakrad = 1;

    int peakradcheck = peakrad * peakrad + peakrad;
    float peakbaseelev = (float)centreelev / 10.0f;

    for (const CellId c : terrain::cellswithinkm(world, centre, sqrt(peakradcheck + 0.5) * pixelkm))
    {
        int thischeck = pixeldistancesquared(world, centre, c, pixelkm);

        if (thischeck <= peakradcheck)
        {
            // Work out how much elevation to add at this point.

            float dist = 1.0f;

            if (thischeck > 1)
                dist = (float)squareroot[thischeck];

            dist = (float)peakrad - dist;

            if (dist < 0.0f)
                dist = 0.0f;

            float mult = dist / (float)peakrad;

            int newelev = (int)(peakbaseelev * mult);

            float thiselev = world.nomfloat(c) + static_cast<float>(newelev);

            world.setnomfloat(c, thiselev);
        }
    }

    // Now blur around the edges a little: the cells within two pixels of the floor that are
    // next to it (the 3 x 3 and 5 x 5 neighbour boxes of the legacy loops).

    int largesize = size + 2;

    for (const CellId c : terrain::cellswithinkm(world, centre, (largesize + 0.5) * pixelkm))
    {
        if (placed[c] == 1)
            continue;

        const bool nexttocrater = worldcells::anyinbox(world, c, 1, [&](CellId n) { return placed[n] == 1; });

        if (nexttocrater)
        {
            float highestcrater = 0.0f;

            worldcells::anyinbox(world, c, 2, [&](CellId n)
            {
                if (placed[n] == 1 && world.nomfloat(n) > highestcrater)
                    highestcrater = world.nomfloat(n);
                return false;
            });

            float newelev = (world.nomfloat(c) + highestcrater) / 2.0f;

            world.setnomfloat(c, newelev);
        }
    }

    for (const CellId c : placedcells)
        placed[c] = 0;
}
}
