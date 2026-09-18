// Crater placement and impact relief used by editing tools.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "classes.hpp"
#include "functions.hpp"
#include "planet.hpp"

using namespace std;

namespace
{
void makecrater(planet& world, std::vector <int> &squareroot, int thiscraterno, int centrex, int centrey, int size);
}

// This function creates craters.

void createcratermap(planet& world, int cratertotal, vector<int>& squareroot, bool custom)
{
    int width = world.width();
    int height = world.height();
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

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j))
                totalsea++;
        }
    }

    int lowestelev = maxelev;
    int highestelev = 0;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int thiselev = world.nom(i, j);

            if (thiselev < lowestelev)
                lowestelev = thiselev;

            if (thiselev > highestelev)
                highestelev = thiselev;
        }
    }

    int boundary = random(lowestelev, highestelev);

    // Create a broad fractal, to vary the density of craters.

    int grain = 2; // Level of detail on this fractal map.
    float valuemod = 0.002f;
    float valuemod2 = 0.002f;

    vector<vector<int>> broad(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(broad, width, height, grain, valuemod, valuemod2, 1, maxelev, 1, 0);

    for (int craterno = oldcraterno; craterno < cratertotal; craterno++)
    {
        int thisx = random(1, width);
        int thisy = random(1, height);

        if (random(1, maxelev) < broad[thisx][thisy])
        {
            bool makethiscrater = 0;

            int thiselev = world.nom(thisx, thisy);

            if (checktype == 0)
                makethiscrater = 1;

            if (checktype == 1)
            {
                int slope = 0;

                for (int i = thisx - 1; i <= thisx + 1; i++)
                {
                    int ii = i;

                    if (ii<0 || ii>width)
                        ii = wrap(ii, width);

                    for (int j = thisy - 1; j <= thisy + 1; j++)
                    {
                        if (j >= 0 && j <= height)
                        {
                            int thisslope = thiselev - world.nom(ii, j);

                            if (thisslope > slope)
                                slope = thisslope;
                        }
                    }
                }

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

            for (int n = 0; n <= craterno; n++) // Don't allow two craters with precisely the same centre.
            {
                if (world.craterx(n) == thisx && world.cratery(n) == thisy)
                    makethiscrater = 0;
            }

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

                makecrater(world, squareroot, craterno, thisx, thisy, (int)thissize);
            }
        }
    }

    // If doing this has created sea when before there was none, probably lower the sea level.

    if (totalsea < 100)
    {
        int newtotalsea = 0;

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
                if (world.sea(i, j))
                    newtotalsea++;
        }

        if (newtotalsea > 100)
        {
            if (random(1, 10) != 1)
            {
                int lowest = maxelev;

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        if (world.nom(i, j) < lowest)
                            lowest = world.nom(i, j);
                    }
                }

                int sealevel = lowest - random(10, 500);

                if (sealevel < 1)
                    sealevel = 1;

                world.setsealevel(sealevel);
            }
        }
    }
}

namespace
{
// This function creates a crater on the world map.

void makecrater(planet& world, vector <int> &squareroot, int thiscraterno, int centrex, int centrey, int size)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    int centreelev = size * 20;
    int rimelev = size * 250;

    float var = ((float)random(80, 120))/100.0f;

    centreelev = (int)((float)centreelev * var);

    var = ((float)random(80, 120)) / 100.0f;

    rimelev = (int)((float)rimelev * var);

    world.setcraterx(thiscraterno, centrex);
    world.setcratery(thiscraterno, centrey);
    world.setcraterelev(thiscraterno, centreelev);
    world.setcraterradius(thiscraterno, size);

    int sizecheck = size * size + size;
    int sizechecksmall = (size - 1) * (size - 1) + size - 1;

    vector<vector<bool>> placed(width + 1, vector<bool>(height + 1, 0)); // This notes where the crater has been placed.

    // First, find the base elevation for this crater.

    int baseelev = world.nom(centrex,centrey);

    // Now make a circular depression, with a rim.

    int depthmult = random(30, 70); // The higher this is, the deeper and more sloping the crater will be.

    int deepestelev = baseelev-(size * depthmult);

    if (deepestelev < 1)
        deepestelev = 1;

    int aveelev = deepestelev + (baseelev - deepestelev) / 2; // Roughly average elevation of the crater floor.

    for (int i = 0 - size; i <= size; i++)
    {
        int ii = centrex + i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);
        
        for (int j = 0 - size; j <= size; j++)
        {
            int jj = centrey + j;

            if (jj >= 0 && jj <= height)
            {
                int thischeck = i * i + j * j;
                
                if (thischeck <= sizecheck)
                {
                    // Work out how deep the crater should be at this point.

                    int depth = squareroot[sizecheck - thischeck];

                    depth = depth * depthmult;

                    int thiselev = baseelev - depth;

                    if (thiselev < 1)
                        thiselev = 1;

                    int oldelev = world.nom(ii, jj);

                    if (oldelev >= deepestelev)
                    {                     
                        placed[ii][jj] = 1;

                        world.setnom(ii, jj, thiselev); // This is the flat crater floor.

                        world.setmountainheight(ii, jj, 0); // Problem here.
                        world.setmountainridge(ii, jj, 0);
                        world.setvolcano(ii, jj, 0);
                        world.setcraterrim(ii, jj, 0);
                        world.setcratercentre(ii, jj, 0);

                        if (thischeck <= sizecheck && thischeck >= sizechecksmall)
                        {
                            world.setcraterrim(ii, jj, rimelev); // This is the ring of mountains around the rim.
                        }
                    }
                }
            }
        }
    }

    world.setcratercentre(centrex, centrey, centreelev); // Do this only now, as we wiped any that the crater covers when making the crater floor just now.

    // Now raise the middle a bit around the central peak.

    int peakrad = (int)((float)size / 4.0f);

    if (peakrad < 1)
        peakrad = 1;

    int peakradcheck = peakrad * peakrad + peakrad;
    float peakbaseelev = (float)centreelev / 10.0f;

    for (int i = 0 - peakrad; i <= peakrad; i++)
    {
        int ii = centrex + i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = 0 - peakrad; j <= peakrad; j++)
        {
            int jj = centrey + j;

            if (jj >= 0 && jj <= height)
            {
                int thischeck = i * i + j * j;

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

                    int thiselev = world.nom(ii, jj) + newelev;

                    world.setnom(ii, jj, thiselev);
                }
            }
        }
    }

    // Now blur around the edges a little.

    int largesize = size + 2;

    for (int i = 0 - largesize; i <= largesize; i++)
    {
        int ii = centrex + i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = 0 - largesize; j <= largesize; j++)
        {
            int jj = centrey + j;

            if (jj >= 0 && jj <= height)
            {                
                if (placed[ii][jj] == 0)
                {
                    bool nexttocrater = 0;
                    
                    for (int k = ii - 1; k <= ii + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = jj - 1; l <= jj + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (placed[kk][l] == 1)
                                {
                                    nexttocrater = 1;
                                    k = ii + 1;
                                    l = jj + 1;
                                }
                            }
                        }
                    }

                    if (nexttocrater)
                    {
                        int highestcrater = 0;

                        for (int k = ii - 2; k <= ii + 2; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = jj - 2; l <= jj + 2; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (placed[kk][l] == 1 && world.nom(kk,l)>highestcrater)
                                        highestcrater = world.nom(kk, l);
                                }
                            }
                        }

                        int newelev = (world.nom(ii, jj) + highestcrater) / 2;

                        world.setnom(ii, jj, newelev);
                    }
                }
            }
        }
    }
}
}
