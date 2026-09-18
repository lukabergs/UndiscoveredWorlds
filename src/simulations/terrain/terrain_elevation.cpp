// Elevation smoothing, relief noise, and depression filling.
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
bool checkdepressiontile(planet& world, std::vector<std::vector<int>>& filledmap, int i, int j, int e, int neighbours[8][2], bool somethingdone, std::vector<std::vector<int>>& noise);
}

// This function smoothes, but without turning any land to sea or vice versa.

void smoothland(planet& world, int amount)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();
    float seaamount = (float)(amount * 3);
    float div = maxelev / seaamount;

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            fractal[i][j] = 0;
    }

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float crount = 0.0f;
            float ave = 0.0f;
            int origland = 0;
            int thisamount = amount;

            if (world.sea(i, j) == 0) // Check to see whether this point is originally land or sea
                origland = 1;
            else
                thisamount = (int)((float)fractal[i][j] / div);

            if (thisamount < 1)
                thisamount = 1;

            bool goahead = 1;
            int seacheck = 5;

            if (goahead == 1)
            {
                for (int k = i - thisamount; k <= i + thisamount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - thisamount; l <= j + thisamount; l++)
                    {
                        if (l > 0 && l < height)
                        {
                            ave = ave + (float)world.nom(kk, l);
                            crount++;
                        }
                    }
                }

                if (crount > 0.0f)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < maxelev)
                    {
                        int newland = 0;

                        if (ave > (float)sealevel) // Now check to see whether the new value is land or sea
                            newland = 1;

                        if (origland == 0 && newland == 0)
                            world.setnom(i, j, (int)ave);
                    }
                }
            }
        }
    }
}

// This function smoothes only the sea.

void smoothonlysea(planet& world, int amount)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float crount = 0.0f;
            float ave = 0.0f;
            int origland = 0;
            int thisamount = amount;

            if (world.nom(i, j) <= sealevel)
            {
                for (int k = i - thisamount; k <= i + thisamount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - thisamount; l <= j + thisamount; l++)
                    {
                        if (l > 0 && l < height && world.nom(kk, l) <= sealevel)
                        {
                            ave = ave + (float)world.nom(kk, l);
                            crount++;
                        }
                    }
                }

                if (crount > 0.0f)
                {
                    ave = ave / crount;

                    if (ave > 0.0f && ave < maxelev)
                        world.setnom(i, j, (int)ave);
                }
            }
        }
    }
}

// This function creates areas of raised elevation around where canyons might form later.

void createextraelev(planet& world)
{
    int width = world.width();
    int height = world.height();
    float gravity = world.gravity();

    vector<vector<int>> tempelev(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    float maxextra = 1000.0f; // Maximum amount that the extraelev array can have.

    float maxelev = (float)world.maxelevation();

    float div = maxelev / (maxextra * 2.0f);

    // First create a fractal for this map.

    int grain = 16; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    float valuemod2 = 0.2f;

    createfractal(tempelev, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    int warpfactor = random(20, 80);
    warp(tempelev, width, height, (int)maxelev, warpfactor, 1);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                world.setextraelev(i, j, tempelev[i][j]);
                tempelev[i][j] = 0;
            }
        }
    });

    int v = random(3, 6);
    valuemod2 = (float)v;

    createfractal(tempelev, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    warpfactor = random(20, 80);
    warp(tempelev, width, height, (int)maxelev, warpfactor, 1);

    // Now alter it a bit.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 1)
                    world.setextraelev(i, j, 0);
                else
                {
                    float e = (float)world.extraelev(i, j);
                    e = e / div;
                    e = e - maxextra;

                    if (e < 0.0f)
                        e = 0.0f;

                    if (e > maxextra)
                        e = maxextra;

                    e = e * gravity;

                    world.setextraelev(i, j, (int)e);
                }
            }
        }
    });

    // That gives us a *lot* of extra elevation across the map. We want to apply it a bit more judiciously.
    // So we use a second fractal to mask the first, so that fewer areas get extra elevation.

    float multfactor = (float)random(1, 10);
    multfactor = multfactor / 10.0f + 1.0f; // Different worlds will have different amounts of extra elevation.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 0)
                {
                    float e = (float)tempelev[i][j];
                    e = e / div;
                    e = e - maxextra;

                    e = e * multfactor; // The higher this is, the more widespread the extra elevation will be.

                    tempelev[i][j] = (int)e;

                    if (tempelev[i][j] < 0)
                        tempelev[i][j] = 0;

                    if (tempelev[i][j] < world.extraelev(i, j))
                        world.setextraelev(i, j, tempelev[i][j]);
                }
            }
        }
    });

    // Now smooth it

    world.smoothextraelev(2);
    world.smoothextraelev(5);
}

// Fill land depressions using Planchon-Darboux directional scans. The noisy
// minimum gradient prevents flat areas from blocking later river routing.

void depressionfill(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int e = 2; // This is the extra amount added to ensure that everywhere slopes.

    int neighbours[8][2];

    neighbours[0][0] = 0;
    neighbours[0][1] = -1;

    neighbours[1][0] = 1;
    neighbours[1][1] = -1;

    neighbours[2][0] = 1;
    neighbours[2][1] = 0;

    neighbours[3][0] = 1;
    neighbours[3][1] = 1;

    neighbours[4][0] = 0;
    neighbours[4][1] = 1;

    neighbours[5][0] = -1;
    neighbours[5][1] = 1;

    neighbours[6][0] = -1;
    neighbours[6][1] = -0;

    neighbours[7][0] = -1;
    neighbours[7][1] = -1;

    vector<vector<int>> noise(ARRAYWIDTH, vector<int>(ARRAYWIDTH, 0)); // This will contain noise that allows us to vary the value of e from tile to tile.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            noise[i][j] = random(0, 5);
    }

    vector<vector<int>> filledmap(ARRAYWIDTH, vector<int>(ARRAYWIDTH, 0)); // This will be the new version of the map.

    for (int i = 0; i <= width; i++) // First, fill the new map to a huge height.
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) <= sealevel)
                filledmap[i][j] = world.nom(i, j); // Sea tiles start with the normal heights.
            else
            {
                if (world.outline(i, j) == 1) // If this is a coastal tile
                    filledmap[i][j] = world.nom(i, j); // Coastal tiles start with the normal heights.
                else
                    filledmap[i][j] = maxelev * 2; // Other tiles start very high

            }
        }
    }

    

    // Now we're ready to start!

    twointegers lowest;

    lowest.x = 0;
    lowest.y = 0;

    bool somethingdone = 0;

    do
    {
        for (int i = 0; i <= width; i++) // Eight mini-loops, because there are eight possible ways to scan the map and we must rotate between them.
        {
            for (int j = 0; j <= height; j++)
            {
                if (filledmap[i][j] > world.nom(i, j))
                    somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int j = height; j >= 0; j--)
            {
                for (int i = width; i >= 0; i--)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int i = width; i >= 0; i--)
            {
                for (int j = 0; j <= height; j++)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int j = 0; j <= height; j++)
            {
                for (int i = width; i >= 0; i--)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int i = 0; i <= width; i++)
            {
                for (int j = height; j >= 0; j--)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int j = height; j >= 0; j--)
            {
                for (int i = 0; i <= width; i++)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int i = width; i > 0; i--)
            {
                for (int j = height; j >= 0; j--)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }

        if (somethingdone == 1)
        {
            somethingdone = 0;

            for (int j = 0; j <= height; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    if (filledmap[i][j] > world.nom(i, j))
                        somethingdone = checkdepressiontile(world, filledmap, i, j, e, neighbours, somethingdone, noise);

                }
            }
        }
    } while (somethingdone == 1);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            world.setnom(i, j, filledmap[i][j]);

    }
}

namespace
{
// This does the actual checking/filling of each tile.

bool checkdepressiontile(planet& world, vector<vector<int>>& filledmap, int i, int j, int e, int neighbours[8][2], bool somethingdone, vector<vector<int>>& noise)
{
    int width = world.width();
    int height = world.height();

    int start = random(0, 7);

    for (int n = start; n <= start + 7; n++)
    {
        int nn = wrap(n, 7);

        int k = i + neighbours[nn][0];

        if (k<0 || k>width)
            k = wrap(k, width);

        int l = j + neighbours[nn][1];

        if (l >= 0 && l <= height)
        {
            int ee = 0;

            if (k == i || l == j)
                ee = e + 1 + noise[i][j];

            else
                ee = e + noise[i][j];

            if (world.nom(i, j) >= filledmap[k][l] + ee)
            {
                filledmap[i][j] = world.nom(i, j);
                somethingdone = 1;

            }
            else
            {
                if (filledmap[i][j] > filledmap[k][l] + ee)
                {
                    filledmap[i][j] = filledmap[k][l] + ee;
                    somethingdone = 1;
                }
            }
        }
    }

    return (somethingdone);
}
}

// This function adds a bit of noise to the land.

void addlandnoise(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int maxadjust = 10;
    int adjustchance = 5; // The lower this is, the more there will be.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) > sealevel && random(1, adjustchance) == 1)
            {
                int adjust = randomsign(random(1, maxadjust));

                int newamount = world.nom(i, j) + adjust;

                if (newamount <= sealevel)
                    newamount = sealevel;

                if (newamount > maxelev)
                    newamount = maxelev;

                world.setnom(i, j, newamount);
            }
        }
    }
}
