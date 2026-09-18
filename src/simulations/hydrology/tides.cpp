#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "tides.hpp"

using namespace std;

namespace
{
int gettidalrange(planet& world, int startx, int starty);

// This function finds the tidal range of a given point.

int gettidalrange(planet& world, int startx, int starty)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    float mult1 = 0.01f;

    // We cast in each of the eight directions from the starting point until we hit land.
    // The longer each cast, the more tidal range it generates.

    int north = 0;
    int northeast = 0;
    int east = 0;
    int southeast = 0;
    int south = 0;
    int southwest = 0;
    int west = 0;
    int northwest = 0;

    bool hitland;
    int x, y;

    // North

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        y--;

        if (y < 0)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                north++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // South

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        y++;

        if (y > height)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                south++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // East

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x++;

        if (x > width)
            x = 0;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                east++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // West

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x--;

        if (x < 0)
            x = width;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                west++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // Northeast

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x++;

        if (x > width)
            x = 0;

        y--;

        if (y < 0)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                northeast++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // Southeast

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x++;

        if (x > width)
            x = 0;

        y++;

        if (y > height)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                southeast++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // Southwest

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x--;

        if (x < 0)
            x = width;

        y++;

        if (y > height)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                southwest++;
        }

        if (hitland == 1)
            n = 1000;
    }

    // Northwest

    hitland = 0;
    x = startx;
    y = starty;

    for (int n = 1; n <= 1000; n++)
    {
        x--;

        if (x < 0)
            x = width;

        y--;

        if (y < 0)
            hitland = 1;
        else
        {
            if (world.nom(x, y) > sealevel)
                hitland = 1;
            else
                northwest++;
        }

        if (hitland == 1)
            n = 1000;
    }

    float total = (float)(north + south + east + west + southeast + southwest + northeast + northwest);

    total = total * mult1;

    total = total * world.lunar(); // Lunar pull multiplies tides.

    int tide = (int)total; // In theory this could be as high as 80, but never in practice assuming lunar = 1.0. The highest in the world is 16 metres, so divide this value by 4 to get it in metres.

    return (tide);
}
}

// This function calculates the tides.

void createtidalmap(planet& world)
{
    int width = world.width();
    int height = world.height();
    int lunar = (int)world.lunar();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            world.settide(i, j, -1);
    }
    
    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1)
            {
                bool landfound = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (world.sea(kk, l) == 0)
                            {
                                landfound = 1;
                                k = i + 1;
                                l = j + 1;
                            }
                        }
                    }
                }

                if (landfound == 1)
                {
                    int tide = gettidalrange(world, i, j);

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 2; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.tide(kk, l) == -1)
                                    world.settide(kk, l, tide);
                            }
                        }
                    }

                    if (tide < lunar * 2)
                        tide = lunar * 2;

                    world.settide(i, j, tide);
                }
            }
        }
    }
}
