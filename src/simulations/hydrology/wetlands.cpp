#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "wetlands.hpp"

using namespace std;

namespace
{
void pastewetland(planet& world, int centrex, int centrey, int shapenumber, boolshapetemplate smalllake[]);

// This actually pastes a patch of wetlands onto the map.

void pastewetland(planet& world, int centrex, int centrey, int shapenumber, boolshapetemplate smalllake[])
{
    int width = world.width();
    int height = world.height();

    if (centrey == 0 || centrey == height)
        return;

    int imheight = smalllake[shapenumber].ysize() - 1;
    int imwidth = smalllake[shapenumber].xsize() - 1;

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

    int leftx = centrex;
    int lefty = centrey;
    int rightx = centrex;
    int righty = centrey; // These are the coordinates of the furthermost pixels of the swamp.
    bool wrapped = 0; // If this is 1, the swamp wraps over the edge of the map.

    int leftr = random(0, 1); // If it's 1 then we reverse it left-right
    int downr = random(0, 1); // If it's 1 then we reverse it top-bottom

    int istart = 0, desti = imwidth + 1, istep = 1;
    int jstart = 0, destj = imheight + 1, jstep = 1;

    if (leftr == 1)
    {
        istart = imwidth;
        desti = -1;
        istep = -1;
    }

    if (downr == 1)
    {
        jstart = imwidth;
        destj = -1;
        jstep = -1;
    }

    for (int i = istart; i != desti; i = i + istep)
    {
        for (int j = jstart; j != destj; j = j + jstep)
        {
            if (smalllake[shapenumber].point(i, j) == 1)
            {
                int xx = x + i;
                int yy = y + j;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                    {
                        xx = wrap(xx, width);
                        wrapped = 1;
                    }

                    if (world.sea(xx, yy) == 0 && world.mountainheight(xx, yy) == 0) // Don't try to put wetlands on top of sea or mountains...
                    {
                        world.setspecial(xx, yy, 130);

                        if (xx < leftx)
                            leftx = xx;

                        if (xx > rightx)
                            rightx = xx;

                        if (yy < lefty)
                            lefty = yy;

                        if (yy > righty)
                            righty = yy;

                    }
                }
            }
        }
    }

    if (wrapped == 1)
    {
        leftx = 0;
        lefty = 0;
        rightx = width;
        righty = height;
    }
}
}

// This lays down wetlands.

void createwetlands(planet& world, boolshapetemplate smalllake[])
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    int minflow = 100; // Average flow must be at least this to have wetlands.
    int maxflatness = 5; // Must be no less flat than this.
    float factor = (float)maxelev / 5.0f;

    // First, create a fractal to be the drainage map.

    int grain = 8; // This is the level of detail on this map.
    float valuemod = 0.2f;

    int v = random(3, 6);
    float valuemod2 = (float)v;
    int shapenumber;

    vector<vector<int>> drainage(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(drainage, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    // Now go through the map and place wetlands.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 0 && world.special(i, j) == 0 && world.lakesurface(i, j) == 0 && world.climate(i, j) != 31 && world.climate(i, j) != 5)
            {
                if (world.deltadir(i, j) != 0)
                {
                    shapenumber = random(0, 5);
                    pastewetland(world, i, j, shapenumber, smalllake);
                }
                else
                {
                    if (world.riveraveflow(i, j) >= minflow)
                    {
                        int flatness = getflatness(world, i, j);

                        if (flatness <= maxflatness)
                        {
                            float drain = (float)drainage[i][j];

                            drain = drain / factor;

                            int d = (int)drain; // This should be from 1 to 5.
                            d = d * 2;

                            int prob = flatness + d; // The flatter the land, the more likely wetlands are.

                            int flow = world.riveraveflow(i, j) / 1000;

                            prob = prob - flow; // The bigger the flow here, the more likely wetlands are.

                            for (int k = i - 1; k <= i + 1; k++) // The more lakes/sea there are around here, the more likely wetlands are.
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - 1; l <= j + 1; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.lakesurface(kk, l) != 0)
                                            prob = prob - 3;

                                        if (world.sea(kk, l) == 1)
                                            prob = prob - 4;
                                    }
                                }
                            }

                            //prob=prob+(7-(world.averain(i,j)/60));
                            prob = prob + static_cast<int>(
                                40.0f - world.averainfloat(i, j) / 10.0f);

                            if (world.climate(i, j) == 30 && world.maxtemp(i, j) >= 5) // Much likelier in tundra with warmish summers
                                prob = prob / 4;

                            if (prob < 1)
                                prob = 1;

                            if (random(1, prob) == 1) // Put some wetlands here.
                            {
                                shapenumber = random(0, 5);
                                pastewetland(world, i, j, shapenumber, smalllake);
                            }
                        }
                    }
                }
            }
        }
    }

    int deltaswampchance = 6; // The BIGGER this is, the more probable wetlands will be around delta branches.

    for (int i = 0; i <= width; i++) // Add wetlands to the actual delta branches
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.deltadir(i, j) > 0 && world.sea(i, j) == 0 && random(1, deltaswampchance) != 1)
                world.setspecial(i, j, 130);

        }
    }

    // Now fill in extra wetlands if necessary to reach the sea or other wetlands.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 0 && world.special(i, j) == 0 && world.climate(i, j) != 31)
            {
                int westx = i - 1;
                if (westx < 0)
                    westx = width;

                int eastx = i + 1;
                if (eastx > width)
                    eastx = 0;

                int northy = j - 1;
                int southy = j + 1;

                if (world.special(i, northy) == 130)
                {
                    if (world.special(i, southy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(i, southy) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(i, southy) == 130)
                {
                    if (world.special(i, northy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(i, northy) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(westx, j) == 130)
                {
                    if (world.special(eastx, j) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(eastx, j) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(eastx, j) == 130)
                {
                    if (world.special(westx, j) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(westx, j) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(westx, northy) == 130)
                {
                    if (world.special(eastx, southy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(eastx, southy) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(eastx, northy) == 130)
                {
                    if (world.special(westx, southy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(westx, southy) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(westx, southy) == 130)
                {
                    if (world.special(eastx, northy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(eastx, northy) == 1)
                        world.setspecial(i, j, 130);
                }

                if (world.special(eastx, southy) == 130)
                {
                    if (world.special(westx, northy) == 130)
                        world.setspecial(i, j, 130);

                    if (world.sea(westx, northy) == 1)
                        world.setspecial(i, j, 130);
                }
            }
        }
    }

    for (int i = 0; i <= width; i++) // For some reason it keeps putting wetlands all over the sea, so we remove them.
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1)
                world.setspecial(i, j, 0);

        }
    }

    // Now make the wetlands salty/brackish if necessary.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.special(i, j) == 130)
            {
                int salt = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (world.sea(kk, l) == 1)
                                salt = 1;

                            if (world.special(kk, l) == 100 || world.special(kk, l) == 110)
                                salt = 1;
                        }
                    }
                }

                if (salt == 1)
                    world.setspecial(i, j, 132);
                else
                {
                    for (int k = i - 2; k <= i + 2; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 2; l <= j + 2; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.sea(kk, l) == 1)
                                    salt = 1;

                                if (world.special(kk, l) == 100 || world.special(kk, l) == 110)
                                    salt = 1;
                            }
                        }
                    }

                    if (salt == 1)
                        world.setspecial(i, j, 131);
                }
            }
        }
    }
}

// This removes any wetlands tiles that border large lakes.

void removeexcesswetlands(planet& world)
{
    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.special(i, j) == 130)
            {
                bool nearlake = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (world.truelake(kk, l) == 1)
                            {
                                nearlake = 1;
                                k = i + 1;
                                l = j + 1;
                            }
                        }
                    }
                }

                if (nearlake == 1)
                {
                    world.setlakesurface(i, j, 0);
                    world.setspecial(i, j, 0);
                }
            }
        }
    }
}
