// Isolated volcanoes and mountainous islands.
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

// This makes a small, mountainous island.

void terrain::detail::makemountainisland(planet& world, int x, int y, int peakheight)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int iterations = random(1, 3);

    for (int n = 1; n <= iterations; n++)
    {
        int landheight = sealevel + random(5, 50);

        int xx = x;
        int yy = y;

        do
        {
            xx = x - 2 + random(1, 3);
            yy = y - 2 + random(1, 3);

        } while (x == xx && y == yy);

        if (yy<0 || yy>height)
            return;

        if (xx<0 || xx>width)
            xx = wrap(xx, width);

        if (world.mountainheight(x, y) < peakheight)
            world.setmountainheight(x, y, peakheight);

        if (world.mountainheight(xx, yy) < peakheight)
            world.setmountainheight(xx, yy, peakheight);

        for (int i = x - 1; i <= x + 1; i++) // Mark it all out as mountain islands!
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = y - 1; j <= y + 1; j++)
            {
                if (j >= 0 && j <= height)
                    world.setmountainisland(i, j, 1);
            }
        }

        for (int i = xx - 1; i <= xx + 1; i++)
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = yy - 1; j <= yy + 1; j++)
            {
                if (j >= 0 && j <= height)
                    world.setmountainisland(i, j, 1);
            }
        }

        if (world.nom(x, y) < landheight)
            world.setnom(x, y, landheight);

        if (world.nom(xx, yy) < landheight)
            world.setnom(xx, yy, landheight);

        // If there isn't a ridge between them, make one.

        int dir = getdir(x, y, xx, yy);

        if (getridge(world, x, y, dir) == 0)
        {
            int code = terrain::detail::getcode(dir);

            world.setmountainridge(x, y, world.mountainridge(x, y) + code);

            dir = dir + 8;

            if (dir > 8)
                dir = dir - 8;

            code = terrain::detail::getcode(dir);

            world.setmountainridge(xx, yy, world.mountainridge(xx, yy) + code);
        }

        //world.settest(x,y,1);

        x = xx;
        y = yy;

        peakheight = peakheight + randomsign(random(peakheight / 4, peakheight / 2));
    }
}

// This makes an isolated volcano, together with a line of extinct neighbours.

void createisolatedvolcano(planet& world, int x, int y, vector<vector<bool>>& shelves, vector<vector<int>>& volcanodirection, int peakheight, bool strato)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    float maxshift = 6;
    float shiftdiv = maxshift / (maxelev / 2);

    int shelfcheck = 5;

    if (world.sea(x, y) == 1) // Don't do any submarine volcanoes near continental shelves.
    {
        for (int i = x - shelfcheck; i <= x + shelfcheck; i++)
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = y - shelfcheck; j <= y + shelfcheck; j++)
            {
                if (j >= 0 && j <= height)
                {
                    if (shelves[ii][j] == 1)
                        return;
                }
            }
        }
    }

    bool active = 1;

    int total = random(1, 6);

    for (int n = 1; n != total; n++)
    {
        if (world.sea(x, y) == 1)
        {
            for (int i = x - 1; i <= x + 1; i++) // No shading around undersea volcanoes.
            {
                int ii = i;

                if (ii<0 || ii>width)
                    ii = wrap(ii, width);

                for (int j = y - 1; j <= y + 1; j++)
                {
                    if (j >= 0 && j <= height)
                        world.setnoshade(i, j, 1);
                }
            }
        }

        float thispeakheight = (float)peakheight;

        thispeakheight = thispeakheight / 100.0f;

        thispeakheight = thispeakheight * (float)random(80, 120);

        if (active == 0)
            thispeakheight = (float)(0 - peakheight);

        for (int i = x - 1; i <= x + 1; i++) // Remove any volcanoes immediately around
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = y - 1; j <= y + 1; j++)
            {
                if (j >= 0 && j <= height)
                {
                    world.setvolcano(ii, j, 0);
                    //world.settest(ii,j,0);
                }
            }
        }

        world.setvolcano(x, y, (int)thispeakheight);
        world.setstrato(x, y, strato);

        active = 0;

        int xx = x + width / 2;

        if (xx > width)
            xx = xx - width;

        float xshift = (float)(volcanodirection[x][y] - maxelev / 2);
        float yshift = (float)(volcanodirection[xx][y] - maxelev / 2);

        xshift = xshift * shiftdiv + (float)randomsign(random(0, 2));
        yshift = yshift * shiftdiv + (float)randomsign(random(0, 2));

        x = x + (int)xshift;

        if (x<0 || x>width)
            x = wrap(x, width);

        y = y + (int)yshift;

        if (y<0 || y>height)
            return;

        if (world.sea(x, y) == 1) // Don't do any submarine volcanoes near continental shelves.
        {
            for (int i = x - shelfcheck; i <= x + shelfcheck; i++)
            {
                int ii = i;

                if (ii<0 || ii>width)
                    ii = wrap(ii, width);

                for (int j = y - shelfcheck; j <= y + shelfcheck; j++)
                {
                    if (j >= 0 && j <= height)
                    {
                        if (shelves[ii][j] == 1)
                            return;
                    }
                }
            }
        }
    }
}
