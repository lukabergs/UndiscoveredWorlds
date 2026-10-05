// Questions about coasts on the world grid: distance to land and sea, coastal cells, land in one quadrant.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "simulations/terrain/coast_queries.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <queue>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"

using namespace std;

// This finds the distance to the nearest land in the direction of the equator.

int landdistance(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();
    int equator = height / 2;

    if (x<0 || x>width || y<0 || y>height)
        return (-1);

    int dist = 0;
    int pets = -1;

    int target = equator;

    if (y < equator)
        pets = 1;

    // Reverse the direction of search if the equator is colder than the poles!

    if (y < equator) // Northern hemisphere
    {
        if (world.northpolartemp()>world.eqtemp())
        {
            pets = -1;
            target = 0;
        }
    }
    else // Southern hemisphere
    {
        if (world.southpolartemp() > world.eqtemp())
        {
            pets = 1;
            target = height;
        }
    }

    for (int j = y; j != target; j = j + pets)
    {
        if (world.sea(x, j) == 0)
        {
            if (world.sea(x, j + pets) == 0 && world.sea(x, j + pets * 2) == 0) // To avoid little islands showing up
            {
                int xleft = x - 1;

                if (xleft == -1)
                    xleft = width;

                int xright = x + 1;
                if (xright == width + 1)
                    xright = 0;

                if (world.sea(xleft, j + pets) == 0 && world.sea(xright, j + pets) == 0)
                {
                    return (dist);
                }
            }
        }

        dist++;
    }

    dist = -1; // No land found at all.

    return(dist);
}

// This finds the coordinates of the nearest sea (or land, if land=1) to the point specified. Note that it does not wrap the x coordinate, i.e. it could return one lower than 0 or higher than width.

twointegers nearestsea(planet& world, int i, int j, bool land, int limit, int grain)
{
    int width = world.width();
    int height = world.height();

    twointegers seapoint;

    seapoint.x = -1;
    seapoint.y = -1;

    int distance = -1; // Distance to nearest sea.
    int seax = -1; // x coordinate of nearest sea.
    int seay = -1; // y coordinate of nearest sea.
    int checkdist = 1; // Distance that we're checking to see if it contains sea.

    while (distance == -1)
    {
        checkdist = checkdist + grain; // Increase the size of circle that we're checking

        int rr = checkdist * checkdist;

        for (int x = i - checkdist; x <= i + checkdist; x++)
        {
            int xx = x;

            if (xx<0 || xx>width)
                xx = wrap(xx, width);

            for (int y = j - checkdist; y <= j + checkdist; y++)
            {
                if (y >= 0 && y <= height)
                {
                    if ((x-i)*(x-i)+(y-j)*(y-j)==rr) // if this point is on the circle we're looking at
                    {
                        if (world.sea(xx, y) == 1 && land == 0) // And if this point is sea
                        {
                            distance = checkdist;
                            seax = x;
                            seay = y;
                        }

                        if (world.sea(xx, y) == 0 && land == 1) // And if this point is land
                        {
                            distance = checkdist;
                            seax = x;
                            seay = y;
                        }
                    }
                }
            }
        }
        if (checkdist > limit) // We are surely not going to find any now...
        {
            distance = width * 2;
            seax = -1;
            seay = -1;
        }
    }

    seapoint.x = seax;
    seapoint.y = seay;

    return (seapoint);
}

// This function checks whether a global tile is either sea next to land or land next to sea.

bool vaguelycoastal(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    bool sea = 0;

    if (world.sea(x, y) == 1)
        sea = 1;

    bool found = 0;

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (j >= 0 && j <= height)
            {
                if (sea == 0 && world.sea(ii, j) == 1)
                {
                    found = 1;
                    i = x + 1;
                    j = y + 1;
                }

                if (sea == 1 && world.sea(ii, j) == 0)
                {
                    found = 1;
                    i = x + 1;
                    j = y + 1;
                }
            }
        }
    }

    return found;
}

// This function checks to see if a global tile is land, bordered by sea to the east, south, and southeast.

bool northwestlandonly(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    if (world.sea(x, y) == 1)
        return 0;

    int xx = x + 1;

    if (xx > width)
        xx = 0;

    int yy = y + 1;

    if (yy > height)
        return 0;

    if (world.sea(xx, y) == 0)
        return 0;

    if (world.sea(x, yy) == 0)
        return 0;

    if (world.sea(xx, yy) == 0)
        return 0;

    xx = x - 1;

    if (xx < 0)
        x = width;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (world.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// This function checks to see if a global tile is land, bordered by sea to the west, south, and southwest.

bool northeastlandonly(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    if (world.sea(x, y) == 1)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        xx = width;

    int yy = y + 1;

    if (yy > height)
        return 0;

    if (world.sea(xx, y) == 0)
        return 0;

    if (world.sea(x, yy) == 0)
        return 0;

    if (world.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > width)
        xx = 0;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (world.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// This function checks to see if a global tile is land, bordered by sea to the east, north, and northeast.

bool southwestlandonly(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    if (world.sea(x, y) == 1)
        return 0;

    int xx = x + 1;

    if (xx > width)
        xx = 0;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (world.sea(xx, y) == 0)
        return 0;

    if (world.sea(x, yy) == 0)
        return 0;

    if (world.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > width)
        x = 0;

    yy = y + 1;

    if (yy > height)
        return 0;

    if (world.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// This function checks to see if a global tile is land, bordered by sea to the west, north, and northwest.

bool southeastlandonly(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    if (world.sea(x, y) == 1)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        xx = width;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (world.sea(xx, y) == 0)
        return 0;

    if (world.sea(x, yy) == 0)
        return 0;

    if (world.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > width)
        xx = 0;

    yy = y + 1;

    if (yy > height)
        return 0;

    if (world.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// This does the same thing, for the continental shelves.

void makeshelvesvoronoi(planet& world, vector<vector<short>>& voronoi, vector<vector<bool>>& outline, int pointdist)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int midelev = maxelev / 2;

    int grain = 16; // Level of detail on this fractal map.
    float valuemod = 0.4f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    int warpdist = 1;
    int warpdiv = midelev / warpdist;

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

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int closest = 0;
            int dist = 10000000;

            int warpdistx = 0;
            int warpdisty = 0;

            int warpx = i;
            int warpy = j;

            if (fractal[warpx][warpy] > maxelev / 2)
                warpdistx = (fractal[warpx][warpy] - midelev) / warpdiv;
            else
                warpdistx = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

            int checkx = i + warpdistx;

            if (checkx<0 || checkx>width)
                wrap(checkx, width);

            warpx = warpx + width / 2;

            if (warpx > width)
                warpx = warpx - width;

            if (fractal[warpx][warpy] > maxelev / 2)
                warpdisty = (fractal[warpx][warpy] - midelev) / warpdiv;
            else
                warpdisty = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

            int checky = j + warpdisty;

            if (checky < 0)
                checky = 0;

            if (checky > height)
                checky = height;

            checkx = i;
            checky = j;

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
                                int pointdist = xdist * xdist + ydist * ydist;

                                if (pointdist < dist)
                                {
                                    closest = pointsmap[xx][y];
                                    dist = pointdist;
                                }
                            }
                        }
                    }
                }

                distcheck = distcheck + 10;

            } while (closest == 0);

            voronoi[i][j] = closest;
        }
    }
}
