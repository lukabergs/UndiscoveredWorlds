// Questions about drainage on the world grid: slopes, flow destinations, inflows, lakes nearby.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "simulations/hydrology/flow_queries.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <queue>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"

using namespace std;

// This function returns the direction from one tile to another.

// This function finds the direction of the lowest neighbouring tile.

int findlowestdir(planet& world, int neighbours[8][2], int x, int y)
{
    int width = world.width();
    int height = world.height();
    const std::uint64_t directionseed = deterministiccontextseed(world.seed(), 0x6101);

    float lowest = static_cast<float>(world.maxelevation() * 2);
    int nnn = -1;

    int start = deterministicrandom(directionseed, 0, 7, x, y);

    for (int n = start; n <= start + 7; n++)
    {
        int nn = wrap(n, 7);

        int i = x + neighbours[nn][0];

        if (i<0 || i>width)
            i = wrap(i, width);

        int j = y + neighbours[nn][1];

        if (j >= 0 && j <= height)
        {
            float pointelev = world.nomfloat(i, j);

            if (pointelev < lowest)
            {
                lowest = pointelev;
                nnn = nn;
            }
        }
    }

    int dir = nnn + 1;

    return (dir);
}

// This gets the amount the land is sloping between two points.

float getslope(planet& world, int x, int y, int xx, int yy)
{
    int width = world.width();
    int height = world.height();

    if (x<0 || x>width)
        x = wrap(x, width);

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (y<0 || y>height)
        return (-1.0f);

    if (yy<0 || yy>height)
        return (-1.0f);

    float slope = world.mapfloat(xx, yy) - world.mapfloat(x, y);

    return (slope);
}

// This function checks how flat the land is at a given point.

int getflatness(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    int xminus = x - 1;

    if (xminus < 0)
        xminus = width;

    int xplus = x + 1;

    if (xplus > width)
        xplus = 0;

    const float thiselev = getflatelevation(world, x, y);

    float flatness = 0.0f;

    float thatelev = getflatelevation(world, xminus, y);
    flatness = flatness + std::fabs(thiselev - thatelev);

    float crount = 2.0f;

    if (y > 0)
    {
        thatelev = getflatelevation(world, xminus, y - 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        thatelev = getflatelevation(world, x, y - 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        thatelev = getflatelevation(world, xplus, y - 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        crount = crount + 3.0f;
    }

    if (y < height)
    {
        thatelev = getflatelevation(world, xminus, y + 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        thatelev = getflatelevation(world, x, y + 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        thatelev = getflatelevation(world, xplus, y + 1);
        flatness = flatness + std::fabs(thiselev - thatelev);

        crount = crount + 3.0f;
    }

    flatness = flatness / crount;

    return (int)flatness;
}

// This function gets the elevation of a point, incorporating water level.

float getflatelevation(planet& world, int x, int y)
{
    float elev;

    float surface = world.lakesurfacefloat(x, y);

    if (surface > 0)
        elev = surface;
    else
    {
        if (world.sea(x, y) == 1)
            elev = static_cast<float>(world.sealevel());
        else
            elev = world.nomfloat(x, y);
    }

    return (elev);
}

// This function finds a neighbouring sea tile that isn't the one we just came from.

twointegers findseatile(planet& world, int x, int y, int dir)
{
    int width = world.width();
    int height = world.height();

    twointegers newtile;

    newtile.x = -1;
    newtile.y = -1;

    int oldx = x;
    int oldy = y;

    if (dir == 8 || dir == 1 || dir == 2)
        oldy = y + 1;

    if (dir == 4 || dir == 5 || dir == 6)
        oldy = y - 1;

    if (dir == 2 || dir == 3 || dir == 4)
        oldx = x - 1;

    if (dir == 6 || dir == 7 || dir == 8)
        oldx = x + 1;

    if (oldx<0 || oldx>width)
        oldx = wrap(oldx, width);

    // oldx and oldy are the coordinates of where we are coming to this tile from. We want to avoid going back there!

    float lowest = world.mapfloat(x, y);

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (ii != x && j != y && ii != oldx && j != oldy && j >= 0 && j <= height && world.sea(ii, j) == 1 && world.mapfloat(ii, j) < lowest)
            {
                newtile.x = ii;
                newtile.y = j;
                lowest = world.mapfloat(ii, j);
            }
        }
    }

    if (newtile.x != -1)
        return (newtile);

    // If it didn't find a neighbouring sea tile that's lower, just find one at all.

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (ii != x && j != y && ii != oldx && j != oldy && j >= 0 && j <= height && world.sea(ii, j) == 1)
            {
                newtile.x = ii;
                newtile.y = j;
            }
        }
    }
    return (newtile);
}

// This returns the coordinate that the river in the given tile is flowing into.

twointegers getflowdestination(planet& world, int x, int y, int dir)
{
    int width = world.width();
    int height = world.height();

    twointegers destpoint;

    if (dir == 0)
        dir = world.riverdir(x, y);

    if (dir == 8 || dir == 1 || dir == 2)
        y--;

    if (dir == 4 || dir == 5 || dir == 6)
        y++;

    if (dir == 2 || dir == 3 || dir == 4)
        x++;

    if (dir == 6 || dir == 7 || dir == 8)
        x--;

    if (x<0 || x>width)
        x = wrap(x, width);

    if (y < 0)
        y = 0;

    if (y > height)
        y = height;

    destpoint.x = x;
    destpoint.y = y;

    return (destpoint);
}

// This function finds the cell with the largest flow into the current one.

twointegers getupstreamcell(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    twointegers upcell, destpoint;

    upcell.x = -1;
    upcell.y = -1;

    int largest = 0;

    int riverjan, riverjul;

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (j >= 0 && y <= height)
            {
                riverjan = world.riverjan(ii, j);
                riverjul = world.riverjul(ii, j);

                if (riverjan + riverjul > largest) // If this is larger than the current largest inflow
                {
                    destpoint = getflowdestination(world, ii, j, 0);

                    if (destpoint.x == x && destpoint.y == y) // If this is actually flowing into our cell
                    {
                        upcell.x = ii;
                        upcell.y = j;

                        largest = riverjan + riverjul;
                    }

                }
            }
        }
    }
    return (upcell);
}

// This function checks to see whether a tile on the global map has any water flowing into it.

int checkwaterinflow(planet& world, int x, int y)
{
    int i, j;

    int width = world.width();
    int height = world.height();

    // From the north

    i = x;
    j = y - 1;

    if (j >= 0 && world.riverdir(i, j) == 5)
        return 1;

    // From the northeast

    i = x + 1;
    if (i > width)
        i = 0;
    j = y - 1;

    if (j >= 0 && world.riverdir(i, j) == 6)
        return 1;

    // From the east

    i = x + 1;
    if (i > width)
        i = 0;
    j = y;

    if (j >= 0 && world.riverdir(i, j) == 7)
        return 1;

    // From the southeast

    i = x + 1;
    if (i > width)
        i = 0;
    j = y + 1;

    if (j <= height && world.riverdir(i, j) == 8)
        return 1;

    // From the south

    i = x;
    j = y + 1;

    if (j <= height && world.riverdir(i, j) == 1)
        return 1;

    // From the southwest

    i = x - 1;
    if (i < 0)
        i = width;
    j = y + 1;

    if (j <= height && world.riverdir(i, j) == 2)
        return 1;

    // From the west

    i = x - 1;
    if (i < 0)
        i = width;
    j = y;

    if (j <= height && world.riverdir(i, j) == 3)
        return 1;

    // From the northwest

    i = x - 1;
    if (i < 0)
        i = width;
    j = y - 1;

    if (j >= 0 && world.riverdir(i, j) == 4)
        return 1;

    return 0;
}

// This function gets the total water inflow into a given tile on the global map.

twointegers gettotalinflow(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    twointegers total;
    twointegers dest;

    total.x = 0;
    total.y = 0;

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (j >= 0 && j <= height)
            {
                dest = getflowdestination(world, i, j, 0);

                if (dest.x == x && dest.y == y)
                {
                    total.x = total.x + world.riverjan(i, j);
                    total.y = total.y + world.riverjul(i, j);
                }
            }
        }
    }
    return total;
}

// This function tells whether a given point on the world map is in or next to a lake.

float nearlake(planet& world, int x, int y, int dist, bool rift)
{
    int width = world.width();
    int height = world.height();

    for (int i = x - dist; i <= x + dist; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - dist; j <= y + dist; j++)
        {
            if (j >= 0 && j <= height)
            {
                if (world.lakesurfacefloat(ii, j) != 0)
                    return world.lakesurfacefloat(ii, j);

                if (rift == 1 && world.riftlakesurfacefloat(ii, j) != 0)
                    return world.riftlakesurfacefloat(ii, j);
            }
        }
    }
    return 0.0f;
}

// This function tells us whether a tile is a lake tile on the edge.

int getlakeedge(planet& world, int x, int y)
{
    if (world.lakesurfacefloat(x, y) == 0)
        return 0;

    int width = world.width();
    int height = world.height();

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (j >= 0 && j <= height)
            {
                if (world.lakesurfacefloat(ii, j) == 0)
                    return 1;
            }
        }
    }
    return 0;
}
