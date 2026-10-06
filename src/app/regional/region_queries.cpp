// The same questions on a regional map (coasts, flow, lakes, rivers), and temperature adjusted for regional elevation.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "app/regional/region_queries.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <queue>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "app/regional/region.hpp"
#include "app/regional/regional_map.hpp"

using namespace std;

// This function adds the elevation element to a temperature.

// Same thing, but at the regional level.

int tempelevadd(planet& world, region& region, int temp, int i, int j)
{
    int elevation = region.map(i, j) - region.sealevel();

    if (elevation < 0)
        elevation = 0;

    float elevationadjust = (float)elevation;
    elevationadjust = elevationadjust / 1000.0f;
    elevationadjust = elevationadjust * world.tempdecrease();

    temp = temp - (int)elevationadjust;

    return temp;
}

// This function removes the elevation element of a temperature.

// Same thing, but at the regional level.

int tempelevremove(planet& world, region& region, int temp, int i, int j)
{
    int elevation = region.map(i, j) - region.sealevel();

    if (elevation < 0)
        elevation = 0;

    float elevationadjust = (float)elevation;
    elevationadjust = elevationadjust / 1000.0f;
    elevationadjust = elevationadjust * world.tempdecrease();

    temp = temp + (int)elevationadjust;

    return temp;
}

// The same on the global map, in whole numbers as the regional generator computes (the float form is
// tempelevremove in simulations/climate/climate_elevation.hpp).

int tempelevremove(planet& world, int temp, int i, int j)
{
    int elevation = world.map(i, j) - world.sealevel();

    if (elevation < 0)
        elevation = 0;

    float elevationadjust = (float)elevation;
    elevationadjust = elevationadjust / 1000.0f;
    elevationadjust = elevationadjust * world.tempdecrease();

    temp = temp + (int)elevationadjust;

    return temp;
}

// Same thing on the regional map.

int getslope(region& region, int x, int y, int xx, int yy)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (x<0 || x>rwidth)
        return (-1);

    if (xx<0 || xx>rwidth)
        return (-1);

    if (y<0 || y>rheight)
        return (-1);

    if (yy<0 || yy>rheight)
        return (-1);

    int slope = region.map(xx, yy) - region.map(x, y);

    return (slope);
}

// This does the same thing, but for the regional map.

twointegers nearestsea(region& region, int leftx, int lefty, int rightx, int righty, int i, int j)
{
    twointegers seapoint;

    seapoint.x = -1;
    seapoint.y = -1;

    int distance = -1; // Distance to nearest sea.
    int seax = -1; // x coordinate of nearest sea.
    int seay = -1; // y coordinate of nearest sea.
    int checkdist = 1; // Distance that we're checking to see if it contains sea.

    while (distance == -1)
    {
        checkdist++; // Increase the size of circle that we're checking

        int rr = checkdist * checkdist;

        for (int x = i - checkdist; x <= i + checkdist; x++)
        {
            for (int y = j - checkdist; y <= j + checkdist; y++)
            {
                if ((x - i) * (x - i) + (y - j) * (y - j) == rr) // if this point is on the circle we're looking at
                {
                    if (x >= leftx && x <= rightx && y >= lefty && y <= righty)
                    {
                        if (region.sea(x, y) == 1) // And if this point is sea
                        {
                            distance = checkdist;
                            seax = x;
                            seay = y;
                        }
                        else
                        {
                            if (region.lakesurface(x, y) != 0) // Lake counts as sea for our purposes.
                            {
                                distance = checkdist;
                                seax = x;
                                seay = y;
                            }
                        }
                    }
                }
            }
        }
        if (checkdist > 100) // We are surely not going to find any now...
        {
            distance = region.rwidth() * 2;
            seax = -1;
            seay = -1;
        }
    }

    seapoint.x = seax;
    seapoint.y = seay;

    return seapoint;
}

// Same thing, for a regional cell.

bool vaguelycoastal(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    bool sea = 0;

    if (region.sea(x, y) == 1)
        sea = 1;

    bool found = 0;

    for (int i = x - 1; i <= x + 1; i++)
    {
        for (int j = y - 1; j <= y + 1; j++)
        {
            if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight)
            {
                if (sea == 0 && region.sea(i, j) == 1)
                {
                    found = 1;
                    i = x + 1;
                    j = y + 1;
                }

                if (sea == 1 && region.sea(i, j) == 0)
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

// Same thing, at the regional level.

bool northwestlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.sea(x, y) == 1)
        return 0;

    int xx = x + 1;

    if (xx > rwidth)
        return 0;

    int yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.sea(xx, y) == 0)
        return 0;

    if (region.sea(x, yy) == 0)
        return 0;

    if (region.sea(xx, yy) == 0)
        return 0;

    xx = x - 1;

    if (xx < 0)
        return 0;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// Same again, but for lakes rather than sea.

bool lakenorthwestlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.lakesurface(x, y) != 0)
        return 0;

    int xx = x + 1;

    if (xx > rwidth)
        return 0;

    int yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.lakesurface(xx, y) == 0)
        return 0;

    if (region.lakesurface(x, yy) == 0)
        return 0;

    if (region.lakesurface(xx, yy) == 0)
        return 0;

    xx = x - 1;

    if (xx < 0)
        return 0;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.lakesurface(xx, yy) != 0)
        return 0;

    return 1;
}

// Same thing, at the regional level.

bool northeastlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.sea(x, y) == 1)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        return 0;

    int yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.sea(xx, y) == 0)
        return 0;

    if (region.sea(x, yy) == 0)
        return 0;

    if (region.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// Same again, but for lakes rather than sea.

bool lakenortheastlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.lakesurface(x, y) != 0)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        return 0;

    int yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.lakesurface(xx, y) == 0)
        return 0;

    if (region.lakesurface(x, yy) == 0)
        return 0;

    if (region.lakesurface(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.lakesurface(xx, yy) != 0)
        return 0;

    return 1;
}

// Same thing, at the regional level.

bool southwestlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.sea(x, y) == 1)
        return 0;

    int xx = x + 1;

    if (xx > rwidth)
        return 0;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.sea(xx, y) == 0)
        return 0;

    if (region.sea(x, yy) == 0)
        return 0;

    if (region.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// Same thing, but for lakes rather than sea.

bool lakesouthwestlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.lakesurface(x, y) != 0)
        return 0;

    int xx = x + 1;

    if (xx > rwidth)
        return 0;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.lakesurface(xx, y) == 0)
        return 0;

    if (region.lakesurface(x, yy) == 0)
        return 0;

    if (region.lakesurface(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.lakesurface(xx, yy) != 0)
        return 0;

    return 1;
}

// Same thing, at the regional level.

bool southeastlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.sea(x, y) == 1)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        return 0;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.sea(xx, y) == 0)
        return 0;

    if (region.sea(x, yy) == 0)
        return 0;

    if (region.sea(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.sea(xx, yy) == 1)
        return 0;

    return 1;
}

// Same thing, but for lakes rather than sea..

bool lakesoutheastlandonly(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    if (region.lakesurface(x, y) != 0)
        return 0;

    int xx = x - 1;

    if (xx < 0)
        return 0;

    int yy = y - 1;

    if (yy < 0)
        return 0;

    if (region.lakesurface(xx, y) == 0)
        return 0;

    if (region.lakesurface(x, yy) == 0)
        return 0;

    if (region.lakesurface(xx, yy) == 0)
        return 0;

    xx = x + 1;

    if (xx > rwidth)
        return 0;

    yy = y + 1;

    if (yy > rheight)
        return 0;

    if (region.lakesurface(xx, yy) != 0)
        return 0;

    return 1;
}

// The same thing, but for the regional map.

twointegers getregionalflowdestination(region& region, int x, int y, int dir)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    twointegers destpoint;

    if (dir == 0)
        dir = region.riverdir(x, y);

    if (dir == 8 || dir == 1 || dir == 2)
        y--;

    if (dir == 4 || dir == 5 || dir == 6)
        y++;

    if (dir == 2 || dir == 3 || dir == 4)
        x++;

    if (dir == 6 || dir == 7 || dir == 8)
        x--;

    if (x<0 || x>rwidth || y<0 || y>rheight)
    {
        destpoint.x = -1;
        destpoint.y = -1;
    }
    else
    {
        destpoint.x = x;
        destpoint.y = y;
    }

    return destpoint;
}

// This does the same thing for the regional map.

int checkregionalwaterinflow(region& region, int x, int y)
{
    int i, j;

    // From the north

    i = x;
    j = y - 1;

    if (region.riverdir(i, j) == 5)
        return 1;

    // From the northeast

    i = x + 1;
    j = y - 1;

    if (region.riverdir(i, j) == 6)
        return 1;

    // From the east

    i = x + 1;
    j = y;

    if (region.riverdir(i, j) == 7)
        return 1;

    // From the southeast

    i = x + 1;
    j = y + 1;

    if (region.riverdir(i, j) == 8)
        return 1;

    // From the south

    i = x;
    j = y + 1;

    if (region.riverdir(i, j) == 1)
        return 1;

    // From the southwest

    i = x - 1;
    j = y + 1;

    if (region.riverdir(i, j) == 2)
        return 1;

    // From the west

    i = x - 1;
    j = y;

    if (region.riverdir(i, j) == 3)
        return 1;

    // From the northwest

    i = x - 1;
    j = y - 1;

    if (region.riverdir(i, j) == 4)
        return 1;

    return 0;
}

// This does the same thing on the regional map.

twointegers gettotalinflow(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    twointegers total;
    twointegers dest;

    total.x = 0;
    total.y = 0;

    for (int i = x - 1; i <= x + 1; i++)
    {
        for (int j = y - 1; j <= y + 1; j++)
        {
            if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight)
            {
                dest = getregionalflowdestination(region, i, j, 0);

                if (dest.x == x && dest.y == y)
                {
                    total.x = total.x + region.riverjan(i, j);
                    total.y = total.y + region.riverjul(i, j);
                }
            }
        }
    }
    return total;
}

// This function finds the lowest neighbouring point that's still higher, on the regional map.

twointegers findlowesthigher(region& region, int dx, int dy, int x, int y, int janload, int julload, int crount, vector<vector<bool>>& mainriver)
{
    twointegers nextpoint;
    int origelev = region.map(x, y);

    nextpoint.x = -1;
    nextpoint.y = -1;

    int elev = 100000000;

    for (int i = x - 1; i <= x + 1; i++)
    {
        for (int j = y - 1; j <= y + 1; j++)
        {
            if (region.map(i, j) < elev && region.map(i, j) >= origelev && region.sea(i, j) == 0 && region.riverjan(i, j) == 0 && region.riverjul(i, j) == 0)
            {
                bool keepgoing = 1;

                if (crount > 1)
                {
                    int neighbouringrivers = 0;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (region.riverjan(k, l) != 0 && region.riverjan(k, l) != janload && region.riverjul(k, l) != 0 && region.riverjul(k, l) != julload)
                                keepgoing = 0; // After the first point, avoid going next to the main river in this tile.

                            if (region.riverjan(k, l) != 0 || region.riverjul(k, l) != 0)
                                neighbouringrivers++;

                            if (region.sea(k, l) == 1) // Avoid going near coastlines after the first point, too.
                                keepgoing = 0;
                        }
                    }

                    if (neighbouringrivers > 2)
                        keepgoing = 0;
                }
                else
                {
                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (region.riverdir(k, l) != 0) // && mainriver[dx-k][dy-l]==0)
                                keepgoing = 0;
                        }
                    }
                }

                // Avoid crossing existing rivers diagonally.

                if (i == x + 1 && j == y + 1)
                {
                    if (region.riverdir(x + 1, y) == 6 || region.riverdir(x, y + 1) == 2)
                        keepgoing = 0;
                }

                if (i == x - 1 && j == y - 1)
                {
                    if (region.riverdir(x, y - 1) == 6 || region.riverdir(x - 1, y) == 2)
                        keepgoing = 0;
                }

                if (i == x + 1 && j == y - 1)
                {
                    if (region.riverdir(x, y - 1) == 4 || region.riverdir(x + 1, y) == 8)
                        keepgoing = 0;
                }

                if (i == x - 1 && j == y + 1)
                {
                    if (region.riverdir(x - 1, y) == 4 || region.riverdir(x, y + 1) == 8)
                        keepgoing = 0;
                }

                if (keepgoing == 1)
                {
                    elev = region.map(i, j);
                    nextpoint.x = i;
                    nextpoint.y = j;
                }
            }
        }
    }
    return nextpoint;
}

// This function tells whether a regional cell is next to a lake.

int nexttolake(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    int dist = 1;

    for (int i = x - dist; i <= x + dist; i++)
    {
        int ii = i;

        for (int j = y - dist; j <= y + dist; j++)
        {
            if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight && region.lakesurface(ii, j) != 0)
                return (region.lakesurface(ii, j));
        }
    }
    return (0);
}

// This function finds the level of the nearest lake to the given point in the regional map.

int getnearestlakelevel(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    for (int n = 1; n <= 200; n++)
    {
        for (int i = x - n; i <= x + n; i++)
        {
            for (int j = y - n; j <= y + n; j++)
            {
                if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight)
                {
                    if (region.lakesurface(i, j) != 0)
                        return region.lakesurface(i, j);
                }
            }
        }
    }
    return 0;
}

// This function finds the special value of the nearest lake to the given point in the regional map.

int getnearestlakespecial(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    for (int n = 1; n <= 200; n++)
    {
        for (int i = x - n; i <= x + n; i++)
        {
            for (int j = y - n; j <= y + n; j++)
            {
                if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight)
                {
                    if (region.lakesurface(i, j) != 0)
                        return region.special(i, j);
                }
            }
        }
    }
    return 0;
}

// This function finds the nearest river on the regional map.

twointegers findclosestriver(region& region, int x, int y, bool delta)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    twointegers nearest;

    nearest.x = -1;;
    nearest.y = -1;

    for (int n = 1; n <= 14; n++)
    {
        int ii = 0;

        for (int i = x - 1; i <= x + n; i++)
        {
            int jj = 0;

            for (int j = y - n; j <= y + n; j++)
            {
                if (ii * jj < (n * n) + n)
                {
                    if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight && region.waterdir(i, j, delta) != 0)
                    {
                        nearest.x = i;
                        nearest.y = j;

                        return nearest;
                    }
                }

                jj++;
            }

            ii++;
        }
    }

    return nearest;
}

// As above, but more rough and ready.

twointegers findclosestriverquickly(region& region, int x, int y)
{
    twointegers nearest;

    if (region.riverdir(x, y) != 0)
    {
        nearest.x = x;
        nearest.y = y;
        return nearest;
    }

    nearest.x = -1;;
    nearest.y = -1;

    for (int n = 1; n <= 30; n++)
    {
        for (int i = x - n; i <= x + n; i++)
        {
            int j = y - n;

            if (region.riverdir(i, j) != 0)
            {
                nearest.x = i;
                nearest.y = j;
                return nearest;
            }

            j = y + n;

            if (region.riverdir(i, j) != 0)
            {
                nearest.x = i;
                nearest.y = j;
                return nearest;
            }
        }

        for (int j = y - n; j <= y + n; j++)
        {
            int i = x - n;

            if (region.riverdir(i, j) != 0)
            {
                nearest.x = i;
                nearest.y = j;
                return nearest;
            }

            i = x + n;

            if (region.riverdir(i, j) != 0)
            {
                nearest.x = i;
                nearest.y = j;
                return nearest;
            }
        }
    }

    return nearest;
}

// This function finds how many inflows there are to a cell on the regional map.

int countinflows(region& region, int x, int y)
{
    int rwidth = region.rwidth();
    int rheight = region.rheight();

    int inflows = 0;
    twointegers dest;

    for (int i = x - 1; i <= x + 1; i++)
    {
        for (int j = y - 1; j <= y + 1; j++)
        {
            if (i >= 0 && i <= rwidth && j >= 0 && j <= rheight)
            {
                dest = getregionalflowdestination(region, i, j, 0);

                if (dest.x == x && dest.y == y)
                    inflows++;
            }
        }
    }

    return inflows;
}
