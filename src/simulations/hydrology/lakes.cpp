#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "lakes.hpp"
#include "drainage.hpp"
#include "lake_effects.hpp"

using namespace std;

namespace
{
void createlakedepression(planet& world, int centrex, int centrey, int origlevel, int steepness, vector<vector<int>>& basins, int limit, bool up, vector<vector<int>>& avoid);
void drawlake(planet& world, int shapenumber, int centrex, int centrey, vector<vector<int>>& thislake, int lakeno, vector<vector<int>>& checked, vector<vector<int>>& nolake, int templatesize, int minrain, int maxtemp, boolshapetemplate smalllake[], boolshapetemplate largelake[]);
int lakeoutline(planet& world, vector<vector<int>>& thislake, int lakeno, int x, int y);
void sortlakerivers(planet& world, int leftx, int lefty, int rightx, int righty, int centrex, int centrey, int outflowx, int outflowy, vector<vector<int>>& thislake, int lakeno);
void cleanlakes(planet& world);
twointegers createriftlake(planet& world, int startx, int starty, int lakelength, vector<vector<int>>& nolake);

// This function lowers land around a given point to put it at the centre of a large depression.

void createlakedepression(planet& world, int centrex, int centrey, int origlevel, int steepness, vector<vector<int>>& basins, int limit, bool up, vector<vector<int>>& avoid)
{
    int width = world.width();
    int height = world.height();

    if (limit == 0)
        limit = 1000;

    int maxmountains = 200; // Mountains larger than this will be left alone.

    for (int radius = 1; radius <= limit; radius++)
    {
        int numberchanged = 0;

        int currentlevel = origlevel + (steepness * (radius - 1)); // Cells checked on this round must be no higher than this.

        for (int i = centrex - radius; i <= centrex + radius; i++)
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = centrey - radius; j <= centrey + radius; j++)
            {
                if (j >= 0 && j <= height)
                {
                    if (up == 1 || world.nom(ii, j) > currentlevel) // If this cell is higher than the current limit (or if we're raising cells as well as lowering them)
                    {
                        if (basins[ii][j] == 0 && avoid[ii][j] == 0 && world.sea(ii, j) == 0 && world.truelake(ii, j) == 0 && world.mountainheight(ii, j) < maxmountains)
                        {
                            if ((ii - centrex) * (ii - centrex) + (j - centrey) * (j - centrey) <= radius * radius) // If this cell is within the current radius
                            {
                                int thislevel = currentlevel; // The level of this particular cell, which will be a variation on the current level of the whole radius.

                                if (radius > 1 && up == 0) // Add a bit of variation to the slope.
                                    thislevel = thislevel - random(0, steepness);

                                world.setnom(ii, j, thislevel);
                                numberchanged++;
                                basins[ii][j] = 1;
                            }
                        }
                    }
                }
            }
        }

        if (numberchanged == 0)
            return;
    }
}

// Puts a lake template onto the lake map.

void drawlake(planet& world, int shapenumber, int centrex, int centrey, vector<vector<int>>& thislake, int lakeno, vector<vector<int>>& checked, vector<vector<int>>& nolake, int templatesize, int minrain, int maxtemp, boolshapetemplate smalllake[], boolshapetemplate largelake[])
{
    int width = world.width();
    int height = world.height();
    int riverlandreduce = world.riverlandreduce();

    if (centrey == 0 || centrey == height)
    {
        world.setlakesurface(centrex, centrey, 0);
        return;
    }

    twointegers nearestsea, flow;

    int minlakedistance = 10; // Minimum distance between lakes.
    int minmountaindistance = 2; // Minimum distance from mountains.
    int maxmountainlakeheight = 1000; // Ignore mountains smaller than this.
    float mindepth = 5; // Minimum depth.

    int imheight, imwidth;

    if (templatesize == 1)
    {
        imheight = smalllake[shapenumber].ysize() - 1;
        imwidth = smalllake[shapenumber].xsize() - 1;
    }
    else
    {
        imheight = largelake[shapenumber].ysize() - 1;
        imwidth = largelake[shapenumber].xsize() - 1;
    }

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

    int leftx = centrex;
    int lefty = centrey;
    int rightx = centrex;
    int righty = centrey; // These are the coordinates of the furthermost pixels of the lake.
    bool wrapped = 0; // If this is 1, the lake wraps over the edge of the map.

    int leftr = random(1, 2); // If it's 1 then we reverse it left-right
    int downr = random(1, 2); // If it's 1 then we reverse it top-bottom

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

    // First, check to see if another lake is already too close.

    bool tooclose = 0;

    int mapi = -1;
    int mapj = -1; // These are the offsets on the actual map, as opposed to on the lake template.

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            bool point;

            if (templatesize == 1)
                point = smalllake[shapenumber].point(i, j);
            else
                point = largelake[shapenumber].point(i, j);

            if (point == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                    {
                        xx = wrap(xx, width);
                        wrapped = 1;
                    }

                    if (world.sea(xx, yy) == 0)
                    {
                        if (nolake[xx][yy] == 1)
                            tooclose = 1;

                        if (world.volcano(xx, yy) != 0)
                            tooclose = 1;

                        if (tooclose == 0) // Only do this next check if we haven't already failed
                        {
                            for (int xxx = xx - minlakedistance; xxx <= xx + minlakedistance; xxx++)
                            {
                                int xxxx = xxx;

                                if (xxxx<0 || xxxx>width)
                                    xxxx = wrap(xxxx, width);

                                for (int yyy = yy - minlakedistance; yyy <= yy + minlakedistance; yyy++)
                                {
                                    if (yyy >= 0 && yyy <= height)
                                    {
                                        if (thislake[xxxx][yyy] != 0 || world.volcano(xxxx, yyy) != 0) // There's another lake here already! Or a volcano.
                                        {
                                            tooclose = 1;
                                            xxx = xx + minlakedistance;
                                            yyy = yy + minlakedistance;
                                        }
                                    }
                                }
                            }
                        }

                        if (tooclose == 0) // Only do this next check if we haven't already failed the last!
                        {
                            for (int xxx = xx - minmountaindistance; xxx <= xx + minmountaindistance; xxx++)
                            {
                                int xxxx = xxx;

                                if (xxxx<0 || xxxx>width)
                                    xxxx = wrap(xxxx, width);

                                for (int yyy = yy - minmountaindistance; yyy <= yy + minmountaindistance; yyy++)
                                {
                                    if (yyy >= 0 && yyy <= height)
                                    {
                                        if (world.mountainheight(xxxx, yyy) > maxmountainlakeheight) // There are mountains here.
                                        {
                                            tooclose = 1;
                                            xxx = xx + minmountaindistance;
                                            yyy = yy + minmountaindistance;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                        tooclose = 1;
                }
            }
        }
    }

    if (tooclose == 1)
        return;

    if (wrapped == 1)
    {
        leftx = 0;
        lefty = 0;
        rightx = width;
        righty = height;
    }

    // Now put down the actual lake.

    mapi = -1;
    mapj = -1;

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            bool point;

            if (templatesize == 1)
                point = smalllake[shapenumber].point(i, j);
            else
                point = largelake[shapenumber].point(i, j);

            if (point == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                    {
                        xx = wrap(xx, width);
                        wrapped = 1;
                    }

                    thislake[xx][yy] = lakeno; // Mark it on the keeping-track array.

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

    // Now we find out what the lowest point on the edge of the lake is. This will be the point from which the outflowing river emerges, and it will determine the surface level of the whole lake.

    int smallest = world.maxelevation();
    int outflowx = -1;
    int outflowy = -1;

    for (int i = leftx; i <= rightx; i++)
    {
        for (int j = lefty; j <= righty; j++)
        {
            if (thislake[i][j] == lakeno)
            {
                float border = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (thislake[kk][l] != lakeno)
                            {
                                border = 1;
                                k = i + 1;
                                l = j + 1;
                            }
                        }
                    }
                }

                if (border == 1)
                {
                    if (world.nom(i, j) < smallest)
                    {
                        smallest = world.nom(i, j);
                        outflowx = i;
                        outflowy = j;
                    }
                }
            }
        }
    }

    int surfaceheight = smallest - riverlandreduce;

    // Now make the lake itself that height.

    for (int i = leftx; i <= rightx; i++)
    {
        for (int j = lefty; j <= righty; j++)
        {
            if (thislake[i][j] == lakeno)
                world.setlakesurface(i, j, surfaceheight);
        }
    }

    // Now we must alter the terrain height underneath the lake.

    int maxdepth = imwidth / 8;

    if (maxdepth < 10)
        maxdepth = 10;

    if (maxdepth > 40)
        maxdepth = 40;

    float depthmult = (float)random((int)mindepth, maxdepth); // The higher this is, the deeper the lake will be.

    float maxmaxdepth = (float)random(500, 1000); // Absolute maximum depth for lakes in this world.

    for (int i = leftx; i <= rightx; i++)
    {
        for (int j = lefty; j <= righty; j++)
        {
            if (thislake[i][j] == lakeno)
            {
                float nom = (float)world.nom(i, j);

                float depth = nom - (float)surfaceheight;

                if (depth<0.0f) // If the elevation here is higher than the lake surface
                    depth = (0.0f - depth) * depthmult;

                if (depth > maxmaxdepth)
                    depth = maxmaxdepth;

                if (depth < mindepth)
                    depth = mindepth;

                int newnom = surfaceheight - (int)depth;

                if (newnom < 1)
                    newnom = 1;

                world.setnom(i, j, newnom);
            }
        }
    }

    // Now do the rivers.

    sortlakerivers(world, leftx, lefty, rightx, righty, centrex, centrey, outflowx, outflowy, thislake, lakeno);

    // Now mark a "start point" somewhere on the edge of the lake, which we will use when it comes to the regional map.

    makelakestartpoint(world, thislake, lakeno, leftx, lefty, rightx, righty);
}

// This checks whether a given tile is on the edge of a given lake.

int lakeoutline(planet& world, vector<vector<int>>& thislake, int lakeno, int x, int y)
{
    int width = world.width();
    int height = world.height();

    int xx = x - 1;
    int yy = y;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x + 1;
    yy = y;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x;
    yy = y - 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x;
    yy = y + 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x - 1;
    yy = y - 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x + 1;
    yy = y + 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x + 1;
    yy = y - 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    xx = x - 1;
    yy = y + 1;

    if (xx<0 || xx>width)
        xx = wrap(xx, width);

    if (yy >= 0 && yy <= height)
    {
        if (thislake[xx][yy] != lakeno)
            return(1);
    }

    return (0);
}

// Deals with the rivers going into and out of a lake.

void sortlakerivers(planet& world, int leftx, int lefty, int rightx, int righty, int centrex, int centrey, int outflowx, int outflowy, vector<vector<int>>& thislake, int lakeno)
{
    int width = world.width();
    int height = world.height();

    // Mark the route of the outflowing river (if any) on this array. We will avoid messing with any of it.

    vector<vector<int>> outflow(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    markriver(world, outflowx, outflowy, outflow, 0);

    // Now go over the lake area. Every river tile that is under or next to the lake should be pointing to the centre of the lake.

    vector<vector<int>> checkedriverlaketiles(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> removedrivers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int riverno = 0;

    for (int i = leftx; i <= rightx; i++)
    {
        for (int j = lefty; j <= righty; j++)
        {
            for (int k = i - 1; k <= i + 1; k++)
            {
                int kk = k;

                if (kk<0 || k>width)
                    kk = wrap(kk, width);

                for (int l = j - 1; l <= j + 1; l++)
                {
                    if (l >= 0 && l <= height)
                    {
                        // kk and l are the points around the current one.

                        if (outflow[kk][l] == 1)
                            checkedriverlaketiles[kk][l] = 1;

                        if (checkedriverlaketiles[kk][l] == 0)
                        {
                            bool mustdivert = 0;

                            if (world.lakesurface(kk, l) != 0)
                                mustdivert = 1;
                            else
                            {
                                if (nexttolake(world, kk, l) != 0) // If it's next to our lake
                                    mustdivert = 1;
                            }

                            if (mustdivert == 1)
                            {
                                riverno++;

                                riverno = divertlakeriver(world, kk, l, centrex, centrey, removedrivers, riverno, outflowx, outflowy, outflow);
                            }

                            checkedriverlaketiles[kk][l] = 1;
                        }
                    }
                }
            }
        }
    }

    // At this point, the outflowing river is reduced to nothing or almost nothing.
    // Now we need to work out how big it should actually be.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            checkedriverlaketiles[i][j] = 0;
    }

    leftx--;
    lefty--;
    rightx++;
    righty++;

    if (leftx<0 || rightx>width)
    {
        leftx = 0;
        rightx = width;
    }

    if (lefty < 0)
        lefty = 0;

    if (righty > height)
        righty = height;

    int janload = 0;
    int julload = 0;

    int origjanoutflow = 0;
    int origjuloutflow = 0;

    for (int i = leftx; i <= rightx; i++)
    {
        for (int j = lefty; j <= righty; j++)
        {
            if (thislake[i][j] != lakeno)
            {
                float border = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (thislake[kk][l] == lakeno)
                            {
                                border = 1;
                                k = i + 1;
                                l = j + 1;
                            }
                        }
                    }
                }

                if (border == 1)
                {
                    janload = janload + world.riverjan(i, j);
                    julload = julload + world.riverjul(i, j);
                }
            }
        }
    }

    if (janload < 0)
        janload = 0;

    if (julload < 0)
        julload = 0;

    janload = janload / 2;
    julload = julload / 2; // Because as it stands, it contains the original outflows as well as the inflows.

    if (janload == 0 && julload == 0) // If somehow we've ended up with zero flow for the outflowing river
    {
        janload = origjanoutflow;
        julload = origjuloutflow;
    }

    // Now we apply the new load to the outflowing river.

    vector<vector<bool>> loadadded(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int x = outflowx;
    int y = outflowy;
    bool keepgoing = 1;

    do
    {
        if (loadadded[x][y] == 0) // Just to make sure it's not somehow going round in circles.
        {
            world.setriverjan(x, y, world.riverjan(x, y) + janload);
            world.setriverjul(x, y, world.riverjul(x, y) + julload);

            loadadded[x][y] = 1;

            int dir = world.riverdir(x, y);

            if (dir == 8 || dir == 1 || dir == 2)
                y--;

            if (dir == 4 || dir == 5 || dir == 6)
                y++;

            if (dir == 2 || dir == 3 || dir == 4)
                x++;

            if (dir == 6 || dir == 7 || dir == 8)
                x--;

            if (x < 0)
                x = width;

            if (x > width)
                x = 0;

            if (thislake[x][y] == lakeno && (x != outflowx || y != outflowy)) // If it's somehow flowed back into the lake!
                keepgoing = 0;

            if (y<0 || y>height)
                keepgoing = 0;

            if (world.sea(x, y) == 1)
                keepgoing = 0;

            if (world.truelake(x, y) != 0 && thislake[x][y] != lakeno) // If it's gone into another lake
                keepgoing = 0;
        }
        else
            keepgoing = 0;

    } while (keepgoing == 1);
}

// This function cleans up the lake map, removing any absurd values.

void cleanlakes(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.lakesurface(i, j) > maxelev)
                    world.setlakesurface(i, j, 0);

                if (world.truelake(i, j) == 1)
                    world.setextraelev(i, j, 0);
            }
        }
    });
}

// This function puts a particular rift lake onto the rift lake map.

twointegers createriftlake(planet& world, int startx, int starty, int lakelength, vector<vector<int>>& nolake)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int riverlandreduce = world.riverlandreduce();

    int maxdepth = 1800; // Maximum depth
    int mindepth = 100; // Minimum depth

    vector<vector<int>> removedrivers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> currentriver(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    // First we need to work out the surface level of this lake. It will be at the level of the river flowing out of it.

    int x = startx;
    int y = starty;

    for (int n = 1; n <= lakelength; n++)
    {
        int dir = world.riverdir(x, y);

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

        if (y == 0 || y == height)
            n = lakelength;
    }

    int surfacelevel = world.nom(x, y) - riverlandreduce;

    world.setlakestart(startx, starty, 1); // This marks this as the beginning of a rift lake, which we'll need to know when it comes to the regional map.

    markriver(world, startx, starty, currentriver, 0);

    twointegers nextpoint;

    int depth = random(mindepth, maxdepth);

    x = startx;
    y = starty;
    int done = 0;

    for (int n = 1; n <= lakelength; n++)
    {
        bool keepgoing = 1;

        if (nolake[x][y] != 0)
        {
            keepgoing = 0;
            n = lakelength;
        }

        for (int i = x - 2; i <= x + 2; i++)
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = y - 2; j <= y + 2; j++)
            {
                if (j >= 0 && j <= height && world.lakesurface(ii, j) != 0)
                {
                    keepgoing = 0;
                    i = x + 2;
                    j = y + 2;
                    n = lakelength;
                }
                if (j > 0 && j <= height)
                {
                    if (world.riftlakesurface(ii, j) != 0 && world.riftlakesurface(ii, j) != surfacelevel)
                    {
                        keepgoing = 0;
                        i = x + 2;
                        j = y + 2;
                        n = lakelength;
                    }
                }
            }
        }

        if (keepgoing == 0)
        {
            if (done < 3) // If we weren't able to make this lake long enough, get rid of it again.
            {
                x = startx;
                y = starty;

                for (int n = 1; n <= done; n++)
                {
                    world.setriftlakesurface(x, y, 0);
                    world.setriftlakebed(x, y, 0);

                    int dir = world.riverdir(x, y);

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

                    if (y == 0 || y == height)
                        n = done;
                }
            }

            nextpoint.x = -1;
            nextpoint.y = -1;
            return (nextpoint);
        }

        world.setriftlakesurface(x, y, surfacelevel);

        done++;

        int thislevel = surfacelevel - depth + randomsign(random(0, 200));

        if (thislevel > surfacelevel - 50)
            thislevel = surfacelevel - 50;

        if (thislevel < 1)
            thislevel = 1;

        world.setriftlakebed(x, y, thislevel);

        int dir = world.riverdir(x, y);

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

        if (y == 0 || y == height)
            n = lakelength;
    }

    if (done < 3) // If we weren't able to make this lake long enough, get rid of it again.
    {
        x = startx;
        y = starty;

        for (int n = 1; n <= done; n++)
        {
            world.setriftlakesurface(x, y, 0);
            world.setriftlakebed(x, y, 0);

            int dir = world.riverdir(x, y);

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

            if (y == 0 || y == height)
                n = done;
        }
        nextpoint.x = x;
        nextpoint.y = y;

        return nextpoint;
    }

    // Now we need to make sure that all rivers bordering this rift lake flow into it.

    x = startx;
    y = starty;

    int riverno = 1;

    twointegers dest;

    for (int n = 1; n <= lakelength; n++)
    {
        for (int i = x - 1; i <= x + 1; i++)
        {
            int ii = i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = y - 1; j <= y + 1; j++)
            {
                if (j >= 0 && j <= height)
                {
                    if (world.riftlakesurface(ii, j) == 0) // If this is not in the rift lake
                    {
                        if (currentriver[i][j] == 0)
                        {
                            int lowest = maxelev;

                            int newdestx = -1;
                            int newdesty = -1;

                            for (int k = i - 1; k <= i + 1; k++) // Find the deepest neighbouring rift lake tile
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - 1; l <= j + 1; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.riftlakesurface(kk, l) == surfacelevel)
                                        {
                                            if (world.riftlakesurface(kk, l) < lowest)
                                            {
                                                lowest = world.riftlakesurface(kk, l);

                                                newdestx = kk;
                                                newdesty = l;
                                            }
                                        }
                                    }
                                }
                            }

                            if (newdestx != -1) // Now divert the flow towards that tile
                            {
                                divertriver(world, ii, j, newdestx, newdesty, removedrivers, riverno);

                                riverno++;
                            }
                        }
                    }
                }
            }
        }

        int dir = world.riverdir(x, y);

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

        if (y == 0 || y == height)
            n = lakelength;
    }

    nextpoint.x = x;
    nextpoint.y = y;

    return (nextpoint);
}
}

// This function creates areas of sea that will later become salt lakes.

void createsaltlakes(planet& world, int& lakesplaced, vector<vector<vector<int>>>& saltlakemap, vector<vector<int>>& nolake, vector<vector<int>>& basins, boolshapetemplate smalllake[])
{
    fast_srand(static_cast<long>(deterministiccontextseed(world.seed(), 0x3001) & 0x7fffffffull));

    int width = world.width();
    int height = world.height();

    int maxsaltrain = tuning::climate::saltlakes::maxRain;
    int minsalttemp = tuning::climate::saltlakes::minimumTemperature;
    int minflow = tuning::climate::saltlakes::minimumRiverFlow;

    int saltlakechance = tuning::climate::saltlakes::lakeChance;

    bool dodepressions = 0;

    if (random(1, tuning::climate::saltlakes::depressionChance) == 1)
        dodepressions = 1;

    vector<vector<int>> avoid(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This is for cells *not* to alter. This will remain empty (it's just here so we can use the depression-creating routine for normal lakes as well, where we don't want to mess with the path of the outflowing river.)

    for (int i = tuning::climate::saltlakes::edgeMargin; i <= width - tuning::climate::saltlakes::edgeMargin; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (basins[i][j] == 0 && nolake[i][j] == 0 && world.sea(i, j) == 0 && world.riverjan(i, j) + world.riverjul(i, j) >= minflow)
            {
                if (random(1, saltlakechance) == 1)
                {
                    int averain = (world.winterrain(i, j) + world.summerrain(i, j)) / 2;
                    int avetemp = (world.maxtemp(i, j) + world.mintemp(i, j)) / 2;

                    if (averain<maxsaltrain && avetemp>minsalttemp)
                    {
                        lakesplaced++;
                        placesaltlake(world, i, j, 0, dodepressions, saltlakemap, basins, avoid, nolake, smalllake);
                    }
                }
            }
        }
    }
}

// This function puts a patch of sea on the map that will later be turned into a salt lake.

void placesaltlake(planet& world, int centrex, int centrey, bool large, bool dodepressions, vector<vector<vector<int>>>& saltlakemap, vector<vector<int>>& basins, vector<vector<int>>& avoid, vector<vector<int>>& nolake, boolshapetemplate smalllake[])
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int riverlandreduce = world.riverlandreduce();

    int shapenumber;

    if (large || random(1, tuning::climate::saltlakes::largeShapeChance) == 1)
        shapenumber = random(tuning::climate::saltlakes::largeShapeMin, tuning::climate::saltlakes::largeShapeMax);
    else
        shapenumber = random(tuning::climate::saltlakes::smallShapeMin, tuning::climate::saltlakes::smallShapeMax);

    int depth = random(tuning::climate::saltlakes::depthMin, tuning::climate::saltlakes::depthMax);

    // The surfaceheight will be a lot lower than this point originally is.

    int origheight = world.nom(centrex, centrey);
    int surfaceheight = (sealevel + (origheight - sealevel) / 2) - riverlandreduce;

    if (surfaceheight <= sealevel)
        surfaceheight = sealevel + 1;

    int imheight = smalllake[shapenumber].ysize() - 1;
    int imwidth = smalllake[shapenumber].xsize() - 1;

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

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

    int leftx = centrex;
    int lefty = centrey;
    int rightx = centrex;
    int righty = centrey; // These are the coordinates of the furthermost pixels of the lake.
    bool wrapped = 0; // If this is 1, the lake wraps over the edge of the map.

    bool tooclose = 0;

    int mapi = -1;
    int mapj = -1;

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            if (smalllake[shapenumber].point(i, j) == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    if (nolake[xx][yy] == 1 || world.nom(xx, yy) <= sealevel) // Don't try to put lakes on top of sea...
                        tooclose = 1;
                    else
                    {
                        if (surfaceheight > world.nom(xx, yy) - riverlandreduce) // Lakes can't extend onto areas where local rivers will be lower than the surface of the lake.
                            tooclose = 1;
                        else
                        {

                            for (int k = xx - 2; k <= xx + 2; k++) // Check nearby cells for other lakes.
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = yy - 2; l <= yy + 2; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (saltlakemap[kk][l][0] != 0 && saltlakemap[kk][l][0] != surfaceheight)
                                            tooclose = 1;

                                        if (world.nom(kk, l) <= sealevel && saltlakemap[kk][l][0] != surfaceheight)
                                            tooclose = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (tooclose == 1)
        return;

    // Now actually do the lake.

    mapi = -1;
    mapj = -1;

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            if (smalllake[shapenumber].point(i, j) == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                    {
                        wrapped = 1;
                        xx = wrap(xx, width);
                    }

                    world.setnom(xx, yy, sealevel - 1);
                    saltlakemap[xx][yy][0] = surfaceheight;
                    saltlakemap[xx][yy][1] = surfaceheight - (depth + randomsign(random(0, 4)));

                    if (saltlakemap[xx][yy][1] < 1)
                        saltlakemap[xx][yy][1] = 1;

                    if (xx < leftx)
                        leftx = xx;

                    if (xx > rightx)
                        rightx = xx;

                    if (yy < lefty)
                        lefty = yy;
                    if (yy > righty)
                        righty = yy;

                    if (world.mountainheight(xx, yy) != 0) // Remove any mountains that might be here.
                    {
                        for (int dir = 1; dir <= 8; dir++)
                            deleteridge(world, xx, yy, dir);
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

    // Now we need to lower the land around this sea/lake, to ensure that all rivers for a wide radius will flow into it.

    int steepness = 20; // Steepness of the basin around the lake.
    int limit = 0; // size*3 ` Maximum extent of the basin.

    if (dodepressions)
        createlakedepression(world, centrex, centrey, surfaceheight, steepness, basins, limit, 0, avoid);

    // Now we need to create a "start point" on the sea/lake.

    twointegers startpoint;
    startpoint.x = -1;
    startpoint.y = -1;

    int crount = 0;

    do
    {
        for (int i = leftx; i <= rightx; i++)
        {
            for (int j = lefty; j <= righty; j++)
            {
                if (saltlakemap[i][j][0] == surfaceheight)
                {
                    if (vaguelycoastal(world, i, j) == 1) // If this is on the edge of the sea/lake
                    {
                        if (random(1, 3) == 1)
                        {
                            startpoint.x = i;
                            startpoint.y = j;
                        }
                    }
                }
            }
        }

        crount++;

        if (crount > 1000) // Something's gone horribly wrong
            return;

    } while (startpoint.x == -1);

    world.setlakestart(startpoint.x, startpoint.y, 1);
}

// This function turns the inland seas into proper salt lakes.

void convertsaltlakes(planet& world, vector<vector<vector<int>>>& saltlakemap)
{
    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (saltlakemap[i][j][0] != 0)
            {
                world.setlakesurface(i, j, saltlakemap[i][j][0]);
                world.setnom(i, j, saltlakemap[i][j][1]);

                world.setspecial(i, j, 100);
            }
        }
    }
}

// This function creates the major freshwater lakes.

void createlakemap(planet& world, vector<vector<int>>& nolake, boolshapetemplate smalllake[], boolshapetemplate largelake[])
{
    fast_srand(static_cast<long>(deterministiccontextseed(world.seed(), 0x3002) & 0x7fffffffull));

    int width = world.width();
    int height = world.height();
    int glacialtemp = world.glacialtemp();

    int minlake = 1000; // Flows of this size or higher might have lakes on them.
    int maxlake = 12000; // Flows larger than this can't have lakes.
    int lakechance = random(8, 16); //2; //10; //random(80,150); // The higher this is, the less likely lakes are.
    int minrain = 200;
    int maxtemp = 25; // Lakes won't appear where there is less rain *and* higher temperature than this.

    vector<vector<int>> thislake(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> checked(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            thislake[i][j] = 0; // Clear the thislake array, where we will keep track of the lakes as we do them.
        }
    }

    // This array is used to make the lakes cluster in particular parts of the world.

    vector<vector<bool>> lakeprobability(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            lakeprobability[i][j] = 0;
    }

    int clusternumber = 50; // Number of possible clusters of lakes on the map.
    int minrad = 10;
    int maxrad = 40; // Possible sizes of the clusters

    for (int n = 0; n < clusternumber; n++) // Put some circles on this array, to allow for lake clusters.
    {
        int centrex = random(0, width);
        int centrey = random(0, height);

        int radius = random(minrad, maxrad);

        for (int i = -radius; i <= radius; i++)
        {
            int ii = centrex + i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = -radius; j <= radius; j++)
            {
                int jj = centrey + j;

                if (jj >= 0 && jj <= height)
                {
                    if (i * i + j * j < radius * radius + radius)
                        lakeprobability[ii][jj] = 1;
                }
            }
        }
    }

    // Now make the lakes

    twointegers nearestsea;
    int lakeno = 0;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int flow = world.riveraveflow(i, j);

            if (world.averainfloat(i, j) >= static_cast<float>(minrain) ||
                world.avetemp(i, j) <= maxtemp)
            {
                if (flow >= minlake && flow <= maxlake && lakeprobability[i][j] == 1 && random(1, lakechance) == 1)
                {
                    if (nolake[i][j] == 0)
                    {
                        lakeno++;

                        int templatesize = 1;
                        int shapenumber;

                        if (world.avetemp(i, j) <= glacialtemp)
                            shapenumber = random(0, 3);
                        else
                        {
                            shapenumber = random(2, 11);

                            if (random(1, 6) != 1)
                                shapenumber = random(2, 11);
                            else
                            {
                                shapenumber = random(0, 9);
                                templatesize = 2;
                            }
                        }
                        drawlake(world, shapenumber, i, j, thislake, lakeno, checked, nolake, templatesize, minrain, maxtemp, smalllake, largelake);
                    }
                }
            }
        }
    }

    // Now clean up a little.

    cleanlakes(world);

    // Now we need to add precipitation from the lakes.

    vector<vector<int>> lakewinterrainmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> lakesummerrainmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    lakerain(world, lakewinterrainmap, lakesummerrainmap);

    // Now we add new rivers from that extra precipitation.

    lakerivers(world, lakewinterrainmap, lakesummerrainmap);

    // Now we check that no rivers shrink inexplicably.

    checkglobalflows(world);
}

// This creates the rift lake map.

void createriftlakemap(planet& world, vector<vector<int>>& nolake)
{
    fast_srand(static_cast<long>(deterministiccontextseed(world.seed(), 0x3003) & 0x7fffffffull));

    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int riftlakechance = random(500, 4000); // The higher this is, the fewer rift lakes there will be.
    int minlake = 500; // Flows of this size or higher might have lakes on them.
    int maxlake = 2000; // Flows larger than this can't have lakes.
    int minlakelength = 4;
    int maxlakelength = 10;

    twointegers nextpoint;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) >= sealevel && world.lakesurface(i, j) == 0 && nolake[i][j] == 0 && world.riftlakesurface(i, j) == 0)
            {
                if (world.riveraveflow(i, j) > minlake && world.riveraveflow(i, j) < maxlake && random(1, riftlakechance) == 1)
                {
                    // Check that no rift lakes currently exist downstream of here.

                    bool keepgoing = 1;
                    bool found = 0;

                    int x = i;
                    int y = j;

                    int crount = 0;

                    do
                    {
                        crount++;

                        int dir = world.riverdir(x, y);

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

                        if (y<0 || y>height)
                            keepgoing = 0;

                        if (dir == 0 || world.sea(x, y) == 1)
                            keepgoing = 0;

                        if (world.riftlakesurface(x, y) != 0 || world.lakesurface(x, y) != 0)
                        {
                            found = 1;
                            keepgoing = 0;
                        }

                        if (crount > 100)
                            keepgoing = 0;

                    } while (keepgoing == 1);

                    if (found == 0)
                    {
                        x = i;
                        y = j;
                        keepgoing = 1;

                        do
                        {
                            int lakelength = random(minlakelength, maxlakelength);

                            nextpoint = createriftlake(world, x, y, lakelength, nolake); // Make a lake.

                            x = nextpoint.x;
                            y = nextpoint.y;

                            if (x == -1)
                                keepgoing = 0;
                            else
                            {
                                int gaplength = random(2, 10);

                                for (int n = 1; n <= gaplength; n++) // Move down the river a bit.
                                {
                                    if (x<0 || x>width)
                                        x = wrap(x, width);

                                    if (y >= 0 && y <= height)
                                    {
                                        int dir = world.riverdir(x, y);

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

                                        if (y == 0 || y == height)
                                        {
                                            keepgoing = 0;
                                            n = gaplength;
                                        }
                                        else
                                        {
                                            if (y >= 0 && y <= height)
                                            {
                                                if (world.lakesurface(x, y) != 0)
                                                {
                                                    keepgoing = 0;
                                                    n = gaplength;
                                                }

                                                if (world.nom(x, y) <= sealevel)
                                                {
                                                    keepgoing = 0;
                                                    n = gaplength;
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        keepgoing = 0;
                                        n = gaplength;
                                    }
                                }

                                if (random(1, 4) == 1)
                                    keepgoing = 0;
                            }

                        } while (keepgoing == 1);
                    }
                }
            }
        }
    }

    // Now we need to add precipitation from the lakes.

    vector<vector<int>> lakewinterrainmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> lakesummerrainmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    riftlakerain(world, lakewinterrainmap, lakesummerrainmap);

    // Now we add new rivers from that extra precipitation.

    lakerivers(world, lakewinterrainmap, lakesummerrainmap); // We can use the same function as before, thankfully.
}

// This function tells us whether a tile is next to a lake.

int nexttolake(planet& world, int x, int y)
{
    int width = world.width();
    int height = world.height();

    int dist = 1;

    for (int i = x - dist; i <= x + dist; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = y - dist; j <= y + dist; j++)
        {
            if (j >= 0 && j <= height && world.lakesurface(ii, j) != 0)
                return (world.lakesurface(ii, j));
        }
    }
    return (0);
}

// This removes any cells of sea that are next to lakes.

void removesealakes(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int surface = world.lakesurface(i, j);

            if (surface != 0)
            {
                for (int k = i - 1; k <= i + 1; k++)
                {
                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (world.nom(k, l) <= sealevel && world.lakesurface(k, l) == 0)
                            world.setlakesurface(k, l, surface);
                    }
                }
            }
        }
    }
}

// This ensures that lakes don't have fragments of themselves nearby.

void connectlakes(planet& world)
{
    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 2; j < height - 1; j++)
        {
            if (world.lakesurface(i, j) != 0)
            {
                int left = i - 1;
                int right = i + 1;
                int up = j - 1;
                int down = j + 1;

                if (left < 0)
                    left = width;

                if (right > width)
                    right = 0;

                int lleft = left - 1;
                int rright = right + 1;
                int uup = up - 1;
                int ddown = down + 1;

                if (lleft < 0)
                    lleft = width;

                if (rright > width)
                    rright = 0;

                int surface = world.lakesurface(i, j);
                int elev = world.nom(i, j);
                int special = world.special(i, j);

                if (world.lakesurface(right, j) == 0 && world.lakesurface(rright, j) == surface)
                {
                    world.setnom(right, j, elev);
                    world.setspecial(right, j, special);
                    world.setlakesurface(right, j, surface);
                }

                if (world.lakesurface(left, j) == 0 && world.lakesurface(lleft, j) == surface)
                {
                    world.setnom(left, j, elev);
                    world.setspecial(left, j, special);
                    world.setlakesurface(left, j, surface);
                }

                if (world.lakesurface(i, up) == 0 && world.lakesurface(i, uup) == surface)
                {
                    world.setnom(i, up, elev);
                    world.setspecial(i, up, special);
                    world.setlakesurface(i, uup, surface);
                }

                if (world.lakesurface(i, down) == 0 && world.lakesurface(i, ddown) == surface)
                {
                    world.setnom(i, down, elev);
                    world.setspecial(i, up, special);
                    world.setlakesurface(i, ddown, surface);
                }
            }
        }
    }
}

// Select the river source cell for a lake or terrain water feature.
void makelakestartpoint(planet& world, vector<vector<int>>& thislake, int lakeno, int leftx, int lefty, int rightx, int righty)
{
    twointegers startpoint;
    startpoint.x = -1;
    startpoint.y = -1;

    int crount = 0;

    do
    {
        for (int i = leftx; i <= rightx; i++)
        {
            for (int j = lefty; j <= righty; j++)
            {
                if (thislake[i][j] == lakeno)
                {
                    if (lakeoutline(world, thislake, lakeno, i, j) == 1) // If this is on the edge of the lake
                    {
                        if (random(1, 3) == 1)
                        {
                            startpoint.x = i;
                            startpoint.y = j;
                        }
                    }
                }
            }
        }

        crount++;

        if (crount > 1000) // Something's gone horribly wrong
            return;

    } while (startpoint.x == -1);

    world.setlakestart(startpoint.x, startpoint.y, 1);
}
