// Ridge connectivity, mountain bases, and fjord relief.
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
int getdestinationland(planet& world, int x, int y, int dir);
}

// This function removes any mountains that aren't over land.

void terrain::detail::removefloatingmountains(planet& world)
{
    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1)
            {
                if (world.mountainridge(i, j) != 0)
                {
                    for (int dir = 1; dir <= 8; dir++)
                        deleteridge(world, i, j, dir);
                }

                world.setcraterrim(i, j, 0);
            }
        }
    }
}

// This function ensures that all mountain ridges connect up properly, removing rogue ones.

void cleanmountainridges(planet& world)
{
    int width = world.width();
    int height = world.height();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.mountainheight(i, j) != 0 && world.sea(i, j) == 0)
            {
                for (int dir = 1; dir <= 8; dir++)
                {
                    if (getridge(world, i, j, dir) == 1 && getdestinationland(world, i, j, dir) == 0) // If there is a ridge going in this direction from this point
                        deleteridge(world, i, j, dir);
                }
            }
        }
    }
}

// This function translates a mountain direction into a binary code.

int terrain::detail::getcode(int dir)
{
    if (dir == 1)
        return(1);

    if (dir == 2)
        return(2);

    if (dir == 3)
        return(4);

    if (dir == 4)
        return(8);

    if (dir == 5)
        return(16);

    if (dir == 6)
        return(32);

    if (dir == 7)
        return(64);

    if (dir == 8)
        return(128);

    return(0);
}

// This function tells whether a ridge goes from the specified tile in the specified direction.

int getridge(planet& world, int x, int y, int dir)
{
    bool check = (world.mountainridge(x, y) & (1 << (dir - 1))) != 0;

    if (check == 1)
        return(1);

    return(0);
}

// Same thing, but on a specified vector.

int getridge(vector<vector<int>>& arr, int x, int y, int dir)
{
    bool check = (arr[x][y] & (1 << (dir - 1))) != 0;

    if (check == 1)
        return(1);

    return(0);
}

// Same thing, but for an ocean ridge.

int getoceanridge(planet& world, int x, int y, int dir)
{
    bool check = (world.oceanridges(x, y) & (1 << (dir - 1))) != 0;

    if (check == 1)
        return(1);

    return(0);
}

// This function deletes a ridge going from the specified tile in the specified direction.

void deleteridge(planet& world, int x, int y, int dir)
{
    int width = world.width();
    int height = world.height();

    if (x<0 || x>width)
        x = wrap(x, width);

    if (y<0 || y>height)
        return;

    if (getridge(world, x, y, dir) != 1) // If there isn't a ridge going this way
        return;

    int code = terrain::detail::getcode(dir);

    int currentridge = world.mountainridge(x, y);

    currentridge = currentridge - code;

    world.setmountainridge(x, y, currentridge);

    if (currentridge == 0)
        world.setmountainheight(x, y, 0);

    // Now remove one going the other way.

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
        return;

    dir = dir + 4;

    if (dir > 8)
        dir = dir - 8;

    if (getridge(world, x, y, dir) != 1) // If there isn't a ridge going in this direction
        return;

    code = terrain::detail::getcode(dir);

    currentridge = world.mountainridge(x, y);

    currentridge = currentridge - code;

    world.setmountainridge(x, y, currentridge);

    if (currentridge == 0)
        world.setmountainheight(x, y, 0);
}

// Same thing, but on a specified vector.

void deleteridge(planet& world, vector<vector<int>>& ridgesarr, vector<vector<int>>& heightsarr, int x, int y, int dir)
{
    int width = world.width();
    int height = world.height();

    if (x<0 || x>width)
        x = wrap(x, width);

    if (y<0 || y>height)
        return;

    if (getridge(ridgesarr, x, y, dir) != 1) // If there isn't a ridge going this way
        return;

    int code = terrain::detail::getcode(dir);

    int currentridge = ridgesarr[x][y];

    currentridge = currentridge - code;

    ridgesarr[x][y] = currentridge;

    if (currentridge == 0)
        heightsarr[x][y] = 0;

    // Now remove one going the other way.

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
        return;

    dir = dir + 4;

    if (dir > 8)
        dir = dir - 8;

    if (getridge(ridgesarr, x, y, dir) != 1) // If there isn't a ridge going in this direction
        return;

    code = terrain::detail::getcode(dir);

    currentridge = ridgesarr[x][y];

    currentridge = currentridge - code;

    ridgesarr[x][y] = currentridge;

    if (currentridge == 0)
        heightsarr[x][y] = 0;
}

// Same thing, but for an ocean ridge.

void terrain::detail::deleteoceanridge(planet& world, int x, int y, int dir)
{
    int width = world.width();
    int height = world.height();

    if (x<0 || x>width)
        x = wrap(x, width);

    if (y<0 || y>height)
        return;

    if (getoceanridge(world, x, y, dir) != 1) // If there isn't a ridge going this way
        return;

    int code = terrain::detail::getcode(dir);

    int currentridge = world.oceanridges(x, y);

    currentridge = currentridge - code;

    world.setoceanridges(x, y, currentridge);

    // Now remove one going the other way.

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
        return;

    dir = dir + 4;

    if (dir > 8)
        dir = dir - 8;

    if (getoceanridge(world, x, y, dir) != 1) // If there isn't a ridge going in this direction
        return;

    code = terrain::detail::getcode(dir);

    currentridge = world.oceanridges(x, y);

    currentridge = currentridge - code;

    world.setoceanridges(x, y, currentridge);
}

namespace
{
// This function tells whether there is land in the specified direction from this point.

int getdestinationland(planet& world, int x, int y, int dir)
{
    int newx = x;
    int newy = y;

    if (dir == 8 || dir == 1 || dir == 2)
        newy--;

    if (dir == 4 || dir == 5 || dir == 6)
        newy++;

    if (dir == 2 || dir == 3 || dir == 4)
        newx++;

    if (dir == 6 || dir == 7 || dir == 8)
        newx--;

    if (newy<0 || newy>world.height())
        return(0);

    if (newx<0 || newx>world.width())
        newx = wrap(newx, world.width());

    // newx and newy are the coordinates of where this ridge is pointing to.

    if (world.sea(newx, newy) == 1)
        return(0);

    return(1);
}
}

// This function raises the land beneath mountains.

void raisemountainbases(planet& world, vector<vector<int>>& mountaindrainage, vector<vector<bool>>& OKmountains)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();
    float gravity = world.gravity();

    int maxradius = 6;
    
    float heightreduce = 0.5f; // Multiply the extra height by this each time.
    float extraheightreduce = 0.045f; // Each time, reduce the heightreduce by this.

    if (gravity < 1.0f) // Low gravity means wider areas of higher ground.
    {
        float diff = 1.0f - gravity;
        diff = diff * 10.0f;

        maxradius = maxradius + (int)diff;

        if (maxradius > 12)
            maxradius = 12;

        heightreduce = heightreduce + diff / 10.0f;

        if (heightreduce > 0.8f)
            heightreduce = 0.8f;      
    }

    if (gravity > 1.0f) // High gravity means smaller areas of higher ground.
    {
        float diff = gravity - 1.0f;

        maxradius = maxradius - (int)(diff * 10.0f);

        if (maxradius < 2)
            maxradius = 2;

        heightreduce = heightreduce - diff;

        if (heightreduce < 0.2f)
            heightreduce = 0.2f;
    }

    int volcanomaxradius = maxradius / 3;

    // First, adjust the heights of the mountains themselves to take account of gravity.

    float mountainheightfactor = 1.0f;

    if (gravity < 1.0f) // Low gravity means higher mountains.
    {
        float diff = 1.0f - gravity;

        mountainheightfactor = mountainheightfactor + diff * 2.0f;
    }

    if (gravity > 1.0f) // High gravity means lower mountains.
    {
        float diff = gravity - 1.0f;

        mountainheightfactor = mountainheightfactor - diff / 2.0f;

        if (mountainheightfactor < 0.2f)
            mountainheightfactor = 0.2f;
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (OKmountains[i][j] == 0) // Don't change the heights of any mountains that the user has imported.
            {
                float thismountain = (float)world.mountainheight(i, j);

                thismountain = thismountain * mountainheightfactor;

                world.setmountainheight(i, j, (int)thismountain);
            }
        }
    }

    // Adjust the maximum elevation for the world.

    float abovesea = (float)maxelev - (float)sealevel;

    abovesea = abovesea * mountainheightfactor;

    if (abovesea > 20000.0f) // Hard maximum, as weirdness happens above this height.
        abovesea = 20000.0f;

    world.setmaxelevation(sealevel + (int)abovesea);

    // Now do the mountain bases themselves.

    bool volcano = 0;

    float lowrandom = 20.0f;
    float highrandom = 50.0f;

    if (gravity > 1.0f)
    {
        lowrandom = lowrandom - gravity * 2.0f;

        if (lowrandom < 2.0f)
            lowrandom = 2.0f;

        highrandom = lowrandom * 2.0f;
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            volcano = 0;

            float mheights = (float) world.mountainheight(i, j);

            if (mheights == 0 && world.strato(i, j) == 0)
            {
                mheights = (float)world.volcano(i, j);
                volcano = 1;
            }

            if (mheights < 0.0f)
                mheights = 0.0f - mheights;

            if (mheights != 0.0f)
            {
                heightreduce = 0.7f; //0.5f;

                float baseheightfraction = (float)random((int)lowrandom, (int)highrandom);

                if (volcano == 1)
                    baseheightfraction = (float)random((int)lowrandom * 5, (int)highrandom * 4);

                baseheightfraction = baseheightfraction / 100.0f;
                float baseheight = mheights * baseheightfraction;

                int thismaxradius = maxradius;

                if (volcano == 1)
                    thismaxradius = volcanomaxradius;

                for (int radius = 1; radius <= thismaxradius; radius++) // Go through ever-increasing circles
                {
                    for (int k = -radius; k <= radius; k++)
                    {
                        for (int l = -radius; l <= radius; l++)
                        {
                            if (k * k + l * l < radius * radius + radius) // If we're within the current circle
                            {
                                int ll = j + l;

                                if (ll >= 0 && ll <= height)
                                {
                                    int kk = i + k;

                                    if (kk<0 || kk>width)
                                        kk = wrap(kk, width);

                                    if (world.sea(kk, ll) == 0 && mountaindrainage[kk][ll] < baseheight)
                                        mountaindrainage[kk][ll] = (int)baseheight;
                                }
                            }
                        }
                    }

                    baseheight = baseheight * heightreduce;
                }
            }
        }
    }

    // Now add the extra heights onto the actual map.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 0)
                world.setnom(i, j, world.nom(i, j) + mountaindrainage[i][j]);
        }
    }
}

// Extend coastal mountains into the sea to form fjords in glaciated regions
// and rocky peninsulas elsewhere. Climate stages call this after temperatures settle.

void addfjordmountains(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int glacialtemp = world.glacialtemp();
    int glacialmountainheight = 20; //200; // In glacial regions, mountains this height or more are guaranteed to spawn fjords.

    int minmountheight = 0; // Ignore mountains lower than this.
    int noglacchance = 8; //3; //15; // Chance of making sea mountain ridges in non-glaciated areas (to add occasional rocky peninsulas).

    vector<vector<int>> justadded(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> previouslyadded(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int n = 1; n <= 2; n++)
    {
        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                if (world.sea(i, j)) // && random(1,maxelev)<world.roughness(i,j))
                {
                    int avetemp = (world.mintemp(i, j) + world.maxtemp(i, j)) / 2;

                    bool fjordreason1 = 0;
                    bool fjordreason2 = 0;

                    if (avetemp < glacialtemp)
                        fjordreason1 = 1; // Because it's cold enough. But this will depend upon the mountains being high enough.

                    if (random(1, noglacchance) == 1 && random(1, maxelev) < (int)world.roughness(i, j))
                        fjordreason2 = 1; // Because it's just one of those ones that appear throughout the map.

                    if (fjordreason1 == 1 || fjordreason2 == 1)
                    {
                        int nexttoland = 0;

                        int roughquotient = random(1, maxelev);

                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (n == 1 && world.sea(kk, l) == 0 && world.mountainisland(kk, l) == 0) // Don't allow these from mountain islands, as that causes weirdness.
                                    {
                                        nexttoland = 1;
                                        k = i + 1;
                                        l = j + 1;
                                    }

                                    if (n != 1 && l >= 0 && l <= height)
                                    {
                                        if (previouslyadded[kk][l] == 1)
                                        {
                                            int roughness = (int)world.roughness(i, j);

                                            if (avetemp < glacialtemp)
                                                roughness = roughness * 5; // This is more likely in glacial areas.

                                            if (roughquotient < roughness) // The rougher the area, the more likely this is
                                            {
                                                nexttoland = 1;
                                                k = i + 1;
                                                l = j + 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (nexttoland == 1)
                        {
                            int highest = 0;

                            int landx = -1;
                            int landy = -1;

                            for (int k = i - 1; k <= i + 1; k++) // Find the highest mountain nearby, if there is one
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - 1; l <= j + 1; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.mountainheight(kk, l) > highest && justadded[kk][l] == 0)
                                        {
                                            highest = world.mountainheight(kk, l);

                                            landx = kk;
                                            landy = l;
                                        }
                                    }
                                }
                            }

                            if (fjordreason1 == 1 && fjordreason2 == 0 && highest < glacialmountainheight) // If these are glacial mountain fjords, the mountains have to be high enough.
                                landx = -1;

                            if (landx != -1 && highest >= minmountheight) // If we found a nearby mountain, add a ridge to it.
                            {
                                int dir1 = getdir(i, j, landx, landy);
                                int dir2 = getdir(landx, landy, i, j);

                                int code1 = terrain::detail::getcode(dir1);
                                int code2 = terrain::detail::getcode(dir2);

                                world.setmountainridge(i, j, world.mountainridge(i, j) + code1);
                                world.setmountainridge(landx, landy, world.mountainridge(landx, landy) + code2);

                                world.setmountainheight(i, j, world.mountainheight(landx, landy));

                                world.setsummerrain(i, j, world.summerrain(landx, landy));
                                world.setwinterrain(i, j, world.winterrain(landx, landy));

                                justadded[i][j] = 1;

                                //world.settest(i,j,world.mountainheight(i,j));
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                if (previouslyadded[i][j] == 0)
                    previouslyadded[i][j] = justadded[i][j];

                justadded[i][j] = 0;
            }
        }
    }
}
