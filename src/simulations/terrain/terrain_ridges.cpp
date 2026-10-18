// Ridge connectivity and mountain bases.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/terrain/terrain_ridges.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;

namespace
{
int getdestinationland(const planet& world, mesh::CellId c, int dir);
}

// This function removes any mountains that aren't over land.

void terrain::detail::removefloatingmountains(planet& world)
{
    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (world.sea(c) == 1)
        {
            if (world.mountainridge(c) != 0)
            {
                for (int dir = 1; dir <= 8; dir++)
                    terrain::deleteridge(world, c, dir);
            }

            world.setcraterrim(c, 0);
        }
    }
}

// This function ensures that all mountain ridges connect up properly, removing rogue ones.

void cleanmountainridges(planet& world)
{
    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (world.mountainheightfloat(c) != 0 && world.sea(c) == 0)
        {
            for (int dir = 1; dir <= 8; dir++)
            {
                if (terrain::getridge(world, c, dir) == 1 && getdestinationland(world, c, dir) == 0) // If there is a ridge going in this direction from this point
                    terrain::deleteridge(world, c, dir);
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

// The same on the cell grid.

int terrain::getridge(const planet& world, mesh::CellId c, int slot)
{
    return (world.mountainridge(c) & (1 << (slot - 1))) != 0 ? 1 : 0;
}

int terrain::getoceanridge(const planet& world, mesh::CellId c, int slot)
{
    return (world.oceanridges(c) & (1 << (slot - 1))) != 0 ? 1 : 0;
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

    terrain::deleteridge(world, world.grid().cellof(x, y), dir);
}

void terrain::deleteridge(planet& world, mesh::CellId c, int slot)
{
    if (getridge(world, c, slot) != 1) // If there isn't a ridge going this way
        return;

    int currentridge = world.mountainridge(c) - terrain::detail::getcode(slot);

    world.setmountainridge(c, currentridge);

    if (currentridge == 0)
        world.setmountainheight(c, 0.0f);

    // Now remove the other end, at the neighbour in its own frame.

    const std::int32_t neighbour = world.grid().neighbour(c, slot);

    if (neighbour < 0)
        return;

    const mesh::CellId n = static_cast<mesh::CellId>(neighbour);
    const int back = world.grid().reverse(c, slot);

    if (getridge(world, n, back) != 1) // If there isn't a ridge going in this direction
        return;

    currentridge = world.mountainridge(n) - terrain::detail::getcode(back);

    world.setmountainridge(n, currentridge);

    if (currentridge == 0)
        world.setmountainheight(n, 0.0f);
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

    terrain::deleteoceanridge(world, world.grid().cellof(x, y), dir);
}

void terrain::deleteoceanridge(planet& world, mesh::CellId c, int slot)
{
    if (getoceanridge(world, c, slot) != 1) // If there isn't a ridge going this way
        return;

    world.setoceanridges(c, world.oceanridges(c) - terrain::detail::getcode(slot));

    // Now remove the other end, at the neighbour in its own frame.

    const std::int32_t neighbour = world.grid().neighbour(c, slot);

    if (neighbour < 0)
        return;

    const mesh::CellId n = static_cast<mesh::CellId>(neighbour);
    const int back = world.grid().reverse(c, slot);

    if (getoceanridge(world, n, back) != 1) // If there isn't a ridge going in this direction
        return;

    world.setoceanridges(n, world.oceanridges(n) - terrain::detail::getcode(back));
}

namespace
{
// This function tells whether there is land in the specified direction from this point.

int getdestinationland(const planet& world, mesh::CellId c, int dir)
{
    const std::int32_t neighbour = world.grid().neighbour(c, dir);

    if (neighbour < 0)
        return(0);

    // The neighbour is where this ridge is pointing to.

    if (world.sea(static_cast<mesh::CellId>(neighbour)) == 1)
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
                float thismountain = world.mountainheightfloat(i, j);

                thismountain = thismountain * mountainheightfactor;

                world.setmountainheight(i, j, thismountain);
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

            float mheights = world.mountainheightfloat(i, j);

            if (mheights == 0 && world.strato(i, j) == 0)
            {
                mheights = world.volcanofloat(i, j);
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
                world.setnomfloat(i, j, world.nomfloat(i, j) + mountaindrainage[i][j]);
        }
    }
}
