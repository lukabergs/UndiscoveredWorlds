#include <vector>

#include "core/classes.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"
#include "simulations/hydrology/lakes.hpp"

using namespace std;

namespace
{
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
}

// This function tells us whether a tile is next to a lake.

float nexttolake(planet& world, int x, int y)
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
            if (j >= 0 && j <= height && world.lakesurfacefloat(ii, j) != 0)
                return world.lakesurfacefloat(ii, j);
        }
    }
    return 0.0f;
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
