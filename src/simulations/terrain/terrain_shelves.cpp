// Continental shelf masks and their connected gaps.
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
void removeshelfgaps(planet& world, std::vector<std::vector<bool>>& shelves);
int nonshelfareacheck(planet& world, std::vector<std::vector<bool>>& shelves, std::vector<std::vector<bool>>& checked, int startx, int starty);
}

// This function makes continental shelves.

void makecontinentalshelves(planet& world, vector<vector<bool>>& shelves, int pointdist)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int midelev = maxelev / 2;

    // First, sort out the outline.

    vector<vector<bool>> outline(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f; //0.4f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    int maxwarp = 60;
    int warpdiv = maxelev / maxwarp;
    int maxradius = 20;
    int raddiv = maxelev / maxradius;

    twofloats pt, mm1, mm2, mm3;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.outline(i, j) == 1)
            {
                // First, work out where the centre of this circle will be. It's getting offset to make the shelves more interesting.

                int warpdistx = 0;
                int warpdisty = 0;

                int warpx = i + width / 4;

                if (warpx > width)
                    warpx = warpx - width;

                int warpy = j;

                if (fractal[warpx][warpy] > maxelev / 2)
                    warpdistx = (fractal[warpx][warpy] - midelev) / warpdiv;
                else
                    warpdistx = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

                int ii = i + warpdistx;

                if (ii<0 || ii>width)
                    wrap(ii, width);

                warpx = warpx + width / 2;

                if (warpx > width)
                    warpx = warpx - width;

                if (fractal[warpx][warpy] > maxelev / 2)
                    warpdisty = (fractal[warpx][warpy] - midelev) / warpdiv;
                else
                    warpdisty = 0 - (midelev - fractal[warpx][warpy]) / warpdiv;

                int jj = j + warpdisty;

                // ii and jj are the centre of the circle we're drawing. But we want to draw a line from the actual outline to that centre, to ensure that there isn't a gap in the shelf.

                mm1.x = (float) i;
                mm1.y = (float) j;

                mm2.x = (float) (i + ii) / 2;
                mm2.y = (float) (j + jj) / 2;

                mm3.x = (float) ii;
                mm3.y = (float) jj;

                for (int n = 1; n <= 2; n++) // Two halves of the curve.
                {
                    for (float t = 0.0f; t <= 1.0f; t = t + 0.01f)
                    {
                        if (n == 1)
                            pt = curvepos(mm1, mm1, mm2, mm3, t);
                        else
                            pt = curvepos(mm1, mm2, mm3, mm3, t);

                        int x = (int)pt.x;
                        int y = (int)pt.y;

                        if (x<0 || x>width)
                            x = wrap(x, width);

                        if (y < 0)
                            y = 0;

                        if (y > height)
                            y = height;

                        outline[x][y] = 1;
                    }
                }

                // Now we can draw a circle.

                if (ii<0 || ii>width)
                    wrap(ii, width);

                if (jj < 0)
                    jj = 0;

                if (jj > height)
                    jj = height;

                int radius = fractal[i][j] / raddiv;

                for (int x = -radius; x <= radius; x++)
                {
                    int xx = ii + x;

                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    for (int y = -radius; y <= radius; y++)
                    {
                        int yy = jj + y;

                        if (yy >= 0 && yy <= height)
                        {
                            if (x * x + y * y < radius * radius + radius)
                                outline[xx][yy] = 1;
                        }
                    }
                }
            }
        }
    }

    // Now, we need a voronoi map.

    vector<vector<short>> voronoi(width + 1, vector<short>(height + 1, 0));

    makeshelvesvoronoi(world, voronoi, outline, pointdist);

    // Now, every panel of the voronoi map that has any outline in it is continental shelf.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (outline[i][j] == 1 && shelves[i][j] == 0)
            {
                int thiscell = voronoi[i][j];

                for (int k = i - 50; k <= i + 50; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 50; l <= j + 50; l++)
                    {
                        if (l >= 0 && l < height)
                        {
                            if (voronoi[kk][l] == thiscell)
                                shelves[kk][l] = 1;
                        }
                    }
                }
            }
        }
    }

    // And anywhere that's land is continental shelf too, technically.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (world.sea(i, j) == 0)
                shelves[i][j] = 1;
        }
    }

    // Now add a circles of shelf to the shelves, to make the edges a bit more varied.

    vector<vector<bool>> circles(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int circlechance = 40;
    int mindist = 5; // Minimum distance to non-shelf sea.
    int radius = 8;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (shelves[i][j] == 1 && random(1, circlechance) == 1)
            {
                bool nearsea = 0;

                for (int k = -mindist; k <= mindist; k++)
                {
                    for (int l = -mindist; l <= mindist; l++)
                    {
                        int ii = i + k;
                        int jj = j + l;

                        if (ii<0 || ii>width)
                            ii = wrap(ii, width);

                        if (jj >= 0 && jj <= height)
                        {
                            if (shelves[ii][jj] == 0)
                            {
                                nearsea = 1;
                                k = mindist;
                                l = mindist;
                            }
                        }
                    }
                }

                if (nearsea == 0)
                {

                    for (int x = -radius; x <= radius; x++)
                    {
                        int xx = i + x;

                        if (xx<0 || xx>width)
                            xx = wrap(xx, width);

                        for (int y = -radius; y <= radius; y++)
                        {
                            int yy = j + y;

                            if (yy >= 0 && yy <= height)
                            {
                                if (x * x + y * y < radius * radius + radius)
                                    circles[xx][yy] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (circles[i][j] == 1)
                shelves[i][j] = 1;
        }
    }

    // Now we need to remove any holes in the continental shelves!

    removeshelfgaps(world, shelves);
}

namespace
{
// This function removes any gaps in the continental shelves

void removeshelfgaps(planet& world, vector<vector<bool>>& shelves)
{
    int width = world.width();
    int height = world.height();

    int minseasize = (width * height) / 4;

    vector<vector<bool>> checked(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0)); // This array marks which sea areas have already been scanned.

    vector<vector<int>> area(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This array marks the area of each bit of sea.

    int seax = -1;
    int seay = -1;

    // First, we find the area of each section of non-shelf and mark them on the area array.

    for (int i = width - 1; i > 0; i--)
    {
        for (int j = height - 1; j > 0; j--)
        {
            if (shelves[i][j] == 0 && checked[i][j] == 0 && area[i][j] == 0)
            {
                int size = nonshelfareacheck(world, shelves, checked, i, j);

                for (int k = 0; k <= width; k++) // Clear the checked array and mark all adjoining points with this area.
                {
                    for (int l = 0; l <= height; l++)
                    {
                        if (checked[k][l] == 1)
                        {
                            area[k][l] = size;
                            checked[k][l] = 0;
                        }
                    }
                }
            }
        }
    }

    // Now we just find a point that's in the largest area of non-shelf.

    int largest = 0;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (area[i][j] > largest)
            {
                largest = area[i][j];
                seax = i;
                seay = j;
            }
        }
    }

    // Now we simply remove any sea tiles that aren't part of the largest area of sea.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (area[i][j] != largest)
                shelves[i][j] = 1;
        }
    }
}
}

namespace
{
// This tells how large a given area of non-shelf sea is.

int nonshelfareacheck(planet& world, vector<vector<bool>>& shelves, vector<vector<bool>>& checked, int startx, int starty)
{
    int width = world.width();
    int height = world.height();

    int total = 0;

    int row[] = { -1,0,0,1 };
    int col[] = { 0,-1,1,0 };

    twointegers node;

    node.x = startx;
    node.y = starty;

    queue<twointegers> q; // Create a queue
    q.push(node); // Put our starting node into the queue

    while (q.empty() != 1) // Keep going while there's anything left in the queue
    {
        node = q.front(); // Take the node at the front of the queue
        q.pop(); // And then pop it out of the queue

        if (node.x >= 0 && node.x <= width && node.y >= 0 && node.y <= height && shelves[node.x][node.y] == 0)
        {
            total++;
            checked[node.x][node.y] = 1;

            for (int k = 0; k < 4; k++) // Look at the four neighbouring nodes in turn
            {
                twointegers nextnode;

                nextnode.x = node.x + row[k];
                nextnode.y = node.y + col[k];

                if (nextnode.x > width)
                    nextnode.x = 0;

                if (nextnode.x < 0)
                    nextnode.x = width;

                if (nextnode.y >= 0 && nextnode.y < height)
                {
                    if (shelves[nextnode.x][nextnode.y] == 0 && checked[nextnode.x][nextnode.y] == 0) // If this node is sea
                    {
                        checked[nextnode.x][nextnode.y] = 1;
                        q.push(nextnode); // Put that node onto the queue
                    }
                }
            }
        }
    }
    return total;
}
}
