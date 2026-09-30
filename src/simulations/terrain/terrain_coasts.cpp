// Sea connectivity, coastal elevation, and island cleanup.
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
int areacheck(planet& world, std::vector<std::vector<bool>>& checked, int startx, int starty);
}

// This function removes any small seas from the map. level is level to fill them in to.

void removesmallseas(planet& world, int minseasize, int level)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int j = 0; j <= height; j++)
    {
        int amount = (world.nom(width - 1, j) + world.nom(0, j)) / 2;
        world.setnom(width, j, amount);
    }

    vector<vector<bool>> checked(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0)); // This array marks which sea areas have already been scanned.

    vector<vector<int>> area(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This array marks the area of each bit of sea.

    int seax = -1;
    int seay = -1;

    // First, we find the area of each section of sea and mark them on the area array.

    for (int i = width - 1; i > 0; i--)
    {
        for (int j = height - 1; j > 0; j--)
        {
            if (world.nom(i, j) <= sealevel && checked[i][j] == 0 && area[i][j] == 0)
            {
                int size = areacheck(world, checked, i, j);

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

    // Now turn all the points whose area is too small into land.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) <= sealevel && area[i][j] < minseasize)
                world.setnom(i, j, level);
        }
    }
}

namespace
{
// This tells how large a given area of sea is.

int areacheck(planet& world, vector<vector<bool>>& checked, int startx, int starty)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

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

        if (node.x >= 0 && node.x <= width && node.y >= 0 && node.y <= height && world.nom(node.x, node.y) <= sealevel)
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
                    if (world.nom(nextnode.x, nextnode.y) <= sealevel && checked[nextnode.x][nextnode.y] == 0) // If this node is sea
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

// This function removes straight edges on the coastlines.

void removestraights(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int maxlength = 4; // Any straights longer than this will be disrupted.

    int minsize = 2;
    int maxsize = 4; // Range of sizes for the chunks to be cut out of the coasts.

    for (int i = 0; i <= width; i++) // First do vertical/horizontal lines.
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.map(i, j) <= sealevel)
            {
                if (j > 0)
                {
                    if (world.map(i, j - 1) > sealevel) // Sea to the south, land to the north
                    {
                        bool keepgoing = 1;
                        int x = i;
                        int xx = x;
                        int y = j;

                        do
                        {
                            if (world.map(x, y) > sealevel)
                                keepgoing = 0;

                            if (world.map(x, y - 1) <= sealevel)
                                keepgoing = 0;

                            x++;

                            xx = x;

                            if (xx > width)
                                xx = wrap(xx, width);

                        } while (keepgoing == 1);

                        if (x > i + maxlength) // If this line is too long
                        {
                            x = i + (x - i) / 2; // Move to the middle of the line

                            int size = random(minsize, maxsize);
                            disruptseacoastline(world, x, y, world.map(i, j), 0, size);
                        }
                    }
                }

                if (j < height)
                {
                    if (world.map(i, j + 1) > sealevel) // Sea to the north, land to the south
                    {
                        bool keepgoing = 1;
                        int x = i;
                        int xx = x;
                        int y = j;

                        do
                        {
                            if (world.map(x, y) > sealevel)
                                keepgoing = 0;

                            if (world.map(x, y + 1) <= sealevel)
                                keepgoing = 0;

                            x++;

                            xx = x;

                            if (xx > width)
                                xx = wrap(xx, width);

                        } while (keepgoing == 1);

                        if (x > i + maxlength) // If this line is too long
                        {
                            if (x > i + maxlength) // If this line is too long
                            {
                                x = i + (x - i) / 2; // Move to the middle of the line

                                int size = random(minsize, maxsize);
                                disruptseacoastline(world, x, y, world.map(i, j), 0, size);
                            }
                        }
                    }
                }

                if (world.map(i + 1, j) > sealevel) // Sea to the west, land to the east
                {
                    bool keepgoing = 1;
                    int x = i;
                    int y = j;
                    int xx = x + 1;

                    if (xx > width)
                        xx = 0;

                    do
                    {
                        if (world.map(x, y) > sealevel)
                            keepgoing = 0;

                        if (world.map(xx, y) <= sealevel)
                            keepgoing = 0;

                        y++;

                        if (y > height - 1)
                            keepgoing = 0;

                    } while (keepgoing == 1);

                    if (y > j + maxlength) // If this line is too long
                    {
                        y = j + (y - j) / 2; // Move to the middle of the line

                        int size = random(minsize, maxsize);
                        disruptseacoastline(world, x, y, world.map(i, j), 0, size);
                    }
                }

                if (world.map(i - 1, j) > sealevel) // Sea to the east, land to the west
                {
                    bool keepgoing = 1;
                    int x = i;
                    int y = j;
                    int xx = i - 1;

                    if (xx < 0)
                        xx = width;

                    do
                    {
                        if (world.map(x, y) > sealevel)
                            keepgoing = 0;

                        if (world.map(x - 1, y) <= sealevel)
                            keepgoing = 0;

                        y++;

                        if (y > height - 1)
                            keepgoing = 0;

                    } while (keepgoing == 1);

                    if (y > j + maxlength) // If this line is too long
                    {
                        y = j + (y - j) / 2; // Move to the middle of the line

                        int size = random(minsize, maxsize);
                        disruptseacoastline(world, x, y, world.map(i, j), 0, size);
                    }
                }
            }
        }
    }

    for (int i = 0; i <= width; i++) // Now do diagonals.
    {
        for (int j = 1; j < height; j++)
        {
            if (world.map(i, j) > sealevel)
            {
                if (northwestlandonly(world, i, j) == 1)
                {
                    int ii = i + 1;

                    if (ii > width)
                        ii = 0;

                    int jj = j - 1;

                    if (northwestlandonly(world, ii, jj) == 1)
                    {
                        int iii = i - 1;

                        if (iii < 0)
                            iii = width;

                        int jjj = j + 1;

                        if (northwestlandonly(world, iii, jjj) == 1)
                        {
                            int choose = random(1, 3);

                            if (choose == 1)
                                world.setnom(i, j, world.nom(ii, jjj));

                            if (choose == 2)
                                world.setnom(ii, j, world.nom(i, j));

                            if (choose == 3)
                                world.setnom(i, jjj, world.nom(i, j));
                        }
                    }
                }

                if (northeastlandonly(world, i, j) == 1)
                {
                    int ii = i + 1;

                    if (ii > width)
                        ii = 0;

                    int jj = j + 1;

                    if (northeastlandonly(world, ii, jj) == 1)
                    {
                        int iii = i - 1;

                        if (iii < 0)
                            iii = width;

                        int jjj = j - 1;

                        if (northeastlandonly(world, iii, jjj) == 1)
                        {
                            int choose = random(1, 3);

                            if (choose == 1)
                                world.setnom(i, j, world.nom(iii, jj));

                            if (choose == 2)
                                world.setnom(iii, j, world.nom(i, j));

                            if (choose == 3)
                                world.setnom(i, jj, world.nom(i, j));
                        }
                    }
                }

                if (southwestlandonly(world, i, j) == 1)
                {
                    int ii = i + 1;

                    if (ii > width)
                        ii = 0;

                    int jj = j + 1;

                    if (southwestlandonly(world, ii, jj) == 1)
                    {
                        int iii = i - 1;

                        if (iii < 0)
                            iii = width;

                        int jjj = j - 1;

                        if (southwestlandonly(world, iii, jjj) == 1)
                        {
                            int choose = random(1, 3);

                            if (choose == 1)
                                world.setnom(i, j, world.nom(ii, jjj));

                            if (choose == 2)
                                world.setnom(ii, j, world.nom(i, j));

                            if (choose == 3)
                                world.setnom(i, jjj, world.nom(i, j));
                        }
                    }



                }

                if (southeastlandonly(world, i, j) == 1)
                {
                    int ii = i + 1;

                    if (ii > width)
                        ii = 0;

                    int jj = j - 1;

                    if (southeastlandonly(world, ii, jj) == 1)
                    {
                        int iii = i - 1;

                        if (iii < 0)
                            iii = width;

                        int jjj = j + 1;

                        if (southeastlandonly(world, iii, jjj) == 1)
                        {
                            int choose = random(1, 3);

                            if (choose == 1)
                                world.setnom(i, j, world.nom(iii, jj));

                            if (choose == 2)
                                world.setnom(iii, j, world.nom(i, j));

                            if (choose == 3)
                                world.setnom(i, jj, world.nom(i, j));
                        }
                    }
                }
            }
        }
    }
}

// This function cuts a slice out of the coastline, to try to make it more interesting.

void disruptseacoastline(planet& world, int centrex, int centrey, int avedepth, bool raise, int size)
{
    int width = world.width();
    int height = world.height();

    int half = size / 2;

    for (int i = -half; i <= half; i++)
    {
        for (int j = -half; j <= half; j++)
        {
            int ii = centrex + i;
            int jj = centrey + j;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            if (jj >= 0 && jj <= height)
            {
                if (world.nom(ii, jj) > avedepth)
                {
                    if (i * i + j * j < half * half + half)
                        world.setnom(ii, jj, avedepth);
                }
            }
        }
    }
}

// This function ensures that channels of sea are at least two pixels wide.

void widenchannels(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int i = 1; i < width; i++) // First remove any one-pixel channels.
    {
        for (int j = 1; j < height; j++)
        {
            if (world.sea(i, j))
            {
                if (world.sea(i, j - 1) == 0 && world.sea(i, j + 1) == 0)
                    world.setnom(i, j - 1, sealevel - 15);

                if (world.sea(i - 1, j) == 0 && world.sea(i + 1, j) == 0)
                    world.setnom(i - 1, j, sealevel - 15);

                if (world.sea(i - 1, j - 1) == 0 && world.sea(i + 1, j + 1) == 0)
                    world.setnom(i - 1, j - 1, sealevel - 15);

                if (world.sea(i - 1, j + 1) == 0 && world.sea(i + 1, j - 1) == 0)
                    world.setnom(i - 1, j + 1, sealevel - 15);
            }
        }
    }

    for (int i = 2; i < width - 1; i++) // Now remove any two-pixel channels.
    {
        for (int j = 2; j < height - 1; j++)
        {
            if (world.sea(i, j))
            {
                if (world.sea(i, j - 1) == 0 && world.sea(i, j + 1) == 1 && world.sea(i, j + 2) == 0)
                    world.setnom(i, j, sealevel - 15);

                if (world.sea(i - 1, j) == 0 && world.sea(i + 1, j) == 1 && world.sea(i + 2, j) == 0)
                    world.setnom(i + 1, j, sealevel - 15);
            }
        }
    }
}

void loweroceans(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) && world.nom(i, j) > sealevel - 170)
            {
                bool nearland = 0;

                for (int ii = i - 2; ii <= i + 2; ii++)
                {
                    int iii = ii;

                    if (iii<0 || iii>width)
                        iii = wrap(iii, width);

                    for (int jj = j - 2; jj <= j + 2; jj++)
                    {
                        if (jj >= 0 && jj <= height)
                        {
                            if (world.sea(iii, jj) == 0)
                            {
                                nearland = 1;

                                ii = i + 2;
                                jj = j + 2;
                            }
                        }
                    }
                }

                if (nearland == 0)
                    world.setnom(i, j, sealevel - 170 + randomsign(random(0, 10)));
            }
        }
    }
}

// This function forces land and sea at coastlines to normalise towards certain values, to improve the appearance of the coastline. Note that landheight and seadepth are relative to the sea level, not absolute.

void normalisecoasts(planet& world, int landheight, int seadepth, int severity)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int depth, elev;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.volcano(i, j) == 0 && world.mountainisland(i, j) == 0)
            {
                if (world.coast(i, j))
                {
                    depth = sealevel - world.nom(i, j);
                    int newdepth = (depth + seadepth * severity) / (severity + 1);

                    world.setnom(i, j, sealevel - newdepth);
                }

                if (world.outline(i, j))
                {
                    elev = world.nom(i, j) - sealevel;
                    int newelev = (elev + landheight * severity) / (severity + 1);

                    world.setnom(i, j, sealevel + newelev);
                }
            }
        }
    }

    // Now we need to smooth the seabeds around them.

    vector<vector<bool>> done(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    int amount = 2;
    int amount2 = 2;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.coast(i, j) == 1)
            {
                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (done[kk][l] == 0 && world.sea(kk, l) == 1 && world.coast(kk, l) == 0)
                            {
                                int total = 0;
                                int crount = 0;

                                for (int m = kk - amount2; m <= kk + amount2; m++)
                                {
                                    int mm = m;

                                    if (mm<0 || mm>width)
                                        mm = wrap(mm, width);

                                    for (int n = l - amount2; n <= l + amount2; n++)
                                    {
                                        if (n >= 0 && n <= height && world.sea(mm, n) == 1)
                                        {
                                            crount++;
                                            total = total + world.nom(mm, n);
                                        }
                                    }
                                }

                                if (total != 0)
                                {
                                    int newnom = total / crount;
                                    world.setnom(kk, l, newnom);
                                    done[kk][l] = 1;

                                    world.setnoshade(kk, l, 1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Do it again, going the other way.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            done[i][j] = 0;
    }

    for (int i = width; i >= 0; i--)
    {
        for (int j = height; j >= 0; j--)
        {
            if (world.coast(i, j) == 1)
            {
                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (done[kk][l] == 0 && world.sea(kk, l) == 1 && world.coast(kk, l) == 0)
                            {
                                int total = 0;
                                int crount = 0;

                                for (int m = kk - amount2; m <= kk + amount2; m++)
                                {
                                    int mm = m;

                                    if (mm<0 || mm>width)
                                        mm = wrap(mm, width);

                                    for (int n = l - amount2; n <= l + amount2; n++)
                                    {
                                        if (n >= 0 && n <= height && world.sea(mm, n) == 1)
                                        {
                                            crount++;
                                            total = total + world.nom(mm, n);
                                        }
                                    }
                                }

                                if (total != 0)
                                {
                                    int newnom = total / crount;
                                    world.setnom(kk, l, newnom);
                                    done[kk][l] = 1;
                                    world.setnoshade(kk, l, 1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// This notes down all the very small islands on the map.

void checkislands(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 1; j < height; j++)
        {
            if (world.sea(i, j) == 0)
            {
                int crount = -1; // Because there will definitely be one positive, i.e. the current cell.

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (world.sea(kk, l) == 0)
                            crount++;
                    }
                }

                if (crount == 0) // This is a one-tile island!
                {
                    world.setisland(i, j, 1);

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (world.sea(kk, l) == 1)
                                world.setnoshade(kk, l, 1);
                        }
                    }
                }
            }
            else
            {
                int volcano = world.volcano(i, j);

                if (volcano != 0)
                {
                    bool extinct = 0;

                    if (volcano < 0)
                    {
                        extinct = 1;
                        volcano = 0 - volcano;
                    }

                    int totalheight = world.nom(i, j) + volcano;

                    if (totalheight > sealevel)
                    {
                        volcano = totalheight - sealevel;

                        terrain::detail::makemountainisland(world, i, j, volcano);

                        if (extinct == 1)
                            volcano = 0 - volcano;

                        world.setvolcano(i, j, volcano);
                    }
                }
            }
        }
    }
}

// This extends the noshade areas by one cell in each direction.

void extendnoshade(planet& world)
{
    int width = world.width();
    int height = world.height();

    vector<vector<bool>> extra(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0)); // This is to mark all cells touching the no-shade ones as themselves no-shade.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.noshade(i, j) == 1)
            {
                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height && l >= 0 && l <= height && world.sea(kk, l) == 1)
                            extra[kk][l] = 1;
                    }
                }
            }
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (extra[i][j] == 1)
                world.setnoshade(i, j, 1);
        }
    }
}
