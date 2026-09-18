#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <unordered_set>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "drainage.hpp"

using namespace std;

namespace
{
int seaflowdir(int x, int y, int newx, int newy);
void addriverflowlocked(planet& world, int x, int y, int janadd, int juladd, vector<mutex>& rowlocks);
void tracedropparallel(planet& world, int x, int y, int minimum, vector<mutex>& rowlocks);
void removediagonalrivers(planet& world);
void adddiagonalriverjunctions(planet& world);
void avoidvolcanoes(planet& world);
void removeparallelrivers(planet& world);
void removecrossingrivers(planet& world);

int seaflowdir(int x, int y, int newx, int newy)
{
    if (newx == x && newy < y)
        return 1;

    if (newx > x && newy < y)
        return 2;

    if (newx > x && newy == y)
        return 3;

    if (newx > x && newy > y)
        return 4;

    if (newx == x && newy > y)
        return 5;

    if (newx < x && newy > y)
        return 6;

    if (newx < x && newy == y)
        return 7;

    if (newx < x && newy < y)
        return 8;

    return 0;
}

void addriverflowlocked(planet& world, int x, int y, int janadd, int juladd, vector<mutex>& rowlocks)
{
    lock_guard<mutex> lock(rowlocks[y]);
    world.setriverjan(x, y, world.riverjan(x, y) + janadd);
    world.setriverjul(x, y, world.riverjul(x, y) + juladd);
}

void tracedropparallel(planet& world, int x, int y, int minimum, vector<mutex>& rowlocks)
{
    const int startx = x;
    const int starty = y;
    const int width = world.width();
    const int height = world.height();
    const int maxelevation = world.maxelevation();
    const float riverfactor = world.riverfactor();

    float janload = world.seasonalrainfloat(seasonjanuary, startx, starty);
    float julload = world.seasonalrainfloat(seasonjuly, startx, starty);

    if ((janload + julload) / 2.0f < static_cast<float>(minimum))
        return;

    janload = janload / riverfactor;
    julload = julload / riverfactor;

    float amount = 10.0f;

    if (world.mintemp(startx, starty) < 0)
        amount = 3.0f;

    if (starty < height / 2)
    {
        const float diff = janload / amount;
        janload = janload - diff;
        julload = julload + diff;
    }
    else
    {
        const float diff = julload / amount;
        julload = julload - diff;
        janload = janload + diff;
    }

    if (janload < 0.0f)
        janload = 0.0f;

    if (julload < 0.0f)
        julload = 0.0f;

    const int janadd = static_cast<int>(janload);
    const int juladd = static_cast<int>(julload);

    if (janadd == 0 && juladd == 0)
        return;

    vector<pair<int, int>> path;
    path.reserve(256);

    unordered_set<int> visited;
    visited.reserve(256);

    bool loopdetected = false;
    int loopx = -1;
    int loopy = -1;

    for (;;)
    {
        const int currentkey = y * (width + 1) + x;
        visited.insert(currentkey);
        path.emplace_back(x, y);

        int dir = world.riverdir(x, y);

        if (dir == 8 || dir == 1 || dir == 2)
            y--;

        if (dir == 4 || dir == 5 || dir == 6)
            y++;

        if (dir == 2 || dir == 3 || dir == 4)
            x++;

        if (dir == 6 || dir == 7 || dir == 8)
            x--;

        if (x < 0 || x > width)
            x = wrap(x, width);

        if (y < 0)
            y = 0;

        if (y > height)
            y = height;

        if (world.sea(x, y) == 1)
        {
            int seax = x;
            int seay = y;
            int seadir = dir;

            for (int n = 1; n <= 3; n++)
            {
                path.emplace_back(seax, seay);

                const twointegers newseatile = findseatile(world, seax, seay, seadir);

                if (newseatile.x != -1)
                {
                    seadir = seaflowdir(seax, seay, newseatile.x, newseatile.y);

                    {
                        lock_guard<mutex> lock(rowlocks[seay]);
                        world.setriverdir(seax, seay, seadir);
                    }

                    seax = newseatile.x;
                    seay = newseatile.y;
                }
            }

            break;
        }

        if (y == 0 || y == height)
            break;

        const int nextkey = y * (width + 1) + x;

        if (visited.find(nextkey) != visited.end())
        {
            loopdetected = true;
            loopx = x;
            loopy = y;
            break;
        }
    }

    for (const pair<int, int>& cell : path)
        addriverflowlocked(world, cell.first, cell.second, janadd, juladd, rowlocks);

    if (loopdetected == true)
    {
        lock_guard<mutex> lock(rowlocks[loopy]);
        world.setlakesurface(loopx, loopy, maxelevation * 2);
    }
}

// This function forces rivers to avoid diagonals, as they don't look so good on the regional map.

void removediagonalrivers(planet& world)
{
    int width = world.width();
    int height = world.height();

    int adjustchance = 1; // Probability of adjusting these bits. The higher it is, the less likely.

    int desti, destj, adji, adjj, newdir1, newdir2, newheight;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 2; j <= height - 2; j++)
        {
            if (random(1, adjustchance) == 1)
            {
                bool candoit = 1;

                // Going northeast

                if (world.riverdir(i, j) == 2)
                {
                    desti = i + 1;
                    destj = j - 1;

                    if (desti > width)
                        desti = 0;

                    if (random(0, 1) == 1)
                    {
                        adji = i;
                        adjj = j - 1;

                        newdir1 = 1;
                        newdir2 = 3;

                        if (world.riverdir(adji, adjj) == 5 || world.riverdir(desti, destj) == 7)
                            candoit = 0;
                    }
                    else
                    {
                        adji = i + 1;
                        adjj = j;

                        if (adji > width)
                            adji = 0;

                        newdir1 = 3;
                        newdir2 = 1;

                        if (world.riverdir(adji, adjj) == 7 || world.riverdir(desti, destj) == 5)
                            candoit = 0;

                    }

                    newheight = (world.nom(i, j) + world.nom(desti, destj)) / 2;

                    if (newheight > world.nom(adji, adjj))
                        candoit = 0;

                    if (candoit == 1)
                    {
                        world.setnom(adji, adjj, newheight);
                        world.setriverdir(i, j, newdir1);
                        world.setriverdir(adji, adjj, newdir2);
                    }
                }

                // Going northwest

                if (world.riverdir(i, j) == 8)
                {
                    desti = i - 1;
                    destj = j - 1;

                    if (desti < 0)
                        desti = width;

                    if (random(0, 1) == 1)
                    {
                        adji = i;
                        adjj = j - 1;

                        newdir1 = 1;
                        newdir2 = 7;

                        if (world.riverdir(adji, adjj) == 5 || world.riverdir(desti, destj) == 3)
                            candoit = 0;
                    }
                    else
                    {
                        adji = i - 1;
                        adjj = j;

                        if (adji < 0)
                            adji = width;

                        newdir1 = 7;
                        newdir2 = 1;

                        if (world.riverdir(adji, adjj) == 5 || world.riverdir(desti, destj) == 3)
                            candoit = 0;

                    }

                    newheight = (world.nom(i, j) + world.nom(desti, destj)) / 2;

                    if (newheight > world.nom(adji, adjj))
                        candoit = 0;

                    if (candoit == 1)
                    {
                        world.setnom(adji, adjj, newheight);
                        world.setriverdir(i, j, newdir1);
                        world.setriverdir(adji, adjj, newdir2);
                    }
                }

                // Going southeast

                if (world.riverdir(i, j) == 4)
                {
                    desti = i + 1;
                    destj = j + 1;

                    if (desti > width)
                        desti = 0;

                    if (random(0, 1) == 1)
                    {
                        adji = i;
                        adjj = j + 1;

                        newdir1 = 5;
                        newdir2 = 3;

                        if (world.riverdir(adji, adjj) == 1 || world.riverdir(desti, destj) == 7)
                            candoit = 0;
                    }
                    else
                    {
                        adji = i + 1;
                        adjj = j;

                        if (adji > width)
                            adji = 0;

                        newdir1 = 3;
                        newdir2 = 5;

                        if (world.riverdir(adji, adjj) == 7 || world.riverdir(desti, destj) == 1)
                            candoit = 0;

                    }

                    newheight = (world.nom(i, j) + world.nom(desti, destj)) / 2;

                    if (newheight > world.nom(adji, adjj))
                        candoit = 0;

                    if (candoit == 1)
                    {
                        world.setnom(adji, adjj, newheight);
                        world.setriverdir(i, j, newdir1);
                        world.setriverdir(adji, adjj, newdir2);
                    }
                }

                // Going southwest

                if (world.riverdir(i, j) == 6)
                {
                    desti = i - 1;
                    destj = j + 1;

                    if (desti < 0)
                        desti = width;

                    if (random(0, 1) == 1)
                    {
                        adji = i;
                        adjj = j + 1;

                        newdir1 = 5;
                        newdir2 = 7;

                        if (world.riverdir(adji, adjj) == 1 || world.riverdir(desti, destj) == 3)
                            candoit = 0;
                    }
                    else
                    {
                        adji = i - 1;
                        adjj = j;

                        if (adji < 0)
                            adji = width;

                        newdir1 = 7;
                        newdir2 = 5;

                        if (world.riverdir(adji, adjj) == 3 || world.riverdir(desti, destj) == 1)
                            candoit = 0;

                    }

                    newheight = (world.nom(i, j) + world.nom(desti, destj)) / 2;

                    if (newheight > world.nom(adji, adjj))
                        candoit = 0;

                    if (candoit == 1)
                    {
                        world.setnom(adji, adjj, newheight);
                        world.setriverdir(i, j, newdir1);
                        world.setriverdir(adji, adjj, newdir2);
                    }
                }
            }
        }
    }
}

// This function makes rivers meet each other at diagonals, which looks a little more realistic.

void adddiagonalriverjunctions(planet& world)
{
    int width = world.width();
    int height = world.height();

    // First the N/S and E/W ones

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            // Going east

            if (world.riverdir(i, j) == 3 && world.riverdir(i, j + 1) == 1)
                world.setriverdir(i, j + 1, 2);

            if (world.riverdir(i, j + 1) == 3 && world.riverdir(i, j) == 5)
                world.setriverdir(i, j, 4);

            // Going west

            if (world.riverdir(i + 1, j) == 7 && world.riverdir(i + 1, j + 1) == 1)
                world.setriverdir(i + 1, j + 1, 8);

            if (world.riverdir(i + 1, j + 1) == 7 && world.riverdir(i + 1, j) == 5)
                world.setriverdir(i + 1, j, 6);

            // Going north

            if (world.riverdir(i, j + 1) == 1 && world.riverdir(i + 1, j + 1) == 7)
                world.setriverdir(i + 1, j + 1, 8);

            if (world.riverdir(i + 1, j + 1) == 1 && world.riverdir(i, j + 1) == 3)
                world.setriverdir(i, j + 1, 2);

            // Going south

            if (world.riverdir(i, j) == 5 && world.riverdir(i + 1, j) == 7)
                world.setriverdir(i + 1, j, 6);

            if (world.riverdir(i + 1, j) == 5 && world.riverdir(i, j) == 3)
                world.setriverdir(i, j, 4);
        }
    }

    // Now the diagonal ones.

    for (int i = 0; i < width - 1; i++)
    {
        for (int j = 0; j < height - 1; j++)
        {
            // Going northeast

            if (world.riverdir(i, j) == 4 && world.riverdir(i + 1, j + 1) == 2)
            {
                if (world.nom(i + 1, j) > world.nom(i + 2, j))
                {
                    world.setriverdir(i, j, 3);
                    world.setriverdir(i + 1, j, 3);

                    if (world.nom(i + 1, j) > world.nom(i, j))
                        world.setnom(i + 1, j, world.nom(i, j) - 1);
                }
            }

            if (world.riverdir(i, j + 1) == 2 && world.riverdir(i + 1, j + 2) == 8)
            {
                if (world.nom(i + 1, j + 1) > world.nom(i + 1, j))
                {
                    world.setriverdir(i + 1, j + 2, 1);
                    world.setriverdir(i + 1, j + 1, 1);

                    if (world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2))
                        world.setnom(i + 1, j + 1, world.nom(i + 1, j + 1) - 1);
                }
            }

            // Going southwest

            if (world.riverdir(i, j) == 4 && world.riverdir(i + 1, j + 1) == 6)
            {
                if (world.nom(i, j + 1) > world.nom(i, j + 1))
                {
                    world.setriverdir(i, j, 5);
                    world.setriverdir(i, j + 1, 5);

                    if (world.nom(i, j + 1) > world.nom(i, j))
                        world.setnom(i, j + 1, world.nom(i, j) - 1);
                }
            }

            if (world.riverdir(i + 1, j) == 6 && world.riverdir(i + 2, j + 1) == 8)
            {
                if (world.nom(i + 1, j + 1) > world.nom(i, j + 1))
                {
                    world.setriverdir(i + 2, j + 1, 7);
                    world.setriverdir(i + 1, j + 1, 7);

                    if (world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1))
                        world.setnom(i + 1, j + 1, world.nom(i + 1, j + 1) - 1);
                }
            }

            // Going southeast

            if (world.riverdir(i + 1, j) == 4 && world.riverdir(i, j + 1) == 2)
            {
                if (world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1))
                {
                    world.setriverdir(i, j + 1, 3);
                    world.setriverdir(i + 1, j + 1, 3);

                    if (world.nom(i + 1, j + 1) > world.nom(i, j + 1))
                        world.setnom(i + 1, j + 1, world.nom(i, j + 1) - 1);
                }
            }

            if (world.riverdir(i, j + 1) == 4 && world.riverdir(i + 1, j) == 6)
            {
                if (world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2))
                {
                    world.setriverdir(i + 1, j, 5);
                    world.setriverdir(i + 1, j + 1, 5);

                    if (world.nom(i + 1, j + 1) > world.nom(i + 1, j))
                        world.setnom(i + 1, j + 1, world.nom(i + 1, j) - 1);
                }
            }

            // Going northwest

            if (world.riverdir(i + 1, j + 1) == 8 && world.riverdir(i + 2, j) == 6)
            {
                if (world.nom(i + 1, j) > world.nom(i, j))
                {
                    world.setriverdir(i + 2, j, 7);
                    world.setriverdir(i + 1, j, 7);

                    if (world.nom(i + 1, j) > world.nom(i + 2, j))
                        world.setnom(i + 1, j, world.nom(i + 2, j) - 1);
                }
            }

            if (world.riverdir(i + 1, j + 1) == 8 && world.riverdir(i, j + 2) == 2)
            {
                if (world.nom(i, j + 1) > world.nom(i, j))
                {
                    world.setriverdir(i, j + 2, 1);
                    world.setriverdir(i, j + 1, 1);

                    if (world.nom(i, j + 1) > world.nom(i, j + 2))
                        world.setnom(i, j + 1, world.nom(i, j + 2) - 1);
                }
            }
        }
    }
}

// This function makes rivers avoid volcanoes where possible.

void avoidvolcanoes(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    vector<vector<short>> newdirs(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, -1));

    parallelforrows(0, height - 1, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i < width; i++)
            {
                int dir = world.riverdir(i, j);

                if (dir == 0)
                    continue;

                twointegers dest = getdestination(i, j, dir);

                if (world.volcano(dest.x, dest.y) != 0)
                {
                    int currentelev = world.nom(i, j);

                    int lowest = maxelev;
                    int newx = dest.x;
                    int newy = dest.y;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk < 0 || kk > width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height && (kk != i || l != j) && world.volcano(kk, l) == 0)
                            {
                                int nom = world.nom(kk, l);

                                if (nom < currentelev && nom < lowest)
                                {
                                    lowest = nom;
                                    newx = kk;
                                    newy = l;
                                }
                            }
                        }
                    }

                    newdirs[i][j] = getdir(i, j, newx, newy);
                }
            }
        }
    });

    parallelforrows(0, height - 1, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i < width; i++)
            {
                if (newdirs[i][j] != -1)
                    world.setriverdir(i, j, newdirs[i][j]);
            }
        }
    });
}

// This function makes parallel rivers run into each other.

void removeparallelrivers(planet& world)
{
    int width = world.width();
    int height = world.height();

    int removechance = 1;

    vector<vector<int>> altered(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This will show any points that have been involved in a change

    // First, N/S and E/W

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            if (random(1, removechance) == 1)
            {
                // Going east

                if (world.riverdir(i, j) == 3 && world.riverdir(i + 1, j) == 3 && world.riverdir(i, j + 1) == 3 && world.riverdir(i + 1, j + 1) == 3)
                {
                    if (world.nom(i, j) > world.nom(i + 1, j + 1) && altered[i][j] == 0)
                    {
                        world.setriverdir(i, j, 4);

                        altered[i][j] = 1;
                        altered[i + 1][j] = 1;
                        altered[i][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;
                    }
                    else
                    {
                        if (world.nom(i, j + 1) > world.nom(i + 1, j) && altered[i][j + 1] == 0)
                        {
                            world.setriverdir(i, j + 1, 2);

                            altered[i][j] = 1;
                            altered[i + 1][j] = 1;
                            altered[i][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;
                        }
                    }
                }

                // Going west

                if (world.riverdir(i, j) == 7 && world.riverdir(i + 1, j) == 7 && world.riverdir(i, j + 1) == 7 && world.riverdir(i + 1, j + 1) == 7)
                {
                    if (world.nom(i + 1, j) > world.nom(i, j + 1) && altered[i + 1][j] == 0)
                    {
                        world.setriverdir(i + 1, j, 6);

                        altered[i][j] = 1;
                        altered[i + 1][j] = 1;
                        altered[i][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;
                    }
                    else
                    {
                        if (world.nom(i + 1, j + 1) > world.nom(i, j) && altered[i + 1][j + 1] == 0)
                        {
                            world.setriverdir(i + 1, j + 1, 8);

                            altered[i][j] = 1;
                            altered[i + 1][j] = 1;
                            altered[i][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;
                        }
                    }
                }

                // Going north

                if (world.riverdir(i, j) == 1 && world.riverdir(i, j + 1) == 1 && world.riverdir(i + 1, j) == 1 && world.riverdir(i + 1, j + 1) == 7)
                {
                    if (world.nom(i, j + 1) > world.nom(i + 1, j) && altered[i][j + 1] == 0)
                    {
                        world.setriverdir(i, j + 1, 2);

                        altered[i][j] = 1;
                        altered[i + 1][j] = 1;
                        altered[i][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;
                    }
                    else
                    {
                        if (world.nom(i + 1, j + 1) > world.nom(i, j) && altered[i + 1][j + 1] == 0)
                        {
                            world.setriverdir(i + 1, j + 1, 8);

                            altered[i][j] = 1;
                            altered[i + 1][j] = 1;
                            altered[i][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;
                        }
                    }
                }

                // Going south

                if (world.riverdir(i, j) == 5 && world.riverdir(i, j + 1) == 5 && world.riverdir(i + 1, j) == 5 && world.riverdir(i + 1, j + 1) == 5)
                {
                    if (world.nom(i, j) > world.nom(i + 1, j + 1) && altered[i][j] == 0)
                    {
                        world.setriverdir(i, j, 4);

                        altered[i][j] = 1;
                        altered[i + 1][j] = 1;
                        altered[i][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;
                    }
                    else
                    {
                        if (world.nom(i + 1, j) > world.nom(i, j + 1) && altered[i + 1][j] == 0)
                        {
                            world.setriverdir(i + 1, j, 6);

                            altered[i][j] = 1;
                            altered[i + 1][j] = 1;
                            altered[i][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;
                        }
                    }
                }
            }
        }
    }

    // Now the diagonals.

    for (int i = 0; i <= width - 2; i++)
    {
        for (int j = 0; j <= height - 2; j++)
        {
            if (random(1, removechance) == 1)
            {
                // Going northeast

                if (world.riverdir(i, j + 1) == 2 && world.riverdir(i + 1, j) == 2 && world.riverdir(i + 1, j + 2) == 2 && world.riverdir(i + 2, j + 1) == 2)
                {
                    if (world.nom(i, j + 1) > world.nom(i + 2, j + 1) && world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1) && altered[i][j + 1] == 0)
                    {
                        world.setriverdir(i, j + 1, 3);
                        world.setriverdir(i + 1, j + 1, 3);

                        altered[i][j + 1] = 1;
                        altered[i + 1][j] = 1;
                        altered[i + 2][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;

                        if (world.nom(i + 1, j + 1) > world.nom(i, j + 1))
                            world.setnom(i + 1, j + 1, world.nom(i, j + 1) - 1);

                    }
                    else
                    {
                        if (world.nom(i + 1, j + 2) > world.nom(i + 1, j) && world.nom(i + 1, j + 1) > world.nom(i + 1, j) && altered[i + 1][j + 2] == 0)
                        {
                            world.setriverdir(i + 1, j + 2, 1);
                            world.setriverdir(i + 1, j + 1, 1);

                            altered[i][j + 1] = 1;
                            altered[i + 1][j] = 1;
                            altered[i + 2][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;

                            if (world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2))
                                world.setnom(i + 1, j + 1, world.nom(i + 1, j + 2) - 1);
                        }
                    }
                }

                // Going southwest

                if (world.riverdir(i, j + 1) == 6 && world.riverdir(i + 1, j) == 6 && world.riverdir(i + 1, j + 2) == 6 && world.riverdir(i + 2, j + 1) == 6)
                {
                    if (world.nom(i + 1, j) > world.nom(i + 1, j + 2) && world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2) && altered[i + 1][j] == 0)
                    {
                        world.setriverdir(i + 1, j, 5);
                        world.setriverdir(i + 1, j + 1, 5);

                        altered[i][j + 1] = 1;
                        altered[i + 1][j] = 1;
                        altered[i + 2][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;

                        if (world.nom(i + 1, j + 1) > world.nom(i + 1, j))
                            world.setnom(i + 1, j + 1, world.nom(i + 1, j) - 1);

                    }
                    else
                    {
                        if (world.nom(i + 2, j + 1) > world.nom(i, j + 1) && world.nom(i + 1, j + 1) > world.nom(i, j + 1) && altered[i + 2][j + 1] == 0)
                        {
                            world.setriverdir(i + 2, j + 1, 7);
                            world.setriverdir(i + 1, j + 1, 7);

                            altered[i][j + 1] = 1;
                            altered[i + 1][j] = 1;
                            altered[i + 2][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;

                            if (world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1))
                                world.setnom(i + 1, j + 1, world.nom(i + 2, j + 1) - 1);
                        }
                    }
                }

                // Going southeast

                if (world.riverdir(i, j + 1) == 4 && world.riverdir(i + 1, j) == 4 && world.riverdir(i + 1, j + 2) == 4 && world.riverdir(i + 2, j + 1) == 4)
                {
                    if (world.nom(i + 1, j) > world.nom(i + 1, j + 2) && world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2) && altered[i + 1][j] == 0)
                    {
                        world.setriverdir(i + 1, j, 5);
                        world.setriverdir(i + 1, j + 1, 5);

                        altered[i][j + 1] = 1;
                        altered[i + 1][j] = 1;
                        altered[i + 2][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;

                        if (world.nom(i + 1, j + 1) > world.nom(i + 1, j))
                            world.setnom(i + 1, j + 1, world.nom(i + 1, j) - 1);

                    }
                    else
                    {
                        if (world.nom(i, j + 1) > world.nom(i + 2, j + 1) && world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1) && altered[i][j + 1] == 0)
                        {
                            world.setriverdir(i, j + 1, 3);
                            world.setriverdir(i + 1, j + 1, 3);

                            altered[i][j + 1] = 1;
                            altered[i + 1][j] = 1;
                            altered[i + 2][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;

                            if (world.nom(i + 1, j + 1) > world.nom(i, j + 1))
                                world.setnom(i + 1, j + 1, world.nom(i, j + 1) - 1);
                        }
                    }
                }

                // Going northwest

                if (world.riverdir(i, j + 1) == 8 && world.riverdir(i + 1, j) == 8 && world.riverdir(i + 1, j + 2) == 8 && world.riverdir(i + 2, j + 1) == 8)
                {
                    if (world.nom(i + 2, j + 1) > world.nom(i, j + 1) && world.nom(i + 1, j + 1) > world.nom(i, j + 1) && altered[i + 2][j + 1] == 0)
                    {
                        world.setriverdir(i + 2, j + 1, 7);
                        world.setriverdir(i + 1, j + 1, 7);

                        altered[i][j + 1] = 1;
                        altered[i + 1][j] = 1;
                        altered[i + 2][j + 1] = 1;
                        altered[i + 1][j + 1] = 1;

                        if (world.nom(i + 1, j + 1) > world.nom(i + 2, j + 1))
                            world.setnom(i + 1, j + 1, world.nom(i + 2, j + 1) - 1);

                    }
                    else
                    {
                        if (world.nom(i + 1, j + 2) > world.nom(i + 1, j) && world.nom(i + 1, j + 1) > world.nom(i + 1, j) && altered[i + 1][j + 2] == 0)
                        {
                            world.setriverdir(i + 1, j + 2, 1);
                            world.setriverdir(i + 1, j + 1, 1);

                            altered[i][j + 1] = 1;
                            altered[i + 1][j] = 1;
                            altered[i + 2][j + 1] = 1;
                            altered[i + 1][j + 1] = 1;

                            if (world.nom(i + 1, j + 1) > world.nom(i + 1, j + 2))
                                world.setnom(i + 1, j + 1, world.nom(i + 1, j + 2) - 1);
                        }
                    }
                }
            }
        }
    }
}

// This function removes any rivers that cross each other.

void removecrossingrivers(planet& world)
{
    int width = world.width();
    int height = world.height();

    vector<vector<short>> newdirs(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, -1));

    parallelforrows(0, height - 1, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                int currentdir = world.riverdir(i, j);

                if (currentdir == 2)
                {
                    int ii = i;
                    int jj = j - 1;

                    if (jj >= 0 && jj <= height && world.riverdir(ii, jj) == 4 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 1;

                    ii = i + 1;
                    jj = j;

                    if (ii > width)
                        ii = 0;

                    if (world.riverdir(ii, jj) == 8 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 3;
                }

                if (currentdir == 4)
                {
                    int ii = i;
                    int jj = j + 1;

                    if (jj >= 0 && jj <= height && world.riverdir(ii, jj) == 2 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 5;

                    ii = i + 1;
                    jj = j;

                    if (ii > width)
                        ii = 0;

                    if (world.riverdir(ii, jj) == 6 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 3;
                }

                if (currentdir == 6)
                {
                    int ii = i;
                    int jj = j + 1;

                    if (jj >= 0 && jj <= height && world.riverdir(ii, jj) == 8 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 5;

                    ii = i - 1;
                    jj = j;

                    if (ii < 0)
                        ii = width;

                    if (world.riverdir(ii, jj) == 4 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 7;
                }

                if (currentdir == 8)
                {
                    int ii = i;
                    int jj = j - 1;

                    if (jj >= 0 && jj <= height && world.riverdir(ii, jj) == 6 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 1;

                    ii = i - 1;
                    jj = j;

                    if (ii < 0)
                        ii = width;

                    if (world.riverdir(ii, jj) == 2 && world.nom(ii, jj) <= world.nom(i, j))
                        newdirs[i][j] = 7;
                }
            }
        }
    });

    parallelforrows(0, height - 1, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (newdirs[i][j] != -1)
                    world.setriverdir(i, j, newdirs[i][j]);
            }
        }
    });
}

void planriverdirections(planet& world, vector<vector<int>>& mountaindrainage, std::uint64_t riverseed)
{
    const int width = world.width();
    const int height = world.height();
    int neighbours[8][2];

    neighbours[0][0] = 0;
    neighbours[0][1] = -1;

    neighbours[1][0] = 1;
    neighbours[1][1] = -1;

    neighbours[2][0] = 1;
    neighbours[2][1] = 0;

    neighbours[3][0] = 1;
    neighbours[3][1] = 1;

    neighbours[4][0] = 0;
    neighbours[4][1] = 1;

    neighbours[5][0] = -1;
    neighbours[5][1] = 1;

    neighbours[6][0] = -1;
    neighbours[6][1] = -0;

    neighbours[7][0] = -1;
    neighbours[7][1] = -1;

    // First, we go through the map and mark the direction of water flow on every land tile.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (j == 0)
                    world.setriverdir(i, j, 1);
                else if (j == height)
                    world.setriverdir(i, j, 5);
                else if (world.sea(i, j) == 0)
                {
                    int dir = findlowestdirriver(world, neighbours, i, j, mountaindrainage);

                    if (dir == -1)
                        dir = deterministicrandom(riverseed ^ 0x2001ull, 1, 8, i, j);

                    world.setriverdir(i, j, dir);
                }
                else
                    world.setriverdir(i, j, 0);
            }
        }
    });

}

void accumulaterivercatchments(planet& world, int minimum, int mountainheightlimit)
{
    const int width = world.width();
    const int height = world.height();
    vector<mutex> riverrowlocks(height + 1);
    // Now, we go through the map tile by tile. Take the rainfall in each tile and add it to every downstream tile.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 0 && world.strato(i, j) == 0)
                {
                    bool goahead = 1;

                    if (world.mountainheight(i, j) > mountainheightlimit)
                    {
                        int dir = world.riverdir(i, j);
                        int x = i;
                        int y = j;

                        if (dir == 8 || dir == 1 || dir == 2)
                            y--;

                        if (dir == 4 || dir == 5 || dir == 6)
                            y++;

                        if (dir == 2 || dir == 3 || dir == 4)
                            x++;

                        if (dir == 6 || dir == 7 || dir == 8)
                            x--;

                        if (x < 0 || x > width)
                            x = wrap(x, width);

                        if (y < 0)
                            y = 0;

                        if (y > height)
                            y = height;

                        if (world.mountainheight(x, y) > mountainheightlimit)
                            goahead = 0;
                    }

                    if (goahead == 1)
                        tracedropparallel(world, i, j, minimum, riverrowlocks);
                }
            }
        }
    });

}

void reduceseasonalrivercontrast(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const float tilt = world.tilt();
    if (tilt < 10.f) // For worlds with low axial tilt, reduce any seasonal difference in river flow.
    {
        int adjustfactor = (int)tilt;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    int aveflow = (world.riverjan(i, j) + world.riverjul(i, j)) / 2;

                    if (aveflow > 0)
                    {
                        int thisjanflow = world.riverjan(i, j) * adjustfactor + aveflow * (10 - adjustfactor);
                        int thisjulflow = world.riverjul(i, j) * adjustfactor + aveflow * (10 - adjustfactor);

                        world.setriverjan(i, j, thisjanflow);
                        world.setriverjul(i, j, thisjulflow);
                    }
                }
            }
        });
    }
}
}

// This function creates the rivers.

void createrivermap(planet& world, vector<vector<int>>& mountaindrainage)
{
    const std::uint64_t riverseed = deterministiccontextseed(world.seed(), 0x25b3f2a1);
    fast_srand(deterministicfastseed(riverseed));

    planriverdirections(world, mountaindrainage, riverseed);
    removediagonalrivers(world);
    adddiagonalriverjunctions(world);
    removeparallelrivers(world);
    avoidvolcanoes(world);
    removecrossingrivers(world);

    accumulaterivercatchments(world, tuning::climate::rivers::minimumFlow,
        tuning::climate::rivers::mountainHeightLimit);
    reduceseasonalrivercontrast(world);
}

// This function removes mountains that have rivers running over them.

void removerivermountains(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    int amount = 1; // Distance around mountain rivers to remove mountains.
    int mountainremovechance = 16; // The higher this is, the fewer mountains around rivers will be removed.

    int minremoveheight = 100; // Only mountains higher than this will be removed.

    vector<vector<bool>> toremove(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));
    vector<vector<int>> pathseen(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    int pathmark = 1;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 0 && world.mountainheight(i, j) == 0)
            {
                twointegers dest = getflowdestination(world, i, j, 0);

                if (world.mountainheight(dest.x, dest.y) > minremoveheight && world.sea(dest.x, dest.y) == 0) // We've got flow from a non-mountain cell into a mountain cell!
                {
                    bool keepgoing = 1;

                    int x = i;
                    int y = j;
                    int thispathmark = pathmark++;

                    do
                    {
                        if (pathseen[x][y] == thispathmark)
                        {
                            world.setlakesurface(x, y, maxelev * 2);
                            keepgoing = 0;
                            break;
                        }

                        pathseen[x][y] = thispathmark;

                        for (int k = x - amount; k <= x + amount; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = y - amount; l <= y + amount; l++)
                            {
                                if (l >= 0 && l <= height && random(1, mountainremovechance) == 1)
                                    toremove[kk][l] = 1;
                            }
                        }

                        toremove[x][y] = 1;

                        dest = getflowdestination(world, x, y, 0);

                        x = dest.x;
                        y = dest.y;

                        if (world.riverdir(x, y) == 0)
                            keepgoing = 0;

                        if (y == 0 || y == height)
                            keepgoing = 0;

                        if (world.sea(x, y) == 1)
                            keepgoing = 0;

                        if (world.mountainheight(x, y) == 0)
                            keepgoing = 0;

                    } while (keepgoing == 1);
                }
            }
        }
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (toremove[i][j] == 1 && world.sea(i, j) == 0)
            {
                world.setmintemp(i, j, tempelevremove(world, world.mintemp(i, j), i, j));
                world.setmaxtemp(i, j, tempelevremove(world, world.maxtemp(i, j), i, j));

                int totalrain = world.winterrain(i, j) + world.summerrain(i, j);

                short crount = 0;

                int winterraintotal = 0;
                int summerraintotal = 0;

                for (int x = i - 1; x <= i + 1; x++)
                {
                    int xx = x;

                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    for (int y = j - 1; y <= j + 1; y++)
                    {
                        if (y >= 0 && y <= height)
                        {
                            if (xx != y || y != j)
                            {
                                if (world.mountainheight(xx, y) == 0 && toremove[xx][y] == 0)
                                {
                                    crount++;

                                    winterraintotal = winterraintotal + world.winterrain(xx, y);
                                    summerraintotal = summerraintotal + world.summerrain(xx, y);
                                }
                            }
                        }
                    }
                }

                if (crount > 0)
                {
                    world.setwinterrain(i, j, winterraintotal / crount);
                    world.setsummerrain(i, j, summerraintotal / crount);
                }

                world.setmountainridge(i, j, 0);
                world.setmountainheight(i, j, 0);
            }
        }
    }

    cleanmountainridges(world);
}

// This function checks to make sure that river flows don't decrease.

void checkglobalflows(planet& world)
{
    int width = world.width();
    int height = world.height();

    bool found;

    do
    {
        vector<vector<int>> nextjan(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
        vector<vector<int>> nextjul(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
        vector<mutex> rowlocks(height + 1);
        atomic<bool> anychanged(false);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    if (world.sea(i, j) == 0 && (world.riverjan(i, j) != 0 || world.riverjul(i, j) != 0))
                    {
                        twointegers destpoint = getflowdestination(world, i, j, 0);

                        if (destpoint.x != -1)
                        {
                            int dex = destpoint.x;
                            int dey = destpoint.y;

                            int jan = world.riverjan(i, j);
                            int jul = world.riverjul(i, j);
                            int jand = world.riverjan(dex, dey);
                            int juld = world.riverjul(dex, dey);

                            if (jand < jan || juld < jul)
                            {
                                lock_guard<mutex> lock(rowlocks[dey]);

                                if (jan > nextjan[dex][dey])
                                    nextjan[dex][dey] = jan;

                                if (jul > nextjul[dex][dey])
                                    nextjul[dex][dey] = jul;

                                anychanged.store(true, memory_order_relaxed);
                            }
                        }
                    }
                }
            }
        });

        found = anychanged.load(memory_order_relaxed);

        if (found)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int j = startrow; j <= endrow; j++)
                {
                    for (int i = 0; i <= width; i++)
                    {
                        if (nextjan[i][j] > world.riverjan(i, j))
                            world.setriverjan(i, j, nextjan[i][j]);

                        if (nextjul[i][j] > world.riverjul(i, j))
                            world.setriverjul(i, j, nextjul[i][j]);
                    }
                }
            });
        }
    } while (found == true);
}

// This marks the route of a river on an array, from the given point.

void markriver(planet& world, int x, int y, vector<vector<int>>& markedarray, int extra)
{
    int width = world.width();
    int height = world.height();

    if (y<0 || y>height || x<0 || x>width)
        return;

    bool keepgoing = 1;

    do
    {
        if (markedarray[x][y] == 1) // If we're somehow going round in a loop)
            return;

        markedarray[x][y] = 1;

        if (extra == 1)
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
                        if (markedarray[ii][j] == 0)
                            markedarray[ii][j] = 2;
                    }
                }
            }
        }

        int dir = world.riverdir(x, y);

        if (dir == 0)
            return;

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

        if (world.sea(x, y) == 1)
            keepgoing = 0;

        if (y == 0 || y == height)
            keepgoing = 0;

    } while (keepgoing == 1);
}

// This diverts a river from the given point to the destination point.

int divertriver(planet& world, int x, int y, int destx, int desty, vector<vector<int>>& removedrivers, int riverno)
{
    twointegers destpoint;

    while (1 == 1)
    {
        if (x == destx && y == desty)
            return (riverno);

        int dir = getdir(x, y, destx, desty);

        if (dir == 0)
            return (riverno);

        if (world.riverdir(x, y) == dir)
        {
            destpoint = getflowdestination(world, x, y, 0);
            x = destpoint.x;
            y = destpoint.y;
        }
        else
        {
            // Delete the existing river

            riverno++;

            int janflow = world.riverjan(x, y);
            int julflow = world.riverjul(x, y);

            removeriver(world, removedrivers, riverno, x, y);

            // Now reinstate the river, pointing in the new direction

            world.setriverdir(x, y, dir);
            world.setriverjan(x, y, janflow);
            world.setriverjul(x, y, julflow);

            // Now move to the next point in the new river

            destpoint = getflowdestination(world, x, y, 0);

            x = destpoint.x;
            y = destpoint.y;

            addtoriver(world, x, y, janflow, julflow);
        }
    }
    return (riverno);
}

// This diverts a river that's next to a lake, to go into it.

int divertlakeriver(planet& world, int x, int y, int destx, int desty, vector<vector<int>>& removedrivers, int riverno, int outflowx, int outflowy, vector<vector<int>>& avoidarray)
{
    twointegers destpoint;

    for (int n = 1; n < 1000; n++) // It very occasionally gets stuck in an endless loop, so limit the number of times it can go around.
    {
        if (x == destx && y == desty)
            return (riverno);

        if (avoidarray[x][y] == 1)
            return (riverno);

        // First, get the direction to the destination.

        int dir = getdir(x, y, destx, desty);

        if (dir == 0)
            return (riverno);

        // Now we have to check to see whether that goes onto dry land, which we don't want.

        destpoint = getflowdestination(world, x, y, dir);

        int newdir = 0;

        if (world.lakesurface(destpoint.x, destpoint.y) == 0)
        {
            // First, see if we can swerve right.

            newdir = dir + 1;

            if (newdir == 9)
                newdir = 1;

            destpoint = getflowdestination(world, x, y, newdir);

            if (world.lakesurface(destpoint.x, destpoint.y) == 0)
            {

                // now see if we can swerve left.

                newdir = dir - 1;

                if (newdir == 0)
                    newdir = 8;

                destpoint = getflowdestination(world, x, y, newdir);

                if (world.lakesurface(destpoint.x, destpoint.y) == 0)
                {

                    // If we can't, then we will have to stop.

                    int janflow = world.riverjan(x, y);
                    int julflow = world.riverjul(x, y);

                    world.setriverjan(destx, desty, world.riverjan(destx, desty) + janflow);
                    world.setriverjul(destx, desty, world.riverjul(destx, desty) + julflow);

                    return (riverno);
                }
            }

            dir = newdir;
        }

        if (world.riverdir(x, y) == dir)
        {
            destpoint = getflowdestination(world, x, y, dir);

            x = destpoint.x;
            y = destpoint.y;
        }
        else
        {
            if (avoidarray[x][y] == 1)
                return (riverno);

            // Delete the existing river

            riverno++;

            int janflow = world.riverjan(x, y);
            int julflow = world.riverjul(x, y);

            removeriver(world, removedrivers, riverno, x, y);

            // Now reinstate the river, pointing in the new direction

            world.setriverdir(x, y, dir);
            world.setriverjan(x, y, janflow);
            world.setriverjul(x, y, julflow);

            // Now move to the next point in the new river

            destpoint = getflowdestination(world, x, y, 0);

            x = destpoint.x;
            y = destpoint.y;

            if (avoidarray[x][y] != 1)
                addtoriver(world, x, y, janflow, julflow);

        }
    }
    return (riverno);
}

// This removes a river, starting from the given point.

void removeriver(planet& world, vector<vector<int>>& removedrivers, int riverno, int x, int y)
{
    int janload = world.riverjan(x, y);
    int julload = world.riverjul(x, y);

    if (janload == 0 && julload == 0) // No river here...
        return;

    int width = world.width();
    int height = world.height();

    bool keepgoing = 1;

    do
    {
        int dir = world.riverdir(x, y);

        removedrivers[x][y] = riverno;

        if (world.riverjan(x, y) - janload >= 0)
            world.setriverjan(x, y, world.riverjan(x, y) - janload);
        else
            return;

        if (world.riverjul(x, y) - julload >= 0)
            world.setriverjul(x, y, world.riverjul(x, y) - julload);
        else
            return;

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

        if (world.sea(x, y) == 1)
            keepgoing = 0;

        if (y == 0 || y == height)
            keepgoing = 0;

        if (removedrivers[x][y] == riverno)
            keepgoing = 0;

    } while (keepgoing == 1);

}

// This lowers and reduces a river, starting from the given point.

void reduceriver(planet& world, int janreduce, int julreduce, vector<vector<int>>& removedrivers, int riverno, int x, int y)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    bool keepgoing = 1;

    do
    {
        int dir = world.riverdir(x, y);

        removedrivers[x][y] = riverno;

        if (world.riverjan(x, y) - janreduce >= 0)
            world.setriverjan(x, y, world.riverjan(x, y) - janreduce);
        else
            return;

        if (world.riverjul(x, y) - julreduce >= 0)
            world.setriverjul(x, y, world.riverjul(x, y) - julreduce);
        else
            return;

        world.setnom(x, y, sealevel + 1);

        if (world.deltadir(x, y) == 0)
            world.setdeltadir(x, y, -1);

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

        if (world.sea(x, y) == 1)
            keepgoing = 0;

        if (y == 0 || y == height)
            keepgoing = 0;

        if (removedrivers[x][y] == riverno)
            keepgoing = 0;

    } while (keepgoing == 1);
}

// This adds water to an existing river, starting from the given point.

void addtoriver(planet& world, int x, int y, int janload, int julload)
{
    int width = world.width();
    int height = world.height();

    if (y<0 || y>height || x<0 || x>width)
        return;

    vector<vector<int>> thisdrop(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    bool keepgoing = 1;

    do
    {
        int dir = world.riverdir(x, y);

        thisdrop[x][y] = 1;

        world.setriverjan(x, y, world.riverjan(x, y) + janload);
        world.setriverjul(x, y, world.riverjul(x, y) + julload);

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

        if (world.sea(x, y) == 1)
            keepgoing = 0;

        if (thisdrop[x][y] == 1)
            keepgoing = 0;

    } while (keepgoing == 1);

}
