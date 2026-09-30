// Procedural ocean ridges, faults, and trenches.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "classes.hpp"
#include "functions.hpp"
#include "planet.hpp"
#include "generation_tuning.hpp"

using namespace std;

namespace
{
bool shelfedge(planet& world, std::vector<std::vector<bool>>& shelves, int x, int y);
void drawoceanridgeline(planet& world, int fromx, int fromy, int tox, int toy, std::vector<std::vector<int>>& array, int value);
void createoceanfault(planet& world, int midx, int midy, int mindist, int maxdist, std::vector<std::vector<int>>& ridgesmap, std::vector<std::vector<bool>>& checked, int masksize);
}

// This creates mid-ocean ridges.

void createoceanridges(planet& world, vector<vector<bool>>& shelves)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    const std::uint64_t ridgeseed = deterministiccontextseed(world.seed(), 0x0cea6e51);
    fast_srand(deterministicfastseed(ridgeseed));

    vector<vector<int>> nearestshelfdist(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> nearestshelfx(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
    vector<vector<int>> nearestshelfy(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
    vector<vector<bool>> edgepoints(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));
    vector<vector<bool>> boundaries(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));
    vector<vector<int>> grid(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> gridnumbers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> ridges(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> ridgesmap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> ridgedistances(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    // We need a grid of edge points - points where the coastal shelves meet ocean.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (shelves[i][j] == 1 && shelfedge(world, shelves, i, j) == 1)
                    edgepoints[i][j] = 1;
            }
        }
    });

    // Now we need to go over all the ocean points and work out their closest edge points.

    vector<twointegers> frontier;
    frontier.reserve((width + 1) * 4);

    auto addfrontierseed = [&](int x, int y, int sourcex, int sourcey)
    {
        if (nearestshelfdist[x][y] != 0)
            return;

        nearestshelfdist[x][y] = 1;
        nearestshelfx[x][y] = sourcex;
        nearestshelfy[x][y] = sourcey;

        twointegers point;
        point.x = x;
        point.y = y;
        frontier.push_back(point);
    };

    for (int i = 0; i <= width; i++) // Every cell that *is* an edge point is closest to itself.
    {
        for (int j = 0; j <= height; j++)
        {
            if (edgepoints[i][j] == 1)
                addfrontierseed(i, j, i, j);
        }
    }

    for (int i = 0; i <= width; i++) // The northern and southern edges of the map count as shelf edges for this purpose.
    {
        addfrontierseed(i, 0, -1, -1);
        addfrontierseed(i, 1, -1, -1);
        addfrontierseed(i, height, -1, -1);
        addfrontierseed(i, height - 1, -1, -1);
    }

    int sweep = 2; // Because the first sweep was just setting up the edge points themselves

    while (frontier.empty() == false) // Now fill out the zones using a true frontier expansion.
    {
        vector<twointegers> nextfrontier;
        nextfrontier.reserve(frontier.size() * 2);

        for (const twointegers& point : frontier)
        {
            const int i = point.x;
            const int j = point.y;

            for (int k = i - 1; k <= i + 1; k++)
            {
                int kk = k;

                if (kk < 0 || kk > width)
                    kk = wrap(kk, width);

                for (int l = j - 1; l <= j + 1; l++)
                {
                    if (l >= 0 && l <= height)
                    {
                        if (k == i || l == j || deterministicrandom(ridgeseed ^ 0x1001ull, 1, 2, i, j, kk, l, sweep) == 1) // Only sometimes do diagonals, otherwise the end result looks too angular.
                        {
                            if (shelves[kk][l] == 0 && nearestshelfdist[kk][l] == 0)
                            {
                                nearestshelfdist[kk][l] = sweep;
                                nearestshelfx[kk][l] = nearestshelfx[i][j];
                                nearestshelfy[kk][l] = nearestshelfy[i][j];

                                twointegers nextpoint;
                                nextpoint.x = kk;
                                nextpoint.y = l;
                                nextfrontier.push_back(nextpoint);
                            }
                        }
                    }
                }
            }
        }

        frontier.swap(nextfrontier);
        sweep++;
    }

    // Now we need to find where the zones meet.

    int maxdiff = tuning::terrain::oceanridges::boundaryMaxSourceDifference;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (nearestshelfdist[i][j] > 0)
                {
                    bool found = 0;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk < 0 || kk > width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (shelves[i][j] == 0 && shelves[kk][l] == 0)
                                {
                                    if (boundaries[kk][l] == 0 && (nearestshelfx[i][j] - nearestshelfx[kk][l] > maxdiff || nearestshelfx[kk][l] - nearestshelfx[i][j] > maxdiff || nearestshelfy[i][j] - nearestshelfy[kk][l] > maxdiff || nearestshelfy[kk][l] - nearestshelfy[i][j] > maxdiff))
                                    {
                                        boundaries[i][j] = 1;
                                        found = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    });

    // We've got the boundaries marked out. Now we need to make a grid.

    int gridsize = tuning::terrain::oceanridges::gridSize;
    int halfgrid = gridsize / 2 - 1;

    for (int i = 0; i <= width; i = i + gridsize)
    {
        for (int j = 0; j <= height; j = j + gridsize)
        {
            for (int k = i - halfgrid; k <= i + halfgrid; k++)
            {
                int kk = k;

                if (kk<0 || kk>width)
                    kk = wrap(kk, width);

                for (int l = j - halfgrid; l <= j + halfgrid; l++)
                {
                    if (l >= 0 && l <= height)
                    {
                        if (boundaries[kk][l] == 1)
                        {
                            grid[i][j] = nearestshelfdist[kk][l];
                            k = i + halfgrid;
                            l = j + halfgrid;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i <= width; i = i + gridsize) // This one is to get rid of extra points around the diagonals.
    {
        int righti = i + gridsize;

        if (righti > width)
            righti = righti - width;

        for (int j = 0; j <= height; j = j + gridsize)
        {
            if (grid[i][j] != 0)
            {
                int upj = j - gridsize;
                int downj = j + gridsize;

                if (downj <= height)
                {
                    if (grid[righti][j] != 0 && grid[i][downj] != 0)
                        grid[i][j] = 0;
                }

                if (upj >= 0)
                {
                    if (grid[righti][j] != 0 && grid[i][upj] != 0)
                        grid[i][j] = 0;
                }
            }
        }
    }

    // Now we have a grid of points. We need to assign each one a unique number.

    int gridno = 1;

    for (int i = 0; i <= width; i = i + gridsize)
    {
        for (int j = 0; j <= height; j = j + gridsize)
        {
            if (grid[i][j] != 0)
            {
                gridnumbers[i][j] = gridno;
                gridno++;
            }
        }
    }

    int totalpoints = gridno;

    // Now we store which points neighbour each one.

    vector<int> tonorth(totalpoints + 1);
    vector<int> tonortheast(totalpoints + 1);
    vector<int> toeast(totalpoints + 1);
    vector<int> tosoutheast(totalpoints + 1);
    vector<int> distance(totalpoints + 1);
    vector<int> pointx(totalpoints + 1);
    vector<int> pointy(totalpoints + 1);

    for (int n = 0; n <= totalpoints; n++)
    {
        tonorth[n] = 0;
        tonortheast[n] = 0;
        toeast[n] = 0;
        tosoutheast[n] = 0;
        distance[n] = 0;
    }

    int checkdist = 1;

    for (int i = 0; i <= width; i = i + gridsize)
    {
        int righti = i + gridsize;

        if (righti > width)
            righti = righti - width;

        for (int j = 0; j <= height; j = j + gridsize)
        {
            if (grid[i][j] != 0)
            {
                int thispoint = gridnumbers[i][j];

                distance[thispoint] = nearestshelfdist[i][j];
                pointx[thispoint] = i;
                pointy[thispoint] = j;

                int upj = j - gridsize;
                int downj = j + gridsize;

                if (upj >= 0)
                {
                    for (int k = i - halfgrid; k <= i + halfgrid; k++)
                    {
                        int kk = k;

                        if (kk < -0 || kk >= width)
                            kk = wrap(kk, width);

                        if (grid[kk][upj] != 0)
                        {
                            tonorth[thispoint] = gridnumbers[kk][upj];
                            k = i + halfgrid;
                        }
                    }

                    for (int k = righti - halfgrid; k <= righti + halfgrid; k++)
                    {
                        int kk = k;

                        if (kk < -0 || kk >= width)
                            kk = wrap(kk, width);

                        if (grid[kk][upj] != 0)
                        {
                            tonortheast[thispoint] = gridnumbers[kk][upj];
                            k = righti + halfgrid;
                        }
                    }
                }

                for (int k = righti - halfgrid; k <= righti + halfgrid; k++)
                {
                    int kk = k;

                    if (kk < -0 || kk >= width)
                        kk = wrap(kk, width);

                    if (grid[kk][j] != 0)
                    {
                        toeast[thispoint] = gridnumbers[kk][j];
                        k = righti + halfgrid;
                    }
                }


                if (downj <= height)
                {
                    for (int k = righti - halfgrid; k <= righti + halfgrid; k++)
                    {
                        int kk = k;

                        if (kk < -0 || kk >= width)
                            kk = wrap(kk, width);

                        if (grid[kk][downj] != 0)
                        {
                            tosoutheast[thispoint] = gridnumbers[kk][downj];
                            k = righti + halfgrid;
                        }
                    }
                }
            }
        }
    }

    // Now for each point, we know its location, the distance to the closest shelf, and also which points (if any) border it to the N, NE, E, and SE. Now we just offset it a bit.

    // First, we want a fractal to offset all points.

    int grain = tuning::terrain::oceanridges::pointShiftFractalGrain;
    float valuemod = tuning::terrain::oceanridges::pointShiftFractalValueMod;
    int v = deterministicrandom(ridgeseed ^ 0x1002ull, tuning::terrain::oceanridges::pointShiftFractalValueMod2Min, tuning::terrain::oceanridges::pointShiftFractalValueMod2Max, width, height);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    int maxshift = tuning::terrain::oceanridges::maxShift;

    int maxadditionalshift = tuning::terrain::oceanridges::maxAdditionalShift;
    int minadditionalshift = tuning::terrain::oceanridges::minAdditionalShift;

    float div = (float)maxelev / (float)(maxshift * 2);

    for (int thispoint = 1; thispoint <= totalpoints; thispoint++)
    {
        float xshift = (float)fractal[pointx[thispoint]][pointy[thispoint]];

        xshift = xshift / div;
        xshift = xshift - (float)maxshift;

        int yy = pointx[thispoint] + width;

        if (yy > width)
            yy = yy - width;

        float yshift = (float)fractal[yy][pointy[thispoint]];

        yshift = yshift / div;
        yshift = yshift - (float)maxshift;

        if (xshift < 0.0f-(float)maxshift)
            xshift = 0.0f-(float)maxshift;

        if (xshift > (float)maxshift)
            xshift = (float)maxshift;

        if (yshift < 0.0f-(float)maxshift)
            yshift = 0.0f-(float)maxshift;

        if (yshift > (float)maxshift)
            yshift = (float)maxshift;

        int xextrashift = deterministicsignedrandom(ridgeseed ^ 0x1003ull, minadditionalshift, maxadditionalshift, thispoint, pointx[thispoint], pointy[thispoint], 0);
        int yextrashift = deterministicsignedrandom(ridgeseed ^ 0x1004ull, minadditionalshift, maxadditionalshift, thispoint, pointx[thispoint], pointy[thispoint], 1);

        int newx = pointx[thispoint] + (int)xshift + xextrashift;

        if (newx<0 || newx>width)
            newx = wrap(newx, width);

        int newy = pointy[thispoint] + (int)yshift + yextrashift;

        if (newy < 0)
            newy = 0;

        if (newy > height)
            newy = height;

        pointx[thispoint] = newx;
        pointy[thispoint] = newy;
    }

    // Now join up the points.

    for (int thispoint = 1; thispoint <= totalpoints; thispoint++)
    {
        if (tonorth[thispoint] != 0)
            drawoceanridgeline(world, pointx[thispoint], pointy[thispoint], pointx[tonorth[thispoint]], pointy[tonorth[thispoint]], ridges, distance[thispoint]);

        if (tonortheast[thispoint] != 0)
            drawoceanridgeline(world, pointx[thispoint], pointy[thispoint], pointx[tonortheast[thispoint]], pointy[tonortheast[thispoint]], ridges, distance[thispoint]);

        if (toeast[thispoint] != 0)
            drawoceanridgeline(world, pointx[thispoint], pointy[thispoint], pointx[toeast[thispoint]], pointy[toeast[thispoint]], ridges, distance[thispoint]);

        if (tosoutheast[thispoint] != 0)
            drawoceanridgeline(world, pointx[thispoint], pointy[thispoint], pointx[tosoutheast[thispoint]], pointy[tosoutheast[thispoint]], ridges, distance[thispoint]);
    }

    // Now get rid of excess points around the diagonals.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (ridges[i][j] != 0)
            {
                int ii = i + 1;

                if (ii > width)
                    ii = 0;

                if (j < height)
                {
                    if (ridges[ii][j + 1] != 0)
                    {
                        ridges[i][j + 1] = 0;
                        ridges[ii][j] = 0;
                    }
                }

                if (j > 0)
                {
                    if (ridges[ii][j - 1] != 0)
                    {
                        ridges[i][j - 1] = 0;
                        ridges[ii][j] = 0;
                    }
                }
            }
        }
    }

    // Now make sure each ridge point is adjacent to no more than one other ridge point.

    bool found = 0;

    do
    {
        found = 0;

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                if (ridges[i][j] != 0)
                {
                    int neighbours = -1; // Because we'll add one for itself

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (ridges[kk][l] != 0)
                                    neighbours++;

                            }
                        }
                    }

                    if (neighbours > 2)
                    {
                        found = 1;

                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (ridges[kk][l] != 0)
                                    {
                                        ridges[kk][l] = 0;
                                        k = i + 1;
                                        l = j + 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    } while (found == 1);

    vector<twointegers> ridgecells;
    ridgecells.reserve((width + 1) * (height + 1) / 32);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (ridges[i][j] != 0)
            {
                twointegers point;
                point.x = i;
                point.y = j;
                ridgecells.push_back(point);
            }
        }
    }

    auto buildcircleoffsets = [](int maxradius)
    {
        vector<vector<twointegers>> offsets(maxradius + 1);

        for (int radius = 1; radius <= maxradius; radius++)
        {
            vector<twointegers>& current = offsets[radius];

            for (int k = -radius; k <= radius; k++)
            {
                for (int l = -radius; l <= radius; l++)
                {
                    if (k * k + l * l < radius * radius + radius)
                    {
                        twointegers offset;
                        offset.x = k;
                        offset.y = l;
                        current.push_back(offset);
                    }
                }
            }
        }

        return offsets;
    };

    // Now raise the land around the ridges.

    const int firstmaxradius = tuning::terrain::oceanridges::firstPassMaxRadius;
    int heightmult = tuning::terrain::oceanridges::firstPassHeightMultiplier;
    int maxvolcanoradius = tuning::terrain::oceanridges::maxVolcanoRadius;
    const vector<vector<twointegers>> firstcircleoffsets = buildcircleoffsets(firstmaxradius);
    const int widthplusone = width + 1;
    const size_t ridgecellcount = static_cast<size_t>(widthplusone) * static_cast<size_t>(height + 1);

    auto flatridgeindex = [&](int x, int y)
    {
        return static_cast<size_t>(y) * static_cast<size_t>(widthplusone) + static_cast<size_t>(x);
    };

    auto ridgeworkerstouse = [](int totalitems, int minitemsperworker)
    {
        unsigned int workerstouse = std::thread::hardware_concurrency();

        if (workerstouse == 0)
            workerstouse = 4;

        if (workerstouse <= 1 || totalitems <= minitemsperworker)
            return 1u;

        const int maxworkers = std::max(1, totalitems / minitemsperworker);
        workerstouse = std::min(workerstouse, static_cast<unsigned int>(maxworkers));

        if (workerstouse == 0)
            workerstouse = 1;

        return workerstouse;
    };

    float thisheightperthousand[firstmaxradius + 1]; // Create a lookup table for these values, for speed.

    thisheightperthousand[1] = 1.1037f;

    for (int n = 2; n <= firstmaxradius; n++)
        thisheightperthousand[n] = thisheightperthousand[n - 1] * thisheightperthousand[1];

    const unsigned int firstpaintworkers = ridgeworkerstouse(static_cast<int>(ridgecells.size()), 16);
    vector<vector<int>> localridgeheights(firstpaintworkers, vector<int>(ridgecellcount, 0));
    vector<vector<unsigned char>> localridgedists(firstpaintworkers, vector<unsigned char>(ridgecellcount, 0));
    vector<std::thread> firstworkers;
    firstworkers.reserve(firstpaintworkers > 0 ? firstpaintworkers - 1 : 0);

    int firstcurrent = 0;
    const int firstbasechunk = static_cast<int>(ridgecells.size()) / static_cast<int>(firstpaintworkers);
    const int firstremainder = static_cast<int>(ridgecells.size()) % static_cast<int>(firstpaintworkers);

    for (unsigned int worker = 0; worker < firstpaintworkers; worker++)
    {
        const int chunksize = firstbasechunk + (worker < static_cast<unsigned int>(firstremainder) ? 1 : 0);
        const int chunkstart = firstcurrent;
        const int chunkend = firstcurrent + chunksize - 1;
        firstcurrent = chunkend + 1;

        auto paintfirstpass = [&](unsigned int workerindex, int startindex, int endindex)
        {
            vector<int>& workerheights = localridgeheights[workerindex];
            vector<unsigned char>& workerdists = localridgedists[workerindex];

            for (int ridgeindex = startindex; ridgeindex <= endindex; ridgeindex++)
            {
                const twointegers& ridgepoint = ridgecells[ridgeindex];
                const int i = ridgepoint.x;
                const int j = ridgepoint.y;

                float riftheight = (float)(ridges[i][j] * heightmult);
                float riftheightdiv = riftheight / 1000.0f;

                int mult = 0;

                for (int radius = firstmaxradius; radius >= 1; radius--)
                {
                    mult++;
                    const int thisheight = (int)(riftheightdiv * thisheightperthousand[mult]);

                    for (const twointegers& offset : firstcircleoffsets[radius])
                    {
                        int kk = i + offset.x;

                        if (kk < 0 || kk > width)
                            kk = wrap(kk, width);

                        int ll = j + offset.y;

                        if (ll >= 0 && ll <= height)
                        {
                            const size_t flatindex = flatridgeindex(kk, ll);

                            if (workerheights[flatindex] < thisheight)
                            {
                                workerheights[flatindex] = thisheight;
                                workerdists[flatindex] = static_cast<unsigned char>(mult);
                            }
                        }
                    }
                }
            }
        };

        if (chunkstart > chunkend)
            continue;

        if (worker + 1 == firstpaintworkers)
            paintfirstpass(worker, chunkstart, chunkend);
        else
            firstworkers.emplace_back(paintfirstpass, worker, chunkstart, chunkend);
    }

    for (std::thread& worker : firstworkers)
        worker.join();

    vector<int> ridgeheightflat(ridgecellcount, 0);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                const size_t flatindex = flatridgeindex(i, j);
                int bestheight = 0;
                unsigned char bestdist = 0;

                for (unsigned int worker = 0; worker < firstpaintworkers; worker++)
                {
                    const int candidate = localridgeheights[worker][flatindex];

                    if (candidate > bestheight)
                    {
                        bestheight = candidate;
                        bestdist = localridgedists[worker][flatindex];
                    }
                }

                if (bestheight > 0)
                {
                    ridgesmap[i][j] = bestheight;
                    ridgedistances[i][j] = bestdist;
                    ridgeheightflat[flatindex] = bestheight;
                    world.setoceanridgeheights(i, j, bestheight);
                }
            }
        }
    });

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (ridgedistances[i][j] != 0 && ridgedistances[i][j] <= maxvolcanoradius)
                {
                    const int thisheight = ridgeheightflat[flatridgeindex(i, j)];

                    if (thisheight > 0)
                    {
                        unsigned int hash = static_cast<unsigned int>(world.seed());
                        hash ^= static_cast<unsigned int>(i) + 0x9e3779b9u + (hash << 6) + (hash >> 2);
                        hash ^= static_cast<unsigned int>(j) + 0x9e3779b9u + (hash << 6) + (hash >> 2);
                        hash ^= static_cast<unsigned int>(thisheight) + 0x9e3779b9u + (hash << 6) + (hash >> 2);

                        if ((hash % 4000000u) < static_cast<unsigned int>(std::max(1, thisheight / 10)))
                        {
                            const int minvolcano = std::max(1, thisheight / 4);
                            const int maxvolcano = std::max(minvolcano, (thisheight / 4) * 3);
                            const unsigned int span = static_cast<unsigned int>(maxvolcano - minvolcano + 1);

                            hash ^= 0x85ebca6bu + (hash << 6) + (hash >> 2);
                            world.setvolcano(i, j, minvolcano + static_cast<int>(hash % span));
                        }
                    }
                }
            }
        }
    });

    int maxradius = tuning::terrain::oceanridges::secondPassMaxRadius;
    heightmult = tuning::terrain::oceanridges::secondPassHeightMultiplier;
    const vector<vector<twointegers>> secondcircleoffsets = buildcircleoffsets(maxradius);
    const unsigned int secondpaintworkers = ridgeworkerstouse(static_cast<int>(ridgecells.size()), 16);
    vector<vector<int>> localridgeadds(secondpaintworkers, vector<int>(ridgecellcount, 0));
    vector<std::thread> secondworkers;
    secondworkers.reserve(secondpaintworkers > 0 ? secondpaintworkers - 1 : 0);

    int secondcurrent = 0;
    const int secondbasechunk = static_cast<int>(ridgecells.size()) / static_cast<int>(secondpaintworkers);
    const int secondremainder = static_cast<int>(ridgecells.size()) % static_cast<int>(secondpaintworkers);

    for (unsigned int worker = 0; worker < secondpaintworkers; worker++)
    {
        const int chunksize = secondbasechunk + (worker < static_cast<unsigned int>(secondremainder) ? 1 : 0);
        const int chunkstart = secondcurrent;
        const int chunkend = secondcurrent + chunksize - 1;
        secondcurrent = chunkend + 1;

        auto paintsecondpass = [&](unsigned int workerindex, int startindex, int endindex)
        {
            vector<int>& workerheights = localridgeadds[workerindex];

            for (int ridgeindex = startindex; ridgeindex <= endindex; ridgeindex++)
            {
                const twointegers& ridgepoint = ridgecells[ridgeindex];
                const int i = ridgepoint.x;
                const int j = ridgepoint.y;
                const int ridgeheight = ridges[i][j] * heightmult;
                const int heightdivs = ridgeheight / maxradius;
                int mult = 0;

                for (int radius = maxradius; radius >= 1; radius--)
                {
                    mult++;

                    int thisheight = heightdivs * mult;
                    thisheight = (thisheight * 3 + ridgeheight / radius) / 4; // 3, 4

                    for (const twointegers& offset : secondcircleoffsets[radius])
                    {
                        int kk = i + offset.x;

                        if (kk < 0 || kk > width)
                            kk = wrap(kk, width);

                        int ll = j + offset.y;

                        if (ll >= 0 && ll <= height)
                        {
                            const size_t flatindex = flatridgeindex(kk, ll);

                            if (workerheights[flatindex] < thisheight)
                                workerheights[flatindex] = thisheight;
                        }
                    }
                }
            }
        };

        if (chunkstart > chunkend)
            continue;

        if (worker + 1 == secondpaintworkers)
            paintsecondpass(worker, chunkstart, chunkend);
        else
            secondworkers.emplace_back(paintsecondpass, worker, chunkstart, chunkend);
    }

    for (std::thread& worker : secondworkers)
        worker.join();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                const size_t flatindex = flatridgeindex(i, j);
                int bestheight = ridgesmap[i][j];

                for (unsigned int worker = 0; worker < secondpaintworkers; worker++)
                {
                    const int candidate = localridgeadds[worker][flatindex];

                    if (candidate > bestheight)
                        bestheight = candidate;
                }

                ridgesmap[i][j] = bestheight;
            }
        }
    });

    // Add the rift in the middle

    vector<twointegers> riftcells;
    riftcells.reserve(ridgecells.size());

    for (const twointegers& ridgepoint : ridgecells)
    {
        const int i = ridgepoint.x;
        const int j = ridgepoint.y;
        ridgesmap[i][j] = ridgesmap[i][j] / 2;
        world.setoceanrifts(i, j, ridgesmap[i][j]);

        if (ridgesmap[i][j] != 0)
            riftcells.push_back(ridgepoint);
    }

    // Now we have to work out the angles of lines crossing the rifts.

    // First, make another fractal.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            fractal[i][j] = 0;
    }

    grain = 4; // Level of detail on this fractal map.
    valuemod = 0.01f;
    v = 1; //random(3,6);
    valuemod2 = 0.1f;

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 360, 0, 1); // The values in this fractal wrap!

    // Now we need to go through every rift cell, and work out the angle of the line that would cross it at right angles.

    vector<vector<bool>> checked(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    short dist = tuning::terrain::oceanridges::ridgeAngleSearchDistance;

    int x1, y1, x2, y2;

    for (const twointegers& ridgepoint : riftcells)
    {
        const int i = ridgepoint.x;
        const int j = ridgepoint.y;

        // Clear this area in the checked array.

        for (int k = i - dist - 1; k <= i + dist + 1; k++)
        {
            int kk = k;

            if (kk < 0 || kk > width)
                kk = wrap(kk, width);

            for (int l = j - dist - 1; l <= j + dist + 1; l++)
            {
                if (l >= 0 && l <= height)
                    checked[kk][l] = 0;
            }
        }

        // Now go along the line twice - once in each direction.

        for (int dir = 1; dir <= 2; dir++)
        {
            int x = i;
            int y = j;

            for (int n = 1; n <= dist; n++)
            {
                checked[x][y] = 1;

                for (int k = x - 1; k <= x + 1; k++)
                {
                    int kk = k;
                    if (kk < 0 || kk > width)
                        kk = wrap(kk, width);

                    for (int l = y - 1; l <= y + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (checked[kk][l] == 0 && world.oceanrifts(kk, l) != 0)
                            {
                                x = kk;
                                y = l;

                                k = x + 20;
                                l = y + 20;
                            }
                        }
                    }
                }
            }

            if (dir == 1)
            {
                x1 = x;
                y1 = y;
            }
            else
            {
                x2 = x;
                y2 = y;
            }
        }

        // Now we have the two nearby points, we just find the angle between them.

        float angle = (float)atan2((float)y2 - (float)y1, (float)x2 - (float)x1) * 180.0f / 3.14159265358979323846f; // This gives us the angle of the ridge

        angle = angle + 180.0f; // Because we want the line that crosses the ridge

        while (angle > 360.0f)
            angle = angle - 360.0f;

        // Now we just have to blend that with the fractal angle for this point.

        int angle2 = wrappedaverage((int)angle, fractal[i][j], 360);

        while (angle2 > 180)
            angle2 = angle2 - 180;

        world.setoceanridgeangle(i, j, angle2);
    }

    // Mark out the ridge mountains, following the contours of the raised land

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (ridgedistances[i][j] != 0)
                {
                    int ii = i; // Looking north
                    int jj = j - 1;
                    int dir = 1;

                    if (jj >= 0)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }

                    ii = i + 1; // Looking northeast
                    jj = j - 1;
                    dir = 2;

                    if (ii > width)
                        ii = 0;

                    if (jj >= 0)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }

                    ii = i + 1; // Looking east
                    jj = j;
                    dir = 3;

                    if (ii > width)
                        ii = 0;

                    if (ridgedistances[ii][jj] == ridgedistances[i][j])
                    {
                        if (getoceanridge(world, i, j, dir) == 0)
                        {
                            int code = terrain::detail::getcode(dir);
                            world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                        }
                    }

                    ii = i + 1; // Looking southeast
                    jj = j + 1;
                    dir = 4;

                    if (ii > width)
                        ii = 0;

                    if (jj <= height)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }

                    ii = i; // Looking south
                    jj = j + 1;
                    dir = 5;

                    if (jj <= height)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }

                    ii = i - 1; // Looking southwest
                    jj = j + 1;
                    dir = 6;

                    if (ii < 0)
                        ii = width;

                    if (jj <= height)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }

                    ii = i - 1; // Looking west
                    jj = j;
                    dir = 7;

                    if (ii < 0)
                        ii = width;

                    if (ridgedistances[ii][jj] == ridgedistances[i][j])
                    {
                        if (getoceanridge(world, i, j, dir) == 0)
                        {
                            int code = terrain::detail::getcode(dir);
                            world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                        }
                    }

                    ii = i - 1; // Looking northwest
                    jj = j - 1;
                    dir = 8;

                    if (ii < 0)
                        ii = width;

                    if (jj >= 0)
                    {
                        if (ridgedistances[ii][jj] == ridgedistances[i][j])
                        {
                            if (getoceanridge(world, i, j, dir) == 0)
                            {
                                int code = terrain::detail::getcode(dir);
                                world.setoceanridges(i, j, world.oceanridges(i, j) + code);
                            }
                        }
                    }
                }
            }
        }
    });

    // Now remove extraneous ridges

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.oceanridges(i, j) != 0)
            {
                short total = 0;

                for (int dir = 1; dir <= 8; dir++)
                {
                    if (getoceanridge(world, i, j, dir) == 1)
                        total++;
                }

                if (total > 2)
                {
                    int extra = total - 2;

                    for (int n = 1; n <= extra; n++)
                    {
                        int a = 1;
                        int b = 8;
                        int c = 1;

                        if (deterministicrandom(ridgeseed ^ 0x1005ull, 1, 2, i, j, n, extra) == 1)
                        {
                            a = 8;
                            b = 1;
                            c = -1;
                        }

                        for (int dir = a; dir != b; dir = dir + c)
                        {
                            if (getoceanridge(world, i, j, dir) == 1)
                                terrain::detail::deleteoceanridge(world, i, j, dir);
                        }
                    }
                }
            }
        }
    }

    // Now we'll make a fractal. This will be used to displace the ridges in the regional map.

    int maxdisplace = tuning::terrain::oceanridges::regionalDisplacement;

    grain = tuning::terrain::oceanridges::regionalDisplacementFractalGrain;
    valuemod = tuning::terrain::oceanridges::regionalDisplacementValueMod;
    valuemod2 = tuning::terrain::oceanridges::regionalDisplacementValueMod2;

    vector<vector<int>> ridgefractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(ridgefractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    div = (float)maxelev / ((float)maxdisplace * 2.0f);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                float amount = (float)ridgefractal[i][j];

                amount = amount / div;
                amount = amount - (float)maxdisplace;

                if (amount < 0.0f - (float)maxdisplace)
                    amount = 0.0f - (float)maxdisplace;

                if (amount > (float)maxdisplace)
                    amount = (float)maxdisplace;

                world.setoceanridgeoffset(i, j, (int)amount);
            }
        }
    });

    // Now we need to disrupt the rifts by adding some faults.

    for (int n = 1; n <= tuning::terrain::oceanridges::faultPasses; n++)
    {
        int faultstep = tuning::terrain::oceanridges::faultStep;
        int faultvar = tuning::terrain::oceanridges::faultVariation;
        int lookdist = tuning::terrain::oceanridges::faultLookDistance;

        int faultstotal = 0;

        for (int i = 0; i < ARRAYWIDTH; i = i + faultstep)
        {
            for (int j = 0; j < ARRAYHEIGHT; j = j + faultstep)
            {
                int ii = i + deterministicsignedrandom(ridgeseed ^ 0x1006ull, 1, faultvar, n, i, j, 0);
                int jj = j + deterministicsignedrandom(ridgeseed ^ 0x1007ull, 1, faultvar, n, i, j, 1);

                if (ii<0 || ii>width)
                    ii = wrap(ii, width);

                if (jj < 0)
                    jj = 0;

                if (jj > height)
                    jj = height;

                if (world.sea(ii, jj) == 1)
                {
                    for (int x = ii - lookdist; x <= ii + lookdist; x++)
                    {
                        int xx = x;

                        if (xx<0 || xx>width)
                            xx = wrap(xx, width);

                        for (int y = jj - lookdist; y <= jj + lookdist; y++)
                        {
                            if (y >= 0 && y <= height)
                            {
                                if (world.oceanrifts(xx, y) != 0)
                                {
                                    int maxwidth = 30;
                                    int minwidth = 6;

                                    int masksize = 200 + maxwidth;

                                    createoceanfault(world, xx, y, minwidth, maxwidth, ridgesmap, checked, masksize);

                                    x = ii + lookdist;
                                    y = jj + lookdist;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Remove any ridge-related stuff that's on continental plates

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (shelves[i][j] == 1 || world.sea(i, j) == 0)
                {
                    ridgesmap[i][j] = 0;
                    world.setoceanrifts(i, j, 0);
                    world.setoceanridges(i, j, 0);
                    world.setoceanridgeheights(i, j, 0);
                    world.setoceanridgeangle(i, j, 0);
                    world.setoceanridgeoffset(i, j, 0);
                }
            }
        }
    });

    // Draw the raised land onto the actual sea bed

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
                world.setnom(i, j, world.nom(i, j) + ridgesmap[i][j]);
        }
    });
}

namespace
{
// This checks to see whether the given point is on the edge of a continental shelf.

bool shelfedge(planet& world, vector<vector<bool>>& shelves, int x, int y)
{
    if (shelves[x][y] == 0)
        return 0;

    int width = world.width();
    int height = world.height();

    for (int i = x - 1; i <= x + 1; i++)
    {
        int ii = i;

        if (ii < 0)
            ii = width;

        if (ii > width)
            ii = 0;

        for (int j = y - 1; j <= y + 1; j++)
        {
            if (j >= 0 && j <= height)
            {
                if (shelves[ii][j] == 0)
                    return 1;
            }
        }
    }
    return 0;
}
}

namespace
{
// This draws a curvy line between two points on the ocean ridge map.

void drawoceanridgeline(planet& world, int fromx, int fromy, int tox, int toy, vector<vector<int>>& array, int value)
{
    int width = world.width();
    int height = world.height();
    const std::uint64_t lineseed = deterministiccontextseed(world.seed(), 0x51a9e245);

    if (tox < fromx - width / 2)
        tox = tox + width;

    if (fromx < tox - width / 2)
        fromx = fromx + width;

    int distance = abs(fromx - tox) + abs(fromy - toy);

    int minvar = 2;
    int maxvar = distance / 8;

    if (maxvar < 3)
        maxvar = 3;

    twofloats pt, mm1, mm2, mm3;

    mm1.x = (float) fromx;
    mm1.y = (float) fromy;

    mm2.x = (float)(fromx + tox) / 2 + deterministicsignedrandom(lineseed ^ 0x1008ull, minvar, maxvar, fromx, fromy, tox, toy, value, 0);
    mm2.y = (float)(fromy + toy) / 2 + deterministicsignedrandom(lineseed ^ 0x1009ull, minvar, maxvar, fromx, fromy, tox, toy, value, 1);

    mm3.x = (float) tox;
    mm3.y = (float) toy;

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

            array[x][y] = value;
        }
    }

}
}

namespace
{
// This creates a fault in an oceanic rift.

void createoceanfault(planet& world, int midx, int midy, int mindist, int maxdist, vector<vector<int>>& ridgesmap, vector<vector<bool>>& checked, int masksize)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    const std::uint64_t faultseed = deterministiccontextseed(world.seed(), 0x7f4a7c15);

    int minshiftdist = 2;
    int maxshiftdist = 6;

    int dist = deterministicrandom(faultseed ^ 0x100aull, mindist, maxdist, midx, midy, masksize); // Distance from the midpoint that the edge points should be.

    int x1 = 0, y1 = 0, x2 = 0, y2 = 0; // Coordinates of the rift points at the ends of the area to be moved.

    // Clear this area in the checked array.

    for (int i = midx - 20; i <= midx + 20; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = midy - 20; j <= midy + 20; j++)
        {
            if (j >= 0 && j <= height)
                checked[ii][j] = 0;
        }
    }

    // Now go along the line twice - once in each direction.

    for (int dir = 1; dir <= 2; dir++)
    {
        int x = midx;
        int y = midy;

        for (int n = 1; n <= dist; n++)
        {
            checked[x][y] = 1;

            bool found = 0;

            for (int k = x - 1; k <= x + 1; k++)
            {
                int kk = k;
                if (kk<0 || kk>width)
                    kk = wrap(kk, width);

                for (int l = y - 1; l <= y + 1; l++)
                {
                    if (l >= 0 && l <= height)
                    {
                        if (checked[kk][l] == 0 && world.oceanrifts(kk, l) != 0)
                        {
                            x = kk;
                            y = l;

                            found = 1;
                            k = x + 20;
                            l = y + 20;
                        }
                    }
                }
            }

            if (found == 0) // This means we reached the end of the rift line, so we can't do it.
                return;
        }

        if (dir == 1)
        {
            x1 = x;
            y1 = y;
        }
        else
        {
            x2 = x;
            y2 = y;
        }
    }

    // Now we need to find the coordinates of the four corners of the area to be moved.

    int length = 70; //100; // Distance to go out from the central rift.

    int cornerx1 = 0, cornery1 = 0, cornerx2 = 0, cornery2 = 0, cornerx3 = 0, cornery3 = 0, cornerx4 = 0, cornery4 = 0;

    for (int n = 1; n <= 2; n++) // Do it for each of the two end points
    {
        int x, y;

        if (n == 1)
        {
            x = x1;
            y = y1;
        }
        else
        {
            x = x2;
            y = y2;
        }

        float angle = (float)(0 - world.oceanridgeangle(x, y));

        while (angle < 0.0f)
            angle = angle + 360.0f;

        for (int nn = 1; nn <= 2; nn++) // Two halves - one on each side of the rift.
        {
            float fangle;

            if (nn == 1)
                fangle = angle * 0.01745329f;
            else
            {
                fangle = angle + 180.0f;

                if (fangle > 360.0f)
                    fangle = fangle - 360.0f;

                fangle = fangle * 0.01745329f;
            }

            int endx = x + (int)((float)length * (float)sin(fangle));
            int endy = y + (int)((float)length * (float)cos(fangle));

            if (n == 1)
            {
                if (nn == 1)
                {
                    cornerx1 = endx;
                    cornery1 = endy;
                }
                else
                {
                    cornerx2 = endx;
                    cornery2 = endy;
                }
            }
            else
            {
                if (nn == 1)
                {
                    cornerx4 = endx;
                    cornery4 = endy;
                }
                else
                {
                    cornerx3 = endx;
                    cornery3 = endy;
                }
            }
        }
    }

    if (cornery1<0 || cornery1>height || cornery2<0 || cornery2>height || cornery3<0 || cornery3>height || cornery4<0 || cornery4>height)
        return;

    // Now make sure that all of our points are equally wrapped/unwrapped.

    normalise(midx, cornerx1, width);
    normalise(midx, cornerx2, width);
    normalise(midx, cornerx3, width);
    normalise(midx, cornerx4, width);

    // We now have the four corners of the whole area to be shifted.

    // Check that they don't cross over each other:

    if (x1 < x2)
    {
        if (cornerx1 > cornerx4)
            return;

        if (cornerx2 > cornerx3)
            return;
    }
    else
    {
        if (cornerx1 < cornerx4)
            return;

        if (cornerx2 < cornerx3)
            return;
    }

    if (y1 < y2)
    {
        if (cornery1 > cornery4)
            return;

        if (cornery2 > cornery3)
            return;
    }
    else
    {
        if (cornery1 < cornery4)
            return;

        if (cornery2 < cornery3)
            return;
    }

    // Now we need to create a mask of this area.

    vector<vector<bool>> mask(masksize + 1, vector<bool>(masksize + 1, 0));

    // To convert coordinates from the map to the mask, add the offsets.
    // To convert them from the mask to the map, subtract the offsets.

    int offsetx = masksize / 2 - midx;
    int offsety = masksize / 2 - midy;

    cornerx1 = cornerx1 + offsetx;
    cornerx2 = cornerx2 + offsetx;
    cornerx3 = cornerx3 + offsetx;
    cornerx4 = cornerx4 + offsetx;

    cornery1 = cornery1 + offsety;
    cornery2 = cornery2 + offsety;
    cornery3 = cornery3 + offsety;
    cornery4 = cornery4 + offsety;

    if (cornerx1<0 || cornerx1>width)
        cornerx1 = wrap(cornerx1, width);

    if (cornerx2<0 || cornerx2>width)
        cornerx2 = wrap(cornerx2, width);

    if (cornerx3<0 || cornerx3>width)
        cornerx3 = wrap(cornerx3, width);

    if (cornerx4<0 || cornerx4>width)
        cornerx4 = wrap(cornerx4, width);

    // Draw the outline of the mask.

    drawline(mask, cornerx1, cornery1, cornerx2, cornery2);
    drawline(mask, cornerx2, cornery2, cornerx3, cornery3);
    drawline(mask, cornerx3, cornery3, cornerx4, cornery4);
    drawline(mask, cornerx4, cornery4, cornerx1, cornery1);

    // Fill in the mask.

    fill(mask, masksize, masksize, masksize / 2, masksize / 2, 1);

    if (mask[0][0] == 1) // If somehow we've filled the outside of the shape rather than the inside - too weird to deal with!
        return;

    // Now we need to remove any links to ridges outside this area.

    for (int i = 0; i <= masksize; i++)
    {
        int ii = i - offsetx;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = 0; j <= masksize; j++)
        {
            int jj = j - offsety;

            if (jj >= 0 && jj <= height)
            {
                if (mask[i][j] == 1)
                {
                    // Looking north

                    int iii = i;
                    int jjj = j - 1;

                    bool docheck = 0;

                    if (jjj < 0)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 1);

                    // Looking northeast

                    iii = i + 1;
                    jjj = j - 1;

                    docheck = 0;

                    if (jjj<0 || iii>masksize)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 2);

                    // Looking east

                    iii = i + 1;
                    jjj = j;

                    docheck = 0;

                    if (iii > masksize)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 3);

                    // Looking southeast

                    iii = i + 1;
                    jjj = j + 1;

                    docheck = 0;

                    if (iii > masksize || jjj > masksize)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 4);

                    // Looking south

                    iii = i;
                    jjj = j + 1;

                    docheck = 0;

                    if (jjj > masksize)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 5);

                    // Looking southwest

                    iii = i - 1;
                    jjj = j + 1;

                    docheck = 0;

                    if (iii<0 || jjj>masksize)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 6);

                    // Looking west

                    iii = i - 1;
                    jjj = j;

                    docheck = 0;

                    if (iii < 0)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 7);

                    // Looking northwest

                    iii = i - 1;
                    jjj = j - 1;

                    docheck = 0;

                    if (iii < 0 || jjj < 0)
                        docheck = 1;

                    if (mask[iii][jjj] == 0)
                        docheck = 1;

                    if (docheck == 1)
                        terrain::detail::deleteoceanridge(world, ii, jj, 8);
                }
            }
        }
    }

    // Now we need to record all of the ridge information for this area

    vector<vector<int>> oldridgesmap(masksize + 1, vector<int>(masksize + 1, 0));
    vector<vector<int>> oldoceanrifts(masksize + 1, vector<int>(masksize + 1, 0));
    vector<vector<int>> oldoceanridges(masksize + 1, vector<int>(masksize + 1, 0));
    vector<vector<int>> oldoceanridgeheights(masksize + 1, vector<int>(masksize + 1, 0));
    vector<vector<int>> oldoceanridgeangle(masksize + 1, vector<int>(masksize + 1, 0));
    vector<vector<int>> oldoceanridgeoffset(masksize + 1, vector<int>(masksize + 1, 0));

    for (int i = 0; i <= masksize; i++)
    {
        int ii = i - offsetx;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = 0; j <= masksize; j++)
        {
            int jj = j - offsety;

            if (jj >= 0 && jj <= height)
            {
                if (mask[i][j] == 1)
                {
                    oldridgesmap[i][j] = ridgesmap[ii][jj];
                    oldoceanrifts[i][j] = world.oceanrifts(ii, jj);
                    oldoceanridges[i][j] = world.oceanridges(ii, jj);
                    oldoceanridgeheights[i][j] = world.oceanridgeheights(ii, jj);
                    oldoceanridgeangle[i][j] = world.oceanridgeangle(ii, jj);
                    oldoceanridgeoffset[i][j] = world.oceanridgeoffset(ii, jj);
                }
            }
        }
    }

    // Now we need to work out the direction and distance that all of this is going to be shifted.

    float angle = (float)(0 - (world.oceanridgeangle(x1, y1) + world.oceanridgeangle(x2, y1)) / 2);

    while (angle < 0.0f)
        angle = angle + 360.0f;

    if (deterministicrandom(faultseed ^ 0x100bull, 1, 2, midx, midy, x1, y1, x2, y2) == 1)
    {
        angle = angle + 180.0f;

        while (angle > 360.0f)
            angle = angle - 360.0f;
    }

    float fangle = angle * 0.01745329f;

    float shiftdist = (float)deterministicrandom(faultseed ^ 0x100cull, minshiftdist, maxshiftdist, midx, midy, x1, y1, x2, y2);

    int xshift = (int)(shiftdist * (float)sin(fangle));
    int yshift = (int)(shiftdist * (float)cos(fangle));

    // Now we rewrite all that information, offset by the new shift values.

    for (int i = 0; i <= masksize; i++)
    {
        int iii = i - offsetx + xshift;

        if (iii<0 || iii>width)
            iii = wrap(iii, width);

        for (int j = 0; j <= masksize; j++)
        {
            int jjj = j - offsety + yshift;

            if (jjj >= 0 && jjj <= height)
            {
                if (mask[i][j] == 1)
                {
                    ridgesmap[iii][jjj] = oldridgesmap[i][j];
                    world.setoceanrifts(iii, jjj, oldoceanrifts[i][j]);
                    world.setoceanridges(iii, jjj, oldoceanridges[i][j]);
                    world.setoceanridgeheights(iii, jjj, oldoceanridgeheights[i][j]);
                    world.setoceanridgeangle(iii, jjj, oldoceanridgeangle[i][j]);
                    world.setoceanridgeoffset(iii, jjj, oldoceanridgeoffset[i][j]);
                }
            }
        }
    }
}
}

// This creates the ocean trenches.

void createoceantrenches(planet& world, vector<vector<bool>>& shelves)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();

    int trenchmin = maxelev / 2; // Values in the fractal higher than this will spawn trenches.

    int maxradius = 20;

    int div = (maxelev - trenchmin) / maxradius;

    vector<vector<bool>> trenchmap(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));
    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int grain = 4; // Level of detail on this fractal map.
    float valuemod = 0.01f;
    int v = 1; //random(3,6);
    float valuemod2 = 0.1f;

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 2; j < height - 1; j++)
        {
            if (shelves[i][j] == 1 && fractal[i][j] > trenchmin)
            {
                if (shelfedge(world, shelves, i, j) == 1)
                {
                    int radius = (fractal[i][j] - trenchmin) / div;

                    for (int k = -radius; k <= radius; k++)
                    {
                        for (int l = -radius; l <= radius; l++)
                        {
                            if (k * k + l * l < radius * radius + radius)
                            {
                                int ll = j + l;

                                if (ll >= 0 && ll <= height)
                                {
                                    int kk = i + k;

                                    if (kk<0 || kk>width)
                                        kk = wrap(kk, width);

                                    trenchmap[kk][ll] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Now just add the trenches to the world map.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (trenchmap[i][j] == 1)
            {
                if (world.sea(i, j) == 1 && shelves[i][j] == 0 && world.oceanridges(i, j) == 0)
                {
                    int newval = world.nom(i, j) - random(4900, 5100);

                    if (newval < 1)
                        newval = 1;

                    world.setnom(i, j, newval);
                }
            }
        }
    }
}
