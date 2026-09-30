#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "lake_effects.hpp"

using namespace std;

namespace
{
void tracelakedrop(planet& world, int x, int y, int minimum, int dropno, vector<vector<int>>& thisdrop, int maxrepeat, int neighbours[8][2], vector<vector<int>>& lakesummerrainmap, vector<vector<int>>& lakewinterrainmap);
int findnearbyriver(planet& world, int x, int y, int dropno, vector<vector<int>>& thisdrop, int neighbours[8][2]);
bool swervecheck(planet& world, int x, int y, int origx, int origy, int lowered, int neighbours[8][2]);

// This function traces a drop from lake rainfall downhill to the sea (or a lake), depositing water as it goes.

void tracelakedrop(planet& world, int x, int y, int minimum, int dropno, vector<vector<int>>& thisdrop, int maxrepeat, int neighbours[8][2], vector<vector<int>>& lakesummerrainmap, vector<vector<int>>& lakewinterrainmap)
{
    twointegers newseatile;
    twointegers dest;

    int width = world.width();
    int height = world.height();
    float riverfactor = world.riverfactor();

    float janload = (float)lakewinterrainmap[x][y];
    float julload = (float)lakesummerrainmap[x][y];

    if (y >= height / 2)
    {
        janload = (float)lakesummerrainmap[x][y];
        julload = (float)lakewinterrainmap[x][y];
    }

    if ((janload + julload) / 2.0f < (float)minimum)
        return;

    janload = janload / riverfactor;
    julload = julload / riverfactor;

    bool keepgoing = 1;

    // This bit moves some water from the winter load to the summer. This is because winter precipitation feeds into the river more slowly than summer precipitation does.

    float amount = 10.0f;

    if (world.mintemp(x, y) < 0) // In cold areas, more winter water gets held over until summer because it's snow.
        amount = 3.0f;

    if (y < height / 2) // Northern hemisphere
    {
        float diff = janload / amount;

        janload = janload - diff;
        julload = julload + diff;
    }
    else // Southern hemisphere
    {
        float diff = julload / amount;

        julload = julload - diff;
        janload = janload + diff;
    }

    if (janload < 0.0f)
        janload = 0.0f;

    if (julload < 0.0f)
        julload = 0.0f;

    // Now we just trace the drop through the map, increasing the water flow wherever it goes.

    int dir = -1;
    int repeated = 0;

    do
    {
        // First check to see if there's another river nearby, if we haven't been here before.

        if (world.riverjan(x, y) == 0 && world.riverjul(x, y) == 0)
        {
            int riverdir = findnearbyriver(world, x, y, dropno, thisdrop, neighbours);

            if (riverdir != 0) // If there is a suitable river nearby
                world.setriverdir(x, y, riverdir); // Just divert the flow of this tile into it.
        }

        int newdir = world.riverdir(x, y);

        if (newdir == dir)
            repeated++;

        // Now we have a check to see if we've been going in a straight line for too long.

        if (repeated > maxrepeat && world.riverjan(x, y) == 0 && world.riverjul(x, y) == 0 && y > 0 && y < height) // If we've been going in the same direction for too long and we're not at the northern/southern edge.
        {
            int leftdir = dir - 1;
            if (leftdir < 1)
                leftdir = 8;

            int rightdir = dir + 1;
            if (rightdir > 8)
                rightdir = 1;

            int lx = x;
            int ly = y;
            int rx = x;
            int ry = y;

            // leftdir and rightdir are now possible directions we might move in.

            if (leftdir == 8 || leftdir == 1 || leftdir == 2)
                ly--;

            if (leftdir == 4 || leftdir == 5 || leftdir == 6)
                ly++;

            if (leftdir == 2 || leftdir == 3 || leftdir == 4)
                lx++;

            if (leftdir == 6 || leftdir == 7 || leftdir == 8)
                lx--;

            if (rightdir == 8 || rightdir == 1 || rightdir == 2)
                ry--;

            if (rightdir == 4 || rightdir == 5 || rightdir == 6)
                ry++;

            if (rightdir == 2 || rightdir == 3 || rightdir == 4)
                rx++;

            if (rightdir == 6 || rightdir == 7 || rightdir == 8)
                rx--;

            if (rx<0 || rx>width)
                rx = wrap(rx, width);

            if (lx<0 || lx>width)
                lx = wrap(lx, width);

            // Now we have the coordinates of the two possible alternative points to go to.

            int lowered = world.nom(x, y) - 1; // This is what we would lower the alternative point to.

            bool leftposs = swervecheck(world, lx, ly, x, y, lowered, neighbours);
            bool rightposs = swervecheck(world, rx, ry, x, y, lowered, neighbours);

            // Now we just see whether we can do the swerve, and do it.

            bool swerved = 0;

            if (leftposs == 1 && world.nom(lx, ly) < world.nom(rx, ry))
            {
                world.setriverdir(x, y, leftdir);
                newdir = leftdir;
                world.setnom(lx, ly, lowered);
                swerved = 1;
            }

            if (rightposs == 1 && world.nom(rx, ry) < world.nom(lx, ly))
            {
                world.setriverdir(x, y, rightdir);
                newdir = rightdir;
                world.setnom(rx, ry, lowered);
                swerved = 1;
            }

            if (leftposs == 1 && swerved == 0)
            {
                world.setriverdir(x, y, leftdir);
                newdir = leftdir;
                world.setnom(lx, ly, lowered);
                swerved = 1;
            }

            if (rightposs == 1 && swerved == 0)
            {
                world.setriverdir(x, y, rightdir);
                newdir = rightdir;
                world.setnom(rx, ry, lowered);
                swerved = 1;
            }
        }

        if (newdir != dir)
            repeated = 0;

        thisdrop[x][y] = dropno;

        world.setriverjan(x, y, world.riverjan(x, y) + (int)janload);
        world.setriverjul(x, y, world.riverjul(x, y) + (int)julload);

        dir = newdir;

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

        if (world.sea(x, y) == 1) // We'll continue the river into the sea for a couple of tiles, to help with the regional map generation later.
        {
            for (int n = 1; n <= 3; n++)
            {

                world.setriverjan(x, y, world.riverjan(x, y) + (int)janload);
                world.setriverjul(x, y, world.riverjul(x, y) + (int)julload);

                newseatile = findseatile(world, x, y, dir);

                int newx = newseatile.x;
                int newy = newseatile.y;

                if (newx != -1)
                {
                    if (newx == x && newy < y)
                        dir = 1;

                    if (newx > x && newy < y)
                        dir = 2;

                    if (newx > x && newy == y)
                        dir = 3;

                    if (newx > x && newy > y)
                        dir = 4;

                    if (newx == x && newy > y)
                        dir = 5;

                    if (newx<x && newy>y)
                        dir = 6;

                    if (newx < x && newy == y)
                        dir = 7;

                    if (newx < x && newy < y)
                        dir = 8;

                    world.setriverdir(x, y, dir);

                    x = newx;
                    y = newy;
                }
            }

            return;
        }

        if (world.truelake(x, y) == 1) // These rivers will stop at lakes.
            keepgoing = 0;

        if (world.riftlakesurface(x, y) != 0)
            keepgoing = 0;

        if (y == 0 || y == height)
            keepgoing = 0;

        if (thisdrop[x][y] == dropno)
        {
            keepgoing = 0;
        }

    } while (keepgoing == 1);
}

// This function sees if there's a river nearby and returns its direction.

int findnearbyriver(planet& world, int x, int y, int dropno, vector<vector<int>>& thisdrop, int neighbours[8][2])
{
    int width = world.width();
    int height = world.height();

    int start = random(0, 7);

    for (int n = start; n < start + 8; n++)
    {
        int nn = wrap(n, y);

        int i = x + neighbours[nn][0];

        if (i<0 || i>width)
            i = wrap(i, width);

        int j = y + neighbours[nn][1];

        if (j >= 0 && j <= height)
        {
            if (thisdrop[i][j] != 0 && thisdrop[i][j] != dropno) // If theres a river there, and it's not the current one
            {
                if (world.nom(i, j) < world.nom(x, y))
                    return (nn + 1);
            }
        }
    }
    return (0);
}

// This function works out whether a river could divert into the given square if its height were lowered.

bool swervecheck(planet& world, int x, int y, int origx, int origy, int lowered, int neighbours[8][2])
{
    if (world.riverjan(x, y) != 0 || world.riverjul(x, y) != 0)
        return (0);

    int width = world.width();
    int height = world.height();

    int dir = findlowestdir(world, neighbours, x, y);

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

    if (newx<0 || newx>width)
        newx = wrap(newx, width);

    if (y<0 || y>height)
        return (0);

    if (newx == origx && newy == origy) // If it's just pointing back to the original tile
        return (0);

    if (world.nom(newx, newy) >= lowered) // If it would be higher than this tile when lowered
        return (0);

    if (world.lakesurface(newx, newy) != 0) // If there's a lake here
        return (0);

    return (1);
}
}

// This adds new precipitation from the lakes.

void lakerain(planet& world, vector<vector<int>>& lakewinterrainmap, vector<vector<int>>& lakesummerrainmap)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int lakemult = tuning::climate::lakerain::lakeMultiplier;
    float tempfactor = tuning::climate::lakerain::temperatureFactor;
    float mintemp = tuning::climate::lakerain::minimumTemperatureMultiplier;
    int dumprate = tuning::climate::lakerain::dumpRate;
    int pickuprate = tuning::climate::lakerain::pickupRate;
    int swervechance = tuning::climate::lakerain::swerveChance;
    int splashsize = tuning::climate::lakerain::splashSize;
    float slopefactor = tuning::climate::lakerain::slopeFactor;
    float elevationfactor = tuning::climate::lakerain::elevationFactor;
    int slopemin = tuning::climate::lakerain::slopeMinimum;
    float seasonalvar = tuning::climate::lakerain::seasonalVariation;
    int maxrain = tuning::climate::lakerain::maxRain;
    float capfactor = tuning::climate::lakerain::capFactor;

    // First the westerly winds.

    for (int j = 0; j <= height; j++)
    {
        for (int i = 1; i <= width; i++)
        {
            int currentwind = world.wind(i, j);

            if (currentwind > 0 && currentwind != 99)
            {
                if (currentwind == 101)
                    currentwind = 4;

                if (world.map(i, j) > sealevel && world.lakesurface(i, j) == 0 && world.lakesurface(i - 1, j) != 0)// If this is next to a lake
                {
                    int crount = currentwind * lakemult;
                    int waterlog = 0; // This will hold the amount of water being carried onto this land
                    bool lakeice = 0; // This will be 1 if any seasonal ice is found when picking up moisture.

                    // First we pick up water from the lake to the west.

                    while (crount > 0)
                    {
                        int ii = i - crount;

                        if (ii < 0)
                            ii = wrap(ii, width);

                        if (world.lakesurface(ii, j) != 0 && world.maxtemp(ii, j) > 0) // It picks up water from unfrozen lake
                        {
                            if (world.mintemp(ii, j) < 5)
                                lakeice = 1;

                            float temp = world.avetemp(ii, j) / tempfactor; // Less water is picked up from colder lakes.

                            if (temp < mintemp)
                                temp = mintemp;

                            temp = (temp / 100) + 1;

                            waterlog = waterlog + (int)((float)pickuprate * temp);
                        }

                        crount--;
                    }

                    if (lakeice == 1)
                        waterlog = waterlog / 2;

                    crount = 0;
                    int jj = j;

                    // Now we deposit water onto the land to the east.

                    while (crount < currentwind * lakemult)
                    {
                        int ii = i + crount;

                        if (random(1, swervechance) == 1)
                        {
                            jj = jj + randomsign(1);

                            if (jj < 0)
                                jj = 0;

                            if (jj > height)
                                jj = height;
                        }

                        if (ii > width)
                            ii = wrap(ii, width);

                        float waterdumped = (float)waterlog / (float)dumprate;

                        if (world.map(ii, j) - sealevel > slopemin)
                        {
                            float slope = (float)getslope(world, ii - 1, j, ii, j);

                            slope = slope / slopefactor;

                            if (slope > 1.0f) // If it's going uphill
                            {
                                waterdumped = waterdumped * slope;

                                float elevation = (float(world.map(ii, j) - sealevel));

                                waterdumped = waterdumped * elevationfactor * elevation;

                                if (waterdumped > (float)waterlog)
                                    waterdumped = (float)waterlog;
                            }
                        }

                        waterlog = waterlog - (int)waterdumped;

                        if (world.map(ii, jj) > sealevel && world.wind(ii, jj) < 50 && world.lakesurface(ii, jj) == 0)
                        {
                            for (int iii = ii - splashsize; iii <= ii + splashsize; iii++)
                            {
                                int iiii = iii;

                                if (iiii<0 || iiii>width)
                                    iiii = wrap(iiii, width);

                                for (int jjj = jj - splashsize; jjj <= jj + splashsize; jjj++)
                                {
                                    int jjjj = jjj;

                                    if (jjjj < 0)
                                        jjjj = 0;

                                    if (jjjj > height)
                                        jjjj = height;

                                    if (lakesummerrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakesummerrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                    if (lakewinterrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakewinterrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                }
                            }

                            if (lakesummerrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;

                            if (lakewinterrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;
                        }

                        crount++;

                        if (waterlog < 50)
                            crount = currentwind * lakemult;
                    }

                }
            }
        }
    }

    // Now the easterly winds.

    for (int j = 0; j <= height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int currentwind = world.wind(i, j);

            if (currentwind < 0 || currentwind == 99)
            {
                if (currentwind == 99)
                    currentwind = -4;

                if (world.map(i, j) > sealevel && world.lakesurface(i, j) == 0 && world.lakesurface(i + 1, j) != 0)// If this is next to a lake
                {
                    int crount = 0 - (currentwind * lakemult);
                    int waterlog = 0; // This will hold the amount of water being carried onto this land
                    bool lakeice = 0; // This will be 1 if any seasonal ice is found when picking up moisture.

                    // First we pick up water from the lake to the east.

                    while (crount > 0)
                    {
                        int ii = i + crount;

                        if (ii > width)
                            ii = wrap(ii, width);

                        if (world.lakesurface(ii, j) != 0 && world.maxtemp(ii, j) > 0) // It picks up water from unfrozen lake
                        {
                            if (world.mintemp(ii, j) < 5)
                                lakeice = 1;

                            float temp = world.avetemp(ii, j) / tempfactor; // Less water is picked up from colder lakes.

                            if (temp < mintemp)
                                temp = mintemp;

                            temp = (temp / 100) + 1;

                            waterlog = waterlog + (int)((float)pickuprate * temp);
                        }

                        crount--;
                    }

                    if (lakeice == 1)
                        waterlog = waterlog / 2;

                    crount = 0;
                    int jj = j;

                    // Now we deposit water onto the land to the west.

                    while (crount < 0 - (currentwind * lakemult))
                    {
                        int ii = i - crount;

                        if (random(1, swervechance) == 1)
                        {
                            jj = jj + randomsign(1);

                            if (jj < 0)
                                jj = 0;

                            if (jj > height)
                                jj = height;
                        }

                        if (ii < 0)
                            ii = wrap(ii, width);

                        float waterdumped = (float)waterlog / (float)dumprate;

                        if (world.map(ii, j) - sealevel > slopemin)
                        {
                            float slope = (float)getslope(world, ii - 1, j, ii, j);

                            slope = slope / slopefactor;

                            if (slope > 1.0f) // If it's going uphill
                            {
                                waterdumped = waterdumped * slope;

                                float elevation = (float(world.map(ii, j) - sealevel));

                                waterdumped = waterdumped * elevationfactor * elevation;

                                if (waterdumped > (float)waterlog)
                                    waterdumped = (float)waterlog;
                            }
                        }

                        waterlog = waterlog - (int)waterdumped;

                        if (world.map(ii, jj) > sealevel && world.wind(ii, jj) < 50 && world.lakesurface(ii, jj) == 0)
                        {
                            for (int iii = ii - splashsize; iii <= ii + splashsize; iii++)
                            {
                                int iiii = iii;

                                if (iiii<0 || iiii>width)
                                    iiii = wrap(iiii, width);

                                for (int jjj = jj - splashsize; jjj <= jj + splashsize; jjj++)
                                {
                                    int jjjj = jjj;

                                    if (jjjj < 0)
                                        jjjj = 0;

                                    if (jjjj > height)
                                        jjjj = height;

                                    if (lakesummerrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakesummerrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                    if (lakewinterrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakewinterrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);
                                }
                            }

                            if (lakesummerrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;

                            if (lakewinterrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;
                        }

                        crount++;

                        if (waterlog < 50)
                            crount = 0 - (currentwind * lakemult);
                    }
                }
            }
        }
    }

    // Now adjust for seasonal variation on land.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.map(i, j) > sealevel)
                {
                    float tempdiff = (float)(world.maxtemp(i, j) - world.mintemp(i, j));

                    if (tempdiff < 0.0f)
                        tempdiff = 0.0f;

                    float winterrain = (float)lakewinterrainmap[i][j];

                    tempdiff = tempdiff * seasonalvar * winterrain;

                    lakewinterrainmap[i][j] = lakewinterrainmap[i][j] + (int)tempdiff;
                    lakesummerrainmap[i][j] = lakesummerrainmap[i][j] - (int)tempdiff;

                    if (lakesummerrainmap[i][j] < 0)
                        lakesummerrainmap[i][j] = 0;
                }
            }
        }
    });

    // Now cap excessive rainfall.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                float rain[2];
                rain[0] = (float)lakewinterrainmap[i][j];
                rain[1] = (float)lakesummerrainmap[i][j];

                for (int n = 0; n <= 1; n++)
                {
                    if (rain[n] > (float)maxrain)
                        rain[n] = ((rain[n] - (float)maxrain) * capfactor) + (float)maxrain;

                    if (rain[n] < 0.0f)
                        rain[n] = 0.0f;
                }

                lakewinterrainmap[i][j] = (int)rain[0];
                lakesummerrainmap[i][j] = (int)rain[1];
            }
        }
    });

    // Now blur the lake rain.

    smooth(lakewinterrainmap, width, height, maxelev, 1, 0);
    smooth(lakesummerrainmap, width, height, maxelev, 1, 0);

    // Now add the lake rain to the existing rain maps.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                world.setwinterrain(i, j, world.winterrain(i, j) + lakewinterrainmap[i][j]);
                world.setsummerrain(i, j, world.summerrain(i, j) + lakesummerrainmap[i][j]);
            }
        }
    });
}

// This does the same thing, but for rift lakes.

void riftlakerain(planet& world, vector<vector<int>>& lakewinterrainmap, vector<vector<int>>& lakesummerrainmap)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int lakemult = tuning::climate::riftlakerain::lakeMultiplier;
    float tempfactor = tuning::climate::riftlakerain::temperatureFactor;
    float mintemp = tuning::climate::riftlakerain::minimumTemperatureMultiplier;
    int dumprate = tuning::climate::riftlakerain::dumpRate;
    int pickuprate = tuning::climate::riftlakerain::pickupRate;
    int swervechance = tuning::climate::riftlakerain::swerveChance;
    int splashsize = tuning::climate::riftlakerain::splashSize;
    int slopefactor = tuning::climate::riftlakerain::slopeFactor;
    float elevationfactor = tuning::climate::riftlakerain::elevationFactor;
    int slopemin = tuning::climate::riftlakerain::slopeMinimum;
    float seasonalvar = tuning::climate::riftlakerain::seasonalVariation;
    int maxrain = tuning::climate::riftlakerain::maxRain;
    float capfactor = tuning::climate::riftlakerain::capFactor;

    // First the westerly winds.

    for (int j = 0; j <= height; j++)
    {
        for (int i = 1; i <= width; i++)
        {
            int currentwind = world.wind(i, j);

            if (currentwind > 0 && currentwind != 99)
            {
                if (currentwind == 101)
                    currentwind = 4;

                if (world.map(i, j) > sealevel && world.riftlakesurface(i, j) == 0 && world.riftlakesurface(i - 1, j) != 0)// If this is next to a lake
                {
                    int crount = currentwind * lakemult;
                    int waterlog = 0; // This will hold the amount of water being carried onto this land
                    bool lakeice = 0; // This will be 1 if any seasonal ice is found when picking up moisture.

                    // First we pick up water from the lake to the west.

                    while (crount > 0)
                    {
                        int ii = i - crount;

                        if (ii < 0)
                            ii = wrap(ii, width);

                        if (world.riftlakesurface(ii, j) != 0 && world.maxtemp(ii, j) > 0) // It picks up water from unfrozen lake
                        {
                            if (world.mintemp(ii, j) < 5)
                                lakeice = 1;

                            float temp = world.avetemp(ii, j) / tempfactor; // Less water is picked up from colder lakes.

                            if (temp < mintemp)
                                temp = mintemp;

                            temp = (temp / 100) + 1;

                            waterlog = waterlog + (int)((float)pickuprate * temp);
                        }

                        crount--;
                    }

                    if (lakeice == 1)
                        waterlog = waterlog / 2;

                    crount = 0;
                    int jj = j;

                    // Now we deposit water onto the land to the east.

                    while (crount < currentwind * lakemult)
                    {
                        int ii = i + crount;

                        if (random(1, swervechance) == 1)
                        {
                            jj = jj + randomsign(1);

                            if (jj < 0)
                                jj = 0;

                            if (jj > height)
                                jj = height;
                        }

                        if (ii > width)
                            ii = wrap(ii, width);

                        float waterdumped = (float)waterlog / (float)dumprate;

                        if (world.map(ii, j) - sealevel > slopemin)
                        {
                            float slope = (float)getslope(world, ii - 1, j, ii, j);

                            slope = slope / slopefactor;

                            if (slope > 1.0f) // If it's going uphill
                            {
                                waterdumped = waterdumped * slope;

                                float elevation = (float)(world.map(ii, j) - sealevel);

                                waterdumped = waterdumped * elevationfactor * elevation;

                                if (waterdumped > (float)waterlog)
                                    waterdumped = (float)waterlog;
                            }
                        }

                        waterlog = waterlog - (int)waterdumped;

                        if (world.map(ii, jj) > sealevel && world.wind(ii, jj) < 50 && world.riftlakesurface(ii, jj) == 0)
                        {
                            for (int iii = ii - splashsize; iii <= ii + splashsize; iii++)
                            {

                                int iiii = iii;

                                if (iiii<0 || iiii>width)
                                    iiii = wrap(iiii, width);

                                for (int jjj = jj - splashsize; jjj <= jj + splashsize; jjj++)
                                {
                                    int jjjj = jjj;

                                    if (jjjj < 0)
                                        jjjj = 0;

                                    if (jjjj > height)
                                        jjjj = height;

                                    if (lakesummerrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakesummerrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                    if (lakewinterrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakewinterrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                }
                            }

                            if (lakesummerrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;

                            if (lakewinterrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;
                        }

                        crount++;

                        if (waterlog < 50)
                            crount = currentwind * lakemult;
                    }

                }
            }
        }
    }

    // Now the easterly winds.

    for (int j = 0; j <= height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            int currentwind = world.wind(i, j);

            if (currentwind < 0 || currentwind == 99)
            {
                if (currentwind == 99)
                    currentwind = -4;

                if (world.map(i, j) > sealevel && world.riftlakesurface(i, j) == 0 && world.riftlakesurface(i + 1, j) != 0)// If this is next to a lake
                {
                    int crount = 0 - (currentwind * lakemult);
                    int waterlog = 0; // This will hold the amount of water being carried onto this land
                    bool lakeice = 0; // This will be 1 if any seasonal ice is found when picking up moisture.

                    // First we pick up water from the lake to the east.

                    while (crount > 0)
                    {
                        int ii = i + crount;

                        if (ii > width)
                            ii = wrap(ii, width);

                        if (world.riftlakesurface(ii, j) != 0 && world.maxtemp(ii, j) > 0) // It picks up water from unfrozen lake
                        {
                            if (world.mintemp(ii, j) < 5)
                                lakeice = 1;

                            float temp = world.avetemp(ii, j) / tempfactor; // Less water is picked up from colder lakes.

                            if (temp < mintemp)
                                temp = mintemp;

                            temp = (temp / 100) + 1;

                            waterlog = waterlog + (int)((float)pickuprate * temp);
                        }

                        crount--;
                    }

                    if (lakeice == 1)
                        waterlog = waterlog / 2;

                    crount = 0;
                    int jj = j;

                    // Now we deposit water onto the land to the west.

                    while (crount < 0 - (currentwind * lakemult))
                    {
                        int ii = i - crount;

                        if (random(1, swervechance) == 1)
                        {
                            jj = jj + randomsign(1);

                            if (jj < 0)
                                jj = 0;

                            if (jj > height)
                                jj = height;
                        }

                        if (ii < 0)
                            ii = wrap(ii, width);

                        float waterdumped = (float)waterlog / (float)dumprate;

                        if (world.map(ii, j) - sealevel > slopemin)
                        {
                            float slope = (float)getslope(world, ii - 1, j, ii, j);

                            slope = slope / slopefactor;

                            if (slope > 1.0f) // If it's going uphill
                            {
                                waterdumped = waterdumped * slope;

                                float elevation = (float)(world.map(ii, j) - sealevel);

                                waterdumped = waterdumped * elevationfactor * elevation;

                                if (waterdumped > (float)waterlog)
                                    waterdumped = (float)waterlog;
                            }
                        }

                        waterlog = waterlog - (int)waterdumped;

                        if (world.map(ii, jj) > sealevel && world.wind(ii, jj) < 50 && world.riftlakesurface(ii, jj) == 0)
                        {
                            for (int iii = ii - splashsize; iii <= ii + splashsize; iii++)
                            {
                                int iiii = iii;

                                if (iiii<0 || iiii>width)
                                    iiii = wrap(iiii, width);

                                for (int jjj = jj - splashsize; jjj <= jj + splashsize; jjj++)
                                {
                                    int jjjj = jjj;

                                    if (jjjj < 0)
                                        jjjj = 0;

                                    if (jjjj > height)
                                        jjjj = height;

                                    if (lakesummerrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakesummerrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);

                                    if (lakewinterrainmap[iiii][jjjj] < (int)(waterdumped / 2.0f))
                                        lakewinterrainmap[iiii][jjjj] = (int)(waterdumped / 2.0f);
                                }
                            }

                            if (lakesummerrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;

                            if (lakewinterrainmap[ii][jj] < (int)waterdumped)
                                lakesummerrainmap[ii][jj] = (int)waterdumped;
                        }

                        crount++;

                        if (waterlog < 50)
                            crount = 0 - (currentwind * lakemult);
                    }
                }
            }
        }
    }

    // Now adjust for seasonal variation on land.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.map(i, j) > sealevel)
                {
                    float tempdiff = (float)(world.maxtemp(i, j) - world.mintemp(i, j));

                    if (tempdiff < 0.0f)
                        tempdiff = 0.0f;

                    float winterrain = (float)lakewinterrainmap[i][j];

                    tempdiff = tempdiff * seasonalvar * winterrain;

                    lakewinterrainmap[i][j] = lakewinterrainmap[i][j] + (int)tempdiff;
                    lakesummerrainmap[i][j] = lakesummerrainmap[i][j] - (int)tempdiff;

                    if (lakesummerrainmap[i][j] < 0)
                        lakesummerrainmap[i][j] = 0;
                }
            }
        }
    });

    // Now cap excessive rainfall.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                float rain[2];
                rain[0] = (float)lakewinterrainmap[i][j];
                rain[1] = (float)lakesummerrainmap[i][j];

                for (int n = 0; n <= 1; n++)
                {
                    if (rain[n] > (float)maxrain)
                        rain[n] = ((rain[n] - (float)maxrain) * capfactor) + (float)maxrain;

                    if (rain[n] < 0.0f)
                        rain[n] = 0.0f;
                }

                lakewinterrainmap[i][j] = (int)rain[0];
                lakesummerrainmap[i][j] = (int)rain[1];
            }
        }
    });

    // Now blur the lake rain.

    smooth(lakewinterrainmap, width, height, maxelev, 1, 0);
    smooth(lakesummerrainmap, width, height, maxelev, 1, 0);

    // Now add the lake rain to the existing rain maps.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                world.setwinterrain(i, j, world.winterrain(i, j) + lakewinterrainmap[i][j]);
                world.setsummerrain(i, j, world.summerrain(i, j) + lakesummerrainmap[i][j]);
            }
        }
    });
}

// This adds new rivers from the lake-related rainfall.

void lakerivers(planet& world, vector<vector<int>>& lakewinterrainmap, vector<vector<int>>& lakesummerrainmap)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int minimum = 40; // We won't bother with tiles with lower average flow than this.
    int maxrepeat = 2; // If a river goes in the same direction for this long, it will try to change course.

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

    // Now, we go through the map tile by tile. Take the lake rainfall in each tile and add it to every downstream tile.

    vector<vector<int>> thisdrop(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int dropno = 1;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.map(i, j) > sealevel)
            {
                tracelakedrop(world, i, j, minimum, dropno, thisdrop, maxrepeat, neighbours, lakesummerrainmap, lakewinterrainmap);
                dropno++;
            }
        }
    }
}
