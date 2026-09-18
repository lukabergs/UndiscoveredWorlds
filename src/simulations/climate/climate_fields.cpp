#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "world_generation_debug.hpp"
#include "climate_fields.hpp"
#include "climate_coupling.hpp"
#include "climate_energy.hpp"
#include "land_distance.hpp"

using namespace std;

namespace
{
struct maritimesample
{
    float influence = 0.0f;
    float thermalanomaly = 0.0f;
    int fetchdistance = 0;
};

int coldseasonindex(planet& world, int x, int y);
int warmseasonindex(planet& world, int x, int y);
int coldoceansst(planet& world, int x, int y);
int warmoceansst(planet& world, int x, int y);
maritimesample sampleupwindmaritimeinfluence(planet& world, int season, int x, int y);
maritimesample storedmaritimeinfluence(planet& world, int season, int x, int y);
void storescaledmaritimeinfluence(planet& world, int season, int x, int y, const maritimesample& sample);
void createcoastalclimateinfluence(planet& world);
void populateinlanddistances(planet& world, vector<vector<int>>& inland);
void createdesertworldrain(planet& world);
void adjustseasonalrainfall(planet& world, vector<vector<int>>& inland);
void smoothrainfall(planet& world, int maxmountainheight);
void caprainfall(planet& world);
void applycoastalclimates(planet& world);
void smoothtemperatures(planet& world);

int coldseasonindex(planet& world, int x, int y)
{
    return (world.jantemp(x, y) <= world.jultemp(x, y)) ? seasonjanuary : seasonjuly;
}

int warmseasonindex(planet& world, int x, int y)
{
    return (world.jantemp(x, y) >= world.jultemp(x, y)) ? seasonjanuary : seasonjuly;
}

int coldoceansst(planet& world, int x, int y)
{
    return world.seasonalsst(coldseasonindex(world, x, y), x, y);
}

int warmoceansst(planet& world, int x, int y)
{
    return world.seasonalsst(warmseasonindex(world, x, y), x, y);
}

maritimesample sampleupwindmaritimeinfluence(planet& world, int season, int x, int y)
{
    maritimesample result;

    if (world.sea(x, y) == 1)
        return result;

    const int width = world.width();
    const int height = world.height();
    const int maxsearchdistance = tuning::climateresolution::scaleDistance(
        tuning::climate::maritime::maxSearchDistance, width, height);
    const int oceansamplecount = tuning::climateresolution::scaleDistance(
        tuning::climate::maritime::oceanSampleCount, width, height);

    const float u = static_cast<float>(world.seasonaluwind(season, x, y));
    const float v = static_cast<float>(world.seasonalvwind(season, x, y));
    const float magnitude = std::sqrt(u * u + v * v);

    if (magnitude < tuning::climate::maritime::minimumWindStrength)
    {
        result.fetchdistance = maxsearchdistance;
        return result;
    }

    const float dirx = -u / magnitude;
    const float diry = -v / magnitude;

    float posx = static_cast<float>(x);
    float posy = static_cast<float>(y);
    int landsteps = 0;
    int seacells = 0;
    float ssttotal = 0.0f;
    float evaporationtotal = 0.0f;
    bool foundsea = false;

    for (int step = 1; step <= maxsearchdistance + oceansamplecount; step++)
    {
        posx = posx + dirx;
        posy = posy + diry;

        int xx = static_cast<int>(std::round(posx));
        int yy = static_cast<int>(std::round(posy));

        if (xx < 0 || xx > width)
            xx = wrap(xx, width);

        yy = std::clamp(yy, 0, height);

        if (foundsea == false)
        {
            if (world.sea(xx, yy) == 1 && world.seaice(xx, yy) != 2)
            {
                foundsea = true;
                result.fetchdistance = landsteps;
            }
            else
            {
                landsteps++;
                continue;
            }
        }

        if (world.sea(xx, yy) == 0)
            break;

        ssttotal = ssttotal + static_cast<float>(world.seasonalsst(season, xx, yy));
        evaporationtotal = evaporationtotal + static_cast<float>(world.seasonalevaporation(season, xx, yy));
        seacells++;

        if (seacells >= oceansamplecount)
            break;
    }

    if (foundsea == false || seacells == 0)
    {
        result.fetchdistance = maxsearchdistance;
        return result;
    }

    const float averagedsstonsea = ssttotal / static_cast<float>(seacells);
    const float averageevaporation = evaporationtotal / static_cast<float>(seacells);
    const float fetchfactor = 1.0f - std::clamp(static_cast<float>(result.fetchdistance) / static_cast<float>(maxsearchdistance), 0.0f, 1.0f);
    const float windfactor = std::clamp(magnitude / tuning::climate::maritime::windScale, 0.0f, 1.0f);
    const float continuityfactor = std::clamp(static_cast<float>(seacells) / static_cast<float>(oceansamplecount), 0.0f, 1.0f);
    const float moisturefactor = 0.5f + std::clamp(averageevaporation / tuning::climate::maritime::evaporationScale, 0.0f, 0.5f);

    result.influence = fetchfactor * windfactor * continuityfactor * moisturefactor;
    result.thermalanomaly = averagedsstonsea - static_cast<float>(world.seasonaltemp(season, x, y));
    result.thermalanomaly = std::clamp(result.thermalanomaly, -tuning::climate::maritime::maxThermalAnomaly, tuning::climate::maritime::maxThermalAnomaly);

    return result;
}

maritimesample storedmaritimeinfluence(planet& world, int season, int x, int y)
{
    maritimesample result;
    result.influence = static_cast<float>(world.seasonalmaritimeinfluence(season, x, y)) / tuning::climate::maritime::influenceStorageScale;
    result.thermalanomaly = static_cast<float>(world.seasonalmaritimethermalanomaly(season, x, y)) / tuning::climate::maritime::thermalAnomalyStorageScale;
    result.fetchdistance = world.seasonalmaritimefetch(season, x, y);
    return result;
}

void storescaledmaritimeinfluence(planet& world, int season, int x, int y, const maritimesample& sample)
{
    world.setseasonalmaritimeinfluence(season, x, y, static_cast<int>(roundf(sample.influence * tuning::climate::maritime::influenceStorageScale)));
    world.setseasonalmaritimethermalanomaly(season, x, y, static_cast<int>(roundf(sample.thermalanomaly * tuning::climate::maritime::thermalAnomalyStorageScale)));
    world.setseasonalmaritimefetch(season, x, y, sample.fetchdistance);
}

void createcoastalclimateinfluence(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.sea(x, y) == 1)
                    {
                        world.setseasonalmaritimeinfluence(season, x, y, 0);
                        world.setseasonalmaritimethermalanomaly(season, x, y, 0);
                        world.setseasonalmaritimefetch(season, x, y, 0);
                    }
                    else
                    {
                        storescaledmaritimeinfluence(world, season, x, y, sampleupwindmaritimeinfluence(world, season, x, y));
                    }
                }
            }
        });
    }
}

void populateinlanddistances(planet& world, vector<vector<int>>& inland)
{
    const int width = world.width();
    const int height = world.height();
    const int maximumdistance = min(30000, width + height);
    vector<vector<short>> inlanddistances(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, static_cast<short>(maximumdistance + 1)));

    const bool hascoast = terraindetail::buildlanddistancefield(world, inlanddistances, maximumdistance, [&](int x, int y)
    {
        return world.outline(x, y);
    });

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1)
            {
                inland[i][j] = 0;
                continue;
            }

            if (hascoast == false)
                inland[i][j] = maximumdistance;
            else
                inland[i][j] = inlanddistances[i][j];
        }
    }
}

// This adds a bit of rainfall to worlds with no sea.

void createdesertworldrain(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    float slopefactor = tuning::climate::desertworld::slopeFactor;
    float idealtemp = tuning::climate::desertworld::idealTemperature;
    float maxdiff = tuning::climate::desertworld::maxTemperatureDifference;

    float maxrain = (float)(random(tuning::climate::desertworld::maxRainMin, tuning::climate::desertworld::maxRainMax));

    int grain = tuning::climate::desertworld::fractalGrain;
    float valuemod = tuning::climate::desertworld::fractalValueMod;
    int v = random(tuning::climate::desertworld::fractalValueMod2Min, tuning::climate::desertworld::fractalValueMod2Max);
    float valuemod2 = float(v);

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    int warpfactor = tuning::climate::desertworld::warpFactor;
    warp(fractal, width, height, (int)maxelev, warpfactor, 0);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                int elev = world.map(i, j);
                int biggestdiff = 0;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;
                    if (kk<0 || k>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            int thisdiff = elev - world.map(kk, l);

                            if (thisdiff > biggestdiff)
                                biggestdiff = thisdiff;
                        }
                    }
                }

                float thisslopemult = (float)biggestdiff / slopefactor;

                float thisfractalmult = (float)fractal[i][j] / maxelev;

                float thisrain = maxrain * thisslopemult * thisfractalmult;

                world.setjanrain(i, j, (int)thisrain);
                world.setjulrain(i, j, (int)thisrain);
            }
        }
    });

    // Now adjust for seasonal variation.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                float oldjanrain = (float)world.janrain(i, j);
                float oldjulrain = (float)world.julrain(i, j);
                
                float jandiff = (float)world.jantemp(i, j) - idealtemp;
                float juldiff = (float)world.jultemp(i, j) - idealtemp;

                if (jandiff < 0.0f)
                    jandiff = 0.0f - jandiff;

                if (juldiff < 0.0f)
                    juldiff = 0.0f - juldiff;

                jandiff = maxdiff - jandiff;
                juldiff = maxdiff - juldiff;

                if (jandiff < 0.0f)
                    jandiff = 0.0f;

                if (juldiff < 0.0f)
                    juldiff = 0.0f;

                float janmult = jandiff / maxdiff;
                float julmult = juldiff / maxdiff;

                float newjanrain = oldjanrain * janmult;
                float newjulrain = oldjulrain * julmult;
                
                world.setjanrain(i, j, (int)newjanrain);
                world.setjulrain(i, j, (int)newjulrain);
            }
        }
    });

    // Now blur.

    int dist = tuning::climate::desertworld::blurDistance;

    vector<vector<int>> finaljanrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> finaljulrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                int crount = 0;
                int jantotal = 0;
                int jultotal = 0;

                for (int k = i - dist; k <= i + dist; k++)
                {
                    int kk = k;
                    if (kk<0 || k>width)
                        kk = wrap(kk, width);

                    for (int l = j - dist; l <= j + dist; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            jantotal = jantotal + world.janrain(i, j);
                            jultotal = jultotal + world.julrain(i, j);

                            crount++;
                        }
                    }
                }

                float newjanrain = (float)jantotal / (float)crount;
                float newjulrain = (float)jultotal / (float)crount;

                finaljanrain[i][j] = (int)newjanrain;
                finaljulrain[i][j] = (int)newjulrain;
            }
        }
    });

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                world.setjanrain(i, j, finaljanrain[i][j]);
                world.setjulrain(i, j, finaljulrain[i][j]);
            }
        }
    });
}

// This adjusts the seasonal rainfall (to ensure certain climates)

void adjustseasonalrainfall(planet& world, vector<vector<int>>& inland)
{
    (void)world;
    (void)inland;

    // Seasonal rainfall is now expected to emerge from the pressure, wind, and moisture-advection passes.
    return;

    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();
    float tilt = world.tilt();

    // First, do some tinkering to encourage Mediterranean climates.

    vector<vector<float>> mediterranean(ARRAYWIDTH, vector<float>(ARRAYHEIGHT, 0));

    float medstrength = 1.8f - (abs(tuning::climate::mediterranean::strengthCenterTilt - tilt) / tuning::climate::mediterranean::strengthTiltDivisor);

    if (medstrength > 0.0f)
    {
        int avemedmaxtemp = tuning::climate::mediterranean::targetMaxTemperature;
        float medmaxtempdifffactor = tuning::climate::mediterranean::maxTemperatureDifferenceFactor;
        float minmedcoldtemp = tuning::climate::mediterranean::minimumColdTemperature;
        float medmintempdifffactor = tuning::climate::mediterranean::minimumColdTemperatureFactor;
        int maxmedrain = tuning::climate::mediterranean::maxRain;
        float medmaxraindifffactor = tuning::climate::mediterranean::maxRainDifferenceFactor;
        int maxmedinland = tuning::climate::mediterranean::maxInlandDistance;
        float medmaxinlanddifffactor = tuning::climate::mediterranean::maxInlandDifferenceFactor;
        float medhorsedifffactor = tuning::climate::mediterranean::horseLatitudeDifferenceFactor;

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                if (world.sea(i, j) == 0)
                {
                    float med = 1.0f;

                    float diff = (float)(abs(world.maxtemp(i, j) - avemedmaxtemp));

                    med = med - diff * medmaxtempdifffactor;

                    diff = minmedcoldtemp - world.mintemp(i, j);

                    if (diff > 0.0f)
                        med = med - diff * medmintempdifffactor;

                    diff = world.averainfloat(i, j) - static_cast<float>(maxmedrain);

                    if (diff > 0.0f)
                        med = med - diff * medmaxraindifffactor;

                    diff = (float)(inland[i][j] - maxmedinland);

                    if (diff > 0.0f)
                        med = med - diff * medmaxinlanddifffactor;

                    if (j < height / 2)
                    {
                        diff = (float)abs(world.horse(i, 1) - j);

                        if (j > world.horse(i, 1))
                            diff = diff * 2.0f;
                    }
                    else
                    {
                        diff = (float)abs(world.horse(i, 4) - j);

                        if (j < world.horse(i, 4))
                            diff = diff * 2.0f;
                    }

                    med = med - diff * medhorsedifffactor;

                    if (med > 0.0f)
                        mediterranean[i][j] = med;

                }
            }
        }

        int dist = tuning::climate::mediterranean::smoothDistance;

        for (int i = 0; i <= width; i++) // Smooth it
        {
            for (int j = 0; j <= height; j++)
            {
                float total = 0;
                short crount = 0;

                for (int k = i - dist; k <= i + dist; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - dist; l <= j + dist; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            total = total + mediterranean[kk][l];
                            crount++;
                        }
                    }
                }

                if (crount > 0)
                    mediterranean[i][j] = total / crount;
            }
        }

        for (int i = 0; i <= width; i++) // Now use it to alter the seasonal variation in rainfall.
        {
            for (int j = 0; j <= height; j++)
            {
                if (mediterranean[i][j] != 0)
                {
                    float winterrain = (float)world.winterrain(i, j);
                    float summerrain = (float)world.summerrain(i, j);

                    float rainshift = summerrain * medstrength * mediterranean[i][j];

                    winterrain = winterrain + rainshift;
                    summerrain = summerrain - rainshift;

                    world.setwinterrain(i, j, (int)winterrain);
                    world.setsummerrain(i, j, (int)summerrain);
                }
            }
        }
    }

    // Now some tinkering to encourage rainforest near the equator.

    float rainstrength = 0.0f; // Seasonal equatorial enhancement now comes from the moisture-advection solver.

    if (rainstrength > 0.0f)
    {
        float equator = (float)(height / 2);

        float winteradd = tuning::climate::equatorialrain::winterAdditionFactor * rainstrength;
        float summeradd = tuning::climate::equatorialrain::summerAdditionFactor * rainstrength;

        for (int i = 0; i <= width; i++)
        {
            float tropicnorth = (float)world.horse(i, 2);
            float tropicnorthheight = equator - tropicnorth;
            float winterstep = winteradd / tropicnorthheight;
            float summerstep = summeradd / tropicnorthheight;

            float currentwinteradd = 0.0f;
            float currentsummeradd = 0.0f;

            for (int j = (int)tropicnorth; j <= (int)equator; j++)
            {
                currentwinteradd = currentwinteradd + winterstep;
                currentsummeradd = currentsummeradd + summerstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentwinteradd;
                currentjul = currentjul + currentjul * currentsummeradd;

                world.setjanrain(i, j, (int)currentjan);
                world.setjulrain(i, j, (int)currentjul);
            }

            float tropicsouth = (float)world.horse(i, 3);
            float tropicsouthheight = tropicsouth - equator;
            winterstep = winteradd / tropicsouthheight;
            summerstep = summeradd / tropicsouthheight;

            currentwinteradd = 0.0f;
            currentsummeradd = 0.0f;

            for (int j = (int)tropicsouth; j > (int)equator; j--)
            {
                currentwinteradd = currentwinteradd + winterstep;
                currentsummeradd = currentsummeradd + summerstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentsummeradd;
                currentjul = currentjul + currentjul * currentwinteradd;

                world.setjanrain(i, j, (int)currentjan);
                world.setjulrain(i, j, (int)currentjul);
            }
        }

        /*

        // And now some more tinkering, to encourage savannah nearer the edge of the tropics. (Turned off now as it isn't really needed.)

        float janadd = -0.8 * rainstrength; // 0.8f;
        float juladd = 0 * rainstrength; // 1.0f; // 0.6f;

        float bandwidth = 30; // The strip affected will extend by this much to the north and south of the boundary of the horse latitudes.

        float janstep = janadd / bandwidth;
        float julstep = juladd / bandwidth;

        for (int i = 0; i <= width; i++)
        {
            float tropicnorth = (world.horse(i, 2) + equator) / 2;
            float tropicsouth = (world.horse(i, 3) + equator) / 2;

            float currentjanadd = 0;
            float currentjuladd = 0;

            for (int j = tropicnorth - bandwidth; j <= tropicnorth; j++)
            {
                currentjanadd = currentjanadd + janstep;
                currentjuladd = currentjuladd + julstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentjanadd;
                currentjul = currentjul + currentjul * currentjuladd;

                int newjan = (int)currentjan;
                int newjul = (int)currentjul;

                world.setjanrain(i, j, newjan);
                world.setjulrain(i, j, newjul);
            }

            currentjanadd = 0;
            currentjuladd = 0;

            for (int j = tropicnorth + bandwidth; j > tropicnorth; j--)
            {
                currentjanadd = currentjanadd + janstep;
                currentjuladd = currentjuladd + julstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentjanadd;
                currentjul = currentjul + currentjul * currentjuladd;

                int newjan = (int)currentjan;
                int newjul = (int)currentjul;

                world.setjanrain(i, j, newjan);
                world.setjulrain(i, j, newjul);
            }

            currentjanadd = 0;
            currentjuladd = 0;

            for (int j = tropicsouth - bandwidth; j <= tropicsouth; j++)
            {
                currentjanadd = currentjanadd + janstep;
                currentjuladd = currentjuladd + julstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentjuladd; // Other way round for southern hemisphere!
                currentjul = currentjul + currentjul * currentjanadd;

                int newjan = (int)currentjan;
                int newjul = (int)currentjul;

                world.setjanrain(i, j, newjan);
                world.setjulrain(i, j, newjul);
            }

            currentjanadd = 0;
            currentjuladd = 0;

            for (int j = tropicsouth + bandwidth; j > tropicsouth; j--)
            {
                currentjanadd = currentjanadd + janstep;
                currentjuladd = currentjuladd + julstep;

                float currentjan = float(world.janrain(i, j));
                float currentjul = float(world.julrain(i, j));

                currentjan = currentjan + currentjan * currentjuladd; // Other way round for southern hemisphere!
                currentjul = currentjul + currentjul * currentjanadd;

                int newjan = (int)currentjan;
                int newjul = (int)currentjul;

                world.setjanrain(i, j, newjan);
                world.setjulrain(i, j, newjul);
            }
        }
        */
    }
}

// This smooths the rainfall.

void smoothrainfall(planet& world, int maxmountainheight)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    vector<vector<int>> smoothedjanrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> smoothedjulrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int n = 0; n < 5; n++) // Quick basic smooth, which will include smoothing from one mountain tile to the next.
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    int crount = 0;
                    int jantotal = 0;
                    int jultotal = 0;
                    int centreheight = world.map(i, j);

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (abs(centreheight - world.map(kk, l)) < 150)
                                {
                                    jantotal = jantotal + world.janrain(kk, l);
                                    jultotal = jultotal + world.julrain(kk, l);
                                    crount++;
                                }
                            }
                        }
                    }

                    jantotal = jantotal / crount;
                    jultotal = jultotal / crount;

                    smoothedjanrain[i][j] = jantotal;
                    smoothedjulrain[i][j] = jultotal;
                }
            }
        });

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    world.setjanrain(i, j, smoothedjanrain[i][j]);
                    world.setjulrain(i, j, smoothedjulrain[i][j]);
                }
            }
        });
    }

    for (int n = 0; n < 5; n++) // Smooth with wider scope, avoiding mountains altogether.
    {
        for (int i = 0; i <= width; i++) // Normal blur.
        {
            int im = i - 1;
            int imm = i - 2;
            int ip = i + 1;
            int ipp = i + 2;

            if (im < 0)
                im = width;

            if (imm < 0)
                imm = wrap(imm, width);

            if (ip > width)
                ip = 0;

            if (ipp > width)
                ipp = wrap(ipp, width);

            for (int j = 0; j <= height; j++)
            {
                if (world.sea(i, j) == 0)
                {
                    if (world.mountainheight(i, j) < maxmountainheight)
                    {
                        float crount = 0;
                        float jantotal = 0;
                        float jultotal = 0;

                        for (int k = i - 1; k <= i + 1; k++) // First, check the cells bordering the central one.
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (world.sea(kk, l) == 0 && world.mountainheight(kk, l) < maxmountainheight)
                                    {
                                        jantotal = jantotal + world.janrain(kk, l);
                                        jultotal = jultotal + world.julrain(kk, l);
                                        crount++;
                                    }
                                }
                            }
                        }

                        int jm = j - 1;
                        int jmm = j - 2;
                        int jp = j + 1;
                        int jpp = j + 2;

                        if (jmm >= 0 && jpp <= height) // Now, check the cells bordering those ones.
                        {
                            // Cells to the north

                            if (world.sea(im, jmm) == 0 && world.mountainheight(im, jmm) < maxmountainheight && (world.mountainheight(im, jm) < maxmountainheight || world.mountainheight(i, jm) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(im, jmm);
                                jultotal = jultotal + world.julrain(im, jmm);
                                crount++;
                            }

                            if (world.sea(i, jmm) == 0 && world.mountainheight(i, jmm) < maxmountainheight && (world.mountainheight(im, jm) < maxmountainheight || world.mountainheight(i, jm) < maxmountainheight || world.mountainheight(ip, jm) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(i, jmm);
                                jultotal = jultotal + world.julrain(i, jmm);
                                crount++;
                            }

                            if (world.sea(ip, jmm) == 0 && world.mountainheight(ip, jmm) < maxmountainheight && (world.mountainheight(i, jm) < maxmountainheight || world.mountainheight(ip, jm) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(ip, jmm);
                                jultotal = jultotal + world.julrain(ip, jmm);
                                crount++;
                            }

                            // Cells to the south

                            if (world.sea(im, jpp) == 0 && world.mountainheight(im, jpp) < maxmountainheight && (world.mountainheight(im, jp) < maxmountainheight || world.mountainheight(i, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(im, jpp);
                                jultotal = jultotal + world.julrain(im, jpp);
                                crount++;
                            }

                            if (world.sea(i, jpp) == 0 && world.mountainheight(i, jpp) < maxmountainheight && (world.mountainheight(im, jp) < maxmountainheight || world.mountainheight(i, jp) < maxmountainheight || world.mountainheight(ip, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(i, jpp);
                                jultotal = jultotal + world.julrain(i, jpp);
                                crount++;
                            }

                            if (world.sea(ip, jpp) == 0 && world.mountainheight(ip, jpp) < maxmountainheight && (world.mountainheight(i, jp) < maxmountainheight || world.mountainheight(ip, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(ip, jpp);
                                jultotal = jultotal + world.julrain(ip, jpp);
                                crount++;
                            }

                            // Cells to the west

                            if (world.sea(imm, jm) == 0 && world.mountainheight(imm, jm) < maxmountainheight && (world.mountainheight(im, jm) < maxmountainheight || world.mountainheight(im, j) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(imm, jm);
                                jultotal = jultotal + world.julrain(imm, jm);
                                crount++;
                            }

                            if (world.sea(imm, j) == 0 && world.mountainheight(imm, j) < maxmountainheight && (world.mountainheight(im, jm) < maxmountainheight || world.mountainheight(im, j) < maxmountainheight || world.mountainheight(im, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(imm, j);
                                jultotal = jultotal + world.julrain(imm, j);
                                crount++;
                            }

                            if (world.sea(imm, jp) == 0 && world.mountainheight(imm, jp) < maxmountainheight && (world.mountainheight(im, j) < maxmountainheight || world.mountainheight(im, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(imm, jp);
                                jultotal = jultotal + world.julrain(imm, jp);
                                crount++;
                            }

                            // Cells to the east

                            if (world.sea(ipp, jm) == 0 && world.mountainheight(ipp, jm) < maxmountainheight && (world.mountainheight(ip, jm) < maxmountainheight || world.mountainheight(ip, j) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(ipp, jm);
                                jultotal = jultotal + world.julrain(ipp, jm);
                                crount++;
                            }

                            if (world.sea(ipp, j) == 0 && world.mountainheight(ipp, j) < maxmountainheight && (world.mountainheight(ip, jm) < maxmountainheight || world.mountainheight(ip, j) < maxmountainheight || world.mountainheight(ip, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(ipp, j);
                                jultotal = jultotal + world.julrain(ipp, j);
                                crount++;
                            }

                            if (world.sea(ipp, jp) == 0 && world.mountainheight(ipp, jp) < maxmountainheight && (world.mountainheight(ip, j) < maxmountainheight || world.mountainheight(ip, jp) < maxmountainheight))
                            {
                                jantotal = jantotal + world.janrain(ipp, jp);
                                jultotal = jultotal + world.julrain(ipp, jp);
                                crount++;
                            }
                        }

                        if (crount > 0)
                        {
                            jantotal = jantotal / crount;
                            jultotal = jultotal / crount;

                            world.setjanrain(i, j, (int)jantotal);
                            world.setjulrain(i, j, (int)jultotal);
                        }
                    }
                }
            }
        }
    }

}

// This caps excessive rainfall.

void caprainfall(planet& world)
{
    int width = world.width();
    int height = world.height();
    float tilt = world.tilt();
    float eccentricity = world.eccentricity();

    float maxrain = 1000.0f; // Any rainfall over this will be greatly reduced.
    float capfactor = 0.1f; // Amount to multiply excessive rain by.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    float seasonalrain = static_cast<float>(world.seasonalrain(season, i, j));

                    if (seasonalrain > maxrain)
                        seasonalrain = (seasonalrain - maxrain) * capfactor + maxrain;

                    world.setseasonalrain(season, i, j, static_cast<int>(std::max(0.0f, seasonalrain)));
                }

                float janrain = (float)world.janrain(i, j);
                float julrain = (float)world.julrain(i, j);

                if (janrain > maxrain)
                {
                    janrain = janrain - maxrain;
                    janrain = janrain * capfactor;
                    janrain = janrain + maxrain;
                }

                if (julrain > maxrain)
                {
                    julrain = julrain - maxrain;
                    julrain = julrain * capfactor;
                    julrain = julrain + maxrain;
                }

                if (janrain < 0.0f)
                    janrain = 0.0f;

                if (julrain < 0.0f)
                    julrain = 0.0f;

                world.setjanrain(i, j, (int)janrain);
                world.setjulrain(i, j, (int)julrain);
                world.setseasonalrain(seasonjanuary, i, j, (int)janrain);
                world.setseasonalrain(seasonjuly, i, j, (int)julrain);
            }
        }
    });

    if (tilt < 10.f && eccentricity < 0.1f) // For worlds with low obliquity and low eccentricity, reduce any seasonal difference in rainfall.
    {
        float adjustfactor = tilt;

        float reducefactor = 0.5f + tilt / 20.0f;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    float averain = (float)((world.janrain(i, j) + world.julrain(i, j)) / 2);

                    if (averain > 0.0f)
                    {
                        float thisjanrain = (float)world.janrain(i, j) * adjustfactor + averain * (10.0f - adjustfactor);
                        float thisjulrain = (float)world.julrain(i, j) * adjustfactor + averain * (10.0f - adjustfactor);

                        thisjanrain = thisjanrain * reducefactor;
                        thisjulrain = thisjulrain * reducefactor;

                        world.setjanrain(i, j, (int)thisjanrain);
                        world.setjulrain(i, j, (int)thisjulrain);
                    }
                }
            }
        });
    }
}

void applycoastalclimates(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const float maxsearchdistance = tuning::climateresolution::scaleDistance(
        static_cast<float>(tuning::climate::maritime::maxSearchDistance), width, height);

    auto moderateseason = [maxsearchdistance](float temperature, float meantemperature, const maritimesample& sample, float thermalfactor)
    {
        const float influence = std::clamp(sample.influence, 0.0f, 1.0f);

        if (influence < tuning::climate::coastalclimate::minimumInfluence)
            return temperature;

        const float fetchfactor = 1.0f - std::clamp(
            static_cast<float>(sample.fetchdistance) / maxsearchdistance, 0.0f, 1.0f);
        const float moderation = std::clamp(
            influence * (tuning::climate::coastalclimate::rangeModerationFactor +
                fetchfactor * tuning::climate::coastalclimate::fetchModerationFactor),
            0.0f, 0.95f);
        float adjusted = meantemperature + (temperature - meantemperature) * (1.0f - moderation);
        adjusted = adjusted + sample.thermalanomaly * influence * thermalfactor;

        const float maxshift = tuning::climate::coastalclimate::maximumSeasonalShift;
        return temperature + std::clamp(adjusted - temperature, -maxshift, maxshift);
    };

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) == 1)
                    continue;

                const float januarytemp = static_cast<float>(world.jantemp(x, y));
                const float julytemp = static_cast<float>(world.jultemp(x, y));
                const float meantemperature = (januarytemp + julytemp) / 2.0f;
                const bool januaryiscold = januarytemp <= julytemp;
                const maritimesample januarymaritime = storedmaritimeinfluence(world, seasonjanuary, x, y);
                const maritimesample julymaritime = storedmaritimeinfluence(world, seasonjuly, x, y);

                const float januaryfactor = januaryiscold ?
                    tuning::climate::coastalclimate::coldSeasonThermalFactor :
                    tuning::climate::coastalclimate::warmSeasonThermalFactor;
                const float julyfactor = januaryiscold ?
                    tuning::climate::coastalclimate::warmSeasonThermalFactor :
                    tuning::climate::coastalclimate::coldSeasonThermalFactor;

                world.setjantemp(x, y, static_cast<int>(roundf(moderateseason(januarytemp, meantemperature, januarymaritime, januaryfactor))));
                world.setjultemp(x, y, static_cast<int>(roundf(moderateseason(julytemp, meantemperature, julymaritime, julyfactor))));
            }
        }
    });

    world.syncseasonalclimatefromlegacy();
}

// This smooths the temperatures in light of rainfall.

void smoothtemperatures(planet& world)
{
    int width = world.width();
    int height = world.height();

    vector<vector<int>> jantemp(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> jultemp(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 0)
                {
                    jantemp[i][j] = tempelevremove(world, world.jantemp(i, j), i, j);
                    jultemp[i][j] = tempelevremove(world, world.jultemp(i, j), i, j);
                }
            }
        }
    });

    for (int i = 0; i <= width; i++)
    {
        for (int j = 1; j < height; j++)
        {
            if (world.sea(i, j) == 0)
            {
                if (world.sea(i, j - 1) == 0 && world.sea(i, j + 1) == 0)
                {
                    if (jantemp[i][j] < jantemp[i][j - 1] && jantemp[i][j] < jantemp[i][j + 1])
                    {
                        jantemp[i][j] = (jantemp[i][j - 1] + jantemp[i][j + 1]) / 2;
                    }

                    if (jultemp[i][j] < jultemp[i][j - 1] && jultemp[i][j] < jultemp[i][j + 1])
                    {
                        jultemp[i][j] = (jultemp[i][j - 1] + jultemp[i][j + 1]) / 2;
                    }
                }
            }
        }
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            if (j == 0 || j == height)
                continue;

            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 0)
                {
                    world.setjantemp(i, j, tempelevadd(world, jantemp[i][j], i, j));
                    world.setjultemp(i, j, tempelevadd(world, jultemp[i][j], i, j));
                }
            }
        }
    });
}
}

// Rain map creator.

void createrainmap(planet& world, vector<vector<int>>& fractal,  int landtotal, int seatotal, boolshapetemplate smalllake[], boolshapetemplate shape[])
{
    int width = world.width();
    int height = world.height();
    
    int slopewaterreduce = 20; // The higher this is, the less extra rain falls on slopes.
    int maxmountainheight = 100;

    vector<vector<int>> inland(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    if (seatotal > 0)
    {
        populateinlanddistances(world, inland);

        if (beginworldgenstep("Advecting moisture and rainfall"))
            createadvectedrainfall(world, inland, fractal);

    }
    else if (landtotal > 0)
    {
        if (beginworldgenstep("Calculating rainfall"))
        {
            createdesertworldrain(world);

            for (int i = 0; i <= width; i++)
            {
                for (int j = 0; j <= height; j++)
                    inland[i][j] = 10000;
            }
        }
    }

    if (seatotal == 0)
    {
        if (beginworldgenstep("Calculating seasonal rainfall"))
            adjustseasonalrainfall(world, inland);

        if (beginworldgenstep("Smoothing rainfall"))
            smoothrainfall(world, maxmountainheight);
    }

    // The conservative hydrology already bounds water physically. Retain the
    // legacy presentation cap only for the fallback all-land rainfall path.
    if (seatotal == 0 && beginworldgenstep("Capping rainfall"))
        caprainfall(world);

    if (beginworldgenstep("Calculating coastal climate influence"))
        createcoastalclimateinfluence(world);

    if (beginworldgenstep("Applying coastal climate influence"))
        applycoastalclimates(world);

    // Land/ocean heat capacity and air-sea coupling now establish seasonal range.
    // Rainfall-derived warming/cooling and distance-to-coast amplification are no longer applied.

    if (beginworldgenstep("Smoothing temperatures"))
        smoothtemperatures(world);

    // Now we just sort out the mountain precipitation arrays, which will be used at the regional level for ensuring that higher mountain precipitation isn't splashed too far.

    if (beginworldgenstep("Calculating mountain rainfall"))
        createmountainprecipitation(world);
}

// Wind map creator.
void refreshadvectedrainfall(planet& world, vector<vector<int>>& fractal)
{
    vector<vector<int>> inland(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    populateinlanddistances(world, inland);
    createadvectedrainfall(world, inland, fractal);
    createmountainprecipitation(world);
}

// Temperature map creator.

void createtemperaturemap(planet& world, vector<vector<int>>& fractal)
{
    climateenergy::createSurfaceEnergyBalanceTemperatureMap(world, fractal);
}

// Work out the sea ice.

void createseaicemap(planet& world, vector<vector<int>>& fractal)
{
    // The coupled ocean already classified seasonal/permanent cover from ice
    // enthalpy. Liquid SST is pinned at freezing and cannot classify ice age.
    if (!tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
        return;
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int permice = -3;
    int seasice = -1;
    int maxadjust = 4;
    int adjustfactor = maxelev / maxadjust;
    int tempdist = 10;
    bool foundcoldwater = 0;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
                world.setseaice(i, j, 0);
        }
    });

    for (int i = 0; i <= width && foundcoldwater == 0; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) <= sealevel && coldoceansst(world, i, j) <= seasice)
            {
                foundcoldwater = 1;
                break;
            }
        }
    }

    if (foundcoldwater == 0)
        return;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.nom(i, j) > sealevel)
                    continue;

                int warmtotal = 0;
                int coldtotal = 0;
                int crount = 0;

                for (int k = -tempdist; k <= tempdist; k++)
                {
                    int kk = i + k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = -tempdist; l <= tempdist; l++)
                    {
                        if (k * k + l * l < tempdist * tempdist + tempdist)
                        {
                            int ll = j + l;

                            if (ll >= 0 && ll <= height && world.nom(kk, ll) <= sealevel)
                            {
                                warmtotal = warmtotal + warmoceansst(world, kk, ll);
                                coldtotal = coldtotal + coldoceansst(world, kk, ll);
                                crount++;
                            }
                        }
                    }
                }

                if (crount == 0)
                    continue;

                int warmsst = warmtotal / crount + (fractal[i][j] / adjustfactor) - maxadjust / 2;
                int coldsst = coldtotal / crount + (fractal[i][j] / adjustfactor) - maxadjust / 2;

                if (warmsst <= permice)
                    world.setseaice(i, j, 2);
                else if (coldsst <= seasice)
                    world.setseaice(i, j, 1);
            }
        }
    });

    // Remove odd bits

    for (int i = 0; i <= width; i++)
    {
        for (int j = 1; j < height; j++)
        {
            if (world.seaice(i, j) == 1)
            {
                if (world.seaice(i, j - 1) == 0 && world.seaice(i, j + 1) == 0)
                    world.setseaice(i, j, 0);
            }

            if (world.seaice(i, j) == 2)
            {
                if (world.seaice(i, j - 1) == 1 && world.seaice(i, j + 1) == 1)
                    world.setseaice(i, j, 2);
            }

        }
    }

    // Now we clean up the edge.

    for (int j = 0; j <= height; j++)
        world.setseaice(0, j, world.seaice(width, j));
}

// This creates mountain precipitation arrays (which are used at the regional level).

void createmountainprecipitation(planet& world)
{
    int width = world.width();
    int height = world.height();

    float totalmountaineffect = 1000; // Mountains higher than this will get the full effect.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.mountainheight(i, j) == 0 && world.craterrim(i, j) == 0)
            {
                world.setwintermountainraindir(i, j, 0);
                world.setsummermountainraindir(i, j, 0);
            }
        }
    }

    for (int n = 0; n < 1; n++)
    {
        // Using average neighbouring rainfall

        for (int i = 0; i <= width; i++) // Winter precipitation
        {
            for (int j = 0; j <= height; j++)
            {
                if ((world.mountainheight(i, j) != 0 || world.craterrim(i, j) != 0) && world.wintermountainraindir(i, j) == 0)
                {
                    // First, find the cell that the wind's coming from. It's the non-mountain cell with the highest precipitation.

                    int windfromcellx = -1;
                    int windfromcelly = -1;
                    int highestamount = -1;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.sea(kk, l) == 0 && world.mountainheight(kk, l) == 0 && world.craterrim(kk, l) == 0)
                                {
                                    if (world.winterrain(kk, l) > highestamount)
                                    {
                                        highestamount = world.winterrain(kk, l);
                                        windfromcellx = kk;
                                        windfromcelly = l;
                                    }
                                }
                            }
                        }
                    }

                    if (windfromcellx == -1) // We didn't find one!
                    {
                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (world.wintermountainraindir(kk, l) != 0)
                                    {
                                        if (world.winterrain(kk, l) > highestamount)
                                        {
                                            highestamount = world.winterrain(kk, l);
                                            windfromcellx = kk;
                                            windfromcelly = l;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Now find the average precipitation of neighbouring cells.

                    short crount = 0;
                    int total = 0;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.sea(kk, l) == 0 && world.mountainheight(kk, l) == 0 && world.craterrim(kk, l) == 0)
                                {
                                    crount++;
                                    total = total + world.winterrain(kk, l);
                                }
                            }
                        }
                    }

                    if (crount == 0) // We didn't find one!
                    {
                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (world.wintermountainraindir(kk, l) != 0)
                                    {
                                        crount++;
                                        total = total + world.winterrain(kk, l);
                                    }
                                }
                            }
                        }
                    }

                    if (windfromcellx != -1 && crount != 0) // If we found one!
                    {
                        float thisheight = (float)world.mountainheight(i, j);
                        float craterheight = (float)world.craterrim(i, j);

                        if (craterheight > thisheight)
                            thisheight = craterheight;

                        float effect = 1.0f;

                        if (thisheight < totalmountaineffect) // Reduce the strength of this.
                            effect = thisheight / totalmountaineffect;

                        float rainamount = (float)total / (float)crount;
                        float rainamount2 = (float)world.winterrain(i, j);

                        float newrainamount = (rainamount * effect) + (rainamount2 * (1.0f - effect));

                        world.setwintermountainrain(i, j, (int)newrainamount);

                        int dir = getdir(windfromcellx, windfromcelly, i, j);
                        world.setwintermountainraindir(i, j, dir);
                    }
                }
            }
        }

        for (int i = 0; i <= width; i++) // Summer precipitation
        {
            for (int j = 0; j <= height; j++)
            {
                if ((world.mountainheight(i, j) != 0 || world.craterrim(i, j) != 0) && world.summermountainraindir(i, j) == 0)
                {
                    // First, find the cell that the wind's coming from. It's the non-mountain cell with the highest precipitation.

                    int windfromcellx = -1;
                    int windfromcelly = -1;
                    int highestamount = -1;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.sea(kk, l) == 0 && world.mountainheight(kk, l) == 0 && world.craterrim(kk, l) == 0)
                                {
                                    if (world.summerrain(kk, l) > highestamount)
                                    {
                                        highestamount = world.summerrain(kk, l);
                                        windfromcellx = kk;
                                        windfromcelly = l;
                                    }
                                }
                            }
                        }
                    }

                    if (windfromcellx == -1) // We didn't find one!
                    {
                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (world.summermountainraindir(kk, l) != 0)
                                    {
                                        if (world.summerrain(kk, l) > highestamount)
                                        {
                                            highestamount = world.summerrain(kk, l);
                                            windfromcellx = kk;
                                            windfromcelly = l;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Now find the average precipitation of neighbouring cells.

                    short crount = 0;
                    int total = 0;

                    for (int k = i - 1; k <= i + 1; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - 1; l <= j + 1; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.sea(kk, l) == 0 && world.mountainheight(kk, l) == 0 && world.craterrim(kk, l) == 0)
                                {
                                    crount++;
                                    total = total + world.summerrain(kk, l);
                                }
                            }
                        }
                    }

                    if (crount == 0) // We didn't find one!
                    {
                        for (int k = i - 1; k <= i + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = j - 1; l <= j + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if (world.summermountainraindir(kk, l) != 0)
                                    {
                                        crount++;
                                        total = total + world.summerrain(kk, l);
                                    }
                                }
                            }
                        }
                    }

                    if (windfromcellx != -1 && crount != 0) // If we found one!
                    {
                        float thisheight = (float)world.mountainheight(i, j);
                        float craterheight = (float)world.craterrim(i, j);

                        if (craterheight > thisheight)
                            thisheight = craterheight;

                        float effect = 1.0f;

                        if (thisheight < totalmountaineffect) // Reduce the strength of this.
                            effect = thisheight / totalmountaineffect;

                        float rainamount = (float)total / (float)crount;
                        float rainamount2 = (float)world.summerrain(i, j);

                        float newrainamount = (rainamount * effect) + (rainamount2 * (1.0f - effect));

                        world.setsummermountainrain(i, j, (int)newrainamount);

                        int dir = getdir(windfromcellx, windfromcelly, i, j);
                        world.setsummermountainraindir(i, j, dir);
                    }
                }
            }
        }
    }
}

// This ensures that the climates at the bottom of the map are correct.

void checkpoleclimates(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    for (int i = 0; i <= width; i++)
    {
        world.setwinterrain(i, 0, world.winterrain(i, 1));
        world.setsummerrain(i, 0, world.summerrain(i, 1));
        world.setmaxtemp(i, 0, world.maxtemp(i, 1));
        world.setmintemp(i, 0, world.mintemp(i, 1));
        world.setwinterrain(i, height, world.winterrain(i, height - 1));
        world.setsummerrain(i, height, world.summerrain(i, height - 1));
        world.setmaxtemp(i, height, world.maxtemp(i, height - 1));
        world.setmintemp(i, height, world.mintemp(i, height - 1));
    }
}
