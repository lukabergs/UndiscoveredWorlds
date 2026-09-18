#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "world_generation_debug.hpp"
#include "climate_generation.hpp"
#include "climate_coupling.hpp"
#include "climate_fields.hpp"
#include "climate_classification.hpp"
#include "drainage.hpp"
#include "lakes.hpp"
#include "deltas.hpp"
#include "wetlands.hpp"
#include "tides.hpp"
#include "climate_landforms.hpp"
#include "map_imports.hpp"

using namespace std;

namespace
{
void reseedglobalclimatepass(planet& world, int salt);
bool hasmatchingimportedclimatemapdimensions(const ImportedClimateMaps& importedclimate, const planet& world);
bool hasimportedtemperaturemap(const ImportedClimateMaps* importedclimate, const planet& world);
bool hasimportedprecipitationmap(const ImportedClimateMaps* importedclimate, const planet& world);
void synclegacyrainfallfromseasonal(planet& world);
void synclegacytemperaturesfromseasonal(planet& world);
void applyimportedtemperaturemap(planet& world, const ImportedClimateMaps& importedclimate);
void applyimportedprecipitationmap(planet& world, const ImportedClimateMaps& importedclimate);

void reseedglobalclimatepass(planet& world, int salt)
{
    fast_srand(deterministicfastseed(deterministiccontextseed(world.seed(), salt)));
}

bool hasmatchingimportedclimatemapdimensions(const ImportedClimateMaps& importedclimate, const planet& world)
{
    const int cellcount = (world.width() + 1) * (world.height() + 1);
    return importedclimate.width == world.width()
        && importedclimate.height == world.height()
        && static_cast<int>(importedclimate.annualTemperature.size()) == cellcount
        && static_cast<int>(importedclimate.annualPrecipitation.size()) == cellcount;
}

bool hasimportedtemperaturemap(const ImportedClimateMaps* importedclimate, const planet& world)
{
    return importedclimate != nullptr
        && importedclimate->hasTemperature
        && hasmatchingimportedclimatemapdimensions(*importedclimate, world);
}

bool hasimportedprecipitationmap(const ImportedClimateMaps* importedclimate, const planet& world)
{
    return importedclimate != nullptr
        && importedclimate->hasPrecipitation
        && hasmatchingimportedclimatemapdimensions(*importedclimate, world);
}

void synclegacyrainfallfromseasonal(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                world.setjanrain(x, y, world.seasonalrain(seasonjanuary, x, y));
                world.setjulrain(x, y, world.seasonalrain(seasonjuly, x, y));
            }
        }
    });
}

void synclegacytemperaturesfromseasonal(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                world.setjantemp(x, y, world.seasonaltemp(seasonjanuary, x, y));
                world.setjultemp(x, y, world.seasonaltemp(seasonjuly, x, y));
            }
        }
    });
}

void applyimportedtemperaturemap(planet& world, const ImportedClimateMaps& importedclimate)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const int index = y * (width + 1) + x;
                const float targetmean = static_cast<float>(importedclimate.annualTemperature[index]);
                float currentmean = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    currentmean += static_cast<float>(world.seasonaltemp(season, x, y));

                currentmean = currentmean / static_cast<float>(CLIMATESEASONCOUNT);
                const float delta = targetmean - currentmean;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    const float adjusted = static_cast<float>(world.seasonaltemp(season, x, y)) + delta;
                    world.setseasonaltemp(season, x, y, static_cast<int>(roundf(adjusted)));
                }
            }
        }
    });

    synclegacytemperaturesfromseasonal(world);
}

void applyimportedprecipitationmap(planet& world, const ImportedClimateMaps& importedclimate)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const int index = y * (width + 1) + x;
                const float targetmean = static_cast<float>(importedclimate.annualPrecipitation[index]);
                float currentmean = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    currentmean += static_cast<float>(world.seasonalrain(season, x, y));

                currentmean = currentmean / static_cast<float>(CLIMATESEASONCOUNT);

                if (currentmean > 0.0f)
                {
                    const float factor = targetmean / currentmean;

                    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    {
                        const float adjusted = max(0.0f, static_cast<float>(world.seasonalrain(season, x, y)) * factor);
                        world.setseasonalrain(season, x, y, static_cast<int>(roundf(adjusted)));
                    }
                }
                else
                {
                    const int fallbackrain = static_cast<int>(roundf(targetmean));

                    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                        world.setseasonalrain(season, x, y, fallbackrain);
                }
            }
        }
    });

    synclegacyrainfallfromseasonal(world);
}

// Establish a drainage level below the lowest land, consuming the same seeded draws.
int preparelandworldsealevel(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int seatotal = world.seatotal();
    int sealevel = world.sealevel();
    int raisesealevelchance = 6; // The higher this is, the *more* likely we are to try this.

    if (seatotal == 0 && random(1, raisesealevelchance) != 1)
    {
        int lowest = maxelev;

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                int thiselev = world.nom(i, j);

                if (thiselev < lowest)
                    lowest = thiselev;
            }
        }

        sealevel = lowest - random(5, 10);

        if (sealevel < 1)
            sealevel = 1;

        world.setsealevel(sealevel);
    }
    return sealevel;
}

// Value 1 excludes all lake tiles; value 2 excludes only lake centres.
void markcoastallakeexclusions(planet& world, vector<vector<int>>& nolake)
{
    const int width = world.width();
    const int height = world.height();
    const int minseadistance = 15;
    const int minseadistance2 = 8;
    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.outline(i, j) == 1)
            {
                for (int k = i - minseadistance2; k <= i + minseadistance2; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - minseadistance2; l <= j + minseadistance2; l++)
                    {
                        if (l >= 0 && l <= height)
                            nolake[kk][l] = 1;
                    }
                }

                for (int k = i - minseadistance; k <= i + minseadistance; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - minseadistance; l <= j + minseadistance; l++)
                    {
                        if (l >= 0 && l <= height && nolake[kk][l] == 0)
                            nolake[kk][l] = 2;
                    }
                }
            }
        }
    }
}

// Seed endorheic sinks before river planning on an all-land world.
bool seeddesertworldlakes(planet& world, int lakeattempts,
    vector<vector<vector<int>>>& saltlakemap, vector<vector<int>>& basins,
    vector<vector<int>>& nolake, boolshapetemplate smalllake[])
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int sealevel = world.sealevel();
    int highestelev = 0;
    int lowestelev = maxelev;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int thiselev = world.nom(i, j);

            if (thiselev < lowestelev)
                lowestelev = thiselev;

            if (thiselev > highestelev)
                highestelev = thiselev;
        }
    }

    int elevdiff = highestelev - lowestelev;

    int maxlakeelev = lowestelev + elevdiff / 4;

    vector<vector<int>> avoid(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This is for cells *not* to alter. This will remain empty (it's just here so we can use the depression-creating routine for normal lakes as well, where we don't want to mess with the path of the outflowing river.)

    for (int n = 0; n < lakeattempts; n++)
    {
        int x = random(1, width);
        int y = random(10, height - 10);

        for (int m = 0; m < 1000; m++)
        {
            if (world.nom(x, y) > maxlakeelev)
            {
                x = random(1, width);
                y = random(10, height - 10);
            }
            else
                m = 1000;
        }

        placesaltlake(world, x, y, 1, 0, saltlakemap, basins, avoid, nolake, smalllake);
    }

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) <= sealevel)
            {
                return true;
            }
        }
    }
    return false;
}

// Atmosphere supplies winds to the ocean; accepted SST then drives evaporation.
void generateinitialatmosphereandocean(planet& world)
{
    if (beginworldgenstep("Generating pressure map"))
    {
        reseedglobalclimatepass(world, 0x50023);
        createpressuremap(world);
    }

    if (beginworldgenstep("Generating vector wind map"))
    {
        reseedglobalclimatepass(world, 0x50024);
        createvectorwindmap(world);
        updatehorsebeltsfrompressure(world);
    }

    if (beginworldgenstep("Coupling wind-driven ocean and SST"))
    {
        reseedglobalclimatepass(world, 0x50021);
        createoceancurrentmap(world);
    }

    if (beginworldgenstep("Diagnosing ocean evaporation"))
    {
        reseedglobalclimatepass(world, 0x50022);
        createsurfacetemperaturemap(world);
    }
}
}

// This function creates the global climate.

void generateglobalclimate(planet& world, bool dorivers, bool dolakes, bool dodeltas, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate landshape[], vector<vector<int>>& mountaindrainage, vector<vector<bool>>& shelves, const ImportedClimateMaps* importedClimate)
{
    long seed = world.seed();
    fast_srand(seed);

    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();
    int seatotal = world.seatotal();
    int landtotal = world.landtotal();
    const bool importedtemperature = hasimportedtemperaturemap(importedClimate, world);
    const bool importedprecipitation = hasimportedprecipitationmap(importedClimate, world);

    // If there is no sea on a world, it may still have rain, provided we can find somewhere to put some salt lakes for the rivers to run into.

    int desertrainchance = 4; // On worlds with no sea, chance of trying to create rain anyway.
    int lakeattempts = random(1,8); // On worlds with no sea where there may be rain, make this many attempts to place a salt lake.

    int saltlakesplaced = 0;

    sealevel = preparelandworldsealevel(world);

    // If there is no sea but there is rain, we need to try to create some small bits of sea that will later become saltwater lakes, and fill depressions, to ensure that rivers run towards them.

    // First we need to prepare a no-lake template, marking out areas too close to the coasts, where lakes can't go. (We'll use this later whether or not the world has sea, so may as well do it now.)

    vector<vector<vector<int>>> saltlakemap(ARRAYWIDTH, vector<vector<int>>(ARRAYHEIGHT, vector<int>(2)));
    vector<vector<int>> nolake(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> basins(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // This will store where endorheic basins have already been carved.

    if (dolakes)
        markcoastallakeexclusions(world, nolake);

    bool desertworldrain = 0; // If this is 1, then this is a world with no sea but which will have rain.

    if (dorivers && seatotal == 0 && (importedprecipitation || random(1, desertrainchance) == 1))
        desertworldrain = seeddesertworldlakes(world, lakeattempts, saltlakemap, basins, nolake, smalllake);

    if (desertworldrain && beginworldgenstep("Filling depressions"))
    {
        depressionfill(world);

        addlandnoise(world); // Add a bit of noise, then do remove depressions again. This is to add variety to the river courses.

        depressionfill(world);
    }

    // Now, set the river land reduce factor.

    float riverlandreduce = 20.0f * world.gravity() * world.gravity();
    world.setriverlandreduce((int)riverlandreduce);

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = 0;
    float valuemod2 = 0.0f;
    int warpfactor = 0;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    if (beginworldgenstep("Generating global temperature map"))
    {
        reseedglobalclimatepass(world, 0x5002);

        // Start by generating a new fractal map.

        v = random(1, 4);
        valuemod2 = float(v);

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

        warpfactor = random(40, 80);
        warp(fractal, width, height, maxelev, warpfactor, 0);

        createtemperaturemap(world, fractal);

        if (importedtemperature)
            applyimportedtemperaturemap(world, *importedClimate);
    }

    if (seatotal > 0)
        generateinitialatmosphereandocean(world);

    // Now do the sea ice.

    if (seatotal > 0)
    {
        if (beginworldgenstep("Generating sea ice map"))
        {
            reseedglobalclimatepass(world, 0x5003);

            for (int i = 0; i <= width; i++)
            {
                for (int j = 0; j <= height; j++)
                    fractal[i][j] = 0;
            }

            createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

            warpfactor = random(40, 80);
            warp(fractal, width, height, maxelev, warpfactor, 0);

            createseaicemap(world, fractal);
        }

        if (beginworldgenstep("Calculating tides"))
            createtidalmap(world);

    }

    // Now do rainfall.

    if (seatotal > 0 || desertworldrain || importedprecipitation)
    {
        reseedglobalclimatepass(world, 0x5004);

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
                fractal[i][j] = 0;
        }

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

        warpfactor = random(40, 80);
        warp(fractal, width, height, maxelev, warpfactor, 0);
        warp(fractal, width, height, maxelev, warpfactor, 0);

        createrainmap(world, fractal, landtotal, seatotal, smalllake, landshape);

        if (seatotal > 0 &&
            tuning::climate::circulation::enableLaggedDiabaticCoupling)
        {
            convergeclimatecoupling(world, fractal);
        }
    }

    if (importedprecipitation)
        applyimportedprecipitationmap(world, *importedClimate);

    // Now add fjord mountains.

    if (seatotal > 0)
    {
        if (beginworldgenstep("Carving fjords"))
        {
            reseedglobalclimatepass(world, 0x5005);
            addfjordmountains(world);
        }
    }

    if (dorivers && (seatotal > 0 || desertworldrain || importedprecipitation))
    {
        // Now work out the rivers initially. We do this the first time so that after the first time we can place the salt lakes in appropriate places, and then we work out the rivers again.

        if (beginworldgenstep("Planning river courses"))
            createrivermap(world, mountaindrainage);

        if (dolakes || seatotal==0) // If there's no sea we need at least one salt lake.
        {
            // Now create salt lakes.

            if (beginworldgenstep("Placing hydrological basins"))
            {
                reseedglobalclimatepass(world, 0x5006);

                createsaltlakes(world, saltlakesplaced, saltlakemap, nolake, basins, smalllake);

                addlandnoise(world);
                depressionfill(world);

                for (int i = 0; i <= width; i++)
                {
                    for (int j = 0; j <= height; j++)
                    {
                        world.setriverdir(i, j, 0);
                        world.setriverjan(i, j, 0);
                        world.setriverjul(i, j, 0);
                    }
                }
            }

            if (beginworldgenstep("Generating rivers"))
                createrivermap(world, mountaindrainage);
        }

        // Now check river valleys in mountains.

        if (beginworldgenstep("Checking mountain river valleys"))
        {
            reseedglobalclimatepass(world, 0x5007);
            removerivermountains(world);
        }

        if (dolakes)
        {
            // Now create the lakes.

            if (beginworldgenstep("Generating lakes"))
            {
                reseedglobalclimatepass(world, 0x5008);

                convertsaltlakes(world, saltlakemap);

                createlakemap(world, nolake, smalllake, largelake);

                createriftlakemap(world, nolake);
            }
        }
        else
        {
            if (seatotal==0)
                convertsaltlakes(world, saltlakemap);
        }      
    }

    if (beginworldgenstep("Broadening FastLEM terrain from rivers"))
        broadenfastlemterrainfromrivers(world);

    world.setmaxriverflow();

    if (beginworldgenstep("Applying mountain temperature lapse"))
    {
        reseedglobalclimatepass(world, 0x5009);
        checkpoleclimates(world);
        world.syncseasonalclimatefromlegacy();

        if (importedtemperature)
            applyimportedtemperaturemap(world, *importedClimate);
    }

    if (beginworldgenstep("Calculating Koppen climates"))
        createclimatemap(world);

    if (beginworldgenstep("Calculating Holdridge biomes"))
        createbiomemap(world);

    // Now specials.

    if (beginworldgenstep("Generating sand dunes"))
    {
        reseedglobalclimatepass(world, 0x500a);
        createergs(world, smalllake, largelake, landshape);
    }

    if (beginworldgenstep("Generating salt pans"))
    {
        reseedglobalclimatepass(world, 0x500b);
        createsaltpans(world, smalllake, largelake);
    }

    // Add river deltas.

    if (dodeltas)
    {
        if (beginworldgenstep("Generating river deltas"))
        {
            reseedglobalclimatepass(world, 0x500c);
            createriverdeltas(world);
            checkrivers(world);
        }
    }

    // Now wetlands.

    if (beginworldgenstep("Generating wetlands"))
    {
        reseedglobalclimatepass(world, 0x500d);
        createwetlands(world, smalllake);
        removeexcesswetlands(world);
    }

    // Now it's time to finesse the roughness map.

    if (beginworldgenstep("Refining roughness map"))
    {
        reseedglobalclimatepass(world, 0x500e);
        refineroughnessmap(world);
    }

    // Check the rift lake map too.

    for (int i = 0; i < ARRAYWIDTH; i++)
    {
        for (int j = 0; j < ARRAYHEIGHT; j++)
        {
            if (world.lakestart(i, j) == 1 && world.riftlakesurface(i, j) == 0 && world.lakesurface(i, j) == 0)
                world.setlakestart(i, j, 0);
        }
    }

    // Check that the climates at the edges of the map are correct.

    if (beginworldgenstep("Checking poles"))
        checkpoleclimates(world);

    world.syncseasonalclimatefromlegacy();

    if (importedtemperature)
        applyimportedtemperaturemap(world, *importedClimate);

    if (importedprecipitation)
        applyimportedprecipitationmap(world, *importedClimate);

    if (importedprecipitation)
        createmountainprecipitation(world);

    createclimatemap(world);
    createbiomemap(world);
    exportclimatevalidationreport(world);

    if (dolakes == 0)
    {
        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                if (world.special(i, j) != 110)
                {
                    world.setlakesurface(i, j, 0);
                    world.setlakestart(i, j, 0);
                    world.setriftlakesurface(i, j, 0);
                    world.setriftlakebed(i, j, 0);
                }
            }
        }
    }

    removesealakes(world); // Also, make sure there are no weird bits of sea next to lakes.

    connectlakes(world); // Make sure lakes aren't fragmented.

    for (int i = 0; i <= width; i++) // Check erg/salt pans are the right depth.
    {
        for (int j = 0; j <= height; j++)
        {
            int special = world.special(i, j);
            
            if (special == 110 || special == 120)
            {
                int level = world.lakesurface(i, j);

                if (level <= sealevel)
                {
                    level = sealevel + 1;
                    world.setlakesurface(i, j, level);
                }

                world.setnom(i, j, level);
            }
        }
    }
}
