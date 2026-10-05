#include "simulations/terrain/import_elements.hpp"

#include "world/functions.hpp"

void generateimportedshelves(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    worldgrid<bool>& shelves = world.continentalshelves();
    shelves.fill(false);

    makecontinentalshelves(world, shelves, 4);
    createoceantrenches(world, shelves);

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> seafractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(seafractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    float coastalvarreduce = (float)maxelev / 3000.0f;
    float oceanvarreduce = (float)maxelev / 100.0f;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1 && shelves[i][j] == 1)
                world.setnom(i, j, sealevel - 100);
        }
    }
}

void generateimportedoceanridges(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    worldgrid<bool>& shelves = world.continentalshelves();
    shelves.fill(false);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) > sealevel - 400)
                shelves[i][j] = 1;
        }
    }

    createoceanridges(world, shelves);
}

void generateimportedseabed(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> seafractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(seafractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    int warpfactor = random(20, 80);
    warp(seafractal, width, height, maxelev, warpfactor, 1);

    float coastalvarreduce = (float)maxelev / 3000.0f;
    float oceanvarreduce = (float)maxelev / 1000.0f;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 1)
            {
                float var = (float)(seafractal[i][j] - maxelev / 2);
                var = var / coastalvarreduce;

                int newval = world.nom(i, j) + (int)var;

                if (newval > sealevel - 10)
                    newval = sealevel - 10;

                if (newval < 1)
                    newval = 1;

                world.setnom(i, j, newval);
            }
        }
    }

    // Smooth the seabed.

    smoothonlysea(world, 2);
}

void generateimportedlandelevation(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    // First, make a fractal map.

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

    int warpfactor = random(20, 80);
    warp(fractal, width, height, maxelev, warpfactor, 1);

    int fractaladd = sealevel - 2500;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            fractal[i][j] = fractal[i][j] + fractaladd;
    }

    // Now use it to change the land heights.

    fractaladdland(world, fractal);

    // Smooth the land.

    //smoothonlyland(world, 2);

    // Also, create extra elevation.

    createextraelev(world);
}

void generateimportedmountainchains(planet& world, boolshapetemplate landshape[], boolshapetemplate chainland[])
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    int baseheight = sealevel - 4500;
    if (baseheight < 1)
        baseheight = 1;
    int conheight = sealevel + 50;

    vector<vector<int>> plateaumap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    // First, make a fractal map.

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

    int fractaladd = sealevel - 2500;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            fractal[i][j] = fractal[i][j] + fractaladd;
    }

    twointegers dummy[1];

    createchains(world, baseheight, conheight, fractal, plateaumap, landshape, chainland, dummy, 0, 0, 5);
}

void generateimportedcraters(planet& world, std::vector<int>& squareroot)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    vector<vector<int>> oldterrain(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // Copy the terrain as it is before adding craters. This is so that we can add some variation to the craters afterwards, if there's sea on this world, so the depression filling doesn't fill them completely.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            oldterrain[i][j] = world.nom(i, j);
    }

    int cratertotal = random(500, 10000);

    createcratermap(world, cratertotal, squareroot, 1);

    float minseaproportion = (float)random(1, 500);
    minseaproportion = minseaproportion / 1000.0f;

    int minseasize = (int)(((float)width * (float)height) * minseaproportion);

    if (random(1, 4) != 1) // That may have produced seas inside craters, so probably remove those now.
        removesmallseas(world, minseasize, sealevel + 1);
    else // Definitely remove any really little bits of sea.
        removesmallseas(world, 20, sealevel + 1);

    int totalsea = 0;

    for (int i = 0; i <= -width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j))
                totalsea++;
        }
    }

    if (totalsea > 40) // If there's sea, add some variation to how prominent the craters are. This will create some gaps in the craters so they don't get entirely filled up by the depression filling.
    {
        vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

        int grain = 8; // Level of detail on this fractal map.
        float valuemod = 0.2f;
        int v = random(3, 6);
        float valuemod2 = (float)v;

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 1, 0);

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
            {
                float oldmult = (float)fractal[i][j] / (float)maxelev;
                float newmult = 1.0f - oldmult;

                float thiselev = (float)oldterrain[i][j] * oldmult + (float)world.nom(i, j) * newmult;

                world.setnom(i, j, (int)thiselev);
            }
        }
    }
}
