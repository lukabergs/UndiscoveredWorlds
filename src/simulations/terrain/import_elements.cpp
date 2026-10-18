#include "simulations/terrain/import_elements.hpp"

// On the cell grid since T8 W2-T2: the lat-lon grid keeps the legacy fractals; every other grid
// takes sphere noise (terrain::coastgrid::spherefractal) and hands the functions that still take
// legacy arrays their values in the legacy slots (terrain::coastgrid::legacyslot).

#include "simulations/terrain/terrain_coasts.hpp"
#include "world/functions.hpp"

using mesh::CellId;
namespace cg = terrain::coastgrid;

namespace
{
bool latlon(const planet& world)
{
    return world.grid().kind() == mesh::GridKind::latlon;
}

// Per-cell values in a legacy ARRAYWIDTH x ARRAYHEIGHT array, at the cells' legacy slots.
vector<vector<int>> legacyarray(const planet& world, const vector<int>& values)
{
    vector<vector<int>> array(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    const mesh::CellGrid& grid = world.grid();
    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        int x = 0, y = 0;
        cg::legacyslot(grid, c, x, y);
        array[static_cast<size_t>(x)][static_cast<size_t>(y)] = values[c];
    }
    return array;
}
}

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

    if (latlon(world))
    {
        int grain = 8; // Level of detail on this fractal map.
        float valuemod = 0.2f;
        int v = random(3, 6);
        float valuemod2 = (float)v;

        vector<vector<int>> seafractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

        createfractal(seafractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0); // Not used; the call keeps the draws.
    }

    for (CellId c = 0; c < world.cellcount(); c++)
    {
        if (world.sea(c) == 1 && shelves[c] == 1)
            world.setnom(c, sealevel - 100);
    }
}

void generateimportedoceanridges(planet& world)
{
    int sealevel = world.sealevel();

    worldgrid<bool>& shelves = world.continentalshelves();
    shelves.fill(false);

    for (CellId c = 0; c < world.cellcount(); c++)
    {
        if (world.nom(c) > sealevel - 400)
            shelves[c] = 1;
    }

    createoceanridges(world, shelves);
}

void generateimportedseabed(planet& world)
{
    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    vector<int> seafractal;

    if (latlon(world))
    {
        int grain = 8; // Level of detail on this fractal map.
        float valuemod = 0.2f;
        int v = random(3, 6);
        float valuemod2 = (float)v;

        vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

        int warpfactor = random(20, 80);
        warp(fractal, width, height, maxelev, warpfactor, 1);

        seafractal.assign(world.cellcount(), 0);
        for (CellId c = 0; c < world.cellcount(); c++)
            seafractal[c] = fractal[static_cast<size_t>(world.grid().column(c))][static_cast<size_t>(world.grid().row(c))];
    }
    else
        seafractal = cg::spherefractal(world, 0x5eab0001, cg::grain8noise(), maxelev);

    float coastalvarreduce = (float)maxelev / 3000.0f;

    for (const CellId c : cg::columnorder(world.grid()))
    {
        if (world.sea(c) == 1)
        {
            float var = (float)(seafractal[c] - maxelev / 2);
            var = var / coastalvarreduce;

            int newval = world.nom(c) + (int)var;

            if (newval > sealevel - 10)
                newval = sealevel - 10;

            if (newval < 1)
                newval = 1;

            world.setnom(c, newval);
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

    int fractaladd = sealevel - 2500;

    vector<vector<int>> fractal;

    // First, make a fractal map.

    if (latlon(world))
    {
        int grain = 8; // Level of detail on this fractal map.
        float valuemod = 0.2f;
        int v = random(3, 6);
        float valuemod2 = (float)v;

        fractal.assign(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

        int warpfactor = random(20, 80);
        warp(fractal, width, height, maxelev, warpfactor, 1);

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
                fractal[i][j] = fractal[i][j] + fractaladd;
        }
    }
    else
    {
        vector<int> values = cg::spherefractal(world, 0x1a4d0001, cg::grain8noise(), 12750);
        for (int& value : values)
            value += fractaladd;
        fractal = legacyarray(world, values);
    }

    // Now use it to change the land heights.

    fractaladdland(world, fractal);

    // Also, create extra elevation.

    createextraelev(world);
}

void generateimportedmountainchains(planet& world, boolshapetemplate landshape[], boolshapetemplate chainland[])
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int baseheight = sealevel - 4500;
    if (baseheight < 1)
        baseheight = 1;
    int conheight = sealevel + 50;

    vector<vector<int>> plateaumap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int fractaladd = sealevel - 2500;

    vector<vector<int>> fractal;

    // First, make a fractal map.

    if (latlon(world))
    {
        int grain = 8; // Level of detail on this fractal map.
        float valuemod = 0.2f;
        int v = random(3, 6);
        float valuemod2 = (float)v;

        fractal.assign(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

        createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, 12750, 0, 0);

        for (int i = 0; i <= width; i++)
        {
            for (int j = 0; j <= height; j++)
                fractal[i][j] = fractal[i][j] + fractaladd;
        }
    }
    else
    {
        vector<int> values = cg::spherefractal(world, 0xc4a10001, cg::grain8noise(), 12750);
        for (int& value : values)
            value += fractaladd;
        fractal = legacyarray(world, values);
    }

    twointegers dummy[1];

    createchains(world, baseheight, conheight, fractal, plateaumap, landshape, chainland, dummy, 0, 0, 5);
}

void generateimportedcraters(planet& world, std::vector<int>& squareroot)
{
    const mesh::CellGrid& grid = world.grid();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    vector<int> oldterrain(world.cellcount(), 0); // Copy the terrain as it is before adding craters. This is so that we can add some variation to the craters afterwards, if there's sea on this world, so the depression filling doesn't fill them completely.

    for (CellId c = 0; c < world.cellcount(); c++)
        oldterrain[c] = world.nom(c);

    int cratertotal = random(500, 10000);

    createcratermap(world, cratertotal, squareroot, 1);

    float minseaproportion = (float)random(1, 500);
    minseaproportion = minseaproportion / 1000.0f;

    // The legacy count of cells (width x height, the largest indices) on lat-lon; the cells on the cube.
    const float cells = latlon(world) ? (float)world.width() * (float)world.height() : (float)world.cellcount();
    int minseasize = (int)(cells * minseaproportion);

    if (random(1, 4) != 1) // That may have produced seas inside craters, so probably remove those now.
        removesmallseas(world, minseasize, sealevel + 1);
    else // Definitely remove any really little bits of sea.
        removesmallseas(world, 20, sealevel + 1);

    // The legacy count of sea cells ran over no column (its loop ended at -width), so the
    // variation below never ran; kept as it was (T8 W2-T2 reports it).
    int totalsea = 0;

    if (totalsea > 40) // If there's sea, add some variation to how prominent the craters are. This will create some gaps in the craters so they don't get entirely filled up by the depression filling.
    {
        vector<int> fractal;

        if (latlon(world))
        {
            int grain = 8; // Level of detail on this fractal map.
            float valuemod = 0.2f;
            int v = random(3, 6);
            float valuemod2 = (float)v;

            vector<vector<int>> legacy(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
            createfractal(legacy, world.width(), world.height(), grain, valuemod, valuemod2, 1, maxelev, 1, 0);
            fractal.assign(world.cellcount(), 0);
            for (CellId c = 0; c < world.cellcount(); c++)
                fractal[c] = legacy[static_cast<size_t>(grid.column(c))][static_cast<size_t>(grid.row(c))];
        }
        else
            fractal = cg::spherefractal(world, 0xc4a70001, cg::grain8noise(), maxelev);

        for (CellId c = 0; c < world.cellcount(); c++)
        {
            float oldmult = (float)fractal[c] / (float)maxelev;
            float newmult = 1.0f - oldmult;

            float thiselev = (float)oldterrain[c] * oldmult + (float)world.nom(c) * newmult;

            world.setnom(c, (int)thiselev);
        }
    }
}
