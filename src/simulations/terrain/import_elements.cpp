#include "simulations/terrain/import_elements.hpp"

// On the cell grid since T8 W2-T2: sphere noise (terrain::spherefractal) replaces the legacy
// fractals; the functions that still take legacy arrays get their values in the legacy slots
// (legacyslot, core/grid_shapes.hpp).

#include "core/deterministic_random.hpp"
#include "core/grid_shapes.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "world/functions.hpp"

using mesh::CellId;

namespace
{
// The legacy fractals' whole values.
vector<int> legacyvalues(const vector<float>& values)
{
    vector<int> whole(values.size());
    for (size_t c = 0; c < values.size(); c++)
        whole[c] = static_cast<int>(values[c]);
    return whole;
}

// Per-cell values in a legacy ARRAYWIDTH x ARRAYHEIGHT array, at the cells' legacy slots.
vector<vector<int>> legacyarray(const planet& world, const vector<int>& values)
{
    vector<vector<int>> array(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    const mesh::CellGrid& grid = world.grid();
    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        int x = 0, y = 0;
        legacyslot(grid, c, x, y);
        array[static_cast<size_t>(x)][static_cast<size_t>(y)] = values[c];
    }
    return array;
}
}

void generateimportedshelves(planet& world)
{
    int sealevel = world.sealevel();

    worldgrid<bool>& shelves = world.continentalshelves();
    shelves.fill(false);

    makecontinentalshelves(world, shelves, 4);
    createoceantrenches(world, shelves);

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
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    // A sphere noise for the legacy warped fractal (its warp only moved the field).
    const vector<float> seafractal = terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x5eab0001), terrain::terrainfractal(), (float)maxelev);

    float coastalvarreduce = (float)maxelev / 3000.0f;

    for (const CellId c : columnmajorcells(world.grid()))
    {
        if (world.sea(c) == 1)
        {
            float var = (float)((int)seafractal[c] - maxelev / 2);
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
    int sealevel = world.sealevel();

    int fractaladd = sealevel - 2500;

    // First, make a fractal map (a sphere noise for the legacy warped fractal).

    vector<int> values = legacyvalues(terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0x1a4d0001), terrain::terrainfractal(), 12750.0f));
    for (int& value : values)
        value += fractaladd;
    vector<vector<int>> fractal = legacyarray(world, values);

    // Now use it to change the land heights.

    fractaladdland(world, fractal);

    // Also, create extra elevation.

    createextraelev(world);
}

void generateimportedmountainchains(planet& world, boolshapetemplate landshape[], boolshapetemplate chainland[])
{
    int sealevel = world.sealevel();

    int baseheight = sealevel - 4500;
    if (baseheight < 1)
        baseheight = 1;
    int conheight = sealevel + 50;

    vector<vector<int>> plateaumap(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int fractaladd = sealevel - 2500;

    // First, make a fractal map.

    vector<int> values = legacyvalues(terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0xc4a10001), terrain::terrainfractal(), 12750.0f));
    for (int& value : values)
        value += fractaladd;
    vector<vector<int>> fractal = legacyarray(world, values);

    twointegers dummy[1];

    createchains(world, baseheight, conheight, fractal, plateaumap, landshape, chainland, dummy, 0, 0, 5);
}

void generateimportedcraters(planet& world, std::vector<int>& squareroot)
{
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    vector<int> oldterrain(world.cellcount(), 0); // Copy the terrain as it is before adding craters. This is so that we can add some variation to the craters afterwards, if there's sea on this world, so the depression filling doesn't fill them completely.

    for (CellId c = 0; c < world.cellcount(); c++)
        oldterrain[c] = world.nom(c);

    int cratertotal = random(500, 10000);

    createcratermap(world, cratertotal, squareroot, 1);

    float minseaproportion = (float)random(1, 500);
    minseaproportion = minseaproportion / 1000.0f;

    const float cells = (float)world.cellcount();
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
        const vector<int> fractal = legacyvalues(terrain::spherefractal(world, deterministiccontextseed(world.seed(), 0xc4a70001), terrain::terrainfractal(), (float)maxelev));

        for (CellId c = 0; c < world.cellcount(); c++)
        {
            float oldmult = (float)fractal[c] / (float)maxelev;
            float newmult = 1.0f - oldmult;

            float thiselev = (float)oldterrain[c] * oldmult + (float)world.nom(c) * newmult;

            world.setnom(c, (int)thiselev);
        }
    }
}
