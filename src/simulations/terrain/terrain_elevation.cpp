// Elevation smoothing, relief noise, and depression filling.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;

namespace
{
bool checkdepressiontile(planet& world, std::vector<float>& filledmap, mesh::CellId c, int e, bool somethingdone, const std::vector<int>& noise);
}

// This function smoothes, but without turning any land to sea or vice versa.

void smoothland(planet& world, int amount)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    int grain = 8; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    int v = random(3, 6);
    float valuemod2 = (float)v;

    vector<vector<int>> fractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            fractal[i][j] = 0;
    }

    createfractal(fractal, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    // A sea cell's radius is set in km (tuning::terrain::ocean_floor) and chosen between the two
    // radii by the fractal across its own range; in cells it is the km over the cell size.
    int fractalmin = fractal[0][0];
    int fractalmax = fractal[0][0];

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            fractalmin = min(fractalmin, fractal[i][j]);
            fractalmax = max(fractalmax, fractal[i][j]);
        }
    }

    const float cellkm = world.radius() / 1000.0f * 3.14159265358979f / static_cast<float>(height + 1);
    const float minimumradiuskm = tuning::terrain::ocean_floor::smoothingMinimumRadiusKm;
    const float maximumradiuskm = tuning::terrain::ocean_floor::smoothingMaximumRadiusKm;
    const float fractalrange = static_cast<float>(max(1, fractalmax - fractalmin));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float crount = 0.0f;
            float ave = 0.0f;
            int origland = 0;
            int thisamount = amount;

            if (world.sea(i, j) == 0) // Check to see whether this point is originally land or sea
                origland = 1;
            else
            {
                const float position = static_cast<float>(fractal[i][j] - fractalmin) / fractalrange;
                const float radiuskm = minimumradiuskm + position * (maximumradiuskm - minimumradiuskm);
                thisamount = cellkm > 0.0f ? static_cast<int>(lroundf(radiuskm / cellkm)) : 0;
            }

            // Land is never written; a sea cell whose radius is under half a cell keeps its value.
            if (origland == 0 && thisamount < 1)
                continue;

            if (thisamount < 1)
                thisamount = 1;

            bool goahead = 1;
            int seacheck = 5;

            if (goahead == 1)
            {
                for (int k = i - thisamount; k <= i + thisamount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - thisamount; l <= j + thisamount; l++)
                    {
                        if (l > 0 && l < height)
                        {
                            ave = ave + world.nomfloat(kk, l);
                            crount++;
                        }
                    }
                }

                if (crount > 0.0f)
                {
                    ave = ave / crount;

                    if (ave > 0 && ave < maxelev)
                    {
                        int newland = 0;

                        if (wholemetres(ave) > sealevel) // Now check to see whether the new value is land or sea
                            newland = 1;

                        if (origland == 0 && newland == 0)
                            world.setnomfloat(i, j, ave);
                    }
                }
            }
        }
    }
}

// This function smoothes only the sea.

void smoothonlysea(planet& world, int amount)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float crount = 0.0f;
            float ave = 0.0f;
            int origland = 0;
            int thisamount = amount;

            if (world.basebelowsea(i, j))
            {
                for (int k = i - thisamount; k <= i + thisamount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - thisamount; l <= j + thisamount; l++)
                    {
                        if (l > 0 && l < height && world.basebelowsea(kk, l))
                        {
                            ave = ave + world.nomfloat(kk, l);
                            crount++;
                        }
                    }
                }

                if (crount > 0.0f)
                {
                    ave = ave / crount;

                    if (ave > 0.0f && ave < maxelev)
                        world.setnomfloat(i, j, ave);
                }
            }
        }
    }
}

// This function creates areas of raised elevation around where canyons might form later.

void createextraelev(planet& world)
{
    int width = world.width();
    int height = world.height();
    float gravity = world.gravity();

    vector<vector<int>> tempelev(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    float maxextra = 1000.0f; // Maximum amount that the extraelev array can have.

    float maxelev = (float)world.maxelevation();

    float div = maxelev / (maxextra * 2.0f);

    // First create a fractal for this map.

    int grain = 16; // Level of detail on this fractal map.
    float valuemod = 0.2f;
    float valuemod2 = 0.2f;

    createfractal(tempelev, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    int warpfactor = random(20, 80);
    warp(tempelev, width, height, (int)maxelev, warpfactor, 1);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                world.setextraelev(i, j, tempelev[i][j]);
                tempelev[i][j] = 0;
            }
        }
    });

    int v = random(3, 6);
    valuemod2 = (float)v;

    createfractal(tempelev, width, height, grain, valuemod, valuemod2, 1, (int)maxelev, 0, 0);

    warpfactor = random(20, 80);
    warp(tempelev, width, height, (int)maxelev, warpfactor, 1);

    // Now alter it a bit.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 1)
                    world.setextraelev(i, j, 0);
                else
                {
                    float e = world.extraelevfloat(i, j);
                    e = e / div;
                    e = e - maxextra;

                    if (e < 0.0f)
                        e = 0.0f;

                    if (e > maxextra)
                        e = maxextra;

                    e = e * gravity;

                    world.setextraelev(i, j, e);
                }
            }
        }
    });

    // That gives us a *lot* of extra elevation across the map. We want to apply it a bit more judiciously.
    // So we use a second fractal to mask the first, so that fewer areas get extra elevation.

    float multfactor = (float)random(1, 10);
    multfactor = multfactor / 10.0f + 1.0f; // Different worlds will have different amounts of extra elevation.

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 0)
                {
                    float e = (float)tempelev[i][j];
                    e = e / div;
                    e = e - maxextra;

                    e = e * multfactor; // The higher this is, the more widespread the extra elevation will be.

                    if (e < 0.0f)
                        e = 0.0f;

                    if (e < world.extraelevfloat(i, j))
                        world.setextraelev(i, j, e);
                }
            }
        }
    });

    // Now smooth it

    world.smoothextraelev(2);
    world.smoothextraelev(5);
}

namespace
{
// The cells of a grid as panels of columns x rows, id = panel x columns x rows + row x columns
// + column: one panel on lat-lon (W x W/2), the six faces on the cube (N x N).
struct ScanPanels
{
    int panels = 0;
    int columns = 0;
    int rows = 0;
};

ScanPanels scanpanels(const mesh::CellGrid& grid)
{
    if (grid.kind() == mesh::GridKind::cube)
        return { 6, grid.facecells(), grid.facecells() };
    return { 1, grid.columns(), grid.rows() };
}

// One of the eight directional scans of the fill: column by column (columnmajor) or row by row,
// columns and rows ascending or descending; the panels in turn, descending when both are.
// On lat-lon these are the legacy (x, y) scans; skipfirstcolumn keeps the seventh scan's
// "x > 0" (column 0 is left to the other seven).
vector<mesh::CellId> scanorder(const ScanPanels& panels, bool columnmajor, bool columnsdescending, bool rowsdescending,
                               bool skipfirstcolumn = false)
{
    vector<mesh::CellId> cells;
    cells.reserve(static_cast<size_t>(panels.panels) * panels.columns * panels.rows);
    const int firstcolumn = skipfirstcolumn ? 1 : 0;
    for (int p = 0; p < panels.panels; p++)
    {
        const int panel = columnsdescending && rowsdescending ? panels.panels - 1 - p : p;
        const mesh::CellId base = static_cast<mesh::CellId>(panel) * panels.columns * panels.rows;
        const int outercount = columnmajor ? panels.columns : panels.rows;
        const int innercount = columnmajor ? panels.rows : panels.columns;
        const bool outerdescending = columnmajor ? columnsdescending : rowsdescending;
        const bool innerdescending = columnmajor ? rowsdescending : columnsdescending;
        for (int o = 0; o < outercount; o++)
        {
            const int outer = outerdescending ? outercount - 1 - o : o;
            for (int n = 0; n < innercount; n++)
            {
                const int inner = innerdescending ? innercount - 1 - n : n;
                const int column = columnmajor ? outer : inner;
                const int row = columnmajor ? inner : outer;
                if (column < firstcolumn)
                    continue;
                cells.push_back(base + static_cast<mesh::CellId>(row * panels.columns + column));
            }
        }
    }
    return cells;
}

// A land cell with a sea cell among its four cardinal neighbours (planet::outline on lat-lon:
// a cell missing a cardinal neighbour, a polar row there, is not).
bool coastalcell(const planet& world, mesh::CellId c)
{
    if (world.sea(c) != 0)
        return false;
    const mesh::CellGrid& grid = world.grid();
    bool coastal = false;
    for (int slot = 1; slot <= mesh::slotcount; slot += 2)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);
        if (neighbour < 0)
            return false;
        if (world.sea(static_cast<mesh::CellId>(neighbour)) == 1)
            coastal = true;
    }
    return coastal;
}
}

// Fill land depressions using Planchon-Darboux directional scans. The noisy
// minimum gradient prevents flat areas from blocking later river routing.

void depressionfill(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = world.cellcount();
    int maxelev = world.maxelevation();

    int e = 2; // This is the extra amount added to ensure that everywhere slopes.

    vector<int> noise(count, 0); // This will contain noise that allows us to vary the value of e from tile to tile.

    for (const mesh::CellId c : columnmajorcells(grid))
        noise[c] = random(0, 5);

    vector<float> filledmap(count, 0.0f); // This will be the new version of the map.

    for (mesh::CellId c = 0; c < count; c++) // First, fill the new map to a huge height.
    {
        if (world.basebelowsea(c))
            filledmap[c] = world.nomfloat(c); // Sea tiles start with the normal heights.
        else
        {
            if (coastalcell(world, c)) // If this is a coastal tile
                filledmap[c] = world.nomfloat(c); // Coastal tiles start with the normal heights.
            else
                filledmap[c] = static_cast<float>(maxelev * 2); // Other tiles start very high
        }
    }

    // Eight scans, because there are eight possible ways to scan the map and we must rotate between them.
    const ScanPanels panels = scanpanels(grid);
    const vector<mesh::CellId> scans[8] =
    {
        scanorder(panels, true, false, false),
        scanorder(panels, false, true, true),
        scanorder(panels, true, true, false),
        scanorder(panels, false, true, false),
        scanorder(panels, true, false, true),
        scanorder(panels, false, false, true),
        scanorder(panels, true, true, true, true),
        scanorder(panels, false, false, false),
    };

    // Now we're ready to start!

    bool somethingdone = 0;

    do
    {
        for (int scan = 0; scan < 8; scan++)
        {
            if (scan > 0)
            {
                if (somethingdone == 0)
                    break;
                somethingdone = 0;
            }

            for (const mesh::CellId c : scans[scan])
            {
                if (filledmap[c] > world.nomfloat(c))
                    somethingdone = checkdepressiontile(world, filledmap, c, e, somethingdone, noise);
            }
        }
    } while (somethingdone == 1);

    for (mesh::CellId c = 0; c < count; c++)
        world.setnomfloat(c, filledmap[c]);
}

namespace
{
// This does the actual checking/filling of each tile: its eight neighbours from a random
// one, a cardinal step one metre steeper than a diagonal.

bool checkdepressiontile(planet& world, vector<float>& filledmap, mesh::CellId c, int e, bool somethingdone, const vector<int>& noise)
{
    const mesh::CellGrid& grid = world.grid();

    int start = random(0, 7);

    for (int n = start; n <= start + 7; n++)
    {
        int nn = wrap(n, 7);

        const std::int32_t neighbour = grid.neighbour(c, nn + 1);

        if (neighbour >= 0)
        {
            const mesh::CellId k = static_cast<mesh::CellId>(neighbour);
            int ee = 0;

            if ((nn + 1) % 2 == 1)
                ee = e + 1 + noise[c];

            else
                ee = e + noise[c];

            if (world.nomfloat(c) >= filledmap[k] + ee)
            {
                filledmap[c] = world.nomfloat(c);
                somethingdone = 1;

            }
            else
            {
                if (filledmap[c] > filledmap[k] + ee)
                {
                    filledmap[c] = filledmap[k] + ee;
                    somethingdone = 1;
                }
            }
        }
    }

    return (somethingdone);
}
}

// This function adds a bit of noise to the land.

void addlandnoise(planet& world)
{
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int maxadjust = 10;
    int adjustchance = 5; // The lower this is, the more there will be.

    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (!world.basebelowsea(c) && random(1, adjustchance) == 1)
        {
            int adjust = randomsign(random(1, maxadjust));

            float newamount = world.nomfloat(c) + static_cast<float>(adjust);

            if (newamount <= static_cast<float>(sealevel))
                newamount = static_cast<float>(sealevel);

            if (newamount > static_cast<float>(maxelev))
                newamount = static_cast<float>(maxelev);

            world.setnomfloat(c, newamount);
        }
    }
}
