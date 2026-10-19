// Template mountain chains and land shapes used by editing tools.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "simulations/terrain/terrain_coasts.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_shapes.hpp"
#include "simulations/terrain/terrain_volcanoes.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;

namespace
{
// Salt of the archipelagos' sphere noise (it replaces a fractal drawn from the stage's stream).
constexpr int archipelagosalt = 0x5304;

// The side of the mountain range templates' arrays (legacy pixels; the 2048-column world's rows).
constexpr int rangetemplatesize = ARRAYHEIGHT;

// The point east_px legacy pixels east and south_px south of p, through p's tangent frame.
mesh::Vec3 pixelpoint(const planet& world, const mesh::Vec3& p, double east_px, double south_px)
{
    if (east_px == 0.0 && south_px == 0.0)
        return p;
    const double pixelkm = terrain::legacycellkm(world);
    return mesh::fromtangent(mesh::tangentframe(p, world.grid().radius_m()), east_px * pixelkm, -south_px * pixelkm);
}

// A template's ridge direction (1 north, clockwise to 8 north-west) as a slot of the cell: the
// code itself on lat-lon, elsewhere the slot whose bearing is nearest the direction's.
int ridgeslot(const mesh::CellGrid& grid, mesh::CellId c, int dir)
{
    if (dir < 1 || dir > 8 || grid.kind() == mesh::GridKind::latlon)
        return dir;
    const double pi = 3.14159265358979323846;
    const double target = static_cast<double>(dir - 1) * pi / 4.0;
    int best = dir;
    double bestdifference = 10.0;
    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        if (grid.neighbour(c, slot) < 0)
            continue;
        double difference = fabs(grid.bearing_rad(c, slot) - target);
        difference = min(difference, 2.0 * pi - difference);
        if (difference < bestdifference)
        {
            bestdifference = difference;
            best = slot;
        }
    }
    return best;
}

void drawshape(planet& world, int shapenumber, const mesh::Vec3& centre, bool land, int baseheight, int conheight, boolshapetemplate shape[]);
void drawplateaushape(planet& world, int shapenumber, const mesh::Vec3& centre, int platheight, std::vector<int>& plateaumap, boolshapetemplate chainland[]);
}

// This function creates chains of mountains, with associated land where appropriate. On the
// cell grid (T8 W2-T3): positions are points on the sphere, moved in legacy pixels
// (terrain.sphere_noise.legacy_columns) through their tangent frame; a template pixel lies
// one legacy pixel from the next, so a cell wider than a pixel takes the first pixel of a range
// that reaches it.

void createchains(planet& world, int baseheight, int conheight, const vector<int>& fractal, vector<int>& plateaumap, boolshapetemplate landshape[], boolshapetemplate chainland[], const vector<mesh::CellId>& focuspoints, int focaldistance, int mode)
{
    // mode=0: continental chains (with associated land).
    // mode=1: smaller mountains (no associated land, can form peninsulas).
    // mode=2: smaller mountains (no associated land, can't form peninsulas).
    // mode=3: island chains.
    // mode=5: like mode 2, but fewer.
    // mode=6: like 3, but fewer.
    // mode=7: hills.

    (void)landshape;

    bool lower = 0;

    if (mode == 7)
    {
        mode = 2;
        lower = 1;
    }

    bool fewer = 0;

    if (mode == 5)
    {
        mode = 2;
        fewer = 1;
    }

    if (mode == 6)
    {
        mode = 3;
        fewer = 1;
    }

    const int focustotal = static_cast<int>(focuspoints.size());

    std::uint64_t chainseed = deterministiccontextseed(world.seed(), 0x4102);
    chainseed = deterministiccombine(chainseed, baseheight, conheight, focustotal, focaldistance, mode, lower ? 1 : 0, fewer ? 1 : 0);

    for (const mesh::CellId focus : focuspoints)
        chainseed = deterministiccombine(chainseed, static_cast<std::int64_t>(focus));

    fast_srand(deterministicfastseed(chainseed));

    int rangestartvar = 2; //6; // Possible distance between previous range and next one in a chain
    int changedirchance = 3; // Chance of changing direction (the lower it is, the more likely)
    int landchance = 10; // Chance of any given mountain pixel generating a nearby lump of land (the lower it is, the more likely)
    int minlanddist = 5;
    int maxlanddist = 10; // Range for how far lumps of land can be from the previous one (or the originating mountain)
    int stringmax = 20; //10; // Maximum number of lumps in a string of land
    int landvar = 4; //8; // Maximum distance a lump of land can be displaced from its proper position
    int platchance = 4; // Chance of any range generating a plateau
    int platmaxheight = 30; // Percentage of the range height for the plateau maximum height
    int platminheight = 5; // Percentage of the range height for the plateau minimum height
    int minplatreduce = 80;
    int maxplatreduce = 95; // Range of percentages to reduce the height of the plateau with each step
    int mode1landmin = 30; // Minimum distance to land from the starting point of a mode 1 chain.

    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);

    vector<vector<unsigned char>> rangeheighttemplate(rangetemplatesize, vector<unsigned char>(rangetemplatesize, 0));
    vector<vector<unsigned char>> rangeridgetemplate(rangetemplatesize, vector<unsigned char>(rangetemplatesize, 0)); // These will hold the range we're currently putting on.

    auto intemplate = [&](int i, int j) { return i >= 0 && i < rangetemplatesize && j >= 0 && j < rangetemplatesize; };

    vector<int> rangestamp(world.cellcount(), 0); // The range that set a cell's mountain (0: none).
    int rangeid = 0;

    int chainno = 0;
    int minchainlength = 0;
    int maxchainlength = 0;

    switch (mode)
    {
    case 0:

        chainno = random(10, 40); //(5,20);
        maxchainlength = 15;
        minchainlength = 4;
        break;

    case 1:

        chainno = random(8, 60); //(2,15);
        maxchainlength = 8;
        minchainlength = 1;
        break;

    case 2:

        chainno = random(4000, 32000); //(1000,8000);

        if (fewer == 1)
            chainno = 1000;

        maxchainlength = 8;
        minchainlength = 1;
        break;

    case 3:

        chainno = random(20, 200); //(80,800);

        if (fewer == 1)
            chainno = 10;

        maxchainlength = 6;
        minchainlength = 1;
        break;
    }

    // Now do the chains.

    for (int chain = 1; chain <= chainno; chain++)
    {
        mesh::Vec3 oldp = grid.centre(terrain::randomcellbyarea(world)); // The point of the current blob.

        if ((mode == 0 || mode == 3) && focustotal > 0) // Continental ranges should begin near the continental focus points, if possible. Islands should begin near their focus points too.
        {
            int thisfocus = random(0, focustotal - 1);

            float xdiff = (float)random(1, 100);
            float ydiff = (float)random(1, 100);

            xdiff = xdiff / 100.0f;
            ydiff = ydiff / 100.0f;

            xdiff = xdiff * xdiff;
            ydiff = ydiff * ydiff;

            xdiff = xdiff * (float)focaldistance;
            ydiff = ydiff * (float)focaldistance;

            const int movex = (int)randomsign(xdiff);
            const int movey = (int)randomsign(ydiff);

            oldp = pixelpoint(world, grid.centre(focuspoints[static_cast<size_t>(thisfocus)]), movex, movey);
        }

        int chaindir = random(1, 8); // Chaindir is the direction of the chain, from 1-8 (sort of on the diagonals!).
        int maxranges = random(minchainlength, maxchainlength); // The chain won't have more ranges in it than this.
        int rangeno = 1;
        bool goahead = 1;

        int tempno = -1;

        if (!(mode == 2 && world.sea(grid.locate(oldp)) == 1)) // If we're placing land-only mountains and we're in the sea, just skip the whole thing.
        {
            do
            {
                int lasttemp = tempno;

                do // Get a new template, making sure it's not the same as the last one.
                {
                    tempno = random(1, MOUNTAINTEMPLATESTOTAL);
                } while (tempno == lasttemp);

                int dampner = 1800; //maxelev/7; // This is to reduce the maximum height a bit.
                int rangeheight = 0;

                switch (mode)
                {
                case 0:

                    rangeheight = random((maxelev - sealevel) / 2, maxelev - sealevel); // Approximate height of the range (above sea level)
                    rangeheight = rangeheight - dampner;

                    break;

                case 1:

                    rangeheight = random(200, (maxelev - sealevel) / 2);
                    rangeheight = rangeheight - dampner;

                    break;

                case 2:

                    rangeheight = random(100, (maxelev - sealevel) / 3);
                    rangeheight = rangeheight - dampner;

                    if (lower == 1)
                        rangeheight = rangeheight / 4;

                    break;

                case 3:

                    rangeheight = random(sealevel / 2, sealevel + 3000);

                    break;
                }

                if (rangeheight < 200)
                    rangeheight = 200;

                bool plateau = 0;
                int platheight = 0;

                if (mode == 0 && random(1, platchance) == 1) // See if we might do a plateau.
                {
                    plateau = 1;

                    float platheightfactor = (float)random(platminheight, platmaxheight);
                    platheightfactor = platheightfactor / 100.0f;

                    platheight = (int)((float)rangeheight * platheightfactor);
                }

                int rangeinc = rangeheight / 10; // Difference between the different template heights
                int rangelittleinc = rangeinc / 100; // Amount to vary the template heights by

                // Get the template.

                int thisdir = chaindir;

                if (random(1, 2) == 1)
                    thisdir = thisdir + 4;

                if (thisdir > 8)
                    thisdir = thisdir - 8;

                int span = createmountainrangetemplate(rangeridgetemplate, rangeheighttemplate, tempno, thisdir);

                // Now we need to find the start of the actual range within the template.

                int startx = 0;
                int starty = 0;

                if (chaindir == 1 || chaindir == 2)
                {
                    for (int i = 0; i <= span; i++)
                    {
                        for (int j = 0; j <= span; j++)
                        {
                            if (intemplate(i, j) && rangeridgetemplate[i][j] != 0)
                            {
                                startx = i;
                                starty = j;

                                i = span;
                                j = span;
                            }
                        }
                    }
                }

                if (chaindir == 3 || chaindir == 4)
                {
                    for (int j = 0; j <= span; j++)
                    {
                        for (int i = 0; i <= span; i++)
                        {
                            if (intemplate(i, j) && rangeridgetemplate[i][j] != 0)
                            {
                                startx = i;
                                starty = j;

                                i = span;
                                j = span;
                            }
                        }
                    }
                }

                if (chaindir == 5 || chaindir == 6)
                {
                    for (int i = span; i >= 0; i--)
                    {
                        for (int j = 0; j <= span; j++)
                        {
                            if (intemplate(i, j) && rangeridgetemplate[i][j] != 0)
                            {
                                startx = i;
                                starty = j;

                                i = 0;
                                j = span;
                            }
                        }
                    }
                }

                if (chaindir == 7 || chaindir == 8)
                {
                    for (int j = span; j >= 0; j--)
                    {
                        for (int i = 0; i <= span; i++)
                        {
                            if (intemplate(i, j) && rangeridgetemplate[i][j] != 0)
                            {
                                startx = i;
                                starty = j;

                                i = span;
                                j = 0;
                            }
                        }
                    }
                }

                int movex = 0;
                int movey = 0;

                if (chaindir == 1 || chaindir == 2)
                {
                    movex = random(0, rangestartvar);
                    movey = randomsign(random(0, rangestartvar));
                }

                if (chaindir == 3 || chaindir == 4)
                {
                    movex = randomsign(random(0, rangestartvar));
                    movey = random(0, rangestartvar);
                }

                if (chaindir == 5 || chaindir == 6)
                {
                    movex = 0 - random(0, rangestartvar);
                    movey = randomsign(random(0, rangestartvar));
                }

                if (chaindir == 7 || chaindir == 8)
                {
                    movex = randomsign(random(0, rangestartvar));
                    movey = 0 - random(0, rangestartvar);
                }

                oldp = pixelpoint(world, oldp, movex, movey);

                // oldp is now the point where we want the new range to start: template pixel
                // (i, j) lies i - startx pixels east and j - starty south of it.

                const mesh::CellId oldcell = grid.locate(oldp);
                const mesh::TangentFrame frame = mesh::tangentframe(oldp, grid.radius_m());

                auto pixelofrange = [&](int i, int j)
                {
                    return mesh::fromtangent(frame, static_cast<double>(i - startx) * pixelkm, -static_cast<double>(j - starty) * pixelkm);
                };

                if (mode == 1 && world.sea(oldcell) == 1) // Peninsular mountains have to start within a certain distance of land.
                {
                    bool foundone = 0;

                    for (const mesh::CellId n : terrain::cellswithinkm(world, oldcell, (mode1landmin + 0.5) * pixelkm))
                    {
                        if (world.sea(n) == 0)
                        {
                            foundone = 1;
                            break;
                        }
                    }

                    if (foundone == 0)
                        goahead = 0;
                }

                if (mode == 2) // Non-continental, non-peninsular mountains can't start in the sea.
                {
                    if (world.sea(oldcell) == 1)
                        goahead = 0;
                }

                if (mode == 3) // Island chains can't start on the land.
                {
                    if (world.sea(oldcell) == 0)
                        goahead = 0;
                }

                // Now we need to see whether there is a clear area to paste it.

                if (mode != 2)
                {
                    for (int i = 0; i <= span && goahead == 1; i++)
                    {
                        for (int j = 0; j <= span && goahead == 1; j++)
                        {
                            const mesh::CellId c = grid.locate(pixelofrange(i, j));

                            if (c != mesh::nocell && world.mountainheightfloat(c) != 0)
                                goahead = 0;
                        }
                    }
                }

                int hislandadd = random(1000, 5000);

                if (goahead == 1) // If the area is clear
                {
                    rangeid++;

                    for (int i = 0; i <= span; i++)
                    {
                        for (int j = 0; j <= span; j++)
                        {
                            const bool templatepixel = intemplate(i, j) && rangeheighttemplate[i][j] != 0;

                            if (templatepixel == false && goahead == 0) // Nothing more to learn from this pixel.
                                continue;

                            const mesh::Vec3 p = pixelofrange(i, j);
                            const mesh::CellId c = grid.locate(p);

                            if (c == mesh::nocell || rangestamp[c] == rangeid) // This range has already set this cell.
                                continue;

                            if (world.mountainheightfloat(c) == 0)
                            {
                                if (templatepixel)
                                {
                                    int h = rangeheighttemplate[i][j] * rangeinc;
                                    h = h + randomsign(rangelittleinc * random(0, 100));

                                    if (mode == 3)
                                        h = h + hislandadd;

                                    float frac = (float)fractal[c];
                                    frac = frac / (float)maxelev;
                                    const float scaled = (float)h * frac;
                                    const int ridge = ridgeslot(grid, c, rangeridgetemplate[i][j]);

                                    if (mode == 0 || mode == 1)
                                    {
                                        world.setmountainheight(c, scaled);
                                        world.setmountainridge(c, ridge);
                                        world.setnomfloat(c, conheight);
                                    }

                                    if (mode == 2) // Mountains of this kind cannot create their own land, so must stop at the sea.
                                    {
                                        if (world.mapfloat(c) >= conheight)
                                        {
                                            world.setmountainheight(c, scaled);
                                            world.setmountainridge(c, ridge);

                                            if (fewer == 0)
                                                world.setnomfloat(c, conheight);
                                        }
                                    }

                                    if (mode == 3) // Islands are just strange and different.
                                    {
                                        if (scaled > (float)sealevel)
                                        {
                                            world.setmountainheight(c, scaled - (float)sealevel);
                                            world.setmountainridge(c, ridge);
                                            world.setnomfloat(c, conheight);
                                        }
                                    }

                                    if (world.mountainheightfloat(c) != 0)
                                        rangestamp[c] = rangeid;

                                    if (mode == 0) // Only create substantial land if these are continental mountains.
                                    {
                                        if (random(1, landchance) == 1) // See if we'll paste a chunk of land nearby
                                        {
                                            int dist = random(minlanddist, maxlanddist);

                                            int landmovex = 0;
                                            int landmovey = 0;

                                            switch (chaindir)
                                            {
                                            case 1:
                                                landmovex = 0 - dist / 2;
                                                landmovey = 0 - dist;
                                                break;

                                            case 2:
                                                landmovex = dist / 2;
                                                landmovey = 0 - dist;
                                                break;

                                            case 3:
                                                landmovex = dist;
                                                landmovey = 0 - dist / 2;
                                                break;

                                            case 4:
                                                landmovex = dist;
                                                landmovey = dist / 2;
                                                break;

                                            case 5:
                                                landmovex = dist / 2;
                                                landmovey = dist;
                                                break;

                                            case 6:
                                                landmovex = 0 - dist / 2;
                                                landmovey = dist;
                                                break;

                                            case 7:
                                                landmovex = 0 - dist;
                                                landmovey = dist / 2;
                                                break;

                                            case 8:
                                                landmovex = 0 - dist;
                                                landmovey = 0 - dist / 2;
                                                break;
                                            }

                                            mesh::Vec3 landp = p;

                                            int max = rangeno;
                                            int max2 = maxranges - rangeno;

                                            if (max2 < max)
                                                max = max2;

                                            int mult = 4;

                                            max = max * mult;

                                            if (max > stringmax)
                                                max = stringmax;

                                            int totallumps = random((max / 4) * 3, max);

                                            for (int lump = 1; lump <= totallumps; lump++) // Now we will paste a string of land lumps leading away from the mountain.
                                            {
                                                const int jitterx = randomsign(random(0, landvar));
                                                const int jittery = randomsign(random(0, landvar));

                                                landp = pixelpoint(world, landp, landmovex + jitterx, landmovey + jittery);

                                                int shapenumber = random(0, 1);
                                                drawshape(world, shapenumber, landp, 1, baseheight, conheight, chainland);

                                                if (plateau == 1) // Draw a bit of plateau, and reduce the height for the next pass.
                                                {
                                                    drawplateaushape(world, shapenumber, landp, platheight, plateaumap, chainland);

                                                    float platreduce = (float)random(minplatreduce, maxplatreduce);
                                                    platreduce = platreduce / 100.0f;
                                                    platheight = (int)((float)platheight * platreduce);
                                                }
                                            }
                                        }
                                    }

                                    if (mode == 1) // Mountains of this kind can create land just immediately around themselves.
                                    {
                                        if (world.nomfloat(c) < conheight)
                                            world.setnomfloat(c, conheight);
                                    }

                                    if (mode == 3) // Island mountains can create land just immediately around themselves.
                                    {
                                        if (h > sealevel)
                                        {
                                            if (world.nomfloat(c) < conheight)
                                                world.setnomfloat(c, conheight);
                                        }
                                    }
                                }
                            }
                            else // If this range runs into other mountains, end the chain.
                                goahead = 0;
                        }
                    }

                    // Now we need to find the new starting point.

                    int endx = -1;
                    int endy = -1;

                    if (chaindir == 1 || chaindir == 2)
                    {
                        for (int i = span; i >= 0 && endx < 0; i--)
                        {
                            for (int j = 0; j <= span; j++)
                            {
                                if (intemplate(i, j) && rangeheighttemplate[i][j] != 0)
                                {
                                    endx = i;
                                    endy = j;
                                    break;
                                }
                            }
                        }
                    }

                    if (chaindir == 3 || chaindir == 4)
                    {
                        for (int j = span; j >= 0 && endx < 0; j--)
                        {
                            for (int i = 0; i <= span; i++)
                            {
                                if (intemplate(i, j) && rangeheighttemplate[i][j] != 0)
                                {
                                    endx = i;
                                    endy = j;
                                    break;
                                }
                            }
                        }
                    }

                    if (chaindir == 5 || chaindir == 6)
                    {
                        for (int i = 0; i <= span && endx < 0; i++)
                        {
                            for (int j = 0; j <= span; j++)
                            {
                                if (intemplate(i, j) && rangeheighttemplate[i][j] != 0)
                                {
                                    endx = i;
                                    endy = j;
                                    break;
                                }
                            }
                        }
                    }

                    if (chaindir == 7 || chaindir == 8)
                    {
                        for (int j = 0; j <= span && endx < 0; j++)
                        {
                            for (int i = 0; i <= span; i++)
                            {
                                if (intemplate(i, j) && rangeheighttemplate[i][j] != 0)
                                {
                                    endx = i;
                                    endy = j;
                                    break;
                                }
                            }
                        }
                    }

                    if (endx >= 0)
                        oldp = pixelofrange(endx, endy);

                    // Now change direction, possibly.

                    if (random(1, changedirchance) == 1)
                    {
                        chaindir = chaindir + randomsign(1);

                        if (chaindir > 8)
                            chaindir = 1;

                        if (chaindir < 1)
                            chaindir = 8;
                    }
                }

                rangeno++;

                if (rangeno > maxranges)
                    goahead = 0;

            } while (goahead == 1);
        }
    }
}

namespace
{
// Draws a shape from an image onto the map, centred at a point, one legacy pixel per pixel
// (terrain::coastgrid::stamptemplate), mirrored at random.

void drawshape(planet& world, int shapenumber, const mesh::Vec3& centre, bool land, int baseheight, int conheight, boolshapetemplate shape[])
{
    const bool leftr = random(0, 1) == 1; // If it's 1 then we reverse it left-right
    const bool downr = random(0, 1) == 1; // If it's 1 then we reverse it top-bottom

    for (const mesh::CellId c : terrain::coastgrid::stamptemplate(world, centre, shape[shapenumber], leftr, downr))
    {
        if (land == 1)
            world.setnomfloat(c, conheight);
        else
        {
            world.setnomfloat(c, baseheight);
            world.setmountainheight(c, 0);
            world.setmountainridge(c, 0);
        }
    }
}
}

namespace
{
// Draws a shape from an image onto the plateau map. As in the legacy function, every set pixel
// of the shape writes the one cell a pixel east and south of the shape's top-left corner (its
// mirror draws are made but do not matter).

void drawplateaushape(planet& world, int shapenumber, const mesh::Vec3& centre, int platheight, vector<int>& plateaumap, boolshapetemplate chainland[])
{
    const boolshapetemplate& shape = chainland[shapenumber];

    int imheight = shape.ysize() - 1;
    int imwidth = shape.xsize() - 1;

    (void)random(0, 1); // The legacy left-right and top-bottom mirror draws.
    (void)random(0, 1);

    bool anypixel = false;

    for (int i = 0; i <= imwidth && !anypixel; i++)
    {
        for (int j = 0; j <= imheight && !anypixel; j++)
            anypixel = shape.point(i, j) == 1;
    }

    if (!anypixel)
        return;

    const mesh::CellId c = world.grid().locate(pixelpoint(world, centre, 1 - imwidth / 2, 1 - imheight / 2));

    if (c != mesh::nocell && world.sea(c) == 0)
        plateaumap[c] = platheight;
}
}

// This function does something similar, but adds the fractal map onto the land map.

void fractaladdland(planet& world, const vector<int>& fractal, const vector<int>& divfield)
{
    int div = 300; //150; // The higher this number, the more highlands there will be.
    int maxamount = 25; // The higher this is, the flatter the flat areas will be.
    int minamount = 10; // The lower this is, the higher the high areas will be.

    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    // divfield (a second field of the fractal's values: the legacy flipped and shifted copy of
    // it) varies the divisor, to create highlands and lowlands.

    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (!world.basebelowsea(c))
        {
            int valtoadd = fractal[c] - sealevel;

            int divamount = divfield[c];
            divamount = divamount / div;
            if (divamount > maxamount)
                divamount = maxamount;
            if (divamount < minamount)
                divamount = minamount;

            valtoadd = valtoadd / divamount;

            float newval = world.nomfloat(c) + static_cast<float>(valtoadd);

            if (newval > static_cast<float>(maxelev))
                newval = static_cast<float>(maxelev);

            if (newval > world.nomfloat(c))
                world.setnomfloat(c, newval);
        }
    }
}

// This scatters islands across areas where land has been removed.

void terrain::makearchipelagos(planet& world, const vector<std::uint8_t>& removedland, boolshapetemplate landshape[])
{
    int maxelev = world.maxelevation();

    // Two sphere noises for the islands (the legacy grain-8 fractal, and its value half a turn
    // of longitude away, which acted as a second, independent one).
    const std::uint64_t seed = deterministiccontextseed(world.seed(), archipelagosalt);
    const vector<float> islandfractal = terrain::spherefractal(world, seed, terrain::archipelagofractal(), static_cast<float>(maxelev));
    const vector<float> otherfractal =
        terrain::spherefractal(world, deterministiccombine(seed, 1), terrain::archipelagofractal(), static_cast<float>(maxelev));

    const vector<mesh::CellId> cells = columnmajorcells(world.grid());

    for (const mesh::CellId c : cells)
    {
        int thisrandom = random(1, maxelev * 100);
        (void)thisrandom;

        int islandchance = 0;

        if (islandfractal[c] < maxelev / 5)
            islandchance = 70;

        if (islandfractal[c] < maxelev / 6)
            islandchance = 40;

        if (islandfractal[c] < maxelev / 7)
            islandchance = 30;

        if (islandchance != 0 && otherfractal[c] > maxelev / 3 && removedland[c] == 1 && world.sea(c) == 1 && random(1, islandchance) == 1) // Put an island here.
        {
            int peakheight = random(500, 3000);
            terrain::makemountainisland(world, c, static_cast<float>(peakheight));
        }
    }

    // No mountain islands within two legacy pixels of existing land (at least not marked as such).
    const double radius_km = 2.5 * terrain::legacycellkm(world);

    for (const mesh::CellId c : cells)
    {
        if (world.mountainisland(c) == 1)
        {
            bool nearland = 0;

            for (const mesh::CellId n : terrain::cellswithinkm(world, c, radius_km))
            {
                if (world.sea(n) == 0 && world.mountainisland(n) == 0)
                {
                    nearland = 1;
                    break;
                }
            }

            if (nearland == 1)
                world.setmountainisland(c, 0);
        }
    }
}
