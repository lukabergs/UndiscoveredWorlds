// Template mountain chains and land shapes used by editing tools.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp.

#include "terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "classes.hpp"
#include "functions.hpp"
#include "planet.hpp"

using namespace std;

namespace
{
void drawshape(planet& world, int shapenumber, int centrex, int centrey, bool land, int baseheight, int conheight, boolshapetemplate shape[]);
void drawplateaushape(planet& world, int shapenumber, int centrex, int centrey, int baseheight, int platheight, std::vector<std::vector<int>>& plateaumap, boolshapetemplate chainland[]);
}

// This function creates chains of mountains, with associated land where appropriate.

void createchains(planet& world, int baseheight, int conheight, vector<vector<int>>& fractal, vector<vector<int>>& plateaumap, boolshapetemplate landshape[], boolshapetemplate chainland[], twointegers focuspoints[], int focustotal, int focaldistance, int mode)
{
    // mode=0: continental chains (with associated land).
    // mode=1: smaller mountains (no associated land, can form peninsulas).
    // mode=2: smaller mountains (no associated land, can't form peninsulas).
    // mode=3: island chains.
    // mode=5: like mode 2, but fewer.
    // mode=6: like 3, but fewer.
    // mode=7: hills.

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

    std::uint64_t chainseed = deterministiccontextseed(world.seed(), 0x4102);
    chainseed = deterministiccombine(chainseed, baseheight, conheight, focustotal, focaldistance, mode, lower ? 1 : 0, fewer ? 1 : 0);

    for (int n = 0; n < focustotal; n++)
        chainseed = deterministiccombine(chainseed, focuspoints[n].x, focuspoints[n].y);

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

    int width = world.width();
    int height = world.height();
    int maxelev = world.maxelevation();
    int sealevel = world.sealevel();

    vector<vector<unsigned char>> rangeheighttemplate(height + 1, vector<unsigned char>(height + 1, 0));
    vector<vector<unsigned char>> rangeridgetemplate(height + 1, vector<unsigned char>(height + 1, 0)); // These will hold the range we're currently putting on.

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
        int oldx = random(0, width);
        int oldy = random(0, height);

        if (mode == 0 || mode == 3) // Continental ranges should begin near the continental focus points, if possible. Islands should begin near their focus points too.
        {
            int thisfocus = random(0, focustotal);

            oldx = focuspoints[thisfocus].x;
            oldy = focuspoints[thisfocus].y; // Coordinates of the current blob.

            float xdiff = (float)random(1, 100);
            float ydiff = (float)random(1, 100);

            xdiff = xdiff / 100.0f;
            ydiff = ydiff / 100.0f;

            xdiff = xdiff * xdiff;
            ydiff = ydiff * ydiff;

            xdiff = xdiff * (float)focaldistance;
            ydiff = ydiff * (float)focaldistance;

            oldx = oldx + (int)randomsign(xdiff);
            oldy = oldy + (int)randomsign(ydiff);

            if (oldy < 0)
                oldy = 0;

            if (oldy > height)
                oldy = height;

            if (oldx<0 || oldx>width)
                oldx = wrap(oldx, width);
        }

        int chaindir = random(1, 8); // Chaindir is the direction of the chain, from 1-8 (sort of on the diagonals!).
        int maxranges = random(minchainlength, maxchainlength); // The chain won't have more ranges in it than this.
        int rangeno = 1;
        bool goahead = 1;

        int tempno = -1;

        if (!(mode == 2 && world.seawrap(oldx, oldy) == 1)) // If we're placing land-only mountains and we're in the sea, just skip the whole thing.
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
                            if (rangeridgetemplate[i][j] != 0)
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
                            if (rangeridgetemplate[i][j] != 0)
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
                            if (i >= 0 && i <= height && j >= 0 && j <= height && rangeridgetemplate[i][j] != 0)
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
                            if (rangeridgetemplate[i][j] != 0)
                            {
                                startx = i;
                                starty = j;

                                i = span;
                                j = 0;
                            }
                        }
                    }
                }

                if (chaindir == 1 || chaindir == 2)
                {
                    oldx = oldx + random(0, rangestartvar);
                    oldy = oldy + randomsign(random(0, rangestartvar));
                }

                if (chaindir == 3 || chaindir == 4)
                {
                    oldx = oldx + randomsign(random(0, rangestartvar));
                    oldy = oldy + random(0, rangestartvar);
                }

                if (chaindir == 5 || chaindir == 6)
                {
                    oldx = oldx - random(0, rangestartvar);
                    oldy = oldy + randomsign(random(0, rangestartvar));
                }

                if (chaindir == 7 || chaindir == 8)
                {
                    oldx = oldx + randomsign(random(0, rangestartvar));
                    oldy = oldy - random(0, rangestartvar);
                }

                // oldx and oldy are now the point where we want the new range to start.

                int x = oldx - startx;
                int y = oldy - starty;

                if (x <= 0 || x > width)
                    x = wrap(x, width);

                if (mode == 1 && oldx >= 0 && oldx <= width && oldy >= 0 && oldy <= height) // Peninsular mountains have to start within a certain distance of land.
                {
                    if (world.sea(oldx, oldy) == 1)
                    {
                        bool foundone = 0;

                        for (int i = oldx - mode1landmin; i <= oldx + mode1landmin; i++)
                        {
                            for (int j = oldy - mode1landmin; j <= oldy + mode1landmin; j++)
                            {
                                if (world.seawrap(i, j) == 0)
                                {
                                    foundone = 1;

                                    i = oldx + mode1landmin;
                                    j = oldy + mode1landmin;
                                }
                            }
                        }

                        if (foundone == 0)
                            goahead = 0;
                    }
                }

                if (mode == 2) // Non-continental, non-peninsular mountains can't start in the sea.
                {
                    if (world.seawrap(oldx, oldy) == 1)
                        goahead = 0;
                }

                if (mode == 3) // Island chains can't start on the land.
                {
                    if (oldx<0 || oldx>width || oldy<0 || oldy>height)
                        goahead = 0;
                    else
                    {
                        if (world.sea(oldx, oldy) == 0)
                            goahead = 0;
                    }
                }

                // x and y are now the coordinates of the top left of the new range template.
                // Now we need to see whether there is a clear area to paste it.

                if (mode != 2)
                {
                    for (int i = x; i <= x + span; i++)
                    {
                        for (int j = y; j <= y + span; j++)
                        {
                            if (world.mountainheightwrap(i, j) != 0)
                                goahead = 0;
                        }
                    }
                }

                int hislandadd = random(1000, 5000);

                if (goahead == 1) // If the area is clear
                {
                    for (int i = 0; i <= span; i++)
                    {
                        int ii = x + i;
                        if (ii<0 || ii>width)
                            ii = wrap(ii, width);

                        for (int j = 0; j <= span; j++)
                        {
                            int jj = y + j;

                            if (jj >= 0 && jj <= height)
                            {
                                if (world.mountainheight(ii, jj) == 0)
                                {
                                    if (i >= 0 && i <= height && j >= 0 && j <= height && rangeheighttemplate[i][j] != 0)
                                    {
                                        int h = rangeheighttemplate[i][j] * rangeinc;
                                        h = h + randomsign(rangelittleinc * random(0, 100));

                                        if (mode == 3)
                                            h = h + hislandadd;

                                        float frac = (float)fractal[ii][jj];
                                        frac = frac / (float)maxelev;
                                        h = (int)((float)h * frac);

                                        if (mode == 0 || mode == 1)
                                        {
                                            world.setmountainheight(ii, jj, h);
                                            world.setmountainridge(ii, jj, rangeridgetemplate[i][j]);
                                            world.setnom(ii, jj, conheight);
                                        }

                                        if (mode == 2) // Mountains of this kind cannot create their own land, so must stop at the sea.
                                        {
                                            if (world.map(ii, jj) >= conheight)
                                            {
                                                world.setmountainheight(ii, jj, h);
                                                world.setmountainridge(ii, jj, rangeridgetemplate[i][j]);

                                                if (fewer == 0)
                                                    world.setnom(ii, jj, conheight);
                                            }
                                        }

                                        if (mode == 3) // Islands are just strange and different.
                                        {
                                            if (h > sealevel)
                                            {
                                                world.setmountainheight(ii, jj, h - sealevel);
                                                world.setmountainridge(ii, jj, rangeridgetemplate[i][j]);
                                                world.setnom(ii, jj, conheight);
                                            }
                                        }

                                        if (mode == 0) // Only create substantial land if these are continental mountains.
                                        {
                                            if (random(1, landchance) == 1) // See if we'll paste a chunk of land nearby
                                            {
                                                int dist = random(minlanddist, maxlanddist);

                                                int movex = 0;
                                                int movey = 0;

                                                switch (chaindir)
                                                {
                                                case 1:
                                                    movex = 0 - dist / 2;
                                                    movey = 0 - dist;
                                                    break;

                                                case 2:
                                                    movex = dist / 2;
                                                    movey = 0 - dist;
                                                    break;

                                                case 3:
                                                    movex = dist;
                                                    movey = 0 - dist / 2;
                                                    break;

                                                case 4:
                                                    movex = dist;
                                                    movey = dist / 2;
                                                    break;

                                                case 5:
                                                    movex = dist / 2;
                                                    movey = dist;
                                                    break;

                                                case 6:
                                                    movex = 0 - dist / 2;
                                                    movey = dist;
                                                    break;

                                                case 7:
                                                    movex = 0 - dist;
                                                    movey = dist / 2;
                                                    break;

                                                case 8:
                                                    movex = 0 - dist;
                                                    movey = 0 - dist / 2;
                                                    break;
                                                }

                                                int landx = ii;
                                                int landy = jj;

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
                                                    landx = landx + movex;
                                                    landy = landy + movey;

                                                    landx = landx + randomsign(random(0, landvar));
                                                    landy = landy + randomsign(random(0, landvar));

                                                    if (landx<0 || landx>width)
                                                        landx = wrap(landx, width);

                                                    int shapenumber = random(0, 1);
                                                    drawshape(world, shapenumber, landx, landy, 1, baseheight, conheight, chainland);

                                                    if (plateau == 1) // Draw a bit of plateau, and reduce the height for the next pass.
                                                    {
                                                        drawplateaushape(world, shapenumber, landx, landy, baseheight, platheight, plateaumap, chainland);

                                                        float platreduce = (float)random(minplatreduce, maxplatreduce);
                                                        platreduce = platreduce / 100.0f;
                                                        platheight = (int)((float)platheight * platreduce);
                                                    }
                                                }
                                            }
                                        }

                                        if (mode == 1) // Mountains of this kind can create land just immediately around themselves.
                                        {
                                            if (world.nom(ii, jj) < conheight)
                                                world.setnom(ii, jj, conheight);
                                        }

                                        if (mode == 3) // Island mountains can create land just immediately around themselves.
                                        {
                                            if (h > sealevel)
                                            {
                                                if (world.nom(ii, jj) < conheight)
                                                    world.setnom(ii, jj, conheight);
                                            }
                                        }
                                    }
                                }
                                else // If this range goes off the map, and the chain.
                                    goahead = 0;

                            }
                        }
                    }

                    // Now we need to find the new oldx and oldy.

                    if (chaindir == 1 || chaindir == 2)
                    {
                        for (int i = span; i >= 0; i--)
                        {
                            for (int j = 0; j <= span; j++)
                            {
                                if (i >= 0 && i <= height && j >= 0 && j <= height && rangeheighttemplate[i][j] != 0)
                                {
                                    oldx = x + i;
                                    oldy = y + j;

                                    if (oldx<0 || oldx>width)
                                        oldx = wrap(oldx, width);

                                    i = 0;
                                    j = span;
                                }
                            }
                        }
                    }

                    if (chaindir == 3 || chaindir == 4)
                    {
                        for (int j = span; j >= 0; j--)
                        {
                            for (int i = 0; i <= span; i++)
                            {
                                if (rangeheighttemplate[i][j] != 0)
                                {
                                    oldx = x + i;
                                    oldy = y + j;

                                    if (oldx<0 || oldx>width)
                                        oldx = wrap(oldx, width);

                                    j = 0;
                                    i = span;
                                }
                            }
                        }
                    }

                    if (chaindir == 5 || chaindir == 6)
                    {
                        for (int i = 0; i <= span; i++)
                        {
                            for (int j = 0; j <= span; j++)
                            {
                                if (rangeheighttemplate[i][j] != 0)
                                {
                                    oldx = x + i;
                                    oldy = y + j;

                                    if (oldx<0 || oldx>width)
                                        oldx = wrap(oldx, width);

                                    i = span;
                                    j = span;
                                }
                            }
                        }
                    }

                    if (chaindir == 7 || chaindir == 8)
                    {
                        for (int j = 0; j <= span; j++)
                        {
                            for (int i = 0; i <= span; i++)
                            {
                                if (rangeheighttemplate[i][j] != 0)
                                {
                                    oldx = x + i;
                                    oldy = y + j;

                                    if (oldx<0 || oldx>width)
                                        oldx = wrap(oldx, width);

                                    j = span;
                                    i = span;
                                }
                            }
                        }
                    }

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
// Draws a shape from an image onto the map.

void drawshape(planet& world, int shapenumber, int centrex, int centrey, bool land, int baseheight, int conheight, boolshapetemplate shape[])
{
    int width = world.width();
    int height = world.height();

    int imheight = shape[shapenumber].ysize() - 1;
    int imwidth = shape[shapenumber].xsize() - 1;

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

    int leftr = random(0, 1); // If it's 1 then we reverse it left-right
    int downr = random(0, 1); // If it's 1 then we reverse it top-bottom

    int istart = 0, desti = imwidth + 1, istep = 1;
    int jstart = 0, destj = imheight + 1, jstep = 1;

    if (leftr == 1)
    {
        istart = imwidth;
        desti = -1;
        istep = -1;
    }

    if (downr == 1)
    {
        jstart = imwidth;
        destj = -1;
        jstep = -1;
    }

    int imap = -1;
    int jmap = -1;

    for (int i = istart; i != desti; i = i + istep)
    {
        imap++;
        jmap = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            jmap++;

            if (shape[shapenumber].point(i, j) == 1)
            {
                int xx = x + imap;
                int yy = y + jmap;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    if (land == 1)
                        world.setnom(xx, yy, conheight);
                    else
                    {
                        world.setnom(xx, yy, baseheight);
                        world.setmountainheight(xx, yy, 0);
                        world.setmountainridge(xx, yy, 0);
                    }
                }
            }
        }
    }
}
}

namespace
{
// Draws a shape from an image onto the plateau map.

void drawplateaushape(planet& world, int shapenumber, int centrex, int centrey, int baseheight, int platheight, vector<vector<int>>& plateaumap, boolshapetemplate chainland[])
{

    int width = world.width();
    int height = world.height();

    int imheight = chainland[shapenumber].ysize() - 1;
    int imwidth = chainland[shapenumber].xsize() - 1;

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

    int leftr = random(0, 1); // If it's 1 then we reverse it left-right
    int downr = random(0, 1); // If it's 1 then we reverse it top-bottom

    int istart = 0, desti = imwidth + 1, istep = 1;
    int jstart = 0, destj = imheight + 1, jstep = 1;

    if (leftr == 1)
    {
        istart = imwidth;
        desti = -1;
        istep = -1;
    }

    if (downr == 1)
    {
        jstart = imwidth;
        destj = -1;
        jstep = -1;
    }

    for (int i = istart; i != desti; i = i + istep)
    {
        for (int j = jstart; j != destj; j = j + jstep)
        {
            if (chainland[shapenumber].point(i, j) == 1)
            {
                int xx = x + 1;
                int yy = y + 1;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    if (world.sea(xx, yy) == 0)
                        plateaumap[xx][yy] = platheight;
                }
            }
        }
    }
}
}

// This function does something similar, but adds the fractal map onto the land map.

void fractaladdland(planet& world, vector<vector<int>>& fractal)
{
    int div = 300; //150; // The higher this number, the more highlands there will be.
    int maxamount = 25; // The higher this is, the flatter the flat areas will be.
    int minamount = 10; // The lower this is, the higher the high areas will be.

    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    // First we create a template map, which we will use to vary divamount.

    vector<vector<int>> templ(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
            templ[i][j] = fractal[i][j];
    }

    flip(templ, width, height, 1, 1);
    int offset = random(1, width);
    shift(templ, width, height, offset);

    // Now we use that template to apply the fractal to the map, to create highlands and lowlands.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) > sealevel)
            {
                int valtoadd = fractal[i][j] - sealevel;

                int divamount = templ[i][j];
                divamount = divamount / div;
                if (divamount > maxamount)
                    divamount = maxamount;
                if (divamount < minamount)
                    divamount = minamount;

                valtoadd = valtoadd / divamount;

                int newval = world.nom(i, j) + valtoadd;

                if (newval > maxelev)
                    newval = maxelev;

                if (newval > world.nom(i, j))
                    world.setnom(i, j, newval);
            }
        }
    }
}

// This scatters islands across areas where land has been removed.

void makearchipelagos(planet& world, vector<vector<bool>>& removedland, boolshapetemplate landshape[])
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();
    int conheight = sealevel + 50;

    int baseheight = sealevel - 4500;
    if (baseheight < 1)
        baseheight = 1;

    int grain = 8; // This is a fractal, for the islands
    float valuemod = 0.05f;
    float valuemod2 = valuemod;

    vector<vector<int>> islandfractal(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    createfractal(islandfractal, width, height, grain, valuemod, valuemod2, 1, maxelev, 0, 0);

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int thisrandom = random(1, maxelev * 100);

            int islandchance = 0;

            if (islandfractal[i][j] < maxelev / 5)
                islandchance = 70;

            if (islandfractal[i][j] < maxelev / 6)
                islandchance = 40;

            if (islandfractal[i][j] < maxelev / 7)
                islandchance = 30;

            int ii = i + width / 2;

            if (ii > width)
                ii = ii - width;

            if (islandchance != 0 && islandfractal[ii][j] > maxelev / 3 && removedland[i][j] == 1 && world.sea(i, j) == 1 && random(1, islandchance) == 1) // Put an island here.
            {
                int peakheight = random(500, 3000);
                terrain::detail::makemountainisland(world, i, j, peakheight);
            }
        }
    }

    int amount = 2; // No mountain islands this close to existing land (at least not marked as such).

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.mountainisland(i, j) == 1)
            {
                bool nearland = 0;

                for (int k = i - amount; k <= i + amount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (world.sea(kk, l) == 0 && world.mountainisland(kk, l) == 0)
                            {
                                nearland = 1;
                                k = i + amount;
                                l = j + amount;
                            }
                        }
                    }
                }

                if (nearland == 1)
                    world.setmountainisland(i, j, 0);
            }
        }
    }
}
