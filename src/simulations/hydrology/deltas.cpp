#include <algorithm>
#include <atomic>
#include <cmath>
#include <mutex>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "deltas.hpp"
#include "drainage.hpp"

using namespace std;

namespace
{
void placedelta(planet& world, int centrex, int centrey, int upriver, vector<vector<int>>& deltarivers);
void divertdeltarivers(planet& world, vector<vector<int>>& deltarivers);

// This creates a river delta.

void placedelta(planet& world, int centrex, int centrey, int upriver, vector<vector<int>>& deltarivers)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    int dir = world.riverdir(centrex, centrey);
    int maxbranches = 40;
    int branchestotal = 0;
    int branches = maxbranches; //20; //random(upriver*2,upriver*4); // Number of branches this delta will (hopefully) have.
    int size = upriver; // Possible distance from the centre that the branches can end.
    int addition = 2; // Amount to push branch ends out to sea.
    int x, y;

    if (branches > maxbranches)
        branches = maxbranches;

    twointegers destpoint;

    vector<vector<int>> branchends(maxbranches + 1, vector<int>(2));

    for (int i = 0; i < maxbranches + 1; i++)
    {
        for (int j = 0; j < 2; j++)
            branchends[i][j] = -1;
    }

    for (int n = 1; n <= branches; n++)
    {
        x = centrex;
        y = centrey;

        if (dir == 3 || dir == 7)
            y = y + randomsign(random(1, size));

        if (dir == 1 || dir == 5)
            x = x + randomsign(random(1, size));

        if (dir == 2 || dir == 6)
        {
            int amount = randomsign(random(1, size));
            x = x + amount;
            y = y + amount;
        }

        if (dir == 4 || dir == 8)
        {
            int amount = randomsign(random(1, size));
            x = x + amount;
            y = y - amount;
        }

        if (y < 0)
            y = 0;

        if (y > height)
            y = height;

        if (x<0 || x>width)
            x = wrap(x, width);

        if (world.sea(x, y) == 0)
        {
            destpoint = nearestsea(world, x, y, 0, 10, 1);

            x = destpoint.x;
            y = destpoint.y;
        }

        if (x != -1 && y != -1)
        {
            bool goahead = 1;

            if (world.deltadir(x, y) != 0 && random(1, 2) == 1)
                goahead = 0;

            if (goahead == 1)
            {
                branchestotal++;

                branchends[branchestotal][0] = x;
                branchends[branchestotal][1] = y;
            }
        }
    }

    x = centrex;
    y = centrey;

    // Now move upstream the required amount.

    for (int n = 1; n <= upriver; n++)
    {
        deltarivers[x][y] = 1; // Mark this bit of river to show that it's part of a river that's being turned into a delta.
        destpoint = getupstreamcell(world, x, y);

        x = destpoint.x;
        y = destpoint.y;

        if (x == -1 || y == -1)
            return;
    }

    deltarivers[x][y] = 1;

    int targetx = x;
    int targety = y; // These are the coordinates of the point where the delta branches are aiming for.

    // Now push those branch ends futher out to sea.

    for (int branch = 1; branch <= branchestotal; branch++)
    {
        x = branchends[branch][0];
        y = branchends[branch][1];

        if (x > targetx)
            x = x + addition;

        if (x < targetx)
            x = x - addition;

        if (y > targety)
            y = y + addition;

        if (y < targety)
            y = y + addition;

        if (x<0 || x>width)
            x = wrap(x, width);

        if (y < 0)
            y = 0;

        if (y > height)
            y = height;

        if (world.sea(x, y) == 0)
        {
            destpoint = nearestsea(world, x, y, 0, 10, 1);

            x = destpoint.x;
            y = destpoint.y;
        }

        branchends[branch][0] = x;
        branchends[branch][1] = y;
    }

    int leftx = targetx;
    int rightx = targetx;
    int lefty = targety;
    int righty = targety; // These will mark the corners of the whole area containing the delta.

    // Now work out the flow for each branch.

    float janflow = (float)world.riverjan(targetx, targety);
    float julflow = (float)world.riverjul(targetx, targety);

    float branchjanflow = janflow / (float)(branchestotal + 1);
    float branchjulflow = julflow / (float)(branchestotal + 1);

    // Now do each branch in turn.

    for (int branch = 1; branch <= branchestotal; branch++)
    {
        x = branchends[branch][0];
        y = branchends[branch][1];

        world.setdeltajan(x, y, (int)branchjanflow);
        world.setdeltajul(x, y, (int)branchjulflow);

        for (int n = 1; n <= 100; n++)
        {
            if (world.deltadir(x, y) == 0) // If we're drawing a new branch route
            {
                int rshift = targetx - x;
                int dshift = targety - y;

                if (abs(rshift) < 2 && abs(dshift) < 2) // If we're next to the target go straight towards it
                {
                    if (targetx == x && targety == y - 1)
                        dir = 1;

                    if (targetx == x + 1 && targety == y - 1)
                        dir = 2;

                    if (targetx == x + 1 && targety == y)
                        dir = 3;

                    if (targetx == x + 1 && targety == y + 1)
                        dir = 5;

                    if (targetx == x - 1 && targety == y + 1)
                        dir = 6;

                    if (targetx == x - 1 && targety == y)
                        dir = 7;

                    if (targetx == x - 1 && targety == y - 1)
                        dir = 8;
                }
                else // Find a direction going roughly in the right direction
                {
                    if (rshift >= 0 && dshift >= 0) // Going roughly southeast
                    {
                        if (rshift > dshift)
                        {
                            if (random(1, 2) == 1)
                                dir = 3;
                            else
                                dir = 4;
                        }
                        else
                        {
                            if (random(1, 2) == 1)
                                dir = 4;
                            else
                                dir = 5;
                        }
                    }

                    if (rshift >= 0 && dshift < 0) // Going roughly northeast
                    {
                        if (rshift > 0 - dshift)
                        {
                            if (random(1, 2) == 1)
                                dir = 3;
                            else
                                dir = 2;
                        }
                        else
                        {
                            if (random(1, 2) == 1)
                                dir = 2;
                            else
                                dir = 1;
                        }
                    }

                    if (rshift < 0 && dshift >= 0) // Going roughly southwest
                    {
                        if (0 - rshift > dshift)
                        {
                            if (random(1, 2) == 1)
                                dir = 7;
                            else
                                dir = 6;
                        }
                        else
                        {
                            if (random(1, 2) == 1)
                                dir = 6;
                            else
                                dir = 5;
                        }
                    }

                    if (rshift < 0 && dshift < 0) // Going roughly northwest
                    {
                        if (rshift < dshift)
                        {
                            if (random(1, 2) == 1)
                                dir = 7;
                            else
                                dir = 8;
                        }
                        else
                        {
                            if (random(1, 2) == 1)
                                dir = 8;
                            else
                                dir = 1;
                        }
                    }
                    world.setdeltadir(x, y, dir);
                }
            }
            else
                dir = world.deltadir(x, y);

            if (dir == 8 || dir == 1 || dir == 2)
                y--;

            if (dir == 4 || dir == 5 || dir == 6)
                y++;

            if (dir == 2 || dir == 3 || dir == 4)
                x++;

            if (dir == 6 || dir == 7 || dir == 8)
                x--;

            if (x < 0)
                x = width;

            if (x > width)
                x = 0;

            if (y >= 0 && y <= height)
            {
                world.setdeltajan(x, y, world.deltajan(x, y) + (int)branchjanflow);
                world.setdeltajul(x, y, world.deltajul(x, y) + (int)branchjulflow);

                if (world.nom(x, y) > sealevel)
                    world.setnom(x, y, sealevel + 1);

                if (x < leftx)
                    leftx = x;

                if (x > rightx)
                    rightx = x;

                if (y < lefty)
                    lefty = y;

                if (y > righty)
                    righty = y;

            }
            else
            {
                x = targetx;
                y = targety;
            }

            if (x == targetx && y == targety)
            {
                n = 100;
            }
        }
    }

    leftx--;
    rightx++;
    lefty--;
    righty++;

    if (leftx > rightx)
    {
        leftx = 0;
        rightx = width;
    }

    // Now fill in missing tiles.

    for (int i = targetx - 1; i <= targetx + 1; i++)
    {
        int ii = i;

        if (ii<0 || ii>width)
            ii = wrap(ii, width);

        for (int j = targety - 1; j <= targety + 1; j++)
        {
            if (j >= 0 && j <= height)
            {
                // We're looking at each tile around the destination one.
                // For each of these, see whether any delta branch is pointing into it and add those amounts to its flow

                int janflow = 0;
                int julflow = 0;

                for (int k = ii - 1; k <= ii + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (kk != targetx || l != targety) // Don't do the actual destination tile
                            {
                                int dir = world.deltadir(kk, l);

                                if (dir > 0) // If there's a delta branch here
                                {
                                    twointegers destpoint;

                                    destpoint = getflowdestination(world, kk, l, dir);

                                    if (destpoint.x == i && destpoint.y == j)
                                    {
                                        janflow = janflow + world.deltajan(kk, l);
                                        julflow = julflow + world.deltajul(kk, l);
                                    }
                                }
                            }
                        }
                    }
                }

                if (janflow > 0 || julflow > 0) // If there are delta branches flowing into this tile
                {
                    // First find the direction to the destination tile

                    int dir = getdir(ii, j, x, y);

                    // Now put a branch on our tile.

                    world.setdeltadir(ii, j, dir);
                    world.setdeltajan(ii, j, janflow);
                    world.setdeltajul(ii, j, julflow);
                }
            }
        }
    }

    // Now reduce the river below that point.

    vector<vector<int>> removedrivers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int janreduce = world.riverjan(targetx, targety) - (int)(branchjanflow * 4.0f); // Shouldn't really be multiplied, but this is to make it more visible on the map.
    int julreduce = world.riverjul(targetx, targety) - (int)(branchjulflow * 4.0f);

    reduceriver(world, janreduce, julreduce, removedrivers, 1, targetx, targety);

    //Now we simply make one final delta cell in the target cell, bringing the branches back together to meet the main river.

    int ux = -1;
    int uy = -1;

    x = targetx;
    y = targety;

    destpoint = getupstreamcell(world, x, y);

    ux = destpoint.x;
    uy = destpoint.y;

    if (ux == x && uy == y - 1)
        dir = 1;

    if (ux == x + 1 && uy == y - 1)
        dir = 2;

    if (ux == x + 1 && uy == y)
        dir = 3;

    if (ux == x + 1 && uy == y + 1)
        dir = 4;

    if (ux == x && uy == y + 1)
        dir = 5;

    if (ux == x - 1 && uy == y + 1)
        dir = 6;

    if (ux == x - 1 && uy == y)
        dir = 7;

    if (ux == x - 1 && uy == y - 1)
        dir = 8;

    world.setdeltadir(x, y, dir);

    world.setdeltajan(x, y, 0 - (int)(branchjanflow * (float)branchestotal));
    world.setdeltajul(x, y, 0 - (int)(branchjulflow * (float)branchestotal));

    deltarivers[x][y] = 1;

}

// This function diverts rivers that run through deltas, so that they follow the pattern of the delta branches more closely.

void divertdeltarivers(planet& world, vector<vector<int>>& deltarivers)
{
    int width = world.width();
    int height = world.height();

    vector<vector<int>> removedrivers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));
    vector<vector<int>> checked(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    int riverno = 0;

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (checked[i][j] == 0 && world.sea(i, j) == 0 && world.riverdir(i, j) != 0)
            {
                int x = i;
                int y = j;

                bool keepgoing = 1;
                bool found = 0;
                bool alreadydeleted = 0;

                do
                {                    
                    checked[x][y] = 1; // Mark this as checked so we won't try to follow any rivers from this point again.

                    int dir = world.riverdir(x, y);

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

                    if (world.sea(x, y) == 1) // If it reaches the sea
                        keepgoing = 0;

                    if (checked[x][y] == 1) // If it reaches a river that's already been done
                    {
                        keepgoing = 0;
                    }

                    if (deltarivers[x][y] == 1) // If it reaches a river that's part of the delta
                    {
                        keepgoing = 0;
                    }

                    if (keepgoing == 1 && (world.deltajan(x, y) != 0 || world.deltajul(x, y) != 0)) // If it reaches an actual delta
                    {
                        keepgoing = 0;
                        found = 1;
                    }

                    // Now check whether it's about to cross over a delta diagonally. If it is, we will need to divert it into one of the neighbouring delta cells.

                    dir = world.riverdir(x, y);

                    if (dir == 2 || dir == 4 || dir == 6 || dir == 8)
                    {
                        int dx = -1;
                        int dy = -1; // Coordinates of the cell that this river is about to head into

                        int mx = -1;
                        int my = -1;
                        int nx = -1;
                        int ny = -1; // These are the coordinates of neighbouring cells that this river isn't actually heading into, to check for delta branches.

                        bool mfound = 0;
                        bool nfound = 0; // To tell whether we have found a delta branch in either of these cells.

                        if (dir == 2)
                        {
                            dx = x + 1;
                            dy = y - 1;

                            if (dx > width)
                                dx = 0;

                            if (dy >= 0 && (world.deltajan(dx, dy) == 0 || world.deltajul(dx, dy) == 0)) // If this river isn't about to run into a delta branch
                            {
                                // Check the ones to the sides to see if they contain a delta branch.

                                mx = x;
                                my = y - 1;

                                nx = x + 1;
                                ny = y;

                                if (nx > width)
                                    nx = 0;

                                if (my >= 0 && (world.deltajan(mx, my) != 0 || world.deltajul(mx, my) != 0))
                                    mfound = 1;

                                if (world.deltajan(nx, ny) != 0 || world.deltajul(nx, ny) != 0)
                                    nfound = 1;
                            }
                        }

                        if (dir == 4)
                        {
                            dx = x + 1;
                            dy = y + 1;

                            if (dx > width)
                                dx = 0;

                            if (dy <= height && (world.deltajan(dx, dy) == 0 || world.deltajul(dx, dy) == 0)) // If this river isn't about to run into a delta branch
                            {
                                // Check the ones to the sides to see if they contain a delta branch.

                                mx = x;
                                my = y + 1;

                                nx = x + 1;
                                ny = y;

                                if (nx > width)
                                    nx = 0;

                                if (my <= height && (world.deltajan(mx, my) != 0 || world.deltajul(mx, my) != 0))
                                    mfound = 1;

                                if (world.deltajan(nx, ny) != 0 || world.deltajul(nx, ny) != 0)
                                    nfound = 1;
                            }
                        }

                        if (dir == 6)
                        {
                            dx = x - 1;
                            dy = y + 1;

                            if (dx < 0)
                                dx = width;

                            if (dy <= height && (world.deltajan(dx, dy) == 0 || world.deltajul(dx, dy) == 0)) // If this river isn't about to run into a delta branch
                            {
                                // Check the ones to the sides to see if they contain a delta branch.

                                mx = x - 1;
                                my = y;

                                nx = x;
                                ny = y + 1;

                                if (mx < 0)
                                    mx = width;

                                if (world.deltajan(mx, my) != 0 || world.deltajul(mx, my) != 0)
                                    mfound = 1;

                                if (ny <= height && (world.deltajan(nx, ny) != 0 || world.deltajul(nx, ny) != 0))
                                    nfound = 1;
                            }
                        }

                        if (dir == 8)
                        {
                            dx = x - 1;
                            dy = y - 1;

                            if (dx < 0)
                                dx = width;

                            if (dy >= 0 && (world.deltajan(dx, dy) == 0 || world.deltajul(dx, dy) == 0)) // If this river isn't about to run into a delta branch
                            {
                                // Check the ones to the sides to see if they contain a delta branch.

                                mx = x;
                                my = y - 1;

                                nx = x - 1;
                                ny = y;

                                if (nx < 0)
                                    nx = width;

                                if (my >= 0 && (world.deltajan(mx, my) != 0 || world.deltajul(mx, my) != 0))
                                    mfound = 1;

                                if (world.deltajan(nx, ny) != 0 || world.deltajul(nx, ny) != 0)
                                    nfound = 1;
                            }
                        }

                        if (mfound == 1 || nfound == 1) // We found one!
                        {
                            keepgoing = 1;

                            int tx = nx;
                            int ty = ny; // The new target to divert our river into. By default it's n.

                            if (mfound == 1 && nfound == 1) // If both possibilities contain a delta branch, we need the downstream one. (Remember it will show as the upstream one because delta branches are calculated in reverse.)
                            {
                                twointegers dest = getflowdestination(world, mx, my, world.deltadir(mx, my));

                                if (dest.x == nx && dest.y == ny) // m is flowing into n. That means we want to divert into m, as it's really downstream.
                                {
                                    tx = mx;
                                    ty = my;
                                }
                            }

                            if (nfound == 0)
                            {
                                tx = mx;
                                ty = my;
                            }
                            // tx and ty now contains the coordinates of the tile we want to divert our river into.

                            // Delete the river from this point onwards.

                            int janload = world.riverjan(x, y);
                            int julload = world.riverjul(x, y);

                            riverno++;

                            removeriver(world, removedrivers, riverno, x, y);

                            // Now recreate it in this cell, pointing to the new destination cell.

                            world.setriverjan(x, y, janload);
                            world.setriverjul(x, y, julload);
                            world.setriverdir(x, y, getdir(x, y, tx, ty));

                            alreadydeleted = 1;

                            // Now on the next pass of this loop, the program will move x and y in this new direction and discover that it's on a delta branch!
                        }
                    }
                } while (keepgoing == 1);

                if (found == 1) // If we've actually hit a delta, time to sort out this river!
                {
                    int janload = world.riverjan(x, y);
                    int julload = world.riverjul(x, y);

                    // First, delete the river from this point onwards.

                    if (alreadydeleted == 0)
                    {
                        riverno++;

                        removeriver(world, removedrivers, riverno, x, y);
                    }

                    // Now, recreate the river, following the line of the delta branch.

                    keepgoing = 1;
                    twointegers dest, check;

                    int tally = 0;

                    do
                    {
                        tally++;
                        
                        // We need to find a delta branch that's flowing into this tile.

                        int dir;
                        dest.x = -1;
                        dest.y = -1;

                        for (int k = x - 1; k <= x + 1; k++)
                        {
                            int kk = k;

                            if (kk<0 || kk>width)
                                kk = wrap(kk, width);

                            for (int l = y - 1; l <= y + 1; l++)
                            {
                                if (l >= 0 && l <= height)
                                {
                                    if ((k != x || l != y) && world.deltadir(kk, l) != 0)
                                    {
                                        check = getflowdestination(world, kk, l, world.deltadir(kk, l));

                                        if (check.x == x && check.y == y)
                                        {
                                            dest.x = kk;
                                            dest.y = l;

                                        }
                                    }
                                }
                            }
                        }

                        if (dest.x == -1) // We didn't find one
                            keepgoing = 0;
                        else
                        {
                            dir = getdir(x, y, dest.x, dest.y);

                            world.setriverdir(x, y, dir);
                            world.setriverjan(x, y, world.riverjan(x, y) + janload);
                            world.setriverjul(x, y, world.riverjul(x, y) + julload);

                            if (dir == 8 || dir == 1 || dir == 2)
                                y--;

                            if (dir == 4 || dir == 5 || dir == 6)
                                y++;

                            if (dir == 2 || dir == 3 || dir == 4)
                                x++;

                            if (dir == 6 || dir == 7 || dir == 8)
                                x--;

                            if (x < 0)
                                x = width;

                            if (x > width)
                                x = 0;
                        }

                        if (tally > 10000) // In case of infinite loops...
                            keepgoing = 0;

                    } while (keepgoing == 1);
                }
            }
        }
    }
}
}

// This function puts river deltas on the global map.

void createriverdeltas(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int glaciertemp = world.glaciertemp();
    float gravity = world.gravity();

    vector<vector<int>> deltarivers(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0)); // Marks all the rivers that have been turned into deltas.

    twointegers frompoint, destpoint;

    int margin = 10; // Don't do any closer to the east/west edges of the map than this.
    int deltachance = 200000; //65000; // The lower this is, the more deltas there will be.

    if (gravity > 1.0) // Higher gravity means fewer deltas.
    {
        float fchance = (float)deltachance;
        fchance = fchance * gravity * gravity;

        deltachance = (int)fchance;
    }

    int tidefactor = 8; // The higher this is, the more deltas there will be.
    int range = 10; // Minimum gap between deltas.

    for (int i = margin; i <= width - margin; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            int riverjan = world.riverjan(i, j);
            int riverjul = world.riverjul(i, j);

            if ((riverjan != 0 || riverjul != 0) && world.avetemp(i, j) > glaciertemp) // There's a river here and it's not too cold
            {
                int flowtotal = riverjan + riverjul;

                if (random(1, deltachance) < flowtotal)
                {
                    if (world.sea(i, j) == 0)
                    {
                        destpoint = getflowdestination(world, i, j, 0);

                        if (destpoint.x != -1 && destpoint.y != -1)
                        {
                            if (world.sea(destpoint.x, destpoint.y) == 1)
                            {
                                int tide = world.tide(i, j);

                                bool found = 0;

                                if (random(1, tidefactor) > tide)
                                {
                                    for (int k = i - range; k <= i + range; k++)
                                    {
                                        int kk = k;

                                        if (kk<0 || kk>width)
                                            kk = wrap(kk, width);

                                        for (int l = j - range; l <= j + range; l++)
                                        {
                                            if (l >= 0 && l <= height)
                                            {
                                                if (world.deltadir(kk, l) != 0)
                                                {
                                                    found = 1;
                                                    k = i + range;
                                                    l = j + range;
                                                }
                                            }
                                        }
                                    }
                                    if (found == 0)
                                    {
                                        int x = destpoint.x;
                                        int y = destpoint.y;

                                        int upriver = (world.riveraveflow(x, y) / 10000) + randomsign(random(0, 3));

                                        if (upriver < 2)
                                            upriver = 2;

                                        if (upriver > 10)
                                            upriver = 10;

                                        //upriver=upriver*2; // Just to make them bigger to check more easily!

                                        // Upriver is how far upstream from the coast we'll start the delta.

                                        placedelta(world, x, y, upriver, deltarivers);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Now we flatten all land around the delta branches.

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.deltadir(i, j) > 0)
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
                            if (world.sea(kk, l) == 0)
                            {
                                world.setnom(kk, l, sealevel + 1);

                                if (world.deltadir(kk, l) == 0)
                                    world.setdeltadir(kk, l, -1);
                            }
                        }
                    }
                }
            }
        }
    }

    // Now we lower any land that borders flat land.

    vector<vector<int>> donethese(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) == sealevel + 1)
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
                            if (donethese[kk][l] == 0 && world.nom(kk, l) > sealevel + 1)
                            {
                                int elev = world.nom(kk, l) - sealevel;
                                elev = elev / 2;

                                world.setnom(kk, l, sealevel + elev);

                                donethese[kk][l] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Now we divert other rivers that flow through the delta areas, so they flow more naturally in the delta shape.

    divertdeltarivers(world, deltarivers);
}

// This goes over the map and ensures that rivers don't inexplicably grow too much within delta regions.

void checkrivers(planet& world)
{
    //highres_timer_t timer("Check Rivers"); // 666: 4156, 999: 6442 => 666: 2926, 999: 2907
    int width = world.width();
    int height = world.height();

    int maxsource = 50; // Largest size a river source can be (jan and jul each)

    // First, go through and make sure there aren't any absurd river sources (these may be accidentally created in the delta regions)

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.deltadir(i, j) != 0 && checkwaterinflow(world, i, j) == 0)
                {
                    if (world.riverjan(i, j) > maxsource)
                        world.setriverjan(i, j, maxsource);

                    if (world.riverjul(i, j) > maxsource)
                        world.setriverjul(i, j, maxsource);
                }
            }
        }
    });

    // Now go through and ensure that no river grows weirdly for no apparent reason.

    int amount = 4;
    vector<vector<unsigned char>> neardelta(ARRAYWIDTH, vector<unsigned char>(ARRAYHEIGHT, 0));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                bool founddelta = 0;

                for (int k = i - amount; k <= i + amount && founddelta == 0; k++)
                {
                    int kk = k;

                    if (kk < 0 || kk > width)
                        kk = wrap(kk, width);

                    for (int l = j - amount; l <= j + amount; l++)
                    {
                        if (l >= 0 && l <= height && world.deltadir_no_bounds_check(kk, l) != 0)
                        {
                            founddelta = 1;
                            break;
                        }
                    }
                }

                neardelta[i][j] = founddelta ? 1 : 0;
            }
        }
    });

    bool found = 0;

    do
    {
        vector<vector<int>> nextjan(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
        vector<vector<int>> nextjul(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, -1));
        atomic<bool> anychanged(false);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int j = startrow; j <= endrow; j++)
            {
                for (int i = 0; i <= width; i++)
                {
                    if (neardelta[i][j] != 0)
                    {
                        twointegers inflow = gettotalinflow(world, i, j);

                        if (world.riverjan(i, j) > inflow.x)
                        {
                            nextjan[i][j] = inflow.x;
                            anychanged.store(true, memory_order_relaxed);
                        }

                        if (world.riverjul(i, j) > inflow.y)
                        {
                            nextjul[i][j] = inflow.y;
                            anychanged.store(true, memory_order_relaxed);
                        }
                    }
                }
            }
        });

        found = anychanged.load(memory_order_relaxed);

        if (found)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int j = startrow; j <= endrow; j++)
                {
                    for (int i = 0; i <= width; i++)
                    {
                        if (nextjan[i][j] != -1)
                            world.setriverjan(i, j, nextjan[i][j]);

                        if (nextjul[i][j] != -1)
                            world.setriverjul(i, j, nextjul[i][j]);
                    }
                }
            });
        }
    } while (found == 1);
}
