#include <algorithm>
#include <cmath>
#include <mutex>
#include <queue>
#include <vector>

#include "classes.hpp"
#include "planet.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "climate_landforms.hpp"
#include "lakes.hpp"
#include "land_distance.hpp"

using namespace std;

namespace
{
void checkergelevation(planet& world);
int lowestergelevation(planet& world, int startx, int starty, vector<vector<bool>>& thiserg);
void drawspeciallake(planet& world, int shapenumber, int centrex, int centrey, int lakeno, vector<vector<int>>& thislake, boolshapetemplate laketemplate[], bool canoverlap, int special);

// This ensures that the elevation of each erg remains constant.

void checkergelevation(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();

    vector<vector<bool>> checked(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.special(i, j) == 120 && checked[i][j] == 0)
            {
                vector<vector<bool>> thiserg(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, 0));
                
                // First, find the lowest elevation within this erg.

                int lowestelev = lowestergelevation(world, i, j, thiserg);

                if (lowestelev < 2)
                    lowestelev = 2;

                if (lowestelev <= sealevel)
                    lowestelev = sealevel + 1;

                // Now apply that to the whole erg.

                for (int k = 0; k <= width; k++)
                {
                    for (int l = 0; l <= height; l++)
                    {
                        if (thiserg[k][l] == 1)
                        {
                            checked[k][l] = 1;
                            world.setlakesurface(k, l, lowestelev);
                        }
                    }
                }
            }
        }
    }
}

// This finds the lowest elevation in an erg.

int lowestergelevation(planet& world, int startx, int starty, vector<vector<bool>>& thiserg)
{
    int width = world.width();
    int height = world.height();

    int lowest = 1000000;

    int row[] = { -1,0,0,1 };
    int col[] = { 0,-1,1,0 };

    twointegers node;

    node.x = startx;
    node.y = starty;

    queue<twointegers> q; // Create a queue
    q.push(node); // Put our starting node into the queue

    while (q.empty() != 1) // Keep going while there's anything left in the queue
    {
        node = q.front(); // Take the node at the front of the queue
        q.pop(); // And then pop it out of the queue

        if (node.x >= 0 && node.x <= width && node.y >= 0 && node.y <= height && world.special(node.x, node.y) == 120)
        {
            int thiselev = world.nom(node.x, node.y);

            if (thiselev < lowest)
                lowest = thiselev;

            thiserg[node.x][node.y] = 1;

            for (int k = 0; k < 4; k++) // Look at the four neighbouring nodes in turn
            {
                twointegers nextnode;

                nextnode.x = node.x + row[k];
                nextnode.y = node.y + col[k];

                if (nextnode.x > width)
                    nextnode.x = 0;

                if (nextnode.x < 0)
                    nextnode.x = width;

                if (nextnode.y >= 0 && nextnode.y < height)
                {
                    if (world.special(nextnode.x, nextnode.y) == 120 && thiserg[nextnode.x][nextnode.y] == 0) // If this node is erg
                    {
                        thiserg[nextnode.x][nextnode.y] = 1;
                        q.push(nextnode); // Put that node onto the queue
                    }
                }
            }
        }
    }

    return lowest;
}

// Puts a lake template onto the lakemap, but as a "special", e.g. salt pans etc.

void drawspeciallake(planet& world, int shapenumber, int centrex, int centrey, int lakeno, vector<vector<int>>& thislake, boolshapetemplate laketemplate[], bool canoverlap, int special)
{
    int width = world.width();
    int height = world.height();

    if (centrey == 0 || centrey == height)
    {
        world.setlakesurface(centrex, centrey, 0);
        return;
    }

    twointegers nearestsea, flow;

    int minlakedistance = 4; // Minimum distance between lakes.
    int maxmountainlakeheight = 1000; // Ignore mountains smaller than this.
    int mindepth = 5; // Minimum depth.

    int surfaceheight;

    int imheight = laketemplate[shapenumber].ysize() - 1;
    int imwidth = laketemplate[shapenumber].xsize() - 1;

    int x = centrex - imwidth / 2; // Coordinates of the top left corner.
    int y = centrey - imheight / 2;

    if (x<0 || x>width)
        x = wrap(x, width);

    int leftx = centrex;
    int lefty = centrey;
    int rightx = centrex;
    int righty = centrey; // These are the coordinates of the furthermost pixels of the lake.
    bool wrapped = 0; // If this is 1, the lake wraps over the edge of the map.

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

    // First, check that the location is clear. Also, get the height of the lowest point in it, and make that the surface height.

    int mapi = -1;
    int mapj = -1;

    surfaceheight = 1000000000;

    bool tooclose = 0;

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            if (laketemplate[shapenumber].point(i, j) == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (xx != 0 && xx <= width && world.nom(xx, yy) < surfaceheight)
                    surfaceheight = world.nom(xx, yy);

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                        xx = wrap(xx, width);

                    if (world.sea(xx, yy) != 0 || world.mountainheight(xx, yy) >= maxmountainlakeheight) // Don't try to put them on top of sea or mountains...
                        tooclose = 1;
                    else
                    {
                        int clim = 1; // This marks whether the climate is appropriate.

                        if (special == 110 || special == 120)
                        {
                            for (int k = xx - 1; k <= xx + 1; k++)
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = yy - 1; l <= yy + 1; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.climate(kk, l) != 5 || world.riverjan(kk, l) > 0 || world.riverjul(kk, l) > 0)
                                        {
                                            clim = 0;

                                            k = xx + 1;
                                            l = yy + 1;
                                        }
                                    }
                                }
                            }
                        }

                        if (clim != 1) // Don't put lakes onto inappropriate climates.
                            tooclose = 1;
                        else
                        {
                            if (world.outline(xx, yy) == 1)
                                tooclose = 1;
                            else
                            {
                                for (int xxx = xx - minlakedistance; xxx <= xx + minlakedistance; xxx++)
                                {
                                    int xxxx = xxx;

                                    if (xxxx<0 || xxxx>width)
                                        xxxx = wrap(xxxx, width);

                                    for (int yyy = yy - minlakedistance; yyy <= yy + minlakedistance; yyy++)
                                    {
                                        if (yyy >= 0 && yyy <= height)
                                        {
                                            if (thislake[xxxx][yyy] != 0 && thislake[xxxx][yyy] != lakeno && canoverlap == 0) // There's another lake here already!
                                                tooclose = 1;

                                            if (world.lakesurface(xxxx, yyy) != 0 && thislake[xxxx][yyy] != lakeno && (canoverlap == 0 || (special != 120 || world.special(xxxx, yyy) != 120))) // There's another kind of lake here already! (Doesn't apply when both are ergs.)
                                                tooclose = 1;

                                            if (world.sea(xxxx, yyy) == 1) // && special != 120) // Sea here (doesn't apply to ergs).
                                                tooclose = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (tooclose == 1)
        return;

    // Now, do the actual thing.

    mapi = -1;
    mapj = -1;

    for (int i = istart; i != desti; i = i + istep)
    {
        mapi++;
        mapj = -1;

        for (int j = jstart; j != destj; j = j + jstep)
        {
            mapj++;

            if (laketemplate[shapenumber].point(i, j) == 1)
            {
                int xx = x + mapi;
                int yy = y + mapj;

                if (yy >= 0 && yy <= height)
                {
                    if (xx<0 || xx>width)
                    {
                        xx = wrap(xx, width);
                        wrapped = 1;
                    }

                    world.setlakesurface(xx, yy, surfaceheight); // Put the lake on this bit of the lake map.
                    world.setnom(xx, yy, surfaceheight);
                    world.setspecial(xx, yy, special);

                    thislake[xx][yy] = lakeno; // Mark it on the keeping-track array too.

                    if (world.mountainheight(xx, yy) != 0) // Remove any mountains that might be here.
                    {
                        for (int dir = 1; dir <= 8; dir++)
                            deleteridge(world, xx, yy, dir);
                    }

                    if (xx < leftx)
                        leftx = xx;

                    if (xx > rightx)
                        rightx = xx;

                    if (yy < lefty)
                        lefty = yy;
                    if (yy > righty)
                        righty = yy;
                }
            }
        }
    }

    if (wrapped == 1)
    {
        leftx = 0;
        lefty = 0;
        rightx = width;
        righty = height;
    }

    // Now mark a "start point" somewhere on the edge of the "lake", which we will use when it comes to the regional map.

    makelakestartpoint(world, thislake, lakeno, leftx, lefty, rightx, righty);
}
}

void broadenfastlemterrainfromrivers(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int maxelev = world.maxelevation();
    const int maxriverdistance = tuning::terrain::fastlem::postRiverMaximumDistance;
    const int maxmountaindistance = tuning::terrain::fastlem::postRiverMountainDistance;

    vector<vector<short>> riverdistance(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, static_cast<short>(maxriverdistance + 1)));
    vector<vector<short>> mountaindistance(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, static_cast<short>(maxmountaindistance + 1)));

    const bool hasrivers = terraindetail::buildlanddistancefield(world, riverdistance, maxriverdistance, [&](int x, int y)
    {
        return world.riverdir(x, y) != 0 || world.riverjan(x, y) != 0 || world.riverjul(x, y) != 0 || world.lakesurface(x, y) != 0 || world.riftlakesurface(x, y) != 0;
    });

    const bool hasmountains = terraindetail::buildlanddistancefield(world, mountaindistance, maxmountaindistance, [&](int x, int y)
    {
        return world.mountainheight(x, y) >= tuning::terrain::fastlem::postRiverMinimumMountainHeight;
    });

    if (hasrivers == false || hasmountains == false)
        return;

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int j = startrow; j <= endrow; j++)
        {
            for (int i = 0; i <= width; i++)
            {
                if (world.sea(i, j) == 1 || world.outline(i, j) == 1)
                    continue;

                if (world.lakesurface(i, j) != 0 || world.riftlakesurface(i, j) != 0)
                    continue;

                const int thisriverdistance = riverdistance[i][j];
                const int thismountaindistance = mountaindistance[i][j];

                if (thisriverdistance <= 0 || thisriverdistance > maxriverdistance || thismountaindistance > maxmountaindistance)
                    continue;

                float riverfactor = static_cast<float>(thisriverdistance) / static_cast<float>(maxriverdistance);
                riverfactor = powf(clamp(riverfactor, 0.0f, 1.0f), tuning::terrain::fastlem::postRiverDistanceExponent);

                float mountainfactor = 1.0f - static_cast<float>(thismountaindistance) / static_cast<float>(maxmountaindistance);
                mountainfactor = powf(clamp(mountainfactor, 0.0f, 1.0f), tuning::terrain::fastlem::postRiverMountainExponent);

                float elevationfactor = static_cast<float>(world.nom(i, j) - sealevel) / static_cast<float>(max(1, tuning::terrain::fastlem::postRiverElevationRange));
                elevationfactor = clamp(elevationfactor, 0.0f, 1.0f);

                const float upliftstrength = riverfactor * mountainfactor * (0.35f + elevationfactor * 0.65f);
                const int uplift = static_cast<int>(roundf(static_cast<float>(tuning::terrain::fastlem::postRiverMaximumUplift) * upliftstrength));

                if (uplift < tuning::terrain::fastlem::postRiverMinimumUplift)
                    continue;

                int newheight = world.nom(i, j) + uplift;

                if (newheight > maxelev - 1)
                    newheight = maxelev - 1;

                world.setnom(i, j, newheight);
            }
        }
    });
}

// This function puts ergs in deserts.

void createergs(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate shape[])
{
    int width = world.width();
    int height = world.height();
    int worldsize = world.size();

    vector<vector<short>> ergprobability(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, 0)); // This is the probability of making an erg here.
    vector<vector<short>> overlapprobability(ARRAYWIDTH, vector<short>(ARRAYHEIGHT, 0)); // This is the probability of this erg being able to overlap its neighbours.

    int clusternumber = random(2, 50); // 50; Number of possible clusters of erg on the map.

    if (worldsize == 1)
        clusternumber = clusternumber * 4;

    if (worldsize == 2)
        clusternumber = clusternumber * 16;


    int minrad = 10;
    int maxrad = 40; // Possible sizes of the clusters

    for (int n = 0; n < clusternumber; n++) // Draw some shapes on this array, to allow for erg clusters.
    {
        int shapenumber = random(1, 11);

        int probability = random(20, 200);

        int overlap = random(1, 100);
        
        int centrex = random(0, width);
        int centrey = random(0, height);

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

                        ergprobability[xx][yy] = probability;
                        overlapprobability[xx][yy] = overlap;

                    }
                }
            }
        }

        /*
        int radius = random(minrad, maxrad);

        for (int i = -radius; i <= radius; i++)
        {
            int ii = centrex + i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = -radius; j <= radius; j++)
            {
                int jj = centrey + j;

                if (jj >= 0 && jj <= height)
                {
                    if (i * i + j * j < radius * radius + radius)
                        ergprobability[ii][jj] = 2;
                }
            }
        }

        radius = (int)((float)radius * 1.5f);

        for (int i = -radius; i <= radius; i++)
        {
            int ii = centrex + i;

            if (ii<0 || ii>width)
                ii = wrap(ii, width);

            for (int j = -radius; j <= radius; j++)
            {
                int jj = centrey + j;

                if (jj >= 0 && jj <= height)
                {
                    if (ergprobability[ii][jj] == 0 && i * i + j * j < radius * radius + radius)
                        ergprobability[ii][jj] = 1;
                }
            }
        }
        */
    }

    int mindist = 2; // Minimum distance to existing lakes etc.

    int lakeno = 0;

    vector<vector<int>> thislake(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.climate(i, j) == 5 && ergprobability[i][j] != 0) // If this is a hot desert
            {
                if (random(1, ergprobability[i][j]) == 1)
                {
                    bool goahead = 1;

                    for (int k = i - mindist; k <= i + mindist; k++) // Check that there aren't any lakes too close.
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - mindist; l <= j + mindist; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.truelake(kk, l) != 0)
                                    goahead = 0;
                            }
                        }
                    }

                    if (goahead == 1)
                    {
                        bool canoverlap = 0;

                        if (random(1, 100) < overlapprobability[i][j])
                            canoverlap = 1;
                        
                        int shapenumber;

                        if (random(1, 4) == 1)
                        {
                            shapenumber = random(0, 11);
                            drawspeciallake(world, shapenumber, i, j, lakeno, thislake, smalllake, canoverlap, 120); // 120 is the code for ergs.

                        }
                        else
                        {
                            if (random(1, 2) == 1)
                                shapenumber = random(0, 2);
                            else
                                shapenumber = random(0, 9);
                            drawspeciallake(world, shapenumber, i, j, lakeno, thislake, largelake, canoverlap, 120); // 120 is the code for ergs.
                        }

                        lakeno++;
                    }
                }
            }
        }
    }

    // Some ergs may have overlapped, so now we need to make sure their elevation remains constant.

    checkergelevation(world);
}

// This function puts salt pans in deserts.

void createsaltpans(planet& world, boolshapetemplate smalllake[], boolshapetemplate largelake[])
{
    int width = world.width();
    int height = world.height();

    int saltchance = 1000; //250; // Probability of salt pans in any given desert tile.
    int mindist = 2; // Minimum distance to existing lakes etc.

    int lakeno = 0;

    vector<vector<int>> thislake(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.climate(i, j) == 5) // If this is a hot desert
            {
                if (random(1, saltchance) == 1)
                {
                    bool goahead = 1;

                    for (int k = i - mindist; k <= i + mindist; k++) // Check that there aren't any other lakes too close.
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - mindist; l <= j + mindist; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.lakesurface(kk, l) != 0)
                                    goahead = 0;
                            }
                        }
                    }

                    if (goahead == 1)
                    {
                        int shapenumber;

                        if (random(1, 4) != 1)
                        {
                            shapenumber = random(0, 3);
                            drawspeciallake(world, shapenumber, i, j, lakeno, thislake, smalllake, 0, 110); // 110 is the code for salt pans.
                        }
                        else
                        {
                            shapenumber = random(3, 11);
                            drawspeciallake(world, shapenumber, i, j, lakeno, thislake, smalllake, 0, 110); // 110 is the code for salt pans.
                        }

                        lakeno++;
                    }
                }
            }
        }
    }
}

// This refines the roughness map so that it takes into account more factors, producing a roughness factor that will be used for the fractal terrain in the regional map.

void refineroughnessmap(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();
    float maxvaluemod = 1.0f; // Maximum that a valuemod can be.
    float maxcoastalvaluemod = 0.1f; //0.08; // Maximum that it can at the coasts.

    // This array will hold the valuemod for every tile on the map. The valuemod is used in the regional map to determine how rough the diamond-square routine makes the terrain.

    vector<vector<float>> valuemod(ARRAYWIDTH, vector<float>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.sea(i, j) == 0)
            {
                int eax = i + 1;
                int eay = j;
                if (eax > width)
                    eax = 0;

                int wex = i - 1;
                int wey = j;
                if (wex < 0)
                    wex = width;

                int nox = i;
                int noy = j - 1;
                if (noy < 0)
                    noy = 0;

                int sox = i;
                int soy = j + 1;
                if (soy > height)
                    soy = height;

                int sex = eax;
                int sey = soy;
                /*
                int aveheight = (world.nom(i, j) + world.nom(eax, eay) + world.nom(sox, soy) + world.nom(sex, sey)) / 4;
                */
                if (world.nom(i, j) <= sealevel && world.nom(nox, noy) <= sealevel && world.nom(sox, soy) <= sealevel && world.nom(eax, eay) <= sealevel && world.nom(wex, wey) <= sealevel) // Sea beds will be very smooth.
                {
                    float m = 0.08f / (float)maxelev;
                    float n = world.roughness(i, j);
                    float extrabit = 0.01f * (float)random(1, 7);

                    valuemod[i][j] = (m * n) + extrabit;
                }
                else // On land, the higher it is, the rougher it is.
                {
                    valuemod[i][j] = world.roughness(i, j) / (float)maxelev;
                    /*
                    float mult = 3.0f / (float)(maxelev - (sealevel - 50)); // 6.0f
                    float rough = world.roughness(i, j);
                    rough = rough / (float)maxelev;

                    valuemod[i][j] = ((float)aveheight - ((float)sealevel - 50.0f)) * mult * rough;

                    // Alter it according to how flat this whole area is

                    float diff = 0.0f;

                    diff = diff + (float)abs(world.map(i, j) - world.map(eax, eay));
                    diff = diff + (float)abs(world.map(sox, soy) - world.map(sex, sey));
                    diff = diff + (float)abs(world.map(i, j) - world.map(sox, soy));
                    diff = diff + (float)abs(world.map(eax, eay) - world.map(sex, sey));

                    diff = diff / 20.0f;

                    // Alter it according to how mountainous the area is.

                    float mountain = (float)world.mountainheight(i, j);

                    if (mountain == 0.0f)
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
                                    if (world.mountainheight(kk, l) > 0)
                                        mountain = (float)world.mountainheight(kk, l) * 0.7f;
                                }
                            }
                        }

                        if (mountain == 0.0f)
                        {
                            for (int k = i - 3; k <= i + 3; k++)
                            {
                                int kk = k;

                                if (kk<0 || kk>width)
                                    kk = wrap(kk, width);

                                for (int l = j - 3; l <= j + 3; l++)
                                {
                                    if (l >= 0 && l <= height)
                                    {
                                        if (world.mountainheight(kk, l) > 0)
                                            mountain = (float)world.mountainheight(kk, l) * 0.3f;
                                    }
                                }
                            }
                        }
                    }

                    mountain = mountain / 100.0f; // 50.0f;

                    diff = diff + mountain;

                    if (diff > 8.0f)
                        diff = 8.0f;

                    valuemod[i][j] = valuemod[i][j] * diff;
                    */
                }
            }
            else // Sea beds - very smooth.
            {

                float m = 0.08f / (float)maxelev;
                float n = world.roughness(i, j);
                float extrabit = 0.01f * (float)random(1, 7);

                valuemod[i][j] = (m * n) + extrabit;
            }
        }
    }

    // Now we blur that, so that there aren't sharp transitions in roughness on the regional map.

    /*
    int amount = 2; // Amount to blur by.

    vector<vector<float>> valuemod2(ARRAYWIDTH, vector<float>(ARRAYHEIGHT, 0));

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (world.nom(i, j) > sealevel && world.truelake(i, j) == 0)
            {
                bool goahead = 1;

                for (int k = i - 1; k <= i + 1; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        if (l >= 0 && l <= height)
                        {
                            if (world.sea(kk, l) == 1)
                            {
                                goahead = 0;
                                k = i + 1;
                                l = j + 1;
                            }
                        }
                    }
                }

                if (goahead == 1)
                {
                    float total = 0.0f;
                    float crount = 0.0f;

                    for (int k = i - amount; k <= i + amount; k++)
                    {
                        int kk = k;

                        if (kk<0 || kk>width)
                            kk = wrap(kk, width);

                        for (int l = j - amount; l <= j + amount; l++)
                        {
                            if (l >= 0 && l <= height)
                            {
                                if (world.nom(kk, l) > sealevel && world.truelake(kk, l) == 0)
                                {
                                    total = total + valuemod[kk][l];
                                    crount++;
                                }
                            }
                        }
                    }
                    valuemod2[i][j] = total / crount;

                    if (valuemod2[i][j] > maxvaluemod)
                        valuemod2[i][j] = maxvaluemod;
                }
                else
                {
                    valuemod2[i][j] = valuemod[i][j];

                    if (valuemod2[i][j] > maxcoastalvaluemod)
                        valuemod2[i][j] = maxcoastalvaluemod;
                }
            }
            else
                valuemod2[i][j] = 0.0f;
        }
    }
    */

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            if (valuemod[i][j] < 0.0f)
                world.setroughness(i, j, 0.0f);
            else
                world.setroughness(i, j, valuemod[i][j]);
        }
    }
}
