// Convert raw ridge masks from imports to mountain relief.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the sphere since T8 W2-T2.

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"
#include "simulations/terrain/terrain_coasts.hpp"

using namespace std;
using mesh::CellId;

// createmountainsfromraw on the cell grid (T8 W2-T2): the raw ridges read from their legacy slots
// (terrain::coastgrid::legacyslot), the ridge clean-ups on the neighbour table, the distances
// from the central peaks within the legacy pixels' km, the variation from sphere noise.

namespace
{
namespace cg = terrain::coastgrid;

// Removes the ridge leaving c in slot from the arrays, and its other end.
void deletearrayridge(const mesh::CellGrid& grid, vector<int>& ridges, vector<int>& heights, CellId c, int slot)
{
    const int code = 1 << (slot - 1);
    if ((ridges[c] & code) == 0) // If there isn't a ridge going this way
        return;
    ridges[c] -= code;
    if (ridges[c] == 0)
        heights[c] = 0;

    // Now remove one going the other way.
    const std::int32_t n = grid.neighbour(c, slot);
    if (n < 0)
        return;
    const int back = 1 << (grid.reverse(c, slot) - 1);
    if ((ridges[n] & back) == 0)
        return;
    ridges[n] -= back;
    if (ridges[n] == 0)
        heights[n] = 0;
}

void spheremountainsfromraw(planet& world, vector<vector<int>>& rawmountainsslots, vector<vector<bool>>& OKmountains)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    const double pixelkm = cg::legacypixelkm(world);
    int maxelev = world.maxelevation();

    vector<int> rawmountains(count, 0);
    for (CellId c = 0; c < count; c++)
    {
        int x = 0, y = 0;
        cg::legacyslot(grid, c, x, y);
        rawmountains[c] = rawmountainsslots[static_cast<size_t>(x)][static_cast<size_t>(y)];
    }

    vector<int> extraraw(count, 0); // This will hold extra raw ridges that we add near the start.
    vector<int> mountaindist(count, 0); // This holds the proximity to the central peaks (the higher the value, the closer it is).
    vector<int> mountainbaseheight(count, 0); // The height of the central peak it is measured from.
    vector<char> timesten(count, 0); // Whether the ridge distance has been multiplied by ten (so that the ridges on either side of the main ridge don't join up).
    vector<int> mountainridges(count, 0);
    vector<int> mountainheights(count, 0);

    const vector<int> heightsfractal = cg::spherefractal(world, 0x6d0e0001, cg::finenoise(), maxelev); // This will vary the heights of the peaks in the subsidiary ridges.
    const vector<int> radiusfractal = cg::spherefractal(world, 0x6d0e0002, cg::grain16noise(), maxelev); // This will vary the effective radius of each mountain point.

    float maxradius = 6; // The bigger this is, the wider the mountain ranges will be (pixels).

    const vector<CellId> order = cg::columnorder(grid);
    auto nb = [&](CellId c, int slot) { return grid.neighbour(c, slot); };
    auto raw = [&](std::int32_t n) { return n >= 0 ? rawmountains[n] : 0; };
    auto clear = [&](std::int32_t n) { if (n >= 0) rawmountains[n] = 0; };
    auto interior = [&](CellId c) { return nb(c, 1) >= 0 && nb(c, 5) >= 0; };

    // First, ensure that the raw mountains array has ridges that are only one cell wide.

    for (const CellId c : order)
    {
        if (interior(c) && rawmountains[c] != 0 && raw(nb(c, 1)) != 0 && raw(nb(c, 5)) != 0)
        {
            clear(nb(c, 3));
            clear(nb(c, 7));
        }
    }

    for (const CellId c : order)
    {
        if (interior(c) && rawmountains[c] != 0 && raw(nb(c, 3)) != 0 && raw(nb(c, 7)) != 0)
        {
            clear(nb(c, 5));
            clear(nb(c, 1));
        }
    }

    for (const CellId c : order)
    {
        if (!interior(c))
            continue;

        if (rawmountains[c] != 0 && raw(nb(c, 4)) != 0)
        {
            clear(nb(c, 3));
            clear(nb(c, 5));
        }

        if (raw(nb(c, 3)) != 0 && raw(nb(c, 5)) != 0)
        {
            rawmountains[c] = 0;
            clear(nb(c, 4));
        }
    }

    // Now add some extra ridges! We do this on another array to start with so we don't add them onto each other.

    int extrachance = 10; // The higher this is, the fewer extra ridges there will be.
    int swervechance = 4; // The higher this is, the straighter they will be.
    int minlength = 2;
    int maxlength = 4; // Possible lengths of these extra ridges.

    for (int t = 1; t <= 2; t++) // Do the whole thing twice!
    {
        for (const CellId c : order)
        {
            if (!interior(c) || rawmountains[c] == 0 || random(1, extrachance) != 1)
                continue;

            int dir = 0;
            int x = 0;
            int y = 0; // The start, as an offset (east, south) from c.

            if (raw(nb(c, 1)) != 0)
            {
                x = random(1, 2) == 1 ? -2 : 2;
                y = 0;
                dir = random(1, 2) == 1 ? 1 : 5;
            }

            if (raw(nb(c, 2)) != 0)
            {
                if (random(1, 2) == 1) { x = -1; y = -1; } else { x = 1; y = 1; }
                dir = random(1, 2) == 1 ? 2 : 6;
            }

            if (raw(nb(c, 3)) != 0)
            {
                y = random(1, 2) == 1 ? -2 : 2;
                x = 0;
                dir = random(1, 2) == 1 ? 3 : 7;
            }

            if (raw(nb(c, 4)) != 0)
            {
                if (random(1, 2) == 1) { x = 1; y = -1; } else { x = -1; y = 1; }
                dir = random(1, 2) == 1 ? 4 : 8;
            }

            if (raw(nb(c, 5)) != 0)
            {
                x = random(1, 2) == 1 ? -2 : 2;
                y = 0;
                dir = random(1, 2) == 1 ? 5 : 1;
            }

            if (raw(nb(c, 6)) != 0)
            {
                if (random(1, 2) == 1) { x = -1; y = -1; } else { x = 1; y = 1; }
                dir = random(1, 2) == 1 ? 6 : 2;
            }

            if (raw(nb(c, 7)) != 0)
            {
                y = random(1, 2) == 1 ? -2 : 2;
                x = 0;
                dir = random(1, 2) == 1 ? 7 : 3;
            }

            if (raw(nb(c, 8)) != 0)
            {
                if (random(1, 2) == 1) { x = 1; y = -1; } else { x = -1; y = 1; }
                dir = random(1, 2) == 1 ? 8 : 4;
            }

            if (dir != 0)
            {
                int peakheight = rawmountains[c];
                int origpeakheight = peakheight;

                int length = random(minlength, maxlength);

                std::int32_t at = cg::offset(grid, c, x, y);

                for (int n = 1; n <= length && at >= 0; n++)
                {
                    at = grid.neighbour(static_cast<CellId>(at), dir);

                    if (at < 0)
                        break;

                    extraraw[at] = peakheight;

                    if (random(1, swervechance) == 1)
                    {
                        if (random(1, 2) == 1)
                            dir--;
                        else
                            dir++;

                        if (dir == 0)
                            dir = 8;

                        if (dir == 9)
                            dir = 1;
                    }

                    peakheight = peakheight + randomsign(random(1, 5));

                    if (peakheight < origpeakheight / 2)
                        peakheight = origpeakheight / 2;

                    if (peakheight > origpeakheight)
                        peakheight = origpeakheight;
                }
            }
        }

        // Now just add those extra raw ridges to the main raw array.

        for (CellId c = 0; c < count; c++)
        {
            if (rawmountains[c] < extraraw[c])
                rawmountains[c] = extraraw[c];

            extraraw[c] = 0;
        }
    }

    // Now get rid of straight sections, where possible.

    auto move = [&](std::int32_t to, std::int32_t from)
    {
        if (to < 0 || from < 0)
            return;
        rawmountains[to] = rawmountains[from];
        rawmountains[from] = 0;
    };

    for (int n = 1; n <= 2; n++) // Do this whole thing twice, to be sure.
    {
        for (const CellId c : order)
        {
            if (!interior(c) || rawmountains[c] == 0)
                continue;

            const std::int32_t self = static_cast<std::int32_t>(c);

            // First, east and west.

            if (raw(nb(c, 7)) != 0 && raw(nb(c, 3)) != 0)
            {
                bool up = random(1, 2) == 1;

                move(nb(c, up ? 1 : 5), self);

                if (random(1, 2) == 1) // Maybe move one of the neighbouring ones too
                {
                    const std::int32_t side = random(1, 2) == 1 ? nb(c, 3) : nb(c, 7);

                    if (side >= 0)
                        move(nb(static_cast<CellId>(side), up ? 1 : 5), side);
                }
            }

            // Now, north and south.

            if (raw(nb(c, 1)) != 0 && raw(nb(c, 5)) != 0)
            {
                bool left = random(1, 2) == 1;

                move(nb(c, left ? 7 : 3), self);

                if (random(1, 2) == 1) // Maybe move one of the neighbouring ones too
                {
                    const std::int32_t side = random(1, 2) == 1 ? nb(c, 5) : nb(c, 1);

                    if (side >= 0)
                        move(nb(static_cast<CellId>(side), left ? 7 : 3), side);
                }
            }

            // Now, diagonals from NW to SE.

            if (raw(nb(c, 8)) != 0 && raw(nb(c, 4)) != 0)
            {
                bool left = random(1, 2) == 1;

                const std::int32_t first = nb(c, left ? 7 : 3);
                const std::int32_t second = nb(c, left ? 5 : 1);

                if (first >= 0)
                    rawmountains[first] = rawmountains[c];
                if (second >= 0)
                    rawmountains[second] = rawmountains[c];
                rawmountains[c] = 0;
            }
        }
    }

    // Now, work out the distances from the central peaks: within the disc of maxradius pixels,
    // maxradius minus the squared distance in pixels, times ten on the west and north side.

    for (const CellId c : order)
    {
        if (rawmountains[c] == 0)
            continue;

        const float peakheight = (float)rawmountains[c];
        const mesh::Vec3 p = grid.centre(c);
        const mesh::TangentFrame frame = mesh::tangentframe(p, grid.radius_m());

        for (const CellId n : cg::cellswithinkm(world, p, sqrt(maxradius * maxradius + maxradius) * pixelkm))
        {
            double east = 0.0, north = 0.0;
            cg::totangent(frame, grid.centre(n), east, north);
            const double d2 = (east * east + north * north) / (pixelkm * pixelkm);
            const int dist = (int)maxradius - static_cast<int>(lround(d2));

            if (mountaindist[n] < dist)
            {
                mountaindist[n] = dist;
                mountainbaseheight[n] = (int)peakheight;

                if (east < -0.5 * pixelkm || north > 0.5 * pixelkm)
                {
                    mountaindist[n] = dist * 10;
                    timesten[n] = 1;
                }
            }
        }
    }

    // Now sort out the heights

    float heightmult = 1.0f / maxradius;
    float heightmult2 = 1.0f / maxelev;

    float finalmult = 0.6f; // Reduce it all! Because we'll have raised ground underneath.

    for (CellId c = 0; c < count; c++)
    {
        float thisheight = (float)mountainbaseheight[c];

        float thismountaindist = (float)mountaindist[c];

        if (timesten[c] == 1)
            thismountaindist = thismountaindist / 10.0f;

        if (thismountaindist < maxradius - 1.0f)
            thismountaindist = thismountaindist * 0.8f;

        if (thismountaindist < maxradius - 2.0f)
            thismountaindist = thismountaindist * 0.8f;

        float distmult = (float)radiusfractal[c];
        distmult = distmult * heightmult2;

        if (distmult < 0.95f)
            distmult = 0.95f;

        if (thismountaindist > 2)
            thismountaindist = thismountaindist * distmult;

        thisheight = thisheight * thismountaindist * heightmult;

        float thisheightmult = (float)heightsfractal[c]; // This is to vary it using the fractal
        thisheightmult = thisheightmult * heightmult2;

        if (thisheightmult < 0.95f)
            thisheightmult = 0.95f;

        thisheight = thisheight * thisheightmult;

        thisheight = thisheight * finalmult;

        if (thisheight > (float)mountainbaseheight[c])
            thisheight = (float)mountainbaseheight[c];

        mountainheights[c] = (int)thisheight;
    }

    // Now we draw out the ridges of the mountains, based on similar distances from the central peaks.

    for (const CellId c : order)
    {
        if (mountaindist[c] == 0)
            continue;
        for (int slot = 1; slot <= mesh::slotcount; slot++)
        {
            const std::int32_t n = nb(c, slot);
            if (n >= 0 && mountaindist[n] == mountaindist[c] && (mountainridges[c] & (1 << (slot - 1))) == 0)
                mountainridges[c] += 1 << (slot - 1);
        }
    }

    // Now we want to add some variation to the ridge directions.

    int varchance = 3; // The lower this is, the more variation there will be.
    int deletechance = 2; // The higher this is, the more extra ridges will be deleted.

    auto hasridge = [&](CellId c, int slot) { return (mountainridges[c] & (1 << (slot - 1))) != 0; };

    // Turns a ridge leaving c in slot to dir (one of the two slots beside it).
    auto turnridge = [&](CellId c, int slot, int dir)
    {
        const std::int32_t target = nb(c, dir);
        if (target < 0 || hasridge(c, dir)) // If there isn't one in this direction
            return;

        mountainridges[c] += 1 << (dir - 1);
        const int back = grid.reverse(c, dir);
        if ((mountainridges[target] & (1 << (back - 1))) == 0)
            mountainridges[target] += 1 << (back - 1);

        if (mountainheights[target] == 0)
            mountainheights[target] = mountainheights[c];

        if (random(1, deletechance) != 1)
            deletearrayridge(grid, mountainridges, mountainheights, c, slot);
    };

    for (const CellId c : order)
    {
        if (!interior(c) || mountainridges[c] == 0)
            continue;

        if (hasridge(c, 1) && random(1, varchance) == 1) // Looking north
            turnridge(c, 1, random(1, 2) == 1 ? 2 : 8);
        else if (hasridge(c, 2) && random(1, varchance) == 1) // Looking northeast
            turnridge(c, 2, random(1, 2) == 1 ? 3 : 1);
        else if (hasridge(c, 3) && random(1, varchance) == 1) // Looking east
            turnridge(c, 3, random(1, 2) == 1 ? 4 : 2);
        else if (hasridge(c, 4) && random(1, varchance) == 1) // Looking southeast
            turnridge(c, 4, random(1, 2) == 1 ? 5 : 3);
    }

    // Now apply the new mountain heights and ridges to the world.

    for (CellId c = 0; c < count; c++)
    {
        if (mountainridges[c] != 0)
        {
            world.setmountainridge(c, mountainridges[c]);
            world.setmountainheight(c, static_cast<float>(mountainheights[c]));

            int x = 0, y = 0;
            cg::legacyslot(grid, c, x, y);
            OKmountains[static_cast<size_t>(x)][static_cast<size_t>(y)] = 1; // Mark them as user-added, so we won't tinker with their heights later.
        }
    }
}
}

// This creates mountains from an imported raw mountain map.

void createmountainsfromraw(planet& world, vector<vector<int>>& rawmountains, vector<vector<bool>>& OKmountains)
{
    spheremountainsfromraw(world, rawmountains, OKmountains);

    terrain::detail::removefloatingmountains(world);
    cleanmountainridges(world);
}
