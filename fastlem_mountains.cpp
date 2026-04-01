#include "fastlem_mountains.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include "functions.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"

using namespace std;

namespace
{
struct FastLEMEdge
{
    int to = -1;
    float distance = 1.0f;
    bool boundary = false;
};

struct FastLEMBlock
{
    bool land = false;
    bool coastal = false;
    int x = 0;
    int y = 0;
    int gx = 0;
    int gy = 0;
    int landTiles = 0;
    float meanNomNormalised = 0.0f;
    float meanFractalNormalised = 0.0f;
    int owner = -1;
};

struct FastLEMSite
{
    int x = 0;
    int y = 0;
    int gx = 0;
    int gy = 0;
    float weight = 1.0f;
    int area = 0;
    bool coastalRegion = false;
    bool outlet = false;
    int coastDistance = 0;
    float baseElevation = 0.0f;
    float uplift = 0.0f;
    float erodibility = 1.0f;
    float solvedElevation = 0.0f;
    int peakHeight = 0;
    int divideCount = 0;
    float candidateScore = 0.0f;
    vector<FastLEMEdge> edges;
};

int blockindex(int gx, int gy, int gridwidth)
{
    return gy * gridwidth + gx;
}

uint64_t mixbits(uint64_t value)
{
    value ^= value >> 33;
    value *= 0xff51afd7ed558ccdULL;
    value ^= value >> 33;
    value *= 0xc4ceb9fe1a85ec53ULL;
    value ^= value >> 33;
    return value;
}

float hash01(long seed, int x, int y, int salt)
{
    const uint64_t bits = mixbits(static_cast<uint64_t>(seed) ^ (static_cast<uint64_t>(x) << 21) ^ (static_cast<uint64_t>(y) << 42) ^ static_cast<uint64_t>(salt) * 0x9e3779b97f4a7c15ULL);
    return static_cast<float>(bits & 0xFFFFFF) / static_cast<float>(0x1000000);
}

float wrappedxdistance(int x1, int x2, int width)
{
    const int span = width + 1;
    int delta = abs(x1 - x2);
    delta = min(delta, span - delta);
    return static_cast<float>(delta);
}

float wrappedblockdistance(int gx1, int gx2, int gridwidth)
{
    const int span = gridwidth;
    int delta = abs(gx1 - gx2);
    delta = min(delta, span - delta);
    return static_cast<float>(delta);
}

float worlddistance(int x1, int y1, int x2, int y2, int width)
{
    const float dx = wrappedxdistance(x1, x2, width);
    const float dy = static_cast<float>(y1 - y2);
    return sqrtf(dx * dx + dy * dy);
}

float blockdistance(int gx1, int gy1, int gx2, int gy2, int gridwidth)
{
    const float dx = wrappedblockdistance(gx1, gx2, gridwidth);
    const float dy = static_cast<float>(gy1 - gy2);
    return sqrtf(dx * dx + dy * dy);
}

bool iscoastaltile(planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();

    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            if (dx == 0 && dy == 0)
                continue;

            const int yy = y + dy;

            if (yy < 0 || yy > height)
                return true;

            const int xx = wrap(x + dx, width);

            if (world.sea(xx, yy) == 1)
                return true;
        }
    }

    return false;
}

void addedge(vector<FastLEMSite>& sites, int from, int to, float distance, bool boundary)
{
    if (from == to || from < 0 || to < 0)
        return;

    bool foundforward = false;

    for (FastLEMEdge& edge : sites[from].edges)
    {
        if (edge.to == to)
        {
            edge.boundary = edge.boundary || boundary;
            foundforward = true;
            break;
        }
    }

    bool foundreverse = false;

    for (FastLEMEdge& edge : sites[to].edges)
    {
        if (edge.to == from)
        {
            edge.boundary = edge.boundary || boundary;
            foundreverse = true;
            break;
        }
    }

    if (foundforward == false)
        sites[from].edges.push_back({ to, distance, boundary });

    if (foundreverse == false)
        sites[to].edges.push_back({ from, distance, boundary });
}

void drawrawline(vector<vector<int>>& rawmountains, int width, int height, int x1, int y1, int x2, int y2, int value)
{
    const int span = width + 1;
    int drawx1 = x1;
    int drawx2 = x2;
    int dx = drawx2 - drawx1;

    if (abs(dx) > span / 2)
    {
        if (dx > 0)
            drawx1 += span;
        else
            drawx2 += span;
    }

    const int steps = max(abs(drawx2 - drawx1), abs(y2 - y1));

    if (steps == 0)
    {
        rawmountains[wrap(drawx1, width)][y1] = max(rawmountains[wrap(drawx1, width)][y1], value);
        return;
    }

    for (int step = 0; step <= steps; step++)
    {
        const float t = static_cast<float>(step) / static_cast<float>(steps);
        const int xx = wrap(static_cast<int>(roundf(drawx1 + (drawx2 - drawx1) * t)), width);
        const int yy = static_cast<int>(roundf(y1 + (y2 - y1) * t));

        if (yy >= 0 && yy <= height)
            rawmountains[xx][yy] = max(rawmountains[xx][yy], value);
    }
}

int findnearestlandblock(const vector<FastLEMBlock>& blocks, int gridwidth, int gridheight, int gx, int gy, int maxradius)
{
    gx = wrap(gx, gridwidth - 1);
    gy = clamp(gy, 0, gridheight - 1);

    for (int radius = 0; radius <= maxradius; radius++)
    {
        int bestindex = -1;
        float bestdistance = numeric_limits<float>::max();

        for (int dy = -radius; dy <= radius; dy++)
        {
            for (int dx = -radius; dx <= radius; dx++)
            {
                if (max(abs(dx), abs(dy)) != radius)
                    continue;

                const int ngx = wrap(gx + dx, gridwidth - 1);
                const int ngy = gy + dy;

                if (ngy < 0 || ngy >= gridheight)
                    continue;

                const int index = blockindex(ngx, ngy, gridwidth);

                if (blocks[index].land == false)
                    continue;

                const float distance = blockdistance(gx, gy, ngx, ngy, gridwidth);

                if (distance < bestdistance)
                {
                    bestdistance = distance;
                    bestindex = index;
                }
            }
        }

        if (bestindex != -1)
            return bestindex;
    }

    return -1;
}

bool hassitenear(const vector<FastLEMSite>& sites, int gx, int gy, int gridwidth, float minimumdistance)
{
    for (const FastLEMSite& site : sites)
    {
        if (blockdistance(gx, gy, site.gx, site.gy, gridwidth) < minimumdistance)
            return true;
    }

    return false;
}

void addsitefromblock(vector<FastLEMSite>& sites, const FastLEMBlock& block, long seed)
{
    FastLEMSite site;
    site.x = block.x;
    site.y = block.y;
    site.gx = block.gx;
    site.gy = block.gy;
    site.baseElevation = block.meanNomNormalised;
    site.solvedElevation = block.meanNomNormalised * 0.25f + block.meanFractalNormalised * 0.02f + hash01(seed, block.gx, block.gy, 31) * 0.0001f;

    const float weightnoise = block.meanFractalNormalised * 0.55f + hash01(seed, block.gx, block.gy, 29) * 0.45f;
    site.weight = tuning::terrain::fastlem::siteWeightMin +
        (tuning::terrain::fastlem::siteWeightMax - tuning::terrain::fastlem::siteWeightMin) * clamp(weightnoise, 0.0f, 1.0f);

    sites.push_back(site);
}

int ensureminimumoutlets(vector<FastLEMSite>& sites, vector<int>& coastalids, int width)
{
    sort(coastalids.begin(), coastalids.end(), [&](int left, int right)
    {
        if (sites[left].x != sites[right].x)
            return sites[left].x < sites[right].x;

        return sites[left].y < sites[right].y;
    });

    const float minimumdistance = static_cast<float>(tuning::terrain::fastlem::outletMinimumDistanceBlocks * tuning::terrain::fastlem::cellSize);
    int outletcount = 0;

    for (int attempt = 0; attempt < 2; attempt++)
    {
        const float targetdistance = attempt == 0 ? minimumdistance : minimumdistance * 0.5f;

        for (int siteid : coastalids)
        {
            if (sites[siteid].outlet)
                continue;

            bool fartarget = true;

            for (const FastLEMSite& site : sites)
            {
                if (site.outlet == false)
                    continue;

                if (worlddistance(site.x, site.y, sites[siteid].x, sites[siteid].y, width) < targetdistance)
                {
                    fartarget = false;
                    break;
                }
            }

            if (fartarget)
            {
                sites[siteid].outlet = true;
                outletcount++;
            }
        }

        if (outletcount >= max(3, static_cast<int>(coastalids.size() / 20)))
            break;
    }

    if (outletcount == 0 && coastalids.empty() == false)
    {
        const int step = max(1, static_cast<int>(coastalids.size() / 6));

        for (size_t index = 0; index < coastalids.size(); index += step)
        {
            sites[coastalids[index]].outlet = true;
            outletcount++;
        }
    }

    return outletcount;
}

}

bool generatefastlemmountains(planet& world, const vector<vector<int>>& fractal)
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int sealevel = world.sealevel();
    const long seed = world.seed();

    const int cellsize = tuning::terrain::fastlem::cellSize;
    const int gridwidth = width / cellsize + 1;
    const int gridheight = height / cellsize + 1;
    const int basestep = tuning::terrain::fastlem::baseLatticeStep;

    vector<FastLEMBlock> blocks(gridwidth * gridheight);

    for (int gy = 0; gy < gridheight; gy++)
    {
        const int ystart = gy * cellsize;
        const int yend = min(height, ystart + cellsize - 1);

        for (int gx = 0; gx < gridwidth; gx++)
        {
            const int xstart = gx * cellsize;
            const int xend = min(width, xstart + cellsize - 1);
            FastLEMBlock& block = blocks[blockindex(gx, gy, gridwidth)];
            block.gx = gx;
            block.gy = gy;
            block.x = xstart;
            block.y = ystart;

            int nomsum = 0;
            int fractalsum = 0;
            int bestnom = numeric_limits<int>::min();

            for (int y = ystart; y <= yend; y++)
            {
                for (int x = xstart; x <= xend; x++)
                {
                    if (world.sea(x, y) == 1)
                        continue;

                    block.land = true;
                    block.landTiles++;
                    nomsum += world.nom(x, y);
                    fractalsum += fractal[x][y];

                    if (world.nom(x, y) > bestnom)
                    {
                        bestnom = world.nom(x, y);
                        block.x = x;
                        block.y = y;
                    }

                    if (block.coastal == false && iscoastaltile(world, x, y))
                        block.coastal = true;
                }
            }

            if (block.land)
            {
                block.meanNomNormalised = static_cast<float>(max(0, nomsum / max(1, block.landTiles) - sealevel)) / static_cast<float>(max(1, maxelev - sealevel));
                block.meanFractalNormalised = static_cast<float>(fractalsum / max(1, block.landTiles)) / static_cast<float>(max(1, maxelev));
            }
        }
    }

    vector<FastLEMSite> sites;
    sites.reserve((gridwidth / max(1, basestep) + 1) * (gridheight / max(1, basestep) + 1));

    for (int basegy = 0; basegy < gridheight; basegy += basestep)
    {
        for (int basegx = 0; basegx < gridwidth; basegx += basestep)
        {
            int gx = basegx + static_cast<int>(hash01(seed, basegx, basegy, 1) * 3.0f) - 1;
            int gy = basegy + static_cast<int>(hash01(seed, basegx, basegy, 2) * 3.0f) - 1;

            gx = wrap(gx, gridwidth - 1);
            gy = clamp(gy, 0, gridheight - 1);

            const int nearestblock = findnearestlandblock(blocks, gridwidth, gridheight, gx, gy, tuning::terrain::fastlem::fallbackNeighbourRadius);

            if (nearestblock == -1)
                continue;

            const FastLEMBlock& block = blocks[nearestblock];

            if (block.landTiles < tuning::terrain::fastlem::minimumLandTilesPerSite)
                continue;

            if (hassitenear(sites, block.gx, block.gy, gridwidth, static_cast<float>(basestep) * 0.75f))
                continue;

            addsitefromblock(sites, block, seed);
        }
    }

    for (const FastLEMBlock& block : blocks)
    {
        if (block.land == false || block.coastal)
            continue;

        const float clusterstrength = block.meanFractalNormalised * 0.5f + hash01(seed, block.gx, block.gy, 11) * 0.5f;

        if (clusterstrength < tuning::terrain::fastlem::clusterNoiseThreshold)
            continue;

        if (hash01(seed, block.gx, block.gy, 13) > 0.28f)
            continue;

        if (hassitenear(sites, block.gx, block.gy, gridwidth, static_cast<float>(basestep) * 0.55f))
            continue;

        addsitefromblock(sites, block, seed);
    }

    if (sites.size() < static_cast<size_t>(tuning::terrain::fastlem::minimumCandidateSites))
        return false;

    for (FastLEMBlock& block : blocks)
    {
        if (block.land == false)
            continue;

        int bestowner = -1;
        float bestscore = numeric_limits<float>::max();

        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        {
            const FastLEMSite& site = sites[siteindex];
            const float dx = wrappedxdistance(block.x, site.x, width);
            const float dy = static_cast<float>(block.y - site.y);
            const float distance = sqrtf(dx * dx + dy * dy);
            const float weighted = distance / max(0.1f, site.weight);

            if (weighted < bestscore)
            {
                bestscore = weighted;
                bestowner = siteindex;
            }
        }

        block.owner = bestowner;

        if (bestowner != -1)
        {
            sites[bestowner].area += max(1, block.landTiles);

            if (block.coastal)
                sites[bestowner].coastalRegion = true;
        }
    }

    const int neighbourdirs[8][2] =
    {
        { 1, 0 }, { 0, 1 }, { 1, 1 }, { -1, 1 },
        { -1, 0 }, { 0, -1 }, { -1, -1 }, { 1, -1 }
    };

    for (const FastLEMBlock& block : blocks)
    {
        if (block.land == false || block.owner == -1)
            continue;

        for (int n = 0; n < 4; n++)
        {
            const int ngx = wrap(block.gx + neighbourdirs[n][0], gridwidth - 1);
            const int ngy = block.gy + neighbourdirs[n][1];

            if (ngy < 0 || ngy >= gridheight)
                continue;

            const FastLEMBlock& other = blocks[blockindex(ngx, ngy, gridwidth)];

            if (other.land == false || other.owner == -1 || other.owner == block.owner)
                continue;

            const float distance = worlddistance(block.x, block.y, other.x, other.y, width);
            addedge(sites, block.owner, other.owner, max(1.0f, distance), true);
        }
    }

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        while (static_cast<int>(sites[siteindex].edges.size()) < tuning::terrain::fastlem::minimumConnections)
        {
            int nearest = -1;
            float nearestdistance = numeric_limits<float>::max();

            for (int otherindex = 0; otherindex < static_cast<int>(sites.size()); otherindex++)
            {
                if (otherindex == siteindex)
                    continue;

                bool connected = false;

                for (const FastLEMEdge& edge : sites[siteindex].edges)
                {
                    if (edge.to == otherindex)
                    {
                        connected = true;
                        break;
                    }
                }

                if (connected)
                    continue;

                const float distance = worlddistance(sites[siteindex].x, sites[siteindex].y, sites[otherindex].x, sites[otherindex].y, width);

                if (distance < nearestdistance)
                {
                    nearestdistance = distance;
                    nearest = otherindex;
                }
            }

            if (nearest == -1)
                break;

            addedge(sites, siteindex, nearest, nearestdistance, false);
        }
    }

    vector<int> coastalids;
    coastalids.reserve(sites.size());

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        if (sites[siteindex].coastalRegion)
            coastalids.push_back(siteindex);
    }

    if (ensureminimumoutlets(sites, coastalids, width) == 0)
        return false;

    const int infdistance = numeric_limits<int>::max() / 4;
    queue<int> bfsqueue;
    int maxcoastdistance = 0;

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        sites[siteindex].coastDistance = infdistance;

        if (sites[siteindex].outlet)
        {
            sites[siteindex].coastDistance = 0;
            bfsqueue.push(siteindex);
        }
    }

    while (!bfsqueue.empty())
    {
        const int current = bfsqueue.front();
        bfsqueue.pop();

        for (const FastLEMEdge& edge : sites[current].edges)
        {
            if (sites[edge.to].coastDistance > sites[current].coastDistance + 1)
            {
                sites[edge.to].coastDistance = sites[current].coastDistance + 1;
                maxcoastdistance = max(maxcoastdistance, sites[edge.to].coastDistance);
                bfsqueue.push(edge.to);
            }
        }
    }

    for (FastLEMSite& site : sites)
    {
        if (site.coastDistance == infdistance)
        {
            site.outlet = true;
            site.coastDistance = 0;
        }

        const float inlandnormalised = maxcoastdistance > 0 ? static_cast<float>(site.coastDistance) / static_cast<float>(maxcoastdistance) : 0.0f;
        const float reliefnormalised = site.baseElevation;
        const float fractalnormalised = static_cast<float>(fractal[site.x][site.y]) / static_cast<float>(max(1, maxelev));
        const float signednoise = (fractalnormalised - 0.5f) * 2.0f;
        const float latitude = fabsf(((static_cast<float>(site.y) / static_cast<float>(max(1, height))) * 2.0f) - 1.0f);
        const float equatorialhumidity = 1.0f - latitude;
        const float midlatitudehumidity = 1.0f - min(1.0f, fabsf(latitude - 0.45f) / 0.35f);
        const float precipitationpotential = clamp((1.0f - inlandnormalised) * 0.55f + max(equatorialhumidity, midlatitudehumidity * 0.65f) * 0.35f + max(0.0f, signednoise) * 0.10f, 0.0f, 1.0f);

        site.uplift =
            tuning::terrain::fastlem::baseUplift +
            powf(inlandnormalised, 1.10f) * tuning::terrain::fastlem::inlandUplift +
            reliefnormalised * tuning::terrain::fastlem::reliefUplift +
            max(0.0f, signednoise) * tuning::terrain::fastlem::noiseUplift;

        site.erodibility =
            tuning::terrain::fastlem::baseErodibility +
            precipitationpotential * tuning::terrain::fastlem::coastalErodibility +
            (1.0f - reliefnormalised) * tuning::terrain::fastlem::noiseErodibility;

        site.erodibility = clamp(site.erodibility, tuning::terrain::fastlem::minimumErodibility, tuning::terrain::fastlem::maximumErodibility);
    }

    vector<float> elevations(sites.size(), 0.0f);
    vector<int> nextsite(sites.size(), 0);
    vector<float> nextdistance(sites.size(), 1.0f);
    const float maxslopetan = tanf(tuning::terrain::fastlem::maxSlopeRadians);

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        elevations[siteindex] = sites[siteindex].solvedElevation;

    for (int iteration = 0; iteration < tuning::terrain::fastlem::iterations; iteration++)
    {
        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        {
            FastLEMSite& site = sites[siteindex];

            if (site.outlet || site.edges.empty())
            {
                nextsite[siteindex] = siteindex;
                nextdistance[siteindex] = 1.0f;
                continue;
            }

            int bestdownhill = siteindex;
            float bestdownhillslope = 0.0f;
            float bestdownhilldistance = 1.0f;

            int bestescape = -1;
            int bestescapecoast = numeric_limits<int>::max();
            float bestescapeelevation = numeric_limits<float>::max();
            float bestescapedistance = 1.0f;

            int fallback = -1;
            float fallbackelevation = numeric_limits<float>::max();
            float fallbackdistance = 1.0f;

            for (const FastLEMEdge& edge : site.edges)
            {
                const int other = edge.to;
                const float otherheight = elevations[other];

                if (elevations[siteindex] > otherheight + 0.000001f)
                {
                    const float slope = (elevations[siteindex] - otherheight) / edge.distance;

                    if (slope > bestdownhillslope)
                    {
                        bestdownhillslope = slope;
                        bestdownhill = other;
                        bestdownhilldistance = edge.distance;
                    }
                }

                if (sites[other].coastDistance < site.coastDistance)
                {
                    if (sites[other].coastDistance < bestescapecoast ||
                        (sites[other].coastDistance == bestescapecoast && otherheight < bestescapeelevation))
                    {
                        bestescape = other;
                        bestescapecoast = sites[other].coastDistance;
                        bestescapeelevation = otherheight;
                        bestescapedistance = edge.distance;
                    }
                }

                if (otherheight < fallbackelevation)
                {
                    fallback = other;
                    fallbackelevation = otherheight;
                    fallbackdistance = edge.distance;
                }
            }

            if (bestdownhill != siteindex)
            {
                nextsite[siteindex] = bestdownhill;
                nextdistance[siteindex] = bestdownhilldistance;
            }
            else if (bestescape != -1)
            {
                nextsite[siteindex] = bestescape;
                nextdistance[siteindex] = bestescapedistance;
            }
            else if (fallback != -1)
            {
                nextsite[siteindex] = fallback;
                nextdistance[siteindex] = fallbackdistance;
            }
            else
            {
                nextsite[siteindex] = siteindex;
                nextdistance[siteindex] = 1.0f;
            }
        }

        vector<vector<int>> children(sites.size());

        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        {
            if (nextsite[siteindex] != siteindex)
                children[nextsite[siteindex]].push_back(siteindex);
        }

        vector<int> order;
        order.reserve(sites.size());
        vector<bool> visited(sites.size(), false);
        vector<int> stack;

        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        {
            if (nextsite[siteindex] != siteindex)
                continue;

            stack.push_back(siteindex);

            while (!stack.empty())
            {
                const int current = stack.back();
                stack.pop_back();

                if (visited[current])
                    continue;

                visited[current] = true;
                order.push_back(current);

                for (int child : children[current])
                    stack.push_back(child);
            }
        }

        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        {
            if (visited[siteindex])
                continue;

            stack.push_back(siteindex);

            while (!stack.empty())
            {
                const int current = stack.back();
                stack.pop_back();

                if (visited[current])
                    continue;

                visited[current] = true;
                order.push_back(current);

                for (int child : children[current])
                    stack.push_back(child);
            }
        }

        vector<float> drainageareas(sites.size(), 1.0f);
        vector<float> responsetimes(sites.size(), 0.0f);

        for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
            drainageareas[siteindex] = static_cast<float>(max(1, sites[siteindex].area));

        for (int orderindex = static_cast<int>(order.size()) - 1; orderindex >= 0; orderindex--)
        {
            const int current = order[orderindex];
            const int downstream = nextsite[current];

            if (downstream != current)
                drainageareas[downstream] += drainageareas[current];
        }

        vector<float> newelevations = elevations;
        float maxdelta = 0.0f;

        for (int current : order)
        {
            const int downstream = nextsite[current];

            if (downstream == current)
            {
                newelevations[current] = 0.0f;
                continue;
            }

            const float celerity = max(tuning::terrain::fastlem::minimumErodibility, sites[current].erodibility) * powf(max(1.0f, drainageareas[current]), tuning::terrain::fastlem::mExponent);
            responsetimes[current] = responsetimes[downstream] + nextdistance[current] / celerity;

            float targetelevation = sites[current].baseElevation + sites[current].uplift * responsetimes[current];
            const float maxslopeelevation = newelevations[downstream] + maxslopetan * nextdistance[current];

            if (targetelevation > maxslopeelevation)
                targetelevation = maxslopeelevation;

            if (targetelevation <= newelevations[downstream])
                targetelevation = newelevations[downstream] + 0.0001f;

            maxdelta = max(maxdelta, fabsf(targetelevation - elevations[current]));
            newelevations[current] = targetelevation;
        }

        elevations.swap(newelevations);

        if (maxdelta < 0.0005f)
            break;
    }

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
        sites[siteindex].solvedElevation = elevations[siteindex];

    vector<int> basinroots(sites.size(), -1);

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        int current = siteindex;
        vector<int> trail;

        while (nextsite[current] != current && basinroots[current] == -1)
        {
            trail.push_back(current);
            current = nextsite[current];

            if (trail.size() > sites.size())
                break;
        }

        if (basinroots[current] == -1)
            basinroots[current] = current;

        for (int trailsite : trail)
            basinroots[trailsite] = basinroots[current];
    }

    float maxsolvedelevation = 0.0f;

    for (const FastLEMSite& site : sites)
        maxsolvedelevation = max(maxsolvedelevation, site.solvedElevation);

    if (maxsolvedelevation <= 0.0f)
        return false;

    vector<bool> preliminarycandidate(sites.size(), false);

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        FastLEMSite& site = sites[siteindex];
        const float elevationnormalised = site.solvedElevation / maxsolvedelevation;
        const float coastnormalised = maxcoastdistance > 0 ? static_cast<float>(site.coastDistance) / static_cast<float>(maxcoastdistance) : 0.0f;

        bool localmaximum = true;
        int dividecount = 0;

        for (const FastLEMEdge& edge : site.edges)
        {
            const FastLEMSite& other = sites[edge.to];

            if (other.solvedElevation > site.solvedElevation + maxsolvedelevation * 0.02f)
                localmaximum = false;

            if (basinroots[edge.to] != basinroots[siteindex])
                dividecount++;
        }

        site.divideCount = dividecount;
        const float dividestrength = clamp(static_cast<float>(dividecount) / 3.0f, 0.0f, 1.0f);
        const float score =
            elevationnormalised * tuning::terrain::fastlem::elevationScoreWeight +
            coastnormalised * tuning::terrain::fastlem::coastScoreWeight +
            dividestrength * 0.20f;

        site.candidateScore = score + dividestrength * 0.25f + elevationnormalised * 0.15f;

        if ((dividecount > 0 &&
            site.coastDistance >= tuning::terrain::fastlem::minimumRidgeCoastDistance &&
            elevationnormalised >= tuning::terrain::fastlem::minimumRidgeElevationNormalised &&
            score >= tuning::terrain::fastlem::minimumRidgeScore) ||
            (localmaximum &&
                site.coastDistance >= tuning::terrain::fastlem::minimumPeakCoastDistance &&
                elevationnormalised >= tuning::terrain::fastlem::minimumPeakElevationNormalised))
        {
            preliminarycandidate[siteindex] = true;
        }

        if (preliminarycandidate[siteindex])
        {
            const float heightnormalised = clamp(elevationnormalised * 0.65f + coastnormalised * 0.20f + dividestrength * 0.15f, 0.0f, 1.0f);
            site.peakHeight = tuning::terrain::fastlem::minimumPeakHeight +
                static_cast<int>(powf(heightnormalised, tuning::terrain::fastlem::peakHeightExponent) *
                    static_cast<float>(tuning::terrain::fastlem::maximumPeakHeight - tuning::terrain::fastlem::minimumPeakHeight));
        }
    }

    vector<int> candidateorder;
    candidateorder.reserve(sites.size());

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        if (preliminarycandidate[siteindex])
            candidateorder.push_back(siteindex);
    }

    sort(candidateorder.begin(), candidateorder.end(), [&](int left, int right)
    {
        if (sites[left].candidateScore != sites[right].candidateScore)
            return sites[left].candidateScore > sites[right].candidateScore;

        if (sites[left].peakHeight != sites[right].peakHeight)
            return sites[left].peakHeight > sites[right].peakHeight;

        return left < right;
    });

    vector<bool> candidate(sites.size(), false);
    vector<int> keptcandidates;
    keptcandidates.reserve(candidateorder.size());
    const float minimumcandidatespacing = tuning::terrain::fastlem::minimumCandidateSpacingBlocks;

    for (int siteindex : candidateorder)
    {
        bool tooclose = false;

        for (int keptindex : keptcandidates)
        {
            if (blockdistance(sites[siteindex].gx, sites[siteindex].gy, sites[keptindex].gx, sites[keptindex].gy, gridwidth) < minimumcandidatespacing)
            {
                tooclose = true;
                break;
            }
        }

        if (tooclose)
            continue;

        candidate[siteindex] = true;
        keptcandidates.push_back(siteindex);
    }

    int candidatecount = 0;
    int boundarycandidatecount = 0;

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        if (candidate[siteindex] == false)
            continue;

        candidatecount++;

        if (sites[siteindex].divideCount > 0)
            boundarycandidatecount++;
    }

    const int minimumcandidatecount = min(tuning::terrain::fastlem::minimumCandidateSites, max(6, static_cast<int>(sites.size() / 30)));
    const int minimumboundarycount = min(tuning::terrain::fastlem::minimumBoundaryCandidates, max(3, static_cast<int>(sites.size() / 60)));

    if (candidatecount < minimumcandidatecount || boundarycandidatecount < minimumboundarycount)
        return false;

    vector<vector<int>> rawmountains(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        if (candidate[siteindex] == false)
            continue;

        rawmountains[sites[siteindex].x][sites[siteindex].y] = max(rawmountains[sites[siteindex].x][sites[siteindex].y], sites[siteindex].peakHeight);
    }

    const int boundaryneighbours[4][2] =
    {
        { 1, 0 }, { 0, 1 }, { 1, 1 }, { -1, 1 }
    };

    for (const FastLEMBlock& block : blocks)
    {
        if (block.land == false || block.owner == -1)
            continue;

        const FastLEMSite& source = sites[block.owner];
        const float sourcenormalised = source.solvedElevation / maxsolvedelevation;

        for (int n = 0; n < 4; n++)
        {
            const int ngx = wrap(block.gx + boundaryneighbours[n][0], gridwidth - 1);
            const int ngy = block.gy + boundaryneighbours[n][1];

            if (ngy < 0 || ngy >= gridheight)
                continue;

            const FastLEMBlock& otherblock = blocks[blockindex(ngx, ngy, gridwidth)];

            if (otherblock.land == false || otherblock.owner == -1 || otherblock.owner == block.owner)
                continue;

            const FastLEMSite& target = sites[otherblock.owner];
            const float targetnormalised = target.solvedElevation / maxsolvedelevation;
            const bool differentbasin = basinroots[block.owner] != basinroots[otherblock.owner];
            const bool similarelevation = fabsf(source.solvedElevation - target.solvedElevation) <= maxsolvedelevation * 0.15f;
            const bool strongboundary = max(sourcenormalised, targetnormalised) >= tuning::terrain::fastlem::minimumRidgeElevationNormalised;

            if (strongboundary == false)
                continue;

            if (differentbasin || similarelevation)
            {
                const int peak = max(source.peakHeight, target.peakHeight);
                drawrawline(rawmountains, width, height, block.x, block.y, otherblock.x, otherblock.y, max(peak, tuning::terrain::fastlem::minimumPeakHeight));
            }
        }
    }

    for (int siteindex = 0; siteindex < static_cast<int>(sites.size()); siteindex++)
    {
        if (candidate[siteindex] == false)
            continue;

        for (const FastLEMEdge& edge : sites[siteindex].edges)
        {
            if (edge.to <= siteindex)
                continue;

            if (candidate[edge.to] == false)
                continue;

            if (edge.boundary && basinroots[siteindex] != basinroots[edge.to])
            {
                drawrawline(rawmountains, width, height, sites[siteindex].x, sites[siteindex].y, sites[edge.to].x, sites[edge.to].y, max(sites[siteindex].peakHeight, sites[edge.to].peakHeight));
            }
        }
    }

    vector<vector<bool>> dummyok(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, false));

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            world.setmountainheight(x, y, 0);
            world.setmountainridge(x, y, 0);
        }
    }

    createmountainsfromraw(world, rawmountains, dummyok);
    return true;
}
