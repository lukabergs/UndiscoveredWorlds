#include "fastlem_mountains.hpp"

#include <algorithm>
#include <cmath>
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
};

struct FastLEMSite
{
    int x = 0;
    int y = 0;
    int gx = 0;
    int gy = 0;
    int landTiles = 0;
    int coastDistance = 0;
    bool outlet = false;
    float baseElevation = 0.0f;
    float uplift = 0.0f;
    float erodibility = 1.0f;
    float solvedElevation = 0.0f;
    int peakHeight = 0;
    vector<FastLEMEdge> edges;
};

float wrappedxdistance(int x1, int x2, int width)
{
    const int span = width + 1;
    int delta = abs(x1 - x2);
    delta = min(delta, span - delta);
    return static_cast<float>(delta);
}

float graphdistance(int x1, int y1, int x2, int y2, int width)
{
    const float dx = wrappedxdistance(x1, x2, width);
    const float dy = static_cast<float>(y1 - y2);
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

void addedge(vector<FastLEMSite>& sites, int from, int to, float distance)
{
    if (from == to || from < 0 || to < 0)
        return;

    for (const FastLEMEdge& edge : sites[from].edges)
    {
        if (edge.to == to)
            return;
    }

    sites[from].edges.push_back({ to, distance });
    sites[to].edges.push_back({ from, distance });
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

int findedgepeak(const FastLEMSite& first, const FastLEMSite& second)
{
    return max(first.peakHeight, second.peakHeight);
}
}

bool generatefastlemmountains(planet& world, const vector<vector<int>>& fractal)
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int sealevel = world.sealevel();

    const int cellsize = tuning::terrain::fastlem::cellSize;
    const int gridwidth = width / cellsize + 1;
    const int gridheight = height / cellsize + 1;

    vector<int> blocktosite(gridwidth * gridheight, -1);
    vector<FastLEMSite> sites;
    sites.reserve(gridwidth * gridheight / 2);

    for (int gy = 0; gy < gridheight; gy++)
    {
        const int ystart = gy * cellsize;
        const int yend = min(height, ystart + cellsize - 1);

        for (int gx = 0; gx < gridwidth; gx++)
        {
            const int xstart = gx * cellsize;
            const int xend = min(width, xstart + cellsize - 1);

            int landtiles = 0;
            int bestx = xstart;
            int besty = ystart;
            int bestnom = numeric_limits<int>::min();
            int nomsum = 0;
            int fractalsum = 0;
            bool coastal = false;

            for (int y = ystart; y <= yend; y++)
            {
                for (int x = xstart; x <= xend; x++)
                {
                    if (world.sea(x, y) == 1)
                        continue;

                    landtiles++;
                    nomsum += world.nom(x, y);
                    fractalsum += fractal[x][y];

                    if (world.nom(x, y) > bestnom)
                    {
                        bestnom = world.nom(x, y);
                        bestx = x;
                        besty = y;
                    }

                    if (coastal == false && iscoastaltile(world, x, y))
                        coastal = true;
                }
            }

            if (landtiles == 0)
                continue;

            if (landtiles < tuning::terrain::fastlem::minimumLandTilesPerSite && coastal == false)
                continue;

            FastLEMSite site;
            site.x = bestx;
            site.y = besty;
            site.gx = gx;
            site.gy = gy;
            site.landTiles = landtiles;
            site.outlet = coastal;
            site.baseElevation = static_cast<float>(max(0, nomsum / landtiles - sealevel)) / static_cast<float>(max(1, maxelev - sealevel));

            const float fractalnormalised = static_cast<float>(fractalsum / landtiles) / static_cast<float>(max(1, maxelev));
            site.solvedElevation = site.baseElevation * 0.1f + fractalnormalised * 0.001f + static_cast<float>(sites.size()) * 0.000001f;

            blocktosite[gy * gridwidth + gx] = static_cast<int>(sites.size());
            sites.push_back(site);
        }
    }

    if (sites.size() < static_cast<size_t>(tuning::terrain::fastlem::minimumCandidateSites))
        return false;

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        const int gx = sites[index].gx;
        const int gy = sites[index].gy;

        for (int radius = 1; radius <= tuning::terrain::fastlem::fallbackNeighbourRadius; radius++)
        {
            for (int dy = -radius; dy <= radius; dy++)
            {
                for (int dx = -radius; dx <= radius; dx++)
                {
                    if (dx == 0 && dy == 0)
                        continue;

                    if (max(abs(dx), abs(dy)) != radius)
                        continue;

                    const int ngx = wrap(gx + dx, gridwidth - 1);
                    const int ngy = gy + dy;

                    if (ngy < 0 || ngy >= gridheight)
                        continue;

                    const int other = blocktosite[ngy * gridwidth + ngx];

                    if (other == -1)
                        continue;

                    const float distance = graphdistance(sites[index].x, sites[index].y, sites[other].x, sites[other].y, width);
                    addedge(sites, index, other, max(1.0f, distance));
                }
            }

            if (static_cast<int>(sites[index].edges.size()) >= tuning::terrain::fastlem::minimumConnections)
                break;
        }

        if (sites[index].edges.empty())
            sites[index].outlet = true;
    }

    const int infdistance = numeric_limits<int>::max() / 4;
    queue<int> bfsqueue;
    int maxcoastdistance = 0;

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        sites[index].coastDistance = infdistance;

        if (sites[index].outlet)
        {
            sites[index].coastDistance = 0;
            bfsqueue.push(index);
        }
    }

    if (bfsqueue.empty())
    {
        sites.front().outlet = true;
        sites.front().coastDistance = 0;
        bfsqueue.push(0);
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

        site.uplift =
            tuning::terrain::fastlem::baseUplift +
            powf(inlandnormalised, 1.15f) * tuning::terrain::fastlem::inlandUplift +
            reliefnormalised * tuning::terrain::fastlem::reliefUplift +
            max(0.0f, signednoise) * tuning::terrain::fastlem::noiseUplift;

        site.erodibility =
            tuning::terrain::fastlem::baseErodibility +
            (1.0f - inlandnormalised) * tuning::terrain::fastlem::coastalErodibility +
            (1.0f - reliefnormalised) * tuning::terrain::fastlem::noiseErodibility;

        site.erodibility = clamp(site.erodibility, tuning::terrain::fastlem::minimumErodibility, tuning::terrain::fastlem::maximumErodibility);
    }

    vector<float> elevations(sites.size(), 0.0f);
    vector<int> nextsite(sites.size(), 0);
    vector<float> nextdistance(sites.size(), 1.0f);
    const float maxslopetan = tanf(tuning::terrain::fastlem::maxSlopeRadians);

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
        elevations[index] = sites[index].solvedElevation;

    for (int iteration = 0; iteration < tuning::terrain::fastlem::iterations; iteration++)
    {
        for (int index = 0; index < static_cast<int>(sites.size()); index++)
        {
            FastLEMSite& site = sites[index];

            if (site.outlet || site.edges.empty())
            {
                nextsite[index] = index;
                nextdistance[index] = 1.0f;
                continue;
            }

            int bestdownhill = index;
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

                if (elevations[index] > otherheight + 0.000001f)
                {
                    const float slope = (elevations[index] - otherheight) / edge.distance;

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

            if (bestdownhill != index)
            {
                nextsite[index] = bestdownhill;
                nextdistance[index] = bestdownhilldistance;
            }
            else if (bestescape != -1)
            {
                nextsite[index] = bestescape;
                nextdistance[index] = bestescapedistance;
            }
            else if (fallback != -1)
            {
                nextsite[index] = fallback;
                nextdistance[index] = fallbackdistance;
            }
            else
            {
                nextsite[index] = index;
                nextdistance[index] = 1.0f;
            }
        }

        vector<vector<int>> children(sites.size());

        for (int index = 0; index < static_cast<int>(sites.size()); index++)
        {
            if (nextsite[index] != index)
                children[nextsite[index]].push_back(index);
        }

        vector<int> order;
        order.reserve(sites.size());
        vector<bool> visited(sites.size(), false);
        vector<int> stack;

        for (int index = 0; index < static_cast<int>(sites.size()); index++)
        {
            if (nextsite[index] != index)
                continue;

            stack.push_back(index);

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

        for (int index = 0; index < static_cast<int>(sites.size()); index++)
        {
            if (visited[index])
                continue;

            stack.push_back(index);

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

        for (int index = 0; index < static_cast<int>(sites.size()); index++)
            drainageareas[index] = static_cast<float>(max(1, sites[index].landTiles));

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

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
        sites[index].solvedElevation = elevations[index];

    vector<int> basinroots(sites.size(), -1);

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        int current = index;
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

        for (int trailindex : trail)
            basinroots[trailindex] = basinroots[current];
    }

    float maxsolvedelevation = 0.0f;

    for (const FastLEMSite& site : sites)
        maxsolvedelevation = max(maxsolvedelevation, site.solvedElevation);

    if (maxsolvedelevation <= 0.0f)
        return false;

    vector<bool> candidate(sites.size(), false);
    int candidatecount = 0;
    int boundarycandidatecount = 0;

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        FastLEMSite& site = sites[index];
        const float elevationnormalised = site.solvedElevation / maxsolvedelevation;
        const float coastnormalised = maxcoastdistance > 0 ? static_cast<float>(site.coastDistance) / static_cast<float>(maxcoastdistance) : 0.0f;
        const float score = elevationnormalised * tuning::terrain::fastlem::elevationScoreWeight + coastnormalised * tuning::terrain::fastlem::coastScoreWeight;

        bool localmaximum = true;
        int uniqueroots[8] = {};
        int uniquerootcount = 0;

        for (const FastLEMEdge& edge : site.edges)
        {
            const FastLEMSite& other = sites[edge.to];
            const float otherscore =
                (other.solvedElevation / maxsolvedelevation) * tuning::terrain::fastlem::elevationScoreWeight +
                (maxcoastdistance > 0 ? static_cast<float>(other.coastDistance) / static_cast<float>(maxcoastdistance) : 0.0f) * tuning::terrain::fastlem::coastScoreWeight;

            if (otherscore > score + 0.02f)
                localmaximum = false;

            const int root = basinroots[edge.to];
            bool seenroot = false;

            for (int rootindex = 0; rootindex < uniquerootcount; rootindex++)
            {
                if (uniqueroots[rootindex] == root)
                {
                    seenroot = true;
                    break;
                }
            }

            if (seenroot == false && uniquerootcount < 8)
                uniqueroots[uniquerootcount++] = root;
        }

        const bool boundarycandidate = uniquerootcount >= 2;

        if (site.coastDistance >= tuning::terrain::fastlem::minimumPeakCoastDistance &&
            elevationnormalised >= tuning::terrain::fastlem::minimumPeakElevationNormalised &&
            localmaximum)
        {
            candidate[index] = true;
        }

        if (site.coastDistance >= tuning::terrain::fastlem::minimumRidgeCoastDistance &&
            elevationnormalised >= tuning::terrain::fastlem::minimumRidgeElevationNormalised &&
            score >= tuning::terrain::fastlem::minimumRidgeScore &&
            (boundarycandidate || localmaximum))
        {
            candidate[index] = true;
        }

        if (candidate[index])
        {
            const float heightnormalised = clamp(elevationnormalised * 0.75f + coastnormalised * 0.25f, 0.0f, 1.0f);
            site.peakHeight = tuning::terrain::fastlem::minimumPeakHeight +
                static_cast<int>(powf(heightnormalised, tuning::terrain::fastlem::peakHeightExponent) *
                    static_cast<float>(tuning::terrain::fastlem::maximumPeakHeight - tuning::terrain::fastlem::minimumPeakHeight));

            candidatecount++;

            if (boundarycandidate)
                boundarycandidatecount++;
        }
    }

    const int minimumcandidatecount = min(tuning::terrain::fastlem::minimumCandidateSites, max(6, static_cast<int>(sites.size() / 40)));
    const int minimumboundarycount = min(tuning::terrain::fastlem::minimumBoundaryCandidates, max(3, static_cast<int>(sites.size() / 80)));

    if (candidatecount < minimumcandidatecount || boundarycandidatecount < minimumboundarycount)
        return false;

    vector<vector<int>> rawmountains(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, 0));

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        if (candidate[index] == false)
            continue;

        rawmountains[sites[index].x][sites[index].y] = max(rawmountains[sites[index].x][sites[index].y], sites[index].peakHeight);
    }

    for (int index = 0; index < static_cast<int>(sites.size()); index++)
    {
        if (candidate[index] == false)
            continue;

        for (const FastLEMEdge& edge : sites[index].edges)
        {
            if (edge.to <= index || candidate[edge.to] == false)
                continue;

            const bool samebasin = basinroots[index] == basinroots[edge.to];
            const bool closeband = abs(sites[index].coastDistance - sites[edge.to].coastDistance) <= 1;

            if (samebasin == false || closeband)
            {
                drawrawline(rawmountains, width, height, sites[index].x, sites[index].y, sites[edge.to].x, sites[edge.to].y, findedgepeak(sites[index], sites[edge.to]));
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
