#include <algorithm>
#include <cmath>
#include <deque>
#include <queue>
#include <string>
#include <vector>

#include "core/mesh/mesh_geometry.hpp"
#include "world/planet.hpp"
#include "simulations/society/social_generation_internal.hpp"

using namespace std;

namespace
{
using mesh::CellId;
using socialgen::scanslots;

struct SettlementCandidate
{
    CellId cell = 0;
    int suitability = 0;
    int choke = 0;
    int harbor = 0;
    int river = 0;
    double drawkey = 0.0;
};

// Steps to the nearest sea cell over the four edge neighbours (the former
// diamond rings), maxradius + 1 beyond maxradius.
vector<int> seadistances(const planet& world, int maxradius)
{
    const mesh::CellGrid& grid = world.grid();
    vector<int> distance(grid.cellcount(), maxradius + 1);
    deque<CellId> frontier;

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (world.sea(c))
        {
            distance[c] = 0;
            frontier.push_back(c);
        }
    }

    while (!frontier.empty())
    {
        const CellId c = frontier.front();
        frontier.pop_front();

        if (distance[c] >= maxradius)
            continue;

        for (int slot = 1; slot <= 7; slot += 2)
        {
            const int n = grid.neighbour(c, slot);

            if (n < 0 || distance[n] <= distance[c] + 1)
                continue;

            distance[n] = distance[c] + 1;
            frontier.push_back(static_cast<CellId>(n));
        }
    }

    return distance;
}

int countinflows(const planet& world, CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    int inflows = 0;

    for (int slot : scanslots)
    {
        const int n = grid.neighbour(c, slot);

        if (n < 0)
            continue;

        const int dir = world.riverdir(static_cast<CellId>(n));

        if (dir <= 0)
            continue;

        if (grid.neighbour(static_cast<CellId>(n), dir) == static_cast<int>(c))
            inflows++;
    }

    return inflows;
}

int computeslopepenalty(const planet& world, CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    float totaldelta = 0.0f;
    int samples = 0;
    const float baseelevation = world.mapfloat(c);

    for (int slot : scanslots)
    {
        const int n = grid.neighbour(c, slot);

        if (n < 0)
            continue;

        totaldelta += std::fabs(baseelevation - world.mapfloat(static_cast<CellId>(n)));
        samples++;
    }

    if (samples == 0)
        return 0;

    return std::clamp(static_cast<int>(totaldelta / static_cast<float>(samples) / 16.0f), 0, 100);
}

int computebiomesuitability(int biome)
{
    if (biome == biomeice)
        return 0;

    if (biome == biomepolardesert || biome == biomesubpolardesert || biome == biomeborealdesert || biome == biomecooltemperatedesert || biome == biomewarmtemperatedesert || biome == biomesubtropicaldesert || biome == biometropicaldesert)
        return 22;

    if (biome == biomecooltemperatesteppe || biome == biomewarmtemperatethornsteppe || biome == biomesubtropicalthornsteppe || biome == biometropicalthornsteppe)
        return 74;

    if (biome == biomecooltemperatemoistforest || biome == biomewarmtemperatemoistforest || biome == biomesubtropicalmoistforest || biome == biometropicalmoistforest)
        return 88;

    if (biome == biomecooltemperatewetforest || biome == biomewarmtemperatewetforest || biome == biomesubtropicalwetforest || biome == biometropicalwetforest || biome == biomecooltemperaterainforest || biome == biomewarmtemperaterainforest || biome == biomesubtropicalrainforest || biome == biometropicalrainforest)
        return 84;

    return 56;
}

bool hasphysicalreserve(const planet& world, CellId c)
{
    return world.floodplainfertility(c) > 0 || world.fisheryreserve(c) > 0 || world.metalorereserve(c) > 0 || world.placerreserve(c) > 0 || world.volcanicreserve(c) > 0;
}

bool hasphysicalreserveinputs(const planet& world)
{
    const mesh::CellGrid& grid = world.grid();

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (hasphysicalreserve(world, c))
            return true;
    }

    return false;
}

int computeresourceopportunity(const planet& world, CellId c, bool usephysicalreserves, int coastscore, int riverscore, int slopescore)
{
    if (usephysicalreserves)
    {
        const int reserve =
            std::max(world.metalorereserve(c), world.placerreserve(c)) +
            std::max(world.volcanicreserve(c), world.fisheryreserve(c)) +
            world.floodplainfertility(c);

        return std::clamp(reserve / 2, 0, 100);
    }

    int proxy = 0;
    proxy += coastscore / 4;
    proxy += riverscore / 3;
    proxy += slopescore / 6;
    proxy += world.mountainheightfloat(c) > 0 ? 24 : 0;
    proxy += std::abs(world.volcanofloat(c)) > 0 ? 30 : 0;
    return std::clamp(proxy, 0, 100);
}

// The cells whose centres lie within a radius in cell widths of a cell's centre.
vector<CellId> cellsaround(const mesh::CellGrid& grid, CellId centre, double radiuscells)
{
    return mesh::cellswithin(grid, grid.centre(centre), radiuscells * socialgen::cellwidth_m(grid));
}

string makesettlementname(int id, std::mt19937_64& rng)
{
    static const vector<string> starts = { "Al", "Bel", "Cor", "Dor", "Er", "Fal", "Gal", "Har", "Is", "Kor", "Lan", "Mor", "Nor", "Or", "Pel", "Quor", "Riv", "Sel", "Tor", "Val" };
    static const vector<string> ends = { "a", "ar", "en", "eth", "ford", "haven", "ia", "in", "mere", "on", "or", "os", "port", "stead", "ton", "vale" };
    std::uniform_int_distribution<int> startpick(0, static_cast<int>(starts.size()) - 1);
    std::uniform_int_distribution<int> endpick(0, static_cast<int>(ends.size()) - 1);
    return starts[startpick(rng)] + ends[endpick(rng)] + " " + to_string(id + 1);
}
}

namespace socialgen
{
void computesitesandinfrastructure(planet& world, const SocialProfile& profile, std::mt19937_64& rng)
{
    const mesh::CellGrid& grid = world.grid();
    const int cellcount = static_cast<int>(grid.cellcount());
    const int maxcoastradius = 12;
    const bool usephysicalreserves = hasphysicalreserveinputs(world);

    world.clearsocialstate();

    const vector<int> seadistance = seadistances(world, maxcoastradius);
    const vector<double> drawkeys = sitedrawkeys(grid, world.seed());
    double landarea_m2 = 0.0;
    vector<SettlementCandidate> candidates;
    candidates.reserve(static_cast<size_t>(cellcount) / 8);

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        world.setownersettlementid(c, -1);
        world.setownerpolityid(c, -1);
        world.setinfrastructure(c, 0);
        world.setroutetraffic(c, 0);

        if (world.sea(c))
        {
            world.setsettlementsuitability(c, 0);
            world.setagriculturalcapacity(c, 0);
            world.setharborscore(c, 0);
            world.setriveraccess(c, 0);
            continue;
        }

        landarea_m2 += grid.area_m2(c);

        const int coastdistance = seadistance[c];
        const int coastscore = std::clamp(100 - coastdistance * 8, 0, 100);

        int adjacentsea = 0;
        int adjacentland = 0;

        for (int slot : scanslots)
        {
            const int n = grid.neighbour(c, slot);

            if (n < 0)
                continue;

            if (world.sea(static_cast<CellId>(n)))
                adjacentsea++;
            else
                adjacentland++;
        }

        int harborscore = 0;

        if (coastdistance <= 2)
        {
            const int shelter = std::max(0, adjacentland - adjacentsea);
            harborscore = 42 + shelter * 7 + adjacentsea * 3;

            if (adjacentsea >= 2 && adjacentsea <= 5)
                harborscore += 12;
        }

        harborscore = std::clamp(harborscore, 0, 100);
        world.setharborscore(c, harborscore);

        const int riverflow = world.riveraveflow(c);
        const int riverscore = std::clamp(static_cast<int>(std::log1p(static_cast<float>(riverflow)) * 14.0f), 0, 100);
        const int inflows = countinflows(world, c);
        const int confluencescore = inflows >= 2 ? std::clamp(18 + inflows * 18, 0, 100) : 0;
        const int slopepenalty = computeslopepenalty(world, c);
        const int slopescore = std::clamp(100 - slopepenalty, 0, 100);

        const float averageclimate = (world.jantempfloat(c) + world.jultempfloat(c)) / 2.0f;
        const float averagerain = (world.janrainfloat(c) + world.julrainfloat(c)) / 2.0f;
        const int temperaturescore = std::clamp(100 - static_cast<int>(std::abs(averageclimate - 14.0f) * 4.0f), 0, 100);
        const int rainscore = std::clamp(100 - static_cast<int>(std::abs(averagerain - 80.0f) * 0.8f), 0, 100);
        const int biomescore = computebiomesuitability(world.biome(c));
        const int climatescore = std::clamp((temperaturescore * 2 + rainscore + biomescore) / 4, 0, 100);

        int gradientpenalty = 0;
        const int downstream = grid.neighbour(c, world.riverdir(c));

        if (downstream >= 0)
            gradientpenalty = static_cast<int>(std::max(0.0f, world.mapfloat(c) - world.mapfloat(static_cast<CellId>(downstream))) / 24.0f);

        const int riveraccess = std::clamp(riverscore - gradientpenalty + (riverflow > 300 ? 12 : 0), 0, 100);
        world.setriveraccess(c, riveraccess);

        int fertility = 0;

        if (usephysicalreserves)
        {
            fertility = world.floodplainfertility(c);
        }
        else
        {
            fertility = std::clamp((riverscore * 4 + rainscore * 3 + slopescore * 2 + coastscore) / 10, 0, 100);
        }

        const int agriculture = std::clamp((fertility * 45 + climatescore * 25 + riveraccess * 20 + slopescore * 10) / 100, 0, 100);
        world.setagriculturalcapacity(c, agriculture);

        const int routeaccessproxy = std::clamp((riverscore + coastscore + harborscore) / 3, 0, 100);
        const int resourceopportunity = computeresourceopportunity(world, c, usephysicalreserves, coastscore, riverscore, slopescore);

        const float weightedsum =
            static_cast<float>(coastscore) * profile.coastWeight +
            static_cast<float>(harborscore) * profile.harborWeight +
            static_cast<float>(riverscore) * profile.riverWeight +
            static_cast<float>(confluencescore) * profile.confluenceWeight +
            static_cast<float>(climatescore) * profile.climateWeight +
            static_cast<float>(agriculture) * profile.agricultureWeight +
            static_cast<float>(routeaccessproxy) * profile.routeWeight +
            static_cast<float>(resourceopportunity) * profile.resourceWeight;

        const float weighttotal =
            profile.coastWeight +
            profile.harborWeight +
            profile.riverWeight +
            profile.confluenceWeight +
            profile.climateWeight +
            profile.agricultureWeight +
            profile.routeWeight +
            profile.resourceWeight;

        int suitability = 0;

        if (weighttotal > 0.0f)
            suitability = static_cast<int>(std::round(weightedsum / weighttotal));

        suitability -= static_cast<int>(std::round(profile.terrainPenalty * static_cast<float>(slopepenalty) * 0.30f));
        suitability = std::clamp(suitability, 0, 100);

        world.setsettlementsuitability(c, suitability);

        if (suitability >= 35)
        {
            SettlementCandidate candidate;
            candidate.cell = c;
            candidate.suitability = suitability;
            candidate.harbor = harborscore;
            candidate.river = riveraccess;
            candidate.drawkey = drawkeys[c];
            candidate.choke = std::clamp((harborscore + riveraccess + confluencescore + routeaccessproxy) / 4, 0, 100);
            candidates.push_back(candidate);
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const SettlementCandidate& left, const SettlementCandidate& right)
    {
        if (left.suitability != right.suitability)
            return left.suitability > right.suitability;

        if (left.choke != right.choke)
            return left.choke > right.choke;

        if (left.harbor != right.harbor)
            return left.harbor > right.harbor;

        if (left.river != right.river)
            return left.river > right.river;

        if (left.drawkey != right.drawkey)
            return left.drawkey > right.drawkey;

        return left.cell < right.cell;
    });

    vector<Settlement>& settlements = world.settlements();
    // The land area in export cells (the former count of land cells).
    const int landcells = static_cast<int>(landarea_m2 / exportcellarea_m2(grid));
    const int targetsettlements = std::clamp(landcells / 12000, 8, 180);
    const float worldscalefactor = static_cast<float>(exportcolumns(grid)) / 1024.0f;
    const int spacingradius = std::max(6, static_cast<int>(std::round(static_cast<float>(profile.spacingRadius) * worldscalefactor)));

    for (const SettlementCandidate& candidate : candidates)
    {
        if (static_cast<int>(settlements.size()) >= targetsettlements)
            break;

        bool blocked = false;

        for (const Settlement& existing : settlements)
        {
            const float distance = celldistance(grid, candidate.cell, existing.cell);

            if (distance >= static_cast<float>(spacingradius))
                continue;

            const int existingscore = world.settlementsuitability(existing.cell);
            const bool privilegedchokepoint = candidate.harbor >= 85 || candidate.river >= 85;

            if (!privilegedchokepoint || distance < static_cast<float>(spacingradius) * 0.5f)
            {
                if (candidate.suitability + candidate.choke / 6 < existingscore + 8)
                {
                    blocked = true;
                    break;
                }
            }
        }

        if (blocked)
            continue;

        Settlement settlement;
        settlement.id = static_cast<int>(settlements.size());
        placesettlement(settlement, grid, candidate.cell);
        settlement.name = makesettlementname(settlement.id, rng);
        settlement.harbor = static_cast<float>(candidate.harbor);
        settlement.riverAccess = static_cast<float>(candidate.river);
        settlement.infrastructure = static_cast<float>(std::clamp(35 + candidate.suitability / 2, 20, 100));
        settlements.push_back(std::move(settlement));
    }

    if (settlements.empty() && !candidates.empty())
    {
        Settlement settlement;
        settlement.id = 0;
        placesettlement(settlement, grid, candidates[0].cell);
        settlement.name = makesettlementname(0, rng);
        settlement.harbor = static_cast<float>(candidates[0].harbor);
        settlement.riverAccess = static_cast<float>(candidates[0].river);
        settlement.infrastructure = static_cast<float>(std::clamp(35 + candidates[0].suitability / 2, 20, 100));
        settlements.push_back(std::move(settlement));
    }

    for (Settlement& settlement : settlements)
    {
        const int radius = std::max(5, spacingradius / 2);
        int agriculturetotal = 0;
        int routeaccesstotal = 0;
        int samples = 0;

        for (const CellId n : cellsaround(grid, settlement.cell, static_cast<double>(radius)))
        {
            if (world.sea(n))
                continue;

            agriculturetotal += world.agriculturalcapacity(n);
            routeaccesstotal += std::max({ world.riveraccess(n), world.harborscore(n), world.settlementsuitability(n) });
            samples++;
        }

        const int suitability = world.settlementsuitability(settlement.cell);
        const float averageagriculture = samples > 0 ? static_cast<float>(agriculturetotal) / static_cast<float>(samples) : static_cast<float>(world.agriculturalcapacity(settlement.cell));
        const float averagerouteaccess = samples > 0 ? static_cast<float>(routeaccesstotal) / static_cast<float>(samples) : static_cast<float>(std::max(world.riveraccess(settlement.cell), world.harborscore(settlement.cell)));
        float carryingcapacity = 800.0f + averageagriculture * 120.0f + averagerouteaccess * 65.0f + static_cast<float>(suitability) * 70.0f;

        float crowdingpenalty = 0.0f;

        for (const Settlement& other : settlements)
        {
            if (other.id == settlement.id)
                continue;

            const float distance = celldistance(grid, settlement.cell, other.cell);

            if (distance >= static_cast<float>(spacingradius * 2))
                continue;

            const int otherscore = world.settlementsuitability(other.cell);

            if (otherscore <= suitability)
                continue;

            const float pressure = (static_cast<float>(spacingradius * 2) - distance) / static_cast<float>(spacingradius * 2);
            crowdingpenalty += pressure * 0.18f;
        }

        crowdingpenalty = std::clamp(crowdingpenalty, 0.0f, 0.45f);
        carryingcapacity *= (1.0f - crowdingpenalty);
        carryingcapacity *= 1.0f + std::max(settlement.harbor, settlement.riverAccess) * 0.003f;
        carryingcapacity = std::max(250.0f, carryingcapacity);

        const int totalpopulation = std::max(120, static_cast<int>(std::round(carryingcapacity * 0.72f)));
        const float urbanratio = std::clamp(profile.urbanizationBias + averagerouteaccess * 0.0025f + std::max(settlement.harbor, settlement.riverAccess) * 0.0015f, 0.08f, 0.58f);
        const int urbanpopulation = std::clamp(static_cast<int>(std::round(static_cast<float>(totalpopulation) * urbanratio)), 10, totalpopulation);

        settlement.carryingCapacity = carryingcapacity;
        settlement.urbanPopulation = urbanpopulation;
        settlement.ruralPopulation = totalpopulation - urbanpopulation;
        settlement.marketStrength = std::clamp(averagerouteaccess * 0.55f + static_cast<float>(suitability) * 0.45f, 0.0f, 100.0f);
        settlement.polityId = -1;
    }

    struct Node
    {
        int strength = 0;
        CellId cell = 0;
        int settlementId = -1;

        bool operator<(const Node& other) const
        {
            return strength < other.strength;
        }
    };

    vector<int> beststrength(static_cast<size_t>(cellcount), -1);
    vector<int> bestowner(static_cast<size_t>(cellcount), -1);
    priority_queue<Node> frontier;

    for (const Settlement& settlement : settlements)
    {
        const int startstrength = std::clamp(static_cast<int>(55.0f + settlement.infrastructure * 0.65f), 35, 100);
        beststrength[settlement.cell] = startstrength;
        bestowner[settlement.cell] = settlement.id;
        frontier.push({ startstrength, settlement.cell, settlement.id });
    }

    while (!frontier.empty())
    {
        const Node node = frontier.top();
        frontier.pop();

        if (node.strength < beststrength[node.cell] || node.settlementId != bestowner[node.cell])
            continue;

        for (int slot : scanslots)
        {
            const int neighbour = grid.neighbour(node.cell, slot);

            if (neighbour < 0 || world.sea(static_cast<CellId>(neighbour)))
                continue;

            const CellId n = static_cast<CellId>(neighbour);
            int movecost = 4;

            if (!orthogonalslot(slot))
                movecost++;

            const float elevationdelta = std::fabs(world.mapfloat(n) - world.mapfloat(node.cell));
            movecost += std::min(8, static_cast<int>(elevationdelta / 180.0f));

            const int special = world.special(n);
            if (special >= 130 && special <= 132)
                movecost += 3;

            const bool currentriver = world.riveraveflow(node.cell) > 220;
            const bool nextriver = world.riveraveflow(n) > 220;
            if (currentriver != nextriver)
                movecost += 2;

            if (world.riveraccess(n) > 40)
                movecost = std::max(1, movecost - 1);

            if (world.harborscore(n) > 60)
                movecost = std::max(1, movecost - 1);

            const int newstrength = node.strength - movecost;

            if (newstrength <= 0)
                continue;

            if (newstrength > beststrength[n] || (newstrength == beststrength[n] && node.settlementId < bestowner[n]))
            {
                beststrength[n] = newstrength;
                bestowner[n] = node.settlementId;
                frontier.push({ newstrength, n, node.settlementId });
            }
        }
    }

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (world.sea(c))
        {
            world.setinfrastructure(c, 0);
            world.setownersettlementid(c, -1);
            continue;
        }

        world.setinfrastructure(c, std::max(0, beststrength[c]));
        world.setownersettlementid(c, bestowner[c]);
    }

    for (Settlement& settlement : settlements)
    {
        int total = 0;
        int samples = 0;

        // The disc of the former 9 x 9 square's area.
        const double squarehalfwidth = 4.5;
        const double pi = 3.14159265358979323846;

        for (const CellId n : cellsaround(grid, settlement.cell, squarehalfwidth * 2.0 / std::sqrt(pi)))
        {
            if (world.sea(n))
                continue;

            total += world.infrastructure(n);
            samples++;
        }

        if (samples > 0)
            settlement.infrastructure = static_cast<float>(total) / static_cast<float>(samples);
    }
}
}
