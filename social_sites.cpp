#include <algorithm>
#include <cmath>
#include <queue>
#include <string>
#include <vector>

#include "functions.hpp"
#include "planet.hpp"
#include "social_generation_internal.hpp"

using namespace std;

namespace
{
struct SettlementCandidate
{
    int x = 0;
    int y = 0;
    int suitability = 0;
    int choke = 0;
    int harbor = 0;
    int river = 0;
};

int wrapxlocal(int x, int width)
{
    while (x < 0)
        x += width + 1;

    while (x > width)
        x -= width + 1;

    return x;
}

int wrappedxdistance(int x1, int x2, int width)
{
    const int worldspan = width + 1;
    const int direct = std::abs(x1 - x2);
    return std::min(direct, worldspan - direct);
}

float wrappeddistance(int x1, int y1, int x2, int y2, int width)
{
    const int dx = wrappedxdistance(x1, x2, width);
    const int dy = std::abs(y1 - y2);
    return std::sqrt(static_cast<float>(dx * dx + dy * dy));
}

int distancetosea(planet& world, int x, int y, int maxradius)
{
    if (world.sea(x, y))
        return 0;

    for (int radius = 1; radius <= maxradius; radius++)
    {
        for (int dx = -radius; dx <= radius; dx++)
        {
            const int dy = radius - std::abs(dx);
            const int nx1 = wrapxlocal(x + dx, world.width());
            const int ny1 = y - dy;
            const int ny2 = y + dy;

            if (ny1 >= 0 && ny1 <= world.height() && world.sea(nx1, ny1))
                return radius;

            if (dy != 0 && ny2 >= 0 && ny2 <= world.height() && world.sea(nx1, ny2))
                return radius;
        }
    }

    return maxradius + 1;
}

int countinflows(planet& world, int x, int y)
{
    int inflows = 0;

    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            if (dx == 0 && dy == 0)
                continue;

            const int nx = wrapxlocal(x + dx, world.width());
            const int ny = y + dy;

            if (ny < 0 || ny > world.height())
                continue;

            const int dir = world.riverdir(nx, ny);

            if (dir <= 0)
                continue;

            const twointegers destination = getdestination(nx, ny, dir);
            const int destx = wrapxlocal(destination.x, world.width());
            const int desty = std::clamp(destination.y, 0, world.height());

            if (destx == x && desty == y)
                inflows++;
        }
    }

    return inflows;
}

int computeslopepenalty(planet& world, int x, int y)
{
    int totaldelta = 0;
    int samples = 0;
    const int baseelevation = world.map(x, y);

    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            if (dx == 0 && dy == 0)
                continue;

            const int nx = wrapxlocal(x + dx, world.width());
            const int ny = y + dy;

            if (ny < 0 || ny > world.height())
                continue;

            totaldelta += std::abs(baseelevation - world.map(nx, ny));
            samples++;
        }
    }

    if (samples == 0)
        return 0;

    return std::clamp((totaldelta / samples) / 16, 0, 100);
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

bool hasphysicalreserveinputs(planet& world)
{
    for (int x = 0; x <= world.width(); x += 16)
    {
        for (int y = 0; y <= world.height(); y += 16)
        {
            if (world.floodplainfertility(x, y) > 0 || world.fisheryreserve(x, y) > 0 || world.metalorereserve(x, y) > 0 || world.placerreserve(x, y) > 0 || world.volcanicreserve(x, y) > 0)
                return true;
        }
    }

    return false;
}

int computeresourceopportunity(planet& world, int x, int y, bool usephysicalreserves, int coastscore, int riverscore, int slopescore)
{
    if (usephysicalreserves)
    {
        const int reserve =
            std::max(world.metalorereserve(x, y), world.placerreserve(x, y)) +
            std::max(world.volcanicreserve(x, y), world.fisheryreserve(x, y)) +
            world.floodplainfertility(x, y);

        return std::clamp(reserve / 2, 0, 100);
    }

    int proxy = 0;
    proxy += coastscore / 4;
    proxy += riverscore / 3;
    proxy += slopescore / 6;
    proxy += world.mountainheight(x, y) > 0 ? 24 : 0;
    proxy += std::abs(world.volcano(x, y)) > 0 ? 30 : 0;
    return std::clamp(proxy, 0, 100);
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
    const int width = world.width();
    const int height = world.height();
    const int maxcoastradius = 12;
    const bool usephysicalreserves = hasphysicalreserveinputs(world);

    world.clearsocialstate();

    int landcells = 0;
    vector<SettlementCandidate> candidates;
    candidates.reserve(((width + 1) * (height + 1)) / 8);

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            world.setownersettlementid(x, y, -1);
            world.setownerpolityid(x, y, -1);
            world.setinfrastructure(x, y, 0);
            world.setroutetraffic(x, y, 0);

            if (world.sea(x, y))
            {
                world.setsettlementsuitability(x, y, 0);
                world.setagriculturalcapacity(x, y, 0);
                world.setharborscore(x, y, 0);
                world.setriveraccess(x, y, 0);
                continue;
            }

            landcells++;

            const int coastdistance = distancetosea(world, x, y, maxcoastradius);
            const int coastscore = std::clamp(100 - coastdistance * 8, 0, 100);

            int adjacentsea = 0;
            int adjacentland = 0;

            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    if (dx == 0 && dy == 0)
                        continue;

                    const int nx = wrapxlocal(x + dx, width);
                    const int ny = y + dy;

                    if (ny < 0 || ny > height)
                        continue;

                    if (world.sea(nx, ny))
                        adjacentsea++;
                    else
                        adjacentland++;
                }
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
            world.setharborscore(x, y, harborscore);

            const int riverflow = world.riveraveflow(x, y);
            const int riverscore = std::clamp(static_cast<int>(std::log1p(static_cast<float>(riverflow)) * 14.0f), 0, 100);
            const int inflows = countinflows(world, x, y);
            const int confluencescore = inflows >= 2 ? std::clamp(18 + inflows * 18, 0, 100) : 0;
            const int slopepenalty = computeslopepenalty(world, x, y);
            const int slopescore = std::clamp(100 - slopepenalty, 0, 100);

            const float averageclimate = static_cast<float>(world.jantemp(x, y) + world.jultemp(x, y)) / 2.0f;
            const float averagerain = static_cast<float>(world.janrain(x, y) + world.julrain(x, y)) / 2.0f;
            const int temperaturescore = std::clamp(100 - static_cast<int>(std::abs(averageclimate - 14.0f) * 4.0f), 0, 100);
            const int rainscore = std::clamp(100 - static_cast<int>(std::abs(averagerain - 80.0f) * 0.8f), 0, 100);
            const int biomescore = computebiomesuitability(world.biome(x, y));
            const int climatescore = std::clamp((temperaturescore * 2 + rainscore + biomescore) / 4, 0, 100);

            const int riverdir = world.riverdir(x, y);
            int gradientpenalty = 0;

            if (riverdir > 0)
            {
                const twointegers downstream = getdestination(x, y, riverdir);
                const int downstreamx = wrapxlocal(downstream.x, width);
                const int downstreamy = std::clamp(downstream.y, 0, height);
                gradientpenalty = std::max(0, world.map(x, y) - world.map(downstreamx, downstreamy)) / 24;
            }

            const int riveraccess = std::clamp(riverscore - gradientpenalty + (riverflow > 300 ? 12 : 0), 0, 100);
            world.setriveraccess(x, y, riveraccess);

            int fertility = 0;

            if (usephysicalreserves)
            {
                fertility = world.floodplainfertility(x, y);
            }
            else
            {
                fertility = std::clamp((riverscore * 4 + rainscore * 3 + slopescore * 2 + coastscore) / 10, 0, 100);
            }

            const int agriculture = std::clamp((fertility * 45 + climatescore * 25 + riveraccess * 20 + slopescore * 10) / 100, 0, 100);
            world.setagriculturalcapacity(x, y, agriculture);

            const int routeaccessproxy = std::clamp((riverscore + coastscore + harborscore) / 3, 0, 100);
            const int resourceopportunity = computeresourceopportunity(world, x, y, usephysicalreserves, coastscore, riverscore, slopescore);

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

            world.setsettlementsuitability(x, y, suitability);

            if (suitability >= 35)
            {
                SettlementCandidate candidate;
                candidate.x = x;
                candidate.y = y;
                candidate.suitability = suitability;
                candidate.harbor = harborscore;
                candidate.river = riveraccess;
                candidate.choke = std::clamp((harborscore + riveraccess + confluencescore + routeaccessproxy) / 4, 0, 100);
                candidates.push_back(candidate);
            }
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

        return left.river > right.river;
    });

    vector<Settlement>& settlements = world.settlements();
    const int targetsettlements = std::clamp(landcells / 12000, 8, 180);
    const float worldscalefactor = static_cast<float>(width + 1) / 1025.0f;
    const int spacingradius = std::max(6, static_cast<int>(std::round(static_cast<float>(profile.spacingRadius) * worldscalefactor)));

    for (const SettlementCandidate& candidate : candidates)
    {
        if (static_cast<int>(settlements.size()) >= targetsettlements)
            break;

        bool blocked = false;

        for (const Settlement& existing : settlements)
        {
            const float distance = wrappeddistance(candidate.x, candidate.y, existing.x, existing.y, width);

            if (distance >= static_cast<float>(spacingradius))
                continue;

            const int existingscore = world.settlementsuitability(existing.x, existing.y);
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
        settlement.x = candidate.x;
        settlement.y = candidate.y;
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
        settlement.x = candidates[0].x;
        settlement.y = candidates[0].y;
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

        for (int dx = -radius; dx <= radius; dx++)
        {
            for (int dy = -radius; dy <= radius; dy++)
            {
                if (dx * dx + dy * dy > radius * radius)
                    continue;

                const int nx = wrapxlocal(settlement.x + dx, width);
                const int ny = settlement.y + dy;

                if (ny < 0 || ny > height || world.sea(nx, ny))
                    continue;

                agriculturetotal += world.agriculturalcapacity(nx, ny);
                routeaccesstotal += std::max({ world.riveraccess(nx, ny), world.harborscore(nx, ny), world.settlementsuitability(nx, ny) });
                samples++;
            }
        }

        const int suitability = world.settlementsuitability(settlement.x, settlement.y);
        const float averageagriculture = samples > 0 ? static_cast<float>(agriculturetotal) / static_cast<float>(samples) : static_cast<float>(world.agriculturalcapacity(settlement.x, settlement.y));
        const float averagerouteaccess = samples > 0 ? static_cast<float>(routeaccesstotal) / static_cast<float>(samples) : static_cast<float>(std::max(world.riveraccess(settlement.x, settlement.y), world.harborscore(settlement.x, settlement.y)));
        float carryingcapacity = 800.0f + averageagriculture * 120.0f + averagerouteaccess * 65.0f + static_cast<float>(suitability) * 70.0f;

        float crowdingpenalty = 0.0f;

        for (const Settlement& other : settlements)
        {
            if (other.id == settlement.id)
                continue;

            const float distance = wrappeddistance(settlement.x, settlement.y, other.x, other.y, width);

            if (distance >= static_cast<float>(spacingradius * 2))
                continue;

            const int otherscore = world.settlementsuitability(other.x, other.y);

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
        int x = 0;
        int y = 0;
        int settlementId = -1;

        bool operator<(const Node& other) const
        {
            return strength < other.strength;
        }
    };

    const int cellcount = (width + 1) * (height + 1);
    vector<int> beststrength(static_cast<size_t>(cellcount), -1);
    vector<int> bestowner(static_cast<size_t>(cellcount), -1);
    priority_queue<Node> frontier;

    const auto indexof = [height](int x, int y) -> int
    {
        return x * (height + 1) + y;
    };

    for (const Settlement& settlement : settlements)
    {
        const int startstrength = std::clamp(static_cast<int>(55.0f + settlement.infrastructure * 0.65f), 35, 100);
        const int index = indexof(settlement.x, settlement.y);
        beststrength[index] = startstrength;
        bestowner[index] = settlement.id;
        frontier.push({ startstrength, settlement.x, settlement.y, settlement.id });
    }

    while (!frontier.empty())
    {
        const Node node = frontier.top();
        frontier.pop();

        const int currentindex = indexof(node.x, node.y);

        if (node.strength < beststrength[currentindex] || node.settlementId != bestowner[currentindex])
            continue;

        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                if (dx == 0 && dy == 0)
                    continue;

                const int nx = wrapxlocal(node.x + dx, width);
                const int ny = node.y + dy;

                if (ny < 0 || ny > height || world.sea(nx, ny))
                    continue;

                int movecost = 4;

                if (dx != 0 && dy != 0)
                    movecost++;

                const int elevationdelta = std::abs(world.map(nx, ny) - world.map(node.x, node.y));
                movecost += std::min(8, elevationdelta / 180);

                const int special = world.special(nx, ny);
                if (special >= 130 && special <= 132)
                    movecost += 3;

                const bool currentriver = world.riveraveflow(node.x, node.y) > 220;
                const bool nextriver = world.riveraveflow(nx, ny) > 220;
                if (currentriver != nextriver)
                    movecost += 2;

                if (world.riveraccess(nx, ny) > 40)
                    movecost = std::max(1, movecost - 1);

                if (world.harborscore(nx, ny) > 60)
                    movecost = std::max(1, movecost - 1);

                const int newstrength = node.strength - movecost;

                if (newstrength <= 0)
                    continue;

                const int neighbourindex = indexof(nx, ny);

                if (newstrength > beststrength[neighbourindex] || (newstrength == beststrength[neighbourindex] && node.settlementId < bestowner[neighbourindex]))
                {
                    beststrength[neighbourindex] = newstrength;
                    bestowner[neighbourindex] = node.settlementId;
                    frontier.push({ newstrength, nx, ny, node.settlementId });
                }
            }
        }
    }

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            if (world.sea(x, y))
            {
                world.setinfrastructure(x, y, 0);
                world.setownersettlementid(x, y, -1);
                continue;
            }

            const int index = indexof(x, y);
            world.setinfrastructure(x, y, std::max(0, beststrength[index]));
            world.setownersettlementid(x, y, bestowner[index]);
        }
    }

    for (Settlement& settlement : settlements)
    {
        int total = 0;
        int samples = 0;

        for (int dx = -4; dx <= 4; dx++)
        {
            for (int dy = -4; dy <= 4; dy++)
            {
                const int nx = wrapxlocal(settlement.x + dx, width);
                const int ny = settlement.y + dy;

                if (ny < 0 || ny > height || world.sea(nx, ny))
                    continue;

                total += world.infrastructure(nx, ny);
                samples++;
            }
        }

        if (samples > 0)
            settlement.infrastructure = static_cast<float>(total) / static_cast<float>(samples);
    }
}
}
