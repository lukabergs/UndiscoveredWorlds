#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "planet.hpp"
#include "social_generation_internal.hpp"

using namespace std;

namespace
{
int wrapxlocal(int x, int width)
{
    while (x < 0)
        x += width + 1;

    while (x > width)
        x -= width + 1;

    return x;
}

float settlementdistance(const Settlement& left, const Settlement& right, int width)
{
    const int worldspan = width + 1;
    const int directdx = std::abs(left.x - right.x);
    const int dx = std::min(directdx, worldspan - directdx);
    const int dy = std::abs(left.y - right.y);
    return std::sqrt(static_cast<float>(dx * dx + dy * dy));
}

void recalcpolitystats(vector<Settlement>& settlements, vector<Polity>& polities, int worldwidth)
{
    vector<int> settlementcounts(polities.size(), 0);
    vector<float> distanceaccum(polities.size(), 0.0f);

    for (Polity& polity : polities)
    {
        polity.population = 0;
        polity.infrastructure = 0.0f;
        polity.wealth = 0.0f;
        polity.militaryPressure = 0.0f;
    }

    for (const Settlement& settlement : settlements)
    {
        if (settlement.polityId < 0 || settlement.polityId >= static_cast<int>(polities.size()))
            continue;

        Polity& polity = polities[settlement.polityId];
        polity.population += settlement.urbanPopulation + settlement.ruralPopulation;
        polity.infrastructure += settlement.infrastructure;
        polity.wealth += settlement.marketStrength;
        settlementcounts[settlement.polityId]++;

        const int capitalsettlementid = polity.capitalSettlementId;

        if (capitalsettlementid >= 0 && capitalsettlementid < static_cast<int>(settlements.size()))
            distanceaccum[settlement.polityId] += settlementdistance(settlement, settlements[capitalsettlementid], worldwidth);
    }

    for (Polity& polity : polities)
    {
        if (polity.id < 0 || polity.id >= static_cast<int>(settlementcounts.size()))
            continue;

        const int count = std::max(1, settlementcounts[polity.id]);
        polity.infrastructure /= static_cast<float>(count);
        polity.wealth /= static_cast<float>(count);
        const float averagedistance = distanceaccum[polity.id] / static_cast<float>(count);
        const float overextendpenalty = count > 8 ? static_cast<float>(count - 8) * 1.8f : 0.0f;
        const float frontierpenalty = averagedistance / std::max(1.0f, static_cast<float>(worldwidth)) * 110.0f;
        const float cohesion = 82.0f + polity.infrastructure * 0.18f - frontierpenalty - overextendpenalty;
        polity.cohesion = std::clamp(cohesion, 12.0f, 95.0f);
    }
}
}

namespace socialgen
{
void buildpolities(planet& world, const SocialGenerationOptions&, std::mt19937_64& rng)
{
    vector<Settlement>& settlements = world.settlements();
    vector<Polity>& polities = world.polities();
    vector<HistoryEvent>& history = world.historyevents();
    polities.clear();
    history.clear();

    if (settlements.empty())
        return;

    vector<vector<pair<int, float>>> adjacency(settlements.size());

    for (const RouteEdge& route : world.routeedges())
    {
        if (route.fromSettlementId < 0 || route.toSettlementId < 0 || route.fromSettlementId >= static_cast<int>(settlements.size()) || route.toSettlementId >= static_cast<int>(settlements.size()))
            continue;

        adjacency[route.fromSettlementId].push_back({ route.toSettlementId, route.traffic });
        adjacency[route.toSettlementId].push_back({ route.fromSettlementId, route.traffic });
    }

    vector<int> component(settlements.size(), -1);
    int componentcount = 0;

    for (int i = 0; i < static_cast<int>(settlements.size()); i++)
    {
        if (component[i] != -1)
            continue;

        queue<int> frontier;
        frontier.push(i);
        component[i] = componentcount;

        while (!frontier.empty())
        {
            const int current = frontier.front();
            frontier.pop();

            for (const pair<int, float>& link : adjacency[current])
            {
                if (component[link.first] == -1)
                {
                    component[link.first] = componentcount;
                    frontier.push(link.first);
                }
            }
        }

        componentcount++;
    }

    vector<pair<float, int>> capitalscores;
    capitalscores.reserve(settlements.size());

    for (const Settlement& settlement : settlements)
    {
        const float score = static_cast<float>(settlement.urbanPopulation) * std::max(1.0f, settlement.infrastructure);
        capitalscores.push_back({ score, settlement.id });
    }

    std::sort(capitalscores.begin(), capitalscores.end(), [](const pair<float, int>& left, const pair<float, int>& right)
    {
        return left.first > right.first;
    });

    const int desiredpolities = std::clamp(static_cast<int>(settlements.size()) / 6, 1, 24);
    const int spacingsquared = std::max(64, (world.width() / 18) * (world.width() / 18));
    vector<int> selectedcapitals;

    for (const pair<float, int>& scoredcapital : capitalscores)
    {
        if (static_cast<int>(selectedcapitals.size()) >= desiredpolities)
            break;

        const Settlement& candidate = settlements[scoredcapital.second];
        bool blocked = false;

        for (int selected : selectedcapitals)
        {
            const Settlement& existing = settlements[selected];
            const int dx = std::min(std::abs(existing.x - candidate.x), world.width() + 1 - std::abs(existing.x - candidate.x));
            const int dy = std::abs(existing.y - candidate.y);

            if (dx * dx + dy * dy < spacingsquared && component[existing.id] == component[candidate.id])
            {
                blocked = true;
                break;
            }
        }

        if (!blocked)
            selectedcapitals.push_back(scoredcapital.second);
    }

    if (selectedcapitals.empty())
        selectedcapitals.push_back(capitalscores.front().second);

    std::uniform_int_distribution<int> politynamepick(100, 999);

    for (int capitalid : selectedcapitals)
    {
        Polity polity;
        polity.id = static_cast<int>(polities.size());
        polity.capitalSettlementId = capitalid;
        polity.name = "Polity " + to_string(politynamepick(rng));
        polities.push_back(std::move(polity));
        settlements[capitalid].polityId = polities.back().id;

        HistoryEvent event;
        event.year = 0;
        event.type = "polity_founding";
        event.primarySettlementId = capitalid;
        event.primaryPolityId = polities.back().id;
        event.summary = polities.back().name + " formed around " + settlements[capitalid].name;
        history.push_back(std::move(event));
    }

    for (Settlement& settlement : settlements)
    {
        if (settlement.polityId >= 0)
            continue;

        float bestscore = std::numeric_limits<float>::max();
        int bestpolity = 0;

        for (const Polity& polity : polities)
        {
            const Settlement& capital = settlements[polity.capitalSettlementId];
            const float distance = settlementdistance(settlement, capital, world.width());
            float score = distance;

            if (component[settlement.id] == component[capital.id])
                score -= 35.0f;

            if (settlement.harbor > 65.0f && capital.harbor > 65.0f)
                score -= 8.0f;

            if (score < bestscore)
            {
                bestscore = score;
                bestpolity = polity.id;
            }
        }

        settlement.polityId = bestpolity;
    }

    recalcpolitystats(settlements, polities, world.width());

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
        {
            const int owner = world.ownersettlementid(x, y);

            if (owner >= 0 && owner < static_cast<int>(settlements.size()))
                world.setownerpolityid(x, y, settlements[owner].polityId);
            else
                world.setownerpolityid(x, y, -1);

            world.settest(x, y, 0);
        }
    }

    vector<float> borderpressure(polities.size(), 0.0f);

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 1; y < world.height(); y++)
        {
            if (world.sea(x, y))
                continue;

            const int owner = world.ownerpolityid(x, y);
            int friction = 0;

            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    if (dx == 0 && dy == 0)
                        continue;

                    const int nx = wrapxlocal(x + dx, world.width());
                    const int ny = std::clamp(y + dy, 0, world.height());
                    const int neighbourowner = world.ownerpolityid(nx, ny);

                    if (owner >= 0 && neighbourowner >= 0 && owner != neighbourowner)
                    {
                        friction++;
                        borderpressure[owner] += 0.25f;
                    }
                }
            }

            world.settest(x, y, std::clamp(friction * 18, 0, 100));
        }
    }

    for (Polity& polity : polities)
    {
        if (polity.id >= 0 && polity.id < static_cast<int>(borderpressure.size()))
        {
            polity.militaryPressure = borderpressure[polity.id];
            polity.cohesion = std::clamp(polity.cohesion - polity.militaryPressure * 0.4f, 8.0f, 95.0f);
        }
    }

    for (const Settlement& settlement : settlements)
    {
        HistoryEvent event;
        event.year = 0;
        event.type = "founding";
        event.primarySettlementId = settlement.id;
        event.primaryPolityId = settlement.polityId;
        event.summary = settlement.name + " was founded";
        history.push_back(std::move(event));
    }
}

void runhistoricalmode(planet& world, const SocialGenerationOptions& options, std::mt19937_64& rng)
{
    if (options.mode != SocialGenerationOptions::Mode::historical)
        return;

    vector<Settlement>& settlements = world.settlements();
    vector<Polity>& polities = world.polities();
    vector<HistoryEvent>& history = world.historyevents();

    if (settlements.empty() || polities.empty())
        return;

    if (options.usePrehistory)
    {
        HistoryEvent prehistory;
        prehistory.year = -std::max(10, options.historyYears);
        prehistory.type = "prehistory";
        prehistory.summary = "Prehistory migration pressures shaped early settlement distribution";
        history.push_back(std::move(prehistory));
    }

    const int years = std::max(20, options.historyYears);
    const int tickyears = 10;
    const int ticks = std::max(1, years / tickyears);
    std::uniform_real_distribution<float> randomdrift(-2.5f, 2.5f);
    std::uniform_int_distribution<int> rebellionroll(0, 99);

    int currentyear = -years;

    for (int tick = 0; tick < ticks; tick++, currentyear += tickyears)
    {
        recalcpolitystats(settlements, polities, world.width());

        const size_t initialpolitycount = polities.size();

        for (size_t polityindex = 0; polityindex < initialpolitycount; polityindex++)
        {
            polities[polityindex].cohesion = std::clamp(polities[polityindex].cohesion + randomdrift(rng) - polities[polityindex].militaryPressure * 0.05f, 5.0f, 95.0f);

            if (polities[polityindex].cohesion >= 30.0f)
                continue;

            if (rebellionroll(rng) > 28)
                continue;

            Settlement* frontier = nullptr;
            float bestdistance = -1.0f;

            for (Settlement& settlement : settlements)
            {
                if (settlement.polityId != polities[polityindex].id || settlement.id == polities[polityindex].capitalSettlementId)
                    continue;

                const float distance = settlementdistance(settlement, settlements[polities[polityindex].capitalSettlementId], world.width());

                if (distance > bestdistance)
                {
                    bestdistance = distance;
                    frontier = &settlement;
                }
            }

            if (frontier == nullptr)
                continue;

            Polity secession;
            secession.id = static_cast<int>(polities.size());
            secession.capitalSettlementId = frontier->id;
            secession.name = "Successor " + to_string(secession.id + 1);
            secession.cohesion = std::clamp(45.0f + randomdrift(rng), 20.0f, 75.0f);
            const int oldpolityid = polities[polityindex].id;
            const string oldpolityname = polities[polityindex].name;
            polities.push_back(std::move(secession));
            frontier->polityId = polities.back().id;
            polities[polityindex].cohesion = std::min(95.0f, polities[polityindex].cohesion + 8.0f);

            HistoryEvent rebellion;
            rebellion.year = currentyear;
            rebellion.type = "rebellion";
            rebellion.primarySettlementId = frontier->id;
            rebellion.primaryPolityId = frontier->polityId;
            rebellion.secondaryPolityId = oldpolityid;
            rebellion.summary = frontier->name + " broke away from " + oldpolityname;
            history.push_back(std::move(rebellion));
        }
    }

    recalcpolitystats(settlements, polities, world.width());

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
        {
            const int owner = world.ownersettlementid(x, y);

            if (owner >= 0 && owner < static_cast<int>(settlements.size()))
                world.setownerpolityid(x, y, settlements[owner].polityId);
            else
                world.setownerpolityid(x, y, -1);
        }
    }
}
}
