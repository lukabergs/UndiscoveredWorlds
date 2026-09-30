#include <algorithm>
#include <numeric>
#include <queue>
#include <vector>

#include "planet.hpp"
#include "social_generation_internal.hpp"

using namespace std;

namespace socialgen
{
void runtradeandknowledge(planet& world, const SocialTables& tables, std::mt19937_64&)
{
    vector<Settlement>& settlements = world.settlements();
    vector<Polity>& polities = world.polities();
    const vector<RouteEdge>& routes = world.routeedges();

    if (settlements.empty())
        return;

    vector<vector<pair<int, float>>> adjacency(settlements.size());

    for (const RouteEdge& route : routes)
    {
        if (route.fromSettlementId < 0 || route.toSettlementId < 0 || route.fromSettlementId >= static_cast<int>(settlements.size()) || route.toSettlementId >= static_cast<int>(settlements.size()))
            continue;

        adjacency[route.fromSettlementId].push_back({ route.toSettlementId, route.traffic });
        adjacency[route.toSettlementId].push_back({ route.fromSettlementId, route.traffic });
    }

    vector<int> component(settlements.size(), -1);
    vector<float> componentconnectivity;
    int componentcount = 0;

    for (int i = 0; i < static_cast<int>(settlements.size()); i++)
    {
        if (component[i] != -1)
            continue;

        queue<int> frontier;
        frontier.push(i);
        component[i] = componentcount;
        float connectivitysum = 0.0f;

        while (!frontier.empty())
        {
            const int current = frontier.front();
            frontier.pop();

            for (const pair<int, float>& edge : adjacency[current])
            {
                connectivitysum += edge.second;

                if (component[edge.first] == -1)
                {
                    component[edge.first] = componentcount;
                    frontier.push(edge.first);
                }
            }
        }

        componentconnectivity.push_back(connectivitysum);
        componentcount++;
    }

    float knowledgeproductivitybonus = 0.0f;
    float knowledgecohesionbonus = 0.0f;
    float knowledgespreadfactor = 0.0f;

    for (const KnowledgeDefinition& knowledge : tables.knowledges)
    {
        knowledgespreadfactor += knowledge.spreadFactor;
        knowledgeproductivitybonus += knowledge.productivityBonus;
        knowledgecohesionbonus += knowledge.cohesionBonus;
    }

    const float normalizedknowledgeproductivity = tables.knowledges.empty() ? 0.0f : knowledgeproductivitybonus / static_cast<float>(tables.knowledges.size());
    const float normalizedknowledgecohesion = tables.knowledges.empty() ? 0.0f : knowledgecohesionbonus / static_cast<float>(tables.knowledges.size());
    const float normalizedknowledgespread = tables.knowledges.empty() ? 0.0f : knowledgespreadfactor / static_cast<float>(tables.knowledges.size());

    const float commodityvalueweight = std::accumulate(tables.commodities.begin(), tables.commodities.end(), 0.0f, [](float value, const CommodityDefinition& commodity)
    {
        return value + commodity.baseValue * std::max(0.1f, commodity.weight);
    });

    const float recipedepth = tables.recipes.empty() ? 0.0f : static_cast<float>(tables.recipes.size()) * 0.35f;

    for (Settlement& settlement : settlements)
    {
        float routeconnectivity = 0.0f;

        for (const pair<int, float>& edge : adjacency[settlement.id])
            routeconnectivity += edge.second;

        const int localresources = std::max(world.agriculturalcapacity(settlement.x, settlement.y), world.riveraccess(settlement.x, settlement.y));
        const int reservebonus = std::max(world.fisheryreserve(settlement.x, settlement.y), std::max(world.metalorereserve(settlement.x, settlement.y), world.placerreserve(settlement.x, settlement.y)));
        const int componentid = component[settlement.id];
        const float componenttraffic = componentid >= 0 && componentid < static_cast<int>(componentconnectivity.size()) ? componentconnectivity[componentid] : routeconnectivity;
        const float scarcitypenalty = adjacency[settlement.id].empty() ? 14.0f : 0.0f;
        const float productivityboost = normalizedknowledgeproductivity + normalizedknowledgespread * std::min(1.0f, routeconnectivity / 80.0f);
        const float commodityboost = std::min(20.0f, commodityvalueweight * 0.0025f + recipedepth);

        const float marketstrength =
            settlement.marketStrength * 0.35f +
            routeconnectivity * 0.52f +
            static_cast<float>(localresources) * 0.24f +
            static_cast<float>(reservebonus) * 0.18f +
            componenttraffic * 0.08f +
            productivityboost +
            commodityboost -
            scarcitypenalty;

        settlement.marketStrength = std::clamp(marketstrength, 0.0f, 100.0f);

        const int urbanlift = std::min(settlement.ruralPopulation / 6, static_cast<int>(std::round(settlement.marketStrength * 0.8f)));
        settlement.urbanPopulation += urbanlift;
        settlement.ruralPopulation = std::max(0, settlement.ruralPopulation - urbanlift);
    }

    vector<float> politymarket(polities.size(), 0.0f);
    vector<int> politysettlementcounts(polities.size(), 0);

    for (const Settlement& settlement : settlements)
    {
        if (settlement.polityId < 0 || settlement.polityId >= static_cast<int>(polities.size()))
            continue;

        politymarket[settlement.polityId] += settlement.marketStrength;
        politysettlementcounts[settlement.polityId]++;
    }

    for (Polity& polity : polities)
    {
        if (polity.id < 0 || polity.id >= static_cast<int>(politysettlementcounts.size()))
            continue;

        const int settlementcount = std::max(1, politysettlementcounts[polity.id]);
        const float averagemarket = politymarket[polity.id] / static_cast<float>(settlementcount);
        polity.wealth = std::clamp(averagemarket * (1.15f + commodityvalueweight * 0.0008f), 0.0f, 200.0f);
        polity.infrastructure = std::clamp(polity.infrastructure + averagemarket * 0.08f, 0.0f, 100.0f);
        polity.cohesion = std::clamp(polity.cohesion + normalizedknowledgecohesion * 0.4f + averagemarket * 0.03f - polity.militaryPressure * 0.05f, 5.0f, 95.0f);
    }
}
}
