#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/mesh/mesh_geometry.hpp"
#include "world/planet.hpp"
#include "simulations/society/social_generation_internal.hpp"

using namespace std;

namespace
{
using mesh::CellId;
using socialgen::scanslots;

struct OpenNode
{
    CellId cell = 0;
    float g = 0.0f;
    float f = 0.0f;
};

struct OpenNodeCompare
{
    bool operator()(const OpenNode& left, const OpenNode& right) const
    {
        return left.f > right.f;
    }
};

float heuristiccost(const mesh::CellGrid& grid, CellId cell, CellId end)
{
    return socialgen::celldistance(grid, cell, end) * 10.0f;
}

// A* over the neighbour table; passable and movecost take (from, to, slot).
template <typename PassableFn, typename CostFn>
vector<CellId> findpathastar(const planet& world, CellId start, CellId end, int maxexpanded, PassableFn&& passable, CostFn&& movecost)
{
    const mesh::CellGrid& grid = world.grid();
    priority_queue<OpenNode, vector<OpenNode>, OpenNodeCompare> open;
    unordered_map<CellId, float> bestg;
    unordered_map<CellId, CellId> camefrom;
    unordered_set<CellId> closed;
    bestg.reserve(8192);
    camefrom.reserve(8192);
    closed.reserve(8192);

    open.push({ start, 0.0f, heuristiccost(grid, start, end) });
    bestg.emplace(start, 0.0f);

    int expanded = 0;

    while (!open.empty() && expanded < maxexpanded)
    {
        const OpenNode node = open.top();
        open.pop();

        const auto bestfornode = bestg.find(node.cell);

        if (bestfornode != bestg.end() && node.g > bestfornode->second)
            continue;

        if (closed.find(node.cell) != closed.end())
            continue;

        closed.emplace(node.cell);
        expanded++;

        if (node.cell == end)
            break;

        for (int slot : scanslots)
        {
            const int neighbour = grid.neighbour(node.cell, slot);

            if (neighbour < 0)
                continue;

            const CellId n = static_cast<CellId>(neighbour);

            if (!passable(node.cell, n, slot))
                continue;

            if (closed.find(n) != closed.end())
                continue;

            const float movement = movecost(node.cell, n, slot);
            const float tentativeg = node.g + movement;
            const auto existing = bestg.find(n);

            if (existing == bestg.end() || tentativeg < existing->second)
            {
                bestg[n] = tentativeg;
                camefrom[n] = node.cell;
                const float f = tentativeg + heuristiccost(grid, n, end);
                open.push({ n, tentativeg, f });
            }
        }
    }

    if (camefrom.find(end) == camefrom.end() && start != end)
        return {};

    vector<CellId> reversed;
    reversed.push_back(end);
    CellId current = end;

    while (current != start)
    {
        const auto found = camefrom.find(current);

        if (found == camefrom.end())
            return {};

        current = found->second;
        reversed.push_back(current);
    }

    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

// The cells along the great circle between two cells, when no passable path exists.
vector<CellId> buildstraightpath(const mesh::CellGrid& grid, CellId start, CellId end)
{
    // Antipodal centres have no single arc: no straight path.
    if (mesh::greatcircle_rad(grid.centre(start), grid.centre(end)) > 3.14159265358979323846 - 1e-9)
        return {};

    return mesh::cellsalongarc(grid, grid.centre(start), grid.centre(end));
}

// The compass direction of a step for the sea-current assist: east and, as the
// former row offsets counted it, south.
void stepdirection(const mesh::CellGrid& grid, CellId from, int slot, float& ux, float& uy)
{
    const double bearing = grid.bearing_rad(from, slot);
    ux = static_cast<float>(std::sin(bearing));
    uy = static_cast<float>(-std::cos(bearing));
}

// A step's length in cell widths.
float steplength(const mesh::CellGrid& grid, CellId from, int slot)
{
    return static_cast<float>(grid.distance_m(from, slot) / socialgen::cellwidth_m(grid));
}

// A path's extent in cell widths: centre to centre along the path plus half a
// cell at each end (the former count of its cells on a straight row).
float pathlength(const mesh::CellGrid& grid, const vector<CellId>& path)
{
    double length_m = 0.0;

    for (size_t k = 1; k < path.size(); k++)
        length_m += mesh::greatcircle_m(grid.centre(path[k - 1]), grid.centre(path[k]), grid.radius_m());

    return 1.0f + static_cast<float>(length_m / socialgen::cellwidth_m(grid));
}

float settlementpopulation(const Settlement& settlement)
{
    return static_cast<float>(settlement.urbanPopulation + settlement.ruralPopulation);
}

uint64_t makeedgekey(int first, int second, RouteMode mode)
{
    const int low = std::min(first, second);
    const int high = std::max(first, second);
    return (static_cast<uint64_t>(static_cast<uint32_t>(low)) << 32) ^
        (static_cast<uint64_t>(static_cast<uint32_t>(high)) << 8) ^
        static_cast<uint64_t>(static_cast<uint8_t>(mode));
}
}

namespace socialgen
{
void buildroutes(planet& world, const SocialProfile& profile, std::mt19937_64&)
{
    const mesh::CellGrid& grid = world.grid();
    vector<Settlement>& settlements = world.settlements();
    vector<RouteEdge>& routes = world.routeedges();
    routes.clear();

    if (settlements.size() < 2)
        return;

    for (CellId c = 0; c < grid.cellcount(); c++)
        world.setroutetraffic(c, 0);

    unordered_set<uint64_t> edgekeys;
    edgekeys.reserve(settlements.size() * 8);

    auto applypathtraffic = [&](const vector<CellId>& path, RouteMode mode, float traffic)
    {
        const int trafficadd = std::clamp(static_cast<int>(std::round(traffic * 0.18f)), 1, 22);

        for (const CellId c : path)
        {
            world.setroutetraffic(c, std::min(100, world.routetraffic(c) + trafficadd));

            if (mode != RouteMode::sea && !world.sea(c))
                world.setinfrastructure(c, std::min(100, world.infrastructure(c) + trafficadd / 4));
        }
    };

    auto addroute = [&](Settlement& from, Settlement& to, RouteMode mode, const vector<CellId>& path)
    {
        if (path.size() < 2)
            return;

        const uint64_t edgekey = makeedgekey(from.id, to.id, mode);

        if (edgekeys.find(edgekey) != edgekeys.end())
            return;

        edgekeys.emplace(edgekey);

        const float distance = pathlength(grid, path);
        const float populationfactor = std::sqrt(std::max(1.0f, settlementpopulation(from) * settlementpopulation(to)));
        const float modemultiplier = mode == RouteMode::land ? 1.0f : (mode == RouteMode::river ? 1.22f : 1.35f);
        const float costmultiplier = mode == RouteMode::land ? 1.0f : (mode == RouteMode::river ? 0.76f : 0.64f);
        const float cost = std::max(1.0f, distance * costmultiplier);
        const float capacity = std::clamp(populationfactor / 38.0f * modemultiplier, 2.0f, 120.0f);
        const float traffic = std::clamp(capacity / std::max(1.0f, cost * 0.085f), 1.0f, 100.0f);

        RouteEdge route;
        route.fromSettlementId = from.id;
        route.toSettlementId = to.id;
        route.mode = mode;
        route.cost = cost;
        route.capacity = capacity;
        route.traffic = traffic;
        routes.push_back(route);

        applypathtraffic(path, mode, traffic);
    };

    for (size_t i = 0; i < settlements.size(); i++)
    {
        Settlement& source = settlements[i];
        vector<pair<float, int>> neighbours;
        neighbours.reserve(settlements.size());

        for (size_t j = 0; j < settlements.size(); j++)
        {
            if (i == j)
                continue;

            const Settlement& target = settlements[j];
            const float distance = celldistance(grid, source.cell, target.cell);
            neighbours.push_back({ distance, static_cast<int>(j) });
        }

        std::sort(neighbours.begin(), neighbours.end(), [](const pair<float, int>& left, const pair<float, int>& right)
        {
            return left.first < right.first;
        });

        const int maxneighbours = std::min(5, static_cast<int>(neighbours.size()));

        for (int n = 0; n < maxneighbours; n++)
        {
            Settlement& target = settlements[neighbours[n].second];
            const float directdistance = neighbours[n].first;

            if (directdistance > 460.0f)
                continue;

            const auto landpassable = [&](CellId, CellId next, int) -> bool
            {
                return !world.sea(next);
            };

            const auto landmovecost = [&](CellId current, CellId next, int slot) -> float
            {
                float cost = 10.0f * steplength(grid, current, slot);
                const float elevationdelta = std::fabs(world.mapfloat(current) - world.mapfloat(next));
                cost += std::min(12.0f, elevationdelta / 180.0f);
                const int special = world.special(next);

                if (special >= 130 && special <= 132)
                    cost += 6.0f;

                const bool currentriver = world.riveraveflow(current) > 250;
                const bool nextriver = world.riveraveflow(next) > 250;

                if (currentriver != nextriver)
                    cost += 5.0f;

                if (world.infrastructure(next) > 30)
                    cost = std::max(2.0f, cost - 1.8f);

                return cost;
            };

            vector<CellId> landpath = findpathastar(world, source.cell, target.cell, 85000, landpassable, landmovecost);

            if (landpath.empty())
                landpath = buildstraightpath(grid, source.cell, target.cell);

            addroute(source, target, RouteMode::land, landpath);

            if (source.riverAccess >= 45.0f && target.riverAccess >= 45.0f && directdistance <= 500.0f)
            {
                const auto riverpassable = [&](CellId, CellId next, int) -> bool
                {
                    return !world.sea(next) && world.riveraccess(next) >= 18;
                };

                const auto rivermovecost = [&](CellId current, CellId next, int slot) -> float
                {
                    float cost = 8.0f * steplength(grid, current, slot);
                    cost -= static_cast<float>(std::min(world.riveraccess(next), 90)) * 0.04f;
                    cost = std::max(2.0f, cost);
                    return cost;
                };

                vector<CellId> riverpath = findpathastar(world, source.cell, target.cell, 70000, riverpassable, rivermovecost);

                if (!riverpath.empty())
                    addroute(source, target, RouteMode::river, riverpath);
            }

            if (source.harbor >= 55.0f && target.harbor >= 55.0f && directdistance <= 760.0f)
            {
                const auto iscoast = [&](CellId c) -> bool
                {
                    // planet::coast: a sea cell with land across an edge, polar rows excluded on lat-lon.
                    if (!world.sea(c))
                        return false;

                    if (grid.kind() == mesh::GridKind::latlon && (grid.row(c) == 0 || grid.row(c) == grid.rows() - 1))
                        return false;

                    for (int slot = 1; slot <= 7; slot += 2)
                    {
                        const int n = grid.neighbour(c, slot);

                        if (n >= 0 && !world.sea(static_cast<CellId>(n)))
                            return true;
                    }

                    return false;
                };

                const auto seapassable = [&](CellId, CellId next, int) -> bool
                {
                    if (world.sea(next) || iscoast(next))
                        return true;

                    return celldistance(grid, next, source.cell) <= 3.0f || celldistance(grid, next, target.cell) <= 3.0f;
                };

                const auto seamovecost = [&](CellId current, CellId next, int slot) -> float
                {
                    float cost = 7.0f * steplength(grid, current, slot);

                    if (!world.sea(next))
                        cost += 6.0f;

                    if (iscoast(next))
                        cost = std::max(2.0f, cost - 0.8f);

                    float ux = 0.0f;
                    float uy = 0.0f;
                    stepdirection(grid, current, slot, ux, uy);
                    const float currentu = static_cast<float>(world.seasonalcurrentu(seasonjuly, current));
                    const float currentv = static_cast<float>(world.seasonalcurrentv(seasonjuly, current));
                    const float assist = (currentu * ux + currentv * uy) / 120.0f;
                    cost *= std::clamp(1.0f - assist * 0.08f, 0.72f, 1.40f);
                    return std::max(2.0f, cost);
                };

                vector<CellId> seapath = findpathastar(world, source.cell, target.cell, 95000, seapassable, seamovecost);

                if (!seapath.empty())
                    addroute(source, target, RouteMode::sea, seapath);
            }
        }
    }

    for (Settlement& settlement : settlements)
    {
        float incidenttraffic = 0.0f;

        for (const RouteEdge& route : routes)
        {
            if (route.fromSettlementId == settlement.id || route.toSettlementId == settlement.id)
                incidenttraffic += route.traffic;
        }

        settlement.marketStrength = std::clamp(settlement.marketStrength + incidenttraffic * 0.35f, 0.0f, 100.0f);
        const int urbanboost = std::min(settlement.ruralPopulation / 4, static_cast<int>(std::round(incidenttraffic * 14.0f)));
        settlement.urbanPopulation += urbanboost;
        settlement.ruralPopulation = std::max(0, settlement.ruralPopulation - urbanboost);
    }

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (world.sea(c))
            continue;

        const int traffic = world.routetraffic(c);
        int adjustment = static_cast<int>(std::round(profile.routeWeight * (static_cast<float>(traffic) / 6.0f)));

        if (traffic > 60)
            adjustment += 4;

        if (traffic > 52 && world.riveraccess(c) > 35)
            adjustment += 3;

        if (traffic > 52 && world.harborscore(c) > 35)
            adjustment += 3;

        world.setsettlementsuitability(c, std::clamp(world.settlementsuitability(c) + adjustment, 0, 100));
    }
}
}
