#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

// The straight line between two cells, when no passable path exists (lat-lon only).
vector<CellId> buildstraightpath(const mesh::CellGrid& grid, CellId start, CellId end)
{
    const int startx = grid.column(start);
    const int starty = grid.row(start);
    const int endx = grid.column(end);
    const int endy = grid.row(end);
    const int worldspan = grid.columns();
    int adjustedendx = endx;
    const int directdx = endx - startx;

    if (std::abs(directdx) > worldspan / 2)
    {
        if (directdx > 0)
            adjustedendx -= worldspan;
        else
            adjustedendx += worldspan;
    }

    vector<CellId> path;
    int x0 = startx;
    int y0 = starty;
    const int x1 = adjustedendx;
    const int y1 = endy;
    const int dx = std::abs(x1 - x0);
    const int sx = x0 < x1 ? 1 : -1;
    const int dy = -std::abs(y1 - y0);
    const int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;

    while (true)
    {
        int wrappedx = x0;
        while (wrappedx < 0)
            wrappedx += worldspan;
        while (wrappedx >= worldspan)
            wrappedx -= worldspan;
        const int clampedy = std::clamp(y0, 0, grid.rows() - 1);
        path.push_back(grid.cellof(wrappedx, clampedy));

        if (x0 == x1 && y0 == y1)
            break;

        const int e2 = 2 * error;

        if (e2 >= dy)
        {
            error += dy;
            x0 += sx;
        }

        if (e2 <= dx)
        {
            error += dx;
            y0 += sy;
        }
    }

    return path;
}

// The unit direction of a step for the sea-current assist (lat-lon only: the
// column and row offsets, y growing south, as the former loops computed them).
void stepdirection(const mesh::CellGrid& grid, CellId from, CellId to, float& ux, float& uy)
{
    const float dirx = static_cast<float>(grid.column(to) - grid.column(from));
    const float diry = static_cast<float>(grid.row(to) - grid.row(from));
    const float length = std::max(1.0f, std::sqrt(dirx * dirx + diry * diry));
    ux = dirx / length;
    uy = diry / length;
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

        const float distance = static_cast<float>(path.size());
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
                float cost = orthogonalslot(slot) ? 10.0f : 14.0f;
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

                const auto rivermovecost = [&](CellId, CellId next, int slot) -> float
                {
                    float cost = orthogonalslot(slot) ? 8.0f : 11.0f;
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

                    const int startdist = cellmanhattan(grid, next, source.cell);
                    const int enddist = cellmanhattan(grid, next, target.cell);
                    return startdist <= 3 || enddist <= 3;
                };

                const auto seamovecost = [&](CellId current, CellId next, int slot) -> float
                {
                    float cost = orthogonalslot(slot) ? 7.0f : 9.5f;

                    if (!world.sea(next))
                        cost += 6.0f;

                    if (iscoast(next))
                        cost = std::max(2.0f, cost - 0.8f);

                    float ux = 0.0f;
                    float uy = 0.0f;
                    stepdirection(grid, current, next, ux, uy);
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
