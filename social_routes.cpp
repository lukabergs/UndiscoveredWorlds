#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "functions.hpp"
#include "planet.hpp"
#include "social_generation_internal.hpp"

using namespace std;

namespace
{
struct OpenNode
{
    int x = 0;
    int y = 0;
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

float heuristiccost(int x, int y, int endx, int endy, int width)
{
    const int dx = wrappedxdistance(x, endx, width);
    const int dy = std::abs(y - endy);
    return std::sqrt(static_cast<float>(dx * dx + dy * dy)) * 10.0f;
}

template <typename PassableFn, typename CostFn>
vector<pair<int, int>> findpathastar(planet& world, int startx, int starty, int endx, int endy, int maxexpanded, PassableFn&& passable, CostFn&& movecost)
{
    const int width = world.width();
    const int height = world.height();
    const int cellheight = height + 1;
    const auto indexof = [cellheight](int x, int y) -> int
    {
        return x * cellheight + y;
    };

    const int startindex = indexof(startx, starty);
    const int endindex = indexof(endx, endy);
    priority_queue<OpenNode, vector<OpenNode>, OpenNodeCompare> open;
    unordered_map<int, float> bestg;
    unordered_map<int, int> camefrom;
    unordered_set<int> closed;
    bestg.reserve(8192);
    camefrom.reserve(8192);
    closed.reserve(8192);

    open.push({ startx, starty, 0.0f, heuristiccost(startx, starty, endx, endy, width) });
    bestg.emplace(startindex, 0.0f);

    int expanded = 0;

    while (!open.empty() && expanded < maxexpanded)
    {
        const OpenNode node = open.top();
        open.pop();

        const int nodeindex = indexof(node.x, node.y);
        const auto bestfornode = bestg.find(nodeindex);

        if (bestfornode != bestg.end() && node.g > bestfornode->second)
            continue;

        if (closed.find(nodeindex) != closed.end())
            continue;

        closed.emplace(nodeindex);
        expanded++;

        if (nodeindex == endindex)
            break;

        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                if (dx == 0 && dy == 0)
                    continue;

                const int nx = wrapxlocal(node.x + dx, width);
                const int ny = node.y + dy;

                if (ny < 0 || ny > height)
                    continue;

                if (!passable(node.x, node.y, nx, ny))
                    continue;

                const int neighbourindex = indexof(nx, ny);

                if (closed.find(neighbourindex) != closed.end())
                    continue;

                const float movement = movecost(node.x, node.y, nx, ny);
                const float tentativeg = node.g + movement;
                const auto existing = bestg.find(neighbourindex);

                if (existing == bestg.end() || tentativeg < existing->second)
                {
                    bestg[neighbourindex] = tentativeg;
                    camefrom[neighbourindex] = nodeindex;
                    const float f = tentativeg + heuristiccost(nx, ny, endx, endy, width);
                    open.push({ nx, ny, tentativeg, f });
                }
            }
        }
    }

    if (camefrom.find(endindex) == camefrom.end() && startindex != endindex)
        return {};

    vector<pair<int, int>> reversed;
    reversed.push_back({ endx, endy });
    int current = endindex;

    while (current != startindex)
    {
        const auto found = camefrom.find(current);

        if (found == camefrom.end())
            return {};

        current = found->second;
        const int x = current / cellheight;
        const int y = current % cellheight;
        reversed.push_back({ x, y });
    }

    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

vector<pair<int, int>> buildstraightpath(planet& world, int startx, int starty, int endx, int endy)
{
    const int width = world.width();
    int adjustedendx = endx;
    const int worldspan = width + 1;
    const int directdx = endx - startx;

    if (std::abs(directdx) > worldspan / 2)
    {
        if (directdx > 0)
            adjustedendx -= worldspan;
        else
            adjustedendx += worldspan;
    }

    vector<pair<int, int>> path;
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
        const int wrappedx = wrapxlocal(x0, width);
        const int clampedy = std::clamp(y0, 0, world.height());
        path.push_back({ wrappedx, clampedy });

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
    vector<Settlement>& settlements = world.settlements();
    vector<RouteEdge>& routes = world.routeedges();
    routes.clear();

    if (settlements.size() < 2)
        return;

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
            world.setroutetraffic(x, y, 0);
    }

    unordered_set<uint64_t> edgekeys;
    edgekeys.reserve(settlements.size() * 8);

    auto applypathtraffic = [&](const vector<pair<int, int>>& path, RouteMode mode, float traffic)
    {
        const int trafficadd = std::clamp(static_cast<int>(std::round(traffic * 0.18f)), 1, 22);

        for (const pair<int, int>& point : path)
        {
            const int x = point.first;
            const int y = point.second;
            world.setroutetraffic(x, y, std::min(100, world.routetraffic(x, y) + trafficadd));

            if (mode != RouteMode::sea && !world.sea(x, y))
                world.setinfrastructure(x, y, std::min(100, world.infrastructure(x, y) + trafficadd / 4));
        }
    };

    auto addroute = [&](Settlement& from, Settlement& to, RouteMode mode, const vector<pair<int, int>>& path)
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
            const float distance = std::sqrt(static_cast<float>(wrappedxdistance(source.x, target.x, world.width()) * wrappedxdistance(source.x, target.x, world.width()) + (source.y - target.y) * (source.y - target.y)));
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

            const auto landpassable = [&](int, int, int nx, int ny) -> bool
            {
                return !world.sea(nx, ny);
            };

            const auto landmovecost = [&](int cx, int cy, int nx, int ny) -> float
            {
                float cost = (cx == nx || cy == ny) ? 10.0f : 14.0f;
                const int elevationdelta = std::abs(world.map(cx, cy) - world.map(nx, ny));
                cost += std::min(12, elevationdelta / 180);
                const int special = world.special(nx, ny);

                if (special >= 130 && special <= 132)
                    cost += 6.0f;

                const bool currentriver = world.riveraveflow(cx, cy) > 250;
                const bool nextriver = world.riveraveflow(nx, ny) > 250;

                if (currentriver != nextriver)
                    cost += 5.0f;

                if (world.infrastructure(nx, ny) > 30)
                    cost = std::max(2.0f, cost - 1.8f);

                return cost;
            };

            vector<pair<int, int>> landpath = findpathastar(world, source.x, source.y, target.x, target.y, 85000, landpassable, landmovecost);

            if (landpath.empty())
                landpath = buildstraightpath(world, source.x, source.y, target.x, target.y);

            addroute(source, target, RouteMode::land, landpath);

            if (source.riverAccess >= 45.0f && target.riverAccess >= 45.0f && directdistance <= 500.0f)
            {
                const auto riverpassable = [&](int, int, int nx, int ny) -> bool
                {
                    return !world.sea(nx, ny) && world.riveraccess(nx, ny) >= 18;
                };

                const auto rivermovecost = [&](int cx, int cy, int nx, int ny) -> float
                {
                    float cost = (cx == nx || cy == ny) ? 8.0f : 11.0f;
                    cost -= static_cast<float>(std::min(world.riveraccess(nx, ny), 90)) * 0.04f;
                    cost = std::max(2.0f, cost);
                    return cost;
                };

                vector<pair<int, int>> riverpath = findpathastar(world, source.x, source.y, target.x, target.y, 70000, riverpassable, rivermovecost);

                if (!riverpath.empty())
                    addroute(source, target, RouteMode::river, riverpath);
            }

            if (source.harbor >= 55.0f && target.harbor >= 55.0f && directdistance <= 760.0f)
            {
                const auto seapassable = [&](int, int, int nx, int ny) -> bool
                {
                    if (world.sea(nx, ny) || world.coast(nx, ny))
                        return true;

                    const int startdist = wrappedxdistance(nx, source.x, world.width()) + std::abs(ny - source.y);
                    const int enddist = wrappedxdistance(nx, target.x, world.width()) + std::abs(ny - target.y);
                    return startdist <= 3 || enddist <= 3;
                };

                const auto seamovecost = [&](int cx, int cy, int nx, int ny) -> float
                {
                    float cost = (cx == nx || cy == ny) ? 7.0f : 9.5f;

                    if (!world.sea(nx, ny))
                        cost += 6.0f;

                    if (world.coast(nx, ny))
                        cost = std::max(2.0f, cost - 0.8f);

                    const float dirx = static_cast<float>(nx - cx);
                    const float diry = static_cast<float>(ny - cy);
                    const float length = std::max(1.0f, std::sqrt(dirx * dirx + diry * diry));
                    const float ux = dirx / length;
                    const float uy = diry / length;
                    const float currentu = static_cast<float>(world.seasonalcurrentu(seasonjuly, cx, cy));
                    const float currentv = static_cast<float>(world.seasonalcurrentv(seasonjuly, cx, cy));
                    const float assist = (currentu * ux + currentv * uy) / 120.0f;
                    cost *= std::clamp(1.0f - assist * 0.08f, 0.72f, 1.40f);
                    return std::max(2.0f, cost);
                };

                vector<pair<int, int>> seapath = findpathastar(world, source.x, source.y, target.x, target.y, 95000, seapassable, seamovecost);

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

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
        {
            if (world.sea(x, y))
                continue;

            const int traffic = world.routetraffic(x, y);
            int adjustment = static_cast<int>(std::round(profile.routeWeight * (static_cast<float>(traffic) / 6.0f)));

            if (traffic > 60)
                adjustment += 4;

            if (traffic > 52 && world.riveraccess(x, y) > 35)
                adjustment += 3;

            if (traffic > 52 && world.harborscore(x, y) > 35)
                adjustment += 3;

            world.setsettlementsuitability(x, y, std::clamp(world.settlementsuitability(x, y) + adjustment, 0, 100));
        }
    }
}
}
