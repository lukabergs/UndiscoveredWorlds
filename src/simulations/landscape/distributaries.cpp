#include "simulations/landscape/distributaries.hpp"

#include <cmath>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <string>
#include <vector>

namespace landscape
{
namespace
{
constexpr std::size_t none = static_cast<std::size_t>(-1);

void checksize(std::size_t actual, std::size_t expected, const char* name)
{
    if (actual != expected)
        throw std::invalid_argument(std::string("splitdistributaries: ") + name + " does not have metric.size() elements");
}

struct Link
{
    std::size_t to;
    double weight;
};
}

Distributaries splitdistributaries(const SphereMetric& metric, const std::vector<uint8_t>& sea,
                                   const std::vector<uint8_t>& prograded, const std::vector<float>& surface,
                                   double sealevel, const std::vector<uint8_t>& flowdirection,
                                   const std::vector<int32_t>& januaryflow, const std::vector<int32_t>& julyflow,
                                   double minimumflow)
{
    const std::size_t size = metric.size();
    checksize(sea.size(), size, "sea");
    checksize(prograded.size(), size, "prograded");
    checksize(surface.size(), size, "surface");
    checksize(flowdirection.size(), size, "flowdirection");
    checksize(januaryflow.size(), size, "januaryflow");
    checksize(julyflow.size(), size, "julyflow");

    Distributaries result;
    result.direction.assign(size, 0);
    result.january.assign(size, 0);
    result.july.assign(size, 0);

    auto isdelta = [&](std::size_t cell) { return sea[cell] == 0 && prograded[cell] != 0; };
    auto receiverof = [&](std::size_t cell) {
        const int code = flowdirection[cell];
        if (code < 1 || code > directioncount)
            return none;
        const std::int32_t next = metric.neighbourid(cell, code);
        return next < 0 ? none : static_cast<std::size_t>(next);
    };

    std::vector<double> january(size, 0.0), july(size, 0.0), largest(size, 0.0);
    std::vector<std::size_t> supplier(size, none);
    auto deliver = [&](std::size_t from, std::size_t to, double jan, double jul) {
        if (sea[to] != 0)
        {
            result.seajanuary += jan;
            result.seajuly += jul;
            return;
        }
        january[to] += jan;
        july[to] += jul;
        if (jan + jul > largest[to])
        {
            largest[to] = jan + jul;
            supplier[to] = from;
        }
    };

    // Rivers entering the prograded land.
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (sea[cell] != 0 || isdelta(cell))
            continue;
        const std::size_t next = receiverof(cell);
        if (next == none || !isdelta(next))
            continue;
        const double jan = static_cast<double>(januaryflow[cell]);
        const double jul = static_cast<double>(julyflow[cell]);
        result.inflowjanuary += jan;
        result.inflowjuly += jul;
        deliver(cell, next, jan, jul);
    }

    // The links of every prograded land cell and the in-degrees among them.
    std::vector<std::vector<Link>> links(size);
    std::vector<int> indegree(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (!isdelta(cell))
            continue;
        const double height = static_cast<double>(surface[cell]);
        std::vector<Link>& out = links[cell];
        for (int code = 1; code <= directioncount; code++)
        {
            const std::int32_t neighbour = metric.neighbourid(cell, code);
            if (neighbour < 0)
                continue;
            const std::size_t next = static_cast<std::size_t>(neighbour);
            if (sea[next] == 0 && !isdelta(next))
                continue;
            const double below = sea[next] != 0 ? sealevel : static_cast<double>(surface[next]);
            const double drop = height - below;
            if (drop > 0.0)
                out.push_back({ next, drop });
        }
        if (out.empty())
        {
            const std::size_t next = receiverof(cell);
            if (next != none
                && (sea[next] != 0 || (isdelta(next) && static_cast<double>(surface[next]) <= height)))
                out.push_back({ next, 1.0 });
        }
        for (const Link& link : out)
            if (sea[link.to] == 0)
                indegree[link.to]++;
    }

    // Kahn's order: a cell is split once every cell linking into it has been.
    std::deque<std::size_t> ready;
    std::vector<uint8_t> done(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
        if (isdelta(cell) && indegree[cell] == 0)
            ready.push_back(cell);
    while (!ready.empty())
    {
        const std::size_t cell = ready.front();
        ready.pop_front();
        done[cell] = 1;
        const std::vector<Link>& out = links[cell];
        double total = 0.0;
        for (const Link& link : out)
            total += link.weight;
        if (out.empty())
        {
            result.strandedjanuary += january[cell];
            result.strandedjuly += july[cell];
        }
        for (const Link& link : out)
        {
            const double share = link.weight / total;
            deliver(cell, link.to, january[cell] * share, july[cell] * share);
            if (sea[link.to] == 0 && --indegree[link.to] == 0)
                ready.push_back(link.to);
        }
    }
    // A cell left out by a cycle (none on a filled surface) keeps its flow: stranded.
    for (std::size_t cell = 0; cell < size; cell++)
        if (isdelta(cell) && !done[cell])
        {
            result.strandedjanuary += january[cell];
            result.strandedjuly += july[cell];
        }

    for (std::size_t cell = 0; cell < size; cell++)
    {
        // A cell without a link onward, or left out by a cycle, is no distributary.
        if (!isdelta(cell) || supplier[cell] == none || links[cell].empty() || !done[cell])
            continue;
        const long jan = std::lround(january[cell]);
        const long jul = std::lround(july[cell]);
        if (static_cast<double>(jan + jul) / 2.0 < minimumflow)
            continue;
        int towards = 0;
        for (int code = 1; code <= directioncount && towards == 0; code++)
        {
            const std::int32_t next = metric.neighbourid(cell, code);
            if (next >= 0 && static_cast<std::size_t>(next) == supplier[cell])
                towards = code;
        }
        result.direction[cell] = static_cast<int16_t>(towards);
        result.january[cell] = static_cast<int32_t>(jan);
        result.july[cell] = static_cast<int32_t>(jul);
    }
    return result;
}
}
