#include "simulations/landscape/analytical_relief.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "simulations/landscape/sphere_metric.hpp"

// Erosion plan step E2c (docs/architecture/erosion-hydrology-plan.md, Phase E and
// section 6, "The solve"). With n = 1 and no uplift the stream power law
// dz/dt = -K A^m |grad z| is an advection of the start surface up the flow paths at
// the speed K Q^m (Tzathas et al. 2024, Eqs. 1, 5-9 with u = 0): z(x, t) = z0(D)
// with D the point whose response time is tau(x) - t, or the boundary's elevation
// when tau(x) <= t. Each step below is a pass over the cells in the flood's order
// (receivers before donors) and reads only the cell, its receiver path and its
// donors (in direction-code order), so the result does not depend on how ties in
// that order fall elsewhere in the patch.

namespace landscape
{
namespace
{
constexpr double cubicmetrespercubickilometre = 1.0e9;

struct Routing
{
    std::vector<int32_t> receiver;   // -1 on pinned cells
    std::vector<uint8_t> direction;  // 1..8 towards the receiver, 0 on pinned cells
    std::vector<int32_t> order;      // every cell once, each after its receiver
};

// The direction code from a cell to one of its neighbours, the lowest on a tie; 0 when none.
uint8_t directionto(const PatchMetric& metric, std::size_t cell, std::size_t target)
{
    for (int direction = 1; direction <= directioncount; direction++)
    {
        const std::size_t next = metric.neighbour(cell, direction);
        if (next != cell && next == target)
            return static_cast<uint8_t>(direction);
    }
    return 0;
}

// A diagonal step whose two corner cells (the cells orthogonally adjacent to both
// ends) are pinned passes between two cells of an 8-connected pinned line: refused,
// so that catchments stay on their side of the channel network.
bool crossespinnedcorner(const PatchMetric& metric, const std::vector<uint8_t>& pinned, std::size_t cell,
                         int direction)
{
    const int ox = directionoffsetx(direction), oy = directionoffsety(direction);
    if (ox == 0 || oy == 0)
        return false;
    const int x = static_cast<int>(cell % static_cast<std::size_t>(metric.width));
    const int y = static_cast<int>(cell / static_cast<std::size_t>(metric.width));
    return pinned[metric.index(x + ox, y)] != 0 && pinned[metric.index(x, y + oy)] != 0;
}

// Receivers on `surface`. A Priority-Flood seeded at the pinned cells (Barnes et al.
// 2014, Algorithm 2: a min-queue on (filled elevation, key) and a first-in first-out
// pit queue, as depressions.cpp, which ties by the local index instead) gives every
// cell a filled elevation, the neighbour it was reached from and the order in which
// it was taken; the key (the absolute fine index) makes that order the same in every
// window, and neighbours are visited in direction-code order. A cell drains
// by steepest descent (drop over metric distance, the lowest code on a tie) to a
// neighbour lower on both the surface and the filled surface; with none, to its
// flood receiver. The plain rule (steepest descent to any lower neighbour) loops:
// a rim cell descending into a pit whose flood receiver is that rim cell. Here a
// steepest-descent link lowers the filled elevation strictly and a flood link never
// raises it, while flood links alone form a tree towards the pinned cells, so no
// loop exists; the flood takes cells in non-decreasing filled elevation, so its order
// lists every receiver before its donors. On a surface without depressions the
// filled surface is the surface and every cell with a lower neighbour drains by
// steepest descent. Neither the flood nor the steepest descent takes a diagonal step
// between two pinned corner cells (crossespinnedcorner), so no water crosses an
// 8-connected channel; every cell is still reached, through orthogonal steps.
Routing route(const PatchMetric& metric, const std::vector<float>& surface, const std::vector<uint8_t>& pinned,
              const std::vector<std::uint64_t>& key)
{
    const std::size_t count = metric.size();
    Routing routing;
    routing.receiver.assign(count, -1);
    routing.direction.assign(count, 0);
    routing.order.reserve(count);

    std::vector<float> filled = surface;
    std::vector<int32_t> floodreceiver(count, -1);
    const auto keyof = [&](std::size_t cell) { return key.empty() ? static_cast<std::uint64_t>(cell) : key[cell]; };
    // (level, key, cell): the cell decides only between equal keys, which the contract excludes.
    using Entry = std::tuple<float, std::uint64_t, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    std::vector<uint8_t> closed(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (pinned[cell] == 0)
            continue;
        closed[cell] = 1;
        queue.push(Entry{ surface[cell], keyof(cell), static_cast<int32_t>(cell) });
    }
    std::vector<int32_t> pit;
    std::size_t pithead = 0;
    while (pithead < pit.size() || !queue.empty())
    {
        std::size_t cell;
        if (pithead < pit.size())
            cell = static_cast<std::size_t>(pit[pithead++]);
        else
        {
            cell = static_cast<std::size_t>(std::get<2>(queue.top()));
            queue.pop();
            pit.clear();
            pithead = 0;
        }
        routing.order.push_back(static_cast<int32_t>(cell));
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || closed[next] != 0 || crossespinnedcorner(metric, pinned, cell, direction))
                continue;
            closed[next] = 1;
            floodreceiver[next] = static_cast<int32_t>(cell);
            if (surface[next] < filled[cell])
            {
                filled[next] = filled[cell];
                pit.push_back(static_cast<int32_t>(next));
            }
            else
                queue.push(Entry{ surface[next], keyof(next), static_cast<int32_t>(next) });
        }
    }

    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (pinned[cell] != 0)
            continue;
        int best = 0;
        double steepest = 0.0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || !(surface[next] < surface[cell]) || !(filled[next] < filled[cell])
                || crossespinnedcorner(metric, pinned, cell, direction))
                continue;
            const double gradient = (static_cast<double>(surface[cell]) - static_cast<double>(surface[next]))
                                    / metric.distance_m(cell, direction);
            if (gradient > steepest)
            {
                steepest = gradient;
                best = direction;
            }
        }
        if (best != 0)
        {
            routing.receiver[cell] = static_cast<int32_t>(metric.neighbour(cell, best));
            routing.direction[cell] = static_cast<uint8_t>(best);
        }
        else
        {
            routing.receiver[cell] = floodreceiver[cell];
            routing.direction[cell] = directionto(metric, cell, static_cast<std::size_t>(floodreceiver[cell]));
        }
    }
    return routing;
}

// The solution on fixed receivers: discharge, response time, the start surface
// advected up the flow paths, then the threshold slope.
std::vector<float> solve(const PatchMetric& metric, const AnalyticalReliefInput& input,
                         const AnalyticalReliefParameters& parameters, const Routing& routing,
                         std::vector<double>& discharge)
{
    const std::size_t count = metric.size();
    const std::vector<float>& start = input.start_m;
    const std::vector<uint8_t>& pinned = input.pinned;
    const std::size_t width = static_cast<std::size_t>(metric.width);

    // Discharge: runoff x cell area, summed over each unpinned catchment, donors last
    // first; a cell adds its donors in direction-code order (a canonical sum).
    discharge.assign(count, 0.0);
    for (std::size_t k = routing.order.size(); k-- > 0;)
    {
        const std::size_t cell = static_cast<std::size_t>(routing.order[k]);
        if (pinned[cell] != 0)
            continue;
        double water = static_cast<double>(input.runoff_m_yr[cell]) * metric.rowarea_m2[cell / width];
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next != cell && routing.receiver[next] == static_cast<int32_t>(cell))
                water += discharge[next];
        }
        discharge[cell] = water;
    }

    // Response time (Myr) down to the first pinned cell: tau = tau_receiver + L / (K Q^m),
    // Q in km^3/yr as stream_power.cpp (a negative sum counts as 0); a link without
    // flow (Q^m = 0) never responds (infinite time).
    const double infinity = std::numeric_limits<double>::infinity();
    std::vector<double> tau(count, 0.0);
    std::vector<int32_t> outlet(count, -1);
    for (int32_t entry : routing.order)
    {
        const std::size_t cell = static_cast<std::size_t>(entry);
        if (pinned[cell] != 0)
        {
            outlet[cell] = entry;
            continue;
        }
        const std::size_t receiver = static_cast<std::size_t>(routing.receiver[cell]);
        const double length = metric.distance_m(cell, routing.direction[cell]);
        const double water = std::max(0.0, discharge[cell]) / cubicmetrespercubickilometre;
        const double speed = parameters.erodibility * std::pow(water, parameters.areaexponent);
        tau[cell] = tau[receiver] + (speed > 0.0 ? length / speed : infinity);
        outlet[cell] = outlet[receiver];
    }

    // The characteristic solution (Tzathas et al. 2024, Eqs. 8 and 9 with u = 0): with
    // tau <= t the pinned cell's elevation; otherwise the start elevation of the point
    // at response time s = tau - t downstream, linear in tau along the link a -> r
    // with tau(r) <= s < tau(a) (the paper's locally constant 1 / a). The walk covers
    // the links within time t of the cell. Neither value may exceed the cell's own
    // start: with u = 0 the erosion term K A^m |grad z| never raises a cell, while a
    // receiver path that climbs on the start surface (a depression's flood links, or
    // receivers taken from a later surface) would carry a higher z0 up to it.
    const double elapsed = parameters.elapsed_myr;
    std::vector<double> advected(count, 0.0);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const double own = static_cast<double>(start[cell]);
        if (pinned[cell] != 0)
        {
            advected[cell] = own;
            continue;
        }
        if (tau[cell] <= elapsed)
        {
            advected[cell] = std::min(own, static_cast<double>(start[static_cast<std::size_t>(outlet[cell])]));
            continue;
        }
        const double target = tau[cell] - elapsed;
        std::size_t upper = cell;
        for (;;)
        {
            const std::size_t lower = static_cast<std::size_t>(routing.receiver[upper]);
            if (tau[lower] <= target)
            {
                const double zupper = static_cast<double>(start[upper]);
                const double zlower = static_cast<double>(start[lower]);
                if (target >= tau[upper])
                    advected[cell] = zupper;  // t = 0, or a link that never responds
                else
                    advected[cell] = zlower + (zupper - zlower) * (target - tau[lower]) / (tau[upper] - tau[lower]);
                break;
            }
            upper = lower;
        }
        advected[cell] = std::min(own, advected[cell]);
    }

    // Threshold slope (Tzathas et al. 2024, section 6): from the pinned cells upwards,
    // a cell above its receiver + S x L is lowered to it. The stored float is rounded
    // down when rounding to nearest would put it above the limit.
    std::vector<float> elevation(count, 0.0f);
    for (int32_t entry : routing.order)
    {
        const std::size_t cell = static_cast<std::size_t>(entry);
        if (pinned[cell] != 0)
        {
            elevation[cell] = start[cell];
            continue;
        }
        const std::size_t receiver = static_cast<std::size_t>(routing.receiver[cell]);
        const double limit = static_cast<double>(elevation[receiver])
                             + parameters.thresholdslope * metric.distance_m(cell, routing.direction[cell]);
        const double value = std::min(advected[cell], limit);
        float stored = static_cast<float>(value);
        if (static_cast<double>(stored) > limit)
            stored = std::nextafter(stored, -std::numeric_limits<float>::infinity());
        elevation[cell] = stored;
    }
    return elevation;
}
}

AnalyticalReliefResult solveanalyticalrelief(const PatchMetric& metric, const AnalyticalReliefInput& input,
                                             const AnalyticalReliefParameters& parameters)
{
    const std::size_t count = metric.size();
    const std::size_t rows = static_cast<std::size_t>(metric.height < 0 ? 0 : metric.height);
    if (metric.width < 0 || metric.height < 0 || metric.rowdx_m.size() != rows || metric.rowarea_m2.size() != rows
        || input.start_m.size() != count || input.pinned.size() != count || input.runoff_m_yr.size() != count
        || (!input.key.empty() && input.key.size() != count))
        throw std::invalid_argument("solveanalyticalrelief: input sizes differ from the metric's");
    if (!(parameters.erodibility > 0.0) || !(parameters.areaexponent >= 0.0) || !(parameters.elapsed_myr >= 0.0)
        || !std::isfinite(parameters.elapsed_myr) || !(parameters.thresholdslope > 0.0) || parameters.maxiterations < 1
        || !(parameters.damping > 0.0) || !(parameters.damping <= 1.0))
        throw std::invalid_argument("solveanalyticalrelief: a parameter is out of range");
    if (!input.key.empty())
    {
        std::vector<std::uint64_t> sorted = input.key;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument("solveanalyticalrelief: two cells share a key");
    }
    bool anypinned = false;
    for (uint8_t flag : input.pinned)
        anypinned = anypinned || flag != 0;
    if (count > 0 && !anypinned)
        throw std::invalid_argument("solveanalyticalrelief: no cell is pinned");

    AnalyticalReliefResult result;
    result.elevation_m = input.start_m;
    result.direction.assign(count, 0);
    result.discharge_m3_yr.assign(count, 0.0);
    // The routing surface: the start, then the exponential moving average of the
    // solutions (Tzathas et al. 2024, section 5.1), held in double and read as float.
    std::vector<double> average(input.start_m.begin(), input.start_m.end());
    std::vector<float> routed = input.start_m;
    const double damping = parameters.damping;
    // Always maxiterations passes: stopping when no receiver changes would make the
    // pass count, and so every cell's result, depend on cells anywhere in the window.
    Routing previous;
    for (int pass = 1; pass <= parameters.maxiterations; pass++)
    {
        Routing routing = route(metric, routed, input.pinned, input.key);
        result.iterations = pass;
        result.converged = pass > 1 && routing.receiver == previous.receiver;
        // Every pass solves from the start surface on the receivers of the damped surface.
        result.elevation_m = solve(metric, input, parameters, routing, result.discharge_m3_yr);
        result.direction = routing.direction;
        previous = std::move(routing);
        if (pass == parameters.maxiterations)
            break;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            average[cell] = damping * static_cast<double>(result.elevation_m[cell]) + (1.0 - damping) * average[cell];
            routed[cell] = static_cast<float>(average[cell]);
        }
    }
    return result;
}
}
