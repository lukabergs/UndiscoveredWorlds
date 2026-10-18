#include "simulations/landscape/flow_graph.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/depressions.hpp"

namespace landscape
{
namespace
{
// The direction code from a cell to one of its neighbours, the lowest on a tie; 0 when none.
int directionto(const SphereMetric& metric, std::size_t cell, std::size_t target)
{
    for (int direction = 1; direction <= directioncount; direction++)
    {
        const std::int32_t next = metric.neighbourid(cell, direction);
        if (next >= 0 && static_cast<std::size_t>(next) == target)
            return direction;
    }
    return 0;
}

// A uniform number in [0, 1) from a hash.
double unitdraw(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }

// The randomised receiver direction of a land cell that drains to a lower
// neighbour (`steepest` is the D8 choice, the fallback). D-infinity (Tarboton
// 1997) on the metric: each of the eight triangular facets formed by a cardinal
// neighbour, the adjacent diagonal one and the cell gives a plane through the
// three heights; its descent direction, when it lies between the facet's two
// edges, has the plane's slope, otherwise the facet's slope is that of its
// steeper edge. In the steepest facet, with p1 and p2 the offsets to the two edge
// neighbours (metres) and g the descent direction, the second edge's neighbour is
// drawn with probability |p1 x g| / (|p1 x g| + |p2 x g|), so the expected
// displacement (1 - f) p1 + f p2 is parallel to g (a river's heading is its mean
// displacement). A drawn neighbour that is neither lower nor sea gives way to the
// other edge's, then to `steepest`.
int randomiseddirection(const SphereMetric& metric, std::size_t cell, const std::vector<float>& filled,
                        const std::vector<uint8_t>& sea, float sealevel, std::uint64_t seed, int steepest)
{
    const double centre = static_cast<double>(filled[cell]);
    const auto height = [&](std::size_t next) {
        return sea[next] != 0 ? static_cast<double>(sealevel) : static_cast<double>(filled[next]);
    };
    double best = 0.0, fraction = 0.0;
    int first = 0, second = 0;
    for (int direction = 1; direction <= directioncount; direction++)
    {
        const int other = direction % directioncount + 1;
        const std::int32_t one = metric.neighbourid(cell, direction), two = metric.neighbourid(cell, other);
        if (one < 0 || two < 0)
            continue;
        const double drop1 = centre - height(static_cast<std::size_t>(one));
        const double drop2 = centre - height(static_cast<std::size_t>(two));
        const TangentOffset p1 = metric.tangentoffset(cell, direction), p2 = metric.tangentoffset(cell, other);
        // The descent vector v of the plane: v . p1 = drop1, v . p2 = drop2.
        const double det = p1.east * p2.north - p1.north * p2.east;
        const double ve = (drop1 * p2.north - drop2 * p1.north) / det;
        const double vn = (p1.east * drop2 - p2.east * drop1) / det;
        const double cross1 = p1.east * vn - p1.north * ve;  // p1 x v
        const double cross2 = ve * p2.north - vn * p2.east;  // v x p2
        double slope, share;
        if ((ve != 0.0 || vn != 0.0) && cross1 * det >= 0.0 && cross2 * det >= 0.0)
        {
            // (1 - share) p1 + share p2 is parallel to v: the mean displacement
            // follows the gradient.
            slope = std::hypot(ve, vn);
            share = std::clamp(std::fabs(cross1) / (std::fabs(cross1) + std::fabs(cross2)), 0.0, 1.0);
        }
        else
        {
            const double slope1 = drop1 / std::hypot(p1.east, p1.north);
            const double slope2 = drop2 / std::hypot(p2.east, p2.north);
            slope = std::max(slope1, slope2);
            share = slope2 > slope1 ? 1.0 : 0.0;
        }
        if (slope > best)
        {
            best = slope;
            first = direction;
            second = other;
            fraction = share;
        }
    }
    if (first == 0)
        return steepest;

    const bool drawsecond = unitdraw(deterministiccombine(seed, static_cast<std::int64_t>(cell))) < fraction;
    const auto accepted = [&](int direction) {
        const std::size_t next = static_cast<std::size_t>(metric.neighbourid(cell, direction));  // both exist
        return sea[next] != 0 || filled[next] < filled[cell];
    };
    const int drawn = drawsecond ? second : first, other = drawsecond ? first : second;
    if (accepted(drawn))
        return drawn;
    if (accepted(other))
        return other;
    return steepest;
}

// Sets a land cell's receiver to its neighbour in `direction`; a sea neighbour
// makes the cell an outlet.
void setreceiver(const SphereMetric& metric, const std::vector<uint8_t>& sea, std::size_t cell, int direction,
                 FlowGraph& graph)
{
    const std::size_t target = static_cast<std::size_t>(metric.neighbourid(cell, direction));  // an existing neighbour
    graph.direction[cell] = static_cast<uint8_t>(direction);
    graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, direction));
    graph.receiver[cell] = static_cast<int32_t>(cell);
    graph.outletsea[cell] = -1;
    if (sea[target] != 0)
        graph.outletsea[cell] = static_cast<int32_t>(target);
    else
        graph.receiver[cell] = static_cast<int32_t>(target);
}

// How buildflowgraph routed a land cell, for the flat resolution.
constexpr uint8_t routednone = 0;      // sea, or a land cell the flood never reached
constexpr uint8_t routedflood = 1;     // follows its flood receiver
constexpr uint8_t routedsteepest = 2;  // drains to a lower neighbour or the sea

// Flat resolution (Barnes, Lehman and Mulla 2014) of the cells that follow their
// flood receiver. Two kinds of flat, both at one filled elevation:
// - Raised cells (filled above the terrain): a cell's outlet is the first cell on
//   its flood-receiver chain that is not raised, so a depression's cells have its
//   spill as outlet; the raised cells that share an outlet form one flat (connected
//   through their chains) whose only low edge is the outlet.
// - Cells at their own terrain: a connected region of them at one filled
//   elevation is a flat whose low edges are the adjacent land cells at the same
//   elevation that drain to a lower neighbour or the sea. A region without one
//   (its flood chains lead into a raised flat) keeps its flood receivers.
// In each flat, the gradient away from higher terrain (breadth-first distance from
// the flat cells next to a higher land cell, inverted) and the gradient towards
// the low edges (twice the breadth-first distance from them) add up to an
// increment that decreases towards a low edge; a flat cell drains by steepest
// descent on it (decrease over the metric's centre distance, the lowest code on a
// tie) to a neighbour in its flat or a low edge, whose increment is 2. No cycle forms: a cell at its terrain draining through a
// low edge drops below the flat; the outlet of a raised flat lies before the flat's
// cells on their flood-receiver chains, as does the raised flat an unresolved
// region's chains lead into.
void resolveflatcells(const SphereMetric& metric, const std::vector<float>& elevation, const std::vector<uint8_t>& sea,
                      const Depressions& depressions, const std::vector<uint8_t>& routed, FlowGraph& graph)
{
    const std::size_t count = metric.size();
    const std::vector<float>& filled = depressions.filledelevation;
    const std::vector<int32_t>& floodreceiver = depressions.floodreceiver;
    constexpr int32_t unknown = -1, nooutlet = -2;

    // Kind: 0 not flat, 1 flat at its terrain, 2 flat raised by the fill.
    std::vector<uint8_t> kind(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        if (routed[cell] == routedflood)
            kind[cell] = filled[cell] > elevation[cell] ? 2 : 1;
    // A low edge of a flat cell: a land neighbour at the same filled elevation that
    // drains to a lower neighbour or the sea.
    const auto lowedge = [&](std::size_t cell, std::size_t neighbour) {
        return routed[neighbour] == routedsteepest && filled[neighbour] == filled[cell];
    };

    // The outlets of the raised cells.
    std::vector<int32_t> outlet(count, unknown);
    std::vector<std::size_t> path;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (kind[cell] != 2 || outlet[cell] != unknown)
            continue;
        path.clear();
        std::size_t current = cell;
        int32_t found = nooutlet;
        for (std::size_t step = 0; step <= count; step++)
        {
            path.push_back(current);
            const int32_t next = floodreceiver[current];
            if (next < 0)
                break;
            const std::size_t following = static_cast<std::size_t>(next);
            if (kind[following] != kind[current])
            {
                found = next;
                break;
            }
            if (outlet[following] != unknown)
            {
                found = outlet[following];
                break;
            }
            current = following;
        }
        for (std::size_t member : path)
            outlet[member] = found;
    }

    // Flats numbered in increasing cell index of their first cell: raised ones by
    // outlet, the others by connected region (eight neighbours, one filled elevation).
    std::vector<int32_t> flatid(count, -1);
    std::unordered_map<int32_t, int32_t> raisedids;
    std::vector<std::size_t> stack, region;
    int32_t flats = 0;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (kind[cell] == 2 && outlet[cell] >= 0)
        {
            const auto inserted = raisedids.emplace(outlet[cell], flats);
            if (inserted.second)
                flats++;
            flatid[cell] = inserted.first->second;
            continue;
        }
        if (kind[cell] != 1 || flatid[cell] != -1)
            continue;
        const int32_t id = flats++;
        bool drains = false;
        flatid[cell] = id;
        region.assign(1, cell);
        stack.assign(1, cell);
        while (!stack.empty())
        {
            const std::size_t current = stack.back();
            stack.pop_back();
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::int32_t next = metric.neighbourid(current, direction);
                if (next < 0)
                    continue;
                const std::size_t neighbour = static_cast<std::size_t>(next);
                drains = drains || lowedge(current, neighbour);
                if (kind[neighbour] != 1 || flatid[neighbour] != -1 || filled[neighbour] != filled[current])
                    continue;
                flatid[neighbour] = id;
                region.push_back(neighbour);
                stack.push_back(neighbour);
            }
        }
        if (!drains)
            for (std::size_t member : region)
                flatid[member] = -2;  // keeps its flood receivers; -2 so that it is not labelled again
    }
    for (std::size_t cell = 0; cell < count; cell++)
        if (flatid[cell] == -2)
            flatid[cell] = -1;
    if (flats == 0)
        return;
    // Is `neighbour` a low edge of the flat cell `cell`?
    const auto exit = [&](std::size_t cell, std::size_t neighbour) {
        return kind[cell] == 2 ? static_cast<int32_t>(neighbour) == outlet[cell] : lowedge(cell, neighbour);
    };

    // Breadth-first distances within each flat, from the seeds in `layer` (distance `start`).
    std::vector<int32_t> nextlayer;
    const auto spread = [&](std::vector<int32_t>& layer, std::vector<int32_t>& distance, int32_t start) {
        for (int32_t loop = start + 1; !layer.empty(); loop++)
        {
            nextlayer.clear();
            for (int32_t member : layer)
            {
                const std::size_t cell = static_cast<std::size_t>(member);
                for (int direction = 1; direction <= directioncount; direction++)
                {
                    const std::int32_t next = metric.neighbourid(cell, direction);
                    if (next < 0)
                        continue;
                    const std::size_t neighbour = static_cast<std::size_t>(next);
                    if (flatid[neighbour] != flatid[cell] || distance[neighbour] != 0)
                        continue;
                    distance[neighbour] = loop;
                    nextlayer.push_back(static_cast<int32_t>(neighbour));
                }
            }
            layer.swap(nextlayer);
        }
    };

    // Away from higher terrain: 1 on the cells next to a higher land cell.
    std::vector<int32_t> away(count, 0), towards(count, 0), layer;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (flatid[cell] < 0)
            continue;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::int32_t neighbour = metric.neighbourid(cell, direction);
            if (neighbour >= 0 && sea[static_cast<std::size_t>(neighbour)] == 0
                && filled[static_cast<std::size_t>(neighbour)] > filled[cell])
            {
                away[cell] = 1;
                layer.push_back(static_cast<int32_t>(cell));
                break;
            }
        }
    }
    spread(layer, away, 1);
    std::vector<int32_t> flatheight(static_cast<std::size_t>(flats), 0);
    for (std::size_t cell = 0; cell < count; cell++)
        if (flatid[cell] >= 0)
            flatheight[static_cast<std::size_t>(flatid[cell])]
                = std::max(flatheight[static_cast<std::size_t>(flatid[cell])], away[cell]);

    // Towards the low edges: a low edge is 1, its neighbours in the flat 2.
    layer.clear();
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (flatid[cell] < 0)
            continue;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::int32_t neighbour = metric.neighbourid(cell, direction);
            if (neighbour >= 0 && exit(cell, static_cast<std::size_t>(neighbour)))
            {
                towards[cell] = 2;
                layer.push_back(static_cast<int32_t>(cell));
                break;
            }
        }
    }
    spread(layer, towards, 2);

    // The combined increment (Barnes et al. 2014): a low edge 2, a flat cell twice its
    // distance towards the low edges plus its inverted distance away from higher terrain.
    const auto increment = [&](std::size_t cell) {
        const int32_t height = flatheight[static_cast<std::size_t>(flatid[cell])];
        return 2 * towards[cell] + (away[cell] > 0 ? height - away[cell] : 0);
    };
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (flatid[cell] < 0 || towards[cell] == 0)
            continue;  // not reached from a low edge: keeps its flood receiver
        const int32_t own = increment(cell);
        int best = 0;
        double steepest = 0.0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::int32_t next = metric.neighbourid(cell, direction);
            if (next < 0)
                continue;
            const std::size_t neighbour = static_cast<std::size_t>(next);
            int32_t value;
            if (exit(cell, neighbour))
                value = 2;
            else if (flatid[neighbour] == flatid[cell] && towards[neighbour] != 0)
                value = increment(neighbour);
            else
                continue;
            if (value >= own)
                continue;
            const double gradient = static_cast<double>(own - value) / metric.distance_m(cell, direction);
            if (gradient > steepest)
            {
                steepest = gradient;
                best = direction;
            }
        }
        if (best != 0)
            setreceiver(metric, sea, cell, best, graph);
    }
}
}

FlowGraph buildflowgraph(const SphereMetric& metric, const std::vector<float>& elevation,
                         const std::vector<uint8_t>& sea, float sealevel, const Depressions& depressions,
                         const RoutingOptions& options)
{
    const std::size_t count = metric.size();
    if (elevation.size() != count || sea.size() != count || depressions.filledelevation.size() != count
        || depressions.depressionid.size() != count || depressions.floodreceiver.size() != count)
        throw std::invalid_argument("buildflowgraph: input sizes differ from the metric");

    FlowGraph graph;
    graph.receiver.resize(count);
    graph.direction.assign(count, 0);
    graph.distance_m.assign(count, 0.0f);
    graph.outletsea.assign(count, -1);
    graph.land.resize(count);
    graph.area_m2.assign(count, 0.0);
    graph.discharge_m3_yr.assign(count, 0.0);

    const std::vector<float>& filled = depressions.filledelevation;
    std::vector<uint8_t> routed;  // how each land cell was routed, for resolveflats
    if (options.resolveflats)
        routed.assign(count, routednone);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        graph.receiver[cell] = static_cast<int32_t>(cell);
        if (sea[cell] != 0)
        {
            // Seabed: the lowest strictly lower sea neighbour, the lowest code on a tie.
            graph.land[cell] = 0;
            int best = 0;
            float lowest = elevation[cell];
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::int32_t neighbour = metric.neighbourid(cell, direction);
                if (neighbour < 0)
                    continue;
                const std::size_t next = static_cast<std::size_t>(neighbour);
                if (sea[next] == 0 || !(elevation[next] < lowest))
                    continue;
                lowest = elevation[next];
                best = direction;
            }
            if (best != 0)
            {
                graph.receiver[cell] = metric.neighbourid(cell, best);
                graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, best));
            }
            continue;
        }

        // Land: the steepest descent on the filled surface (equal to the terrain
        // outside depressions), a sea neighbour counting with the sea level; the
        // filled surface keeps a rim cell from draining back into its own depression.
        graph.land[cell] = 1;
        int best = 0;
        double steepest = 0.0;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::int32_t neighbour = metric.neighbourid(cell, direction);
            if (neighbour < 0)
                continue;
            const std::size_t next = static_cast<std::size_t>(neighbour);
            const double height = sea[next] != 0 ? static_cast<double>(sealevel) : static_cast<double>(filled[next]);
            const double gradient = (static_cast<double>(filled[cell]) - height) / metric.distance_m(cell, direction);
            if (gradient > steepest)
            {
                steepest = gradient;
                best = direction;
            }
        }

        std::size_t target = cell;
        if (best != 0
            && (sea[static_cast<std::size_t>(metric.neighbourid(cell, best))] != 0 || depressions.depressionid[cell] == 0))
        {
            if (options.randomised)
                best = randomiseddirection(metric, cell, filled, sea, sealevel, options.seed, best);
            target = static_cast<std::size_t>(metric.neighbourid(cell, best));
            if (options.resolveflats)
                routed[cell] = routedsteepest;
        }
        else if (depressions.floodreceiver[cell] >= 0)
        {
            target = static_cast<std::size_t>(depressions.floodreceiver[cell]);
            best = directionto(metric, cell, target);
            if (options.resolveflats)
                routed[cell] = routedflood;
        }
        else
            best = 0;  // not reached by the flood (no sea in its region): a pit, its own receiver

        if (best == 0)
            continue;
        graph.direction[cell] = static_cast<uint8_t>(best);
        graph.distance_m[cell] = static_cast<float>(metric.distance_m(cell, best));
        if (sea[target] != 0)
            graph.outletsea[cell] = static_cast<int32_t>(target);
        else
            graph.receiver[cell] = static_cast<int32_t>(target);
    }
    if (options.resolveflats)
        resolveflatcells(metric, elevation, sea, depressions, routed, graph);

    // The order: donor lists of the land cells (increasing index), then a
    // breadth-first walk (not Braun and Willett's depth-first stack) from each root
    // (an outlet, or a pit) in increasing index, `order` itself being the queue. It
    // has the properties the kernels need: every receiver before its donors, every
    // land cell once.
    std::vector<int32_t> donorstart(count + 1, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        if (graph.land[cell] != 0 && graph.receiver[cell] != static_cast<int32_t>(cell))
            donorstart[static_cast<std::size_t>(graph.receiver[cell]) + 1]++;
    for (std::size_t cell = 0; cell < count; cell++)
        donorstart[cell + 1] += donorstart[cell];
    std::vector<int32_t> donors(static_cast<std::size_t>(donorstart[count]));
    std::vector<int32_t> fill(donorstart.begin(), donorstart.end() - 1);
    for (std::size_t cell = 0; cell < count; cell++)
        if (graph.land[cell] != 0 && graph.receiver[cell] != static_cast<int32_t>(cell))
            donors[static_cast<std::size_t>(fill[static_cast<std::size_t>(graph.receiver[cell])]++)]
                = static_cast<int32_t>(cell);

    graph.order.reserve(count);
    for (std::size_t root = 0; root < count; root++)
    {
        if (graph.land[root] == 0 || graph.receiver[root] != static_cast<int32_t>(root))
            continue;
        std::size_t head = graph.order.size();
        graph.order.push_back(static_cast<int32_t>(root));
        for (; head < graph.order.size(); head++)
        {
            const std::size_t cell = static_cast<std::size_t>(graph.order[head]);
            for (int32_t k = donorstart[cell]; k < donorstart[cell + 1]; k++)
                graph.order.push_back(donors[static_cast<std::size_t>(k)]);
        }
    }
    // A cycle among the receivers has no root, so its cells would be missing.
    std::size_t landcells = 0;
    for (std::size_t cell = 0; cell < count; cell++)
        landcells += graph.land[cell] != 0 ? 1 : 0;
    if (graph.order.size() != landcells)
        throw std::logic_error("buildflowgraph: the receivers form a cycle (the order misses land cells)");
    return graph;
}

void accumulateflow(const SphereMetric& metric, const std::vector<float>& runoff_m_yr, FlowGraph& graph)
{
    const std::size_t count = graph.size();
    if (metric.size() != count || runoff_m_yr.size() != count)
        throw std::invalid_argument("accumulateflow: input sizes differ from the graph");
    graph.area_m2.assign(count, 0.0);
    graph.discharge_m3_yr.assign(count, 0.0);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (graph.land[cell] == 0)
            continue;
        const double area = metric.area_m2(cell);
        graph.area_m2[cell] = area;
        graph.discharge_m3_yr[cell] = static_cast<double>(runoff_m_yr[cell]) * area;
    }
    for (auto it = graph.order.rbegin(); it != graph.order.rend(); ++it)
    {
        const std::size_t cell = static_cast<std::size_t>(*it);
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver == cell)
            continue;
        graph.area_m2[receiver] += graph.area_m2[cell];
        graph.discharge_m3_yr[receiver] += graph.discharge_m3_yr[cell];
    }
}
}
