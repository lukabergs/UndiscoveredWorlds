#include "simulations/landscape/flow_graph.hpp"

#include <stdexcept>

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
        const std::size_t next = metric.neighbour(cell, direction);
        if (next != cell && next == target)
            return direction;
    }
    return 0;
}
}

FlowGraph buildflowgraph(const SphereMetric& metric, const std::vector<float>& elevation,
                         const std::vector<uint8_t>& sea, float sealevel, const Depressions& depressions)
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
                const std::size_t next = metric.neighbour(cell, direction);
                if (next == cell || sea[next] == 0 || !(elevation[next] < lowest))
                    continue;
                lowest = elevation[next];
                best = direction;
            }
            if (best != 0)
            {
                graph.receiver[cell] = static_cast<int32_t>(metric.neighbour(cell, best));
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
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell)
                continue;
            const double height = sea[next] != 0 ? static_cast<double>(sealevel) : static_cast<double>(filled[next]);
            const double gradient = (static_cast<double>(filled[cell]) - height) / metric.distance_m(cell, direction);
            if (gradient > steepest)
            {
                steepest = gradient;
                best = direction;
            }
        }

        std::size_t target = cell;
        if (best != 0 && (sea[metric.neighbour(cell, best)] != 0 || depressions.depressionid[cell] == 0))
            target = metric.neighbour(cell, best);
        else if (depressions.floodreceiver[cell] >= 0)
        {
            target = static_cast<std::size_t>(depressions.floodreceiver[cell]);
            best = directionto(metric, cell, target);
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
