#include "simulations/landscape/depressions.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <stdexcept>
#include <utility>

namespace landscape
{
int32_t terminalsink(const std::vector<float>& elevation, const std::vector<uint8_t>& sea)
{
    int32_t sink = -1;
    for (std::size_t cell = 0; cell < elevation.size(); cell++)
    {
        if (cell < sea.size() && sea[cell] != 0)
            return -1;
        if (sink < 0 || elevation[cell] < elevation[static_cast<std::size_t>(sink)])
            sink = static_cast<int32_t>(cell);
    }
    return sink;
}

Depressions filldepressions(const SphereMetric& metric, const std::vector<float>& elevation,
                            const std::vector<uint8_t>& sea, float minimumdepth_m)
{
    const std::size_t count = metric.size();
    if (elevation.size() != count || sea.size() != count)
        throw std::invalid_argument("filldepressions: input sizes differ from the metric");

    Depressions result;
    result.filledelevation = elevation;
    result.depressionid.assign(count, 0);
    result.floodreceiver.assign(count, -1);

    // Priority-Flood with flow directions (Barnes et al. 2014): a min-queue on
    // (filled elevation, cell index), so ties go to the lower index.
    using Entry = std::pair<float, int32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    std::vector<uint8_t> closed(count, 0);
    std::vector<int32_t> closeorder;  // land cells in the order they were closed
    closeorder.reserve(count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (sea[cell] == 0)
            continue;
        closed[cell] = 1;
        queue.push({ elevation[cell], static_cast<int32_t>(cell) });
    }
    // No sea cell: the terminal sink is the one seed (a virtual outlet).
    const int32_t sink = terminalsink(elevation, sea);
    if (sink >= 0)
    {
        closed[static_cast<std::size_t>(sink)] = 1;
        queue.push({ elevation[static_cast<std::size_t>(sink)], sink });
    }
    bool raised = false;   // the flood has raised a cell
    float passlevel = 0.0f;  // the filled level at the first raise: the pass out of the sink's basin
    // A raised cell goes to a first-in first-out pit queue that is emptied before
    // the next priority pop (Barnes et al. 2014, Algorithm 2): a depression is then
    // flooded entirely from the cell the flood first entered it from, so every
    // flood-receiver chain in it passes through its spill even when several rim
    // cells share the spill elevation. The filled surface is the same.
    std::vector<int32_t> pit;
    std::size_t pithead = 0;
    while (pithead < pit.size() || !queue.empty())
    {
        std::size_t cell;
        if (pithead < pit.size())
            cell = static_cast<std::size_t>(pit[pithead++]);
        else
        {
            cell = static_cast<std::size_t>(queue.top().second);
            queue.pop();
            pit.clear();
            pithead = 0;
        }
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || closed[next] != 0)
                continue;
            closed[next] = 1;
            closeorder.push_back(static_cast<int32_t>(next));
            result.floodreceiver[next] = static_cast<int32_t>(cell);
            if (elevation[next] < result.filledelevation[cell])
            {
                if (!raised)
                {
                    raised = true;
                    passlevel = result.filledelevation[cell];
                }
                result.filledelevation[next] = result.filledelevation[cell];
                pit.push_back(static_cast<int32_t>(next));
            }
            else
                queue.push({ elevation[next], static_cast<int32_t>(next) });
        }
    }

    // Depressions are the connected regions (eight neighbours) of cells raised by
    // more than minimumdepth, numbered in the order the flood first entered them;
    // the spill is the unraised rim cell that first entry came from. A connected raised region
    // has one filled level, so the spill elevation is that of every cell of it.
    const auto deep = [&](std::size_t cell) {
        return static_cast<double>(result.filledelevation[cell]) - static_cast<double>(elevation[cell])
               > static_cast<double>(minimumdepth_m);
    };
    std::vector<int32_t> stack;
    // The terminal basin (grid without a sea cell): the cells below the pass level
    // connected to the sink, all closed before the first raise (so at their own
    // terrain, with flood receivers leading to the sink), raised to the pass level.
    if (sink >= 0)
    {
        float level = passlevel;
        if (!raised)
            for (float height : elevation)
                level = std::max(level, height);
        if (elevation[static_cast<std::size_t>(sink)] < level)
        {
            Depression depression;
            depression.id = 1;
            depression.spillcell = -1;
            depression.spillelevation_m = level;
            depression.terminal = true;
            result.depressionid[static_cast<std::size_t>(sink)] = depression.id;
            stack.assign(1, sink);
            while (!stack.empty())
            {
                const std::size_t cell = static_cast<std::size_t>(stack.back());
                stack.pop_back();
                result.filledelevation[cell] = level;
                for (int direction = 1; direction <= directioncount; direction++)
                {
                    const std::size_t next = metric.neighbour(cell, direction);
                    if (next == cell || result.depressionid[next] != 0 || !(elevation[next] < level))
                        continue;
                    result.depressionid[next] = depression.id;
                    stack.push_back(static_cast<int32_t>(next));
                }
            }
            result.depressions.push_back(depression);
        }
    }
    for (int32_t entry : closeorder)
    {
        const std::size_t first = static_cast<std::size_t>(entry);
        if (result.depressionid[first] != 0 || !deep(first))
            continue;
        Depression depression;
        depression.id = static_cast<int32_t>(result.depressions.size()) + 1;
        // The spill is the rim: the first unraised cell on the flood-receiver chain
        // (a margin raised by less than minimumdepth lies inside the fill).
        int32_t spill = result.floodreceiver[first];
        while (spill >= 0
               && result.filledelevation[static_cast<std::size_t>(spill)] > elevation[static_cast<std::size_t>(spill)])
            spill = result.floodreceiver[static_cast<std::size_t>(spill)];
        depression.spillcell = spill;
        depression.spillelevation_m = result.filledelevation[first];
        result.depressionid[first] = depression.id;
        stack.assign(1, entry);
        while (!stack.empty())
        {
            const std::size_t cell = static_cast<std::size_t>(stack.back());
            stack.pop_back();
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::size_t next = metric.neighbour(cell, direction);
                if (next == cell || sea[next] != 0 || result.depressionid[next] != 0 || !deep(next))
                    continue;
                result.depressionid[next] = depression.id;
                stack.push_back(static_cast<int32_t>(next));
            }
        }
        result.depressions.push_back(depression);
    }

    for (std::size_t cell = 0; cell < count; cell++)
    {
        const int32_t id = result.depressionid[cell];
        if (id == 0)
            continue;
        Depression& depression = result.depressions[static_cast<std::size_t>(id - 1)];
        const double area = metric.area_m2(cell);
        depression.area_m2 += area;
        depression.volume_m3
            += (static_cast<double>(result.filledelevation[cell]) - static_cast<double>(elevation[cell])) * area;
    }
    return result;
}
}
