#include "simulations/landscape/ocean_mask.hpp"

#include <cstddef>
#include <stdexcept>

namespace landscape
{
std::vector<uint8_t> oceanmask(const SphereMetric& metric, const std::vector<uint8_t>& belowsea,
                               double minimumareafraction)
{
    const std::size_t count = metric.size();
    if (belowsea.size() != count)
        throw std::invalid_argument("oceanmask: the mask's size differs from the metric");

    // Label the 8-connected below-sea regions in cell order and sum their areas; the
    // largest by area is found first on a tie.
    std::vector<int32_t> region(count, -1);
    std::vector<double> areas;
    std::vector<std::size_t> members;
    int32_t largest = -1;
    for (std::size_t start = 0; start < count; start++)
    {
        if (belowsea[start] == 0 || region[start] >= 0)
            continue;
        const int32_t id = static_cast<int32_t>(areas.size());
        region[start] = id;
        members.assign(1, start);
        double area = 0.0;
        for (std::size_t next = 0; next < members.size(); next++)
        {
            const std::size_t cell = members[next];
            area += metric.area_m2(cell);
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::size_t neighbour = metric.neighbour(cell, direction);
                if (neighbour == cell || belowsea[neighbour] == 0 || region[neighbour] >= 0)
                    continue;
                region[neighbour] = id;
                members.push_back(neighbour);
            }
        }
        areas.push_back(area);
        if (largest < 0 || area > areas[static_cast<std::size_t>(largest)])
            largest = id;
    }

    const double minimumarea = minimumareafraction * metric.totalarea_m2();
    std::vector<uint8_t> sea(areas.size(), 0);
    for (std::size_t id = 0; id < areas.size(); id++)
        sea[id] = static_cast<int32_t>(id) == largest || areas[id] >= minimumarea ? 1 : 0;
    std::vector<uint8_t> ocean(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        ocean[cell] = region[cell] >= 0 ? sea[static_cast<std::size_t>(region[cell])] : 0;
    return ocean;
}
}
