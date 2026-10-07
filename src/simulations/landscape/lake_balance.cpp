#include "simulations/landscape/lake_balance.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>

namespace landscape
{
namespace
{
void checksize(std::size_t actual, std::size_t expected, const char* name)
{
    if (actual != expected)
        throw std::invalid_argument(std::string("balancelakes: ") + name + " does not match the metric's size");
}

// A cell downstream of a depression belongs to its flat while its filled level is within
// this of the depression's spill (float rounding of equal filled levels).
constexpr float flatmargin_m = 1.0e-3f;

struct DepressionTotals
{
    std::vector<std::size_t> cells;  // ascending index
    double area_m2 = 0.0;
    double precipitation_m3_yr = 0.0;  // sum of precipitation x area
    double evaporation_m3_yr = 0.0;    // sum of evaporation x area
    double inflow_m3_yr = 0.0;
    float spill_m = 0.0f;
};
}

LakeBalance balancelakes(const LakeBalanceInputs& inputs)
{
    if (inputs.metric == nullptr)
        throw std::invalid_argument("balancelakes: no metric");
    const SphereMetric& metric = *inputs.metric;
    const std::size_t size = metric.size();
    checksize(inputs.elevation.size(), size, "elevation");
    checksize(inputs.filledelevation.size(), size, "filledelevation");
    checksize(inputs.depressionid.size(), size, "depressionid");
    checksize(inputs.receiver.size(), size, "receiver");
    checksize(inputs.discharge_m3_yr.size(), size, "discharge_m3_yr");
    checksize(inputs.runoff_m_yr.size(), size, "runoff_m_yr");
    checksize(inputs.precipitation_m_yr.size(), size, "precipitation_m_yr");
    checksize(inputs.evaporation_m_yr.size(), size, "evaporation_m_yr");

    // Ascending id; cells and sums in ascending cell index, so the result is deterministic.
    std::map<int32_t, DepressionTotals> totals;
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int32_t id = inputs.depressionid[cell];
        if (id == 0)
            continue;
        DepressionTotals& depression = totals[id];
        const double area = metric.area_m2(cell);
        if (depression.cells.empty() || inputs.filledelevation[cell] > depression.spill_m)
            depression.spill_m = inputs.filledelevation[cell];
        depression.cells.push_back(cell);
        depression.area_m2 += area;
        depression.precipitation_m3_yr += static_cast<double>(inputs.precipitation_m_yr[cell]) * area;
        depression.evaporation_m3_yr += static_cast<double>(inputs.evaporation_m_yr[cell]) * area;
        depression.inflow_m3_yr += static_cast<double>(inputs.runoff_m_yr[cell]) * area;
    }
    // Cells of d's flat outside d: shallow flooded cells (raised by less than the routing's
    // minimum depth, so id 0) that lie downstream of d's cells while still at d's level. Their
    // discharge already holds d's own water, so they are not external donors (B2 review).
    std::vector<int32_t> flat(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int32_t id = inputs.depressionid[cell];
        if (id == 0)
            continue;
        const float level = totals[id].spill_m - flatmargin_m;
        int32_t next = inputs.receiver[cell];
        for (std::size_t step = 0; step < size && next >= 0 && static_cast<std::size_t>(next) < size; step++)
        {
            const std::size_t at = static_cast<std::size_t>(next);
            if (inputs.depressionid[at] == id || flat[at] == id || inputs.filledelevation[at] < level)
                break;
            flat[at] = id;
            next = inputs.receiver[at];
        }
    }
    // External donors: a cell outside d and outside d's flat whose receiver is a cell of d
    // brings its discharge.
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int32_t receiver = inputs.receiver[cell];
        if (receiver < 0 || static_cast<std::size_t>(receiver) >= size)
            continue;
        const int32_t target = inputs.depressionid[static_cast<std::size_t>(receiver)];
        if (target == 0 || target == inputs.depressionid[cell] || flat[cell] == target)
            continue;
        totals[target].inflow_m3_yr += inputs.discharge_m3_yr[cell];
    }

    LakeBalance balance;
    balance.surface.assign(size, 0.0f);
    balance.lakes.reserve(totals.size());
    for (auto& [id, depression] : totals)
    {
        std::vector<std::size_t>& cells = depression.cells;
        std::sort(cells.begin(), cells.end(), [&](std::size_t a, std::size_t b) {
            const float za = inputs.elevation[a];
            const float zb = inputs.elevation[b];
            return za != zb ? za < zb : a < b;
        });

        LakeLevel lake;
        lake.id = id;
        lake.spill_m = depression.spill_m;
        lake.lowestcell = static_cast<int32_t>(cells.front());
        lake.depressionarea_m2 = depression.area_m2;
        lake.inflow_m3_yr = depression.inflow_m3_yr;
        const double precipitation = depression.precipitation_m3_yr / depression.area_m2;
        const double evaporation = depression.evaporation_m3_yr / depression.area_m2;
        lake.netloss_m_yr = std::max(0.0, evaporation - precipitation);

        const bool overflowing = lake.netloss_m_yr <= 0.0 ||
                                 lake.inflow_m3_yr / lake.netloss_m_yr >= depression.area_m2;
        if (overflowing)
        {
            lake.overflowing = true;
            lake.level_m = lake.spill_m;
            lake.lakecells = static_cast<int32_t>(cells.size());
            lake.area_m2 = depression.area_m2;
        }
        else
        {
            const double equilibrium = lake.inflow_m3_yr / lake.netloss_m_yr;
            std::size_t k = 0;
            double area = 0.0;
            if (equilibrium >= metric.area_m2(cells.front()))
            {
                while (k < cells.size() && area < equilibrium)
                    area += metric.area_m2(cells[k++]);
                // A lake level cannot separate cells of equal terrain: the whole tie group of
                // the k-th lowest cell is lake (B2 review).
                while (k < cells.size() && inputs.elevation[cells[k]] == inputs.elevation[cells[k - 1]])
                    area += metric.area_m2(cells[k++]);
            }
            lake.lakecells = static_cast<int32_t>(k);
            lake.area_m2 = area;
            if (k > 0)
            {
                // The spill stands for the (k+1)-th level when every cell is lake but the
                // balance still closes the lake.
                const float upper = k < cells.size() ? inputs.elevation[cells[k]] : lake.spill_m;
                const float midpoint = 0.5f * (inputs.elevation[cells[k - 1]] + upper);
                lake.level_m = std::min(midpoint, lake.spill_m);
            }
        }
        for (int32_t index = 0; index < lake.lakecells; index++)
            balance.surface[cells[static_cast<std::size_t>(index)]] = lake.level_m;
        balance.lakes.push_back(lake);
    }
    return balance;
}

std::vector<uint8_t> saltpans(const LakeBalanceInputs& inputs, const LakeBalance& balance,
                              const std::vector<uint8_t>& belowsea)
{
    if (inputs.metric == nullptr)
        throw std::invalid_argument("saltpans: no metric");
    const std::size_t size = inputs.metric->size();
    if (inputs.depressionid.size() != size || balance.surface.size() != size || belowsea.size() != size)
        throw std::invalid_argument("saltpans: a size does not match the metric's");
    std::vector<int32_t> closedids;  // ascending, as balance.lakes
    std::vector<uint8_t> pan(size, 0);
    for (const LakeLevel& lake : balance.lakes)
    {
        if (lake.overflowing)
            continue;
        closedids.push_back(lake.id);
        if (lake.lakecells == 0 && lake.lowestcell >= 0)
            pan[static_cast<std::size_t>(lake.lowestcell)] = 1;
    }
    for (std::size_t cell = 0; cell < size && !closedids.empty(); cell++)
    {
        const int32_t id = inputs.depressionid[cell];
        if (id != 0 && belowsea[cell] != 0 && balance.surface[cell] == 0.0f &&
            std::binary_search(closedids.begin(), closedids.end(), id))
            pan[cell] = 1;
    }
    return pan;
}
}
