#include "simulations/landscape/runoff.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace landscape
{
namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
}

double openwaterevaporation(double temperaturec, double perc, double zeroc)
{
    return perc * std::max(0.0, temperaturec - zeroc);
}

double budykorunoff(double precipitation, double potentialevaporation, double w)
{
    if (!(w > 1.0))
        throw std::invalid_argument("budykorunoff: w must exceed 1");
    if (!(precipitation > 0.0))
        return 0.0;
    if (!(potentialevaporation > 0.0))
        return precipitation;
    // runoff = (P^w + PET^w)^(1 / w) - PET, factored by the larger of the two so that
    // neither power overflows and the difference keeps its precision when PET >> P.
    double runoff;
    if (potentialevaporation >= precipitation)
    {
        const double ratio = precipitation / potentialevaporation;
        runoff = potentialevaporation * std::expm1(std::log1p(std::pow(ratio, w)) / w);
    }
    else
    {
        const double ratio = potentialevaporation / precipitation;
        runoff = precipitation * (std::exp(std::log1p(std::pow(ratio, w)) / w) - ratio);
    }
    return std::clamp(runoff, 0.0, precipitation);
}

double accumulateflowwithlosses(const SphereMetric& metric, const std::vector<float>& runoff_m_yr,
                                const std::vector<float>& elevation, const TransmissionLoss& loss, FlowGraph& graph)
{
    const std::size_t count = graph.size();
    if (metric.size() != count || runoff_m_yr.size() != count || elevation.size() != count ||
        loss.deficit_m_yr.size() != count || (!loss.openwater.empty() && loss.openwater.size() != count) ||
        (!loss.depressionid.empty() && loss.depressionid.size() != count))
        throw std::invalid_argument("accumulateflowwithlosses: input sizes differ from the graph");
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
    double totalloss = 0.0;
    for (auto it = graph.order.rbegin(); it != graph.order.rend(); ++it)
    {
        const std::size_t cell = static_cast<std::size_t>(*it);
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        // Every upstream cell comes earlier in the reversed order, so the flow is complete here.
        const double flow = graph.discharge_m3_yr[cell];
        const double deficit = static_cast<double>(loss.deficit_m_yr[cell]);
        const bool flooded = !loss.depressionid.empty() && loss.depressionid[cell] != 0;
        if (deficit > 0.0 && flow > 0.0 && !flooded)
        {
            double evaporation = 0.0;
            if (!loss.openwater.empty() && loss.openwater[cell] != 0)
                evaporation = deficit * metric.area_m2(cell);
            else
            {
                const double length = static_cast<double>(graph.distance_m[cell]);
                if (length > 0.0 && loss.widthcoefficient > 0.0)
                {
                    double width = loss.widthcoefficient * std::pow(flow / secondsperyear, loss.widthexponent);
                    const double downstream = graph.outletsea[cell] >= 0
                        ? static_cast<double>(loss.sealevel)
                        : static_cast<double>(elevation[receiver]);
                    const double slope = (static_cast<double>(elevation[cell]) - downstream) / length;
                    if (slope < loss.floodplainslope)
                        width *= loss.floodplainfactor;
                    evaporation = width * length * deficit;
                }
            }
            const double removed = std::min(flow, evaporation);
            graph.discharge_m3_yr[cell] = flow - removed;
            totalloss += removed;
        }
        if (receiver == cell)
            continue;
        graph.area_m2[receiver] += graph.area_m2[cell];
        graph.discharge_m3_yr[receiver] += graph.discharge_m3_yr[cell];
    }
    return totalloss;
}
}
