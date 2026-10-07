#include "simulations/landscape/glacial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

// The glacial kernels of erosion plan step G1 (docs/architecture/erosion-hydrology-plan.md,
// section 6, "Decided for Phase G").
namespace landscape
{
namespace
{
constexpr double squaremetrespersquarekilometre = 1.0e6;
// Ohmura, Kasser and Funk 1992: P = a + b T + c T^2 (P mm/yr, T C).
constexpr double ohmuraconstant = 645.0;
constexpr double ohmuralinear = 296.0;
constexpr double ohmuraquadratic = 9.0;

void checksize(std::size_t actual, std::size_t expected, const char* function, const char* name)
{
    if (actual != expected)
        throw std::invalid_argument(std::string(function) + ": " + name + " does not have metric.size() elements");
}
}

double equilibriumlinetemperature(double precipitation_mm_yr)
{
    // c T^2 + b T + (a - P) = 0 on its increasing branch, T = (-b + sqrt(b^2 - 4 c (a - P))) / (2 c),
    // written as 2 (P - a) / (b + sqrt(...)) so that it does not cancel near P = a.
    const double precipitation = std::max(0.0, precipitation_mm_yr);
    const double discriminant =
        ohmuralinear * ohmuralinear - 4.0 * ohmuraquadratic * (ohmuraconstant - precipitation);
    return 2.0 * (precipitation - ohmuraconstant) / (ohmuralinear + std::sqrt(discriminant));
}

void equilibriumline(const SphereMetric& metric, const std::vector<float>& elevation,
                     const std::vector<float>& warmest_c, const std::vector<float>& precipitation_mm_yr,
                     double lapse_c_per_m, double cooling_c, std::vector<float>& ela_m)
{
    const std::size_t cells = metric.size();
    checksize(elevation.size(), cells, "equilibriumline", "elevation");
    checksize(warmest_c.size(), cells, "equilibriumline", "warmest_c");
    checksize(precipitation_mm_yr.size(), cells, "equilibriumline", "precipitation_mm_yr");
    if (!(lapse_c_per_m > 0.0))
        throw std::invalid_argument("equilibriumline: the lapse rate must be positive");
    ela_m.resize(cells);
    for (std::size_t i = 0; i < cells; i++)
    {
        const double temperature = equilibriumlinetemperature(static_cast<double>(precipitation_mm_yr[i]));
        const double excess = static_cast<double>(warmest_c[i]) - cooling_c - temperature;
        ela_m[i] = static_cast<float>(static_cast<double>(elevation[i]) + excess / lapse_c_per_m);
    }
}

void accumulateice(const SphereMetric& metric, const FlowGraph& graph, const std::vector<float>& surface,
                   const std::vector<float>& ela_m, const std::vector<float>& runoff_m_yr, double fullicerange_m,
                   std::vector<double>& iceflux_m3_yr)
{
    const std::size_t cells = metric.size();
    checksize(graph.size(), cells, "accumulateice", "graph.receiver");
    checksize(graph.land.size(), cells, "accumulateice", "graph.land");
    checksize(surface.size(), cells, "accumulateice", "surface");
    checksize(ela_m.size(), cells, "accumulateice", "ela_m");
    checksize(runoff_m_yr.size(), cells, "accumulateice", "runoff_m_yr");
    if (!(fullicerange_m > 0.0))
        throw std::invalid_argument("accumulateice: the full-ice range must be positive");
    for (int32_t entry : graph.order)
    {
        if (entry < 0 || static_cast<std::size_t>(entry) >= cells)
            throw std::invalid_argument("accumulateice: graph.order entry is not a cell index");
        const int32_t receiver = graph.receiver[static_cast<std::size_t>(entry)];
        if (receiver < 0 || static_cast<std::size_t>(receiver) >= cells)
            throw std::invalid_argument("accumulateice: graph.receiver is not a cell index");
    }

    // Each cell of the order starts with its source; over the order reversed every
    // donor has added its flux before the cell is clipped at zero and passes it on.
    iceflux_m3_yr.assign(cells, 0.0);
    for (int32_t entry : graph.order)
    {
        const std::size_t cell = static_cast<std::size_t>(entry);
        const double fraction = std::min(
            (static_cast<double>(surface[cell]) - static_cast<double>(ela_m[cell])) / fullicerange_m, 1.0);
        iceflux_m3_yr[cell] = static_cast<double>(runoff_m_yr[cell]) * metric.area_m2(cell) * fraction;
    }
    for (auto it = graph.order.rbegin(); it != graph.order.rend(); ++it)
    {
        const std::size_t cell = static_cast<std::size_t>(*it);
        iceflux_m3_yr[cell] = std::max(0.0, iceflux_m3_yr[cell]);
        const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
        if (receiver == cell)
            continue;  // an outlet's ice calves; a root keeps it
        iceflux_m3_yr[receiver] += iceflux_m3_yr[cell];
    }
}

void icethickness(const std::vector<double>& iceflux_m3_yr, const GlacialParameters& parameters,
                  std::vector<float>& thickness_m)
{
    thickness_m.resize(iceflux_m3_yr.size());
    for (std::size_t i = 0; i < iceflux_m3_yr.size(); i++)
    {
        const double flux = iceflux_m3_yr[i];
        if (!(flux > 0.0))
        {
            thickness_m[i] = 0.0f;
            continue;
        }
        // A_i = flux / (1 m/yr) in m^2, in km^2 for the width law.
        const double area_km2 = flux / squaremetrespersquarekilometre;
        const double width = parameters.widthreference_m * std::pow(area_km2, parameters.widthexponent);
        thickness_m[i] = static_cast<float>(std::min(parameters.maximumthickness_m, parameters.thicknessratio * width));
    }
}
}
