#include "simulations/landscape/subgrid_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace landscape
{
double characteristiclength_km(const SubgridParameters& parameters, double relativeerodibility, double runoff_mm_yr)
{
    if (!(relativeerodibility > 0.0) || !(runoff_mm_yr > 0.0) || !(parameters.referenceerodibility > 0.0) ||
        !(parameters.referencerunoff_mm_yr > 0.0) || !(parameters.referencedoverk > 0.0))
        return 0.0;
    // K_A(reference) / K_A(cell) with K_A = K x e x r^m_ref (Q = r A in the site's K Q^m_ref).
    const double ratio = (parameters.referenceerodibility / relativeerodibility) *
                         std::pow(parameters.referencerunoff_mm_yr / runoff_mm_yr, parameters.referencem);
    const double length_m = std::pow(parameters.referencedoverk * ratio, 1.0 / (2.0 * parameters.referencem + 1.0));
    return length_m / 1000.0;
}

double hillsloperelief_m(const SubgridParameters& parameters, double erosionrate_m_myr, double halflength_km)
{
    const double length = std::max(0.0, halflength_km);
    if (!(parameters.diffusivity_km2_myr > 0.0))
        return 0.0;
    // m/Myr x km^2 / (km^2/Myr) = m
    return std::max(0.0, erosionrate_m_myr) * length * length / (2.0 * parameters.diffusivity_km2_myr);
}

std::size_t subgridstatistics(const SphereMetric& metric, const SubgridParameters& parameters,
                              const std::vector<float>& elevation_m, const std::vector<std::uint8_t>& sea,
                              const std::vector<float>& erodibility, const std::vector<float>& runoff_mm_yr,
                              const std::vector<float>& erosionrate_m_myr, std::vector<float>& relief_m,
                              std::vector<float>& spacing_km)
{
    const std::size_t size = metric.size();
    if (elevation_m.size() != size || sea.size() != size || erodibility.size() != size ||
        runoff_mm_yr.size() != size || erosionrate_m_myr.size() != size)
        throw std::invalid_argument("subgridstatistics: an input differs in size from the grid");
    relief_m.assign(size, 0.0f);
    spacing_km.assign(size, 0.0f);
    std::size_t capped = 0;
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (sea[cell] != 0)
            continue;
        const double cellsize_km = std::sqrt(metric.area_m2(cell)) / 1000.0;
        double spacing = parameters.spacingcoefficient *
                         characteristiclength_km(parameters, erodibility[cell], runoff_mm_yr[cell]);
        if (spacing > cellsize_km)  // no sub-grid valleys: stored as the cell size
        {
            spacing = cellsize_km;
            capped++;
        }
        float local = 0.0f;
        for (int direction = 1; direction <= directioncount; direction++)
        {
            const std::int32_t other = metric.neighbourid(cell, direction);
            if (other >= 0)
                local = std::max(local, std::fabs(elevation_m[static_cast<std::size_t>(other)] - elevation_m[cell]));
        }
        const double relief = hillsloperelief_m(parameters, erosionrate_m_myr[cell], 0.5 * cellsize_km);
        relief_m[cell] = static_cast<float>(std::min(relief, static_cast<double>(local)));
        spacing_km[cell] = static_cast<float>(spacing);
    }
    return capped;
}
}
