#pragma once

#include "climate_energy.hpp"
#include "climate_seasons.hpp"

#include <algorithm>
#include <array>

namespace climateenergy::detail
{
inline float seasonalPermanentLandIceFraction(
    const std::array<double, CLIMATESEASONCOUNT>& profileTemperatureC,
    double localTemperatureOffsetC)
{
    // Diagnose permanent cover with the existing warm-season survival rule.
    // Elevation already enters the local temperature offset.
    const double warmest = *std::max_element(profileTemperatureC.begin(), profileTemperatureC.end());
    return permanentLandIceFraction(static_cast<float>(warmest + localTemperatureOffsetC));
}
}
