#pragma once

#include <algorithm>
#include <cmath>

namespace climateatmosphere::detail
{
// A bulk two-reservoir approximation, not a resolved humidity profile. The
// diagnosed 850-hPa velocity may carry free-reservoir moisture only where the
// pressure level is above terrain. Missing/buried levels retain the 500 mode.
inline float lowLevelMoistureWeight(float fraction, float availability,
    float terrainMetres, float surfaceTemperatureK, float gravityMps2)
{
    if (!std::isfinite(availability) || !std::isfinite(terrainMetres) ||
        !(surfaceTemperatureK > 0.0f) || !(gravityMps2 > 0.0f)) return 0.0f;
    const float height850 = 287.05f * surfaceTemperatureK / gravityMps2 * std::log(100000.0f / 85000.0f);
    if (terrainMetres >= height850) return 0.0f;
    return std::clamp(fraction, 0.0f, 1.0f) * std::clamp(availability, 0.0f, 1.0f);
}

inline float moistureCarrierMode(float midLevelWind, float lowLevelWind, float weight)
{
    if (!(weight > 0.0f) || !std::isfinite(lowLevelWind)) return midLevelWind;
    return midLevelWind + std::clamp(weight, 0.0f, 1.0f) * (lowLevelWind - midLevelWind);
}
}
