#pragma once

#include <algorithm>
#include <cmath>

namespace climatehydrology::detail
{
// Convert a one-day removal probability to another duration. For a fixed
// coefficient, successive substeps retain the same water as one whole day.
inline float dailyFractionForTimeStep(float dailyFraction, float seconds)
{
    if (!(seconds > 0.0f)) return 0.0f;
    const float fraction = std::clamp(dailyFraction, 0.0f, 1.0f);
    if (seconds == 86400.0f || fraction == 1.0f) return fraction;
    return -std::expm1(std::log1p(-fraction) * (seconds / 86400.0f));
}
}
