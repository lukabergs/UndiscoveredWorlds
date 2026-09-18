#pragma once

#include <algorithm>

namespace climatephysics::detail
{
// Lower-reservoir mean humidity proxy, not a resolved 2 m humidity profile.
// Capacity uses the same lower-layer partition as the prognostic moisture model.
inline float boundarySurfaceRelativeHumidity(float boundaryWaterMm,
    float saturationColumnWaterMm, float boundaryCapacityFraction)
{
    const float capacity = std::max(1.0e-8f,
        saturationColumnWaterMm * std::clamp(boundaryCapacityFraction, 0.0f, 1.0f));
    return std::clamp(std::max(0.0f, boundaryWaterMm) / capacity, 0.0f, 1.0f);
}
}
