#pragma once

#include <algorithm>

namespace climatehydrology::detail
{
// Mean condensable water for a uniform subcell distribution with the same
// cell-mean water. The support is nonnegative; no water is introduced.
inline float subgridCondensableWater(float meanWater, float capacity, float fractionalHalfWidth)
{
    const float water = std::max(0.0f, meanWater);
    const float threshold = std::max(0.0f, capacity);
    const float halfWidth = water * std::clamp(fractionalHalfWidth, 0.0f, 1.0f);
    if (halfWidth == 0.0f || water - halfWidth >= threshold)
        return std::max(0.0f, water - threshold);
    if (water + halfWidth <= threshold) return 0.0f;
    const double wetWidth = static_cast<double>(water) + halfWidth - threshold;
    return static_cast<float>(std::clamp(wetWidth * wetWidth / (4.0 * halfWidth), 0.0, static_cast<double>(water)));
}
}
