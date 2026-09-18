#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
// Calibrated projection of surface temperature gradients onto the upper mode.
// Scale neighbouring differences, then integrate and remove the area mean.
// Multiplying height anomalies by a latitude envelope instead would introduce
// an extra gradient where that envelope changes, potentially reversing winds.
inline std::vector<float> projectUpperThermalGradient(
    const std::vector<float>& temperature, const std::vector<double>& rowAreas,
    float coreFraction, float outerFraction, float coreLatitude, float outerLatitude,
    float polarFraction = 1.0f, float polarStartLatitude = 60.0f, float polarEndLatitude = 80.0f)
{
    if (temperature.empty() || temperature.size() != rowAreas.size()) return {};
    coreFraction = std::clamp(coreFraction, 0.0f, 1.0f);
    outerFraction = std::clamp(outerFraction, 0.0f, 1.0f);
    polarFraction = std::clamp(polarFraction, 0.0f, 1.0f);
    const double inner = std::max(0.0f, coreLatitude);
    const double width = std::max(0.001, static_cast<double>(outerLatitude) - inner);
    const double polarInner = std::max(0.0f, polarStartLatitude);
    const double polarWidth = std::max(0.001, static_cast<double>(polarEndLatitude) - polarInner);
    std::vector<double> integrated(temperature.size(), 0.0);
    double area = 0.0, weighted = 0.0;
    for (std::size_t y = 0; y < temperature.size(); ++y)
    {
        if (!(rowAreas[y] > 0.0) || !std::isfinite(temperature[y])) return {};
        if (y > 0)
        {
            const double latitude = 90.0 - 180.0 * y / temperature.size();
            const double t = std::clamp((std::abs(latitude) - inner) / width, 0.0, 1.0);
            const double polar = std::clamp((std::abs(latitude) - polarInner) / polarWidth, 0.0, 1.0);
            const double fraction = (coreFraction + (outerFraction - coreFraction) * t * t * (3.0 - 2.0 * t)) *
                (1.0 + (polarFraction - 1.0) * polar * polar * (3.0 - 2.0 * polar));
            integrated[y] = integrated[y - 1] + fraction * (temperature[y] - temperature[y - 1]);
        }
        area += rowAreas[y];
        weighted += rowAreas[y] * integrated[y];
    }
    std::vector<float> result(temperature.size());
    for (std::size_t y = 0; y < temperature.size(); ++y)
        result[y] = static_cast<float>(integrated[y] - weighted / area);
    return result;
}
}
