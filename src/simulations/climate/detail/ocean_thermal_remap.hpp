#pragma once

#include "climate_grid.hpp"

#include <cstdint>
#include <vector>

namespace climateocean::detail
{
// Reconstruct an intensive ocean field using only wet donors. The smooth
// numerator and wet support retain the grid remap weights; their quotient
// preserves wet constants and bounds, not an integral over a binary target mask.
inline std::vector<float> remapOceanThermalField(
    const climategrid::ConservativeRemap& plan,
    const std::vector<std::uint8_t>& landMask,
    const std::vector<float>& source,
    const std::vector<float>& destinationFallback)
{
    const auto sourceCount = static_cast<std::size_t>(plan.sourceColumns) * plan.sourceRows;
    const auto columns = plan.longitude.size(), rows = plan.latitude.size();
    if (plan.sourceColumns <= 0 || plan.sourceRows <= 0 || plan.reconstruct ||
        source.size() != sourceCount || landMask.size() != sourceCount ||
        destinationFallback.size() != columns * rows)
        return {};

    auto result = destinationFallback;
    for (std::size_t y = 0; y < rows; ++y)
        for (std::size_t x = 0; x < columns; ++x)
        {
            double weightedValue = 0.0, wetSupport = 0.0;
            for (const auto& yw : plan.latitude[y])
                for (const auto& xw : plan.longitude[x])
                {
                    const auto donor = static_cast<std::size_t>(yw.source) * plan.sourceColumns + xw.source;
                    if (landMask[donor]) continue;
                    const double weight = yw.fraction * xw.fraction;
                    weightedValue += weight * source[donor];
                    wetSupport += weight;
                }
            if (wetSupport > 0.0)
                result[y * columns + x] = static_cast<float>(weightedValue / wetSupport);
        }
    return result;
}
}
