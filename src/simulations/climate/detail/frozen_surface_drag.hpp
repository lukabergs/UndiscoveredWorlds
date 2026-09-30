#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
// Minimum seasonal mean cover is a roughness proxy, not an ice-mass diagnosis.
// Form the minimum before remapping so snow in different subcells/seasons
// cannot manufacture persistent cover. Missing quarters provide no evidence.
template<std::size_t Seasons>
inline std::vector<float> persistentSnowCover(
    const std::array<const std::vector<float>*, Seasons>& quarters, std::size_t cells)
{
    std::vector<float> result(cells, 0.0f);
    for (const auto* quarter : quarters)
        if (!quarter || quarter->size() != cells) return result;
    std::fill(result.begin(), result.end(), 1.0f);
    for (const auto* quarter : quarters)
        for (std::size_t cell = 0; cell < cells; ++cell)
            result[cell] = std::min(result[cell], std::isfinite((*quarter)[cell]) ?
                std::clamp((*quarter)[cell], 0.0f, 1.0f) : 0.0f);
    return result;
}

// Frozen cover is a fraction of the whole grid cell, bounded by land fraction.
// Preserve the legacy drag exactly when cover is absent; blend stress
// coefficients, not velocities. Fully frozen land approaches its own roughness
// rather than retaining an absolute-elevation penalty intended for rough land.
inline float frozenSurfaceDragCoefficient(float landFraction, float elevationMetres,
    float frozenLandFraction, float oceanCd, float landCd, float highReliefCd, float frozenCd)
{
    const float land = std::clamp(landFraction, 0.0f, 1.0f);
    const float base = oceanCd + land * (landCd - oceanCd);
    const float relief = std::clamp(elevationMetres / 4500.0f, 0.0f, 1.0f);
    const float legacy = base + relief * (highReliefCd - base);
    if (!(land > 0.0f) || !(frozenLandFraction > 0.0f)) return legacy;
    const float fraction = std::clamp(frozenLandFraction / land, 0.0f, 1.0f);
    const float frozen = oceanCd + land * (frozenCd - oceanCd);
    return (1.0f - fraction) * legacy + fraction * frozen;
}
}
