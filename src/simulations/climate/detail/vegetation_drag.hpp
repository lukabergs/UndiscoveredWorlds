#pragma once
#include "frozen_surface_drag.hpp"

namespace climateatmosphere::detail
{
inline float vegetationSurfaceDragCoefficient(float landFraction, float elevationMetres,
    float frozenLandFraction, float oceanCd, float landCd, float highReliefCd, float frozenCd,
    float dragContrast, float strength)
{
    const float original = frozenSurfaceDragCoefficient(landFraction, elevationMetres,
        frozenLandFraction, oceanCd, landCd, highReliefCd, frozenCd);
    if (!(strength > 0.0f)) return original;
    const float exposedLand = std::max(0.0f, std::clamp(landFraction, 0.0f, 1.0f) -
        std::clamp(frozenLandFraction, 0.0f, 1.0f));
    // Replace the land component's spatial contrast, retaining relief and ice
    // contributions. Fractions weight stress coefficients, never speeds.
    // Preserve the original calibration path exactly. Above unit strength,
    // amplify protection but keep low-roughness reductions at their input value.
    // This widens the contrast without forcing all exposed terrain to a floor.
    if (strength > 1.0f)
    {
        const float contrast = std::clamp(dragContrast, 0.25f, 4.0f) - 1.0f;
        const float gain = contrast > 0.0f ? std::min(strength, 16.0f) : 1.0f;
        const float excess = gain * contrast;
        return original + exposedLand * landCd * excess;
    }
    return original + exposedLand * landCd * std::clamp(strength, 0.0f, 1.0f) *
        (std::clamp(dragContrast, 0.25f, 4.0f) - 1.0f);
}
}
