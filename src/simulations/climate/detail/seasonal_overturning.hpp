#pragma once
#include <algorithm>
#include <cmath>
#include "../climate_atmosphere.hpp"

namespace climateatmosphere::detail
{
// The winter overturning cell spans the geographical equator to reach the
// summer thermal trough. Optionally anchor its subtropical edge instead of
// translating both entire hemispheric pressure profiles with that trough.
inline float winterAnchoredHadleyWidth(float latitudeDegrees, float troughDegrees,
    float baseWidthDegrees, float anchorFraction)
{
    const bool winterBranch = troughDegrees > 0.0f ? latitudeDegrees < troughDegrees :
        troughDegrees < 0.0f && latitudeDegrees > troughDegrees;
    return baseWidthDegrees + (winterBranch ? std::clamp(anchorFraction, 0.0f, 1.0f) *
        std::abs(troughDegrees) : 0.0f);
}

// Move only the tropical minimum. The old subtropical maxima and every
// extratropical value remain fixed (up to the caller's pressure gauge).
// Piecewise affine coordinates retain zero pressure slope at each extremum.
inline float confinedOverturningPressureHpa(float latitudeDegrees, float controlTrough,
    float baseWidthDegrees, float anchorFraction, float targetTrough, float amplitudeHpa)
{
    if (baseWidthDegrees <= 0.0f || amplitudeHpa <= 0.0f) return 0.0f;
    const auto pressure = [&](float latitude) {
        return axisymmetricOverturningPressureAnomalyHpa(latitude, controlTrough,
            winterAnchoredHadleyWidth(latitude, controlTrough, baseWidthDegrees, anchorFraction), amplitudeHpa);
    };
    const float northEdge = controlTrough + std::min(
        winterAnchoredHadleyWidth(90, controlTrough, baseWidthDegrees, anchorFraction),
        0.8f * (90.0f - controlTrough));
    const float southEdge = controlTrough - std::min(
        winterAnchoredHadleyWidth(-90, controlTrough, baseWidthDegrees, anchorFraction),
        0.8f * (90.0f + controlTrough));
    const float minimumBranch = std::min(2.0f, 0.1f * (northEdge - southEdge));
    targetTrough = std::clamp(targetTrough, southEdge + minimumBranch, northEdge - minimumBranch);
    if (targetTrough == controlTrough || latitudeDegrees <= southEdge || latitudeDegrees >= northEdge)
        return pressure(latitudeDegrees);
    const float edge = latitudeDegrees >= targetTrough ? northEdge : southEdge;
    const float mapped = controlTrough + (latitudeDegrees - targetTrough) *
        (edge - controlTrough) / (edge - targetTrough);
    return pressure(mapped);
}

// Calibrated marine weighting of the negative extratropical template lobe.
// This acts on the UNCENTRED template, before the caller removes its gauge.
// Continental thermal/heating departures and the tropical trough are retained.
// It approximates geography-dependent eddy forcing, not resolved storm tracks.
inline float continentalSubpolarCorrectionHpa(float uncentredPressureHpa,
    float latitudeDegrees, float controlTrough, float baseWidthDegrees,
    float anchorFraction, float landFraction, float strength)
{
    const float pole = latitudeDegrees >= controlTrough ? 90.0f : -90.0f;
    const float edge = std::min(winterAnchoredHadleyWidth(pole, controlTrough,
        baseWidthDegrees, anchorFraction), 0.8f * std::abs(pole - controlTrough));
    if (std::abs(latitudeDegrees - controlTrough) <= edge || uncentredPressureHpa >= 0.0f)
        return 0.0f;
    return -uncentredPressureHpa * std::clamp(landFraction, 0.0f, 1.0f) *
        std::clamp(strength, 0.0f, 1.0f);
}

// Empirical eddy-forcing sensitivity to each hemisphere's thermal contrast.
// Use only the uncentred negative extratropical lobe; a zero exponent recovers
// equal hemispheric amplitudes. This is not a resolved storm-track closure.
inline float hemisphericSubpolarPressureHpa(float uncentredPressureHpa,
    float latitudeDegrees, float controlTrough, float baseWidthDegrees,
    float anchorFraction, float contrastRatio, float exponent)
{
    if (!(exponent > 0.0f)) return uncentredPressureHpa;
    const float depression = continentalSubpolarCorrectionHpa(uncentredPressureHpa,
        latitudeDegrees, controlTrough, baseWidthDegrees, anchorFraction, 1.0f, 1.0f);
    const float scale = std::pow(std::clamp(contrastRatio, 0.25f, 2.0f),
        std::clamp(exponent, 0.0f, 4.0f));
    return uncentredPressureHpa - depression * (scale - 1.0f);
}

// Separate the subpolar depression from the tropical overturning amplitude.
// This smooth, calibrated shape is not a prognostic eddy-momentum closure.
// The caller removes the resulting global area mean, preserving pressure gauge.
inline float deepenSubpolarLowPressureHpa(float pressureHpa, float latitudeDegrees,
    float controlTrough, float baseWidthDegrees, float anchorFraction,
    float amplitudeHpa, float extraFraction, bool preservePolarGradient = false)
{
    if (!(amplitudeHpa > 0.0f) || !(extraFraction > 0.0f)) return pressureHpa;
    const float pole = latitudeDegrees >= controlTrough ? 90.0f : -90.0f;
    const float span = std::abs(pole - controlTrough);
    const float edge = std::min(winterAnchoredHadleyWidth(pole, controlTrough,
        baseWidthDegrees, anchorFraction), 0.8f * span);
    const float distance = std::abs(latitudeDegrees - controlTrough);
    if (distance <= edge || (distance >= span && !preservePolarGradient)) return pressureHpa;
    const float low = edge + 0.65f * (span - edge);
    // A constant offset beyond the low avoids adding a spurious polar slope.
    // Retain the endpoint-pinned alternative for the G6 sensitivity replay.
    const float position = distance <= low ? (distance - edge) / (low - edge) :
        preservePolarGradient ? 1.0f : (span - distance) / (span - low);
    const float wave = std::sin(1.5707963267948966f * position);
    return pressureHpa - amplitudeHpa * std::clamp(extraFraction, 0.0f, 1.0f) * wave * wave;
}

// Adjust only the subpolar-low-to-pole pressure difference. Its low and the
// entire tropical/midlatitude forcing remain fixed before gauge removal.
inline float polarBranchPressureHpa(float pressureHpa, float latitudeDegrees,
    float controlTrough, float baseWidthDegrees, float anchorFraction,
    float amplitudeHpa, float contrastRatio, float strength)
{
    if (!(strength > 0.0f)) return pressureHpa;
    const bool north = latitudeDegrees >= controlTrough;
    const float pole = north ? 90.0f : -90.0f;
    const float span = std::abs(pole - controlTrough);
    const float edge = std::min(winterAnchoredHadleyWidth(pole,controlTrough,baseWidthDegrees,anchorFraction),.8f*span);
    const float lowDistance = edge + .65f * (span-edge);
    if (std::abs(latitudeDegrees-controlTrough) <= lowDistance) return pressureHpa;
    const float scale = 1.0f + std::clamp(strength,0.0f,1.0f) * (std::clamp(contrastRatio,.25f,2.0f)-1.0f);
    return -.85f*amplitudeHpa + (pressureHpa+.85f*amplitudeHpa)*scale;
}
}
