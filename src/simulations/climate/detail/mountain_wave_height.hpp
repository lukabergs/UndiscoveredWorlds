#pragma once

#include <algorithm>
#include <cmath>

namespace climateatmosphere::detail
{
// Experimental launch-amplitude limit: N*h/|U| <= 1. This approximates the
// effective height available to a linear mountain wave when lower flow blocks.
// The caller currently has only a bulk background wind and stability; this is
// not a vertically resolved blocking, drag or wave-breaking scheme.
// Effective-height precedent: https://acp.copernicus.org/articles/20/12483/2020/
inline float effectiveMountainWaveHeightMetres(float heightMetres,
    float backgroundWindMps, float stabilityPerSecond)
{
    if (!std::isfinite(heightMetres) || !std::isfinite(backgroundWindMps) ||
        !std::isfinite(stabilityPerSecond) || stabilityPerSecond <= 0.0f)
        return 0.0f;
    return std::min(std::max(0.0f, heightMetres),
        std::abs(backgroundWindMps) / stabilityPerSecond);
}
}
