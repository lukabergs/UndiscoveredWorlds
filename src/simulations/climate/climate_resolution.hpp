#pragma once

// Distances of the climate rules are written in cells of the grid they were
// tuned on (tuning::climate::resolution); these scale them to another grid.

#include "simulations/climate/climate_parameters.hpp"

namespace tuning::climateresolution
{
inline float gridScale(int width, int height)
{
    const float widthscale = static_cast<float>(width + 1) / static_cast<float>(climate::resolution::referenceWidth + 1);
    const float heightscale = static_cast<float>(height + 1) /
        static_cast<float>(climate::resolution::referenceHeight + 1);
    return (widthscale + heightscale) / 2.0f;
}

inline int scaleDistance(int distance, int width, int height)
{
    const int scaled = static_cast<int>(static_cast<float>(distance) * gridScale(width, height) + 0.5f);
    return scaled > 0 ? scaled : 1;
}

inline float scaleDistance(float distance, int width, int height)
{
    return distance * gridScale(width, height);
}
}
