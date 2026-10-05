#pragma once

#include <array>

namespace tuning
{
namespace worlddefaults
{
inline constexpr int size = 2;
inline constexpr int width = 599;
inline constexpr int height = 299;
inline constexpr bool rotation = true;
inline constexpr float tilt = 22.5f;
inline constexpr float eccentricity = 0.0167f;
inline constexpr short perihelion = 0;
inline constexpr float gravity = 1.0f;
inline constexpr float lunar = 1.0f;
inline constexpr float tempDecrease = 4.54f;
inline constexpr int northPolarAdjust = -4;
inline constexpr int southPolarAdjust = -4;
inline constexpr int averageTemp = 14;
inline constexpr float waterPickup = 1.0f;
inline constexpr float riverFactor = 15.0f;
inline constexpr int riverLandReduce = 20;
inline constexpr int estuaryLimit = 20;
inline constexpr int glacialTemp = 4;
inline constexpr int glacierTemp = -1;
inline constexpr float mountainReduce = 0.75f;
inline constexpr int climateCount = 31;
inline constexpr int maxElevation = 48000;
inline constexpr int seaLevel = 24000;
inline constexpr int craterCount = 0;
}

namespace climateresolution
{
inline constexpr float gridScale(int width, int height)
{
    const float widthscale = static_cast<float>(width + 1) / static_cast<float>(worlddefaults::width + 1);
    const float heightscale = static_cast<float>(height + 1) /
        static_cast<float>(worlddefaults::height + 1);
    return (widthscale + heightscale) / 2.0f;
}

inline constexpr int scaleDistance(int distance, int width, int height)
{
    const int scaled = static_cast<int>(static_cast<float>(distance) * gridScale(width, height) + 0.5f);
    return scaled > 0 ? scaled : 1;
}

inline constexpr float scaleDistance(float distance, int width, int height)
{
    return distance * gridScale(width, height);
}

static_assert(scaleDistance(30, worlddefaults::width, worlddefaults::height) == 30);
static_assert(scaleDistance(30, 2047, 1023) == 102);
}

namespace regional
{
inline constexpr float submarineWarpFactor = 100.0f;
inline constexpr int submarineExtraMargin = 20;
inline constexpr int ridgeRadiationExtraMargin = 20;
inline constexpr int saltPanMaxRain = 30;
inline constexpr float mountainPrecipitationSlopeFactor = 400.0f;
}
}
