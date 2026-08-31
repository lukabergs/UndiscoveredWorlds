#pragma once

namespace climatekoppen
{
inline constexpr float winterDrynessDivisor = 10.0f;

constexpr bool isWinterDry(float driestColdMonthMm, float wettestWarmMonthMm)
{
    return driestColdMonthMm < wettestWarmMonthMm / winterDrynessDivisor;
}
}
