#pragma once
#include <array>

namespace climatekoppen
{
inline constexpr float winterDrynessDivisor = 10.0f;

constexpr bool isWinterDry(float driestColdMonthMm, float wettestWarmMonthMm)
{
    return driestColdMonthMm < wettestWarmMonthMm / winterDrynessDivisor;
}

// Monthly means/totals. Uses the Beck 2018 0 C C/D boundary and standard
// Apr-Sep northern summer (reversed in the southern hemisphere).
short classifyMonthly(const std::array<float, 12>& temperatureC,
    const std::array<float, 12>& precipitationMm, bool northernHemisphere);
}
