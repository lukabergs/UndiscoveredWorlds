#include "climate_koppen.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace climatekoppen
{
short classifyMonthly(const std::array<float, 12>& temperature,
    const std::array<float, 12>& precipitation, bool northern)
{
    float minimum = temperature[0], maximum = temperature[0], mean = 0.0f, annual = 0.0f;
    float summerRain = 0.0f, summerMin = std::numeric_limits<float>::max(), winterMin = summerMin;
    float summerMax = 0.0f, winterMax = 0.0f;
    int warmMonths = 0;
    for (int month = 0; month < 12; ++month)
    {
        if (!std::isfinite(temperature[month]) || !std::isfinite(precipitation[month]) || precipitation[month] < 0.0f) return 0;
        minimum = std::min(minimum, temperature[month]); maximum = std::max(maximum, temperature[month]);
        mean += temperature[month] / 12.0f; annual += precipitation[month];
        warmMonths += temperature[month] > 10.0f;
        const bool summer = (month >= 3 && month <= 8) == northern;
        if (summer)
        {
            summerRain += precipitation[month];
            summerMin = std::min(summerMin, precipitation[month]); summerMax = std::max(summerMax, precipitation[month]);
        }
        else
        {
            winterMin = std::min(winterMin, precipitation[month]); winterMax = std::max(winterMax, precipitation[month]);
        }
    }
    const float summerFraction = annual > 0.0f ? summerRain / annual : 0.0f;
    const float dryThreshold = 20.0f * mean + (summerFraction >= 0.7f ? 280.0f : (summerFraction <= 0.3f ? 0.0f : 140.0f));
    if (annual < dryThreshold) return static_cast<short>((annual < 0.5f * dryThreshold ? 5 : 7) + (mean < 18.0f));
    if (maximum < 10.0f) return maximum < 0.0f ? 31 : 30;
    if (minimum >= 18.0f)
    {
        const float driest = std::min(summerMin, winterMin);
        if (driest >= 60.0f) return 1;
        if (driest >= 100.0f - annual / 25.0f) return 2;
        return summerMin < winterMin ? 4 : 3;
    }
    const int moisture = summerMin < 40.0f && summerMin < winterMax / 3.0f ? 0 :
        (isWinterDry(winterMin, summerMax) ? 1 : 2);
    const int heat = maximum >= 22.0f ? 0 : (warmMonths >= 4 ? 1 : (minimum < -38.0f ? 3 : 2));
    if (minimum > 0.0f) return static_cast<short>(9 + 3 * moisture + std::min(heat, 2));
    return static_cast<short>(18 + 4 * moisture + heat);
}
}
