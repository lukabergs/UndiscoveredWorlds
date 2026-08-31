#include "climate_hydrology.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace climatehydrology
{
namespace
{
constexpr std::array<int, monthCount> daysInMonth = {
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

float relaxedExcess(
    float availableColumnWaterMm,
    float criticalColumnWaterMm,
    float timeStepSeconds,
    float conversionTimeSeconds)
{
    if (timeStepSeconds <= 0.0f || conversionTimeSeconds <= 0.0f)
        return 0.0f;

    const float excess = std::max(
        0.0f,
        availableColumnWaterMm - std::max(0.0f, criticalColumnWaterMm));
    const float convertedFraction = -std::expm1(-timeStepSeconds / conversionTimeSeconds);
    return std::clamp(
        excess * convertedFraction,
        0.0f,
        std::max(0.0f, availableColumnWaterMm));
}

float smoothStep(float lower, float upper, float value)
{
    if (upper <= lower)
        return value >= upper ? 1.0f : 0.0f;

    const float phase = std::clamp((value - lower) / (upper - lower), 0.0f, 1.0f);
    return phase * phase * (3.0f - 2.0f * phase);
}
}

float PrecipitationPartition::totalMm() const
{
    return stratiformMm + orographicMm + convectiveMm;
}

CalendarMonth calendarMonth(int month)
{
    const int wrappedMonth = ((month % monthCount) + monthCount) % monthCount;
    const int firstSeason = wrappedMonth / 3;
    return {
        wrappedMonth,
        firstSeason,
        (firstSeason + 1) % seasonCount,
        static_cast<float>(wrappedMonth % 3) / 3.0f,
        daysInMonth[wrappedMonth]
    };
}

float interpolateSeasonal(float first, float second, float interpolation)
{
    const float phase = std::clamp(interpolation, 0.0f, 1.0f);
    return first + (second - first) * phase;
}

float soilMoistureStress(
    float soilMoistureMm,
    float soilMoistureCapacityMm,
    float criticalCapacityFraction,
    float exponent)
{
    const float criticalStorage = std::max(
        0.001f,
        std::max(0.0f, soilMoistureCapacityMm) *
            std::clamp(criticalCapacityFraction, 0.001f, 1.0f));
    const float availableFraction = std::clamp(
        soilMoistureMm / criticalStorage,
        0.0f,
        1.0f);
    return std::pow(availableFraction, std::max(0.0f, exponent));
}

PrecipitationPartition partitionPrecipitation(
    float availableColumnWaterMm,
    float nonOrographicSaturationCapacityMm,
    float terrainAdjustedSaturationCapacityMm,
    float signedMoistureFluxConvergenceMm,
    float surfaceEvaporationMm,
    float surfaceTemperatureC,
    float timeStepSeconds,
    float stratiformCriticalRelativeHumidity,
    float stratiformConversionTimeSeconds,
    float convectiveResidualRelativeHumidity,
    float convectiveConversionEfficiency,
    float convectiveActivationTemperatureC,
    float convectiveFullStrengthTemperatureC)
{
    const float available = std::max(0.0f, availableColumnWaterMm);
    const float criticalRelativeHumidity = std::clamp(
        stratiformCriticalRelativeHumidity, 0.0f, 1.0f);
    const float nonOrographicThreshold = std::max(
        0.0f,
        nonOrographicSaturationCapacityMm * criticalRelativeHumidity);
    const float terrainAdjustedThreshold = std::max(
        0.0f,
        terrainAdjustedSaturationCapacityMm * criticalRelativeHumidity);
    const float nonOrographicCondensation = relaxedExcess(
        available,
        nonOrographicThreshold,
        timeStepSeconds,
        stratiformConversionTimeSeconds);
    const float terrainAdjustedCondensation = relaxedExcess(
        available,
        terrainAdjustedThreshold,
        timeStepSeconds,
        stratiformConversionTimeSeconds);

    PrecipitationPartition result;
    result.orographicMm = std::max(
        0.0f,
        terrainAdjustedCondensation - nonOrographicCondensation);
    result.stratiformMm = std::max(
        0.0f,
        terrainAdjustedCondensation - result.orographicMm);

    const float remaining = std::max(0.0f, available - terrainAdjustedCondensation);
    const float convectiveFloor = std::max(
        0.0f,
        nonOrographicSaturationCapacityMm *
            std::clamp(convectiveResidualRelativeHumidity, 0.0f, 1.0f));
    const float convectivelyAvailable = std::max(0.0f, remaining - convectiveFloor);
    const float instability = smoothStep(
        convectiveActivationTemperatureC,
        convectiveFullStrengthTemperatureC,
        surfaceTemperatureC);
    const float convergentSupply = std::max(
        0.0f,
        signedMoistureFluxConvergenceMm + std::max(0.0f, surfaceEvaporationMm));
    result.convectiveMm = std::min(
        convectivelyAvailable,
        convergentSupply * instability *
            std::clamp(convectiveConversionEfficiency, 0.0f, 1.0f));

    const float total = result.totalMm();
    if (total > available && total > 0.0f)
    {
        const float scale = available / total;
        result.stratiformMm *= scale;
        result.orographicMm *= scale;
        result.convectiveMm *= scale;
    }

    return result;
}
}
