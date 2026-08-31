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

float FallingPrecipitation::surfaceTotalMm() const
{
    return rainMm + snowMm;
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

float diagnosticCloudFraction(float relativeHumidity, float cloudOnsetRelativeHumidity)
{
    return smoothStep(
        std::clamp(cloudOnsetRelativeHumidity, 0.0f, 1.0f),
        1.0f,
        std::clamp(relativeHumidity, 0.0f, 1.0f));
}

MoistAdjustment moistSaturationAdjustment(
    float availableColumnWaterMm,
    float saturationCapacityMm,
    float temperatureC,
    float timeStepSeconds,
    float conversionTimeSeconds,
    int iterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC)
{
    MoistAdjustment result;
    result.remainingVapourMm = std::max(0.0f, availableColumnWaterMm);
    result.adjustedTemperatureC = temperatureC;
    const float initialCapacity = std::max(0.0f, saturationCapacityMm);

    const int iterationCount = std::max(1, iterations);

    for (int iteration = 0; iteration < iterationCount; iteration++)
    {
        const float warming = result.adjustedTemperatureC - temperatureC;
        const float adjustedCapacity = initialCapacity * std::exp(
            std::clamp(capacityTemperatureSensitivityPerC * warming, -20.0f, 20.0f));
        const float condensed = relaxedExcess(
            result.remainingVapourMm,
            adjustedCapacity,
            timeStepSeconds / static_cast<float>(iterationCount),
            conversionTimeSeconds);

        result.condensedMm += condensed;
        result.remainingVapourMm -= condensed;
        result.adjustedTemperatureC += condensed *
            std::max(0.0f, latentHeatingCPerMillimetre);
    }

    result.condensedMm = std::clamp(
        result.condensedMm,
        0.0f,
        std::max(0.0f, availableColumnWaterMm));
    result.remainingVapourMm = std::max(
        0.0f,
        std::max(0.0f, availableColumnWaterMm) - result.condensedMm);
    return result;
}

MoistureLayerExchange exchangeMoistureLayers(
    float boundaryLayerMm,
    float freeTroposphereMm,
    float upwardFraction,
    float downwardFraction)
{
    const float boundary = std::max(0.0f, boundaryLayerMm);
    const float free = std::max(0.0f, freeTroposphereMm);
    const float upward = boundary * std::clamp(upwardFraction, 0.0f, 1.0f);
    const float downward = free * std::clamp(downwardFraction, 0.0f, 1.0f);

    return {
        boundary - upward + downward,
        free + upward - downward,
        upward,
        downward
    };
}

FallingPrecipitation processFallingPrecipitation(
    float condensateMm,
    float surfaceTemperatureC,
    float boundaryRelativeHumidity,
    float maximumReevaporationFraction,
    float maximumVapourUptakeMm,
    float allSnowTemperatureC,
    float allRainTemperatureC)
{
    const float condensate = std::max(0.0f, condensateMm);
    const float humidityDeficit = 1.0f - std::clamp(boundaryRelativeHumidity, 0.0f, 1.0f);
    const float reevaporation = std::min({
        condensate,
        std::max(0.0f, maximumVapourUptakeMm),
        condensate * std::clamp(maximumReevaporationFraction, 0.0f, 1.0f) *
            humidityDeficit * humidityDeficit
    });
    const float reachingSurface = condensate - reevaporation;
    const float rainFraction = smoothStep(
        allSnowTemperatureC,
        allRainTemperatureC,
        surfaceTemperatureC);

    return {
        reachingSurface * rainFraction,
        reachingSurface * (1.0f - rainFraction),
        reevaporation
    };
}

float snowMeltAmount(
    float snowWaterEquivalentMm,
    float surfaceTemperatureC,
    float timeStepSeconds,
    float degreeDayMeltMmPerDegreeC)
{
    constexpr float secondsPerDay = 86400.0f;
    const float potentialMelt = std::max(0.0f, surfaceTemperatureC) *
        std::max(0.0f, degreeDayMeltMmPerDegreeC) *
        std::max(0.0f, timeStepSeconds) / secondsPerDay;
    return std::min(std::max(0.0f, snowWaterEquivalentMm), potentialMelt);
}

SnowAccumulation accumulateSnowfall(
    float snowWaterEquivalentMm,
    float snowfallMm,
    float maximumSnowStorageMm)
{
    const float combined = std::max(0.0f, snowWaterEquivalentMm) +
        std::max(0.0f, snowfallMm);
    const float storage = std::min(combined, std::max(0.0f, maximumSnowStorageMm));
    return { storage, combined - storage };
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

PrecipitationPartition partitionTwoLayerPrecipitation(
    float boundaryLayerWaterMm,
    float freeTroposphereWaterMm,
    float boundaryLayerSaturationCapacityMm,
    float nonOrographicFreeTroposphereCapacityMm,
    float terrainAdjustedFreeTroposphereCapacityMm,
    float signedBoundaryLayerConvergenceMm,
    float surfaceEvaporationMm,
    float surfaceTemperatureC,
    float freeTroposphereTemperatureC,
    float timeStepSeconds,
    float stratiformConversionTimeSeconds,
    float convectiveResidualRelativeHumidity,
    float convectiveConversionEfficiency,
    float convectiveActivationTemperatureC,
    float convectiveFullStrengthTemperatureC,
    int moistAdjustmentIterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC)
{
    const float boundaryWater = std::max(0.0f, boundaryLayerWaterMm);
    const float freeWater = std::max(0.0f, freeTroposphereWaterMm);
    const MoistAdjustment nonOrographicAdjustment = moistSaturationAdjustment(
        freeWater,
        nonOrographicFreeTroposphereCapacityMm,
        freeTroposphereTemperatureC,
        timeStepSeconds,
        stratiformConversionTimeSeconds,
        moistAdjustmentIterations,
        latentHeatingCPerMillimetre,
        capacityTemperatureSensitivityPerC);
    const MoistAdjustment terrainAdjustment = moistSaturationAdjustment(
        freeWater,
        terrainAdjustedFreeTroposphereCapacityMm,
        freeTroposphereTemperatureC,
        timeStepSeconds,
        stratiformConversionTimeSeconds,
        moistAdjustmentIterations,
        latentHeatingCPerMillimetre,
        capacityTemperatureSensitivityPerC);

    PrecipitationPartition result;
    result.orographicMm = std::max(
        0.0f,
        terrainAdjustment.condensedMm - nonOrographicAdjustment.condensedMm);
    result.stratiformMm = std::max(
        0.0f,
        terrainAdjustment.condensedMm - result.orographicMm);

    const float convectiveFloor = std::max(
        0.0f,
        boundaryLayerSaturationCapacityMm *
            std::clamp(convectiveResidualRelativeHumidity, 0.0f, 1.0f));
    const float convectivelyAvailable = std::max(0.0f, boundaryWater - convectiveFloor);
    const float instability = smoothStep(
        convectiveActivationTemperatureC,
        convectiveFullStrengthTemperatureC,
        surfaceTemperatureC);
    const float convergentSupply = std::max(
        0.0f,
        signedBoundaryLayerConvergenceMm + std::max(0.0f, surfaceEvaporationMm));
    result.convectiveMm = std::min(
        convectivelyAvailable,
        convergentSupply * instability *
            std::clamp(convectiveConversionEfficiency, 0.0f, 1.0f));
    return result;
}
}
