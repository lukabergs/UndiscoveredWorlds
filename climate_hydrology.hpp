#pragma once

namespace climatehydrology
{
inline constexpr int monthCount = 12;
inline constexpr int seasonCount = 4;

struct CalendarMonth
{
    int month = 0;
    int firstSeason = 0;
    int secondSeason = 0;
    float interpolation = 0.0f;
    int days = 31;
};

struct PrecipitationPartition
{
    float stratiformMm = 0.0f;
    float orographicMm = 0.0f;
    float convectiveMm = 0.0f;

    float totalMm() const;
};

struct MoistAdjustment
{
    float condensedMm = 0.0f;
    float remainingVapourMm = 0.0f;
    float adjustedTemperatureC = 0.0f;
};

struct MoistureLayerExchange
{
    float boundaryLayerMm = 0.0f;
    float freeTroposphereMm = 0.0f;
    float upwardTransferMm = 0.0f;
    float downwardTransferMm = 0.0f;
};

struct FallingPrecipitation
{
    float rainMm = 0.0f;
    float snowMm = 0.0f;
    float reevaporatedMm = 0.0f;

    float surfaceTotalMm() const;
};

struct SnowAccumulation
{
    float storageMm = 0.0f;
    float overflowMm = 0.0f;
};

CalendarMonth calendarMonth(int month);
float interpolateSeasonal(float first, float second, float interpolation);
float soilMoistureStress(
    float soilMoistureMm,
    float soilMoistureCapacityMm,
    float criticalCapacityFraction,
    float exponent);
float diagnosticCloudFraction(float relativeHumidity, float cloudOnsetRelativeHumidity);
MoistAdjustment moistSaturationAdjustment(
    float availableColumnWaterMm,
    float saturationCapacityMm,
    float temperatureC,
    float timeStepSeconds,
    float conversionTimeSeconds,
    int iterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC);
MoistureLayerExchange exchangeMoistureLayers(
    float boundaryLayerMm,
    float freeTroposphereMm,
    float upwardFraction,
    float downwardFraction);
FallingPrecipitation processFallingPrecipitation(
    float condensateMm,
    float surfaceTemperatureC,
    float boundaryRelativeHumidity,
    float maximumReevaporationFraction,
    float maximumVapourUptakeMm,
    float allSnowTemperatureC,
    float allRainTemperatureC);
float snowMeltAmount(
    float snowWaterEquivalentMm,
    float surfaceTemperatureC,
    float timeStepSeconds,
    float degreeDayMeltMmPerDegreeC);
SnowAccumulation accumulateSnowfall(
    float snowWaterEquivalentMm,
    float snowfallMm,
    float maximumSnowStorageMm);

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
    float convectiveFullStrengthTemperatureC);
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
    float capacityTemperatureSensitivityPerC);
}
