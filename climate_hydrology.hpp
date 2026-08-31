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

CalendarMonth calendarMonth(int month);
float interpolateSeasonal(float first, float second, float interpolation);
float soilMoistureStress(
    float soilMoistureMm,
    float soilMoistureCapacityMm,
    float criticalCapacityFraction,
    float exponent);

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
}
