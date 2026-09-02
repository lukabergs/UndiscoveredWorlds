#pragma once

#include <vector>

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

struct ClimateGridDimensions
{
    int columns = 1;
    int rows = 1;
};

struct WeatherPhase
{
    float windRotationRadians = 0.0f;
    float synopticPhaseRadians = 0.0f;
    float coastalDirection = 0.0f;
    float landTemperatureAnomalyC = 0.0f;
    float seaTemperatureAnomalyC = 0.0f;
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

struct SphericalTracerTransportDiagnostics
{
    int substeps = 1;
    int correctivePasses = 0;
    float maximumMultidimensionalCourant = 0.0f;
    float maximumZonalCourant = 0.0f;
    float maximumMeridionalCourant = 0.0f;
    double initialAreaWeightedMass = 0.0;
    double finalAreaWeightedMass = 0.0;
    float minimumMixingRatio = 0.0f;
};

struct MpdataOptions
{
    float maximumCourantPerSubstep = 0.80f;
    int correctivePasses = 1;
    bool monotone = true;
    // Optional end-of-interval winds: linearly interpolated at EVERY substep
    // midpoint. Empty means a stationary velocity during this interval.
    std::vector<float> endZonalWindMps;
    std::vector<float> endMeridionalWindMps;
};

CalendarMonth calendarMonth(int month);
ClimateGridDimensions climateGridDimensions(
    int outputColumns,
    int outputRows,
    int targetHorizontalCells);
float climateCellLatitudeDegrees(int row, int rowCount);
float climateCellAreaWeight(int row, int rowCount);
float polarTaperFactor(
    float latitudeDegrees,
    float taperStartDegrees,
    float taperEndDegrees);
int adjacentMeridionalTransportTargetRow(
    int sourceRow,
    float displacementRows,
    int maximumRow);
SphericalTracerTransportDiagnostics advectSphericalTracer(
    int columns,
    int rows,
    const std::vector<float>& source,
    const std::vector<float>& zonalWindMps,
    const std::vector<float>& meridionalWindMps,
    float timeStepSeconds,
    float planetRadiusMetres,
    float maximumMeridionalCourantPerSubstep,
    float maximumDisplacementCells,
    std::vector<float>& destination);
SphericalTracerTransportDiagnostics advectSphericalTracerMpdata(
    int columns,
    int rows,
    const std::vector<float>& source,
    const std::vector<float>& zonalWindMps,
    const std::vector<float>& meridionalWindMps,
    float timeStepSeconds,
    float planetRadiusMetres,
    const MpdataOptions& options,
    std::vector<float>& destination);
WeatherPhase deterministicWeatherPhase(
    int phase,
    int phaseCount,
    float windRotationDegrees,
    float daytimeLandTemperatureAnomalyC,
    float nighttimeLandTemperatureAnomalyC,
    float daytimeSeaTemperatureAnomalyC,
    float nighttimeSeaTemperatureAnomalyC);
float interpolateSeasonal(float first, float second, float interpolation);
float kuoPrecipitationEfficiency(
    float relativeHumidity,
    float criticalRelativeHumidity,
    float humidityExponent);
float convectiveBuoyancyEfficiency(
    float parcelBuoyancyC,
    float activationBuoyancyC,
    float fullStrengthBuoyancyC);
float shallowConvectionExchangeFraction(
    float boundaryRelativeHumidity,
    float freeTroposphereRelativeHumidity,
    float parcelBuoyancyC,
    float verticalWindShearMps,
    float timeStepSeconds,
    float mixingTimeDays,
    float humidityOnset,
    float fullHumidity,
    float fullShearMps,
    float maximumExchangeFraction);
float dryConvectionExchangeFraction(
    float parcelBuoyancyC,
    float timeStepSeconds,
    float mixingTimeDays,
    float activationBuoyancyC,
    float fullStrengthBuoyancyC,
    float maximumExchangeFraction);
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
    float kuoCriticalRelativeHumidity,
    float convectiveConversionEfficiency,
    float convectiveActivationBuoyancyC,
    float convectiveFullStrengthBuoyancyC,
    int moistAdjustmentIterations,
    float latentHeatingCPerMillimetre,
    float capacityTemperatureSensitivityPerC,
    float signedFreeTroposphereConvergenceMm = 0.0f,
    float elevatedMoistureAccessionFraction = 0.0f,
    float kuoHumidityExponent = 3.0f);
}
