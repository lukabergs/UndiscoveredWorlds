#pragma once

#include <array>
#include <vector>

#include "physical_layers.hpp"

namespace climatephysics
{
struct ClimateCouplingDiagnostics
{
    int iteration = 0;
    double windChange = 0.0, sstChange = 0.0, heatingChange = 0.0, rainfallChange = 0.0;
    bool innerSolvesAccepted = false, converged = false;
};
bool climateCouplingConverged(const ClimateCouplingDiagnostics& diagnostics, int minimumIterations, double tolerance);
void clearClimateCouplingDiagnostics();
void appendClimateCouplingDiagnostics(const ClimateCouplingDiagnostics& diagnostics);
const std::vector<ClimateCouplingDiagnostics>& lastClimateCouplingDiagnostics();
struct WaterBudget
{
    double initialAtmosphericStorage = 0.0;
    double initialSoilStorage = 0.0;
    double initialSnowStorage = 0.0;
    double oceanEvaporation = 0.0;
    double landEvaporation = 0.0;
    double oceanPrecipitation = 0.0;
    double landPrecipitation = 0.0;
    double runoff = 0.0;
    double atmosphericStorage = 0.0;
    double soilStorage = 0.0;
    double snowStorage = 0.0;

    double residual() const;
    double relativeResidual() const;
};

struct HydrologySpinupDiagnostics
{
    int cyclesCompleted = 0;
    bool converged = false;
    double relativeStorageChange = 0.0;
    double relativeAtmosphericStorageChange = 0.0;
    double relativeSoilStorageChange = 0.0;
    double relativeSnowStorageChange = 0.0;
    double relativeSnowCoverChange = 0.0;
    double atmosphericStorage = 0.0;
    double soilStorage = 0.0;
    double snowStorage = 0.0;
};

struct PrecipitationDistributionScope
{
    long long cells = 0;
    double areaWeight = 0.0;
    double rawZeroFraction = 0.0;
    double storedZeroFraction = 0.0;
    double belowOneMillimetreFraction = 0.0;
    double areaWeightedRawZeroFraction = 0.0;
    double areaWeightedStoredZeroFraction = 0.0;
    double areaWeightedBelowOneMillimetreFraction = 0.0;
    double meanMonthlyPrecipitationMm = 0.0;
    double areaWeightedMeanMonthlyPrecipitationMm = 0.0;
    double maximumMonthlyPrecipitationMm = 0.0;
    int maximumMonthlyPrecipitationX = -1;
    int maximumMonthlyPrecipitationY = -1;
    int maximumMonthlyPrecipitationMonth = -1;
    int roundedPrecipitationAtMaximumMm = 0;
    double wettestTenPercentShare = 0.0;
};

struct PrecipitationDistributionDiagnostics
{
    PrecipitationDistributionScope land;
    PrecipitationDistributionScope ocean;
};

struct CondensationActivityDiagnostics
{
    double cellStepAreaWeight = 0.0;
    double activeCellStepAreaWeight = 0.0;
    double atmosphericWaterAreaWeighted = 0.0;
    double activeAtmosphericWaterAreaWeighted = 0.0;
    double excessWaterAreaWeighted = 0.0;
    double precipitationAreaWeighted = 0.0;
};

struct PrecipitationProcessDiagnostics
{
    double stratiformPrecipitation = 0.0;
    double orographicPrecipitation = 0.0;
    double convectivePrecipitation = 0.0;
    double reevaporatedPrecipitation = 0.0;
    double snowfall = 0.0;
    double upwardMoistureTransfer = 0.0;
    double downwardMoistureTransfer = 0.0;
    double cloudFractionAreaWeighted = 0.0;
    double positiveMoistureFluxConvergence = 0.0;
    double negativeMoistureFluxConvergence = 0.0;
};

float saturationVapourPressureHpa(float temperatureC);
float surfacePressureHpa(float elevationAboveSeaLevelMetres);
float surfacePressureHpa(float elevationAboveSeaLevelMetres, float gravityMultiplier);
float saturationSpecificHumidity(float temperatureC, float pressureHpa);
float saturationColumnWaterAtPressure(float temperatureC, float surfacePressureHpa);
float saturationColumnWaterAtPressure(
    float temperatureC,
    float surfacePressureHpa,
    float gravityMultiplier);
float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres);
float neutralDragCoefficient(float roughnessLengthMetres, float referenceHeightMetres);
float bulkRichardsonExchangeMultiplier(
    float surfaceTemperatureC,
    float airTemperatureC,
    float windSpeedMetresPerSecond,
    float gravityMultiplier,
    float referenceHeightMetres,
    float minimumMultiplier,
    float maximumMultiplier);
float bulkAerodynamicEvaporationMmAtPressure(
    float temperatureC,
    float surfacePressureHpa,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient);
float bulkAerodynamicEvaporationMm(
    float temperatureC,
    float elevationAboveSeaLevelMetres,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient);
float transientEddyDiffusivityM2S(
    float surfaceUWindMetresPerSecond,
    float surfaceVWindMetresPerSecond,
    float upperUWindMetresPerSecond,
    float upperVWindMetresPerSecond,
    float latitudeDegrees,
    float mixingLengthMetres,
    float maximumDiffusivityM2S,
    float minimumLatitudeDegrees,
    float fullStrengthLatitudeDegrees);
float relaxedExcessCondensationAmount(
    float availableColumnWater,
    float criticalColumnWater,
    float timeStepSeconds,
    float conversionTimeSeconds);

void setLastWaterBudget(int season, const WaterBudget& budget);
const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastWaterBudgets();
void setLastAreaWeightedWaterBudget(int season, const WaterBudget& budget);
const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastAreaWeightedWaterBudgets();
void setLastHydrologySpinupDiagnostics(const HydrologySpinupDiagnostics& diagnostics);
const HydrologySpinupDiagnostics& lastHydrologySpinupDiagnostics();
void setLastPrecipitationDistributionDiagnostics(const PrecipitationDistributionDiagnostics& diagnostics);
const PrecipitationDistributionDiagnostics& lastPrecipitationDistributionDiagnostics();
void setLastCondensationActivityDiagnostics(
    int season,
    const CondensationActivityDiagnostics& diagnostics);
const std::array<CondensationActivityDiagnostics, CLIMATESEASONCOUNT>&
lastCondensationActivityDiagnostics();
void setLastPrecipitationProcessDiagnostics(
    int season,
    const PrecipitationProcessDiagnostics& diagnostics);
const std::array<PrecipitationProcessDiagnostics, CLIMATESEASONCOUNT>&
lastPrecipitationProcessDiagnostics();
}
