#include "climate_physics.hpp"

#include <algorithm>
#include <cmath>

namespace climatephysics
{
namespace
{
constexpr float standardSeaLevelPressureHpa = 1013.25f;
constexpr float pressureScaleHeightMetres = 8434.5f;
constexpr float gravityMetresPerSecondSquared = 9.80665f;
constexpr float dryAirGasConstantJoulesPerKilogramKelvin = 287.05f;
constexpr float activeMoistureColumnFraction = 0.28f;

std::array<WaterBudget, CLIMATESEASONCOUNT> waterBudgets{};
std::array<WaterBudget, CLIMATESEASONCOUNT> areaWeightedWaterBudgets{};
HydrologySpinupDiagnostics hydrologySpinupDiagnostics{};
PrecipitationDistributionDiagnostics precipitationDistributionDiagnostics{};
std::array<CondensationActivityDiagnostics, CLIMATESEASONCOUNT>
    condensationActivityDiagnostics{};
std::array<PrecipitationProcessDiagnostics, CLIMATESEASONCOUNT>
    precipitationProcessDiagnostics{};
}

double WaterBudget::residual() const
{
    return initialAtmosphericStorage + initialSoilStorage + initialSnowStorage +
        oceanEvaporation - oceanPrecipitation - runoff - atmosphericStorage -
        soilStorage - snowStorage;
}

double WaterBudget::relativeResidual() const
{
    const double totalInput = initialAtmosphericStorage + initialSoilStorage +
        initialSnowStorage + oceanEvaporation;
    return residual() / std::max(1.0, std::abs(totalInput));
}

float saturationVapourPressureHpa(float temperatureC)
{
    if (temperatureC >= 0.0f)
    {
        return 6.1121f * std::exp(
            (18.678f - temperatureC / 234.5f) *
            (temperatureC / (257.14f + temperatureC)));
    }

    return 6.1115f * std::exp(
        (23.036f - temperatureC / 333.7f) *
        (temperatureC / (279.82f + temperatureC)));
}

float surfacePressureHpa(float elevationAboveSeaLevelMetres)
{
    return surfacePressureHpa(elevationAboveSeaLevelMetres, 1.0f);
}

float surfacePressureHpa(float elevationAboveSeaLevelMetres, float gravityMultiplier)
{
    const float elevation = std::max(0.0f, elevationAboveSeaLevelMetres);
    const float gravityScale = std::max(0.05f, gravityMultiplier);
    return standardSeaLevelPressureHpa *
        std::exp(-elevation * gravityScale / pressureScaleHeightMetres);
}

float saturationSpecificHumidity(float temperatureC, float pressureHpa)
{
    const float saturationPressureHpa = std::min(
        saturationVapourPressureHpa(temperatureC), pressureHpa * 0.95f);
    return 0.622f * saturationPressureHpa /
        std::max(1.0f, pressureHpa - 0.378f * saturationPressureHpa);
}

float saturationColumnWaterAtPressure(float temperatureC, float surfacePressureHpa)
{
    return saturationColumnWaterAtPressure(temperatureC, surfacePressureHpa, 1.0f);
}

float saturationColumnWaterAtPressure(
    float temperatureC,
    float surfacePressureHpa,
    float gravityMultiplier)
{
    const float pressureHpa = std::max(1.0f, surfacePressureHpa);
    const float specificHumidity = saturationSpecificHumidity(temperatureC, pressureHpa);
    const float gravity = gravityMetresPerSecondSquared * std::max(0.05f, gravityMultiplier);
    const float columnAirMass = pressureHpa * 100.0f / gravity;

    return std::max(
        0.05f,
        specificHumidity * columnAirMass * activeMoistureColumnFraction);
}

float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres)
{
    return saturationColumnWaterAtPressure(
        temperatureC,
        surfacePressureHpa(elevationAboveSeaLevelMetres));
}

float neutralDragCoefficient(float roughnessLengthMetres, float referenceHeightMetres)
{
    constexpr float vonKarmanConstant = 0.4f;
    const float roughness = std::max(1.0e-6f, roughnessLengthMetres);
    const float referenceHeight = std::max(roughness * 1.01f, referenceHeightMetres);
    const float logarithmicProfile = std::log(referenceHeight / roughness);
    return std::clamp(
        vonKarmanConstant * vonKarmanConstant /
            std::max(0.01f, logarithmicProfile * logarithmicProfile),
        0.0002f,
        0.02f);
}

float bulkRichardsonExchangeMultiplier(
    float surfaceTemperatureC,
    float airTemperatureC,
    float windSpeedMetresPerSecond,
    float gravityMultiplier,
    float referenceHeightMetres,
    float minimumMultiplier,
    float maximumMultiplier)
{
    const float absoluteTemperature = std::max(150.0f, surfaceTemperatureC + 273.15f);
    const float windSpeedSquared = std::max(
        0.25f,
        windSpeedMetresPerSecond * windSpeedMetresPerSecond);
    const float richardsonNumber =
        gravityMetresPerSecondSquared * std::max(0.05f, gravityMultiplier) *
        std::max(0.1f, referenceHeightMetres) *
        (airTemperatureC - surfaceTemperatureC) /
        (absoluteTemperature * windSpeedSquared);
    const float multiplier = richardsonNumber >= 0.0f
        ? 1.0f / (1.0f + 10.0f * richardsonNumber)
        : std::sqrt(std::max(0.0f, 1.0f - 16.0f * richardsonNumber));
    return std::clamp(
        multiplier,
        std::max(0.0f, minimumMultiplier),
        std::max(minimumMultiplier, maximumMultiplier));
}

float bulkAerodynamicEvaporationMmAtPressure(
    float temperatureC,
    float surfacePressureHpa,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient)
{
    const float pressureHpa = std::max(1.0f, surfacePressureHpa);
    const float absoluteTemperatureKelvin = std::max(150.0f, temperatureC + 273.15f);
    const float airDensity = pressureHpa * 100.0f /
        (dryAirGasConstantJoulesPerKilogramKelvin * absoluteTemperatureKelvin);
    const float humidityDeficit = saturationSpecificHumidity(temperatureC, pressureHpa) *
        (1.0f - std::clamp(relativeHumidity, 0.0f, 1.0f));

    return std::max(0.0f,
        airDensity * std::max(0.0f, transferCoefficient) *
        std::max(0.0f, windSpeedMetresPerSecond) * humidityDeficit *
        std::max(0.0f, timeStepSeconds));
}

float bulkAerodynamicEvaporationMm(
    float temperatureC,
    float elevationAboveSeaLevelMetres,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient)
{
    return bulkAerodynamicEvaporationMmAtPressure(
        temperatureC,
        surfacePressureHpa(elevationAboveSeaLevelMetres),
        windSpeedMetresPerSecond,
        relativeHumidity,
        timeStepSeconds,
        transferCoefficient);
}

float transientEddyDiffusivityM2S(
    float surfaceUWindMetresPerSecond,
    float surfaceVWindMetresPerSecond,
    float upperUWindMetresPerSecond,
    float upperVWindMetresPerSecond,
    float latitudeDegrees,
    float mixingLengthMetres,
    float maximumDiffusivityM2S,
    float minimumLatitudeDegrees,
    float fullStrengthLatitudeDegrees)
{
    constexpr float pi = 3.14159265358979323846f;
    const float ushear = upperUWindMetresPerSecond - surfaceUWindMetresPerSecond;
    const float vshear = upperVWindMetresPerSecond - surfaceVWindMetresPerSecond;
    const float shearspeed = std::sqrt(ushear * ushear + vshear * vshear);
    const float absolutelatitude = std::abs(latitudeDegrees);
    const float latituderange = std::max(0.001f, fullStrengthLatitudeDegrees - minimumLatitudeDegrees);
    const float latitudephase = std::clamp(
        (absolutelatitude - minimumLatitudeDegrees) / latituderange,
        0.0f,
        1.0f);
    const float baroclinicweight = latitudephase * latitudephase * (3.0f - 2.0f * latitudephase);
    const float rotationweight = std::abs(std::sin(latitudeDegrees * pi / 180.0f)) * baroclinicweight;
    return std::clamp(
        std::max(0.0f, mixingLengthMetres) * shearspeed * rotationweight,
        0.0f,
        std::max(0.0f, maximumDiffusivityM2S));
}

float relaxedExcessCondensationAmount(
    float availableColumnWater,
    float criticalColumnWater,
    float timeStepSeconds,
    float conversionTimeSeconds)
{
    if (timeStepSeconds <= 0.0f || conversionTimeSeconds <= 0.0f)
        return 0.0f;

    const float excessWater = std::max(
        0.0f,
        availableColumnWater - std::max(0.0f, criticalColumnWater));
    const float convertedfraction = -std::expm1(-timeStepSeconds / conversionTimeSeconds);
    return std::clamp(excessWater * convertedfraction, 0.0f, std::max(0.0f, availableColumnWater));
}

void setLastWaterBudget(int season, const WaterBudget& budget)
{
    if (season >= 0 && season < CLIMATESEASONCOUNT)
        waterBudgets[season] = budget;
}

const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastWaterBudgets()
{
    return waterBudgets;
}

void setLastAreaWeightedWaterBudget(int season, const WaterBudget& budget)
{
    if (season >= 0 && season < CLIMATESEASONCOUNT)
        areaWeightedWaterBudgets[season] = budget;
}

const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastAreaWeightedWaterBudgets()
{
    return areaWeightedWaterBudgets;
}

void setLastHydrologySpinupDiagnostics(const HydrologySpinupDiagnostics& diagnostics)
{
    hydrologySpinupDiagnostics = diagnostics;
}

const HydrologySpinupDiagnostics& lastHydrologySpinupDiagnostics()
{
    return hydrologySpinupDiagnostics;
}

void setLastPrecipitationDistributionDiagnostics(const PrecipitationDistributionDiagnostics& diagnostics)
{
    precipitationDistributionDiagnostics = diagnostics;
}

const PrecipitationDistributionDiagnostics& lastPrecipitationDistributionDiagnostics()
{
    return precipitationDistributionDiagnostics;
}

void setLastCondensationActivityDiagnostics(
    int season,
    const CondensationActivityDiagnostics& diagnostics)
{
    if (season >= 0 && season < CLIMATESEASONCOUNT)
        condensationActivityDiagnostics[season] = diagnostics;
}

const std::array<CondensationActivityDiagnostics, CLIMATESEASONCOUNT>&
lastCondensationActivityDiagnostics()
{
    return condensationActivityDiagnostics;
}

void setLastPrecipitationProcessDiagnostics(
    int season,
    const PrecipitationProcessDiagnostics& diagnostics)
{
    if (season >= 0 && season < CLIMATESEASONCOUNT)
        precipitationProcessDiagnostics[season] = diagnostics;
}

const std::array<PrecipitationProcessDiagnostics, CLIMATESEASONCOUNT>&
lastPrecipitationProcessDiagnostics()
{
    return precipitationProcessDiagnostics;
}
}
