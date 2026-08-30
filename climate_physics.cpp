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
HydrologySpinupDiagnostics hydrologySpinupDiagnostics{};
}

double WaterBudget::residual() const
{
    return initialAtmosphericStorage + initialSoilStorage + oceanEvaporation -
        oceanPrecipitation - runoff - atmosphericStorage - soilStorage;
}

double WaterBudget::relativeResidual() const
{
    const double totalInput = initialAtmosphericStorage + initialSoilStorage + oceanEvaporation;
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
    const float elevation = std::max(0.0f, elevationAboveSeaLevelMetres);
    return standardSeaLevelPressureHpa * std::exp(-elevation / pressureScaleHeightMetres);
}

float saturationSpecificHumidity(float temperatureC, float pressureHpa)
{
    const float saturationPressureHpa = std::min(
        saturationVapourPressureHpa(temperatureC), pressureHpa * 0.95f);
    return 0.622f * saturationPressureHpa /
        std::max(1.0f, pressureHpa - 0.378f * saturationPressureHpa);
}

float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres)
{
    const float pressureHpa = surfacePressureHpa(elevationAboveSeaLevelMetres);
    const float specificHumidity = saturationSpecificHumidity(temperatureC, pressureHpa);
    const float columnAirMass = pressureHpa * 100.0f / gravityMetresPerSecondSquared;

    return std::max(
        0.05f,
        specificHumidity * columnAirMass * activeMoistureColumnFraction);
}

float bulkAerodynamicEvaporationMm(
    float temperatureC,
    float elevationAboveSeaLevelMetres,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient)
{
    const float pressureHpa = surfacePressureHpa(elevationAboveSeaLevelMetres);
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

void setLastWaterBudget(int season, const WaterBudget& budget)
{
    if (season >= 0 && season < CLIMATESEASONCOUNT)
        waterBudgets[season] = budget;
}

const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastWaterBudgets()
{
    return waterBudgets;
}

void setLastHydrologySpinupDiagnostics(const HydrologySpinupDiagnostics& diagnostics)
{
    hydrologySpinupDiagnostics = diagnostics;
}

const HydrologySpinupDiagnostics& lastHydrologySpinupDiagnostics()
{
    return hydrologySpinupDiagnostics;
}
}
