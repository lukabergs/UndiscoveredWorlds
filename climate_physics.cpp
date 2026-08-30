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
constexpr float activeMoistureColumnFraction = 0.28f;

std::array<WaterBudget, CLIMATESEASONCOUNT> waterBudgets{};
}

double WaterBudget::residual() const
{
    return oceanEvaporation - oceanPrecipitation - runoff - atmosphericStorage - soilStorage;
}

double WaterBudget::relativeResidual() const
{
    return residual() / std::max(1.0, std::abs(oceanEvaporation));
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

float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres)
{
    const float pressureHpa = surfacePressureHpa(elevationAboveSeaLevelMetres);
    const float saturationPressureHpa = std::min(
        saturationVapourPressureHpa(temperatureC), pressureHpa * 0.95f);
    const float saturationSpecificHumidity =
        0.622f * saturationPressureHpa /
        std::max(1.0f, pressureHpa - 0.378f * saturationPressureHpa);
    const float columnAirMass = pressureHpa * 100.0f / gravityMetresPerSecondSquared;

    return std::max(
        0.05f,
        saturationSpecificHumidity * columnAirMass * activeMoistureColumnFraction);
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
}
