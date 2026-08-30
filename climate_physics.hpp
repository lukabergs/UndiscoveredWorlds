#pragma once

#include <array>

#include "physical_layers.hpp"

namespace climatephysics
{
struct WaterBudget
{
    double initialAtmosphericStorage = 0.0;
    double initialSoilStorage = 0.0;
    double oceanEvaporation = 0.0;
    double landEvaporation = 0.0;
    double oceanPrecipitation = 0.0;
    double landPrecipitation = 0.0;
    double runoff = 0.0;
    double atmosphericStorage = 0.0;
    double soilStorage = 0.0;

    double residual() const;
    double relativeResidual() const;
};

struct HydrologySpinupDiagnostics
{
    int cyclesCompleted = 0;
    bool converged = false;
    double relativeStorageChange = 0.0;
    double atmosphericStorage = 0.0;
    double soilStorage = 0.0;
};

float saturationVapourPressureHpa(float temperatureC);
float surfacePressureHpa(float elevationAboveSeaLevelMetres);
float saturationSpecificHumidity(float temperatureC, float pressureHpa);
float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres);
float bulkAerodynamicEvaporationMm(
    float temperatureC,
    float elevationAboveSeaLevelMetres,
    float windSpeedMetresPerSecond,
    float relativeHumidity,
    float timeStepSeconds,
    float transferCoefficient);

void setLastWaterBudget(int season, const WaterBudget& budget);
const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastWaterBudgets();
void setLastHydrologySpinupDiagnostics(const HydrologySpinupDiagnostics& diagnostics);
const HydrologySpinupDiagnostics& lastHydrologySpinupDiagnostics();
}
