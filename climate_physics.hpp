#pragma once

#include <array>

#include "physical_layers.hpp"

namespace climatephysics
{
struct WaterBudget
{
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

float saturationVapourPressureHpa(float temperatureC);
float surfacePressureHpa(float elevationAboveSeaLevelMetres);
float saturationColumnWater(float temperatureC, float elevationAboveSeaLevelMetres);

void setLastWaterBudget(int season, const WaterBudget& budget);
const std::array<WaterBudget, CLIMATESEASONCOUNT>& lastWaterBudgets();
}
