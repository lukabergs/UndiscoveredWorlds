#include "climate_physics.hpp"

#include <cmath>
#include <iostream>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}
}

int main()
{
    expect(
        climatephysics::saturationVapourPressureHpa(30.0f) >
            climatephysics::saturationVapourPressureHpa(20.0f),
        "saturation vapour pressure must increase with temperature");
    expect(
        climatephysics::saturationVapourPressureHpa(20.0f) >
            climatephysics::saturationVapourPressureHpa(0.0f),
        "warm air must support more vapour than freezing air");
    expect(
        climatephysics::surfacePressureHpa(3000.0f) <
            climatephysics::surfacePressureHpa(0.0f),
        "surface pressure must decrease with elevation");

    const float tropicalColumn = climatephysics::saturationColumnWater(30.0f, 0.0f);
    const float freezingColumn = climatephysics::saturationColumnWater(0.0f, 0.0f);
    expect(tropicalColumn > 60.0f && tropicalColumn < 90.0f, "tropical column water is outside its physical sanity range");
    expect(freezingColumn > 8.0f && freezingColumn < 14.0f, "freezing column water is outside its physical sanity range");

    climatephysics::WaterBudget closed;
    closed.oceanEvaporation = 100.0;
    closed.oceanPrecipitation = 40.0;
    closed.runoff = 20.0;
    closed.atmosphericStorage = 15.0;
    closed.soilStorage = 25.0;
    expect(std::abs(closed.residual()) < 1e-9, "closed water budget must have zero residual");

    closed.soilStorage = 20.0;
    expect(std::abs(closed.relativeResidual() - 0.05) < 1e-9, "water-budget loss must be detectable");

    return failures == 0 ? 0 : 1;
}
