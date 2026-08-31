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
    expect(
        climatephysics::surfacePressureHpa(3000.0f, 2.0f) <
            climatephysics::surfacePressureHpa(3000.0f, 1.0f),
        "stronger gravity must reduce pressure more rapidly with elevation");

    const float tropicalColumn = climatephysics::saturationColumnWater(30.0f, 0.0f);
    const float freezingColumn = climatephysics::saturationColumnWater(0.0f, 0.0f);
    expect(tropicalColumn > 60.0f && tropicalColumn < 90.0f, "tropical column water is outside its physical sanity range");
    expect(freezingColumn > 8.0f && freezingColumn < 14.0f, "freezing column water is outside its physical sanity range");
    expect(
        climatephysics::saturationColumnWaterAtPressure(20.0f, 980.0f) >
            climatephysics::saturationColumnWaterAtPressure(20.0f, 1030.0f),
        "column saturation capacity must respond to surface pressure");
    expect(
        climatephysics::saturationColumnWaterAtPressure(20.0f, 1013.25f, 2.0f) <
            climatephysics::saturationColumnWaterAtPressure(20.0f, 1013.25f, 1.0f),
        "stronger gravity must reduce column water capacity");

    const float oceanDrag = climatephysics::neutralDragCoefficient(0.0002f, 10.0f);
    const float landDrag = climatephysics::neutralDragCoefficient(0.05f, 10.0f);
    expect(landDrag > oceanDrag,
        "rough land must exchange more strongly than a smooth ocean surface");
    expect(
        climatephysics::bulkRichardsonExchangeMultiplier(
            10.0f, 14.0f, 3.0f, 1.0f, 10.0f, 0.25f, 2.0f) < 1.0f,
        "stable surface layers must suppress turbulent exchange");
    expect(
        climatephysics::bulkRichardsonExchangeMultiplier(
            18.0f, 14.0f, 3.0f, 1.0f, 10.0f, 0.25f, 2.0f) > 1.0f,
        "unstable surface layers must enhance turbulent exchange");

    const float dryEvaporation = climatephysics::bulkAerodynamicEvaporationMm(
        20.0f, 0.0f, 5.0f, 0.0f, 86400.0f, 0.0013f);
    const float humidEvaporation = climatephysics::bulkAerodynamicEvaporationMm(
        20.0f, 0.0f, 5.0f, 0.8f, 86400.0f, 0.0013f);
    expect(dryEvaporation > humidEvaporation, "bulk evaporation must fall as relative humidity rises");
    expect(humidEvaporation > 1.0f && humidEvaporation < 5.0f, "humid daily ocean evaporation is outside its physical sanity range");
    expect(
        climatephysics::bulkAerodynamicEvaporationMm(20.0f, 0.0f, 5.0f, 1.0f, 86400.0f, 0.0013f) == 0.0f,
        "saturated air must suppress evaporation");

    const float equatorialDiffusivity = climatephysics::transientEddyDiffusivityM2S(
        0.0f, 0.0f, 20.0f, 0.0f, 0.0f, 25000.0f, 250000.0f, 20.0f, 45.0f);
    const float midlatitudeDiffusivity = climatephysics::transientEddyDiffusivityM2S(
        0.0f, 0.0f, 20.0f, 0.0f, 45.0f, 25000.0f, 250000.0f, 20.0f, 45.0f);
    const float subtropicalDiffusivity = climatephysics::transientEddyDiffusivityM2S(
        0.0f, 0.0f, 20.0f, 0.0f, 15.0f, 25000.0f, 250000.0f, 20.0f, 45.0f);
    expect(equatorialDiffusivity == 0.0f, "rotational eddy closure must vanish at the equator");
    expect(midlatitudeDiffusivity > 0.0f, "vertical shear must generate midlatitude eddy mixing");
    expect(subtropicalDiffusivity == 0.0f, "baroclinic eddy mixing must stay outside the tropical belt");
    expect(midlatitudeDiffusivity <= 250000.0f, "eddy diffusivity must respect its physical cap");
    expect(
        climatephysics::transientEddyDiffusivityM2S(
            5.0f, -2.0f, 5.0f, -2.0f, 45.0f, 25000.0f, 250000.0f, 20.0f, 45.0f) == 0.0f,
        "eddy mixing must vanish without vertical wind shear");

    const float belowCriticalCondensation = climatephysics::relaxedExcessCondensationAmount(
        10.0f, 12.0f, 86400.0f, 172800.0f);
    const float excessCondensation = climatephysics::relaxedExcessCondensationAmount(
        20.0f, 10.0f, 86400.0f, 172800.0f);
    expect(
        belowCriticalCondensation == 0.0f,
        "sub-critical atmospheric water must not produce unconditional desert drizzle");
    expect(
        std::abs(excessCondensation - 3.93469f) < 0.0001f,
        "excess atmospheric water must relax on the configured two-day e-folding timescale");

    climatephysics::WaterBudget closed;
    closed.oceanEvaporation = 100.0;
    closed.oceanPrecipitation = 40.0;
    closed.runoff = 20.0;
    closed.atmosphericStorage = 15.0;
    closed.soilStorage = 25.0;
    expect(std::abs(closed.residual()) < 1e-9, "closed water budget must have zero residual");

    closed.soilStorage = 20.0;
    expect(std::abs(closed.relativeResidual() - 0.05) < 1e-9, "water-budget loss must be detectable");

    climatephysics::WaterBudget persistent;
    persistent.initialAtmosphericStorage = 40.0;
    persistent.initialSoilStorage = 60.0;
    persistent.initialSnowStorage = 10.0;
    persistent.oceanEvaporation = 25.0;
    persistent.oceanPrecipitation = 45.0;
    persistent.runoff = 20.0;
    persistent.atmosphericStorage = 35.0;
    persistent.soilStorage = 25.0;
    persistent.snowStorage = 10.0;
    expect(std::abs(persistent.residual()) < 1e-9,
        "persistent water budget must include initial atmospheric and soil storage");

    return failures == 0 ? 0 : 1;
}
