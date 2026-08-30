#pragma once

namespace climateatmosphere
{
struct CellSpacing
{
    float zonalMetres = 0.0f;
    float meridionalMetres = 0.0f;
};

struct HorizontalWind
{
    float eastMetresPerSecond = 0.0f;
    float southMetresPerSecond = 0.0f;
};

float coriolisParameterPerSecond(
    float latitudeDegrees,
    float rotationRatePerSecond,
    float rotationDirection = 1.0f);
float hypsometricHeightResponseMetresPerKelvin(
    float lowerPressurePa,
    float upperPressurePa,
    float dryAirGasConstant = 287.05f,
    float gravityMetresPerSecondSquared = 9.80665f);
CellSpacing cellSpacingMetres(
    float latitudeDegrees,
    int longitudeCells,
    int latitudeCells,
    float planetRadiusMetres);
HorizontalWind steadyRayleighCoriolisWind(
    float forceEastMetresPerSecondSquared,
    float forceNorthMetresPerSecondSquared,
    float latitudeDegrees,
    float dragTimeSeconds,
    float rotationRatePerSecond,
    float rotationDirection = 1.0f);
}
