#pragma once

#include <array>

#include "physical_layers.hpp"

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

struct CirculationPrecisionStageDiagnostics
{
    double areaWeightedFloatPressureGradientRmsPaPerMetre = 0.0;
    double areaWeightedRoundedPressureGradientRmsPaPerMetre = 0.0;
    double areaWeightedPressureGradientDifferenceRmsPaPerMetre = 0.0;
    double areaWeightedFloatSurfaceWindRmsMetresPerSecond = 0.0;
    double areaWeightedRoundedSurfaceWindRmsMetresPerSecond = 0.0;
    double areaWeightedSurfaceWindDifferenceRmsMetresPerSecond = 0.0;
    double areaWeightedFloatUpperHeightGradientRms = 0.0;
    double areaWeightedRoundedUpperHeightGradientRms = 0.0;
    double areaWeightedUpperHeightGradientDifferenceRms = 0.0;
    double areaWeightedFloatUpperWindRmsMetresPerSecond = 0.0;
    double areaWeightedRoundedUpperWindRmsMetresPerSecond = 0.0;
    double areaWeightedUpperWindDifferenceRmsMetresPerSecond = 0.0;
};

struct CirculationPrecisionDiagnostics
{
    CirculationPrecisionStageDiagnostics base;
    CirculationPrecisionStageDiagnostics final;
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
float heldHouHadleyEdgeLatitudeDegrees(
    float equatorToPoleTemperatureContrastK,
    float troposphereHeightMetres,
    float referenceTemperatureK,
    float gravityMetresPerSecondSquared,
    float rotationRatePerSecond,
    float planetRadiusMetres);
float thermalSurfacePressureAnomalyHpa(
    float temperatureAnomalyK,
    float referencePressureHpa,
    float referenceTemperatureK,
    float massRedistributionEfficiency);
float axisymmetricOverturningPressureAnomalyHpa(
    float latitudeDegrees,
    float thermalEquatorLatitudeDegrees,
    float hadleyHalfWidthDegrees,
    float pressureAmplitudeHpa);
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
void setLastCirculationPrecisionDiagnostics(
    int season,
    const CirculationPrecisionDiagnostics& diagnostics);
const std::array<CirculationPrecisionDiagnostics, CLIMATESEASONCOUNT>&
lastCirculationPrecisionDiagnostics();
}
