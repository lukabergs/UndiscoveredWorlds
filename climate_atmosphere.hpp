#pragma once

#include <array>
#include <vector>

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

struct StationaryWaveResponse
{
    std::vector<float> pressureAnomalyHpa;
    std::vector<float> residualHistory;
    int iterations = 0;
    int restartCycles = 0;
    float relativeResidual = 0.0f;
    bool converged = false;
    bool stagnated = false;
};

struct DiabaticHeatingBudget
{
    std::vector<float> netColumnHeatingWm2;
    std::vector<float> stationaryProjectedHeatingWm2;
    double areaWeightedRadiativeHeatingWm2 = 0.0;
    double areaWeightedSensibleHeatingWm2 = 0.0;
    double areaWeightedLatentHeatingWm2 = 0.0;
    double areaWeightedProjectedHeatingWm2 = 0.0;
    float maximumAbsoluteRowMeanWm2 = 0.0f;
};

struct StationaryParameterDiagnosis
{
    float equivalentDepthMetres = 0.0f;
    float gravityWaveSpeedMps = 0.0f;
    float adjustmentLengthMetres = 0.0f;
    float dampingTimeSeconds = 0.0f;
    int maximumZonalWavenumber = 1;
    int maximumMeridionalWavenumber = 1;
};

struct CirculationModeSwitches
{
    bool zonal = true;
    bool stationary = true;
    bool surface = true;
    bool upper = true;
};

struct ModeSeparatedCirculationConfig
{
    float surfaceHeatingPressureResponseHpaPerWm2 = -0.04f;
    float upperHeatingHeightResponseMetresPerWm2 = 1.5f;
    float surfaceEquivalentPressureDepthHpa = 48.0f;
    float upperEquivalentPressureDepthHpa = 24.0f;
    float surfaceDampingTimeSeconds = 190080.0f;
    float upperDampingTimeSeconds = 345600.0f;
    float surfaceDragTimeSeconds = 43200.0f;
    float upperDragTimeSeconds = 172800.0f;
    float surfaceDragCoefficient = 0.0013f;
    float surfaceBoundaryLayerDepthMetres = 300.0f;
    float surfaceToUpperCoupling = 0.10f;
    float interlayerMomentumCoupling = 0.05f;
    // Upper height is the 500 hPa hydrostatic mode, not a second surface pressure.
    std::vector<float> zonalUpperHeightMetres;
    std::vector<float> upperOrographicHeightMetres;
    std::vector<float> surfaceDragTimesSeconds;
    int maximumZonalWavenumber = 8;
    int maximumMeridionalWavenumber = 24;
    int solverRestartLength = 60;
    float airDensityKgM3 = 1.225f;
    float gravityMetresPerSecondSquared = 9.80665f;
    float planetRadiusMetres = 6371000.0f;
    float rotationRatePerSecond = 7.2921159e-5f;
    float rotationDirection = 1.0f;
    int maximumIterations = 1000;
    float relativeTolerance = 1.0e-4f;
    CirculationModeSwitches enabled;
};

struct ModeSeparatedCirculation
{
    std::vector<float> surfacePressureAnomalyHpa;
    std::vector<float> upperHeightAnomalyMetres;
    std::vector<float> surfaceEastWindMps;
    std::vector<float> surfaceSouthWindMps;
    std::vector<float> upperEastWindMps;
    std::vector<float> upperSouthWindMps;
    StationaryWaveResponse surfaceStationarySolver;
    StationaryWaveResponse upperStationarySolver;
    double areaWeightedKineticEnergyJm2 = 0.0;
    double areaWeightedDragDissipationWm2 = 0.0;
    double areaWeightedMassAnomalyKgM2 = 0.0;
    double maximumMomentumExchangeResidual = 0.0;
    float maximumStationaryRowMeanHpa = 0.0f;
    std::vector<float> ascentHpaPerDay;
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
float thermalModePressureAnomalyHpa(
    float temperatureAnomalyK,
    float airDensityKgM3,
    float lowerPressurePa,
    float upperPressurePa,
    float dryAirGasConstant = 287.05f);
float axisymmetricOverturningPressureAnomalyHpa(
    float latitudeDegrees,
    float thermalEquatorLatitudeDegrees,
    float hadleyHalfWidthDegrees,
    float pressureAmplitudeHpa);
DiabaticHeatingBudget diagnoseDiabaticHeating(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& absorbedShortwaveWm2,
    const std::vector<float>& outgoingLongwaveWm2,
    const std::vector<float>& sensibleHeatingWm2,
    const std::vector<float>& condensationMm,
    float accumulationSeconds,
    float verticalProjection,
    float latentHeatCouplingFraction = 1.0f,
    float latentHeatJkg = 2.5e6f,
    bool removeZonalMean = true);
StationaryParameterDiagnosis diagnoseStationaryParameters(
    float bruntVaisalaFrequencyPerSecond,
    float modeDepthMetres,
    float gravityMetresPerSecondSquared,
    float planetRadiusMetres,
    float rotationRatePerSecond,
    int longitudeCells,
    int latitudeCells,
    float nondimensionalDamping,
    float resolvedForcingScaleMetres);
float diagnoseBruntVaisalaFrequency(
    float temperatureK, float lapseRateKPerMetre, float gravityMps2);
// Linear rotating mountain-wave projection with radiative vertical propagation,
// evanescence and critical-level damping. Zonal terrain modes only (metres).
std::vector<float> upperOrographicHeightForcing(
    int columns, int rows, const std::vector<float>& terrainMetres,
    const std::vector<float>& backgroundEastWindMps,
    float stabilityPerSecond, float levelHeightMetres,
    float dampingTimeSeconds, const ModeSeparatedCirculationConfig& config);
void diagnoseModeWinds(
    int columns, int rows, const ModeSeparatedCirculationConfig& config,
    ModeSeparatedCirculation& circulation);
std::vector<float> nonlocalThermalResponse(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& localTemperatureAnomalyK,
    float tropicalLatitudeDegrees,
    float longitudinalReachDegrees,
    float meridionalReachDegrees,
    float rotationDirection = 1.0f);
std::vector<float> mechanicalTopographicPressureForcingHpa(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& terrainMetres,
    const std::vector<float>& eastWindMps,
    const std::vector<float>& southWindMps,
    float sampleDistanceCells,
    float terrainScaleMetres,
    float maximumAmplitudeHpa,
    float minimumWindMps,
    float fullStrengthWindMps,
    float minimumLatitudeDegrees,
    float fullStrengthLatitudeDegrees);
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
HorizontalWind steadyQuadraticDragCoriolisWind(
    float forceEastMetresPerSecondSquared,
    float forceNorthMetresPerSecondSquared,
    float latitudeDegrees,
    float dragCoefficient,
    float boundaryLayerDepthMetres,
    float rotationRatePerSecond,
    float rotationDirection = 1.0f);
StationaryWaveResponse solveSteadyStationaryWavePressure(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& equilibriumPressureAnomalyHpa,
    const std::vector<float>& dragTimeSeconds,
    float equivalentPressureDepthHpa,
    float pressureDampingTimeSeconds,
    float airDensityKgM3,
    float planetRadiusMetres,
    float rotationRatePerSecond,
    float rotationDirection,
    bool preserveZonalMean,
    int maximumIterations,
    float relativeTolerance,
    int restartLength = 60);
ModeSeparatedCirculation solveModeSeparatedCirculation(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& zonalSurfacePressureHpa,
    const std::vector<float>& stationaryHeatingWm2,
    const std::vector<float>& orographicSurfaceForcingHpa,
    const ModeSeparatedCirculationConfig& config);
void setLastCirculationPrecisionDiagnostics(
    int season,
    const CirculationPrecisionDiagnostics& diagnostics);
const std::array<CirculationPrecisionDiagnostics, CLIMATESEASONCOUNT>&
lastCirculationPrecisionDiagnostics();
}
