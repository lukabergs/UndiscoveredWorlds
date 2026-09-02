#pragma once

#include <cstdint>
#include <vector>

namespace climateocean
{
struct OceanForcing
{
    std::vector<std::uint8_t> landMask;
    std::vector<float> bathymetryMetres;
    std::vector<float> eastWindMps;
    std::vector<float> southWindMps;
    std::vector<float> atmosphericTemperatureC;
    std::vector<float> initialSstC;
    // Optional diagnosed net surface radiative/sensible/latent exchange.
    // When present it replaces the Newtonian surface heat-exchange proxy.
    std::vector<float> surfaceHeatFluxWm2;
};

struct OceanConfig
{
    float planetRadiusMetres = 6371000.0f;
    float rotationRatePerSecond = 7.2921159e-5f;
    float rotationDirection = 1.0f;
    float airDensityKgM3 = 1.225f;
    float waterDensityKgM3 = 1025.0f;
    float dragCoefficient = 0.0013f;
    float barotropicDragPerSecond = 1.5e-6f;
    float minimumCoriolisPerSecond = 2.0e-5f;
    float mixedLayerDepthMetres = 60.0f;
    float heatDiffusivityM2S = 750.0f;
    float surfaceHeatExchangeWm2K = 18.0f;
    float waterHeatCapacityJkgK = 3990.0f;
    float oceanTimeStepSeconds = 86400.0f;
    int streamfunctionIterations = 800;
    int heatStepsPerIteration = 30;
    int couplingIterations = 40;
    float underRelaxation = 0.35f;
    float convergenceTolerance = 1.0e-3f;
    float sstWindFeedbackMpsPerK = 0.08f; // pressure response hPa/K after the 10 m/s reference conversion
    float deepWaterTemperatureContrastK = 8.0f;
    float streamfunctionTolerance = 1.0e-4f;
    float maximumCurrentMps = 2.5f;
    bool oneWay = false;
};

struct OceanState
{
    // Corner streamfunction [columns * (rows + 1)], volume transport m3/s.
    std::vector<float> streamfunctionM3S;
    // East/south face volume transports [columns * rows], m3/s. Exactly closed
    // around land and non-divergent, including spatially varying bathymetry.
    std::vector<double> eastVolumeTransportM3S;
    std::vector<double> southVolumeTransportM3S;
    std::vector<float> eastCurrentMps;
    std::vector<float> southCurrentMps;
    std::vector<float> ekmanUpwellingMps;
    std::vector<float> sstC;
    std::vector<float> coupledEastWindMps;
    std::vector<float> coupledSouthWindMps;
    std::vector<float> coupledPressureAnomalyHpa;
    int couplingIterations = 0;
    float relativeResidual = 0.0f;
    float streamfunctionRelativeResidual = 0.0f;
    std::vector<float> residualHistory;
    double maximumTransportDivergenceMps = 0.0;
    double heatBudgetResidualJ = 0.0;
    double relativeHeatBudgetResidual = 0.0;
    bool converged = false;
    bool finite = true;
};

bool usableOceanState(const OceanState& state, std::size_t cellCount);

OceanState solveWindDrivenOcean(
    int columns,
    int rows,
    const OceanForcing& forcing,
    const OceanConfig& config);
}
