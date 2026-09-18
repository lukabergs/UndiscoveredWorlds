#pragma once

#include <cstdint>
#include <array>
#include <vector>

namespace climateocean
{
enum class CoastalScheme { Legacy = 0, FaceDrag = 1, FaceForcing = 2 };

struct OceanForcing
{
    std::vector<std::uint8_t> landMask;
    std::vector<float> bathymetryMetres;
    std::vector<float> eastWindMps;
    std::vector<float> southWindMps;
    std::vector<float> atmosphericTemperatureC;
    std::vector<float> initialSstC;
    // Optional prescribed seasonal ice reservoir; no ice if omitted.
    std::vector<float> initialIceThicknessMetres;
    // Optional diagnosed net surface radiative/sensible/latent exchange.
    // When present it replaces the Newtonian surface heat-exchange proxy.
    std::vector<float> surfaceHeatFluxWm2;
    // Temperature at which the diagnosed flux was evaluated (otherwise air T).
    std::vector<float> surfaceHeatFluxReferenceTemperatureC;
    // Fixed deep reservoir when continuing the surface state across seasons.
    std::vector<float> deepWaterTemperatureC;
    // Optional carried thermodynamic column; distinct from the fixed deep
    // boundary supplying the circulation's vertical continuity exchange.
    std::vector<float> initialStorageDepthMetres, initialReservoirTemperatureC;
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
    // Linear bottom stress divided by column mass adds a damping rate r_b/H.
    float linearBottomDragMps = 0.0f;
    // Use the transport face's depth for drag and wind acceleration. Keep the
    // corner-depth closure and drag-only intermediate for benchmark replay.
    CoastalScheme coastalScheme = CoastalScheme::FaceForcing;
    float minimumCoriolisPerSecond = 2.0e-5f;
    float mixedLayerDepthMetres = 60.0f;
    // Experimental heat storage only. Momentum/face heat transports retain
    // mixedLayerDepthMetres, so changing capacity does not rescale currents.
    bool variableHeatStorage = false;
    float minimumStorageDepthMetres = 10.0f, maximumStorageDepthMetres = 200.0f;
    float storageColumnDepthMetres = 300.0f;
    float storageStratificationPerSecond2 = 2.0e-5f;
    float storageWindMixingEfficiency = 0.2f;
    float storageMixingMemoryDays = 30.0f, storageAdjustmentDays = 10.0f;
    float heatDiffusivityM2S = 750.0f;
    float surfaceHeatExchangeWm2K = 18.0f;
    float waterHeatCapacityJkgK = 3990.0f;
    float freezingTemperatureC = -1.8f;
    float iceDensityKgM3 = 917.0f;
    float latentHeatFusionJkg = 334000.0f;
    float iceConductivityWmK = 2.0f;
    float oceanTimeStepSeconds = 86400.0f;
    int streamfunctionIterations = 800;
    int heatStepsPerIteration = 30;
    int couplingIterations = 40;
    float underRelaxation = 0.35f;
    float convergenceTolerance = 1.0e-3f;
    float sstWindFeedbackMpsPerK = 0.08f; // pressure response hPa/K after the 10 m/s reference conversion
    float deepWaterTemperatureContrastK = 4.0f;
    float streamfunctionTolerance = 1.0e-4f;
    float maximumCurrentMps = 2.5f;
    bool oneWay = false;
    // Preserve explicit integration for numerical reference/replay. The implicit
    // path removes only zonal diffusion from the explicit stability bound.
    bool implicitZonalHeatDiffusion = true;
    // Used only within one periodic solve; no reuse across changed forcing.
    bool reuseSeasonalCirculation = true;
    bool krylovCirculation = true; // False retains the historical point relaxation.
    bool parallelZonalHeatDiffusion = true; // Independent rows; enabled at >=256 columns.
    bool gatherHeatFluxes = true; // False retains the serial face-scatter reference.
    int heatWorkers = 0; // Zero selects up to sixteen workers on grids >=256 wide; one is serial.
};

struct OceanHeatBudget
{
    double depthMetres = 0.0, deepTemperatureC = 0.0;
    // Finite-interval integration before outer under-relaxation; positive warms.
    double horizontalWm2 = 0.0, surfaceWm2 = 0.0, verticalWm2 = 0.0, storageWm2 = 0.0;
    // Internal transfer into the mixed layer; opposite change in the finite
    // subsurface reservoir. Whole-column storage excludes this internal flux.
    double entrainmentWm2 = 0.0, columnStorageWm2 = 0.0;
};

struct OceanStabilityLimit
{
    int cell = -1;
    double totalRate = 0.0, zonalAdvectionRate = 0.0, meridionalAdvectionRate = 0.0;
    double zonalDiffusionRate = 0.0, meridionalDiffusionRate = 0.0, verticalRate = 0.0, surfaceRate = 0.0;
    double depthMetres = 0.0;
    // E, W, S, N faces; velocities point east on E/W and south on S/N.
    std::array<double, 4> faceVelocityMps{}, barotropicFaceVelocityMps{}, faceDepthMetres{};
};

struct OceanStabilityDistribution
{
    int wetCells = 0;
    double wetAreaM2 = 0.0, maximumRateWithoutAdvection = 0.0;
    // Cells/area whose retained explicit rate forbids a step of 1, 3, 6, 12, 24 hours.
    std::array<int, 5> restrictedCells{};
    std::array<double, 5> restrictedAreaM2{};
    std::array<double, 3> rateQuantiles{}; // Unweighted wet-cell p50, p95, p99.
    std::array<OceanStabilityLimit, 8> highestRates{};
};

struct OceanState
{
    std::vector<float> storageDepthMetres, meanStorageDepthMetres, reservoirTemperatureC;
    // Corner streamfunction [columns * (rows + 1)], volume transport m3/s.
    std::vector<float> streamfunctionM3S;
    // East/south face volume transports [columns * rows], m3/s. Exactly closed
    // around land and non-divergent, including spatially varying bathymetry.
    std::vector<double> eastVolumeTransportM3S;
    std::vector<double> southVolumeTransportM3S;
    std::vector<float> eastCurrentMps;
    std::vector<float> southCurrentMps;
    std::vector<float> ekmanUpwellingMps;
    // Components use the same wet-face averaging as the total current.
    // Upwelling is mixed-layer divergence, not full-column mass divergence.
    struct TransportComponents
    {
        double barotropicEastMps = 0.0, barotropicSouthMps = 0.0;
        double ekmanEastMps = 0.0, ekmanSouthMps = 0.0;
        double barotropicUpwellingMps = 0.0, ekmanUpwellingMps = 0.0;
    };
    std::vector<TransportComponents> transportComponents;
    std::vector<float> sstC;
    std::vector<float> surfaceSkinTemperatureC;
    std::vector<float> iceThicknessMetres;
    // Time means over the finite integration interval; end states above are
    // retained separately for continuation and enthalpy accounting.
    std::vector<float> meanSstC, meanSurfaceSkinTemperatureC, meanIceCover;
    std::vector<OceanHeatBudget> heatBudget;
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
    int heatSubsteps = 0, explicitHeatSubsteps = 0;
    OceanStabilityLimit stabilityLimit;
    OceanStabilityDistribution stabilityDistribution;
    double actualHeatStepSeconds = 0.0;
    int maximumZonalDiffusionIterations = 0;
    double maximumZonalDiffusionResidualK = 0.0;
    double circulationSeconds = 0.0, heatSeconds = 0.0, zonalHeatSeconds = 0.0;
    double heatFluxSeconds = 0.0, heatUpdateSeconds = 0.0;
    int heatWorkersUsed = 1;
    int circulationSolves = 0;
    bool converged = false;
    bool finite = true;
};

bool usableOceanState(const OceanState& state, std::size_t cellCount);

OceanState solveWindDrivenOcean(
    int columns,
    int rows,
    const OceanForcing& forcing,
    const OceanConfig& config);

struct PeriodicOceanState
{
    std::array<OceanState, 4> seasons;
    int years = 0;
    double annualEnthalpyDriftK = 0.0;
    double circulationSeconds = 0.0, heatSeconds = 0.0, zonalHeatSeconds = 0.0;
    double heatFluxSeconds = 0.0, heatUpdateSeconds = 0.0;
    int circulationSolves = 0;
    std::uint64_t heatSubsteps = 0;
    bool accepted = false, converged = false;
};
// Continue end-of-season liquid/ice energy, then compare matching year-end
// states. Forcing[0] supplies the initial reservoir; seasonal intervals are
// supplied explicitly in the configs. Only one-way ocean solves are used.
// Bathymetry, land mask and enthalpy constants must match across all seasons.
PeriodicOceanState solvePeriodicOcean(int columns, int rows,
    std::array<OceanForcing, 4> forcing, std::array<OceanConfig, 4> configs,
    int maximumYears = 24, double toleranceK = 0.05);
}
