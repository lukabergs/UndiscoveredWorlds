#include "climate_ocean_dynamics.hpp"
#include "detail/ocean_thermal_remap.hpp"
#include "detail/ocean_replay.hpp"
#include "detail/ekman_cell_mobility.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        failures++;
    }
}

struct FaceWork { double wind = 0.0, drag = 0.0, maximumSpeed = 0.0; };
FaceWork faceWork(int columns, int rows, const climateocean::OceanForcing& f,
    const climateocean::OceanConfig& config, const climateocean::OceanState& state)
{
    FaceWork result;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    const double dy = config.planetRadiusMetres * grid.latitudeSpacingRadians;
    const auto wet = [&](std::size_t c) { return !f.landMask[c] && f.bathymetryMetres[c] > 0; };
    const auto stress = [&](std::size_t c, bool east)
    {
        return config.airDensityKgM3 * config.dragCoefficient *
            std::hypot(static_cast<double>(f.eastWindMps[c]), f.southWindMps[c]) *
            (east ? f.eastWindMps[c] : f.southWindMps[c]);
    };
    for (int y = 0; y < rows; ++y) for (int x = 0; x < columns; ++x)
    {
        const auto c = grid.index(x, y);
        if (!wet(c)) continue;
        const auto face = [&](std::size_t neighbour, double volume, double length, double area, bool east)
        {
            if (!wet(neighbour)) return;
            const double h = 0.5 * (std::max(50.0f, f.bathymetryMetres[c]) + std::max(50.0f, f.bathymetryMetres[neighbour]));
            const double velocity = volume / (h * length);
            result.maximumSpeed = std::max(result.maximumSpeed, std::abs(velocity));
            result.wind += 0.5 * (stress(c, east) + stress(neighbour, east)) * velocity * area;
            result.drag += config.waterDensityKgM3 * (config.barotropicDragPerSecond * h + config.linearBottomDragMps) *
                velocity * velocity * area;
        };
        const double eastArea = config.planetRadiusMetres * grid.longitudeSpacingRadians *
            std::cos(grid.latitudeCentresRadians[y]) * dy;
        face(grid.index(x + 1, y), state.eastVolumeTransportM3S[c], dy, eastArea, true);
        if (y + 1 < rows) face(grid.index(x, y + 1), state.southVolumeTransportM3S[c],
            grid.southFaceLengthsMetres[y], grid.southFaceLengthsMetres[y] * dy, false);
    }
    return result;
}

int replay(const std::string& path, int days, int stepHours)
{
    if (days < 0 || days > 90 || stepHours < 1 || stepHours > 24 || 24 % stepHours)
        throw std::runtime_error("Replay expects 0..90 days and a step dividing 24 hours");
    auto fixture = climateocean::detail::OceanReplay::read(path);
    auto config = fixture.config;
    config.oneWay = true;
    config.oceanTimeStepSeconds = days ? stepHours * 3600.0f : 1.0f;
    config.heatStepsPerIteration = days ? days * 24 / stepHours : 1;
    std::array<climateocean::OceanState, 3> states;
    std::cout << std::setprecision(12) << "{\"columns\":" << fixture.columns << ",\"season\":" << fixture.season
        << ",\"days\":" << days << ",\"step_hours\":" << stepHours << ",\"cases\":[";
    for (int scheme = 0; scheme < 3; ++scheme)
    {
        config.coastalScheme = static_cast<climateocean::CoastalScheme>(scheme);
        auto& state = states[scheme];
        state = climateocean::solveWindDrivenOcean(fixture.columns, fixture.rows, fixture.forcing, config);
        if (!climateocean::usableOceanState(state, fixture.forcing.landMask.size())) throw std::runtime_error("Replay solve failed");
        const auto power = faceWork(fixture.columns, fixture.rows, fixture.forcing, config, state);
        const auto& dist = state.stabilityDistribution;
        if (scheme) std::cout << ',';
        std::cout << "{\"scheme\":" << scheme << ",\"circulation_seconds\":" << state.circulationSeconds
            << ",\"heat_seconds\":" << state.heatSeconds << ",\"steps\":" << state.heatSubsteps
            << ",\"residual\":" << state.streamfunctionRelativeResidual << ",\"volume_divergence_mps\":" << state.maximumTransportDivergenceMps
            << ",\"heat_relative_residual\":" << state.relativeHeatBudgetResidual
            << ",\"stable_hours\":" << 0.7 / state.stabilityLimit.totalRate / 3600.0
            << ",\"without_advection_hours\":" << 0.7 / dist.maximumRateWithoutAdvection / 3600.0
            << ",\"maximum_barotropic_face_mps\":" << power.maximumSpeed
            << ",\"wind_power_w\":" << power.wind << ",\"drag_power_w\":" << power.drag
            << ",\"six_hour_cells\":" << dist.restrictedCells[2]
            << ",\"six_hour_area_fraction\":" << dist.restrictedAreaM2[2] / dist.wetAreaM2 << '}';
    }
    std::cout << ']';
    if (days && stepHours > 1)
    {
        config.coastalScheme = climateocean::CoastalScheme::FaceForcing;
        config.oceanTimeStepSeconds = 3600.0f; config.heatStepsPerIteration = days * 24;
        const auto refined = climateocean::solveWindDrivenOcean(fixture.columns, fixture.rows, fixture.forcing, config);
        if (!climateocean::usableOceanState(refined, fixture.forcing.landMask.size())) throw std::runtime_error("Refinement failed");
        double squared = 0.0, area = 0.0, maximum = 0.0, maximumIce = 0.0;
        for (std::size_t c = 0; c < fixture.forcing.landMask.size(); ++c)
        {
            if (fixture.forcing.landMask[c] || fixture.forcing.bathymetryMetres[c] <= 0) continue;
            const double weight = climategrid::latitudeBandMeasure(static_cast<int>(c) / fixture.columns, fixture.rows);
            const double delta = states[2].sstC[c] - refined.sstC[c];
            squared += weight * delta * delta; area += weight; maximum = std::max(maximum, std::abs(delta));
            maximumIce = std::max(maximumIce, static_cast<double>(std::abs(states[2].iceThicknessMetres[c] - refined.iceThicknessMetres[c])));
        }
        std::cout << ",\"face_forcing_vs_one_hour\":{\"sst_area_rmse_k\":" << std::sqrt(squared / area)
            << ",\"maximum_sst_k\":" << maximum << ",\"maximum_ice_m\":" << maximumIce << '}';
    }
    std::cout << "}\n";
    return 0;
}
}

int main(int argc, char** argv)
{
    if (argc == 6 && std::string(argv[1]) == "--momentum-replay")
    {
        auto fixture=climateocean::detail::OceanReplay::read(argv[2]);
        auto config=fixture.config;
        config.minimumCoriolisPerSecond=std::stof(argv[4]);
        config.barotropicDragPerSecond=std::stof(argv[5]);
        config.oneWay=true;config.oceanTimeStepSeconds=1.0f;config.heatStepsPerIteration=1;
        const auto state=climateocean::solveWindDrivenOcean(fixture.columns,fixture.rows,fixture.forcing,config);
        if (!climateocean::usableOceanState(state,fixture.forcing.landMask.size()) || !state.converged)
            throw std::runtime_error("Frozen momentum replay failed");
        std::ofstream csv(argv[3]);
        csv << std::setprecision(12) << "east_mps,north_mps,gyre_east_mps,gyre_north_mps,ekman_east_mps,ekman_north_mps\n";
        for (std::size_t c=0;c<state.eastCurrentMps.size();++c)
        {
            const auto& v=state.transportComponents[c];
            csv << state.eastCurrentMps[c] << ',' << -state.southCurrentMps[c] << ','
                << v.barotropicEastMps << ',' << -v.barotropicSouthMps << ','
                << v.ekmanEastMps << ',' << -v.ekmanSouthMps << '\n';
        }
        if (!csv) throw std::runtime_error("Could not export frozen momentum replay");
        const auto power=faceWork(fixture.columns,fixture.rows,fixture.forcing,config,state);
        std::cout << std::setprecision(12) << "{\"season\":" << fixture.season
            << ",\"ekman_damping_per_s\":" << config.minimumCoriolisPerSecond
            << ",\"gyre_damping_per_s\":" << config.barotropicDragPerSecond
            << ",\"volume_divergence_mps\":" << state.maximumTransportDivergenceMps
            << ",\"heat_relative_residual\":" << state.relativeHeatBudgetResidual
            << ",\"momentum_residual\":" << state.streamfunctionRelativeResidual
            << ",\"wind_power_w\":" << power.wind << ",\"drag_power_w\":" << power.drag << "}\n";
        return 0;
    }
    if (argc >= 4 && std::string(argv[1]) == "--coastal-replay")
        return replay(argv[2], std::stoi(argv[3]), argc > 4 ? std::stoi(argv[4]) : 24);
    {
        constexpr double r=2e-6;
        for (const auto bounds : {std::pair<double,double>{-8e-6,0.0}, {0.0,8e-6},
                 {-8e-6,8e-6}, {1.0e-4,1.1e-4}, {-1.1e-4,-1.0e-4}})
        {
            const double lo=bounds.first,hi=bounds.second;
            const auto mobility=climateocean::detail::ekmanCellMobility(r,hi,lo);
            double along=0.0,cross=0.0;
            constexpr int samples=10000;
            for (int i=0;i<samples;++i)
            {
                const double f=lo+(i+.5)*(hi-lo)/samples;
                along+=r/(r*r+f*f)/samples;cross+=f/(r*r+f*f)/samples;
            }
            expect(mobility.alongStressSeconds>0.0 &&
                std::abs(mobility.alongStressSeconds-along)<1e-6*along &&
                std::abs(mobility.crossStressSeconds-cross)<1e-6*std::max(1.0,std::abs(cross)),
                "cell-averaged Ekman response must match independent quadrature and dissipate wind work");
            const auto half1=climateocean::detail::ekmanCellMobility(r,hi,.5*(hi+lo));
            const auto half2=climateocean::detail::ekmanCellMobility(r,.5*(hi+lo),lo);
            expect(std::abs(mobility.alongStressSeconds-.5*(half1.alongStressSeconds+half2.alongStressSeconds))<1e-7 &&
                std::abs(mobility.crossStressSeconds-.5*(half1.crossStressSeconds+half2.crossStressSeconds))<1e-7,
                "uniform-stress Ekman response must be invariant under area subdivision");
        }
        const auto equator=climateocean::detail::ekmanCellMobility(r,0.0,0.0);
        const auto symmetric=climateocean::detail::ekmanCellMobility(r,8e-6,-8e-6);
        expect(equator.alongStressSeconds==1.0/r && equator.crossStressSeconds==0.0 &&
            symmetric.crossStressSeconds==0.0,
            "nonrotating and symmetric equatorial cells must have no cross-stress response");
    }
    constexpr int columns = 32;
    constexpr int rows = 16;
    constexpr std::size_t cellCount = columns * rows;
    constexpr float pi = 3.14159265358979323846f;
    const auto index = [=](int x, int y)
    {
        return static_cast<std::size_t>(y) * columns + x;
    };
    climateocean::OceanForcing forcing;
    forcing.landMask.assign(cellCount, 0);
    forcing.bathymetryMetres.assign(cellCount, 4000.0f);
    forcing.eastWindMps.assign(cellCount, 0.0f);
    forcing.southWindMps.assign(cellCount, 0.0f);
    forcing.atmosphericTemperatureC.assign(cellCount, 0.0f);
    forcing.initialSstC.assign(cellCount, 0.0f);
    for (int y = 0; y < rows; y++)
    {
        const float latitude = pi * 0.5f - pi * (y + 0.5f) / rows;
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = index(x, y);
            const bool wall = x == 0 || x == columns - 1 || y == 0 || y == rows - 1;
            forcing.landMask[cell] = wall ? 1 : 0;
            forcing.bathymetryMetres[cell] = wall ? 0.0f : 4000.0f;
            forcing.eastWindMps[cell] = 10.0f * std::cos(2.0f * latitude);
            forcing.atmosphericTemperatureC[cell] = 18.0f * std::cos(latitude) - 2.0f;
            forcing.initialSstC[cell] = forcing.atmosphericTemperatureC[cell];
        }
    }

    climateocean::OceanConfig oneWayConfig;
    oneWayConfig.oneWay = true;
    oneWayConfig.streamfunctionIterations = 500;
    oneWayConfig.heatStepsPerIteration = 2;
    const auto oneWay = climateocean::solveWindDrivenOcean(
        columns, rows, forcing, oneWayConfig);
    double currentMagnitude = 0.0;
    double westernBoundaryMagnitude = 0.0;
    double easternBoundaryMagnitude = 0.0;
    bool noNormalFlow = true;
    bool hasUpwelling = false;
    bool hasDownwelling = false;
    for (int y = 1; y < rows - 1; y++)
    {
        westernBoundaryMagnitude += std::hypot(
            oneWay.eastCurrentMps[index(1, y)], oneWay.southCurrentMps[index(1, y)]);
        easternBoundaryMagnitude += std::hypot(
            oneWay.eastCurrentMps[index(columns - 2, y)],
            oneWay.southCurrentMps[index(columns - 2, y)]);
        for (int x = 1; x < columns - 1; x++)
        {
            const std::size_t cell = index(x, y);
            currentMagnitude += std::hypot(
                oneWay.eastCurrentMps[cell], oneWay.southCurrentMps[cell]);
            hasUpwelling = hasUpwelling || oneWay.ekmanUpwellingMps[cell] > 0.0f;
            hasDownwelling = hasDownwelling || oneWay.ekmanUpwellingMps[cell] < 0.0f;
            if (x == 1) noNormalFlow = noNormalFlow && oneWay.eastVolumeTransportM3S[index(x - 1, y)] == 0.0;
            if (x == columns - 2) noNormalFlow = noNormalFlow && oneWay.eastVolumeTransportM3S[cell] == 0.0;
            if (y == 1) noNormalFlow = noNormalFlow && oneWay.southVolumeTransportM3S[index(x, y - 1)] == 0.0;
            if (y == rows - 2) noNormalFlow = noNormalFlow && oneWay.southVolumeTransportM3S[cell] == 0.0;
        }
    }
    expect(oneWay.converged && currentMagnitude > 0.0,
        "wind stress must drive a deterministic basin circulation");
    expect(noNormalFlow,
        "the basin solver must enforce no-normal flow at every coastline");
    expect(westernBoundaryMagnitude > easternBoundaryMagnitude,
        "the beta-plane closure must intensify the western boundary current");
    expect(hasUpwelling && hasDownwelling,
        "wind-stress divergence must diagnose both upwelling and downwelling");
    expect(oneWay.maximumTransportDivergenceMps < 1.0e-14,
        "barotropic face transport must be divergence free to roundoff");

    climateocean::OceanConfig coupledConfig = oneWayConfig;
    coupledConfig.oneWay = false;
    coupledConfig.couplingIterations = 40;
    coupledConfig.convergenceTolerance = 1.0e-4f;
    const auto first = climateocean::solveWindDrivenOcean(
        columns, rows, forcing, coupledConfig);
    const auto second = climateocean::solveWindDrivenOcean(
        columns, rows, forcing, coupledConfig);
    expect(first.converged && first.couplingIterations > 0 &&
            first.eastCurrentMps == second.eastCurrentMps && first.sstC == second.sstC,
        "coupled atmosphere-ocean iterations must converge deterministically");
    expect(first.relativeResidual <= coupledConfig.convergenceTolerance &&
            first.streamfunctionRelativeResidual <= coupledConfig.streamfunctionTolerance &&
            first.residualHistory.back() < first.residualHistory.front(),
        "coupled convergence must use both the fixed-point and physical equation residuals");
    expect(first.coupledEastWindMps != forcing.eastWindMps ||
            first.coupledSouthWindMps != forcing.southWindMps,
        "interactive SST gradients must feed back into the atmospheric wind");
    const double heatScale = 1.0e22;
    expect(std::abs(first.heatBudgetResidualJ) / heatScale < 1.0e-10,
        "mixed-layer heat transport must close its diagnosed energy budget");

    auto uniform = forcing;
    std::fill(uniform.initialSstC.begin(), uniform.initialSstC.end(), 2.0f);
    std::fill(uniform.atmosphericTemperatureC.begin(), uniform.atmosphericTemperatureC.end(), 2.0f);
    oneWayConfig.deepWaterTemperatureContrastK = 0.0f;
    oneWayConfig.surfaceHeatExchangeWm2K = 0.0f;
    for (std::size_t cell = 0; cell < cellCount; ++cell)
        if (!uniform.landMask[cell]) uniform.bathymetryMetres[cell] = 500.0f + (cell % 7) * 400.0f;
    const auto constant = climateocean::solveWindDrivenOcean(columns, rows, uniform, oneWayConfig);
    expect(std::all_of(constant.sstC.begin(), constant.sstC.end(), [](float t) { return std::abs(t - 2.0f) < 1.0e-5f; }) &&
            constant.maximumTransportDivergenceMps < 1.0e-14,
        "variable bathymetry and convergent Ekman flow must preserve uniform temperature with closed vertical exchange");
    // With rotation removed, physical wind work must balance physical drag.
    // Integrate on velocity-face dual areas, independently of the curl stencil.
    auto shelf = forcing;
    for (int y = 0; y < rows; ++y) for (int x = 0; x < columns; ++x)
        if (!shelf.landMask[index(x, y)]) shelf.bathymetryMetres[index(x, y)] = x < 3 ? 50.0f : 4000.0f;
    auto shelfConfig = oneWayConfig;
    shelfConfig.rotationRatePerSecond = 0.0f;
    shelfConfig.linearBottomDragMps = 0.001f;
    shelfConfig.oceanTimeStepSeconds = 1.0f;
    shelfConfig.heatStepsPerIteration = 1;
    shelfConfig.streamfunctionTolerance = 1.0e-7f;
    shelfConfig.streamfunctionIterations = 2000;
    const auto shelfFlow = climateocean::solveWindDrivenOcean(columns, rows, shelf, shelfConfig);
    const auto shelfPower = faceWork(columns, rows, shelf, shelfConfig, shelfFlow);
    expect(shelfFlow.converged && shelfPower.wind > 0 && shelfPower.drag > 0 &&
        std::abs(shelfPower.wind - shelfPower.drag) / shelfPower.drag < 5.0e-6,
        "nonrotating steep-shelf flow must balance face-integrated wind work and drag dissipation");
    expect(shelfFlow.maximumTransportDivergenceMps < 1.0e-14 && std::abs(shelfFlow.relativeHeatBudgetResidual) < 1.0e-8,
        "face-depth coastal closure must conserve full-column volume and mixed-layer heat");
    for (int scheme = 0; scheme < 3; ++scheme)
    {
        auto flatConfig = oneWayConfig;
        flatConfig.coastalScheme = static_cast<climateocean::CoastalScheme>(scheme);
        const auto flat = climateocean::solveWindDrivenOcean(columns, rows, forcing, flatConfig);
        double error = 0.0, scale = 0.0;
        for (std::size_t c = 0; c < cellCount; ++c)
        {
            error = std::max(error, std::abs(flat.eastVolumeTransportM3S[c] - oneWay.eastVolumeTransportM3S[c]));
            scale = std::max(scale, std::abs(oneWay.eastVolumeTransportM3S[c]));
        }
        expect(flat.converged && error / scale < 1.0e-6,
            "coastal discretizations must agree for constant bathymetry");
        for (float rotation : {-1.0f, 1.0f})
        {
            auto rotated = shelfConfig;
            rotated.rotationRatePerSecond = oneWayConfig.rotationRatePerSecond;
            rotated.rotationDirection = rotation;
            rotated.coastalScheme = static_cast<climateocean::CoastalScheme>(scheme);
            const auto flow = climateocean::solveWindDrivenOcean(columns, rows, shelf, rotated);
            expect(flow.converged && flow.maximumTransportDivergenceMps < 1.0e-14,
                "all coastal replay schemes must retain converged closed transport under either rotation");
        }
    }
    const climateocean::detail::OceanReplay capture{columns, rows, 2, shelfConfig, shelf};
    const auto captureDirectory = std::filesystem::temp_directory_path() / "uw-ocean-replay-roundtrip";
    const auto capturePath = capture.write(captureDirectory.string());
    const auto restored = climateocean::detail::OceanReplay::read(capturePath.string());
    const auto restoredFlow = climateocean::solveWindDrivenOcean(columns, rows, restored.forcing, restored.config);
    expect(restored.columns == columns && restored.rows == rows && restored.season == 2 &&
        restored.forcing.bathymetryMetres == shelf.bathymetryMetres && restored.forcing.eastWindMps == shelf.eastWindMps &&
        restored.forcing.deepWaterTemperatureC.empty() && restoredFlow.eastVolumeTransportM3S == shelfFlow.eastVolumeTransportM3S &&
        restoredFlow.sstC == shelfFlow.sstC, "versioned ocean fixture must reproduce exact forcing and solver outputs");
    auto storageCapture = capture;
    storageCapture.config.variableHeatStorage = true;
    storageCapture.forcing.initialStorageDepthMetres.assign(columns * rows, 30.0f);
    storageCapture.forcing.initialReservoirTemperatureC.assign(columns * rows, storageCapture.config.freezingTemperatureC);
    const auto storagePath = storageCapture.write(captureDirectory.string());
    const auto storageRestored = climateocean::detail::OceanReplay::read(storagePath.string());
    const auto storageFlow = climateocean::solveWindDrivenOcean(columns, rows, storageCapture.forcing, storageCapture.config);
    const auto storageReplayFlow = climateocean::solveWindDrivenOcean(columns, rows, storageRestored.forcing, storageRestored.config);
    expect(storageFlow.converged && storageReplayFlow.converged && storageRestored.config.variableHeatStorage &&
        storageFlow.sstC == storageReplayFlow.sstC && storageFlow.storageDepthMetres == storageReplayFlow.storageDepthMetres &&
        storageFlow.reservoirTemperatureC == storageReplayFlow.reservoirTemperatureC,
        "V2 diagnostic replay must retain variable heat storage configuration, depth and subsurface heat");
    std::filesystem::remove(storagePath);
    { std::ofstream invalid(capturePath); invalid << "UW_OCEAN_V1 1000000 500000 0\n"; }
    bool rejected = false;
    try { climateocean::detail::OceanReplay::read(capturePath.string()); } catch (const std::runtime_error&) { rejected = true; }
    expect(rejected, "ocean replay must reject invalid dimensions before allocating fields");
    std::filesystem::remove(capturePath);
    std::filesystem::remove(captureDirectory);
    oneWayConfig.streamfunctionIterations = 1;
    const auto incomplete = climateocean::solveWindDrivenOcean(columns, rows, forcing, oneWayConfig);
    expect(!incomplete.converged && incomplete.streamfunctionRelativeResidual > oneWayConfig.streamfunctionTolerance,
        "one-way mode must not report an unconverged basin solve as converged");
    expect(!climateocean::usableOceanState(incomplete, cellCount) && climateocean::usableOceanState(first, cellCount),
        "production acceptance must reject finite but unconverged ocean solutions");
    uniform.surfaceHeatFluxWm2.assign(cellCount, 25.0f);
    oneWayConfig.streamfunctionIterations = coupledConfig.streamfunctionIterations;
    const auto heated = climateocean::solveWindDrivenOcean(columns, rows, uniform, oneWayConfig);
    expect(heated.sstC != constant.sstC && std::abs(heated.relativeHeatBudgetResidual) < 1.0e-8,
        "diagnosed surface heat exchange must change SST and close a normalized ocean heat budget");
    bool localHeatBudgetsClose = heated.heatBudget.size() == cellCount;
    for (std::size_t cell = 0; cell < cellCount; ++cell)
    {
        if (uniform.landMask[cell]) continue;
        const auto& budget = heated.heatBudget[cell];
        const double sources = budget.horizontalWm2 + budget.surfaceWm2 + budget.verticalWm2;
        localHeatBudgetsClose = localHeatBudgetsClose && std::abs(budget.surfaceWm2 - 25.0) < 1.0e-10 &&
            std::abs(budget.storageWm2 - sources) < 1.0e-8 &&
            budget.depthMetres == uniform.bathymetryMetres[cell];
    }
    expect(localHeatBudgetsClose,
        "each variable-depth cell must expose the applied surface, horizontal and vertical heat budget");
    bool componentsClose = heated.transportComponents.size() == cellCount;
    for (std::size_t cell = 0; cell < cellCount; ++cell)
    {
        const auto& part = heated.transportComponents[cell];
        componentsClose = componentsClose &&
            std::abs(part.barotropicEastMps + part.ekmanEastMps - heated.eastCurrentMps[cell]) < 1.0e-7 &&
            std::abs(part.barotropicSouthMps + part.ekmanSouthMps - heated.southCurrentMps[cell]) < 1.0e-7 &&
            std::abs(part.barotropicUpwellingMps + part.ekmanUpwellingMps - heated.ekmanUpwellingMps[cell]) < 1.0e-10;
    }
    expect(componentsClose, "barotropic and Ekman diagnostics must reconstruct variable-depth currents and mixed-layer divergence");
    auto bottomConfig = oneWayConfig;
    bottomConfig.linearBottomDragMps = 0.001f;
    const auto bottomHeated = climateocean::solveWindDrivenOcean(columns, rows, uniform, bottomConfig);
    expect(bottomHeated.converged && bottomHeated.maximumTransportDivergenceMps < 1.0e-14 &&
        std::abs(bottomHeated.relativeHeatBudgetResidual) < 1.0e-8,
        "depth-dependent bottom stress must retain closed volume transport and heat budget over varying bathymetry");
    auto bottomUniform = uniform;
    bottomUniform.surfaceHeatFluxWm2.clear();
    const auto bottomConstant = climateocean::solveWindDrivenOcean(columns, rows, bottomUniform, bottomConfig);
    expect(bottomConstant.converged && std::all_of(bottomConstant.sstC.begin(), bottomConstant.sstC.end(),
        [](float t) { return std::abs(t - 2.0f) < 1.0e-5f; }),
        "depth-dependent bottom stress must preserve uniform temperature with closed vertical exchange");
    auto constantDepth = bottomUniform;
    std::fill(constantDepth.bathymetryMetres.begin(), constantDepth.bathymetryMetres.end(), 1000.0f);
    const auto bottomFlow = climateocean::solveWindDrivenOcean(columns, rows, constantDepth, bottomConfig);
    auto equivalentConfig = bottomConfig;
    equivalentConfig.barotropicDragPerSecond += bottomConfig.linearBottomDragMps / 1000.0f;
    equivalentConfig.linearBottomDragMps = 0.0f;
    const auto equivalentFlow = climateocean::solveWindDrivenOcean(columns, rows, constantDepth, equivalentConfig);
    double maximumBottomDifference = 0.0, maximumBottomTransport = 0.0;
    for (std::size_t cell = 0; cell < bottomFlow.eastVolumeTransportM3S.size(); ++cell)
    {
        maximumBottomDifference = std::max(maximumBottomDifference, std::abs(
            bottomFlow.eastVolumeTransportM3S[cell] - equivalentFlow.eastVolumeTransportM3S[cell]));
        maximumBottomTransport = std::max(maximumBottomTransport, std::abs(bottomFlow.eastVolumeTransportM3S[cell]));
    }
    expect(bottomFlow.converged && equivalentFlow.converged && maximumBottomDifference / maximumBottomTransport < 1.0e-6,
        "constant-depth bottom stress must equal uniform damping with the same r0 plus r_b/H rate");
    for (float invalid : {-0.001f, std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()})
    {
        bottomConfig.linearBottomDragMps = invalid;
        expect(!climateocean::solveWindDrivenOcean(columns, rows, uniform, bottomConfig).finite,
            "bottom stress must reject negative and non-finite drag");
    }
    auto aquaplanet = forcing;
    std::fill(aquaplanet.landMask.begin(), aquaplanet.landMask.end(), 0);
    std::fill(aquaplanet.bathymetryMetres.begin(), aquaplanet.bathymetryMetres.end(), 4000.0f);
    const auto polarOcean = climateocean::solveWindDrivenOcean(columns, rows, aquaplanet, coupledConfig);
    expect(polarOcean.finite && polarOcean.maximumTransportDivergenceMps < 1.0e-14 &&
        std::all_of(polarOcean.sstC.begin(), polarOcean.sstC.end(), [](float t) { return std::isfinite(t); }),
        "open polar oceans must close faces and clamp meridional pressure neighbours safely");

    // Uniform zero-wind flux fixtures isolate phase change from transport.
    auto frozen = aquaplanet;
    std::fill(frozen.eastWindMps.begin(), frozen.eastWindMps.end(), 0.0f);
    std::fill(frozen.southWindMps.begin(), frozen.southWindMps.end(), 0.0f);
    std::fill(frozen.initialSstC.begin(), frozen.initialSstC.end(), oneWayConfig.freezingTemperatureC);
    std::fill(frozen.atmosphericTemperatureC.begin(), frozen.atmosphericTemperatureC.end(), -20.0f);
    frozen.surfaceHeatFluxWm2.assign(cellCount, -100.0f);
    oneWayConfig.heatStepsPerIteration = 1;
    oneWayConfig.oceanTimeStepSeconds = 86400.0f;
    oneWayConfig.heatDiffusivityM2S = 0.0f;
    const double latent = static_cast<double>(oneWayConfig.iceDensityKgM3) * oneWayConfig.latentHeatFusionJkg;
    const auto freeze = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    expect(freeze.converged && std::abs(freeze.iceThicknessMetres[0] - 100.0 * 86400.0 / latent) < 1.0e-7 &&
        freeze.sstC[0] == oneWayConfig.freezingTemperatureC && std::abs(freeze.relativeHeatBudgetResidual) < 1.0e-10,
        "cooling at freezing must create exactly the latent-equivalent ice and close total enthalpy");
    frozen.initialIceThicknessMetres = freeze.iceThicknessMetres;
    frozen.surfaceHeatFluxWm2.assign(cellCount, 50.0f);
    const auto halfMelt = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    expect(halfMelt.sstC[0] == oneWayConfig.freezingTemperatureC &&
        std::abs(halfMelt.iceThicknessMetres[0] - freeze.iceThicknessMetres[0] * 0.5f) < 1.0e-7 &&
        std::abs(halfMelt.relativeHeatBudgetResidual) < 1.0e-10,
        "positive heat must melt existing ice before warming the liquid");
    frozen.surfaceHeatFluxWm2.assign(cellCount, 200.0f);
    const auto melted = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    const double capacity = static_cast<double>(oneWayConfig.waterDensityKgM3) * oneWayConfig.waterHeatCapacityJkgK * oneWayConfig.mixedLayerDepthMetres;
    expect(melted.iceThicknessMetres[0] == 0.0f &&
        std::abs(melted.sstC[0] - (oneWayConfig.freezingTemperatureC + 100.0 * 86400.0 / capacity)) < 1.0e-6 &&
        std::abs(melted.relativeHeatBudgetResidual) < 1.0e-10,
        "heat beyond complete melting must warm SST without losing latent energy");
    frozen.surfaceHeatFluxWm2.clear();
    frozen.initialIceThicknessMetres.assign(cellCount, 1.0f);
    oneWayConfig.surfaceHeatExchangeWm2K = 18.0f;
    const auto coldSkin = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    expect(coldSkin.sstC[0] == oneWayConfig.freezingTemperatureC && coldSkin.iceThicknessMetres[0] > 1.0f &&
        coldSkin.surfaceSkinTemperatureC[0] < coldSkin.sstC[0] && coldSkin.surfaceSkinTemperatureC[0] > -20.0f &&
        std::abs(coldSkin.relativeHeatBudgetResidual) < 1.0e-10,
        "conductive ice skin must remain colder than liquid SST with conserved freezing energy");
    frozen.surfaceHeatFluxWm2.assign(cellCount, 400.0f);
    frozen.surfaceHeatFluxReferenceTemperatureC.assign(cellCount, -20.0f);
    frozen.initialIceThicknessMetres.assign(cellCount, 0.0f);
    const auto iceFree = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    frozen.initialIceThicknessMetres.assign(cellCount, 1.0e-6f);
    const auto thinIce = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
    const double correctedFlux = 400.0 - 18.0 * (oneWayConfig.freezingTemperatureC + 20.0);
    expect(std::abs(iceFree.sstC[0] - (oneWayConfig.freezingTemperatureC + correctedFlux * 86400.0 / capacity)) < 1.0e-6 &&
        std::abs(thinIce.sstC[0] - iceFree.sstC[0]) < 2.0e-6,
        "the diagnosed surface-flux reference must apply continuously across complete melting");
    frozen.initialIceThicknessMetres[0] = -1.0f;
    expect(!climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig).finite,
        "negative initial ice thickness must be rejected");

    // A full year must carry latent heat through each season, including the
    // year boundary. Constant flux admits an exact end state and time mean.
    frozen.initialIceThicknessMetres.assign(cellCount, 0.0f);
    frozen.surfaceHeatFluxWm2.assign(cellCount, 20.0f);
    oneWayConfig.surfaceHeatExchangeWm2K = 0.0f;
    double elapsed = 0.0;
    for (int days : {90, 91, 92, 92})
    {
        oneWayConfig.heatStepsPerIteration = days;
        const double initialT = frozen.initialSstC[0];
        const auto season = climateocean::solveWindDrivenOcean(columns, rows, frozen, oneWayConfig);
        elapsed += days * 86400.0;
        const double expected = oneWayConfig.freezingTemperatureC + 20.0 * elapsed / capacity;
        expect(season.converged && std::abs(season.sstC[0] - expected) < 2.0e-6 &&
            std::abs(season.meanSstC[0] - 0.5 * (initialT + expected)) < 2.0e-6,
            "continued ocean seasons must preserve cumulative heat and report interval means");
        frozen.initialSstC = season.sstC;
        frozen.initialIceThicknessMetres = season.iceThicknessMetres;
    }

    std::array<climateocean::OceanForcing, 4> seasonalForcing;
    std::array<climateocean::OceanConfig, 4> seasonalConfigs;
    for (int season = 0; season < 4; ++season)
    {
        seasonalForcing[season] = frozen;
        seasonalForcing[season].initialSstC.assign(cellCount, oneWayConfig.freezingTemperatureC);
        seasonalForcing[season].initialIceThicknessMetres.assign(cellCount, 0.0f);
        seasonalForcing[season].surfaceHeatFluxWm2.assign(cellCount, season < 2 ? -25.0f : 25.0f);
        seasonalConfigs[season] = oneWayConfig;
        seasonalConfigs[season].heatStepsPerIteration = 90;
    }
    const auto periodic = climateocean::solvePeriodicOcean(columns, rows, seasonalForcing, seasonalConfigs, 4, 1.0e-5);
    expect(periodic.accepted && periodic.converged && periodic.seasons[1].iceThicknessMetres[0] > 1.0f &&
        periodic.seasons[3].iceThicknessMetres[0] < 1.0e-5f,
        "the production annual solver must carry winter ice through spring and melt it before repeating the year");
    for (auto& forcing : seasonalForcing) forcing.surfaceHeatFluxWm2.assign(cellCount, -25.0f);
    const auto growingIce = climateocean::solvePeriodicOcean(columns, rows, seasonalForcing, seasonalConfigs, 2, 1.0e-5);
    expect(growingIce.accepted && !growingIce.converged && growingIce.years == 2 && growingIce.annualEnthalpyDriftK > 1.0,
        "ice growth must prevent annual convergence even when liquid SST is stationary at freezing");

    // Coarse coastal land cells are absent ocean data, not 0 C water. Include
    // wet cells across the periodic seam and both poles, plus an unsupported sea.
    constexpr int coarseColumns = 8, coarseRows = 4;
    const auto thermalPlan = climategrid::makeSmoothRemap(coarseColumns, coarseRows, columns, rows);
    std::vector<std::uint8_t> coastalLand(coarseColumns * coarseRows, 1);
    std::vector<float> coastalTemperature(coastalLand.size(), 0.0f);
    for (int y = 0; y < coarseRows; ++y)
        for (int x : {0, coarseColumns - 1})
        {
            coastalLand[y * coarseColumns + x] = 0;
            coastalTemperature[y * coarseColumns + x] = 17.25f;
        }
    std::vector<float> thermalFallback(cellCount, -31.0f);
    thermalFallback[index(columns / 2, rows / 2)] = 9.0f;
    const auto coastalConstant = climateocean::detail::remapOceanThermalField(
        thermalPlan, coastalLand, coastalTemperature, thermalFallback);
    expect(std::all_of(coastalConstant.begin(), coastalConstant.end(), [](float t) {
            return t == 17.25f || t == -31.0f || t == 9.0f;
        }) && coastalConstant[index(0, 0)] == 17.25f &&
        coastalConstant[index(columns - 1, rows - 1)] == 17.25f &&
        coastalConstant[index(columns / 2, rows / 2)] == 9.0f,
        "wet thermal constants must survive coasts, the seam and poles with explicit unsupported-cell fallback");
    for (std::size_t cell = 0; cell < coastalLand.size(); ++cell)
        if (coastalLand[cell]) coastalTemperature[cell] = std::numeric_limits<float>::quiet_NaN();
    const auto sameGridThermal = climateocean::detail::remapOceanThermalField(
        climategrid::makeSmoothRemap(coarseColumns, coarseRows, coarseColumns, coarseRows),
        coastalLand, coastalTemperature, std::vector<float>(coastalLand.size(), 9.0f));
    bool exactWetValues = true;
    for (std::size_t cell = 0; cell < coastalLand.size(); ++cell)
        exactWetValues = exactWetValues && sameGridThermal[cell] == (coastalLand[cell] ? 9.0f : coastalTemperature[cell]);
    expect(exactWetValues, "same-grid thermal transfer must preserve wet values exactly and use explicit land fallback");
    expect(climateocean::detail::remapOceanThermalField(
        thermalPlan, coastalLand, coastalTemperature, thermalFallback) == coastalConstant,
        "missing or arbitrary land temperature sentinels must not enter ocean reconstruction");

    for (int y = 0; y < coarseRows; ++y)
    {
        coastalTemperature[y * coarseColumns] = 1.0f;
        coastalTemperature[y * coarseColumns + coarseColumns - 1] = 0.0f;
    }
    const auto coastalIce = climateocean::detail::remapOceanThermalField(
        thermalPlan, coastalLand, coastalTemperature, std::vector<float>(cellCount, 0.4f));
    expect(std::all_of(coastalIce.begin(), coastalIce.end(), [](float cover) {
            return cover >= 0.0f && cover <= 1.0f;
        }) && coastalIce[index(columns - 1, rows / 2)] > 0.0f &&
        coastalIce[index(columns - 1, rows / 2)] < coastalIce[index(0, rows / 2)] &&
        coastalIce[index(0, rows / 2)] < coastalIce[index(1, rows / 2)] &&
        coastalIce[index(1, rows / 2)] < 1.0f,
        "wet-normalized ice cover must remain bounded and continuously refine the coastal seam");

    std::fill(coastalLand.begin(), coastalLand.end(), 0);
    for (std::size_t cell = 0; cell < coastalTemperature.size(); ++cell)
        coastalTemperature[cell] = static_cast<float>(cell % 7) - 2.0f;
    const auto unmaskedThermal = climateocean::detail::remapOceanThermalField(
        thermalPlan, coastalLand, coastalTemperature, thermalFallback);
    const auto smoothThermal = climategrid::remapField(thermalPlan, coastalTemperature, coarseColumns);
    bool sameAllWet = true;
    for (std::size_t cell = 0; cell < cellCount; ++cell)
        sameAllWet = sameAllWet && std::abs(unmaskedThermal[cell] - smoothThermal[cell]) < 1.0e-6f;
    expect(sameAllWet, "all-wet thermal reconstruction must retain the accepted smooth remap");

    // Restrict a mixed land/ocean forcing cell: three water donors have a
    // 100 m column, 20 C flux reference and 40 W/m2 heat input. Changing the
    // fourth, land-only donor must not alter the ocean boundary condition.
    const auto forcingPlan = climategrid::makeConservativeRemap(columns, rows, coarseColumns, coarseRows);
    std::vector<std::uint8_t> forcingLand(cellCount, 0);
    std::vector<float> fineDepth(cellCount, 100.0f), fineHeat(cellCount, 40.0f), fineReference(cellCount, 20.0f);
    for (int y = 0; y < rows; y += 2)
        for (int x = 0; x < columns; x += 2)
        {
            forcingLand[index(x, y)] = 1;
            fineDepth[index(x, y)] = 0.0f;
            fineHeat[index(x, y)] = -500.0f;
            fineReference[index(x, y)] = -30.0f;
        }
    const std::vector<float> forcingFallback(coarseColumns * coarseRows, -999.0f);
    const auto wetDepth = climateocean::detail::remapOceanThermalField(forcingPlan, forcingLand, fineDepth, forcingFallback);
    const auto wetHeat = climateocean::detail::remapOceanThermalField(forcingPlan, forcingLand, fineHeat, forcingFallback);
    const auto wetReference = climateocean::detail::remapOceanThermalField(forcingPlan, forcingLand, fineReference, forcingFallback);
    expect(std::all_of(wetDepth.begin(), wetDepth.end(), [](float d) { return std::abs(d - 100.0f) < 1.0e-5f; }),
        "restriction must preserve water depth instead of diluting it with land zeros");
    bool sameBoundaryFlux = true;
    for (std::size_t cell = 0; cell < wetHeat.size(); ++cell)
        sameBoundaryFlux = sameBoundaryFlux && std::abs(wetHeat[cell] - 18.0f * (22.0f - wetReference[cell]) - 4.0f) < 1.0e-5f;
    for (std::size_t cell = 0; cell < cellCount; ++cell)
        if (forcingLand[cell]) fineHeat[cell] = 900.0f;
    const auto changedLandHeat = climateocean::detail::remapOceanThermalField(forcingPlan, forcingLand, fineHeat, forcingFallback);
    expect(sameBoundaryFlux && changedLandHeat == wetHeat,
        "land heating must not affect the ocean's paired heat-flux/temperature boundary condition");

    // With f=0 the existing damping scale must balance the wind stress;
    // retaining only perpendicular Ekman transport would give zero current.
    auto slabForcing = aquaplanet;
    std::fill(slabForcing.initialSstC.begin(), slabForcing.initialSstC.end(), 10.0f);
    std::fill(slabForcing.atmosphericTemperatureC.begin(), slabForcing.atmosphericTemperatureC.end(), 10.0f);
    std::fill(slabForcing.eastWindMps.begin(), slabForcing.eastWindMps.end(), 8.0f);
    std::fill(slabForcing.southWindMps.begin(), slabForcing.southWindMps.end(), 6.0f);
    climateocean::OceanConfig slabConfig;
    slabConfig.oneWay = true;
    slabConfig.heatStepsPerIteration = 1;
    slabConfig.heatDiffusivityM2S = 0.0f;
    slabConfig.surfaceHeatExchangeWm2K = 0.0f;
    slabConfig.deepWaterTemperatureContrastK = 0.0f;
    slabConfig.rotationRatePerSecond = 0.0f;
    // These all-ocean zonal modes relax much more slowly than closed basins.
    slabConfig.streamfunctionIterations = 16000;
    const auto unrotatingSlab = climateocean::solveWindDrivenOcean(columns, rows, slabForcing, slabConfig);
    const double slabMass = static_cast<double>(slabConfig.waterDensityKgM3) * slabConfig.mixedLayerDepthMetres;
    const double drag = slabConfig.minimumCoriolisPerSecond;
    const double windStressFactor = static_cast<double>(slabConfig.airDensityKgM3) * slabConfig.dragCoefficient;
    const auto slabGrid = climategrid::makeSphericalGrid(columns, rows, slabConfig.planetRadiusMetres);
    bool downwindBalance = unrotatingSlab.converged;
    for (int y = 1; y < rows - 1; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto cell = index(x, y), west = slabGrid.index(x - 1, y), north = index(x, y - 1);
            // Uniform zonal stress has spherical curl even when f=0. Isolate
            // the slab component from the separately solved barotropic flow.
            const double eastBarotropic = (unrotatingSlab.eastVolumeTransportM3S[cell] +
                unrotatingSlab.eastVolumeTransportM3S[west]) / (8000.0 * slabGrid.zonalFaceLengthsMetres[y]);
            const double southBarotropic = (unrotatingSlab.southVolumeTransportM3S[cell] +
                unrotatingSlab.southVolumeTransportM3S[north]) / (4000.0 *
                    (slabGrid.northFaceLengthsMetres[y] + slabGrid.southFaceLengthsMetres[y]));
            downwindBalance = downwindBalance &&
                std::abs(unrotatingSlab.eastCurrentMps[cell] - eastBarotropic - windStressFactor * 10.0 * 8.0 / (slabMass * drag)) < 1.0e-7 &&
                std::abs(unrotatingSlab.southCurrentMps[cell] - southBarotropic - windStressFactor * 10.0 * 6.0 / (slabMass * drag)) < 1.0e-7;
        }
    expect(downwindBalance, "zero Coriolis must retain the finite downwind current that balances linear drag");

    // A fully wet pair of polar rows must obey hemisphere symmetry as well as
    // the interior. This also catches a south-face lookup beyond the last row.
    auto polarStress=slabForcing;
    std::fill(polarStress.southWindMps.begin(),polarStress.southWindMps.end(),0.0f);
    auto polarConfig=slabConfig;
    polarConfig.rotationRatePerSecond=climateocean::OceanConfig{}.rotationRatePerSecond;
    const auto symmetricSlab=climateocean::solveWindDrivenOcean(columns,rows,polarStress,polarConfig);
    bool polarSymmetry=symmetricSlab.converged;
    for (int y=0;y<rows;++y) for (int x=0;x<columns;++x)
    {
        const auto& north=symmetricSlab.transportComponents[index(x,y)];
        const auto& south=symmetricSlab.transportComponents[index(x,rows-1-y)];
        polarSymmetry=polarSymmetry && std::isfinite(north.ekmanEastMps) && std::isfinite(north.ekmanSouthMps) &&
            std::abs(north.ekmanEastMps-south.ekmanEastMps)<1.0e-9 &&
            std::abs(north.ekmanSouthMps+south.ekmanSouthMps)<1.0e-9;
    }
    expect(polarSymmetry,"uniform zonal stress must give symmetric downwind and antisymmetric crosswind Ekman flow through both polar rows");

    // Manufacture constant cell-mean currents using independent quadrature of
    // the subcell momentum equations. Constant means make face averaging exact.
    // Work must balance r * mean(|u|^2), not r * |mean(u)|^2: the latter drops
    // real within-cell velocity variation, most significant near the equator.
    slabConfig.rotationRatePerSecond = climateocean::OceanConfig{}.rotationRatePerSecond;
    // This manufactured zonal-mean stress drives a slowly relaxing global
    // barotropic mode; resolve it without weakening either equation tolerance.
    slabConfig.streamfunctionIterations = 16000;
    constexpr double desiredEast = 0.035, desiredSouth = 0.020;
    std::vector<double> rowStressEast(rows),rowStressSouth(rows),rowDissipation(rows);
    for (int y = 0; y < rows; ++y)
    {
        const double fn=2.0*slabConfig.rotationRatePerSecond*std::sin(slabGrid.latitudeNorthFacesRadians[y]);
        const double fs=2.0*slabConfig.rotationRatePerSecond*std::sin(slabGrid.latitudeSouthFacesRadians[y]);
        constexpr int samples=10000;
        double along=0.0,cross=0.0;
        for (int i=0;i<samples;++i)
        {
            const double f=fs+(i+.5)*(fn-fs)/samples;
            along+=drag/(drag*drag+f*f)/samples;
            cross+=f/(drag*drag+f*f)/samples;
        }
        const double determinant=along*along+cross*cross;
        const double stressEast=slabMass*(along*desiredEast+cross*desiredSouth)/determinant;
        const double stressSouth=slabMass*(along*desiredSouth-cross*desiredEast)/determinant;
        rowStressEast[y]=stressEast;rowStressSouth[y]=stressSouth;
        for (int i=0;i<samples;++i)
        {
            const double f=fs+(i+.5)*(fn-fs)/samples;
            const double east=(drag*stressEast-f*stressSouth)/(slabMass*(drag*drag+f*f));
            const double south=(f*stressEast+drag*stressSouth)/(slabMass*(drag*drag+f*f));
            rowDissipation[y]+=slabMass*drag*(east*east+south*south)/samples;
        }
        const double stressMagnitude = std::hypot(stressEast, stressSouth);
        const double speed = std::sqrt(stressMagnitude / windStressFactor);
        for (int x = 0; x < columns; ++x)
        {
            slabForcing.eastWindMps[index(x, y)] = static_cast<float>(speed * stressEast / stressMagnitude);
            slabForcing.southWindMps[index(x, y)] = static_cast<float>(speed * stressSouth / stressMagnitude);
        }
    }
    const auto rotatingSlab = climateocean::solveWindDrivenOcean(columns, rows, slabForcing, slabConfig);
    bool momentumBalance = rotatingSlab.converged, dissipationBalance = true;
    double maximumMomentumRelativeResidual = 0.0;
    for (int y = 1; y < rows - 1; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto cell = index(x, y), west = slabGrid.index(x - 1, y), north = index(x, y - 1);
            const double eastBarotropic = (rotatingSlab.eastVolumeTransportM3S[cell] +
                rotatingSlab.eastVolumeTransportM3S[west]) / (8000.0 * slabGrid.zonalFaceLengthsMetres[y]);
            const double southBarotropic = (rotatingSlab.southVolumeTransportM3S[cell] +
                rotatingSlab.southVolumeTransportM3S[north]) / (4000.0 *
                    (slabGrid.northFaceLengthsMetres[y] + slabGrid.southFaceLengthsMetres[y]));
            const double east = rotatingSlab.eastCurrentMps[cell] - eastBarotropic;
            const double south = rotatingSlab.southCurrentMps[cell] - southBarotropic;
            const double stressEast=rowStressEast[y],stressSouth=rowStressSouth[y];
            const double residual=std::hypot(east-desiredEast,south-desiredSouth)/
                std::hypot(desiredEast,desiredSouth);
            maximumMomentumRelativeResidual = std::max(maximumMomentumRelativeResidual,
                residual);
            if (residual >= 1.0e-5 && x == 0)
                std::cerr << "Slab row=" << y << " east=" << east << " south=" << south
                    << " ekman_east=" << rotatingSlab.transportComponents[cell].ekmanEastMps
                    << " ekman_south=" << rotatingSlab.transportComponents[cell].ekmanSouthMps << '\n';
            momentumBalance = momentumBalance && residual < 1.0e-5;
            const double windWork = stressEast * east + stressSouth * south;
            const double dissipation=rowDissipation[y];
            dissipationBalance = dissipationBalance && windWork > 0.0 &&
                std::abs(windWork - dissipation) / dissipation < 2.0e-5;
        }
    if (!momentumBalance)
        std::cerr << "Manufactured slab: maximum_momentum_relative_residual=" << maximumMomentumRelativeResidual
            << " basin_relative_residual=" << rotatingSlab.streamfunctionRelativeResidual
            << " converged=" << rotatingSlab.converged << '\n';
    expect(momentumBalance, "damped Ekman cell means must match independently integrated momentum in both hemispheres");
    expect(dissipationBalance, "cell wind work must equal integrated positive slab drag; Coriolis does no work");

    // A polar-land annulus avoids the singular poles of tau_east = A/cos(phi).
    // Its spherical curl vanishes; Cartesian d(tau_east)/dy would force a gyre.
    constexpr int annulusColumns = 64, annulusRows = 32;
    constexpr int annulusCount = annulusColumns * annulusRows;
    climateocean::OceanConfig annulusConfig;
    annulusConfig.oneWay = true;
    annulusConfig.heatStepsPerIteration = 1;
    annulusConfig.heatDiffusivityM2S = 0.0f;
    annulusConfig.surfaceHeatExchangeWm2K = 0.0f;
    annulusConfig.deepWaterTemperatureContrastK = 0.0f;
    annulusConfig.streamfunctionIterations = 16000;
    const auto annulusGrid = climategrid::makeSphericalGrid(
        annulusColumns, annulusRows, annulusConfig.planetRadiusMetres);
    climateocean::OceanForcing annulus;
    annulus.landMask.assign(annulusCount, 0);
    annulus.bathymetryMetres.assign(annulusCount, 4000.0f);
    annulus.eastWindMps.assign(annulusCount, 0.0f);
    annulus.southWindMps.assign(annulusCount, 0.0f);
    annulus.initialSstC.assign(annulusCount, 20.0f);
    annulus.atmosphericTemperatureC.assign(annulusCount, 20.0f);
    for (int y = 0; y < annulusRows; ++y)
        for (int x = 0; x < annulusColumns; ++x)
        {
            const auto cell = annulusGrid.index(x, y);
            annulus.landMask[cell] = y < 2 || y >= annulusRows - 2;
            annulus.eastWindMps[cell] = static_cast<float>(10.0 /
                std::sqrt(std::cos(annulusGrid.latitudeCentresRadians[y])));
        }
    const auto curlFree = climateocean::solveWindDrivenOcean(
        annulusColumns, annulusRows, annulus, annulusConfig);
    double maximumCurlFreeTransport = 0.0;
    for (float psi : curlFree.streamfunctionM3S)
        maximumCurlFreeTransport = std::max(maximumCurlFreeTransport, std::abs(static_cast<double>(psi)));
    expect(curlFree.converged && maximumCurlFreeTransport < 100.0,
        "spherically curl-free stress must not drive barotropic circulation beyond float-input roundoff");

    // Uniform eastward stress has spherical curl, unlike Cartesian stress.
    // With no rotation and constant depth, the closed annulus has the analytic
    // solution psi=A[phi-phi_b asinh(tan(phi))/asinh(tan(phi_b))],
    // A=tau R/(rho r). This checks the sign and matching diffusion metric.
    annulusConfig.rotationRatePerSecond = 0.0f;
    std::fill(annulus.eastWindMps.begin(), annulus.eastWindMps.end(), 10.0f);
    const auto uniformStress = climateocean::solveWindDrivenOcean(
        annulusColumns, annulusRows, annulus, annulusConfig);
    const double boundaryLatitude = annulusGrid.latitudeNorthFacesRadians[2];
    const auto secantPrimitive = [](double latitude) { return std::asinh(std::tan(latitude)); };
    const double stressScale = static_cast<double>(annulusConfig.airDensityKgM3) *
        annulusConfig.dragCoefficient * 100.0 * annulusConfig.planetRadiusMetres /
        (static_cast<double>(annulusConfig.waterDensityKgM3) * annulusConfig.barotropicDragPerSecond);
    double maximumAnalyticError = 0.0, maximumExpectedTransport = 0.0;
    for (int y = 3; y < annulusRows - 2; ++y)
        for (int x = 0; x < annulusColumns; ++x)
        {
            const double latitude = annulusGrid.latitudeNorthFacesRadians[y];
            const double expected = stressScale * (latitude - boundaryLatitude *
                secantPrimitive(latitude) / secantPrimitive(boundaryLatitude));
            maximumExpectedTransport = std::max(maximumExpectedTransport, std::abs(expected));
            maximumAnalyticError = std::max(maximumAnalyticError, std::abs(
                uniformStress.streamfunctionM3S[y * annulusColumns + x] - expected));
        }
    // The 5.625-degree latitude mesh has about 1% truncation error here.
    expect(uniformStress.converged && maximumAnalyticError / maximumExpectedTransport < 0.012 &&
        uniformStress.streamfunctionM3S[(annulusRows / 4) * annulusColumns] > 0.0f &&
        uniformStress.streamfunctionM3S[(3 * annulusRows / 4) * annulusColumns] < 0.0f,
        "uniform zonal stress must reproduce the spherical annulus solution and opposite hemisphere curl signs");
    return failures == 0 ? 0 : 1;
}
