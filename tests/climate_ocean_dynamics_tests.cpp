#include "climate_ocean_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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
}

int main()
{
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

    return failures == 0 ? 0 : 1;
}
