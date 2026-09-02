#include "climate_ocean_dynamics.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace climateocean
{
bool usableOceanState(const OceanState& state, std::size_t cellCount)
{
    const auto valid = [&](const auto& values) { return values.size() == cellCount &&
        std::all_of(values.begin(), values.end(), [](auto v) { return std::isfinite(v); }); };
    return state.finite && state.converged && valid(state.sstC) && valid(state.eastCurrentMps) &&
        valid(state.southCurrentMps) && valid(state.coupledEastWindMps) && valid(state.coupledSouthWindMps);
}
OceanState solveWindDrivenOcean(int columns, int rows,
    const OceanForcing& forcing, const OceanConfig& config)
{
    OceanState state;
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) * std::max(0, rows);
    const auto valid = [&](const std::vector<float>& field)
    {
        return field.size() == count && std::all_of(field.begin(), field.end(), [](float v) { return std::isfinite(v); });
    };
    if (columns < 4 || rows < 3 || forcing.landMask.size() != count ||
        !valid(forcing.bathymetryMetres) || !valid(forcing.eastWindMps) || !valid(forcing.southWindMps) ||
        !valid(forcing.atmosphericTemperatureC) || !valid(forcing.initialSstC) ||
        (!forcing.surfaceHeatFluxWm2.empty() && !valid(forcing.surfaceHeatFluxWm2)) ||
        !(config.planetRadiusMetres > 0.0f) || !(config.waterDensityKgM3 > 0.0f) ||
        !(config.airDensityKgM3 > 0.0f) || !(config.barotropicDragPerSecond > 0.0f) ||
        !(config.mixedLayerDepthMetres > 0.0f) || !(config.waterHeatCapacityJkgK > 0.0f) ||
        !(config.oceanTimeStepSeconds > 0.0f) || config.heatDiffusivityM2S < 0.0f ||
        config.surfaceHeatExchangeWm2K < 0.0f || !(config.minimumCoriolisPerSecond > 0.0f))
    {
        state.finite = false;
        return state;
    }
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    const double dy = config.planetRadiusMetres * grid.latitudeSpacingRadians;
    const double heatCapacity = config.waterDensityKgM3 * config.waterHeatCapacityJkgK * config.mixedLayerDepthMetres;
    const auto ocean = [&](std::size_t cell) { return !forcing.landMask[cell] && forcing.bathymetryMetres[cell] > 0.0f; };
    const auto depth = [&](std::size_t cell) { return std::max(50.0f, forcing.bathymetryMetres[cell]); };
    const auto vertex = [&](int x, int y) { return static_cast<std::size_t>(y) * columns + grid.wrapColumn(x); };
    const std::size_t vertices = static_cast<std::size_t>(columns) * (rows + 1);
    std::vector<bool> active(vertices, false);
    std::vector<double> inverseDepth(vertices, 0.0), q(vertices, 0.0), psi(vertices, 0.0);
    for (int y = 1; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto v = vertex(x, y);
            const auto nw = grid.index(x - 1, y - 1), ne = grid.index(x, y - 1);
            const auto sw = grid.index(x - 1, y), se = grid.index(x, y);
            active[v] = ocean(nw) && ocean(ne) && ocean(sw) && ocean(se);
            inverseDepth[v] = 4.0 / (depth(nw) + depth(ne) + depth(sw) + depth(se));
            q[v] = 2.0 * config.rotationRatePerSecond * config.rotationDirection *
                std::sin(grid.latitudeNorthFacesRadians[y]) * inverseDepth[v];
        }
    state.eastCurrentMps.assign(count, 0.0f);
    state.southCurrentMps.assign(count, 0.0f);
    state.eastVolumeTransportM3S.assign(count, 0.0);
    state.southVolumeTransportM3S.assign(count, 0.0);
    state.ekmanUpwellingMps.assign(count, 0.0f);
    state.sstC = forcing.initialSstC;
    state.coupledEastWindMps = forcing.eastWindMps;
    state.coupledSouthWindMps = forcing.southWindMps;
    state.coupledPressureAnomalyHpa.assign(count, 0.0f);
    std::vector<double> stressU(count), stressV(count), ekmanU(count), ekmanV(count);
    std::vector<double> eastHeatTransport(count), southHeatTransport(count);
    std::vector<double> eastDiffusion(count), southDiffusion(count), divergence(count), tendency(count);
    std::vector<double> aE(vertices), aW(vertices), aN(vertices), aS(vertices), diagonal(vertices), rhs(vertices);
    const int iterations = config.oneWay ? 1 : std::max(1, config.couplingIterations);
    for (int coupling = 0; coupling < iterations; ++coupling)
    {
        state.couplingIterations = coupling + 1;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y);
                if (!ocean(c)) continue;
                const double speed = std::hypot(state.coupledEastWindMps[c], state.coupledSouthWindMps[c]);
                stressU[c] = config.airDensityKgM3 * config.dragCoefficient * speed * state.coupledEastWindMps[c];
                stressV[c] = config.airDensityKgM3 * config.dragCoefficient * speed * state.coupledSouthWindMps[c];
                const double f = 2.0 * config.rotationRatePerSecond * config.rotationDirection * std::sin(grid.latitudeCentresRadians[y]);
                // Rayleigh-regularized Ekman balance: no singular equatorial f clamp.
                const double factor = f / (f * f + config.minimumCoriolisPerSecond * config.minimumCoriolisPerSecond);
                ekmanU[c] = -stressV[c] * factor / config.waterDensityKgM3;
                ekmanV[c] = stressU[c] * factor / config.waterDensityKgM3;
            }
        // Finite-volume Stommel/PV closure at corners. Dirichlet psi=0 on every
        // coast/pole means exactly zero normal transport, including islands.
        double rhsNorm = 0.0;
        for (int y = 1; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto v = vertex(x, y);
                if (!active[v]) continue;
                const double dx = config.planetRadiusMetres * grid.longitudeSpacingRadians * std::cos(grid.latitudeNorthFacesRadians[y]);
                const auto nw = grid.index(x - 1, y - 1), ne = grid.index(x, y - 1);
                const auto sw = grid.index(x - 1, y), se = grid.index(x, y);
                const auto qAt = [&](int xx, int yy) { const auto vv = vertex(xx, yy); return active[vv] ? q[vv] : q[v]; };
                const double betaNorth = (qAt(x, y - 1) - qAt(x, y + 1)) / (2.0 * dy);
                const double betaEast = (qAt(x + 1, y) - qAt(x - 1, y)) / (2.0 * dx);
                const auto diffusion = [&](int xx, int yy, double spacing)
                {
                    const auto vv = vertex(xx, yy);
                    const double inv = active[vv] ? 0.5 * (inverseDepth[v] + inverseDepth[vv]) : inverseDepth[v];
                    return config.barotropicDragPerSecond * inv / (spacing * spacing);
                };
                aE[v] = diffusion(x + 1, y, dx) + std::max(0.0, betaNorth) / dx;
                aW[v] = diffusion(x - 1, y, dx) + std::max(0.0, -betaNorth) / dx;
                aS[v] = diffusion(x, y + 1, dy) + std::max(0.0, betaEast) / dy;
                aN[v] = diffusion(x, y - 1, dy) + std::max(0.0, -betaEast) / dy;
                diagonal[v] = aE[v] + aW[v] + aN[v] + aS[v];
                rhs[v] = -((stressV[ne] / depth(ne) + stressV[se] / depth(se) - stressV[nw] / depth(nw) - stressV[sw] / depth(sw)) / (2.0 * dx) -
                    (stressU[sw] / depth(sw) + stressU[se] / depth(se) - stressU[nw] / depth(nw) - stressU[ne] / depth(ne)) / (2.0 * dy)) / config.waterDensityKgM3;
                rhsNorm += rhs[v] * rhs[v];
            }
        for (int sweep = 0; sweep < std::max(1, config.streamfunctionIterations); ++sweep)
        {
            for (int y = 1; y < rows; ++y)
                for (int x = columns - 1; x >= 0; --x)
                {
                    const auto v = vertex(x, y);
                    if (active[v]) psi[v] = (rhs[v] + aE[v] * psi[vertex(x + 1, y)] + aW[v] * psi[vertex(x - 1, y)] +
                        aN[v] * psi[vertex(x, y - 1)] + aS[v] * psi[vertex(x, y + 1)]) / diagonal[v];
                }
            double residual = 0.0;
            for (int y = 1; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const auto v = vertex(x, y);
                    if (!active[v]) continue;
                    const double r = rhs[v] - diagonal[v] * psi[v] + aE[v] * psi[vertex(x + 1, y)] + aW[v] * psi[vertex(x - 1, y)] +
                        aN[v] * psi[vertex(x, y - 1)] + aS[v] * psi[vertex(x, y + 1)];
                    residual += r * r;
                }
            state.streamfunctionRelativeResidual = rhsNorm > 0.0 ? static_cast<float>(std::sqrt(residual / rhsNorm)) : 0.0f;
            if (state.streamfunctionRelativeResidual <= config.streamfunctionTolerance) break;
        }
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y), e = grid.index(x + 1, y), s = grid.index(x, y + 1);
                state.eastVolumeTransportM3S[c] = -(psi[vertex(x + 1, y + 1)] - psi[vertex(x + 1, y)]);
                state.southVolumeTransportM3S[c] = psi[vertex(x + 1, y + 1)] - psi[vertex(x, y + 1)];
                if (ocean(c) && ocean(e))
                {
                    const double faceDepth = 0.5 * (depth(c) + depth(e));
                    eastHeatTransport[c] = state.eastVolumeTransportM3S[c] / faceDepth +
                        0.5 * (ekmanU[c] + ekmanU[e]) * grid.zonalFaceLengthsMetres[y] / config.mixedLayerDepthMetres;
                    eastDiffusion[c] = config.heatDiffusivityM2S * grid.zonalFaceLengthsMetres[y] * grid.zonalFaceLengthsMetres[y] /
                        grid.cellAreasSquareMetres[y];
                }
                if (y + 1 < rows && ocean(c) && ocean(s))
                {
                    const double faceDepth = 0.5 * (depth(c) + depth(s));
                    southHeatTransport[c] = state.southVolumeTransportM3S[c] / faceDepth +
                        0.5 * (ekmanV[c] + ekmanV[s]) * grid.southFaceLengthsMetres[y] / config.mixedLayerDepthMetres;
                    southDiffusion[c] = config.heatDiffusivityM2S * grid.southFaceLengthsMetres[y] / dy;
                }
            }
        double maximumRate = 0.0;
        state.maximumTransportDivergenceMps = 0.0;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y), w = grid.index(x - 1, y), n = grid.index(x, y - 1);
                if (!ocean(c)) continue;
                const double area = grid.cellAreasSquareMetres[y];
                const double northHeat = y > 0 ? southHeatTransport[n] : 0.0;
                const double northVolume = y > 0 ? state.southVolumeTransportM3S[n] : 0.0;
                const double massResidual = (state.eastVolumeTransportM3S[c] - state.eastVolumeTransportM3S[w] +
                    state.southVolumeTransportM3S[c] - northVolume) / area;
                state.maximumTransportDivergenceMps = std::max(state.maximumTransportDivergenceMps, std::abs(massResidual));
                divergence[c] = (eastHeatTransport[c] - eastHeatTransport[w] + southHeatTransport[c] - northHeat) / area;
                state.ekmanUpwellingMps[c] = static_cast<float>(divergence[c] * config.mixedLayerDepthMetres);
                state.eastCurrentMps[c] = static_cast<float>(0.5 * (eastHeatTransport[c] + eastHeatTransport[w]) / grid.zonalFaceLengthsMetres[y]);
                state.southCurrentMps[c] = static_cast<float>((southHeatTransport[c] + northHeat) /
                    std::max(1.0, grid.northFaceLengthsMetres[y] + grid.southFaceLengthsMetres[y]));
                const double outgoing = std::max(0.0, eastHeatTransport[c]) + std::max(0.0, -eastHeatTransport[w]) +
                    std::max(0.0, southHeatTransport[c]) + std::max(0.0, -northHeat);
                const double diffusionRate = eastDiffusion[c] + eastDiffusion[w] + southDiffusion[c] + (y > 0 ? southDiffusion[n] : 0.0);
                maximumRate = std::max(maximumRate, (outgoing + diffusionRate) / area + std::max(0.0, -divergence[c]) +
                    config.surfaceHeatExchangeWm2K / heatCapacity);
            }
        // Each outer iteration solves the SAME finite seasonal interval from
        // T0, not a progressively longer integration mistaken for convergence.
        std::vector<double> temperature(forcing.initialSstC.begin(), forcing.initialSstC.end());
        const double interval = config.oceanTimeStepSeconds * std::max(1, config.heatStepsPerIteration);
        const int substeps = std::max(1, static_cast<int>(std::ceil(interval * maximumRate / 0.7)));
        const double dt = interval / substeps;
        double expectedHeatChange = 0.0, absoluteHeatExchange = 0.0;
        for (int step = 0; step < substeps; ++step)
        {
            std::fill(tendency.begin(), tendency.end(), 0.0);
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const auto c = grid.index(x, y), e = grid.index(x + 1, y), s = grid.index(x, y + 1);
                    const auto flux = [&](std::size_t neighbour, double transport, double diffusion)
                    {
                        const double exchange = transport * (transport >= 0.0 ? temperature[c] : temperature[neighbour]) +
                            diffusion * (temperature[c] - temperature[neighbour]);
                        tendency[c] -= exchange;
                        tendency[neighbour] += exchange;
                    };
                    flux(e, eastHeatTransport[c], eastDiffusion[c]);
                    if (y + 1 < rows) flux(s, southHeatTransport[c], southDiffusion[c]);
                }
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x)
                {
                    const auto c = grid.index(x, y);
                    if (!ocean(c)) continue;
                    const double area = grid.cellAreasSquareMetres[y];
                    // Fixed-depth continuity closes horizontal convergence via
                    // vertical exchange with the deep reservoir. Uniform T is
                    // invariant when the deep contrast and heat exchange vanish.
                    const double deepT = forcing.initialSstC[c] - config.deepWaterTemperatureContrastK;
                    const double surfaceSource = forcing.surfaceHeatFluxWm2.empty()
                        ? config.surfaceHeatExchangeWm2K * (forcing.atmosphericTemperatureC[c] - temperature[c])
                        : forcing.surfaceHeatFluxWm2[c];
                    const double source = surfaceSource / heatCapacity +
                        divergence[c] * (divergence[c] > 0.0 ? deepT : temperature[c]);
                    expectedHeatChange += dt * area * heatCapacity * source;
                    absoluteHeatExchange += std::abs(dt * area * heatCapacity * source);
                    temperature[c] += dt * (tendency[c] / area + source);
                }
        }
        double heatChange = 0.0, residual = 0.0, totalArea = 0.0;
        const float relaxation = config.oneWay ? 1.0f : std::clamp(config.underRelaxation, 0.01f, 1.0f);
        std::vector<float> targetPressure(count, 0.0f);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y);
                if (!ocean(c)) continue;
                const double area = grid.cellAreasSquareMetres[y];
                heatChange += area * heatCapacity * (temperature[c] - forcing.initialSstC[c]);
                const double delta = temperature[c] - state.sstC[c];
                residual += area * delta * delta; // normalized by 1 K below
                totalArea += area;
                state.sstC[c] += relaxation * static_cast<float>(delta);
                targetPressure[c] = -config.sstWindFeedbackMpsPerK * 10.0f *
                    (state.sstC[c] - forcing.atmosphericTemperatureC[c]);
            }
        state.heatBudgetResidualJ = heatChange - expectedHeatChange;
        state.relativeHeatBudgetResidual = state.heatBudgetResidualJ /
            std::max({1.0, std::abs(heatChange), absoluteHeatExchange});
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y);
                if (!ocean(c) || config.oneWay) continue;
                const auto at = [&](int xx, int yy)
                {
                    const auto other = grid.index(xx, std::clamp(yy, 0, rows - 1));
                    return ocean(other) ? targetPressure[other] : targetPressure[c];
                };
                const double dx = grid.cellAreasSquareMetres[y] / grid.zonalFaceLengthsMetres[y];
                const double eastForce = -(at(x + 1, y) - at(x - 1, y)) * 100.0 / (2.0 * dx * config.airDensityKgM3);
                const double northForce = -(at(x, y - 1) - at(x, y + 1)) * 100.0 / (2.0 * dy * config.airDensityKgM3);
                const double f = 2.0 * config.rotationRatePerSecond * config.rotationDirection * std::sin(grid.latitudeCentresRadians[y]);
                const double drag = 1.0 / 43200.0;
                const double deltaU = forcing.eastWindMps[c] + (drag * eastForce + f * northForce) / (drag * drag + f * f) - state.coupledEastWindMps[c];
                const double deltaV = forcing.southWindMps[c] - (drag * northForce - f * eastForce) / (drag * drag + f * f) - state.coupledSouthWindMps[c];
                residual += grid.cellAreasSquareMetres[y] * (deltaU * deltaU + deltaV * deltaV); // 1 m/s reference
                state.coupledEastWindMps[c] += relaxation * static_cast<float>(deltaU);
                state.coupledSouthWindMps[c] += relaxation * static_cast<float>(deltaV);
                state.coupledPressureAnomalyHpa[c] += relaxation * (targetPressure[c] - state.coupledPressureAnomalyHpa[c]);
            }
        state.relativeResidual = totalArea > 0.0 ? static_cast<float>(std::sqrt(residual / totalArea)) : 0.0f;
        state.residualHistory.push_back(state.relativeResidual);
        state.finite = valid(state.sstC) && valid(state.coupledEastWindMps) && valid(state.coupledSouthWindMps) &&
            std::isfinite(state.streamfunctionRelativeResidual);
        state.converged = state.finite && state.streamfunctionRelativeResidual <= config.streamfunctionTolerance &&
            (config.oneWay || state.relativeResidual <= config.convergenceTolerance);
        if (!state.finite || state.converged) break;
    }
    state.streamfunctionM3S.resize(vertices);
    std::transform(psi.begin(), psi.end(), state.streamfunctionM3S.begin(), [](double v) { return static_cast<float>(v); });
    return state;
}
}
