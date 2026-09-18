#include "climate_ocean_dynamics.hpp"
#include "climate_grid.hpp"
#include "detail/ocean_zonal_diffusion.hpp"
#include "detail/ocean_mixed_layer.hpp"
#include "detail/ocean_circulation_solver.hpp"
#include "legacy/ocean_feedback.hpp"
#include "parallel_rows.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <utility>

namespace climateocean
{
namespace
{
// Owned by one seasonal forcing inside solvePeriodicOcean. Only initial liquid/
// ice energy changes during that call; momentum has no thermal feedback there.
struct PreparedCirculation
{
    bool ready = false;
    OceanState diagnostics;
    std::vector<double> psi, eastHeat, southHeat, eastDiffusion, southDiffusion, divergence;
    double maximumRate = 0.0, maximumExplicitRate = 0.0;
};

void copyCirculation(OceanState& to, const OceanState& from)
{
    to.eastVolumeTransportM3S = from.eastVolumeTransportM3S;
    to.southVolumeTransportM3S = from.southVolumeTransportM3S;
    to.eastCurrentMps = from.eastCurrentMps;
    to.southCurrentMps = from.southCurrentMps;
    to.ekmanUpwellingMps = from.ekmanUpwellingMps;
    to.transportComponents = from.transportComponents;
    to.streamfunctionRelativeResidual = from.streamfunctionRelativeResidual;
    to.maximumTransportDivergenceMps = from.maximumTransportDivergenceMps;
    to.stabilityLimit = from.stabilityLimit;
    to.stabilityDistribution = from.stabilityDistribution;
}

OceanState solveOcean(int columns, int rows, const OceanForcing& forcing,
    const OceanConfig& config, PreparedCirculation* prepared);
}

bool usableOceanState(const OceanState& state, std::size_t cellCount)
{
    const auto valid = [&](const auto& values) { return values.size() == cellCount &&
        std::all_of(values.begin(), values.end(), [](auto v) { return std::isfinite(v); }); };
    return state.finite && state.converged && valid(state.sstC) && valid(state.eastCurrentMps) &&
        valid(state.southCurrentMps) && valid(state.coupledEastWindMps) && valid(state.coupledSouthWindMps) &&
        valid(state.surfaceSkinTemperatureC) && valid(state.iceThicknessMetres);
}
PeriodicOceanState solvePeriodicOcean(int columns, int rows,
    std::array<OceanForcing, 4> forcing, std::array<OceanConfig, 4> configs,
    int maximumYears, double toleranceK)
{
    PeriodicOceanState result;
    if (!climategrid::validGlobalGridDimensions(columns, rows) || maximumYears < 1 ||
        !std::isfinite(toleranceK) || toleranceK < 0.0) return result;
    auto sst = forcing[0].initialSstC, ice = forcing[0].initialIceThicknessMetres;
    const auto count = static_cast<std::size_t>(columns) * rows;
    if (sst.size() != count || forcing[0].landMask.size() != count || forcing[0].bathymetryMetres.size() != count) return result;
    if (ice.empty()) ice.assign(count, 0.0f);
    if (ice.size() != count) return result;
    for (int season = 0; season < 4; ++season)
    {
        const auto& first = configs[0];
        const auto& config = configs[season];
        if (forcing[season].landMask != forcing[0].landMask ||
            forcing[season].bathymetryMetres != forcing[0].bathymetryMetres ||
            config.waterDensityKgM3 != first.waterDensityKgM3 ||
            config.waterHeatCapacityJkgK != first.waterHeatCapacityJkgK ||
            config.mixedLayerDepthMetres != first.mixedLayerDepthMetres ||
            config.variableHeatStorage != first.variableHeatStorage ||
            config.minimumStorageDepthMetres != first.minimumStorageDepthMetres ||
            config.maximumStorageDepthMetres != first.maximumStorageDepthMetres ||
            config.storageColumnDepthMetres != first.storageColumnDepthMetres ||
            config.freezingTemperatureC != first.freezingTemperatureC ||
            config.iceDensityKgM3 != first.iceDensityKgM3 || config.latentHeatFusionJkg != first.latentHeatFusionJkg)
            return result;
        if (forcing[season].deepWaterTemperatureC.empty())
        {
            forcing[season].deepWaterTemperatureC = sst;
            for (auto& value : forcing[season].deepWaterTemperatureC)
                value = std::max(config.freezingTemperatureC, value - config.deepWaterTemperatureContrastK);
        }
    }
    auto storageDepth = forcing[0].initialStorageDepthMetres;
    auto reservoirTemperature = forcing[0].initialReservoirTemperatureC;
    if (configs[0].variableHeatStorage)
    {
        if (storageDepth.empty())
        {
            storageDepth.assign(count, configs[0].mixedLayerDepthMetres);
            for (std::size_t c = 0; c < count; ++c) if (!forcing[0].landMask[c])
            {
                const double column = detail::storageColumnDepth(forcing[0].bathymetryMetres[c], configs[0]);
                storageDepth[c] = static_cast<float>(std::clamp(double(configs[0].mixedLayerDepthMetres),
                    detail::minimumStorageDepth(column, configs[0]), detail::maximumStorageDepth(column, configs[0])));
            }
        }
        if (reservoirTemperature.empty()) reservoirTemperature = forcing[0].deepWaterTemperatureC;
        if (storageDepth.size() != count || reservoirTemperature.size() != count) return result;
    }
    std::array<PreparedCirculation, 4> circulation;
    for (int year = 0; year < maximumYears; ++year)
    {
        const auto initialSst = sst, initialIce = ice;
        const auto initialDepth = storageDepth, initialReservoirTemperature = reservoirTemperature;
        for (int season = 0; season < 4; ++season)
        {
            forcing[season].initialSstC = sst;
            forcing[season].initialIceThicknessMetres = ice;
            if (configs[season].variableHeatStorage)
            {
                forcing[season].initialStorageDepthMetres = storageDepth;
                forcing[season].initialReservoirTemperatureC = reservoirTemperature;
            }
            configs[season].oneWay = true;
            auto& state = result.seasons[season];
            state = solveOcean(columns, rows, forcing[season], configs[season],
                configs[season].reuseSeasonalCirculation ? &circulation[season] : nullptr);
            result.circulationSeconds += state.circulationSeconds;
            result.heatSeconds += state.heatSeconds;
            result.zonalHeatSeconds += state.zonalHeatSeconds;
            result.heatFluxSeconds += state.heatFluxSeconds;
            result.heatUpdateSeconds += state.heatUpdateSeconds;
            result.circulationSolves += state.circulationSolves;
            result.heatSubsteps += state.heatSubsteps;
            if (!usableOceanState(state, count)) { result.accepted = false; return result; }
            sst = state.sstC; ice = state.iceThicknessMetres;
            if (configs[season].variableHeatStorage)
            {
                storageDepth = state.storageDepthMetres;
                reservoirTemperature = state.reservoirTemperatureC;
            }
        }
        const auto& config = configs[0];
        const double latentOverCapacity = static_cast<double>(config.iceDensityKgM3) * config.latentHeatFusionJkg /
            (static_cast<double>(config.waterDensityKgM3) * config.waterHeatCapacityJkgK * config.mixedLayerDepthMetres);
        double change = 0.0, surfaceChange = 0.0, area = 0.0;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = y * columns + x;
                if (forcing[0].landMask[c] || forcing[0].bathymetryMetres[c] <= 0.0f) continue;
                const double weight = climategrid::latitudeBandMeasure(y, rows);
                double delta = sst[c] - initialSst[c] - latentOverCapacity * (ice[c] - initialIce[c]);
                surfaceChange += weight * delta * delta;
                if (config.variableHeatStorage)
                {
                    const double column = detail::storageColumnDepth(forcing[0].bathymetryMetres[c], config);
                    const auto sensible = [&](double h, double t, double r)
                    { return h * (t - config.freezingTemperatureC) + (column - h) * (r - config.freezingTemperatureC); };
                    delta = (sensible(storageDepth[c], sst[c], reservoirTemperature[c]) -
                        sensible(initialDepth[c], initialSst[c], initialReservoirTemperature[c])) / config.mixedLayerDepthMetres -
                        latentOverCapacity * (ice[c] - initialIce[c]);
                }
                change += weight * delta * delta; area += weight;
            }
        result.years = year + 1;
        result.accepted = true;
        result.annualEnthalpyDriftK = std::sqrt(std::max(change, surfaceChange) / std::max(1.0e-30, area));
        result.converged = result.annualEnthalpyDriftK <= toleranceK;
        if (result.converged) break;
    }
    return result;
}

OceanState solveWindDrivenOcean(int columns, int rows,
    const OceanForcing& forcing, const OceanConfig& config)
{
    return solveOcean(columns, rows, forcing, config, nullptr);
}

namespace
{
OceanState solveOcean(int columns, int rows,
    const OceanForcing& forcing, const OceanConfig& config, PreparedCirculation* prepared)
{
    const auto preparationStart = std::chrono::steady_clock::now();
    const bool reuse = prepared && prepared->ready && config.oneWay;
    OceanState state;
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) * std::max(0, rows);
    const auto valid = [&](const std::vector<float>& field)
    {
        return field.size() == count && std::all_of(field.begin(), field.end(), [](float v) { return std::isfinite(v); });
    };
    if (!climategrid::validGlobalGridDimensions(columns, rows) ||
        columns < 4 || forcing.landMask.size() != count ||
        !valid(forcing.bathymetryMetres) || !valid(forcing.eastWindMps) || !valid(forcing.southWindMps) ||
        !valid(forcing.atmosphericTemperatureC) || !valid(forcing.initialSstC) ||
        (!forcing.initialIceThicknessMetres.empty() && (!valid(forcing.initialIceThicknessMetres) ||
            std::any_of(forcing.initialIceThicknessMetres.begin(), forcing.initialIceThicknessMetres.end(), [](float h) { return h < 0.0f; }))) ||
        (!forcing.surfaceHeatFluxReferenceTemperatureC.empty() && !valid(forcing.surfaceHeatFluxReferenceTemperatureC)) ||
        (!forcing.surfaceHeatFluxWm2.empty() && !valid(forcing.surfaceHeatFluxWm2)) ||
        (!forcing.deepWaterTemperatureC.empty() && !valid(forcing.deepWaterTemperatureC)) ||
        (!forcing.initialStorageDepthMetres.empty() && !valid(forcing.initialStorageDepthMetres)) ||
        (!forcing.initialReservoirTemperatureC.empty() && !valid(forcing.initialReservoirTemperatureC)) ||
        (config.variableHeatStorage && (!config.oneWay ||
            !(config.minimumStorageDepthMetres > 0.0f) ||
            !(config.maximumStorageDepthMetres >= config.minimumStorageDepthMetres) ||
            !(config.storageColumnDepthMetres > config.maximumStorageDepthMetres) ||
            !(config.storageStratificationPerSecond2 > 0.0f) ||
            !(config.storageWindMixingEfficiency >= 0.0f) || !std::isfinite(config.storageWindMixingEfficiency) ||
            !(config.storageMixingMemoryDays > 0.0f) || !(config.storageAdjustmentDays > 0.0f))) ||
        !(config.planetRadiusMetres > 0.0f) || !(config.waterDensityKgM3 > 0.0f) ||
        !(config.airDensityKgM3 > 0.0f) || !(config.barotropicDragPerSecond > 0.0f) ||
        !std::isfinite(config.linearBottomDragMps) || config.linearBottomDragMps < 0.0f ||
        !(config.mixedLayerDepthMetres > 0.0f) || !(config.waterHeatCapacityJkgK > 0.0f) ||
        !(config.oceanTimeStepSeconds > 0.0f) || config.heatDiffusivityM2S < 0.0f ||
        !std::isfinite(config.freezingTemperatureC) || !(config.iceDensityKgM3 > 0.0f) ||
        !(config.latentHeatFusionJkg > 0.0f) || !(config.iceConductivityWmK > 0.0f) ||
        config.surfaceHeatExchangeWm2K < 0.0f || !(config.minimumCoriolisPerSecond > 0.0f) ||
        config.heatWorkers < 0 || config.heatWorkers > 32 ||
        (config.coastalScheme != CoastalScheme::Legacy && config.coastalScheme != CoastalScheme::FaceDrag &&
         config.coastalScheme != CoastalScheme::FaceForcing))
    {
        state.finite = false;
        return state;
    }
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    const double dy = config.planetRadiusMetres * grid.latitudeSpacingRadians;
    const double heatCapacity = config.waterDensityKgM3 * config.waterHeatCapacityJkgK * config.mixedLayerDepthMetres;
    const double volumetricCapacity = static_cast<double>(config.waterDensityKgM3) * config.waterHeatCapacityJkgK;
    const double freezing = config.freezingTemperatureC;
    const double latentHeatPerMetre = config.iceDensityKgM3 * config.latentHeatFusionJkg;
    const int requestedWorkers = config.heatWorkers ? config.heatWorkers :
        (columns >= 256 ? (config.gatherHeatFluxes ? 16 : std::min(16, rows / 16)) : 1);
    const int rowsPerWorker = (rows + requestedWorkers - 1) / requestedWorkers;
    state.heatWorkersUsed = std::min<int>(rowworkers::concurrency(), rows / rowsPerWorker);
    const auto heatRows = [&](auto&& work)
    {
        parallelforrows(0, rows - 1, work, rowsPerWorker);
    };
    const auto ocean = [&](std::size_t cell) { return !forcing.landMask[cell] && forcing.bathymetryMetres[cell] > 0.0f; };
    const auto depth = [&](std::size_t cell) { return std::max(50.0f, forcing.bathymetryMetres[cell]); };
    const auto vertex = [&](int x, int y) { return static_cast<std::size_t>(y) * columns + grid.wrapColumn(x); };
    const std::size_t vertices = static_cast<std::size_t>(columns) * (rows + 1);
    std::vector<bool> active(vertices, false);
    std::vector<double> inverseDepth(vertices, 0.0), q(vertices, 0.0), psi(vertices, 0.0);
    if (!reuse) for (int y = 1; y < rows; ++y)
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
    state.transportComponents.assign(count, {});
    state.sstC = forcing.initialSstC;
    state.surfaceSkinTemperatureC = forcing.initialSstC;
    state.iceThicknessMetres.assign(count, 0.0f);
    // Enthalpy is relative to liquid water at freezing. Negative values store
    // latent energy in stationary ice; only liquid sensible heat is advected.
    std::vector<double> initialEnthalpy(count), coupledEnthalpy(count);
    std::vector<double> initialStorageDepth(count, config.mixedLayerDepthMetres), initialReservoirEnthalpy(count, 0.0);
    std::vector<double> capacities(count, heatCapacity), columnDepth(count, config.mixedLayerDepthMetres);
    for (std::size_t c = 0; c < count; ++c) if (config.variableHeatStorage && ocean(c))
    {
        columnDepth[c] = detail::storageColumnDepth(forcing.bathymetryMetres[c], config);
        const double minimum = detail::minimumStorageDepth(columnDepth[c], config);
        const double maximum = detail::maximumStorageDepth(columnDepth[c], config);
        initialStorageDepth[c] = forcing.initialStorageDepthMetres.empty()
            ? std::clamp(double(config.mixedLayerDepthMetres), minimum, maximum) : forcing.initialStorageDepthMetres[c];
        const double lowerT = forcing.initialReservoirTemperatureC.empty()
            ? (forcing.deepWaterTemperatureC.empty() ? std::max(freezing, double(forcing.initialSstC[c]) - config.deepWaterTemperatureContrastK)
                : double(forcing.deepWaterTemperatureC[c])) : double(forcing.initialReservoirTemperatureC[c]);
        if (initialStorageDepth[c] < minimum - 1.0e-5 || initialStorageDepth[c] > maximum + 1.0e-5 || lowerT < freezing)
        { state.finite = false; return state; }
        initialStorageDepth[c] = std::clamp(initialStorageDepth[c], minimum, maximum);
        capacities[c] = volumetricCapacity * initialStorageDepth[c];
        initialReservoirEnthalpy[c] = volumetricCapacity * (columnDepth[c] - initialStorageDepth[c]) * (lowerT - freezing);
    }
    for (std::size_t c = 0; c < count; ++c)
    {
        initialEnthalpy[c] = capacities[c] * (forcing.initialSstC[c] - freezing) -
            (forcing.initialIceThicknessMetres.empty() ? 0.0 : latentHeatPerMetre * forcing.initialIceThicknessMetres[c]);
        coupledEnthalpy[c] = initialEnthalpy[c];
        if (ocean(c)) state.sstC[c] = static_cast<float>(freezing + std::max(0.0, initialEnthalpy[c]) / capacities[c]);
    }
    const auto surfaceExchange = [&](std::size_t c, double enthalpy, double liquidT)
    {
        const double airT = forcing.atmosphericTemperatureC[c];
        const double exchange = config.surfaceHeatExchangeWm2K;
        const double referenceT = forcing.surfaceHeatFluxReferenceTemperatureC.empty()
            ? airT : forcing.surfaceHeatFluxReferenceTemperatureC[c];
        const double referenceFlux = forcing.surfaceHeatFluxWm2.empty() ? 0.0 : forcing.surfaceHeatFluxWm2[c];
        double skin = liquidT;
        if (enthalpy < 0.0)
        {
            const double conductance = config.iceConductivityWmK / std::max(1.0e-6, -enthalpy / latentHeatPerMetre);
            // Zero heat-capacity ice: conduction balances the linearized
            // atmospheric exchange; the surface is capped at melting.
            const double sourceT = forcing.surfaceHeatFluxWm2.empty() ? airT : referenceT;
            skin = std::min(freezing, (referenceFlux + exchange * sourceT + conductance * freezing) / (exchange + conductance));
        }
        const double flux = forcing.surfaceHeatFluxWm2.empty() ? exchange * (airT - skin)
            : referenceFlux - exchange * (skin - referenceT);
        return std::pair<double, double>{skin, flux};
    };
    state.coupledEastWindMps = forcing.eastWindMps;
    state.coupledSouthWindMps = forcing.southWindMps;
    state.coupledPressureAnomalyHpa.assign(count, 0.0f);
    std::vector<double> stressU(count), stressV(count), ekmanU(count), ekmanV(count);
    std::vector<double> eastHeatTransport(count), southHeatTransport(count);
    std::vector<double> eastBarotropic(count), southBarotropic(count);
    std::vector<double> eastDiffusion(count), southDiffusion(count), divergence(count), tendency(count);
    std::vector<double> eastFlux(config.gatherHeatFluxes ? count : 0), southFlux(eastFlux.size());
    std::vector<double> rowExpected(rows, 0.0), rowAbsolute(rows, 0.0), deepTemperature(count);
    for (std::size_t c = 0; c < count; ++c)
        deepTemperature[c] = forcing.deepWaterTemperatureC.empty()
            ? std::max(freezing, static_cast<double>(forcing.initialSstC[c] - config.deepWaterTemperatureContrastK))
            : forcing.deepWaterTemperatureC[c];
    std::vector<int> rowIterations(rows, 0);
    std::vector<double> rowResidual(rows, 0.0);
    std::vector<unsigned char> rowValid(rows, 1);
    std::vector<double> aE(vertices), aW(vertices), aN(vertices), aS(vertices), diagonal(vertices), rhs(vertices);
    const int iterations = config.oneWay ? 1 : std::max(1, config.couplingIterations);
    for (int coupling = 0; coupling < iterations; ++coupling)
    {
        state.couplingIterations = coupling + 1;
        const auto circulationStart = coupling == 0 ? preparationStart : std::chrono::steady_clock::now();
        double maximumRate = 0.0, maximumExplicitRate = 0.0;
        if (reuse)
        {
            copyCirculation(state, prepared->diagnostics);
            psi = prepared->psi;
            eastHeatTransport = prepared->eastHeat; southHeatTransport = prepared->southHeat;
            eastDiffusion = prepared->eastDiffusion; southDiffusion = prepared->southDiffusion;
            divergence = prepared->divergence;
            maximumRate = prepared->maximumRate; maximumExplicitRate = prepared->maximumExplicitRate;
        }
        else
        {
        ++state.circulationSolves;
        state.stabilityLimit = {};
        state.stabilityDistribution = {};
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y);
                if (!ocean(c)) continue;
                const double speed = std::hypot(state.coupledEastWindMps[c], state.coupledSouthWindMps[c]);
                stressU[c] = config.airDensityKgM3 * config.dragCoefficient * speed * state.coupledEastWindMps[c];
                stressV[c] = config.airDensityKgM3 * config.dragCoefficient * speed * state.coupledSouthWindMps[c];
                const double f = 2.0 * config.rotationRatePerSecond * config.rotationDirection * std::sin(grid.latitudeCentresRadians[y]);
                // The regularization scale is the slab's linear damping r.
                // In east/south coordinates: r U + f V = tau_east/rho,
                // r V - f U = tau_south/rho. Retain downwind drag as well as
                // perpendicular transport so the f=0 balance and wind work close.
                const double r = config.minimumCoriolisPerSecond;
                const double inverse = 1.0 / (config.waterDensityKgM3 * (f * f + r * r));
                ekmanU[c] = (r * stressU[c] - f * stressV[c]) * inverse;
                ekmanV[c] = (f * stressU[c] + r * stressV[c]) * inverse;
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
                const auto diffusion = [&](int xx, int yy, double spacing, std::size_t left, std::size_t right)
                {
                    if (config.coastalScheme != CoastalScheme::Legacy)
                    {
                        // This psi edge carries transport through the face
                        // between left/right cells. Use that same H in r/H,
                        // including edges next to a fixed coastal psi vertex.
                        const double inv = 1.0 / (0.5 * (depth(left) + depth(right)));
                        return (config.barotropicDragPerSecond * inv + config.linearBottomDragMps * inv * inv) /
                            (spacing * spacing);
                    }
                    const auto vv = vertex(xx, yy);
                    // curl(r u): with transport streamfunction, the elliptic
                    // coefficient is r/H = r0/H + r_b/H^2. Average that full
                    // positive coefficient at each shared edge.
                    const auto coefficient = [&](std::size_t corner)
                    {
                        const double inv = inverseDepth[corner];
                        return config.barotropicDragPerSecond * inv + config.linearBottomDragMps * inv * inv;
                    };
                    const double value = active[vv] ? 0.5 * (coefficient(v) + coefficient(vv)) : coefficient(v);
                    return value / (spacing * spacing);
                };
                aE[v] = diffusion(x + 1, y, dx, ne, se) + std::max(0.0, betaNorth) / dx;
                aW[v] = diffusion(x - 1, y, dx, nw, sw) + std::max(0.0, -betaNorth) / dx;
                // Curl and drag diffusion share the spherical meridional
                // metric: (1/cos(phi)) d/dphi(cos(phi) ...). Zonal stress
                // proportional to sec(phi) therefore has zero curl.
                const double cosHere = std::cos(grid.latitudeNorthFacesRadians[y]);
                const double cosNorth = std::cos(grid.latitudeCentresRadians[y - 1]);
                const double cosSouth = std::cos(grid.latitudeCentresRadians[y]);
                aS[v] = diffusion(x, y + 1, dy, sw, se) * cosSouth / cosHere + std::max(0.0, betaEast) / dy;
                aN[v] = diffusion(x, y - 1, dy, nw, ne) * cosNorth / cosHere + std::max(0.0, -betaEast) / dy;
                diagonal[v] = aE[v] + aW[v] + aN[v] + aS[v];
                rhs[v] = -((stressV[ne] / depth(ne) + stressV[se] / depth(se) - stressV[nw] / depth(nw) - stressV[sw] / depth(sw)) / (2.0 * dx) -
                    (cosSouth * (stressU[sw] / depth(sw) + stressU[se] / depth(se)) -
                        cosNorth * (stressU[nw] / depth(nw) + stressU[ne] / depth(ne))) / (2.0 * dy * cosHere)) / config.waterDensityKgM3;
                if (config.coastalScheme == CoastalScheme::FaceForcing)
                {
                    const auto acceleration = [&](const auto& stress, std::size_t left, std::size_t right)
                    {
                        return (stress[left] + stress[right]) / (depth(left) + depth(right));
                    };
                    rhs[v] = -((acceleration(stressV, ne, se) - acceleration(stressV, nw, sw)) / dx -
                        (cosSouth * acceleration(stressU, sw, se) - cosNorth * acceleration(stressU, nw, ne)) /
                        (dy * cosHere)) / config.waterDensityKgM3;
                }
                rhsNorm += rhs[v] * rhs[v];
            }
        if (config.krylovCirculation)
        {
            const detail::CirculationStencil matrix{columns, active, diagonal, aE, aW, aN, aS};
            const auto solved = detail::solveCirculationKrylov(matrix, rhs, psi,
                config.streamfunctionIterations, config.streamfunctionTolerance);
            state.streamfunctionRelativeResidual = static_cast<float>(solved.relativeResidual);
        }
        else for (int sweep = 0; sweep < std::max(1, config.streamfunctionIterations); ++sweep)
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
                    eastBarotropic[c] = state.eastVolumeTransportM3S[c] / faceDepth;
                    eastHeatTransport[c] = state.eastVolumeTransportM3S[c] / faceDepth +
                        0.5 * (ekmanU[c] + ekmanU[e]) * grid.zonalFaceLengthsMetres[y] / config.mixedLayerDepthMetres;
                    eastDiffusion[c] = config.heatDiffusivityM2S * grid.zonalFaceLengthsMetres[y] * grid.zonalFaceLengthsMetres[y] /
                        grid.cellAreasSquareMetres[y];
                }
                if (y + 1 < rows && ocean(c) && ocean(s))
                {
                    const double faceDepth = 0.5 * (depth(c) + depth(s));
                    southBarotropic[c] = state.southVolumeTransportM3S[c] / faceDepth;
                    southHeatTransport[c] = state.southVolumeTransportM3S[c] / faceDepth +
                        0.5 * (ekmanV[c] + ekmanV[s]) * grid.southFaceLengthsMetres[y] / config.mixedLayerDepthMetres;
                    southDiffusion[c] = config.heatDiffusivityM2S * grid.southFaceLengthsMetres[y] / dy;
                }
            }
        state.maximumTransportDivergenceMps = 0.0;
        auto& distribution = state.stabilityDistribution;
        std::vector<double> localRates;
        localRates.reserve(count);
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
                auto& components = state.transportComponents[c];
                const double northBarotropic = y > 0 ? southBarotropic[n] : 0.0;
                components.barotropicEastMps = 0.5 * (eastBarotropic[c] + eastBarotropic[w]) / grid.zonalFaceLengthsMetres[y];
                components.barotropicSouthMps = (southBarotropic[c] + northBarotropic) /
                    std::max(1.0, grid.northFaceLengthsMetres[y] + grid.southFaceLengthsMetres[y]);
                components.ekmanEastMps = 0.5 * (eastHeatTransport[c] - eastBarotropic[c] +
                    eastHeatTransport[w] - eastBarotropic[w]) / grid.zonalFaceLengthsMetres[y];
                components.ekmanSouthMps = (southHeatTransport[c] - southBarotropic[c] + northHeat - northBarotropic) /
                    std::max(1.0, grid.northFaceLengthsMetres[y] + grid.southFaceLengthsMetres[y]);
                components.barotropicUpwellingMps = (eastBarotropic[c] - eastBarotropic[w] +
                    southBarotropic[c] - northBarotropic) * config.mixedLayerDepthMetres / area;
                components.ekmanUpwellingMps = divergence[c] * config.mixedLayerDepthMetres - components.barotropicUpwellingMps;
                const double storageRateScale = config.variableHeatStorage
                    ? config.mixedLayerDepthMetres / detail::minimumStorageDepth(columnDepth[c], config) : 1.0;
                const double outgoing = std::max(0.0, eastHeatTransport[c]) + std::max(0.0, -eastHeatTransport[w]) +
                    std::max(0.0, southHeatTransport[c]) + std::max(0.0, -northHeat);
                const double diffusionRate = eastDiffusion[c] + eastDiffusion[w] + southDiffusion[c] + (y > 0 ? southDiffusion[n] : 0.0);
                maximumExplicitRate = std::max(maximumExplicitRate, storageRateScale * ((outgoing + diffusionRate) / area + std::max(0.0, -divergence[c]) +
                    config.surfaceHeatExchangeWm2K / heatCapacity));
                const double explicitDiffusion = config.implicitZonalHeatDiffusion
                    ? southDiffusion[c] + (y > 0 ? southDiffusion[n] : 0.0) : diffusionRate;
                const double rate = storageRateScale * ((outgoing + explicitDiffusion) / area + std::max(0.0, -divergence[c]) +
                    config.surfaceHeatExchangeWm2K / heatCapacity);
                ++distribution.wetCells; distribution.wetAreaM2 += area;
                localRates.push_back(rate);
                constexpr std::array<int, 5> hours{1, 3, 6, 12, 24};
                for (std::size_t bin = 0; bin < hours.size(); ++bin)
                    if (rate * hours[bin] * 3600.0 > 0.7)
                    {
                        ++distribution.restrictedCells[bin];
                        distribution.restrictedAreaM2[bin] += area;
                    }
                distribution.maximumRateWithoutAdvection = std::max(distribution.maximumRateWithoutAdvection,
                    storageRateScale * (explicitDiffusion / area + std::max(0.0, -divergence[c]) + config.surfaceHeatExchangeWm2K / heatCapacity));
                if (rate > distribution.highestRates.back().totalRate || distribution.highestRates.back().cell < 0)
                {
                    OceanStabilityLimit limit;
                    limit.cell = static_cast<int>(c); limit.totalRate = rate;
                    limit.zonalAdvectionRate = storageRateScale * ((std::max(0.0, eastHeatTransport[c]) + std::max(0.0, -eastHeatTransport[w])) / area);
                    limit.meridionalAdvectionRate = storageRateScale * ((std::max(0.0, southHeatTransport[c]) + std::max(0.0, -northHeat)) / area);
                    limit.zonalDiffusionRate = storageRateScale * ((eastDiffusion[c] + eastDiffusion[w]) / area);
                    limit.meridionalDiffusionRate = storageRateScale * ((southDiffusion[c] + (y > 0 ? southDiffusion[n] : 0.0)) / area);
                    limit.verticalRate = storageRateScale * (std::max(0.0, -divergence[c]));
                    limit.surfaceRate = storageRateScale * (config.surfaceHeatExchangeWm2K / heatCapacity);
                    limit.depthMetres = depth(c);
                    const std::array<std::size_t, 4> neighbours{grid.index(x + 1, y), w, grid.index(x, y + 1), n};
                    const std::array<double, 4> lengths{grid.zonalFaceLengthsMetres[y], grid.zonalFaceLengthsMetres[y],
                        grid.southFaceLengthsMetres[y], grid.northFaceLengthsMetres[y]};
                    const std::array<double, 4> transports{eastHeatTransport[c], eastHeatTransport[w], southHeatTransport[c], northHeat};
                    const std::array<double, 4> barotropic{eastBarotropic[c], eastBarotropic[w], southBarotropic[c], northBarotropic};
                    for (std::size_t face = 0; face < 4; ++face)
                    {
                        // index() wraps longitude only. Polar neighbours lie
                        // outside the array; reject closed faces BEFORE reading.
                        if ((face == 2 && y + 1 == rows) || (face == 3 && y == 0)) continue;
                        if (lengths[face] > 0.0 && !forcing.landMask.at(neighbours[face]) &&
                            forcing.bathymetryMetres.at(neighbours[face]) > 0.0f)
                        {
                            limit.faceVelocityMps[face] = transports[face] / lengths[face];
                            limit.barotropicFaceVelocityMps[face] = barotropic[face] / lengths[face];
                            limit.faceDepthMetres[face] = 0.5 * (depth(c) + depth(neighbours[face]));
                        }
                    }
                    if (rate > maximumRate) state.stabilityLimit = limit;
                    auto position = std::find_if(distribution.highestRates.begin(), distribution.highestRates.end(),
                        [&](const auto& value) { return value.cell < 0 || rate > value.totalRate; });
                    if (position != distribution.highestRates.end())
                    {
                        std::move_backward(position, distribution.highestRates.end() - 1, distribution.highestRates.end());
                        *position = limit;
                    }
                }
                maximumRate = std::max(maximumRate, rate);
            }
        std::sort(localRates.begin(), localRates.end());
        if (!localRates.empty())
        {
            constexpr std::array<double, 3> quantiles{0.50, 0.95, 0.99};
            for (std::size_t q = 0; q < quantiles.size(); ++q)
                distribution.rateQuantiles[q] = localRates[static_cast<std::size_t>(quantiles[q] * (localRates.size() - 1))];
        }
        if (prepared && config.oneWay)
        {
            copyCirculation(prepared->diagnostics, state);
            prepared->psi = psi;
            prepared->eastHeat = eastHeatTransport; prepared->southHeat = southHeatTransport;
            prepared->eastDiffusion = eastDiffusion; prepared->southDiffusion = southDiffusion;
            prepared->divergence = divergence;
            prepared->maximumRate = maximumRate; prepared->maximumExplicitRate = maximumExplicitRate;
            prepared->ready = true;
        }
        }
        const auto heatStart = std::chrono::steady_clock::now();
        state.circulationSeconds += std::chrono::duration<double>(heatStart - circulationStart).count();
        // Each outer iteration solves the SAME finite seasonal interval from
        // T0, not a progressively longer integration mistaken for convergence.
        std::vector<double> temperature(forcing.initialSstC.begin(), forcing.initialSstC.end());
        auto enthalpy = initialEnthalpy;
        auto storageDepth = initialStorageDepth, reservoirEnthalpy = initialReservoirEnthalpy;
        for (std::size_t c = 0; c < count; ++c)
            capacities[c] = config.variableHeatStorage && ocean(c) ? volumetricCapacity * storageDepth[c] : heatCapacity;
        state.meanStorageDepthMetres.assign(count, 0.0f);
        for (std::size_t c = 0; c < count; ++c)
            if (ocean(c)) temperature[c] = freezing + std::max(0.0, enthalpy[c]) / capacities[c];
        const double interval = config.oceanTimeStepSeconds * std::max(1, config.heatStepsPerIteration);
        const double explicitSteps = std::max<double>(std::max(1, config.heatStepsPerIteration),
            std::ceil(interval * maximumExplicitRate / 0.7));
        if (!std::isfinite(explicitSteps) || explicitSteps > std::numeric_limits<int>::max())
        {
            state.finite = false; state.converged = false;
            return state;
        }
        state.explicitHeatSubsteps = static_cast<int>(explicitSteps);
        const int substeps = std::max(std::max(1, config.heatStepsPerIteration),
            static_cast<int>(std::ceil(interval * maximumRate / 0.7)));
        state.heatSubsteps = substeps;
        const double dt = interval / substeps;
        state.actualHeatStepSeconds = dt;
        state.meanSstC.assign(count, 0.0f);
        state.meanSurfaceSkinTemperatureC.assign(count, 0.0f);
        state.meanIceCover.assign(count, 0.0f);
        state.heatBudget.assign(count, {});
        for (std::size_t c = 0; c < count; ++c) if (ocean(c))
        {
            state.heatBudget[c].depthMetres = depth(c);
            state.heatBudget[c].deepTemperatureC = deepTemperature[c];
        }
        const auto accumulateMean = [&](std::size_t c, double skin)
        {
            const double weight = 0.5 / substeps;
            state.meanSstC[c] += static_cast<float>(weight * temperature[c]);
            state.meanStorageDepthMetres[c] += static_cast<float>(weight * storageDepth[c]);
            state.meanSurfaceSkinTemperatureC[c] += static_cast<float>(weight * skin);
            state.meanIceCover[c] += static_cast<float>(weight * std::clamp(-enthalpy[c] / (latentHeatPerMetre * 0.10), 0.0, 1.0));
        };
        double expectedHeatChange = 0.0, absoluteHeatExchange = 0.0;
        for (int step = 0; step < substeps; ++step)
        {
            const auto fluxStart = std::chrono::steady_clock::now();
            if (config.variableHeatStorage) heatRows([&](int firstRow, int lastRow)
            {
                for (int y = firstRow; y <= lastRow; ++y) for (int x = 0; x < columns; ++x)
                {
                    const auto c = static_cast<std::size_t>(y) * columns + x;
                    if (!ocean(c)) continue;
                    const double flux = surfaceExchange(c, enthalpy[c], temperature[c]).second;
                    const double wind = std::hypot(forcing.eastWindMps[c], forcing.southWindMps[c]);
                    const double lowerT = freezing + reservoirEnthalpy[c] /
                        (volumetricCapacity * (columnDepth[c] - storageDepth[c]));
                    const double stratification = detail::storageThermalStratification(temperature[c], lowerT, columnDepth[c], config);
                    const double target = detail::storageTargetDepth(wind, flux, columnDepth[c], config, stratification);
                    const double next = storageDepth[c] + (target - storageDepth[c]) *
                        (-std::expm1(-dt / (config.storageAdjustmentDays * 86400.0)));
                    const double exchange = detail::remixStorage(next, columnDepth[c], volumetricCapacity,
                        storageDepth[c], enthalpy[c], reservoirEnthalpy[c]);
                    state.heatBudget[c].entrainmentWm2 += exchange / interval;
                    state.heatBudget[c].verticalWm2 += exchange / interval;
                    capacities[c] = volumetricCapacity * storageDepth[c];
                    temperature[c] = freezing + std::max(0.0, enthalpy[c]) / capacities[c];
                }
            });
            if (!config.gatherHeatFluxes)
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
                    flux(e, eastHeatTransport[c], config.implicitZonalHeatDiffusion ? 0.0 : eastDiffusion[c]);
                    if (y + 1 < rows) flux(s, southHeatTransport[c], southDiffusion[c]);
                }
            }
            else heatRows([&](int firstRow, int lastRow)
            {
                for (int y = firstRow; y <= lastRow; ++y)
                    for (int x = 0; x < columns; ++x)
                    {
                        const auto c = static_cast<std::size_t>(y) * columns + x;
                        const auto e = x + 1 < columns ? c + 1 : c + 1 - columns;
                        const auto flux = [&](std::size_t neighbour, double transport, double diffusion)
                        {
                            return transport * (transport >= 0.0 ? temperature[c] : temperature[neighbour]) +
                                diffusion * (temperature[c] - temperature[neighbour]);
                        };
                        eastFlux[c] = flux(e, eastHeatTransport[c], config.implicitZonalHeatDiffusion ? 0.0 : eastDiffusion[c]);
                        southFlux[c] = y + 1 < rows ? flux(c + columns, southHeatTransport[c], southDiffusion[c]) : 0.0;
                    }
            });
            const auto updateStart = std::chrono::steady_clock::now();
            state.heatFluxSeconds += std::chrono::duration<double>(updateStart - fluxStart).count();
            const auto updateRows = [&](int firstRow, int lastRow)
            {
            for (int y = firstRow; y <= lastRow; ++y)
            {
                rowExpected[y] = rowAbsolute[y] = 0.0;
                for (int x = 0; x < columns; ++x)
                {
                    const auto c = static_cast<std::size_t>(y) * columns + x;
                    if (!ocean(c)) continue;
                    if (config.gatherHeatFluxes)
                    {
                        // Reproduce the original face-scatter addition order,
                        // including the final west-face contribution at x=0.
                        double value = 0.0;
                        if (y > 0) value += southFlux[c - columns];
                        if (x > 0) value += eastFlux[c - 1];
                        value -= eastFlux[c];
                        if (y + 1 < rows) value -= southFlux[c];
                        if (x == 0) value += eastFlux[c + columns - 1];
                        tendency[c] = value;
                    }
                    const double area = grid.cellAreasSquareMetres[y];
                    // Fixed-depth continuity closes horizontal convergence via
                    // vertical exchange with the deep reservoir. Uniform T is
                    // invariant when the deep contrast and heat exchange vanish.
                    const double deepT = deepTemperature[c];
                    const auto surface = surfaceExchange(c, enthalpy[c], temperature[c]);
                    const double surfaceSource = surface.second;
                    const double source = surfaceSource / heatCapacity +
                        divergence[c] * (divergence[c] > 0.0 ? deepT : temperature[c]);
                    auto& budget = state.heatBudget[c];
                    budget.horizontalWm2 += heatCapacity * tendency[c] / (area * substeps);
                    budget.surfaceWm2 += surfaceSource / substeps;
                    budget.verticalWm2 += heatCapacity * divergence[c] *
                        (divergence[c] > 0.0 ? deepT : temperature[c]) / substeps;
                    if (config.gatherHeatFluxes)
                    {
                        rowExpected[y] += dt * area * heatCapacity * source;
                        rowAbsolute[y] += std::abs(dt * area * heatCapacity * source);
                    }
                    else
                    {
                        expectedHeatChange += dt * area * heatCapacity * source;
                        absoluteHeatExchange += std::abs(dt * area * heatCapacity * source);
                    }
                    accumulateMean(c, surface.first);
                    enthalpy[c] += dt * heatCapacity * (tendency[c] / area + source);
                    temperature[c] = freezing + std::max(0.0, enthalpy[c]) / capacities[c];
                    if (!config.implicitZonalHeatDiffusion)
                        accumulateMean(c, surfaceExchange(c, enthalpy[c], temperature[c]).first);
                }
            }
            };
            if (config.gatherHeatFluxes)
            {
                heatRows(updateRows);
                for (int y = 0; y < rows; ++y)
                {
                    expectedHeatChange += rowExpected[y];
                    absoluteHeatExchange += rowAbsolute[y];
                }
            }
            else updateRows(0, rows - 1);
            state.heatUpdateSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - updateStart).count();
            if (config.implicitZonalHeatDiffusion)
            {
                const auto zonalStart = std::chrono::steady_clock::now();
                parallelforrows(0, rows - 1, [&](int firstRow, int lastRow)
                {
                thread_local detail::ZonalDiffusionWorkspace zonalDiffusion(0);
                if (zonalDiffusion.theta.size() != static_cast<std::size_t>(columns))
                    zonalDiffusion = detail::ZonalDiffusionWorkspace(columns);
                for (int y = firstRow; y <= lastRow; ++y)
                {
                    const double area = grid.cellAreasSquareMetres[y];
                    const auto first = grid.index(0, y);
                    if (!zonalDiffusion.solve(enthalpy.data() + first, eastDiffusion.data() + first,
                        heatCapacity, dt / area, rowIterations[y], rowResidual[y],
                        config.variableHeatStorage ? capacities.data() + first : nullptr))
                    {
                        rowValid[y] = 0;
                        return;
                    }
                    // Apply each solved face flux once with opposite signs.
                    // This conserves enthalpy, including latent heat, and records
                    // the same backward-Euler flux used to advance the state.
                    for (int x = 0; x < columns; ++x)
                    {
                        const int east = (x + 1) % columns;
                        const auto c = first + x, e = first + east;
                        const double heat = heatCapacity * eastDiffusion[c] / area *
                            (zonalDiffusion.theta[x] - zonalDiffusion.theta[east]);
                        enthalpy[c] -= dt * heat; enthalpy[e] += dt * heat;
                        state.heatBudget[c].horizontalWm2 -= heat / substeps;
                        state.heatBudget[e].horizontalWm2 += heat / substeps;
                    }
                    for (int x = 0; x < columns; ++x)
                    {
                        const auto c = first + x;
                        if (!ocean(c)) continue;
                        temperature[c] = freezing + std::max(0.0, enthalpy[c]) / capacities[c];
                        accumulateMean(c, surfaceExchange(c, enthalpy[c], temperature[c]).first);
                    }
                }
                }, config.parallelZonalHeatDiffusion ? rowsPerWorker : rows);
                state.zonalHeatSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - zonalStart).count();
                state.maximumZonalDiffusionIterations = *std::max_element(rowIterations.begin(), rowIterations.end());
                state.maximumZonalDiffusionResidualK = *std::max_element(rowResidual.begin(), rowResidual.end());
                if (std::find(rowValid.begin(), rowValid.end(), 0) != rowValid.end())
                {
                    state.finite = false; state.converged = false;
                    return state;
                }
            }
        }
        double heatChange = 0.0, residual = 0.0, totalArea = 0.0;
        const float relaxation = config.oneWay ? 1.0f : std::clamp(config.underRelaxation, 0.01f, 1.0f);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const auto c = grid.index(x, y);
                if (!ocean(c)) continue;
                const double area = grid.cellAreasSquareMetres[y];
                heatChange += area * (enthalpy[c] - initialEnthalpy[c] + reservoirEnthalpy[c] - initialReservoirEnthalpy[c]);
                state.heatBudget[c].columnStorageWm2 = (enthalpy[c] - initialEnthalpy[c] +
                    reservoirEnthalpy[c] - initialReservoirEnthalpy[c]) / interval;
                state.heatBudget[c].storageWm2 = (enthalpy[c] - initialEnthalpy[c]) / interval;
                // Include ice-energy changes in convergence even while SST is
                // pinned at freezing; relax energy, not temperature alone.
                const double delta = (enthalpy[c] - coupledEnthalpy[c]) / capacities[c];
                residual += area * delta * delta; // normalized by 1 K below
                totalArea += area;
                coupledEnthalpy[c] += relaxation * (enthalpy[c] - coupledEnthalpy[c]);
                state.sstC[c] = static_cast<float>(freezing + std::max(0.0, coupledEnthalpy[c]) / capacities[c]);
                state.iceThicknessMetres[c] = static_cast<float>(std::max(0.0, -coupledEnthalpy[c]) / latentHeatPerMetre);
                state.surfaceSkinTemperatureC[c] = static_cast<float>(surfaceExchange(c, coupledEnthalpy[c], state.sstC[c]).first);
            }
        state.storageDepthMetres.resize(count);
        std::transform(storageDepth.begin(), storageDepth.end(), state.storageDepthMetres.begin(),
            [](double h) { return static_cast<float>(h); });
        state.reservoirTemperatureC.assign(count, static_cast<float>(freezing));
        for (std::size_t c = 0; c < count; ++c) if (config.variableHeatStorage && ocean(c))
            state.reservoirTemperatureC[c] = static_cast<float>(freezing + reservoirEnthalpy[c] /
                (volumetricCapacity * (columnDepth[c] - storageDepth[c])));
        state.heatBudgetResidualJ = heatChange - expectedHeatChange;
        state.relativeHeatBudgetResidual = state.heatBudgetResidualJ /
            std::max({1.0, std::abs(heatChange), absoluteHeatExchange});
        if (!config.oneWay)
            legacy::applyEmpiricalSstWindFeedback(grid, forcing, config, relaxation, state, residual);
        state.relativeResidual = totalArea > 0.0 ? static_cast<float>(std::sqrt(residual / totalArea)) : 0.0f;
        state.residualHistory.push_back(state.relativeResidual);
        state.finite = valid(state.sstC) && valid(state.coupledEastWindMps) && valid(state.coupledSouthWindMps) &&
            valid(state.iceThicknessMetres) && valid(state.surfaceSkinTemperatureC) &&
            std::isfinite(state.streamfunctionRelativeResidual) && valid(state.storageDepthMetres) &&
            valid(state.reservoirTemperatureC);
        state.converged = state.finite && state.streamfunctionRelativeResidual <= config.streamfunctionTolerance &&
            (config.oneWay || state.relativeResidual <= config.convergenceTolerance);
        state.heatSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - heatStart).count();
        if (!state.finite || state.converged) break;
    }
    state.streamfunctionM3S.resize(vertices);
    std::transform(psi.begin(), psi.end(), state.streamfunctionM3S.begin(), [](double v) { return static_cast<float>(v); });
    return state;
}
}
}
