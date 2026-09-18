#include "climate_ocean_dynamics.hpp"
#include "climate_grid.hpp"
#include "detail/ocean_mixed_layer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

namespace
{
int failures = 0;
void expect(bool condition, const char* message)
{
    if (!condition) { std::cerr << "FAIL: " << message << '\n'; ++failures; }
}

climateocean::OceanConfig config()
{
    climateocean::OceanConfig c;
    c.oneWay = true;
    c.surfaceHeatExchangeWm2K = 0.0f;
    c.deepWaterTemperatureContrastK = 0.0f;
    c.heatStepsPerIteration = 1;
    return c;
}

climateocean::OceanForcing fixture(int columns, bool polarRings)
{
    const int rows = columns / 2, count = columns * rows;
    climateocean::OceanForcing f;
    f.landMask.assign(count, 0); f.bathymetryMetres.assign(count, 4000.0f);
    f.eastWindMps.assign(count, 0.0f); f.southWindMps.assign(count, 0.0f);
    f.initialSstC.assign(count, 10.0f); f.atmosphericTemperatureC.assign(count, 10.0f);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const int i = y * columns + x;
            if (polarRings && y != 0 && y != rows - 1)
            { f.landMask[i] = 1; f.bathymetryMetres[i] = 0.0f; }
        }
    return f;
}

double difference(const std::vector<float>& a, const std::vector<float>& b)
{
    double value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) value = std::max(value, std::abs(double(a[i]) - b[i]));
    return value;
}

double budgetError(const climateocean::OceanState& s)
{
    double error = 0.0;
    for (const auto& b : s.heatBudget)
        error = std::max(error, std::abs(b.storageWm2 - b.horizontalWm2 - b.surfaceWm2 - b.verticalWm2));
    return error;
}

bool sameHeatState(const climateocean::OceanState& a, const climateocean::OceanState& b)
{
    if (a.sstC != b.sstC || a.meanSstC != b.meanSstC || a.iceThicknessMetres != b.iceThicknessMetres ||
        a.meanIceCover != b.meanIceCover || a.surfaceSkinTemperatureC != b.surfaceSkinTemperatureC ||
        a.meanSurfaceSkinTemperatureC != b.meanSurfaceSkinTemperatureC || a.heatSubsteps != b.heatSubsteps ||
        a.heatBudget.size() != b.heatBudget.size()) return false;
    for (std::size_t i = 0; i < a.heatBudget.size(); ++i)
    {
        const auto& x = a.heatBudget[i]; const auto& y = b.heatBudget[i];
        if (x.horizontalWm2 != y.horizontalWm2 || x.surfaceWm2 != y.surfaceWm2 ||
            x.verticalWm2 != y.verticalWm2 || x.storageWm2 != y.storageWm2 ||
            x.depthMetres != y.depthMetres || x.deepTemperatureC != y.deepTemperatureC) return false;
    }
    return true;
}

climateocean::OceanForcing heatFixture(int width)
{
    auto f = fixture(width, false);
    f.initialIceThicknessMetres.assign(width * width / 2, 0.0f);
    f.surfaceHeatFluxWm2.assign(width * width / 2, 12.0f);
    f.surfaceHeatFluxReferenceTemperatureC.assign(width * width / 2, 8.0f);
    for (int y = 0; y < width / 2; ++y)
        for (int x = 0; x < width; ++x)
        {
            const int i = y * width + x;
            f.landMask[i] = (x > width / 3 && x < width / 2 && y > width / 8 && y < width / 3);
            f.bathymetryMetres[i] = 200.0f + 3000.0f * x / width;
            f.eastWindMps[i] = static_cast<float>(6.0 * std::sin(12.0 * y / width));
            f.southWindMps[i] = static_cast<float>(3.0 * std::cos(9.0 * x / width));
            f.initialSstC[i] = static_cast<float>(8.0 + 5.0 * std::sin(20.0 * x / width));
            if (y < 2 || y >= width / 2 - 2)
            {
                f.initialSstC[i] = -1.8f;
                f.initialIceThicknessMetres[i] = 0.05f + 0.2f * (x % 3);
            }
        }
    return f;
}

void heatKernelTests()
{
    constexpr int width = 64;
    const auto f = heatFixture(width);
    for (bool implicit : {false, true}) for (float rotation : {-1.0f, 1.0f})
    {
        auto c = config();
        c.surfaceHeatExchangeWm2K = 18.0f; c.deepWaterTemperatureContrastK = 4.0f;
        c.oceanTimeStepSeconds = 3600.0f; c.heatStepsPerIteration = 48;
        c.implicitZonalHeatDiffusion = implicit; c.rotationDirection = rotation;
        c.gatherHeatFluxes = false;
        const auto reference = climateocean::solveWindDrivenOcean(width, width / 2, f, c);
        expect(reference.converged && budgetError(reference) < 1.0e-5, "heat kernel reference must conserve heat with coasts, ice and both rotations");
        c.gatherHeatFluxes = true;
        for (int workers : {1, 2, 4, 8})
        {
            c.heatWorkers = workers;
            const auto state = climateocean::solveWindDrivenOcean(width, width / 2, f, c);
            expect(state.converged && sameHeatState(reference, state) && budgetError(state) < 1.0e-5 &&
                std::abs(state.relativeHeatBudgetResidual) < 1.0e-11,
                "gather and parallel kernels must reproduce all heat fields and cell budgets exactly");
            const auto& d = state.stabilityDistribution;
            bool valid = d.wetCells > 0 && d.wetAreaM2 > 0.0 && d.highestRates[0].cell == state.stabilityLimit.cell &&
                d.highestRates[0].totalRate == state.stabilityLimit.totalRate && d.maximumRateWithoutAdvection <= state.stabilityLimit.totalRate;
            for (std::size_t i = 0; i < d.restrictedCells.size(); ++i)
                valid = valid && d.restrictedCells[i] <= d.wetCells && d.restrictedAreaM2[i] <= d.wetAreaM2 &&
                    (!i || (d.restrictedCells[i] >= d.restrictedCells[i - 1] && d.restrictedAreaM2[i] >= d.restrictedAreaM2[i - 1]));
            for (std::size_t i = 1; i < d.highestRates.size(); ++i)
                valid = valid && d.highestRates[i].totalRate <= d.highestRates[i - 1].totalRate;
            for (const auto& hot : d.highestRates)
                if (hot.cell >= 0)
                {
                    valid = valid && !f.landMask[hot.cell] &&
                        std::abs(0.5 * (hot.faceVelocityMps[0] + hot.faceVelocityMps[1]) - state.eastCurrentMps[hot.cell]) < 1.0e-5;
                    if (hot.cell < width) valid = valid && hot.faceVelocityMps[3] == 0.0 && hot.faceDepthMetres[3] == 0.0;
                    if (hot.cell >= width * (width / 2 - 1)) valid = valid && hot.faceVelocityMps[2] == 0.0 && hot.faceDepthMetres[2] == 0.0;
                }
            expect(valid, "ocean rate distribution and face velocities must agree with the maximum and cell current diagnostics");
        }
    }
}

void variableStorageTests()
{
    auto c = config(); c.variableHeatStorage = true;
    const double capacity = double(c.waterDensityKgM3) * c.waterHeatCapacityJkgK;
    double h = 30.0, e = capacity * h * 20.0, r = capacity * (300.0 - h) * 10.0;
    const double total = e + r;
    for (double next : {120.0, 15.0, 80.0, 10.0, 190.0, 30.0})
    {
        const double before = e;
        const double transfer = climateocean::detail::remixStorage(next, 300.0, capacity, h, e, r);
        expect(std::abs((e + r) / total - 1.0) < 1.0e-14 && std::abs(e - before - transfer) < 1.0e-5 &&
            r >= 0.0 && h == next, "entrainment/detrainment must conserve column heat and retain the exchanged water");
    }
    h = 30.0; e = -1.0e8; r = 2.0e8;
    climateocean::detail::remixStorage(60.0, 300.0, capacity, h, e, r);
    expect(e >= 0.0 && std::abs(e + r - 1.0e8) < 1.0e-7,
        "warm entrainment must melt ice through conserved latent enthalpy");
    const auto target = [&](double wind, double flux)
    { return climateocean::detail::storageTargetDepth(wind, flux, 300.0, c); };
    expect(target(5.0, -100.0) > target(5.0, 100.0) && target(10.0, 0.0) > target(5.0, 0.0),
        "cooling and stronger wind must increase diagnosed storage depth");
    const double weak = climateocean::detail::storageThermalStratification(10.0, 10.0, 300.0, c);
    const double strong = climateocean::detail::storageThermalStratification(20.0, 10.0, 300.0, c);
    expect(weak == c.storageStratificationPerSecond2 && strong > weak &&
        climateocean::detail::storageTargetDepth(5.0, -100.0, 300.0, c, weak) >
        climateocean::detail::storageTargetDepth(5.0, -100.0, 300.0, c, strong),
        "a stronger retained thermocline must resist mixing relative to an isothermal winter column");

    auto f = fixture(32, false);
    f.initialStorageDepthMetres.resize(512);
    f.initialReservoirTemperatureC.assign(512, 10.0f);
    for (int i = 0; i < 512; ++i)
    {
        f.initialStorageDepthMetres[i] = 12.0f + (i % 80);
        f.eastWindMps[i] = 3.0f;
    }
    c.heatStepsPerIteration = 4;
    for (bool implicit : {false, true})
    {
        c.implicitZonalHeatDiffusion = implicit;
        const auto s = climateocean::solveWindDrivenOcean(32, 16, f, c);
        expect(s.converged && difference(s.sstC, f.initialSstC) < 1.0e-6 &&
            difference(s.reservoirTemperatureC, f.initialReservoirTemperatureC) < 1.0e-6 && budgetError(s) < 1.0e-5 &&
            std::abs(s.relativeHeatBudgetResidual) < 1.0e-10,
            "variable capacity, moving depth and volume divergence must preserve isothermal water and close heat budgets");
    }
    f = heatFixture(32);
    c.implicitZonalHeatDiffusion = true;
    c.surfaceHeatExchangeWm2K = 18.0f; c.deepWaterTemperatureContrastK = 4.0f;
    const auto implicit = climateocean::solveWindDrivenOcean(32, 16, f, c);
    expect(implicit.converged && budgetError(implicit) < 1.0e-5 && std::abs(implicit.relativeHeatBudgetResidual) < 1.0e-10,
        "variable-capacity implicit diffusion must conserve with ice, coasts, surface flux and subsurface exchange");
    for (const auto& b : implicit.heatBudget)
        expect(std::abs(b.columnStorageWm2 - b.horizontalWm2 - b.surfaceWm2 - b.verticalWm2 + b.entrainmentWm2) < 1.0e-5,
            "whole-column storage must exclude internal entrainment from the external heat budget");
    c.oceanTimeStepSeconds = 600.0f; c.heatStepsPerIteration = 576;
    const auto fineImplicit = climateocean::solveWindDrivenOcean(32, 16, f, c);
    c.implicitZonalHeatDiffusion = false;
    const auto fineExplicit = climateocean::solveWindDrivenOcean(32, 16, f, c);
    expect(fineImplicit.converged && fineExplicit.converged && difference(fineImplicit.sstC, fineExplicit.sstC) < 0.002,
        "explicit and variable-capacity implicit heat solves must agree at a refined step");
    c = config(); c.heatStepsPerIteration = 3;
    f = fixture(32, false); f.atmosphericTemperatureC.assign(512, 20.0f);
    c.surfaceHeatExchangeWm2K = 18.0f; c.deepWaterTemperatureContrastK = 4.0f;
    const auto legacy = climateocean::solveWindDrivenOcean(32, 16, f, c);
    c.variableHeatStorage = true;
    c.minimumStorageDepthMetres = c.maximumStorageDepthMetres = c.mixedLayerDepthMetres;
    const auto fixed = climateocean::solveWindDrivenOcean(32, 16, f, c);
    expect(legacy.sstC == fixed.sstC && legacy.eastCurrentMps == fixed.eastCurrentMps,
        "fixed storage depth must reproduce the existing fixed slab without changing momentum");
    f.initialStorageDepthMetres.assign(512, -1.0f);
    expect(!climateocean::solveWindDrivenOcean(32, 16, f, c).finite, "negative carried storage depths must be rejected");
}

void runTests()
{
    variableStorageTests();
    heatKernelTests();
    auto momentum = fixture(64, false);
    for (int y = 0; y < 32; ++y)
        for (int x = 0; x < 64; ++x)
        {
            const int i = y * 64 + x;
            momentum.landMask[i] = (x > 24 && x < 30 && y > 8 && y < 24);
            momentum.bathymetryMetres[i] = 200.0f + 50.0f * x;
            momentum.initialSstC[i] = 5.0f + 0.1f * x;
            momentum.eastWindMps[i] = static_cast<float>(3.0 * std::sin(0.15 * y));
            momentum.southWindMps[i] = static_cast<float>(std::cos(0.1 * x));
        }
    auto precise = config();
    precise.streamfunctionIterations = 16000;
    precise.streamfunctionTolerance = 1.0e-8f;
    precise.krylovCirculation = false;
    const auto relaxation = climateocean::solveWindDrivenOcean(64, 32, momentum, precise);
    precise.krylovCirculation = true;
    const auto krylov = climateocean::solveWindDrivenOcean(64, 32, momentum, precise);
    expect(relaxation.converged && krylov.converged &&
        difference(relaxation.eastCurrentMps, krylov.eastCurrentMps) < 1.0e-6 &&
        difference(relaxation.southCurrentMps, krylov.southCurrentMps) < 1.0e-6 &&
        difference(relaxation.sstC, krylov.sstC) < 2.0e-6 &&
        krylov.maximumTransportDivergenceMps < 1.0e-12 && budgetError(krylov) < 1.0e-5,
        "Krylov circulation must match tightly converged relaxation with bathymetry, island coasts and a periodic seam");
    // Seasonal circulation reuse must reproduce recomputation exactly through
    // changing thermal/ice state, and must not survive changed outer forcing.
    std::array<climateocean::OceanForcing, 4> seasons;
    std::array<climateocean::OceanConfig, 4> configs;
    for (int season = 0; season < 4; ++season)
    {
        seasons[season] = fixture(32, false);
        auto& forcing = seasons[season];
        forcing.initialIceThicknessMetres.assign(32 * 16, 0.0f);
        for (int i = 0; i < 32 * 16; ++i)
        {
            forcing.eastWindMps[i] = static_cast<float>((season + 1) * std::sin(0.2 * (i / 32)));
            forcing.bathymetryMetres[i] = 1000.0f + 10.0f * (i % 32);
            if (i < 32) { forcing.initialSstC[i] = -1.8f; forcing.initialIceThicknessMetres[i] = 0.5f; }
        }
        configs[season] = config();
        configs[season].heatStepsPerIteration = 3;
        configs[season].surfaceHeatExchangeWm2K = 18.0f;
    }
    for (int trial = 0; trial < 3; ++trial)
    {
        for (auto& cfg : configs) cfg.variableHeatStorage = trial == 2;
        for (auto& cfg : configs) cfg.reuseSeasonalCirculation = true;
        const auto cached = climateocean::solvePeriodicOcean(32, 16, seasons, configs, 3, 0.0);
        for (auto& cfg : configs) cfg.reuseSeasonalCirculation = false;
        const auto uncached = climateocean::solvePeriodicOcean(32, 16, seasons, configs, 3, 0.0);
        bool equal = cached.accepted && uncached.accepted && cached.years == 3 &&
            cached.annualEnthalpyDriftK == uncached.annualEnthalpyDriftK && cached.heatSubsteps == uncached.heatSubsteps;
        for (int season = 0; season < 4; ++season)
        {
            const auto& a = cached.seasons[season]; const auto& b = uncached.seasons[season];
            equal = equal && a.sstC == b.sstC && a.meanSstC == b.meanSstC &&
                a.iceThicknessMetres == b.iceThicknessMetres && a.meanIceCover == b.meanIceCover &&
                a.surfaceSkinTemperatureC == b.surfaceSkinTemperatureC &&
                a.meanSurfaceSkinTemperatureC == b.meanSurfaceSkinTemperatureC &&
                a.streamfunctionM3S == b.streamfunctionM3S &&
                a.eastVolumeTransportM3S == b.eastVolumeTransportM3S &&
                a.southVolumeTransportM3S == b.southVolumeTransportM3S && a.heatBudget.size() == b.heatBudget.size();
            equal = equal && a.storageDepthMetres == b.storageDepthMetres && a.reservoirTemperatureC == b.reservoirTemperatureC;
            for (std::size_t i = 0; i < a.heatBudget.size() && i < b.heatBudget.size(); ++i)
                equal = equal && a.heatBudget[i].horizontalWm2 == b.heatBudget[i].horizontalWm2 &&
                    a.heatBudget[i].surfaceWm2 == b.heatBudget[i].surfaceWm2 &&
                    a.heatBudget[i].verticalWm2 == b.heatBudget[i].verticalWm2 &&
                    a.heatBudget[i].storageWm2 == b.heatBudget[i].storageWm2;
        }
        expect(equal && cached.circulationSolves == 4 && uncached.circulationSolves == 12,
            "seasonal reuse must match recomputation bitwise and solve circulation only once per season");
        if (!equal || cached.circulationSolves != 4 || uncached.circulationSolves != 12)
        {
            std::cerr << std::setprecision(17) << "Reuse trial=" << trial
                << " solves=" << cached.circulationSolves << '/' << uncached.circulationSolves
                << " accepted=" << cached.accepted << '/' << uncached.accepted
                << " years=" << cached.years << '/' << uncached.years
                << " drift=" << cached.annualEnthalpyDriftK << '/' << uncached.annualEnthalpyDriftK
                << " steps=" << cached.heatSubsteps << '/' << uncached.heatSubsteps << '\n';
            for (int season=0;season<4;++season)
            {
                const auto& a=cached.seasons[season];const auto& b=uncached.seasons[season];
                double budget=0.0;std::size_t worst=0;
                for (std::size_t i=0;i<a.heatBudget.size();++i)
                    if (std::abs(a.heatBudget[i].horizontalWm2-b.heatBudget[i].horizontalWm2)>budget)
                    { budget=std::abs(a.heatBudget[i].horizontalWm2-b.heatBudget[i].horizontalWm2);worst=i; }
                std::cerr << " season=" << season << " sst_delta=" << difference(a.sstC,b.sstC)
                    << " mean_sst_delta=" << difference(a.meanSstC,b.meanSstC)
                    << " same_psi=" << (a.streamfunctionM3S==b.streamfunctionM3S)
                    << " horizontal_delta=" << budget << " cell=" << worst
                    << " east_delta=" << difference(a.eastCurrentMps,b.eastCurrentMps)
                    << " south_delta=" << difference(a.southCurrentMps,b.southCurrentMps)
                    << " ice_delta=" << difference(a.iceThicknessMetres,b.iceThicknessMetres)
                    << " budget=" << a.heatBudget[worst].horizontalWm2 << '/' << b.heatBudget[worst].horizontalWm2
                    << " vertical=" << a.heatBudget[worst].verticalWm2 << '/' << b.heatBudget[worst].verticalWm2 << '\n';
            }
        }
        for (auto& forcing : seasons) for (auto& wind : forcing.eastWindMps) wind *= -1.5f;
        for (auto& cfg : configs) cfg.linearBottomDragMps = 0.002f;
    }
    constexpr int n = 64, rows = n / 2;
    auto c = config();
    auto f = fixture(n, true);
    const auto grid = climategrid::makeSphericalGrid(n, rows, c.planetRadiusMetres);
    const double area = grid.cellAreasSquareMetres[0], face = grid.zonalFaceLengthsMetres[0];
    for (int y : {0, rows - 1})
        for (int x = 0; x < n; ++x) f.initialSstC[y * n + x] = x % 2 ? 10.0f : 0.0f;
    const auto implicit = climateocean::solveWindDrivenOcean(n, rows, f, c);
    const double ratio = c.oceanTimeStepSeconds * c.heatDiffusivityM2S * face * face / (area * area);
    double analyticError = 0.0;
    for (int y : {0, rows - 1})
        for (int x = 0; x < n; ++x)
            analyticError = std::max(analyticError, std::abs(implicit.sstC[y * n + x] -
                (5.0 + (x % 2 ? 5.0 : -5.0) / (1.0 + 4.0 * ratio))));
    expect(implicit.converged && analyticError < 1.0e-6 && budgetError(implicit) < 1.0e-7,
        "both polar rings must reproduce backward-Euler eigenmode damping and close every cell heat budget");
    expect(std::abs(implicit.heatBudgetResidualJ) / (2.0 * n * area * 2.5e8 * 10.0) < 1.0e-13,
        "closed zonal exchange must conserve globally integrated heat");

    // Time refinement against the retained explicit implementation at the SAME
    // spatial resolution. End states and interval means must both approach it.
    auto referenceConfig = c;
    referenceConfig.implicitZonalHeatDiffusion = false;
    referenceConfig.heatStepsPerIteration = 2048;
    referenceConfig.oceanTimeStepSeconds /= 2048;
    const auto reference = climateocean::solveWindDrivenOcean(n, rows, f, referenceConfig);
    auto coarseConfig = c;
    coarseConfig.heatStepsPerIteration = 8; coarseConfig.oceanTimeStepSeconds /= 8;
    const auto coarse = climateocean::solveWindDrivenOcean(n, rows, f, coarseConfig);
    coarseConfig.heatStepsPerIteration = 16; coarseConfig.oceanTimeStepSeconds /= 2;
    const auto fine = climateocean::solveWindDrivenOcean(n, rows, f, coarseConfig);
    expect(reference.converged && coarse.converged && fine.converged &&
        difference(fine.sstC, reference.sstC) < 0.6 * difference(coarse.sstC, reference.sstC) &&
        difference(fine.meanSstC, reference.meanSstC) < 0.6 * difference(coarse.meanSstC, reference.meanSstC),
        "implicit end states and seasonal means must approach a fine-step explicit reference under time refinement");

    auto basins = fixture(n, true);
    for (int y : {0, rows - 1})
        for (int x = 0; x < n; ++x)
        {
            const int i = y * n + x;
            basins.landMask[i] = x == 0 || x == n / 2;
            basins.initialSstC[i] = x < n / 2 ? 20.0f : 2.0f;
        }
    const auto isolated = climateocean::solveWindDrivenOcean(n, rows, basins, c);
    expect(isolated.converged && difference(isolated.sstC, basins.initialSstC) < 1.0e-6,
        "land barriers must keep disconnected constant-temperature basins unchanged, including at the longitude seam");

    // Only two wet cells straddle the seam. Their latent deficit is one K of
    // mixed-layer heat; low conductance partly melts ice, high conductance warms
    // both cells after melting. A temperature clamp cannot satisfy this test.
    auto ice = fixture(n, true);
    std::fill(ice.landMask.begin(), ice.landMask.end(), 1);
    ice.landMask[0] = ice.landMask[n - 1] = 0;
    std::fill(ice.initialSstC.begin(), ice.initialSstC.end(), c.freezingTemperatureC);
    ice.initialSstC[n - 1] = c.freezingTemperatureC + 2.0f;
    ice.initialIceThicknessMetres.assign(n * rows, 0.0f);
    const double capacity = c.waterDensityKgM3 * c.waterHeatCapacityJkgK * c.mixedLayerDepthMetres;
    const double latent = c.iceDensityKgM3 * c.latentHeatFusionJkg;
    ice.initialIceThicknessMetres[0] = static_cast<float>(capacity / latent);
    for (double r : {0.25, 4.0})
    {
        auto iceConfig = c;
        iceConfig.heatDiffusivityM2S = static_cast<float>(r * area * area / (c.oceanTimeStepSeconds * face * face));
        const auto result = climateocean::solveWindDrivenOcean(n, rows, ice, iceConfig);
        const double warm = r < 1.0 ? 2.0 / (1.0 + r) : (2.0 + r) / (1.0 + 2.0 * r);
        const double cold = r < 1.0 ? 0.0 : (r - 1.0) / (1.0 + 2.0 * r);
        const double thickness = r < 1.0 ? (1.0 - r * warm) * capacity / latent : 0.0;
        expect(result.converged && std::abs(result.sstC[n - 1] - c.freezingTemperatureC - warm) < 2.0e-6 &&
            std::abs(result.sstC[0] - c.freezingTemperatureC - cold) < 2.0e-6 &&
            std::abs(result.iceThicknessMetres[0] - thickness) < 1.0e-6 && budgetError(result) < 1.0e-6,
            "implicit seam exchange must reproduce analytic partial and complete melting with conserved latent heat");
        if (r > 1.0) expect(result.maximumZonalDiffusionIterations > 1,
            "melting must activate the newly liquid cell in the nonlinear solve");
    }

    auto zero = c; zero.heatDiffusivityM2S = 0.0f;
    const auto zeroImplicit = climateocean::solveWindDrivenOcean(n, rows, f, zero);
    zero.implicitZonalHeatDiffusion = false;
    const auto zeroExplicit = climateocean::solveWindDrivenOcean(n, rows, f, zero);
    expect(zeroImplicit.sstC == zeroExplicit.sstC && zeroImplicit.meanSstC == zeroExplicit.meanSstC,
        "zero diffusion must preserve the explicit integration results exactly");

    // Fixed-seed heterogeneous wet/land/ice rows stress repeated activation and
    // closed faces, rather than only the all-liquid cyclic matrix.
    std::mt19937 rng(20260907);
    for (int trial = 0; trial < 8; ++trial)
    {
        auto mixed = fixture(n, true);
        mixed.initialIceThicknessMetres.assign(n * rows, 0.0f);
        for (int y : {0, rows - 1})
            for (int x = 0; x < n; ++x)
            {
                const int i = y * n + x;
                mixed.landMask[i] = rng() % 13 == 0;
                const double h = (static_cast<int>(rng() % 1201) - 200) / 100.0;
                mixed.initialSstC[i] = static_cast<float>(c.freezingTemperatureC + std::max(0.0, h));
                mixed.initialIceThicknessMetres[i] = static_cast<float>(std::max(0.0, -h) * capacity / latent);
            }
        auto stiff = c;
        stiff.heatDiffusivityM2S *= static_cast<float>(std::pow(2.0, trial));
        const auto result = climateocean::solveWindDrivenOcean(n, rows, mixed, stiff);
        bool bounded = true;
        for (std::size_t i = 0; i < result.sstC.size(); ++i)
            if (!mixed.landMask[i]) bounded = bounded && result.sstC[i] >= c.freezingTemperatureC &&
                result.sstC[i] <= c.freezingTemperatureC + 10.00001f;
        if (!result.converged || !bounded || budgetError(result) >= 1.0e-5)
            std::cerr << "Mixed polar fixture seed=20260907 trial=" << trial << " residual_k="
                << result.maximumZonalDiffusionResidualK << " iterations=" << result.maximumZonalDiffusionIterations << '\n';
        expect(result.converged && bounded && budgetError(result) < 1.0e-5,
            "mixed polar basins must remain bounded and conserve local energy through freezing/melting");
    }

    // Large meridional diffusion still needs explicit substeps after removing
    // the zonal restriction. Constant fields must survive those substeps.
    auto meridional = c; meridional.heatDiffusivityM2S = 1.0e8f;
    const auto meridionalState = climateocean::solveWindDrivenOcean(n, rows, fixture(n, false), meridional);
    expect(meridionalState.converged && meridionalState.heatSubsteps > 1 &&
        meridionalState.explicitHeatSubsteps > meridionalState.heatSubsteps &&
        difference(meridionalState.sstC, fixture(n, false).initialSstC) < 1.0e-6,
        "implicit zonal diffusion must retain the meridional explicit stability restriction");

    auto high = fixture(512, true);
    const auto highState = climateocean::solveWindDrivenOcean(512, 256, high, c);
    const auto& limiter = meridionalState.stabilityLimit;
    expect(limiter.cell >= 0 && limiter.meridionalDiffusionRate > 0.0 &&
        limiter.zonalAdvectionRate == 0.0 && limiter.meridionalAdvectionRate == 0.0 &&
        std::abs(limiter.totalRate - limiter.meridionalDiffusionRate - limiter.verticalRate - limiter.surfaceRate) < 1.0e-14 &&
        meridionalState.actualHeatStepSeconds <= 0.7 / limiter.totalRate,
        "ocean limiter diagnostics must identify and reconstruct the retained explicit stability bound");
    expect(highState.converged && highState.heatSubsteps == 1 && highState.explicitHeatSubsteps >= 800 &&
        difference(highState.sstC, high.initialSstC) < 1.0e-6,
        "512x256 open polar rings must preserve constants in one implicit step despite the explicit polar restriction");
    auto threaded = fixture(256, false);
    threaded.initialIceThicknessMetres.assign(threaded.initialSstC.size(), 0.0f);
    for (std::size_t i = 0; i < threaded.initialSstC.size(); ++i)
    {
        threaded.landMask[i] = i % 19 == 0;
        threaded.initialSstC[i] = i % 7 ? static_cast<float>(i % 11) : c.freezingTemperatureC;
        threaded.initialIceThicknessMetres[i] = i % 7 ? 0.0f : 0.05f;
    }
    auto serialConfig = c; serialConfig.parallelZonalHeatDiffusion = false;
    const auto serial = climateocean::solveWindDrivenOcean(256, 128, threaded, serialConfig);
    const auto parallel = climateocean::solveWindDrivenOcean(256, 128, threaded, c);
    bool same = serial.converged && parallel.converged && serial.sstC == parallel.sstC &&
        serial.meanSstC == parallel.meanSstC && serial.iceThicknessMetres == parallel.iceThicknessMetres &&
        serial.meanIceCover == parallel.meanIceCover && serial.meanSurfaceSkinTemperatureC == parallel.meanSurfaceSkinTemperatureC &&
        serial.maximumZonalDiffusionResidualK == parallel.maximumZonalDiffusionResidualK &&
        serial.maximumZonalDiffusionIterations == parallel.maximumZonalDiffusionIterations;
    for (std::size_t i = 0; i < serial.heatBudget.size(); ++i)
        same = same && serial.heatBudget[i].horizontalWm2 == parallel.heatBudget[i].horizontalWm2 &&
            serial.heatBudget[i].storageWm2 == parallel.heatBudget[i].storageWm2;
    expect(same, "parallel zonal rows preserve serial enthalpy, means, budgets and diagnostics bitwise through coastal melting");
}

int accuracy(int columns)
{
    if (columns < 16 || columns > 512 || columns % 4 != 0) return 2;
    constexpr int days = 4, referenceStepsPerDay = 512;
    std::cout << std::setprecision(12) << "{\"columns\":" << columns << ",\"rows\":" << columns / 2
        << ",\"days\":" << days << ",\"reference\":\"explicit 512 minimum steps/day\",\"cases\":[";
    for (int trial = 0; trial < 2; ++trial)
    {
        auto f = fixture(columns, true);
        auto c = config();
        f.initialIceThicknessMetres.assign(f.initialSstC.size(), 0.0f);
        for (int y : {0, columns / 2 - 1})
            for (int x = 0; x < columns; ++x)
            {
                const int i = y * columns + x;
                f.initialSstC[i] = x < columns / 2 ? 12.0f : (trial ? c.freezingTemperatureC : 0.0f);
                if (trial)
                {
                    f.landMask[i] = x == columns / 4 || x == 3 * columns / 4;
                    f.initialIceThicknessMetres[i] = x < columns / 2 ? 0.0f : 0.25f;
                    f.atmosphericTemperatureC[i] = -5.0f;
                }
            }
        c.surfaceHeatExchangeWm2K = trial ? 18.0f : 0.0f;
        c.implicitZonalHeatDiffusion = false;
        c.oceanTimeStepSeconds = 86400.0f / referenceStepsPerDay;
        c.heatStepsPerIteration = days * referenceStepsPerDay;
        const auto reference = climateocean::solveWindDrivenOcean(columns, columns / 2, f, c);
        if (!reference.converged || budgetError(reference) > 1.0e-5) return 1;
        if (trial) std::cout << ',';
        std::cout << "{\"fixture\":\"" << (trial ? "coastal-polar-ice-v1" : "sharp-polar-front-v1")
            << "\",\"reference_steps\":" << reference.heatSubsteps << ",\"steps\":[";
        bool first = true;
        for (int hours : {24, 12, 6, 3})
        {
            c.implicitZonalHeatDiffusion = true;
            c.oceanTimeStepSeconds = static_cast<float>(hours * 3600);
            c.heatStepsPerIteration = days * (24 / hours);
            const auto start = std::chrono::steady_clock::now();
            const auto state = climateocean::solveWindDrivenOcean(columns, columns / 2, f, c);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (!state.converged || budgetError(state) > 1.0e-5) return 1;
            double fluxError = 0.0;
            for (std::size_t i = 0; i < state.heatBudget.size(); ++i)
                fluxError = std::max(fluxError, std::abs(state.heatBudget[i].horizontalWm2 - reference.heatBudget[i].horizontalWm2));
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"hours\":" << hours << ",\"actual_steps\":" << state.heatSubsteps << ",\"seconds\":" << seconds
                << ",\"max_sst_error_k\":" << difference(state.sstC, reference.sstC)
                << ",\"max_mean_sst_error_k\":" << difference(state.meanSstC, reference.meanSstC)
                << ",\"max_ice_error_m\":" << difference(state.iceThicknessMetres, reference.iceThicknessMetres)
                << ",\"max_mean_ice_cover_error\":" << difference(state.meanIceCover, reference.meanIceCover)
                << ",\"max_horizontal_flux_error_wm2\":" << fluxError
                << ",\"cell_budget_error_wm2\":" << budgetError(state) << '}';
        }
        std::cout << "]}";
    }
    std::cout << "]}\n";
    return 0;
}

int heatKernels(int columns, int days)
{
    if (columns < 8 || columns > 512 || columns % 4 || days < 1 || days > 30) return 2;
    const auto f = heatFixture(columns);
    auto c = config(); c.heatStepsPerIteration = days;
    c.surfaceHeatExchangeWm2K = 18.0f; c.deepWaterTemperatureContrastK = 4.0f;
    c.gatherHeatFluxes = false;
    const auto reference = climateocean::solveWindDrivenOcean(columns, columns / 2, f, c);
    if (!reference.converged) return 1;
    constexpr std::array<int, 6> workers{0, 1, 2, 4, 8, 16};
    std::array<std::vector<double>, 6> times, fluxTimes, updateTimes, zonalTimes;
    for (int repeat = 0; repeat < 3; ++repeat)
        for (std::size_t index = 0; index < workers.size(); ++index)
        {
            const auto which = repeat % 2 ? workers.size() - 1 - index : index;
            c.gatherHeatFluxes = which != 0; c.heatWorkers = workers[which];
            const auto state = climateocean::solveWindDrivenOcean(columns, columns / 2, f, c);
            if (!state.converged || !sameHeatState(reference, state) || budgetError(state) > 1.0e-5)
            {
                std::cerr << "Heat kernel mismatch width=" << columns << " workers=" << workers[which]
                    << " repeat=" << repeat << " sst_delta=" << difference(reference.sstC, state.sstC) << '\n';
                return 1;
            }
            times[which].push_back(state.heatSeconds); fluxTimes[which].push_back(state.heatFluxSeconds);
            updateTimes[which].push_back(state.heatUpdateSeconds); zonalTimes[which].push_back(state.zonalHeatSeconds);
        }
    const auto median = [](auto values) { std::sort(values.begin(), values.end()); return values[1]; };
    std::cout << std::setprecision(12) << "{\"fixture\":\"coastal-variable-depth-ice-v1\",\"columns\":" << columns
        << ",\"days\":" << days << ",\"repeats\":3,\"steps\":" << reference.heatSubsteps << ",\"cases\":[";
    for (std::size_t i = 0; i < workers.size(); ++i)
    {
        if (i) std::cout << ',';
        std::cout << "{\"kernel\":\"" << (i ? "gather" : "legacy") << "\",\"workers\":" << workers[i]
            << ",\"heat_seconds_median\":" << median(times[i]) << ",\"flux_seconds_median\":" << median(fluxTimes[i])
            << ",\"update_seconds_median\":" << median(updateTimes[i]) << ",\"zonal_seconds_median\":" << median(zonalTimes[i])
            << ",\"heat_seconds_min\":" << *std::min_element(times[i].begin(), times[i].end())
            << ",\"heat_seconds_max\":" << *std::max_element(times[i].begin(), times[i].end()) << ",\"exact_fields\":true}";
    }
    std::cout << "]}\n";
    return 0;
}

int benchmark(int columns, int days)
{
    if (columns < 8 || columns > 512 || columns % 4 != 0 || days < 1 || days > 90) return 2;
    auto f = fixture(columns, false);
    auto c = config(); c.heatStepsPerIteration = days;
    const auto grid = climategrid::makeSphericalGrid(columns, columns / 2, c.planetRadiusMetres);
    constexpr double pi = 3.14159265358979323846;
    for (int y = 0; y < columns / 2; ++y)
        for (int x = 0; x < columns; ++x)
            f.initialSstC[y * columns + x] = static_cast<float>(10.0 +
                5.0 * std::cos(grid.latitudeCentresRadians[y]) * std::cos(2.0 * pi * (x + 0.5) / columns));
    std::vector<double> explicitTimes, implicitTimes;
    climateocean::OceanState explicitState, implicitState;
    for (int repeat = 0; repeat < 3; ++repeat)
        for (int method = 0; method < 2; ++method)
        {
            const bool implicit = (method + repeat) % 2 != 0;
            c.implicitZonalHeatDiffusion = implicit;
            const auto start = std::chrono::steady_clock::now();
            auto state = climateocean::solveWindDrivenOcean(columns, columns / 2, f, c);
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (!state.converged || budgetError(state) > 1.0e-5) return 1;
            (implicit ? implicitTimes : explicitTimes).push_back(seconds);
            (implicit ? implicitState : explicitState) = std::move(state);
        }
    std::sort(explicitTimes.begin(), explicitTimes.end()); std::sort(implicitTimes.begin(), implicitTimes.end());
    std::cout << std::setprecision(12) << "{\"fixture\":\"smooth-warm-aquaplanet-zero-wind-v1\",\"columns\":" << columns
        << ",\"rows\":" << columns / 2 << ",\"days\":" << days << ",\"repeats\":3"
        << ",\"explicit_seconds_median\":" << explicitTimes[1] << ",\"implicit_seconds_median\":" << implicitTimes[1]
        << ",\"speedup\":" << explicitTimes[1] / implicitTimes[1]
        << ",\"explicit_substeps\":" << explicitState.heatSubsteps << ",\"implicit_substeps\":" << implicitState.heatSubsteps
        << ",\"max_sst_difference_k\":" << difference(explicitState.sstC, implicitState.sstC)
        << ",\"max_mean_sst_difference_k\":" << difference(explicitState.meanSstC, implicitState.meanSstC)
        << ",\"implicit_cell_budget_error_wm2\":" << budgetError(implicitState)
        << ",\"explicit_cell_budget_error_wm2\":" << budgetError(explicitState)
        << ",\"implicit_nonlinear_residual_k\":" << implicitState.maximumZonalDiffusionResidualK
        << ",\"implicit_max_iterations\":" << implicitState.maximumZonalDiffusionIterations << "}\n";
    return 0;
}
}

int main(int argc, char** argv)
{
    if (argc == 4 && std::string(argv[1]) == "--heat-kernels")
    {
        try { return heatKernels(std::stoi(argv[2]), std::stoi(argv[3])); }
        catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 3 && std::string(argv[1]) == "--accuracy")
    {
        try { return accuracy(std::stoi(argv[2])); }
        catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc == 4 && std::string(argv[1]) == "--benchmark")
    {
        try { return benchmark(std::stoi(argv[2]), std::stoi(argv[3])); }
        catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
    }
    if (argc != 1) { std::cerr << "Usage: climate_ocean_heat_tests [--heat-kernels WIDTH DAYS | --benchmark WIDTH DAYS | --accuracy WIDTH]\n"; return 2; }
    runTests();
    return failures ? 1 : 0;
}
