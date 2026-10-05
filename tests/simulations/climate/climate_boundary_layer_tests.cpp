#include "simulations/climate/climate_atmosphere.hpp"
#include "core/climate_grid.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <gtest/gtest.h>

namespace
{
void expect(bool condition, const char* message)
{
    EXPECT_TRUE(condition) << message;
}
}

namespace
{

// The 32 x 16 grid of the boundary-layer tests: closure enabled, uniform 27 C
// surface, zero heating and mechanical forcing.
class ClimateBoundaryLayer : public ::testing::Test
{
protected:
    static constexpr float pi = 3.14159265358979323846f;
    static constexpr int columns = 32, rows = 16, count = columns * rows;
    climateatmosphere::ModeSeparatedCirculationConfig config;
    std::vector<float> temperature = std::vector<float>(count, 27), heating = std::vector<float>(count, 0),
        mechanical = std::vector<float>(count, 0), zeroRow = std::vector<float>(rows, 0);

    ClimateBoundaryLayer()
    {
        config.boundaryLayer.enabled = true;
        config.relativeTolerance = 1.0e-5f;
        config.maximumIterations = 2000;
    }

    // The regional warm patch of the old body, added to the uniform temperature.
    void warmPatch()
    {
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
                temperature[y * columns + x] += 3.0f * std::cos(2 * pi * x / columns) *
                    std::exp(-std::pow((y + 0.5f - rows / 2) / 2, 2));
    }

    // The circulation with the given config's zonal mode only, before any closure is applied.
    auto circulation(const climateatmosphere::ModeSeparatedCirculationConfig& with) const
    {
        return climateatmosphere::solveModeSeparatedCirculation(columns, rows, zeroRow, heating, mechanical, with);
    }
};

}

TEST_F(ClimateBoundaryLayer, MixedDragClosure)
{
    using namespace climateatmosphere;
    for (float f : {-4.0e-5f, 0.0f, 4.0e-5f})
        for (float q : {0.0f, 0.0013f / 300.0f})
        {
            const float east = 0.00017f, north = -0.00008f, linear = 3.5e-5f;
            const auto wind = steadyMixedDragCoriolisWind(east, north, f, linear, q);
            const float drag = linear + q * std::hypot(wind.eastMetresPerSecond, wind.southMetresPerSecond);
            expect(std::abs(drag * wind.eastMetresPerSecond + f * wind.southMetresPerSecond - east) < 1.0e-10f &&
                std::abs(-drag * wind.southMetresPerSecond + f * wind.eastMetresPerSecond - north) < 1.0e-10f,
                "mixed drag must close both momentum equations in either hemisphere and at the equator");
        }
}

TEST_F(ClimateBoundaryLayer, InversionWindsOffsetConvergence)
{
    using namespace climateatmosphere;
    const float f = 4.0e-5f, upperEast = -6.0f, entrainment = 2.0e-5f, friction = 1.5e-5f;
    const float pressureForceNorth = f * upperEast / 1.225f;
    const auto withoutMixing = steadyMixedDragCoriolisWind(0, pressureForceNorth, f, friction, 0);
    const auto withMixing = steadyMixedDragCoriolisWind(entrainment * upperEast, pressureForceNorth, f, friction + entrainment, 0);
    expect(std::abs(withMixing.southMetresPerSecond) < 0.6f * std::abs(withoutMixing.southMetresPerSecond) &&
        withMixing.eastMetresPerSecond < 0,
        "geostrophic inversion winds offset broad pressure-driven convergence while retaining easterlies");
}

TEST_F(ClimateBoundaryLayer, UniformLayer)
{
    using namespace climateatmosphere;
    auto initial = circulation(config);
    const auto uniform = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    expect(uniform.inversionSolver.converged && applyBoundaryLayer(uniform, config, initial),
        "uniform boundary layer must solve and apply");
    expect(initial.areaWeightedKineticEnergyJm2 == 0, "uniform temperature and zero forcing must not invent wind");
}

TEST_F(ClimateBoundaryLayer, WarmPatchFlow)
{
    using namespace climateatmosphere;
    warmPatch();
    const auto warm = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    auto flow = circulation(config);
    expect(warm.inversionSolver.converged && applyBoundaryLayer(warm, config, flow),
        "a warm regional patch must produce a diagnosed boundary-layer flow");
    expect(flow.surfaceEastWindMps[(rows/2) * columns + 2] < 0 &&
        flow.surfaceEastWindMps[(rows/2) * columns + columns - 2] > 0,
        "regional warm-patch convergence must follow longitude-dependent temperature gradients");
    expect(flow.maximumBoundaryLayerMomentumResidualMps2 < 1.0e-9 &&
        std::abs(flow.areaWeightedMassAnomalyKgM2) < 1.0e-3,
        "applied closure must retain momentum balance and a zero global pressure gauge");
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    double ascent = 0, area = 0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        { ascent += grid.cellAreasSquareMetres[y] * flow.ascentHpaPerDay[y * columns + x]; area += grid.cellAreasSquareMetres[y]; }
    expect(std::abs(ascent / area) < 1.0e-6, "new winds must retain globally closed ascent fluxes");
    {
        const int y = rows/2, x = 0;
        const auto i = grid.index(x, y);
        const auto& u = flow.surfaceEastWindMps;
        const auto& v = flow.surfaceSouthWindMps;
        const double flux = 0.5 * (u[grid.index(x+1,y)] - u[grid.index(x-1,y)]) * grid.zonalFaceLengthsMetres[y] +
            0.5 * (v[i] + v[grid.index(x,y+1)]) * grid.southFaceLengthsMetres[y] -
            0.5 * (v[i] + v[grid.index(x,y-1)]) * grid.northFaceLengthsMetres[y];
        expect(std::abs(flow.ascentHpaPerDay[i] + 86400.0 * 150.0 * flux / grid.cellAreasSquareMetres[y]) < 1.0e-4,
            "inversion ascent must use the actual 150 hPa layer mass, not its wave equivalent depth");
    }
}

TEST_F(ClimateBoundaryLayer, PressureGaugeInvariance)
{
    using namespace climateatmosphere;
    warmPatch();
    const auto warm = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    auto flow = circulation(config);
    applyBoundaryLayer(warm, config, flow);
    // Independently chosen pressure datums cannot create a transition-belt jet.
    auto shifted = warm;
    for (auto& p : shifted.surfacePressureHpa) p += 17.0f;
    auto gaugeFlow = circulation(config);
    applyBoundaryLayer(shifted, config, gaugeFlow);
    double gaugeError = 0;
    for (int i = 0; i < count; ++i) gaugeError = std::max(gaugeError,
        static_cast<double>(std::abs(gaugeFlow.surfaceEastWindMps[i] - flow.surfaceEastWindMps[i])));
    expect(gaugeError < 1.0e-4, "changing the inversion pressure gauge must not change the surface wind");
}

TEST_F(ClimateBoundaryLayer, NonfiniteStateRejected)
{
    using namespace climateatmosphere;
    warmPatch();
    const auto warm = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    auto flow = circulation(config);
    applyBoundaryLayer(warm, config, flow);
    auto invalidState = warm;
    invalidState.inversionEastWindMps[0] = std::numeric_limits<float>::quiet_NaN();
    const auto validPressure = flow.surfacePressureAnomalyHpa;
    expect(!applyBoundaryLayer(invalidState, config, flow) && flow.surfacePressureAnomalyHpa == validPressure,
        "nonfinite accepted states must be rejected before changing circulation");
}

TEST_F(ClimateBoundaryLayer, DisabledClosureKeepsLegacy)
{
    using namespace climateatmosphere;
    warmPatch();
    const auto warm = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    auto legacyConfig = config;
    legacyConfig.boundaryLayer.enabled = false;
    const auto legacy = circulation(legacyConfig);
    auto ignoredState = legacy;
    expect(!applyBoundaryLayer(warm, legacyConfig, ignoredState) && ignoredState.surfacePressureAnomalyHpa == legacy.surfacePressureAnomalyHpa,
        "an explicitly disabled closure must not overwrite the legacy result");
    legacyConfig.inversionEastWindMps.clear(); legacyConfig.inversionSouthWindMps.clear();
    const auto withoutState = circulation(legacyConfig);
    expect(legacy.surfaceEastWindMps == withoutState.surfaceEastWindMps && legacy.ascentHpaPerDay == withoutState.ascentHpaPerDay,
        "stored experimental state must not affect legacy wind or ascent");
}

TEST_F(ClimateBoundaryLayer, HadleyProfileMixing)
{
    using namespace climateatmosphere;
    auto jetConfig = config;
    jetConfig.inversionEquilibriumPressureHpa.resize(count);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            jetConfig.inversionEquilibriumPressureHpa[y * columns + x] =
                -6.0f * std::exp(-std::pow((90.0f - 180.0f * (y + .5f) / rows) / 30.0f, 2));
    const auto jet = solveBoundaryLayer(columns, rows, std::vector<float>(count, 27), heating, mechanical, jetConfig);
    auto mixedJet = circulation(jetConfig);
    applyBoundaryLayer(jet, jetConfig, mixedJet);
    jetConfig.boundaryLayer.entrainmentRatePerSecond = 0;
    auto unmixedJet = circulation(jetConfig);
    applyBoundaryLayer(jet, jetConfig, unmixedJet);
    double mixedInflow = 0, unmixedInflow = 0;
    for (int y : {rows/2-1, rows/2})
        for (int x = 0; x < columns; ++x)
        {
            mixedInflow += std::pow(mixedJet.surfaceSouthWindMps[y * columns + x], 2);
            unmixedInflow += std::pow(unmixedJet.surfaceSouthWindMps[y * columns + x], 2);
        }
    expect(jet.inversionSolver.converged && mixedInflow < 0.8 * unmixedInflow,
        "a mass-adjusted Hadley profile aloft must not impose its full inward surface pull");
}

TEST_F(ClimateBoundaryLayer, HeatingRetainsZonalMode)
{
    using namespace climateatmosphere;
    warmPatch();
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x) heating[y * columns + x] = 70.0f * std::cos(pi * (y + .5f) / rows);
    const auto forced = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    expect(forced.inversionSolver.converged && *std::max_element(forced.inversionPressureHpa.begin(), forced.inversionPressureHpa.end()) > 0.01f,
        "lower-tropospheric heating must retain and solve its zonal mode rather than project it away");
}

TEST_F(ClimateBoundaryLayer, UnacceptedSolveRejected)
{
    using namespace climateatmosphere;
    warmPatch();
    const auto warm = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    auto flow = circulation(config);
    applyBoundaryLayer(warm, config, flow);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x) heating[y * columns + x] = 70.0f * std::cos(pi * (y + .5f) / rows);
    config.maximumIterations = 1;
    config.stationarySolver = StationarySolver::LegacyJacobi; // Force non-convergence for the rejection test.
    config.relativeTolerance = 1.0e-12f;
    const auto rejected = solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config);
    const auto before = flow.surfacePressureAnomalyHpa;
    expect(!rejected.inversionSolver.converged && !applyBoundaryLayer(rejected, config, flow) && before == flow.surfacePressureAnomalyHpa,
        "an unaccepted inversion solve must not partially overwrite circulation");
    temperature[0] = std::numeric_limits<float>::quiet_NaN();
    expect(!solveBoundaryLayer(columns, rows, temperature, heating, mechanical, config).inversionSolver.converged,
        "nonfinite boundary inputs must be rejected");
}
