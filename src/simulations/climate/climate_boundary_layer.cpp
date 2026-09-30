#include "climate_atmosphere.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>

namespace climateatmosphere
{
// Surface momentum/thickness balance: Back & Bretherton (2009), sections 3-4,
// https://aos.wisc.edu/~lback/SSTgradients.pdf. The separate forced 850 hPa
// diagnosis here replaces that paper's prescribed reanalysis boundary data.
namespace
{
constexpr double pi = 3.14159265358979323846;
double mean(const climategrid::SphericalGrid& grid, const std::vector<float>& field)
{
    double sum = 0.0, area = 0.0;
    for (int y = 0; y < grid.rows; ++y)
        for (int x = 0; x < grid.columns; ++x)
        {
            sum += grid.cellAreasSquareMetres[y] * (field[grid.index(x, y)] - field.front());
            area += grid.cellAreasSquareMetres[y];
        }
    return field.front() + sum / area;
}
}

float boundaryLayerWeight(float latitude, const BoundaryLayerConfig& config)
{
    if (!config.enabled) return 0.0f;
    const float s = std::clamp((std::abs(latitude) - config.fullLatitudeDegrees) /
        std::max(0.001f, config.outerLatitudeDegrees - config.fullLatitudeDegrees), 0.0f, 1.0f);
    return 1.0f - s * s * (3.0f - 2.0f * s);
}

HorizontalWind steadyMixedDragCoriolisWind(float east, float north,
    float coriolis, float linear, float quadratic)
{
    if (!std::isfinite(east) || !std::isfinite(north) || !std::isfinite(coriolis) || !std::isfinite(linear) ||
        !std::isfinite(quadratic) || linear < 0.0f || quadratic < 0.0f || (linear == 0.0f && quadratic == 0.0f)) return {};
    const double force = std::hypot(east, north);
    if (force == 0.0) return {};
    double lo = 0.0, hi = quadratic > 0.0f ? std::sqrt(force / quadratic) : force / linear;
    // |F| = speed * hypot(linear + quadratic * speed, f), a monotone scalar equation.
    for (int i = 0; i < 40; ++i)
    {
        const double speed = 0.5 * (lo + hi);
        if (speed * std::hypot(linear + quadratic * speed, coriolis) > force) hi = speed;
        else lo = speed;
    }
    const double drag = linear + quadratic * 0.5 * (lo + hi);
    const double denominator = drag * drag + coriolis * coriolis;
    return {static_cast<float>((drag * east + coriolis * north) / denominator),
        static_cast<float>((coriolis * east - drag * north) / denominator)};
}

BoundaryLayerState solveBoundaryLayer(int columns, int rows,
    const std::vector<float>& temperature, const std::vector<float>& heating,
    const std::vector<float>& terrainPressure, const ModeSeparatedCirculationConfig& config)
{
    BoundaryLayerState result;
    const auto& layer = config.boundaryLayer;
    const std::size_t count = static_cast<std::size_t>(std::max(0, columns)) * std::max(0, rows);
    const auto valid = [&](const auto& field) { return field.size() == count &&
        std::all_of(field.begin(), field.end(), [](float v) { return std::isfinite(v); }); };
    for (float value : {layer.inversionPressurePa, layer.inversionDensityKgM3, layer.inversionDragTimeSeconds,
        layer.entrainmentRatePerSecond, layer.oceanFrictionRatePerSecond, layer.fullLatitudeDegrees,
        layer.outerLatitudeDegrees, config.gravityMetresPerSecondSquared, config.planetRadiusMetres,
        config.airDensityKgM3, config.surfaceDampingTimeSeconds, config.rotationRatePerSecond, config.rotationDirection})
        if (!std::isfinite(value)) return result;
    if (!climategrid::validGlobalGridDimensions(columns, rows) || !valid(temperature) || !valid(heating) ||
        (!config.inversionEquilibriumPressureHpa.empty() && !valid(config.inversionEquilibriumPressureHpa)) ||
        !valid(terrainPressure) || layer.inversionPressurePa <= 50000.0f || layer.inversionPressurePa >= 100000.0f ||
        layer.inversionDensityKgM3 <= 0.0f || layer.inversionDragTimeSeconds <= 0.0f ||
        config.gravityMetresPerSecondSquared <= 0.0f || config.planetRadiusMetres <= 0.0f ||
        config.airDensityKgM3 <= 0.0f || config.surfaceDampingTimeSeconds <= 0.0f ||
        layer.entrainmentRatePerSecond < 0.0f || layer.oceanFrictionRatePerSecond <= 0.0f ||
        layer.fullLatitudeDegrees < 0.0f || layer.outerLatitudeDegrees <= layer.fullLatitudeDegrees ||
        layer.outerLatitudeDegrees > 90.0f) return result;
    result.columns = columns; result.rows = rows;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, config.planetRadiusMetres);
    const double temperatureMean = mean(grid, temperature), heatingMean = mean(grid, heating);
    const double columnMass = 100000.0 / config.gravityMetresPerSecondSquared;
    const double thermalTime = config.surfaceDampingTimeSeconds;
    const double thickness = 287.05 * std::log(layer.inversionPressurePa / 50000.0);
    const double inversionHeight = 287.05 * (temperatureMean + 273.15) / config.gravityMetresPerSecondSquared *
        std::log(100000.0 / layer.inversionPressurePa);
    std::vector<float> equilibrium(count), drag(count, layer.inversionDragTimeSeconds);
    result.inversionTemperatureC.resize(count);
    result.thermalPressureHpa.resize(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        // Reduced lower-free-tropospheric mode: net column heating produces a
        // damped temperature anomaly and a lower-level mass deficit. The broad
        // overturning profile is equilibrium forcing aloft; its adjusted zonal
        // mean and regional modes share the same pressure/divergence solve.
        // This approximates the 850 hPa boundary data used by BB09; it is not
        // the existing 500 hPa geopotential or an observed input.
        const double anomaly = thermalTime * (heating[i] - heatingMean) / (1004.0 * columnMass);
        result.inversionTemperatureC[i] = static_cast<float>(temperatureMean - 0.0065 * inversionHeight + anomaly);
        equilibrium[i] = static_cast<float>(-layer.inversionDensityKgM3 * thickness * anomaly / 100.0 + terrainPressure[i]);
        if (config.inversionEquilibriumPressureHpa.size() == count) equilibrium[i] += config.inversionEquilibriumPressureHpa[i];
        // Hydrostatic thickness with anomalies varying linearly from the
        // surface to the diagnosed inversion (BB09, section 4).
        result.thermalPressureHpa[i] = thermalModePressureAnomalyHpa(
            static_cast<float>(0.5 * (temperature[i] - temperatureMean + anomaly)),
            config.airDensityKgM3, 100000.0f, layer.inversionPressurePa);
    }
    result.inversionSolver = solveSteadyStationaryWavePressure(columns, rows, equilibrium, drag,
        config.surfaceEquivalentPressureDepthHpa, config.surfaceDampingTimeSeconds,
        layer.inversionDensityKgM3, config.planetRadiusMetres, config.rotationRatePerSecond,
        config.rotationDirection, false, config.maximumIterations, config.relativeTolerance, config.solverRestartLength,
        config.stationarySolver, config.stationaryCaptureDirectory);
    if (!result.inversionSolver.converged) return result;
    result.inversionPressureHpa = result.inversionSolver.pressureAnomalyHpa;
    result.surfacePressureHpa.resize(count);
    result.inversionEastWindMps.resize(count); result.inversionSouthWindMps.resize(count);
    const double dy = config.planetRadiusMetres * grid.latitudeSpacingRadians;
    for (int y = 0; y < rows; ++y)
    {
        const int north = std::max(0, y - 1), south = std::min(rows - 1, y + 1);
        const double dx = config.planetRadiusMetres * grid.longitudeSpacingRadians * std::cos(grid.latitudeCentresRadians[y]);
        for (int x = 0; x < columns; ++x)
        {
            const auto i = grid.index(x, y);
            const float east = static_cast<float>(-100.0 * (result.inversionPressureHpa[grid.index(x + 1, y)] -
                result.inversionPressureHpa[grid.index(x - 1, y)]) / (2.0 * dx * layer.inversionDensityKgM3));
            const float northForce = static_cast<float>(-100.0 * (result.inversionPressureHpa[grid.index(x, north)] -
                result.inversionPressureHpa[grid.index(x, south)]) / ((south - north) * dy * layer.inversionDensityKgM3));
            const auto wind = steadyRayleighCoriolisWind(east, northForce,
                static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / pi), layer.inversionDragTimeSeconds,
                config.rotationRatePerSecond, config.rotationDirection);
            result.inversionEastWindMps[i] = wind.eastMetresPerSecond;
            result.inversionSouthWindMps[i] = wind.southMetresPerSecond;
            result.surfacePressureHpa[i] = result.inversionPressureHpa[i] + result.thermalPressureHpa[i];
        }
    }
    return result;
}

bool applyBoundaryLayer(const BoundaryLayerState& state,
    ModeSeparatedCirculationConfig& config, ModeSeparatedCirculation& circulation)
{
    const std::size_t count = circulation.surfacePressureAnomalyHpa.size();
    const auto valid = [&](const auto& field) { return field.size() == count &&
        std::all_of(field.begin(), field.end(), [](float v) { return std::isfinite(v); }); };
    if (!config.boundaryLayer.enabled || !state.inversionSolver.converged ||
        !climategrid::validGlobalGridDimensions(state.columns, state.rows) ||
        count != static_cast<std::size_t>(state.columns) * state.rows ||
        !valid(state.surfacePressureHpa) || !valid(state.inversionEastWindMps) ||
        !valid(state.inversionSouthWindMps) || !valid(circulation.upperHeightAnomalyMetres) || count == 0) return false;
    const auto grid = climategrid::makeSphericalGrid(state.columns, state.rows, config.planetRadiusMetres);
    // Align pressure gauges in the transition belt before blending. Constant
    // pressure offsets must never generate spurious subtropical jets.
    double offset = 0.0, area = 0.0;
    for (int y = 0; y < state.rows; ++y)
    {
        const float w = boundaryLayerWeight(static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / pi), config.boundaryLayer);
        for (int x = 0; x < state.columns; ++x)
        {
            const auto i = grid.index(x, y);
            const double weight = grid.cellAreasSquareMetres[y] * w * (1.0 - w);
            offset += weight * (circulation.surfacePressureAnomalyHpa[i] - state.surfacePressureHpa[i]);
            area += weight;
        }
    }
    offset = area > 0.0 ? offset / area : mean(grid, circulation.surfacePressureAnomalyHpa) - mean(grid, state.surfacePressureHpa);
    for (int y = 0; y < state.rows; ++y)
    {
        const float w = boundaryLayerWeight(static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / pi), config.boundaryLayer);
        for (int x = 0; x < state.columns; ++x)
        {
            const auto i = grid.index(x, y);
            circulation.surfacePressureAnomalyHpa[i] += static_cast<float>(w *
                (state.surfacePressureHpa[i] + offset - circulation.surfacePressureAnomalyHpa[i]));
        }
    }
    const double gauge = mean(grid, circulation.surfacePressureAnomalyHpa);
    for (auto& p : circulation.surfacePressureAnomalyHpa) p -= static_cast<float>(gauge);
    config.inversionEastWindMps = state.inversionEastWindMps;
    config.inversionSouthWindMps = state.inversionSouthWindMps;
    diagnoseModeWinds(state.columns, state.rows, config, circulation);
    return true;
}
}
