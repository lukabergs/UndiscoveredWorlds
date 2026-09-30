#pragma once
#include "climate_atmosphere.hpp"
#include "climate_grid.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
struct SurfaceMomentumAdvection
{
    std::vector<float> east, south;
    int iterations = 0;
    double maximumChangeMps = 0.0;
    double relativeMomentumResidual = 0.0;
    double forcingWorkWm2 = 0.0, dragWorkWm2 = 0.0, advectionWorkWm2 = 0.0;
    double upwindMixingWorkWm2 = 0.0;
};

// Experimental stationary tropical bulk momentum, with first-order upwind
// horizontal advection. Its coefficient is one inside 20 degrees and fades
// to zero at 35 degrees; the midlatitude/polar closure is retained exactly.
// Curvature is a skew rotation and does no work. The implicit local drag solve
// avoids an explicit advective CFL constraint, but iterations must be checked.
// This does not close a vertical momentum or prognostic pressure/energy budget.
inline SurfaceMomentumAdvection advectSurfaceMomentum(
    const climategrid::SphericalGrid& grid, const std::vector<float>& forceEast,
    const std::vector<float>& forceNorth, const std::vector<float>& drag,
    const std::vector<float>& initialEast, const std::vector<float>& initialSouth,
    const ModeSeparatedCirculationConfig& config, int maximumIterations,
    float relaxation = 0.5f)
{
    SurfaceMomentumAdvection result;
    const auto count = static_cast<std::size_t>(grid.columns) * grid.rows;
    const auto valid = [&](const auto& field) {
        return field.size() == count && std::all_of(field.begin(), field.end(), [](float x) { return std::isfinite(x); });
    };
    if (!climategrid::validGlobalGridDimensions(grid.columns, grid.rows) ||
        grid.latitudeCentresRadians.size() != static_cast<std::size_t>(grid.rows) ||
        grid.cellAreasSquareMetres.size() != static_cast<std::size_t>(grid.rows) ||
        !valid(forceEast) || !valid(forceNorth) || !valid(drag) || !valid(initialEast) || !valid(initialSouth) ||
        std::any_of(drag.begin(), drag.end(), [](float cd) { return cd <= 0.0f; }) ||
        !std::isfinite(config.surfaceBoundaryLayerDepthMetres) || config.surfaceBoundaryLayerDepthMetres <= 0.0f ||
        !std::isfinite(config.planetRadiusMetres) || config.planetRadiusMetres <= 0.0f ||
        !std::isfinite(grid.radiusMetres) || std::abs(grid.radiusMetres - config.planetRadiusMetres) > 1e-6 * config.planetRadiusMetres ||
        !std::isfinite(config.airDensityKgM3) || config.airDensityKgM3 <= 0.0f ||
        !std::isfinite(config.rotationRatePerSecond) || !std::isfinite(config.rotationDirection) || maximumIterations <= 0 ||
        !std::isfinite(relaxation) || relaxation <= 0.0f || relaxation > 1.0f)
        return result;
    result.east = initialEast; result.south = initialSouth;
    auto nextEast = initialEast, nextSouth = initialSouth;
    constexpr double pi = 3.141592653589793;
    std::vector<double> weight(grid.rows), dx(grid.rows), f(grid.rows), curvature(grid.rows);
    const double dy = pi * config.planetRadiusMetres / grid.rows;
    for (int y = 0; y < grid.rows; ++y)
    {
        const double latitude = grid.latitudeCentresRadians[y];
        const double t = std::clamp((std::abs(latitude * 180.0 / pi) - 20.0) / 15.0, 0.0, 1.0);
        weight[y] = 1.0 - t * t * (3.0 - 2.0 * t);
        dx[y] = 2.0 * pi * config.planetRadiusMetres * std::cos(latitude) / grid.columns;
        f[y] = coriolisParameterPerSecond(static_cast<float>(latitude * 180.0 / pi),
            config.rotationRatePerSecond, config.rotationDirection);
        curvature[y] = std::tan(latitude) / config.planetRadiusMetres;
    }
    for (int iteration = 0; iteration < maximumIterations; ++iteration)
    {
        double change = 0.0;
        for (int y = 0; y < grid.rows; ++y)
            for (int x = 0; x < grid.columns; ++x)
            {
                if (weight[y] == 0.0) continue;
                const auto i = grid.index(x, y);
                const double u = result.east[i], v = result.south[i];
                const auto ix = grid.index(x + (u >= 0.0 ? -1 : 1), y);
                const auto iy = grid.index(x, std::clamp(y + (v >= 0.0 ? -1 : 1), 0, grid.rows - 1));
                const double rx = weight[y] * std::abs(u) / dx[y], ry = weight[y] * std::abs(v) / dy;
                const auto wind = steadyMixedDragCoriolisWind(
                    static_cast<float>(forceEast[i] + rx * result.east[ix] + ry * result.east[iy]),
                    static_cast<float>(forceNorth[i] - rx * result.south[ix] - ry * result.south[iy]),
                    static_cast<float>(f[y] + weight[y] * u * curvature[y]), static_cast<float>(rx + ry),
                    drag[i] / config.surfaceBoundaryLayerDepthMetres);
                nextEast[i] = static_cast<float>(u + relaxation * (wind.eastMetresPerSecond - u));
                nextSouth[i] = static_cast<float>(v + relaxation * (wind.southMetresPerSecond - v));
                change = std::max(change, std::hypot(nextEast[i] - u, nextSouth[i] - v));
            }
        result.east.swap(nextEast); result.south.swap(nextSouth);
        result.iterations = iteration + 1; result.maximumChangeMps = change;
        if (change < 1e-5) break;
    }
    double residual2 = 0.0, forcing2 = 0.0, area = 0.0;
    const double mass = config.airDensityKgM3 * config.surfaceBoundaryLayerDepthMetres;
    for (int y = 0; y < grid.rows; ++y)
        for (int x = 0; x < grid.columns; ++x)
        {
            if (weight[y] == 0.0) continue;
            const auto i = grid.index(x, y);
            const double u = result.east[i], v = result.south[i], cellArea = grid.cellAreasSquareMetres[y];
            const auto ix = grid.index(x + (u >= 0.0 ? -1 : 1), y);
            const auto iy = grid.index(x, std::clamp(y + (v >= 0.0 ? -1 : 1), 0, grid.rows - 1));
            const double rx = weight[y] * std::abs(u) / dx[y], ry = weight[y] * std::abs(v) / dy;
            const double ae = rx * (u - result.east[ix]) + ry * (u - result.east[iy]);
            const double an = -rx * (v - result.south[ix]) - ry * (v - result.south[iy]);
            const double rate = drag[i] * std::hypot(u, v) / config.surfaceBoundaryLayerDepthMetres;
            const double fc = f[y] + weight[y] * u * curvature[y];
            const double re = rate * u + fc * v + ae - forceEast[i];
            const double rn = fc * u - rate * v + an - forceNorth[i];
            residual2 += cellArea * (re * re + rn * rn);
            forcing2 += cellArea * (forceEast[i] * forceEast[i] + forceNorth[i] * forceNorth[i]);
            area += cellArea;
            result.forcingWorkWm2 += cellArea * mass * (forceEast[i] * u - forceNorth[i] * v);
            result.dragWorkWm2 += cellArea * mass * rate * (u * u + v * v);
            result.advectionWorkWm2 += cellArea * mass * (ae * u - an * v);
            result.upwindMixingWorkWm2 += cellArea * mass * 0.5 * (rx *
                (std::pow(u - result.east[ix], 2) + std::pow(v - result.south[ix], 2)) + ry *
                (std::pow(u - result.east[iy], 2) + std::pow(v - result.south[iy], 2)));
        }
    result.relativeMomentumResidual = std::sqrt(residual2 / std::max(1e-30, forcing2));
    if (area > 0.0)
    {
        result.forcingWorkWm2 /= area; result.dragWorkWm2 /= area;
        result.advectionWorkWm2 /= area; result.upwindMixingWorkWm2 /= area;
    }
    return result;
}
}
