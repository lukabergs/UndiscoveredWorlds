#pragma once

#include "climate_atmosphere.hpp"
#include "climate_grid.hpp"

#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
struct SurfaceMomentumForce
{
    float eastMps2 = 0.0f;
    float northMps2 = 0.0f;
};

// Recover on the native grid before interpolation: quadratic stress does not
// commute with averaging wind and drag separately across a coastline.
inline SurfaceMomentumForce surfaceMomentumForce(
    HorizontalWind wind, float cd, float latitudeDegrees,
    const ModeSeparatedCirculationConfig& config)
{
    if (cd <= 0.0f || config.surfaceBoundaryLayerDepthMetres <= 0.0f)
        return {};
    const double drag = cd * std::hypot(wind.eastMetresPerSecond,
        wind.southMetresPerSecond) / config.surfaceBoundaryLayerDepthMetres;
    const double f = coriolisParameterPerSecond(latitudeDegrees,
        config.rotationRatePerSecond, config.rotationDirection);
    return {static_cast<float>(drag * wind.eastMetresPerSecond + f * wind.southMetresPerSecond),
        static_cast<float>(f * wind.eastMetresPerSecond - drag * wind.southMetresPerSecond)};
}

// Smooth reconstruction spreads coarse-cell drag across fine coastlines.
// Recover the effective acceleration balanced by that reconstructed wind and
// drag, then retain it while applying the receiver's roughness. The force also
// includes the coarse stationary momentum exchange; it is not claimed to be
// a pure pressure gradient. No differentiation of reconstructed pressure or
// assumption about a 10 m versus layer-mean wind conversion is required.
inline HorizontalWind rediagnoseSurfaceWindForDrag(
    HorizontalWind wind, float sourceCd, float recipientCd, float latitudeDegrees,
    const ModeSeparatedCirculationConfig& config)
{
    if (sourceCd <= 0.0f || recipientCd <= 0.0f || config.surfaceBoundaryLayerDepthMetres <= 0.0f)
        return {};
    if (sourceCd == recipientCd ||
        (wind.eastMetresPerSecond == 0.0f && wind.southMetresPerSecond == 0.0f))
        return wind;
    const double drag = sourceCd * std::hypot(wind.eastMetresPerSecond,
        wind.southMetresPerSecond) / config.surfaceBoundaryLayerDepthMetres;
    const double f = coriolisParameterPerSecond(latitudeDegrees,
        config.rotationRatePerSecond, config.rotationDirection);
    const float eastForce = static_cast<float>(drag * wind.eastMetresPerSecond + f * wind.southMetresPerSecond);
    const float northForce = static_cast<float>(f * wind.eastMetresPerSecond - drag * wind.southMetresPerSecond);
    return steadyQuadraticDragCoriolisWind(eastForce, northForce, latitudeDegrees,
        recipientCd, config.surfaceBoundaryLayerDepthMetres,
        config.rotationRatePerSecond, config.rotationDirection);
}

// Delta in ascent from the changed lower-layer pressure-mass flux. Use half
// the surface equivalent pressure depth for the existing legacy two-mode
// ascent diagnosis. The same shared faces as diagnoseModeWinds ensure a zero
// global area integral, including the longitude seam and closed polar faces.
inline std::vector<float> surfaceAscentCorrectionHpaPerDay(
    const climategrid::SphericalGrid& grid,
    const std::vector<float>& deltaEastMps, const std::vector<float>& deltaSouthMps,
    float lowerPressureFluxDepthHpa)
{
    const auto cells = static_cast<std::size_t>(grid.columns) * grid.rows;
    if (!climategrid::validGlobalGridDimensions(grid.columns, grid.rows) ||
        deltaEastMps.size() != cells || deltaSouthMps.size() != cells)
        return {};
    std::vector<float> result(cells);
    for (int y = 0; y < grid.rows; ++y)
        for (int x = 0; x < grid.columns; ++x)
        {
            const auto cell = grid.index(x, y);
            const double east = 0.5 * (deltaEastMps[cell] + deltaEastMps[grid.index(x + 1, y)]) * grid.zonalFaceLengthsMetres[y];
            const double west = 0.5 * (deltaEastMps[cell] + deltaEastMps[grid.index(x - 1, y)]) * grid.zonalFaceLengthsMetres[y];
            const double south = y + 1 < grid.rows ?
                0.5 * (deltaSouthMps[cell] + deltaSouthMps[grid.index(x, y + 1)]) * grid.southFaceLengthsMetres[y] : 0.0;
            const double north = y > 0 ?
                0.5 * (deltaSouthMps[cell] + deltaSouthMps[grid.index(x, y - 1)]) * grid.northFaceLengthsMetres[y] : 0.0;
            result[cell] = static_cast<float>(-86400.0 * lowerPressureFluxDepthHpa *
                (east - west + south - north) / grid.cellAreasSquareMetres[y]);
        }
    return result;
}
}
