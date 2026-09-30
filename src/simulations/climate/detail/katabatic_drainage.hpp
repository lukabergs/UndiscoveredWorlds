#pragma once

#include "surface_wind_remap.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateatmosphere::detail
{
// Diagnostic jet increment: p(z)=z(H-z)/[10(H-10)] within a shallow layer.
// It is one at 10 m, zero at the ground/top, and peaks at H/2. This assumed
// shape is not a resolved turbulence/stability model; the carrier is unchanged.
inline double coldJetProfileWeight(double height, double depth)
{
    if (!(depth > 20.0) || !(height > 0.0) || height >= depth) return 0.0;
    return height * (depth - height) / (10.0 * (depth - 10.0));
}

inline double coldJetMeanToTenMetres(double depth)
{
    return depth > 20.0 ? depth * depth / (60.0 * (depth - 10.0)) : 0.0;
}

inline float coldJetTransportFraction(double depth, double transportDepth)
{
    if (!(transportDepth >= depth) || !(depth > 20.0)) return 0.0f;
    return static_cast<float>(coldJetMeanToTenMetres(depth) * depth / transportDepth);
}

// Subtract the unchanged carrier momentum balance from the shallow balance:
// Cd/H (|U10|U10-|U0|U0) + f J mean(p) (U10-U0) = cold buoyancy.
// The same profile sets Coriolis and column transport. Surface stress uses
// the diagnosed 10 m wind; the thermal capacity uses this same depth H.
// Ambient forcing is anchored to the carrier, not claimed to be resolved
// storm forcing. No cold reservoir advection or entrainment is solved here.
inline HorizontalWind coldJetTenMetreWind(HorizontalWind carrier,
    SurfaceMomentumForce cold, float latitude, float cdTenMetres, float depth,
    float rotationRate, float rotationDirection)
{
    if (cold.eastMps2 == 0.0f && cold.northMps2 == 0.0f) return carrier;
    const double mean = coldJetMeanToTenMetres(depth);
    if (!(mean > 0.0) || !(cdTenMetres > 0.0f)) return carrier;
    const double drag = cdTenMetres * std::hypot(carrier.eastMetresPerSecond,
        carrier.southMetresPerSecond) / depth;
    const double f = mean * coriolisParameterPerSecond(latitude, rotationRate, rotationDirection);
    return steadyQuadraticDragCoriolisWind(
        static_cast<float>(drag * carrier.eastMetresPerSecond + f * carrier.southMetresPerSecond + cold.eastMps2),
        static_cast<float>(f * carrier.eastMetresPerSecond - drag * carrier.southMetresPerSecond + cold.northMps2),
        latitude, cdTenMetres, depth, static_cast<float>(mean * rotationRate), rotationDirection);
}

// Diagnostic velocity increment, normalized at the surface and tapered to zero
// at the top of the represented boundary layer. This is an assumed profile,
// not a resolved jet. Uniform density gives its volume/mass mean below.
inline double katabaticProfileWeight(double height, double decayHeight, double layerDepth)
{
    if (!(decayHeight > 0.0) || !(layerDepth > 0.0)) return 0.0;
    const double z = std::clamp(height, 0.0, layerDepth);
    const double q = layerDepth / decayHeight;
    if (q < 1e-5) return 1.0 - z / layerDepth;
    return std::exp(-z / decayHeight) *
        (-std::expm1(-(layerDepth - z) / decayHeight)) / (-std::expm1(-q));
}

inline float katabaticLayerMeanFraction(double decayHeight, double layerDepth)
{
    if (!(decayHeight > 0.0) || !(layerDepth > 0.0)) return 0.0f;
    const double q = layerDepth / decayHeight;
    const double mean = q < 1e-5 ? 0.5 - q / 12.0 :
        1.0 / q - (q > 700.0 ? 0.0 : 1.0 / std::expm1(q));
    return static_cast<float>(std::clamp(mean, 0.0, 0.5));
}

// Ice occupies its diagnosed fraction; additional seasonal snow covers the
// remaining surface. This avoids counting snow on ice twice.
inline float frozenLandAlbedo(float snow, float ice, float openAlbedo,
    float snowAlbedo, float iceAlbedo)
{
    ice = std::clamp(ice, 0.0f, 1.0f);
    const float snowOnly = std::max(0.0f, std::clamp(snow, 0.0f, 1.0f) - ice);
    return openAlbedo + snowOnly * (snowAlbedo - openAlbedo) + ice * (iceAlbedo - openAlbedo);
}

// Reduced cold-layer closure: the layer-mean potential-temperature deficit
// is a fraction of the co-located air-minus-skin temperature difference.
// Both temperatures refer to the local surface, avoiding a second lapse-rate
// correction. Snow cover gates drainage; this does not model piteraq gusts.
inline float katabaticDeficitK(float airC, float skinC, float snowCover,
    float deficitFraction, float maximumDeficitK)
{
    if (!std::isfinite(airC) || !std::isfinite(skinC) || !std::isfinite(snowCover) ||
        !(deficitFraction > 0.0f) || !(maximumDeficitK > 0.0f)) return 0.0f;
    const float frozen = std::clamp(-skinC / 2.0f, 0.0f, 1.0f);
    return frozen * std::clamp(snowCover, 0.0f, 1.0f) *
        std::clamp(deficitFraction * (airC - skinC), 0.0f, maximumDeficitK);
}

// A mixed cell can contain melting snow while its mean skin is above freezing.
// Estimate the remaining snow patch at min(mean skin, melting point). This is
// a subgrid temperature proxy, not an alteration of the surface energy state.
inline float snowPatchKatabaticDeficitK(float airC, float skinC, float snowCover,
    float deficitFraction, float maximumDeficitK)
{
    if (!std::isfinite(airC) || !std::isfinite(skinC) || !std::isfinite(snowCover) ||
        !(deficitFraction > 0) || !(maximumDeficitK > 0)) return 0;
    return std::clamp(snowCover, 0.0f, 1.0f) *
        std::clamp(deficitFraction * (airC - std::min(skinC, 0.0f)), 0.0f, maximumDeficitK);
}

// Reduced inversion estimate from radiative cooling over a finite renewal time:
// Q * tau / (rho * cp * h) has units K. This diagnoses unresolved cold-layer
// buoyancy; it is not an extra surface heat sink or a prognostic temperature.
inline float radiativeKatabaticDeficitK(float netSurfaceRadiativeWm2, float skinC,
    float snowCover, float densityKgM3, float heatCapacityJkgK, float depthMetres,
    float renewalSeconds, float coolingFraction, float maximumDeficitK)
{
    if (!std::isfinite(netSurfaceRadiativeWm2) || !std::isfinite(skinC) || !std::isfinite(snowCover) ||
        !(densityKgM3 > 0) || !(heatCapacityJkgK > 0) || !(depthMetres > 0) ||
        !(renewalSeconds > 0) || !(coolingFraction > 0) || !(maximumDeficitK > 0)) return 0;
    const float deficit = std::max(0.0f, -netSurfaceRadiativeWm2) * renewalSeconds * coolingFraction /
        (densityKgM3 * heatCapacityJkgK * depthMetres);
    return std::clamp(-skinC / 2.0f, 0.0f, 1.0f) * std::clamp(snowCover, 0.0f, 1.0f) *
        std::clamp(deficit, 0.0f, maximumDeficitK);
}

inline SurfaceMomentumForce katabaticForce(float eastSlope, float northSlope,
    float temperatureK, float deficitK, float gravityMps2, float maximumSlope)
{
    if (!std::isfinite(eastSlope) || !std::isfinite(northSlope) ||
        !std::isfinite(temperatureK) || !std::isfinite(deficitK) ||
        !std::isfinite(gravityMps2) || !(temperatureK > 0.0f) ||
        !(deficitK > 0.0f) || !(gravityMps2 > 0.0f) || !(maximumSlope > 0.0f)) return {};
    const float slope = std::hypot(eastSlope, northSlope);
    const float limit = slope > maximumSlope ? maximumSlope / slope : 1.0f;
    const float acceleration = -gravityMps2 * deficitK / temperatureK * limit;
    return {acceleration * eastSlope, acceleration * northSlope};
}

// Metric height gradients on the solver grid: cyclic longitude, one-sided
// derivatives at the first/last latitude centres; no division by cos(90 deg).
struct TerrainSlope { float east = 0.0f, north = 0.0f; };

inline std::vector<TerrainSlope> terrainSlopes(
    const climategrid::SphericalGrid& grid, const std::vector<float>& heightMetres,
    bool reflectedPolarGradient = false)
{
    if (!climategrid::validGlobalGridDimensions(grid.columns, grid.rows) || grid.rows < 2 ||
        heightMetres.size() != static_cast<std::size_t>(grid.columns * grid.rows) ||
        !(grid.radiusMetres > 0.0)) return {};
    std::vector<TerrainSlope> result(heightMetres.size());
    const double dy = grid.radiusMetres * grid.latitudeSpacingRadians;
    for (int y = 0; y < grid.rows; ++y)
    {
        const int north = std::max(0, y - 1), south = std::min(grid.rows - 1, y + 1);
        const double dx = grid.radiusMetres * grid.longitudeSpacingRadians * std::cos(grid.latitudeCentresRadians[y]);
        for (int x = 0; x < grid.columns; ++x)
        {
            // A scalar continued across a pole is sampled half a turn away.
            // This restores a centred derivative instead of doubling the
            // pole-adjacent slope of a smooth axisymmetric dome.
            const int nx = reflectedPolarGradient && y == 0 ? x + grid.columns / 2 : x;
            const int sx = reflectedPolarGradient && y == grid.rows - 1 ? x + grid.columns / 2 : x;
            const double span = reflectedPolarGradient ? 2.0 : south - north;
            result[grid.index(x,y)] = {
                static_cast<float>((heightMetres[grid.index(x+1,y)] - heightMetres[grid.index(x-1,y)]) / (2.0*dx)),
                static_cast<float>((heightMetres[grid.index(nx,north)] - heightMetres[grid.index(sx,south)]) / (span*dy))};
        }
    }
    return result;
}
}
