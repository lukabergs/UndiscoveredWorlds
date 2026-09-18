#pragma once
#include "climate_atmosphere.hpp"
#include "surface_wind_remap.hpp"
#include <algorithm>
#include <cmath>

namespace climateatmosphere::detail
{
struct Surface850Exchange
{
    HorizontalWind surface, lowLevel;
    double stressEastNm2 = 0.0, stressNorthNm2 = 0.0;
    double mixingLossWm2 = 0.0, maximumResidualMps2 = 0.0;
    bool valid = false;
};

// Stationary two-reservoir momentum with equal/opposite interfacial stress.
// Eliminate the linear low-level wind, then solve the surface quadratic balance.
// The 850 hPa state is a reduced bulk reservoir, not a resolved vertical profile.
inline Surface850Exchange coupleSurface850(SurfaceMomentumForce surfaceForce,
    SurfaceMomentumForce lowForce, float coriolis, float surfaceCd, float surfaceDepth,
    float surfaceDensity, float lowMassKgM2, float lowRayleighRate, float exchangeRate)
{
    Surface850Exchange result;
    for (float x : {surfaceForce.eastMps2, surfaceForce.northMps2, lowForce.eastMps2, lowForce.northMps2,
        coriolis, surfaceCd, surfaceDepth, surfaceDensity, lowMassKgM2, lowRayleighRate, exchangeRate})
        if (!std::isfinite(x)) return result;
    if (surfaceCd <= 0 || surfaceDepth <= 0 || surfaceDensity <= 0 || lowMassKgM2 <= 0 ||
        lowRayleighRate <= 0 || exchangeRate < 0) return result;
    const double ms = surfaceDepth * surfaceDensity, ma = lowMassKgM2;
    const double es = exchangeRate, ea = es * ms / ma;
    const double a = lowRayleighRate + ea, f = coriolis, denominator = a * a + f * f;
    const double linear = es - es * ea * a / denominator;
    const double effectiveF = f + es * ea * f / denominator;
    const double forceEast = surfaceForce.eastMps2 + es *
        (a * lowForce.eastMps2 + f * lowForce.northMps2) / denominator;
    const double forceNorth = surfaceForce.northMps2 + es *
        (a * lowForce.northMps2 - f * lowForce.eastMps2) / denominator;
    result.surface = steadyMixedDragCoriolisWind(static_cast<float>(forceEast), static_cast<float>(forceNorth),
        static_cast<float>(effectiveF), static_cast<float>(std::max(0.0, linear)), surfaceCd / surfaceDepth);
    const double us = result.surface.eastMetresPerSecond, ns = -result.surface.southMetresPerSecond;
    const double ae = lowForce.eastMps2 + ea * us, an = lowForce.northMps2 + ea * ns;
    result.lowLevel = {static_cast<float>((a * ae + f * an) / denominator),
        static_cast<float>((f * ae - a * an) / denominator)};
    const double ua = result.lowLevel.eastMetresPerSecond, na = -result.lowLevel.southMetresPerSecond;
    result.stressEastNm2 = ms * es * (ua - us);
    result.stressNorthNm2 = ms * es * (na - ns);
    result.mixingLossWm2 = ms * es * ((ua - us) * (ua - us) + (na - ns) * (na - ns));
    const double rs = surfaceCd * std::hypot(us, ns) / surfaceDepth;
    result.maximumResidualMps2 = std::max(
        std::hypot(rs * us - f * ns - es * (ua - us) - surfaceForce.eastMps2,
                   rs * ns + f * us - es * (na - ns) - surfaceForce.northMps2),
        std::hypot(lowRayleighRate * ua - f * na + ea * (ua - us) - lowForce.eastMps2,
                   lowRayleighRate * na + f * ua + ea * (na - ns) - lowForce.northMps2));
    result.valid = std::isfinite(result.maximumResidualMps2);
    return result;
}
}
