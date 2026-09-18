#pragma once

#include "../climate_ocean_dynamics.hpp"
#include <algorithm>
#include <cmath>

namespace climateocean::detail
{
// Reduced thermal-only bulk closure, not a density/salinity or turbulence model.
// Wind work scales with water-side u*^3; destabilizing thermal buoyancy flux
// mixes against a specified background N^2. Warming limits the layer by a
// Monin-Obukhov length. All forcing is simulated; there are no geographic masks.
inline double storageColumnDepth(double bathymetry, const OceanConfig& c)
{
    return std::max(1.0, std::min(bathymetry, double(c.storageColumnDepthMetres)));
}
inline double minimumStorageDepth(double column, const OceanConfig& c)
{
    return std::min(double(c.minimumStorageDepthMetres), 0.5 * column);
}
inline double maximumStorageDepth(double column, const OceanConfig& c)
{
    return std::min(double(c.maximumStorageDepthMetres), 0.9 * column);
}
inline double storageThermalStratification(double surfaceTemperature, double reservoirTemperature,
    double column, const OceanConfig& c)
{
    // The two homogeneous layers' centres are D/2 apart, independent of h.
    // A weak unresolved background remains when thermal contrast vanishes;
    // this model has no prognostic salinity or freshwater buoyancy.
    return std::max(double(c.storageStratificationPerSecond2),
        9.81 * 2.0e-4 * std::max(0.0, surfaceTemperature - reservoirTemperature) / (0.5 * column));
}
inline double storageTargetDepth(double windSpeed, double heatFluxWm2, double column, const OceanConfig& c,
    double stratification = -1.0)
{
    const double ustar = std::sqrt(c.airDensityKgM3 * c.dragCoefficient / c.waterDensityKgM3) * windSpeed;
    const double work = c.storageWindMixingEfficiency * ustar * ustar * ustar;
    const double buoyancy = 9.81 * 2.0e-4 * heatFluxWm2 / (c.waterDensityKgM3 * c.waterHeatCapacityJkgK);
    const double memory = c.storageMixingMemoryDays * 86400.0;
    const double n2 = stratification > 0.0 ? stratification : c.storageStratificationPerSecond2;
    const double windDepth = std::cbrt(12.0 * work * memory / n2);
    const double convectiveDepth2 = 2.0 * std::max(0.0, -buoyancy) * memory / n2;
    const double target = buoyancy > 0.0 ? std::min(windDepth, 2.0 * work / buoyancy)
        : std::sqrt(windDepth * windDepth + convectiveDepth2);
    return std::clamp(target, minimumStorageDepth(column, c), maximumStorageDepth(column, c));
}

// E and R are J/m2 relative to water at freezing. The lower reservoir stores
// detrained water; deepening entrains its actual retained heat. Latent ice stays
// at the surface, and warm entrainment can melt it. E+R is exactly conserved
// apart from floating-point roundoff; h+(D-h) conserves column water volume.
inline double remixStorage(double newDepth, double column, double volumetricCapacity,
    double& depth, double& enthalpy, double& reservoirEnthalpy)
{
    const double transfer = newDepth >= depth
        ? (newDepth - depth) * reservoirEnthalpy / (column - depth)
        : (newDepth - depth) * std::max(0.0, enthalpy) / depth;
    enthalpy += transfer;
    reservoirEnthalpy -= transfer;
    depth = newDepth;
    // Convective adjustment also prevents a warm retained reservoir beneath a
    // colder slab, including the heat needed to melt an ice-covered slab.
    const double total = enthalpy + reservoirEnthalpy;
    const double surfaceTheta = std::max(0.0, enthalpy) / (volumetricCapacity * depth);
    const double lowerTheta = reservoirEnthalpy / (volumetricCapacity * (column - depth));
    double adjustment = 0.0;
    if (lowerTheta > surfaceTheta)
    {
        const double adjusted = total > 0.0 ? total * depth / column : total;
        adjustment = adjusted - enthalpy;
        enthalpy = adjusted;
        reservoirEnthalpy = total - adjusted;
    }
    return transfer + adjustment;
}
}
