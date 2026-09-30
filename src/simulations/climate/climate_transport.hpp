#pragma once

#include "climate_hydrology.hpp"

namespace climatehydrology
{
struct FfslOptions
{
    // Limit deformation, not distance travelled. Prevents crossing departure
    // faces even during the intermediate multidimensional volume updates.
    double maximumDeformation = 0.45;
    // Scratch/geometry reuse only; false retains fresh allocation for checks.
    bool reuseWorkspace = true;
    std::vector<float> endZonalWindMps, endMeridionalWindMps;
};

// Conservative column density transport on the spherical finite-volume grid.
// SWIFT density splitting with limited linear subcell reconstruction.
SphericalTracerTransportDiagnostics advectSphericalTracerFfsl(
    int columns, int rows, const std::vector<float>& source,
    const std::vector<float>& zonalWindMps, const std::vector<float>& meridionalWindMps,
    float timeStepSeconds, float planetRadiusMetres, const FfslOptions& options,
    std::vector<float>& destination);
}
