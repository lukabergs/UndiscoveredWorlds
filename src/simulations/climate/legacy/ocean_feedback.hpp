#pragma once

namespace climategrid { struct SphericalGrid; }

namespace climateocean
{
struct OceanForcing;
struct OceanConfig;
struct OceanState;

namespace legacy
{
// Adds wind residuals in raster order so the caller's convergence arithmetic stays stable.
// Inputs and state use the validated grid/masks from solveWindDrivenOcean.
void applyEmpiricalSstWindFeedback(const climategrid::SphericalGrid& grid,
    const OceanForcing& forcing, const OceanConfig& config, float relaxation,
    OceanState& state, double& residual);
}
}
