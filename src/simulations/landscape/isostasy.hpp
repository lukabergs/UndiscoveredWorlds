#pragma once

// Local isostatic response of the erosion run (docs/architecture/erosion-hydrology-plan.md,
// Phase C, step C8): local Airy compensation, no flexure. Removing a thickness e of
// rock of density rho_c lets the column rise by e x rho_c / rho_m; laying a thickness
// d of sediment of density rho_s presses it down by d x rho_s / rho_m. Planet-free;
// the process (terrain/fluvial_erosion.cpp) passes landscape.isostasy.*.

#include <vector>

namespace landscape
{
struct IsostasyParameters
{
    double reboundfraction = 0.0;  // uplift per metre of rock removed; 0 turns the rebound off
    double loadingfraction = 0.0;  // subsidence per metre of sediment laid; 0 turns the loading off
};

// Applies one step's response to every cell (land and sea; a marine deposit loads the
// sea floor): elevation += reboundfraction x erosion_m - loadingfraction x deposition_m.
// `erosion_m` and `deposition_m` are the step's non-negative thicknesses removed and
// laid (landscape::stepstreampower); each cell is independent of its neighbours.
// Throws std::invalid_argument on a size mismatch.
void applyisostasy(const IsostasyParameters& parameters, const std::vector<float>& erosion_m,
                   const std::vector<float>& deposition_m, std::vector<float>& elevation);
}
