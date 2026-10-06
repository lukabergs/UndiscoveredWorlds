#pragma once

// Relief roughness (docs/architecture/erosion-hydrology-plan.md, Phase C, step C3):
// a seeded, deterministic fractal perturbation of the land's elevation at the
// start of the erosion pass. Real plains carry tens of metres of relief that the
// smoothed tectonic terrain lacks; with it the routing branches and wanders and
// the erosion carves the wandering courses in. Planet-free; the process
// (terrain/fluvial_erosion.cpp) passes landscape.roughness.* and a seed derived
// from the world seed.

#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct ReliefRoughness
{
    double amplitude_m = 0.0;    // 0 turns the perturbation off
    double wavelength_km = 0.0;  // of the first octave
    int octaves = 1;
    double persistence = 0.5;    // amplitude ratio of successive octaves (each halves the wavelength)
    std::uint64_t seed = 0;
};

// Adds the perturbation to the land cells (sea == 0) of `elevation`. The noise is a
// function of the seed and of each cell centre's position on the sphere (so it has
// no seam at the date line and no stretch at the poles), zero-mean in expectation,
// with a peak magnitude of about amplitude_m when every octave is resolved. The
// grid resolves an octave whose wavelength is at least twice the north-south cell
// size; the others get weight 0 and the normalisation stays that of the full
// octave set, so the value is independent of the grid width up to the octaves the
// grid resolves (a coarse grid carries less relief). A lowering is scaled by
// min(1, (elevation - sealevel - 1 m) / amplitude_m), and a land cell is never
// lowered to sealevel + 1 m or below: it stays land. Sea cells are not changed.
// Deterministic and thread-count independent.
void addreliefroughness(const SphereMetric& metric, const ReliefRoughness& roughness, const std::vector<uint8_t>& sea,
                        float sealevel, std::vector<float>& elevation);
}
