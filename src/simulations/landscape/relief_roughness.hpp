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
// pure function of the seed and of each cell centre's position on the sphere (so
// it has no seam at the date line and no stretch at the poles, and does not depend
// on the grid width beyond its sampling), zero-mean, with a peak magnitude of about
// amplitude_m. A land cell stays land: it is not lowered to sealevel + 1 m or below
// by the perturbation. Sea cells are not changed. Deterministic and thread-count
// independent.
void addreliefroughness(const SphereMetric& metric, const ReliefRoughness& roughness, const std::vector<uint8_t>& sea,
                        float sealevel, std::vector<float>& elevation);
}
