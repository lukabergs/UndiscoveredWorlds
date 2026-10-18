#pragma once

// The terrain stages' noise and legacy-pixel geometry on the cell grid
// (docs/architecture/sphere-mesh.md sections 3, 11 and 12; terrain.sphere_noise.*):
// every noise value is a function of the seed and the cell centre's unit vector
// (landscape::FractalNoiseField), so it is the same on any grid; a legacy pixel
// constant (a radius, a distance, a warp) is a length in km of the 2048-column world.

#include <cstdint>
#include <vector>

#include "core/mesh/cell_grid.hpp"

class planet;

namespace terrain
{
// One caller's replacement of a legacy diamond-square fractal: the first wavelength (km),
// the persistence, and the legacy field's mean and deviation as shares of its maximum.
struct SphereFractal
{
    double wavelength_km = 0.0;
    double persistence = 0.5;
    double mean = 0.5;
    double deviation = 0.2;
};

// The callers' fractals (terrain.sphere_noise.<caller>_*).
SphereFractal terrainfractal();
SphereFractal smoothingfractal();
SphereFractal canyonbasefractal();
SphereFractal canyonmaskfractal();
SphereFractal archipelagofractal();

// One legacy pixel in km: 2 pi R / terrain.sphere_noise.legacy_columns.
double legacycellkm(const planet& world);
// The number of legacy pixels of the whole world: legacy_columns x legacy_columns / 2.
double legacycellcount();

// The fractal at every cell centre, in id order: maximum x clamp(mean + deviation x z, 0, 1),
// at least 1, z the noise's normalised sum (about unit deviation).
std::vector<float> spherefractal(const planet& world, std::uint64_t seed, const SphereFractal& fractal, float maximum);

// The same, sampled where the legacy warp moved it: each cell's point moved east and north by
// (w - mean) x warpfactor x vary legacy pixels, w and vary the warp and warp_vary fields
// (terrain.sphere_noise.warp_*, as shares of their maximum).
std::vector<float> spherefractalwarped(const planet& world, std::uint64_t seed, const SphereFractal& fractal,
                                       float maximum, int warpfactor);

// The cells within radius_km of a cell's centre (mesh::cellswithin; the cell itself at least).
std::vector<mesh::CellId> cellswithinkm(const planet& world, mesh::CellId c, double radius_km);
// The great-circle distance in km between two cell centres.
double distancekm(const planet& world, mesh::CellId a, mesh::CellId b);
// The cell reached from c by moving east_km and north_km in its tangent frame.
mesh::CellId movedcell(const planet& world, mesh::CellId c, double east_km, double north_km);
// A cell drawn by area from the shared legacy stream (two 30-bit draws of random()).
mesh::CellId randomcellbyarea(const planet& world);
}
