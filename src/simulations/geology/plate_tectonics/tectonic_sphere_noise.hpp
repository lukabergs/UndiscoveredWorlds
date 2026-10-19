#ifndef TECTONIC_SPHERE_NOISE_HPP
#define TECTONIC_SPHERE_NOISE_HPP
// Noise and smooth patterns on the sphere for the engine's position-placed processes (W2-E, E4;
// docs/architecture/sphere-mesh.md §11): the seeding noise, the ocean relief, the noise pass and
// the process variation. Every value is a pure function of a unit vector, the seed and lengths
// in km, so it is continuous across the cube's face edges and corners and independent of the
// grid, its resolution and the visiting order.
#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/spherical_noise.hpp"
#include <array>
#include <cstdint>

namespace tectonics::sphere_noise {
// One octave's standard deviation of the lattice noises the torus engine used, measured over
// 4e6 samples (E4, 2026-10-19): the 4-D simplex noise of the seeding (0.297) and the ocean
// relief's periodic value noise (0.455). An Octaves field scales its octaves (one octave's
// deviation landscape::opensimplexdeviation 0.360) to one of them, so the thresholds and
// amplitudes written for the torus noise keep their meaning. Properties of the noises, not
// tuning values.
inline constexpr double simplex_deviation = 0.297;
inline constexpr double value_deviation = 0.455;

// Fractal noise on the sphere (OpenSimplex2S, each octave rotated by its seed: no lattice axis
// shared across octaves or tied to the pole): octave o has wavelength wavelength_km / 2^o and
// weight persistence^o; an octave shorter than shortest_km is left out (the first is always
// kept, at least shortest_km long). The value is the weighted mean of the kept octaves times
// deviation / opensimplexdeviation, clamped to [-1, 1] (the torus noises' range).
class Octaves {
public:
    Octaves(std::uint64_t seed, double wavelength_km, double radius_km, int octaves, double persistence,
            double deviation = simplex_deviation, double shortest_km = 0);
    double operator()(const mesh::Vec3& p) const;

private:
    landscape::FractalNoiseField field_;
    double factor_ = 1;
};

// p moved along the sphere by the tangent part of amplitude_km x (x(p), y(p), z(p)), as a
// great-circle arc of that length on a sphere of radius_km.
mesh::Vec3 warp(const mesh::Vec3& p, const Octaves& x, const Octaves& y, const Octaves& z, double amplitude_km,
                double radius_km);

// A seeded orthonormal frame (the rows of landscape::octaverotation of the seed): the axes of
// the plane-wave patterns.
std::array<mesh::Vec3, 3> frame(std::uint64_t seed);
}
#endif
