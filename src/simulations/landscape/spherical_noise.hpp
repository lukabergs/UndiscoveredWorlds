#pragma once

// Point-evaluated fractal gradient noise on the sphere (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2a): the noise relief_roughness.cpp has used since step C3, moved
// here so that the regional pass evaluates the same kind of noise at fine cell
// centres. A value is a pure function of the parameters and of the position on the
// unit sphere: no stream, no grid, no visiting order, so two windows that share a
// cell compute the same value there. Planet-free.

#include <cstdint>
#include <vector>

namespace landscape
{
// Gradient noise (Perlin 2002) on the unit lattice of 3D space: zero at the lattice
// points, zero mean, continuous with continuous derivatives; the corner gradients
// are picked by deterministiccombine(seed, corner).
double gradientnoise(std::uint64_t seed, double x, double y, double z);

// The unit vector of a latitude and longitude in radians.
void unitvector(double latitude_rad, double longitude_rad, double& px, double& py, double& pz);

// Standard deviation of one octave of gradientnoise on the sphere (area-weighted
// over the 1024 x 512 cell centres, three seeds and two wavelengths: 0.266 to
// 0.271; a property of the noise, not a tuning value; step C3).
constexpr double octavedeviation = 0.268;

struct FractalNoise
{
    std::uint64_t seed = 0;
    double radius_km = 0.0;      // a position is the unit vector x radius_km / wavelength
    double wavelength_km = 0.0;  // of the first octave; octave o has wavelength_km / 2^o
    int octaves = 1;
    double persistence = 0.5;    // octave o weighs persistence^o
    double shortest_km = 0.0;    // an octave whose wavelength is below this gets weight 0
    double longest_km = 0.0;     // an octave whose wavelength is this or above gets weight 0 (0: no limit),
                                 // so bands [shortest, longest) that meet never share an octave
};

// The octave set of a FractalNoise with its seeds, scales, weights and lattice
// offsets computed once (octave o: seed deterministiccombine(seed, o), offset per
// axis 256 x unit(deterministiccombine(octave seed, axis, -1)), as relief_roughness.cpp
// computed them, so its results stay bit-identical).
class FractalNoiseField
{
public:
    explicit FractalNoiseField(const FractalNoise& noise);

    // The weighted sum of the octaves at a unit vector.
    double sum(double px, double py, double pz) const;
    // octavedeviation x sqrt(sum of the squared full weights, resolved or not): the
    // sum divided by it has a standard deviation of about 1 when every octave has
    // weight; octavedeviation is the measured deviation of one octave (0.268, C3).
    double normalisation() const { return normalisation_; }
    // The sum of the squared full weights, resolved or not: relief_roughness.cpp
    // forms its factor from it in the order it always has (floating-point products
    // are not associative), which keeps its results bit-identical.
    double squaredweights() const { return squaredweights_; }
    // Number of octaves with a non-zero weight.
    int resolved() const { return resolved_; }

private:
    std::vector<std::uint64_t> seeds_;
    std::vector<double> scale_, weight_, offset_;
    double squaredweights_ = 0.0;
    double normalisation_ = 1.0;
    int resolved_ = 0;
};
}
