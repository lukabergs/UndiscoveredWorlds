#pragma once

// Point-evaluated fractal gradient noise on the sphere (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2a): the noise relief_roughness.cpp has used since step C3, moved
// here so that the regional pass evaluates the same kind of noise at fine cell
// centres. A value is a pure function of the parameters and of the position on the
// unit sphere: no stream, no grid, no visiting order, so two windows that share a
// cell compute the same value there. Planet-free, except that a FractalNoise whose
// kernel and rotation are left at "configured" takes them from landscape.noise.kernel
// and landscape.noise.rotate_octaves when its field is built (multi-scale plan S3.16:
// the regional callers do not pass them yet).

#include <array>
#include <cstdint>
#include <vector>

namespace landscape
{
// The lattice noise of the octaves (S3.16): improved Perlin noise (gradientnoise) or
// OpenSimplex2S (opensimplex2s.hpp); configured reads landscape.noise.kernel (0 Perlin,
// 1 OpenSimplex2S).
enum class NoiseKernel
{
    configured = -1,
    perlin = 0,
    opensimplex2s = 1,
};

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

// The same for one octave of opensimplex2s, measured the same way (1024 x 512 cell
// centres, area-weighted, seeds 1 to 3, wavelengths 2000 and 500 km on a 6371 km
// sphere: 0.356 to 0.364, mean 0.3603; Perlin there 0.266 to 0.273;
// SphericalNoiseDeviationsAreTheMeasuredOnes checks both; S3.16).
constexpr double opensimplexdeviation = 0.360;

// The rotation of an octave (S3.16): the orthonormal matrix, row-major, of the unit
// quaternion drawn uniformly from three hashed unit numbers of the octave's seed
// (Shoemake 1992); a pure function of the seed, determinant +1.
std::array<double, 9> octaverotation(std::uint64_t octaveseed);

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
    NoiseKernel kernel = NoiseKernel::configured;
    int rotateoctaves = -1;      // 1: octave o's position is rotated by octaverotation(its seed) before
                                 // its scale and offset; 0: no rotation; -1: landscape.noise.rotate_octaves
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
    // The kernel's one-octave deviation (octavedeviation, or opensimplexdeviation)
    // x sqrt(sum of the squared full weights, resolved or not): the sum divided by it
    // has a standard deviation of about 1 when every octave has weight.
    double normalisation() const { return normalisation_; }
    // The kernel's one-octave deviation (octavedeviation or opensimplexdeviation):
    // relief_roughness.cpp forms its factor from it in the order of step C3.
    double deviation() const { return kernel_ == NoiseKernel::opensimplex2s ? opensimplexdeviation : octavedeviation; }
    // The sum of the squared full weights, resolved or not: relief_roughness.cpp
    // forms its factor from it in the order it always has (floating-point products
    // are not associative), which keeps its results bit-identical.
    double squaredweights() const { return squaredweights_; }
    // Number of octaves with a non-zero weight.
    int resolved() const { return resolved_; }

private:
    std::vector<std::uint64_t> seeds_;
    std::vector<double> scale_, weight_, offset_;
    std::vector<double> rotation_;  // 9 per octave when the octaves are rotated, else empty
    NoiseKernel kernel_ = NoiseKernel::perlin;
    double squaredweights_ = 0.0;
    double normalisation_ = 1.0;
    int resolved_ = 0;
};
}
