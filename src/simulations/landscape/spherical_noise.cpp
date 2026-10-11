#include "simulations/landscape/spherical_noise.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/opensimplex2s.hpp"

// The noise of relief_roughness.cpp (step C3), moved here in erosion plan step E2a
// (docs/architecture/erosion-hydrology-plan.md, Phase E, section 6, "Noise band"):
// every operation is the one it was there, in the same order, so the roughness's
// results stay bit-identical. Multi-scale plan S3.16 adds the OpenSimplex2S kernel
// and the per-octave rotation; with the Perlin kernel and no rotation every
// operation is still the one of C3.

namespace landscape
{
namespace
{
double fade(double t) { return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }

double lerp(double a, double b, double t) { return a + t * (b - a); }

// The dot product of the lattice corner's gradient with the offset from it: one of
// the twelve cube-edge directions (Perlin 2002), picked by a hash of the seed and
// the corner.
double cornerdot(std::uint64_t seed, int ix, int iy, int iz, double dx, double dy, double dz)
{
    switch (deterministiccombine(seed, ix, iy, iz) % 12u)
    {
    case 0: return dx + dy;
    case 1: return -dx + dy;
    case 2: return dx - dy;
    case 3: return -dx - dy;
    case 4: return dx + dz;
    case 5: return -dx + dz;
    case 6: return dx - dz;
    case 7: return -dx - dz;
    case 8: return dy + dz;
    case 9: return -dy + dz;
    case 10: return dy - dz;
    default: return -dy - dz;
    }
}

// A uniform number in [0, 1) from a hash.
double unit(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }
}

double gradientnoise(std::uint64_t seed, double x, double y, double z)
{
    const double fx = std::floor(x), fy = std::floor(y), fz = std::floor(z);
    const int ix = static_cast<int>(fx), iy = static_cast<int>(fy), iz = static_cast<int>(fz);
    const double dx = x - fx, dy = y - fy, dz = z - fz;
    const double u = fade(dx), v = fade(dy), w = fade(dz);
    const double x00 = lerp(cornerdot(seed, ix, iy, iz, dx, dy, dz),
                            cornerdot(seed, ix + 1, iy, iz, dx - 1.0, dy, dz), u);
    const double x10 = lerp(cornerdot(seed, ix, iy + 1, iz, dx, dy - 1.0, dz),
                            cornerdot(seed, ix + 1, iy + 1, iz, dx - 1.0, dy - 1.0, dz), u);
    const double x01 = lerp(cornerdot(seed, ix, iy, iz + 1, dx, dy, dz - 1.0),
                            cornerdot(seed, ix + 1, iy, iz + 1, dx - 1.0, dy, dz - 1.0), u);
    const double x11 = lerp(cornerdot(seed, ix, iy + 1, iz + 1, dx, dy - 1.0, dz - 1.0),
                            cornerdot(seed, ix + 1, iy + 1, iz + 1, dx - 1.0, dy - 1.0, dz - 1.0), u);
    return lerp(lerp(x00, x10, v), lerp(x01, x11, v), w);
}

void unitvector(double latitude_rad, double longitude_rad, double& px, double& py, double& pz)
{
    const double coslat = std::cos(latitude_rad);
    px = coslat * std::cos(longitude_rad);
    py = coslat * std::sin(longitude_rad);
    pz = std::sin(latitude_rad);
}

std::array<double, 9> octaverotation(std::uint64_t octaveseed)
{
    // Shoemake's uniform unit quaternion from three uniform numbers, normalised again
    // against rounding, and its rotation matrix.
    const double twopi = 6.28318530717958647692;
    const double u1 = unit(deterministiccombine(octaveseed, 0, -2));
    const double u2 = unit(deterministiccombine(octaveseed, 1, -2));
    const double u3 = unit(deterministiccombine(octaveseed, 2, -2));
    double w = std::sqrt(1.0 - u1) * std::sin(twopi * u2);
    double x = std::sqrt(1.0 - u1) * std::cos(twopi * u2);
    double y = std::sqrt(u1) * std::sin(twopi * u3);
    double z = std::sqrt(u1) * std::cos(twopi * u3);
    const double length = std::sqrt(w * w + x * x + y * y + z * z);
    w /= length;
    x /= length;
    y /= length;
    z /= length;
    return { 1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),       2.0 * (x * z + y * w),
             2.0 * (x * y + z * w),       1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w),
             2.0 * (x * z - y * w),       2.0 * (y * z + x * w),       1.0 - 2.0 * (x * x + y * y) };
}

FractalNoiseField::FractalNoiseField(const FractalNoise& noise)
{
    if (!(noise.wavelength_km > 0.0) || noise.octaves < 1 || !(noise.radius_km >= 0.0))
        throw std::invalid_argument(
            "FractalNoiseField: wavelength_km must be positive, octaves at least 1 and radius_km not negative");
    // The kernel and the rotation, from landscape.noise when left at configured.
    const int kernel = noise.kernel == NoiseKernel::configured ? tuning::landscape::noise::kernel
                                                               : static_cast<int>(noise.kernel);
    if (kernel != static_cast<int>(NoiseKernel::perlin) && kernel != static_cast<int>(NoiseKernel::opensimplex2s))
        throw std::invalid_argument("FractalNoiseField: the kernel must be 0 (Perlin) or 1 (OpenSimplex2S)");
    kernel_ = static_cast<NoiseKernel>(kernel);
    const int rotate = noise.rotateoctaves < 0 ? tuning::landscape::noise::rotateOctaves : noise.rotateoctaves;
    if (rotate != 0 && rotate != 1)
        throw std::invalid_argument("FractalNoiseField: rotate_octaves must be 0 or 1");

    // Octave k samples the noise at the position on the sphere in units of its
    // wavelength (wavelength_km / 2^k is the lattice spacing), with its own seed and
    // lattice offset; its weight is persistence^k, or 0 outside the band
    // shortest_km <= wavelength < longest_km (the upper end exclusive, so that a
    // regional band ending where the global roughness's begins counts no octave
    // twice; longest_km 0 is no upper end). The normalisation is that of the full octave set,
    // octavedeviation x sqrt(sum of the squared weights), whatever the band.
    const int octaves = noise.octaves;
    const std::size_t count = static_cast<std::size_t>(octaves);
    seeds_.resize(count);
    scale_.resize(count);
    weight_.resize(count);
    offset_.resize(3 * count);
    if (rotate == 1)
        rotation_.resize(9 * count);
    double squaredweights = 0.0;
    int resolved = 0;
    for (int k = 0; k < octaves; k++)
    {
        const std::size_t o = static_cast<std::size_t>(k);
        const double wavelength = noise.wavelength_km / std::ldexp(1.0, k);
        seeds_[o] = deterministiccombine(noise.seed, k);
        scale_[o] = noise.radius_km / wavelength;
        const double full = std::pow(noise.persistence, k);
        squaredweights += full * full;
        const bool inband = wavelength >= noise.shortest_km && (noise.longest_km <= 0.0 || wavelength < noise.longest_km);
        weight_[o] = inband ? full : 0.0;
        resolved += weight_[o] != 0.0 ? 1 : 0;
        for (int axis = 0; axis < 3; axis++)
            offset_[3 * o + static_cast<std::size_t>(axis)] = 256.0 * unit(deterministiccombine(seeds_[o], axis, -1));
        if (!rotation_.empty())
        {
            const std::array<double, 9> matrix = octaverotation(seeds_[o]);
            for (std::size_t i = 0; i < 9; i++)
                rotation_[9 * o + i] = matrix[i];
        }
    }
    squaredweights_ = squaredweights;
    normalisation_ = deviation() * std::sqrt(squaredweights);
    resolved_ = resolved;
}

double FractalNoiseField::sum(double px, double py, double pz) const
{
    // Accumulated octave by octave in order; an octave of weight 0 is skipped.
    double total = 0.0;
    const std::size_t count = weight_.size();
    if (kernel_ == NoiseKernel::perlin && rotation_.empty())
    {
        // The operations of step C3.
        for (std::size_t o = 0; o < count; o++)
        {
            if (weight_[o] == 0.0)
                continue;
            total += weight_[o] * gradientnoise(seeds_[o], px * scale_[o] + offset_[3 * o], py * scale_[o] + offset_[3 * o + 1],
                                                pz * scale_[o] + offset_[3 * o + 2]);
        }
        return total;
    }
    for (std::size_t o = 0; o < count; o++)
    {
        if (weight_[o] == 0.0)
            continue;
        // The position rotated (a unit vector stays one), then scaled and offset.
        double qx = px, qy = py, qz = pz;
        if (!rotation_.empty())
        {
            const double* m = &rotation_[9 * o];
            qx = m[0] * px + m[1] * py + m[2] * pz;
            qy = m[3] * px + m[4] * py + m[5] * pz;
            qz = m[6] * px + m[7] * py + m[8] * pz;
        }
        const double x = qx * scale_[o] + offset_[3 * o], y = qy * scale_[o] + offset_[3 * o + 1],
                     z = qz * scale_[o] + offset_[3 * o + 2];
        total += weight_[o] * (kernel_ == NoiseKernel::opensimplex2s ? opensimplex2s(seeds_[o], x, y, z)
                                                                     : gradientnoise(seeds_[o], x, y, z));
    }
    return total;
}
}
