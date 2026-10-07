#include "simulations/landscape/relief_roughness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "core/deterministic_random.hpp"

namespace landscape
{
namespace
{
// Standard deviation of one octave of the gradient noise below on the sphere
// (area-weighted over the 1024 x 512 cell centres, three seeds and two
// wavelengths: 0.266 to 0.271; a property of the noise, not a tuning value). The
// octave sum is divided by its own standard deviation and by `peakdeviations`,
// so that amplitude_m is the value reached at that many standard deviations: with
// every octave resolved, the largest magnitude on a grid of 64 to 1024 columns was
// 0.9 to 1.2 x amplitude_m (C3 measurement, before the band limit below).
constexpr double octavedeviation = 0.268;
constexpr double peakdeviations = 3.5;

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

// Gradient noise (Perlin 2002) on the unit lattice of 3D space: zero at the lattice
// points, zero mean, continuous with continuous derivatives.
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

// A uniform number in [0, 1) from a hash.
double unit(std::uint64_t hash) { return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); }
}

void addreliefroughness(const SphereMetric& metric, const ReliefRoughness& roughness, const std::vector<uint8_t>& sea,
                        float sealevel, std::vector<float>& elevation)
{
    if (roughness.amplitude_m == 0.0)
        return;
    // The noise lives on the sphere; a torus metric has no position on it (erosion plan D2a).
    if (metric.wrapsrows)
        throw std::invalid_argument("addreliefroughness: a torus metric is not supported");
    const std::size_t count = metric.size();
    if (sea.size() != count || elevation.size() != count)
        throw std::invalid_argument("addreliefroughness: input sizes differ from the metric");
    if (!(roughness.wavelength_km > 0.0) || roughness.octaves < 1)
        throw std::invalid_argument("addreliefroughness: wavelength_km must be positive and octaves at least 1");

    // Octave k samples the noise at the cell centre's position on the sphere in
    // units of its wavelength (wavelength_km / 2^k is the lattice spacing), with
    // its own seed and lattice offset; its weight is persistence^k. The weighted
    // sum is divided by the standard deviation of the full octave set,
    // octavedeviation x sqrt(sum of the squared weights). An octave whose
    // wavelength is below twice the north-south cell size is not resolved by the
    // grid (it would act as per-cell noise whose slope grows with the width): its
    // weight is 0, the normalisation unchanged.
    const int octaves = roughness.octaves;
    std::vector<std::uint64_t> seeds(static_cast<std::size_t>(octaves));
    std::vector<double> scale(static_cast<std::size_t>(octaves)), weight(static_cast<std::size_t>(octaves));
    std::vector<double> offset(3 * static_cast<std::size_t>(octaves));
    const double radius_km = metric.radius_m / 1000.0;
    const double shortest_km = 2.0 * metric.dy_m / 1000.0;
    double squaredweights = 0.0;
    int resolved = 0;
    for (int k = 0; k < octaves; k++)
    {
        const std::size_t o = static_cast<std::size_t>(k);
        const double wavelength = roughness.wavelength_km / std::ldexp(1.0, k);
        seeds[o] = deterministiccombine(roughness.seed, k);
        scale[o] = radius_km / wavelength;
        const double full = std::pow(roughness.persistence, k);
        squaredweights += full * full;
        weight[o] = wavelength >= shortest_km ? full : 0.0;
        resolved += weight[o] != 0.0 ? 1 : 0;
        for (int axis = 0; axis < 3; axis++)
            offset[3 * o + static_cast<std::size_t>(axis)] = 256.0 * unit(deterministiccombine(seeds[o], axis, -1));
    }
    if (resolved == 0)
        return;
    const double factor = roughness.amplitude_m / (peakdeviations * octavedeviation * std::sqrt(squaredweights));

    const double pi = 3.14159265358979323846;
    const double floor = static_cast<double>(sealevel) + 1.0;
    for (int y = 0; y < metric.height; y++)
    {
        const double latitude = metric.rowlatitude_rad[static_cast<std::size_t>(y)];
        const double coslat = std::cos(latitude), sinlat = std::sin(latitude);
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            if (sea[cell] != 0)
                continue;
            const double longitude = 2.0 * pi * (static_cast<double>(x) + 0.5) / metric.width;
            const double px = coslat * std::cos(longitude), py = coslat * std::sin(longitude), pz = sinlat;
            double sum = 0.0;
            for (int k = 0; k < octaves; k++)
            {
                const std::size_t o = static_cast<std::size_t>(k);
                if (weight[o] == 0.0)
                    continue;
                sum += weight[o] * gradientnoise(seeds[o], px * scale[o] + offset[3 * o], py * scale[o] + offset[3 * o + 1],
                                                 pz * scale[o] + offset[3 * o + 2]);
            }
            // A lowering tapers towards the coast: it is scaled by
            // min(1, (old - floor) / amplitude), so low land is not flattened onto
            // floor = sealevel + 1 m. The floor stays as a backstop for the rare
            // perturbation beyond the amplitude: a land cell is not lowered to it or
            // below (one already there is not lowered at all), so it stays land.
            const double old = static_cast<double>(elevation[cell]);
            double delta = factor * sum;
            if (delta < 0.0)
                delta *= std::clamp((old - floor) / std::fabs(roughness.amplitude_m), 0.0, 1.0);
            const double lowest = std::min(old, floor);
            elevation[cell] = static_cast<float>(std::max(old + delta, lowest));
        }
    }
}
}
