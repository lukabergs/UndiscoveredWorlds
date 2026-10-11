#include "simulations/landscape/opensimplex2s.hpp"

#include <cmath>

#include "core/deterministic_random.hpp"

// OpenSimplex2S in 3D after KdotJPG's reference (CC0); the source, the algorithm and
// the differences are in the header.

namespace landscape
{
namespace
{
// The reference's 48 gradients: the 24 sign and position variants of
// (1 + sqrt(3/2), 1 + sqrt(3/2), 1) and the 24 of (3.0862664687972017,
// 1.1721513422464978, 0), all of length 3.30136.
constexpr double gradients[48][3] = {
    { -2.22474487139, -2.22474487139, -1.0 },
    { -2.22474487139, -2.22474487139, 1.0 },
    { -3.0862664687972017, -1.1721513422464978, 0.0 },
    { -1.1721513422464978, -3.0862664687972017, 0.0 },
    { -2.22474487139, -1.0, -2.22474487139 },
    { -2.22474487139, 1.0, -2.22474487139 },
    { -1.1721513422464978, 0.0, -3.0862664687972017 },
    { -3.0862664687972017, 0.0, -1.1721513422464978 },
    { -2.22474487139, -1.0, 2.22474487139 },
    { -2.22474487139, 1.0, 2.22474487139 },
    { -3.0862664687972017, 0.0, 1.1721513422464978 },
    { -1.1721513422464978, 0.0, 3.0862664687972017 },
    { -2.22474487139, 2.22474487139, -1.0 },
    { -2.22474487139, 2.22474487139, 1.0 },
    { -1.1721513422464978, 3.0862664687972017, 0.0 },
    { -3.0862664687972017, 1.1721513422464978, 0.0 },
    { -1.0, -2.22474487139, -2.22474487139 },
    { 1.0, -2.22474487139, -2.22474487139 },
    { 0.0, -3.0862664687972017, -1.1721513422464978 },
    { 0.0, -1.1721513422464978, -3.0862664687972017 },
    { -1.0, -2.22474487139, 2.22474487139 },
    { 1.0, -2.22474487139, 2.22474487139 },
    { 0.0, -1.1721513422464978, 3.0862664687972017 },
    { 0.0, -3.0862664687972017, 1.1721513422464978 },
    { -1.0, 2.22474487139, -2.22474487139 },
    { 1.0, 2.22474487139, -2.22474487139 },
    { 0.0, 1.1721513422464978, -3.0862664687972017 },
    { 0.0, 3.0862664687972017, -1.1721513422464978 },
    { -1.0, 2.22474487139, 2.22474487139 },
    { 1.0, 2.22474487139, 2.22474487139 },
    { 0.0, 3.0862664687972017, 1.1721513422464978 },
    { 0.0, 1.1721513422464978, 3.0862664687972017 },
    { 2.22474487139, -2.22474487139, -1.0 },
    { 2.22474487139, -2.22474487139, 1.0 },
    { 1.1721513422464978, -3.0862664687972017, 0.0 },
    { 3.0862664687972017, -1.1721513422464978, 0.0 },
    { 2.22474487139, -1.0, -2.22474487139 },
    { 2.22474487139, 1.0, -2.22474487139 },
    { 3.0862664687972017, 0.0, -1.1721513422464978 },
    { 1.1721513422464978, 0.0, -3.0862664687972017 },
    { 2.22474487139, -1.0, 2.22474487139 },
    { 2.22474487139, 1.0, 2.22474487139 },
    { 1.1721513422464978, 0.0, 3.0862664687972017 },
    { 3.0862664687972017, 0.0, 1.1721513422464978 },
    { 2.22474487139, 2.22474487139, -1.0 },
    { 2.22474487139, 2.22474487139, 1.0 },
    { 3.0862664687972017, 1.1721513422464978, 0.0 },
    { 1.1721513422464978, 3.0862664687972017, 0.0 },
};

// The reference's NORMALIZER_3D: the supremum of the sum over positions and gradient
// assignments, so that the value divided by it lies in [-1, 1].
constexpr double normaliser = 0.2781926117527186;

// The kernel's squared radius (the reference's RSQUARED_3D).
constexpr double rsquared = 0.75;

// The contribution of the lattice point with hash coordinates (hx, hy, hz) at the
// offset (dx, dy, dz) from it.
double contribution(std::uint64_t seed, std::int64_t hx, std::int64_t hy, std::int64_t hz, double dx, double dy,
                    double dz)
{
    const double a = rsquared - dx * dx - dy * dy - dz * dz;
    if (!(a > 0.0))
        return 0.0;
    const double* g = gradients[deterministiccombine(seed, hx, hy, hz) % 48u];
    const double a2 = a * a;
    return a2 * a2 * (g[0] * dx + g[1] * dy + g[2] * dz);
}

// The contributions of the eight corners of the cube of one cubic lattice that holds
// (x, y, z), the lattice's points at integer + shift (shift 0 or 1/2), hashed at twice
// their coordinates (2 n + parity).
double cube(std::uint64_t seed, double x, double y, double z, double shift, std::int64_t parity)
{
    const double sx = x - shift, sy = y - shift, sz = z - shift;
    const double fx = std::floor(sx), fy = std::floor(sy), fz = std::floor(sz);
    const std::int64_t ix = static_cast<std::int64_t>(fx), iy = static_cast<std::int64_t>(fy),
                       iz = static_cast<std::int64_t>(fz);
    const double dx = sx - fx, dy = sy - fy, dz = sz - fz;
    double value = 0.0;
    for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++)
                value += contribution(seed, 2 * (ix + i) + parity, 2 * (iy + j) + parity, 2 * (iz + k) + parity,
                                      dx - i, dy - j, dz - k);
    return value;
}
}

double opensimplex2s(std::uint64_t seed, double x, double y, double z)
{
    // The reference's fallback orientation: an orthonormal half turn about (1, 1, 1).
    const double r = (2.0 / 3.0) * (x + y + z);
    const double xr = r - x, yr = r - y, zr = r - z;
    // The body-centred cubic lattice as its two cubic lattices.
    const double value = cube(seed, xr, yr, zr, 0.0, 0) + cube(seed, xr, yr, zr, 0.5, 1);
    return value / normaliser;
}
}
