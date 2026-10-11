#pragma once

// OpenSimplex2S noise in 3D (multi-scale plan step S3.16; docs/architecture/noise-functions.md,
// section 4), the second kernel of FractalNoise beside the improved Perlin noise.
//
// Source: KdotJPG, OpenSimplex2 (github.com/KdotJPG/OpenSimplex2), java/OpenSimplex2S.java
// at commit e09a94744b3d69e4c58ce615445f18712cb50104 ("Revamp", 2022-01-16), the 3D
// function in its plain orientation (noise3_Fallback over noise3_UnrotatedBase); licence
// CC0 1.0 (public domain). The algorithm, RSQUARED_3D 3/4, the kernel (a^2)^2 g.d, the
// 48 gradients and NORMALIZER_3D 0.2781926117527186 are the reference's (the normaliser
// was also reproduced from the table and the kernel to 1e-16: the supremum of the sum
// over every gradient assignment, out/scratch/multiscale/s3.16/os2s_bound.py).
//
// The algorithm: the input is rotated by the orthonormal map r = 2/3 (x + y + z),
// (x, y, z) -> (r - x, r - y, r - z) (a half turn about the main diagonal); the lattice
// is the body-centred cubic one, the union of the integer lattice and the integer
// lattice shifted by (1/2, 1/2, 1/2); every lattice point within sqrt(3/4) of the
// position contributes (3/4 - d^2)^4 g.d, d the offset from the point, g one of the 48
// gradients. Differences from the reference, none of them in the function's form: the
// gradient is picked by deterministiccombine(seed, 2 i + s, 2 j + s, 2 k + s) % 48 (s 0
// on the integer lattice, 1 on the shifted one: no permutation table, as the Perlin
// kernel), not by the reference's prime-product hash into a 256-entry table cycling
// the 48 gradients; the second lattice is told apart by that parity s, not by the
// reference's flipped seed (seed ^ -0x52D547B2E96ED629) on a lattice offset by -1/2
// per axis (the same points, hashed differently); the
// arithmetic is double, not float; the contributing points are found as the corners of
// each lattice's containing cube that lie within the radius (the radius is half the
// cube's diagonal, so no other point reaches) instead of the reference's branches.
// A pure function of (seed, position).

#include <cstdint>

namespace landscape
{
// The value at (x, y, z), in [-1, 1] (the sum divided by its supremum over every
// gradient assignment); zero mean, continuous with continuous first derivatives.
double opensimplex2s(std::uint64_t seed, double x, double y, double z);
}
