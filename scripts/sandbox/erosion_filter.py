# Phacelle Noise function copyright (c) 2025 Rune Skovbo Johansen
# Advanced Terrain Erosion Filter copyright (c) 2025 Rune Skovbo Johansen
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# NumPy port of PhacelleNoise and ErosionFilter from Runevision's "Fast and Gorgeous Erosion Filter"
# (Shadertoy wXcfWn; blog.runevision.com/2026/03/fast-and-gorgeous-erosion-filter.html), taken from the
# optolab port the blog links as faithful to the Shadertoy. Changes from the GLSL: every function works
# on arrays of points; the cell hash is an integer hash with a seed (the GLSL float hash degenerates at
# the large coordinates of a planet in kilometres); strength and the fade target may be arrays.
import numpy as np

TAU = 2.0 * np.pi


def cell_hash(ix, iy, seed):
    """Two values in [-1, 1) per integer cell (ix, iy), 16 bits each: the lowbias32 mix of the cell and the seed."""
    x = (np.asarray(ix) & 0xFFFFFFFF).astype(np.uint32)
    y = (np.asarray(iy) & 0xFFFFFFFF).astype(np.uint32)
    h = x * np.uint32(0x9E3779B1) ^ y * np.uint32(0x85EBCA77) ^ np.uint32((seed * 0x27D4EB2F + 0x165667B1) & 0xFFFFFFFF)
    h ^= h >> np.uint32(16)
    h *= np.uint32(0x7FEB352D)
    h ^= h >> np.uint32(15)
    h *= np.uint32(0x846CA68B)
    h ^= h >> np.uint32(16)
    scale = np.float32(1.0 / 32768.0)
    hx = (h & np.uint32(0xFFFF)).astype(np.float32) * scale - np.float32(1.0)
    hy = (h >> np.uint32(16)).astype(np.float32) * scale - np.float32(1.0)
    return hx, hy


HASH_X, HASH_Y, HASH_Z = np.uint32(0x9E3779B1), np.uint32(0x85EBCA77), np.uint32(0xC2B2AE3D)
HASH_SECOND = np.uint32(0x68E31DA4)
_MIX1, _MIX2 = np.uint32(0x7FEB352D), np.uint32(0x846CA68B)
_SHIFT15, _SHIFT16 = np.uint32(15), np.uint32(16)
_LOW16 = np.uint32(0xFFFF)
_SCALE, _ONE = np.float32(1.0 / 32768.0), np.float32(1.0)


def seed_term(seed):
    """cell_hash3's seed word."""
    return np.uint32((seed * 0x27D4EB2F + 0x165667B1) & 0xFFFFFFFF)


def as_word(values):
    """Integer lattice coordinates as their low 32 bits (uint32)."""
    return (np.asarray(values) & 0xFFFFFFFF).astype(np.uint32)


def mix_in_place(h, scratch):
    """The lowbias32 finaliser on uint32 h, in place (scratch: a uint32 array of h's shape)."""
    np.right_shift(h, _SHIFT16, out=scratch)
    np.bitwise_xor(h, scratch, out=h)
    np.multiply(h, _MIX1, out=h)
    np.right_shift(h, _SHIFT15, out=scratch)
    np.bitwise_xor(h, scratch, out=h)
    np.multiply(h, _MIX2, out=h)
    np.right_shift(h, _SHIFT16, out=scratch)
    np.bitwise_xor(h, scratch, out=h)


def hash_values(h, second):
    """The three values in [-1, 1) of a mixed word and its second mix."""
    return ((h & _LOW16).astype(np.float32) * _SCALE - _ONE, (h >> _SHIFT16).astype(np.float32) * _SCALE - _ONE,
            (second & _LOW16).astype(np.float32) * _SCALE - _ONE)


def hash_words(xa, yb, zc, seed_word, h, second, scratch):
    """h and second (in place) of the products x A, y B, z C of a cell: h = mix(x A ^ y B ^ z C ^ seed), second =
    mix(h ^ HASH_SECOND)."""
    np.bitwise_xor(xa, yb, out=h)
    np.bitwise_xor(h, zc, out=h)
    np.bitwise_xor(h, seed_word, out=h)
    mix_in_place(h, scratch)
    np.bitwise_xor(h, HASH_SECOND, out=second)
    mix_in_place(second, scratch)


def cell_hash3(ix, iy, iz, seed):
    """Three values in [-1, 1) per integer cell (ix, iy, iz), 16 bits each: two lowbias32 mixes of the cell and the seed."""
    xa = as_word(ix) * HASH_X
    yb = as_word(iy) * HASH_Y
    zc = as_word(iz) * HASH_Z
    xa, yb, zc = np.broadcast_arrays(xa, yb, zc)
    h, second, scratch = (np.empty(xa.shape, dtype=np.uint32) for _ in range(3))
    hash_words(xa, yb, zc, seed_term(seed), h, second, scratch)
    return hash_values(h, second)


def phacelle_noise3(p, direction, normal, freq, offset, normalization, seed):
    """phacelle_noise on a 3-D lattice (the sphere's positions, docs/architecture/sphere-mesh.md section 11): p (3, n)
    in lattice units, the unit stripe direction and the unit surface normal (3, n); the side vector is normal x
    direction (in the tangent plane), so the stripes lie on the surface. Returns (cos, sin, side (3, n))."""
    side = np.cross(normal, direction, axis=0) * (freq * TAU)
    offset = offset * TAU
    cells = np.floor(p)
    f = p - cells
    cells = cells.astype(np.int64)
    phase_c = np.zeros(p.shape[1:])
    phase_s = np.zeros(p.shape[1:])
    weight_sum = np.zeros(p.shape[1:])
    for i in range(-1, 3):
        for j in range(-1, 3):
            for k in range(-1, 3):
                hx, hy, hz = cell_hash3(cells[0] + i, cells[1] + j, cells[2] + k, seed)
                vx = f[0] - i - hx * 0.5
                vy = f[1] - j - hy * 0.5
                vz = f[2] - k - hz * 0.5
                weight = np.maximum(0.0, np.exp(-(vx * vx + vy * vy + vz * vz) * 2.0) - 0.01111)
                weight_sum += weight
                wave = vx * side[0] + vy * side[1] + vz * side[2] + offset
                phase_c += np.cos(wave) * weight
                phase_s += np.sin(wave) * weight
    weight_sum = np.maximum(weight_sum, 1e-30)  # a point far from every jittered lattice point (rare in 3-D)
    phase_c /= weight_sum
    phase_s /= weight_sum
    magnitude = np.maximum(1.0 - normalization, np.sqrt(phase_c * phase_c + phase_s * phase_s))
    return phase_c / magnitude, phase_s / magnitude, side


def erosion_filter3(position, height, slope, normal, fade_target, strength, gully_weight, detail, rounding, onset,
                    assumed_slope, scale, octaves, lacunarity, gain, cell_scale, normalization, seed,
                    octave_weights=None):
    """erosion_filter at 3-D positions (3, n; km) with the slope as a 3-D tangent vector (3, n) and the unit normal:
    the same octaves with phacelle_noise3. Returns (delta height, magnitude)."""
    strength = strength * scale
    fade_target = np.clip(fade_target, -1.0, 1.0)
    h = np.array(height, dtype=np.float64)
    s = np.array(slope, dtype=np.float64)
    h0 = h.copy()
    freq = 1.0 / (scale * cell_scale)
    slope_length = np.maximum(np.sqrt((s * s).sum(0)), 1e-10)
    magnitude = np.zeros_like(h)
    rounding_mult = 1.0

    rounding_input = (rounding[1] + (rounding[0] - rounding[1]) * _clamp01(fade_target + 0.5)) * rounding[2]
    combi_mask = ease_out(smooth_start(slope_length * onset[0], rounding_input * onset[0]))
    ridge_mask = ease_out(slope_length * onset[2])

    gully = s + (s / slope_length * assumed_slope[0] - s) * assumed_slope[1]

    for octave in range(octaves):
        weight = 1.0 if octave_weights is None else octave_weights[octave]
        length = np.sqrt((gully * gully).sum(0))
        safe = length > 1e-10
        direction = np.where(safe, gully / np.where(safe, length, 1.0), gully)
        c, sine, side = phacelle_noise3(position * freq, direction, normal, cell_scale, 0.25, normalization,
                                        seed + octave)
        side = side * -freq
        sloping = np.abs(sine)
        octave_strength = strength * weight
        gully = gully + np.sign(sine) * side * octave_strength * gully_weight

        gh, gs = c, sine * side
        faded_h = fade_target + (gh * gully_weight - fade_target) * combi_mask
        faded_s = gs * gully_weight * combi_mask
        h += faded_h * octave_strength
        s = s + faded_s * octave_strength
        magnitude += octave_strength

        fade_target = faded_h
        rounding_octave = (rounding[1] + (rounding[0] - rounding[1]) * _clamp01(c + 0.5)) * rounding_mult
        new_mask = ease_out(smooth_start(sloping * onset[1], rounding_octave * onset[1]))
        combi_mask = pow_inv(combi_mask, detail) * new_mask
        ridge_mask = ridge_mask * ease_out(sloping * onset[3])

        strength = strength * gain
        freq *= lacunarity
        rounding_mult *= rounding[3]
    return h - h0, magnitude


def phacelle_noise(px, py, dx, dy, freq, offset, normalization, seed):
    """Stripes along the unit direction (dx, dy) at the points (px, py): (cos, sin, side_x, side_y),
    the normalised cosine and sine waves and the side vector (freq tau across the stripes)."""
    side_x = -dy * freq * TAU
    side_y = dx * freq * TAU
    offset = offset * TAU
    ix, iy = np.floor(px), np.floor(py)
    fx, fy = px - ix, py - iy
    ix, iy = ix.astype(np.int64), iy.astype(np.int64)
    phase_c = np.zeros_like(px)
    phase_s = np.zeros_like(px)
    weight_sum = np.zeros_like(px)
    for i in range(-1, 3):
        for j in range(-1, 3):
            hx, hy = cell_hash(ix + i, iy + j, seed)
            vx = fx - i - hx * 0.5
            vy = fy - j - hy * 0.5
            weight = np.maximum(0.0, np.exp(-(vx * vx + vy * vy) * 2.0) - 0.01111)
            weight_sum += weight
            wave = vx * side_x + vy * side_y + offset
            phase_c += np.cos(wave) * weight
            phase_s += np.sin(wave) * weight
    phase_c /= weight_sum
    phase_s /= weight_sum
    magnitude = np.maximum(1.0 - normalization, np.sqrt(phase_c * phase_c + phase_s * phase_s))
    return phase_c / magnitude, phase_s / magnitude, side_x, side_y


def _clamp01(t):
    return np.clip(t, 0.0, 1.0)


def pow_inv(t, power):
    return 1.0 - np.power(1.0 - _clamp01(t), power)


def ease_out(t):
    v = 1.0 - _clamp01(t)
    return 1.0 - v * v


def smooth_start(t, smoothing):
    smoothing = np.maximum(smoothing, 1e-12)
    return np.where(t >= smoothing, t - 0.5 * smoothing, 0.5 * t * t / smoothing)


def _safe_normalize(x, y):
    length = np.sqrt(x * x + y * y)
    safe = length > 1e-10
    return np.where(safe, x / np.where(safe, length, 1.0), x), np.where(safe, y / np.where(safe, length, 1.0), y)


def erosion_filter(px, py, height, slope_x, slope_y, fade_target, strength, gully_weight, detail, rounding, onset,
                   assumed_slope, scale, octaves, lacunarity, gain, cell_scale, normalization, seed,
                   octave_weights=None):
    """ErosionFilter on arrays. Returns (delta height, delta slope x, delta slope y, magnitude, ridge map).
    octave_weights (one per octave, 0..1) scales each octave's strength: the sandbox fades octaves finer
    than its grid can show, so a coarse grid draws the same large gullies as a fine one."""
    strength = strength * scale
    fade_target = np.clip(fade_target, -1.0, 1.0)
    h, sx, sy = np.array(height, dtype=np.float64), np.array(slope_x, dtype=np.float64), np.array(slope_y, dtype=np.float64)
    h0, sx0, sy0 = h.copy(), sx.copy(), sy.copy()
    freq = 1.0 / (scale * cell_scale)
    slope_length = np.maximum(np.sqrt(sx * sx + sy * sy), 1e-10)
    magnitude = np.zeros_like(h)
    rounding_mult = 1.0

    rounding_input = (rounding[1] + (rounding[0] - rounding[1]) * _clamp01(fade_target + 0.5)) * rounding[2]
    combi_mask = ease_out(smooth_start(slope_length * onset[0], rounding_input * onset[0]))
    ridge_mask = ease_out(slope_length * onset[2])
    ridge_target = fade_target.copy()

    gully_x = sx + (sx / slope_length * assumed_slope[0] - sx) * assumed_slope[1]
    gully_y = sy + (sy / slope_length * assumed_slope[0] - sy) * assumed_slope[1]

    for octave in range(octaves):
        weight = 1.0 if octave_weights is None else octave_weights[octave]
        dx, dy = _safe_normalize(gully_x, gully_y)
        c, s, side_x, side_y = phacelle_noise(px * freq, py * freq, dx, dy, cell_scale, 0.25, normalization, seed + octave)
        side_x = side_x * -freq
        side_y = side_y * -freq
        sloping = np.abs(s)
        octave_strength = strength * weight
        gully_x = gully_x + np.sign(s) * side_x * octave_strength * gully_weight
        gully_y = gully_y + np.sign(s) * side_y * octave_strength * gully_weight

        gh, gsx, gsy = c, s * side_x, s * side_y
        faded_h = fade_target + (gh * gully_weight - fade_target) * combi_mask
        faded_sx = gsx * gully_weight * combi_mask
        faded_sy = gsy * gully_weight * combi_mask
        h += faded_h * octave_strength
        sx += faded_sx * octave_strength
        sy += faded_sy * octave_strength
        magnitude += octave_strength

        fade_target = faded_h
        rounding_octave = (rounding[1] + (rounding[0] - rounding[1]) * _clamp01(c + 0.5)) * rounding_mult
        new_mask = ease_out(smooth_start(sloping * onset[1], rounding_octave * onset[1]))
        combi_mask = pow_inv(combi_mask, detail) * new_mask

        ridge_target = ridge_target + (gh - ridge_target) * ridge_mask
        ridge_mask = ridge_mask * ease_out(sloping * onset[3])

        strength = strength * gain
        freq *= lacunarity
        rounding_mult *= rounding[3]

    ridge_map = ridge_target * (1.0 - ridge_mask)
    return h - h0, sx - sx0, sy - sy0, magnitude, ridge_map
