"""Terrain sandbox: procedural terrain laid on the tectonic world, then the pipeline's drainage and erosion
(terrain_sandbox.exe, src/cli/terrain_sandbox_main.cpp), rendered on the whole world and on windows.
It is for iterating on the look in seconds, not a pipeline stage: nothing is registered and the multiscale
contract (exact block means, pinned rivers, seamless windows) is not enforced. The base is a world bundle
stopped before the climate (out/scratch/sandbox/base: procedural-world, seed 1, 256 columns, 600 Myr,
--steps ..terrain.roughness).

Order per domain (the whole world at a preview resolution, each window at a scale k with a margin):
the base surface (cubic B-spline of the world's surface elevation), a slow position warp, the relief noise
(isotropic, blended into noise stretched along the strike where the crust is shortened), the erosion
filter (erosion_filter.py) on base + relief, a rock-hardness field (relative erodibility, banded along the
strike), then routing and stream power. Every noise is a function of absolute position, so a window is a
zoom of the world preview; octaves finer than a domain's grid fade out instead of aliasing.

  uv run python scripts/sandbox/terrain_sandbox.py --label try1 [--config FILE] [--set group.key=value ...]
      [--domains global,A,B] [--no-physics]
Writes out/scratch/sandbox/<label>/: index.html, the PNGs, config.json (the values used) and timings.
"""
import argparse
import copy
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'evaluation'))
from uwtools.bundle import Bundle  # noqa: E402
from uwtools.map_rendering import colourize, save_rgb  # noqa: E402
from uwtools.palettes import palette_rgb  # noqa: E402
from validation import surface_elevation  # noqa: E402
from maps import hillshade, shade_factor, terrain_base  # noqa: E402

import erosion_filter  # noqa: E402

DEFAULT_CONFIG = ROOT / 'configs' / 'sandbox' / 'terrain.json'
PHYSICS = ROOT / 'out' / 'build' / 'headless' / 'bin' / 'Release' / 'terrain_sandbox.exe'
NOISE_STD = 0.22  # standard deviation of one octave of the gradient noise below, measured


# ---------------------------------------------------------------- configuration

def load_config(path, overrides):
    config = json.loads(Path(path).read_text(encoding='utf-8'))
    for item in overrides:
        key, _, value = item.partition('=')
        node = config
        parts = key.split('.')
        for part in parts[:-1]:
            node = node[part]
        try:
            node[parts[-1]] = json.loads(value)
        except json.JSONDecodeError:
            node[parts[-1]] = value
    return config


# ---------------------------------------------------------------- domains

class Domain:
    """A grid of fine cells in global-cell units: x east from the date line, y south from the north pole;
    cell centres at (column + 0.5) / k. km: kilometres per global cell (east-west at the equator and
    north-south alike on a cell-centred W x W/2 grid)."""

    def __init__(self, name, columns, rows, x0, y0, k, world_columns, world_rows, radius_km, crop=None,
                 exaggeration=10.0, physics=None):
        self.name, self.k = name, k
        self.columns, self.rows = columns, rows
        self.x = x0 + (np.arange(columns) + 0.5) / k
        self.y = y0 + (np.arange(rows) + 0.5) / k
        self.world_columns, self.world_rows = world_columns, world_rows
        self.cell_km = 2.0 * math.pi * radius_km / world_columns
        self.fine_km = self.cell_km / k
        self.radius_km = radius_km
        self.north_deg = 90.0 - y0 * 180.0 / world_rows
        self.pixel_deg = 180.0 / world_rows / k
        latitude = self.north_deg - (np.arange(rows) + 0.5) * self.pixel_deg
        self.coslat = np.maximum(np.cos(np.radians(latitude)), 1e-3)
        self.crop = crop  # (row0, row1, column0, column1) of the kept window, or None
        self.whole = crop is None
        self.exaggeration = exaggeration
        self.physics = physics or {}

    def mesh_km(self):
        """Planar kilometres of every cell (equirectangular: east-west at the equator's scale)."""
        return np.meshgrid(self.x * self.cell_km, self.y * self.cell_km)


def make_domains(config, bundle, names):
    radius_km = float(bundle.manifest.get('parameters', {}).get('world.body.radius_km', 6371.0))
    W, H = bundle.columns, bundle.rows
    domains = []
    for name in names:
        spec = config['domains'][name]
        if name == 'global':
            k = spec['columns'] // W
            domains.append(Domain(name, W * k, H * k, 0.0, 0.0, k, W, H, radius_km,
                                  exaggeration=spec.get('exaggeration', 10.0), physics=config['physics']['global']))
        else:
            c0, r0, c1, r1 = spec['pixels']
            k, m = spec['scale'], spec['margin']
            columns, rows = (c1 - c0 + 1 + 2 * m) * k, (r1 - r0 + 1 + 2 * m) * k
            crop = (m * k, rows - m * k, m * k, columns - m * k)
            domains.append(Domain(name, columns, rows, c0 - m, r0 - m, k, W, H, radius_km, crop=crop,
                                  exaggeration=spec.get('exaggeration', 3.0), physics=config['physics']['window']))
    return domains


# ---------------------------------------------------------------- sampling the world

def bspline_matrix(positions, size, wrap):
    """Rows of cubic B-spline weights from `size` cell-centred samples to `positions` (cell units)."""
    t = positions - 0.5
    i0 = np.floor(t).astype(np.int64)
    f = t - i0
    weights = [(1 - f) ** 3 / 6, (3 * f ** 3 - 6 * f ** 2 + 4) / 6, (-3 * f ** 3 + 3 * f ** 2 + 3 * f + 1) / 6, f ** 3 / 6]
    matrix = np.zeros((positions.size, size))
    for offset, w in zip(range(-1, 3), weights):
        index = i0 + offset
        index = np.mod(index, size) if wrap else np.clip(index, 0, size - 1)
        np.add.at(matrix, (np.arange(positions.size), index), w)
    return matrix


def sample(field, domain, mx=None, my=None):
    mx = bspline_matrix(domain.x, domain.world_columns, True) if mx is None else mx
    my = bspline_matrix(domain.y, domain.world_rows, False) if my is None else my
    return my @ np.asarray(field, dtype=np.float64) @ mx.T


def world_fields(bundle):
    elevation = surface_elevation(bundle)
    sea = bundle.load('geosphere.base_elevation').astype(float) <= float(bundle.scalar('geosphere.sea_level'))
    # The strain as the fold fabric reads it: the shortening (m) over the native cell width 2 pi R / W.
    radius_m = float(bundle.manifest.get('parameters', {}).get('world.body.radius_km', 6371.0)) * 1000.0
    strain = bundle.load('geosphere.shortening').astype(float) / (2.0 * math.pi * radius_m / bundle.columns)
    direction = bundle.load('geosphere.shortening_direction').astype(float)
    uplift = bundle.load('geosphere.rock_uplift_rate').astype(float) if bundle.has('geosphere.rock_uplift_rate') else \
        np.zeros_like(elevation)
    return {'elevation': elevation, 'sea': sea, 'shortening': np.maximum(strain, 0.0), 'direction': direction,
            'uplift': uplift}


# ---------------------------------------------------------------- noise

def gradient_noise(x, y, seed):
    """2D gradient noise with a quintic fade; values about +-0.7, standard deviation NOISE_STD."""
    ix, iy = np.floor(x), np.floor(y)
    fx, fy = (x - ix).astype(np.float32), (y - iy).astype(np.float32)
    ix, iy = ix.astype(np.int64), iy.astype(np.int64)
    ux = fx * fx * fx * (fx * (fx * 6 - 15) + 10)
    uy = fy * fy * fy * (fy * (fy * 6 - 15) + 10)
    corners = []
    for cx, cy in ((0, 0), (1, 0), (0, 1), (1, 1)):
        gx, gy = erosion_filter.cell_hash(ix + cx, iy + cy, seed)
        corners.append(gx * (fx - cx) + gy * (fy - cy))
    a, b, c, d = corners
    return a + ux * (b - a) + uy * (c - a) + ux * uy * (a - b - c + d)


def octave_fade(wavelength_km, domain, min_cells):
    """1 for an octave the grid draws with at least 2 min_cells cells per wavelength, 0 below min_cells,
    linear between: the same octave fades at the same wavelength on every domain."""
    cells = wavelength_km / domain.fine_km
    return float(np.clip((cells - min_cells) / min_cells, 0.0, 1.0))


def fbm(x, y, spec, seed, domain, stretch=1.0):
    """Fractal sum with unit standard deviation (about); x across and y along a stretch of `stretch`.
    Each octave is rotated by a fixed irrational angle so no lattice direction repeats."""
    total = np.zeros_like(x)
    norm = 0.0
    wavelength, amplitude = spec['wavelength_km'], 1.0
    for octave in range(spec['octaves']):
        fade = octave_fade(wavelength, domain, spec.get('min_cells', 2.5))
        if fade > 0.0:
            if stretch == 1.0:  # each octave rotated by the golden angle
                angle = 2.399963 * octave
                ca, sa = math.cos(angle), math.sin(angle)
                u, v = (ca * x - sa * y) / wavelength, (sa * x + ca * y) / wavelength
            else:  # a stretched octave keeps the stretch's axes: shifted, not rotated
                u, v = x / wavelength + 0.37 * octave, y / (wavelength * stretch) + 0.61 * octave
            total += amplitude * fade * gradient_noise(u, v, seed * 131 + octave)
        norm += (amplitude * fade) ** 2
        wavelength /= spec.get('lacunarity', 2.0)
        amplitude *= spec.get('gain', 0.5)
    return total / (NOISE_STD * math.sqrt(norm)) if norm > 0 else total


def oriented_fbm(x, y, theta, alignment, spec, seed, domain):
    """Isotropic fBm blended into a bank of fBm stretched along the strike: bank member b has its across-
    strike axis at phi_b = b pi / B, weight max(0, cos 2(theta - phi_b))^sharpness; `alignment` 0..1 mixes
    the isotropic sum into the bank's. theta: the shortening axis (radians from east toward south)."""
    isotropic = fbm(x, y, spec, seed, domain)
    count = spec.get('orientations', 6)
    if count == 0 or not np.any(alignment > 0):
        return isotropic
    bank = np.zeros_like(x)
    weights = np.zeros_like(x)
    aligned = alignment > 0
    for b in range(count):
        phi = math.pi * b / count
        w = np.maximum(0.0, np.cos(2.0 * (theta - phi))) ** spec.get('sharpness', 4.0)
        active = aligned & (w > 1e-3)  # evaluated only where the member carries weight
        if not np.any(active):
            continue
        xa, ya = x[active], y[active]
        across = math.cos(phi) * xa + math.sin(phi) * ya
        along = -math.sin(phi) * xa + math.cos(phi) * ya
        bank[active] += w[active] * fbm(across, along, spec, seed * 7 + 101 + b, domain, stretch=spec.get('stretch', 4.0))
        weights[active] += w[active] ** 2
    bank = bank / np.sqrt(np.maximum(weights, 1e-12))
    a = np.clip(alignment, 0.0, 1.0)
    return ((1 - a) * isotropic + a * bank) / np.sqrt((1 - a) ** 2 + a ** 2)


def ridged(n):
    """Unit-variance noise to sharp-crested noise of mean about 0 (|n| of a unit normal has mean 0.8)."""
    return 0.8 - np.abs(n)


# ---------------------------------------------------------------- the procedural pass

def procedural(domain, world, config, seed):
    timings = {}
    started = time.perf_counter()
    mx = bspline_matrix(domain.x, domain.world_columns, True)
    my = bspline_matrix(domain.y, domain.world_rows, False)
    base = sample(world['elevation'], domain, mx, my)
    strain = sample(world['shortening'], domain, mx, my)
    doubled = world['shortening'] * np.exp(2j * world['direction'])
    doubled = sample(doubled.real, domain, mx, my) + 1j * sample(doubled.imag, domain, mx, my)
    theta = 0.5 * np.angle(doubled)
    coherence = np.abs(doubled) / np.maximum(strain, 1e-9)
    structure = config['structure']
    alignment = np.clip(strain / structure['strain_full'], 0.0, 1.0) * np.clip(coherence, 0.0, 1.0)
    uplift = sample(world['uplift'], domain, mx, my)
    sea_seed = sample(world['sea'].astype(float), domain, mx, my) > 0.5
    timings['base'] = time.perf_counter() - started

    started = time.perf_counter()
    x, y = domain.mesh_km()
    # The warp reads the strain smoothed over warp.smoothing_cells world cells: a drag that changes over one
    # 156 km cell folds the grid over itself (v4's loops along the belts).
    sigma = config['warp']['smoothing_cells']
    smooth = lambda field: ndimage.gaussian_filter(field, sigma, mode=('nearest', 'wrap'))  # noqa: E731
    warp_strain = sample(smooth(world['shortening']), domain, mx, my)
    warp_doubled = world['shortening'] * np.exp(2j * world['direction'])
    warp_doubled = sample(smooth(warp_doubled.real), domain, mx, my) + 1j * sample(smooth(warp_doubled.imag), domain, mx, my)
    warp_alignment = np.clip(warp_strain / structure['strain_full'], 0.0, 1.0) * \
        np.clip(np.abs(warp_doubled) / np.maximum(warp_strain, 1e-9), 0.0, 1.0)
    x, y = warped(x, y, 0.5 * np.angle(warp_doubled), warp_alignment, config['warp'], seed, domain)
    stages = {'base': base, 'warp': (x, y)}
    relief = config['relief']
    noise = oriented_fbm(x, y, theta, alignment, {**relief, **structure}, seed + 20, domain)
    noise = (1 - relief['ridged']) * noise + relief['ridged'] * ridged(noise) / 0.6
    amplitude = relief['base_m'] + relief['elevation_fraction'] * np.maximum(base, 0.0) + \
        relief['strain_m'] * alignment
    amplitude = np.where(base > 0, amplitude, amplitude * relief['sea_fraction'])
    surface = base + amplitude * noise
    stages['relief'] = surface
    uplands = config['uplands']
    if uplands['lift_m'] > 0:  # continental interiors raised in patches: a non-negative lift, none at the coast
        land_share = sample(ndimage.gaussian_filter((~world['sea']).astype(float), uplands['coast_sigma_cells'],
                                                    mode=('nearest', 'wrap')), domain, mx, my)
        inland = np.clip((land_share - 0.5) * 2.0, 0.0, 1.0)
        large = fbm(x, y, uplands, seed + 50, domain)
        surface = surface + uplands['lift_m'] * inland * np.clip(0.5 + 0.5 * large, 0.0, 1.5) * (base > 0)
    stages['uplands'] = surface
    hardness_spec = config['hardness']
    banding = oriented_fbm(x, y, theta, alignment, {**hardness_spec, **structure}, seed + 30, domain)
    contrast = hardness_spec['contrast']
    erodibility = np.exp(hardness_spec['sigma'] * np.tanh(contrast * banding) / math.tanh(contrast))
    stages['hardness'] = erodibility
    timings['noise'] = time.perf_counter() - started

    started = time.perf_counter()
    filtered = surface
    spec = config['filter']
    if spec['enabled']:
        filtered = surface + erosion_pass(domain, surface, x, y, amplitude, erodibility, spec, seed)
    stages['filter'] = filtered
    timings['filter'] = time.perf_counter() - started

    land_or_low = filtered <= 0.0
    labels, _ = ndimage.label(land_or_low, structure=np.ones((3, 3)))
    sea_labels = np.unique(labels[(labels > 0) & sea_seed])
    sea = np.isin(labels, sea_labels[sea_labels > 0])
    if not domain.whole:  # the outer ring drains the window: sea at or below -1 m
        ring = np.zeros_like(sea)
        ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
        sea |= ring
        filtered = np.where(ring, np.minimum(filtered, -1.0), filtered)
    return {'base': base, 'surface': filtered, 'sea': sea, 'erodibility': erodibility, 'uplift': uplift,
            'strain': strain, 'stages': stages}, timings


def warped(x, y, theta, alignment, spec, seed, domain):
    """Positions moved by the warp. Two noise components v = (n1, n2) give a background displacement
    amplitude_km v everywhere and, in the belts, strain_amplitude_km x alignment x (s (s.v) + across_fraction
    c (c.v)), s along the strike and c across it: rock dragged along the belt, a little squeezed across it.
    The projections s (s.v) and c (c.v) do not change when the axis flips sign, so the displacement stays
    continuous where the shortening angle wraps from pi to 0."""
    if spec['amplitude_km'] <= 0 and spec['strain_amplitude_km'] <= 0:
        return x, y
    n1 = fbm(x, y, spec, seed + 11, domain)
    n2 = fbm(x, y, spec, seed + 12, domain)
    cx, cy = np.cos(theta), np.sin(theta)  # across the strike (the shortening axis)
    sx, sy = -cy, cx                        # along the strike
    along = sx * n1 + sy * n2
    across = cx * n1 + cy * n2
    belt = spec['strain_amplitude_km'] * np.clip(alignment, 0.0, 1.0)
    dx = spec['amplitude_km'] * n1 + belt * (sx * along + spec['across_fraction'] * cx * across)
    dy = spec['amplitude_km'] * n2 + belt * (sy * along + spec['across_fraction'] * cy * across)
    return x + dx, y + dy


def erosion_pass(domain, height_m, x, y, amplitude_m, erodibility, spec, seed):
    """The erosion filter on the domain: horizontal unit km, vertical unit km x vertical_scale; the
    slope by central differences; the fade target the height between the local floor and crest (a
    min/max filter over fade_radius_km, smoothed), -1 at the floor and +1 at the crest; the strength
    scaled by the local relief amplitude over relief_reference_m and by the erodibility^hardness_exponent."""
    vs = spec['vertical_scale']
    h = height_m / 1000.0 * vs
    dy_km = domain.fine_km
    dx_km = domain.fine_km * domain.coslat[:, None]
    gy, gx = np.gradient(h)
    slope_x, slope_y = gx / dx_km, gy / dy_km
    radius = max(1, int(round(spec['fade_radius_km'] / domain.fine_km)))
    mode = ('nearest', 'wrap') if domain.whole else 'nearest'
    low = ndimage.uniform_filter(ndimage.minimum_filter(h, size=2 * radius + 1, mode=mode), size=2 * radius + 1, mode=mode)
    high = ndimage.uniform_filter(ndimage.maximum_filter(h, size=2 * radius + 1, mode=mode), size=2 * radius + 1, mode=mode)
    fade = 2.0 * (h - low) / np.maximum(high - low, 1e-6) - 1.0
    land = np.clip(height_m / spec['land_fade_m'], 0.0, 1.0)  # no gullies at sea, faded in above the coast
    strength = spec['strength'] * np.clip(amplitude_m / spec['relief_reference_m'], 0.0, spec['relief_cap']) * \
        erodibility ** spec['hardness_exponent'] * land
    weights, wavelength = [], spec['scale_km']
    for _ in range(spec['octaves']):
        weights.append(octave_fade(wavelength, domain, spec.get('min_cells', 2.5)))
        wavelength /= spec['lacunarity']
    dh, _, _, magnitude, _ = erosion_filter.erosion_filter(
        x, y, h, slope_x, slope_y, fade, strength, spec['gully_weight'], spec['detail'], spec['rounding'],
        spec['onset'], spec['assumed_slope'], spec['scale_km'], spec['octaves'], spec['lacunarity'], spec['gain'],
        spec['cell_scale'], spec['normalization'], seed + 40, octave_weights=weights)
    offset_x, offset_y = spec['height_offset']
    offset = (offset_x + (-np.clip(fade, -1, 1) - offset_x) * offset_y) * magnitude
    return (dh + offset) / vs * 1000.0


# ---------------------------------------------------------------- physics

def run_physics(domain, fields, config, directory):
    spec = {**config['physics'], **domain.physics}
    stem = directory / f'{domain.name}'
    files = {}
    for key, values in (('elevation', fields['surface']), ('sea', fields['sea'].astype(np.float32)),
                        ('erodibility', fields['erodibility']), ('uplift', fields['uplift'] * spec['uplift_factor'])):
        files[key] = stem.with_name(f'{domain.name}-{key}.f32')
        np.asarray(values, dtype='<f4').tofile(files[key])
    out_elevation = stem.with_name(f'{domain.name}-out-elevation.f32')
    out_discharge = stem.with_name(f'{domain.name}-out-discharge.f32')
    command = [str(PHYSICS), '--columns', str(domain.columns), '--rows', str(domain.rows),
               '--elevation', str(files['elevation']), '--sea', str(files['sea']),
               '--erodibility', str(files['erodibility']), '--uplift', str(files['uplift']),
               '--runoff', str(spec['runoff_m_yr']), '--duration', str(spec['duration_myr']),
               '--step', str(spec['step_myr']), '--k', str(spec['k']), '--m', str(spec['m']),
               '--deposition', str(spec['deposition']), '--diffusivity', str(spec['diffusivity_km2_myr']),
               '--delta-top', str(spec['delta_top_m']),
               '--radius-km', str(domain.radius_km),
               '--out-elevation', str(out_elevation), '--out-discharge', str(out_discharge)]
    if not domain.whole:
        command += ['--north', str(domain.north_deg), '--pixel-deg', str(domain.pixel_deg)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'terrain_sandbox failed on {domain.name}: {result.stderr.strip()}')
    shape = (domain.rows, domain.columns)
    elevation = np.fromfile(out_elevation, dtype='<f4').reshape(shape).astype(np.float64)
    discharge = np.fromfile(out_discharge, dtype='<f4').reshape(shape).astype(np.float64)
    for path in list(files.values()) + [out_elevation, out_discharge]:
        path.unlink()
    return elevation, discharge, result.stdout.strip()


# ---------------------------------------------------------------- rendering

def crop(domain, values):
    if domain.crop is None:
        return values
    r0, r1, c0, c1 = domain.crop
    return values[r0:r1, c0:c1]


def relief_rgb(domain, elevation, sea):
    """Elevation colours times the hillshade, on the kept window (cropped already)."""
    coslat = np.maximum(crop(domain, domain.coslat[:, None] * np.ones((1, domain.columns)))[:, 0], 0.25)  # no black poles
    dy = domain.fine_km * 1000.0
    shade = shade_factor(hillshade(elevation, dy, dy * coslat[:elevation.shape[0]], domain.exaggeration))
    return np.clip(terrain_base(elevation, sea) * shade[..., None], 0, 255)


HARDNESS_STOPS = [(-1.0, (196, 120, 60)), (0.0, (232, 226, 205)), (1.0, (60, 120, 170))]


def render_stages(domain, fields, final_elevation, config, directory):
    """One map per stage of the procedural pass and the physics, each the surface so far."""
    sea = crop(domain, fields['sea'])
    stages = fields['stages']
    images = []
    for index, (key, caption) in enumerate((('base', 'base: B-spline of the tectonic world'),
                                            ('warp', 'warp: a grid carried by the displacement'),
                                            ('relief', 'relief noise'), ('uplands', 'uplands lift'),
                                            ('hardness', 'rock hardness: log erodibility (brown hard, blue soft)'),
                                            ('filter', 'erosion filter'), ('physics', 'drainage and erosion')), 1):
        if key == 'warp':
            rgb = relief_rgb(domain, crop(domain, stages['base']), sea)
            spacing = config['render']['warp_grid_km']
            width = 0.75 * domain.fine_km / spacing
            xw, yw = (crop(domain, values) / spacing for values in stages['warp'])
            lines = (np.abs(xw - np.round(xw)) < width) | (np.abs(yw - np.round(yw)) < width)
            rgb[lines] = (255, 255, 255)
        elif key == 'hardness':
            values = np.log(crop(domain, stages['hardness'])) / config['hardness']['sigma']
            rgb = colourize(values, HARDNESS_STOPS)
        elif key == 'physics':
            if final_elevation is None:
                continue
            rgb = relief_rgb(domain, crop(domain, final_elevation), sea)
        else:
            rgb = relief_rgb(domain, crop(domain, stages[key]), sea)
        name = f'{domain.name}-stage{index}-{key}.png'
        save_rgb(directory / name, rgb)
        images.append((f'{index}. {caption}', name))
    return images


def render(domain, elevation, sea, discharge, directory, suffix=''):
    elevation, sea = crop(domain, elevation), crop(domain, sea)
    base_rgb = terrain_base(elevation, sea)
    names = {}
    names['elevation'] = f'{domain.name}{suffix}-elevation.png'
    save_rgb(directory / names['elevation'], base_rgb)
    names['relief'] = f'{domain.name}{suffix}-relief.png'
    save_rgb(directory / names['relief'], relief_rgb(domain, elevation, sea))
    if discharge is not None:
        discharge = crop(domain, discharge)
        logarithm = np.where(~sea & (discharge > 0), np.log10(np.maximum(discharge, 1e-9)), np.nan)
        # The palette's 0..6 stretched over this map's median to maximum, so a window's rivers are readable.
        low, high = np.nanpercentile(logarithm, 50), np.nanmax(logarithm)
        logarithm = (logarithm - low) / max(high - low, 1e-6) * 6.0
        rgb = palette_rgb(logarithm, 'tectonics:discharge_log').astype(float)
        rgb[sea] = (18, 22, 30)
        names['discharge'] = f'{domain.name}{suffix}-discharge.png'
        save_rgb(directory / names['discharge'], rgb)
    return names


PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Terrain sandbox {label}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}h2{{margin:24px 0 8px}}
.row{{display:flex;gap:8px;flex-wrap:wrap}}figure{{margin:0;flex:1 1 30%;min-width:280px}}
img{{width:100%;image-rendering:pixelated;cursor:zoom-in}}figcaption{{color:#999;margin:4px 0}}
pre{{background:#0d0f12;padding:8px;overflow:auto;font-size:12px}}</style></head><body>
<h1>Terrain sandbox: {label}</h1><p>{note}</p>{sections}<h2>Configuration</h2><pre>{config}</pre>
<h2>Timings and physics</h2><pre>{timings}</pre></body></html>"""


def write_page(directory, label, sections, config, timings):
    html = []
    for title, images in sections:
        figures = ''.join(f'<figure><a href="{path}"><img src="{path}"></a><figcaption>{caption}</figcaption></figure>'
                          for caption, path in images)
        html.append(f'<h2>{title}</h2><div class="row">{figures}</div>')
    note = ('Seed 1; base: procedural-world at 256 columns after 600 Myr of tectonics and the terrain stages, '
            'before the climate and the global erosion. Uniform runoff, no climate. Windows are zooms of the same '
            'noise functions; their physics runs alone (no rivers pinned from the world), on the window plus a '
            'margin with a sea ring, cropped.')
    (directory / 'index.html').write_text(PAGE.format(label=label, note=note, sections=''.join(html),
                                                      config=json.dumps(config, indent=1),
                                                      timings=json.dumps(timings, indent=1)), encoding='utf-8')


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--label', required=True)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    parser.add_argument('--set', action='append', default=[], metavar='GROUP.KEY=VALUE')
    parser.add_argument('--domains', default=None, help="comma-separated, default the config's default_domains")
    parser.add_argument('--no-physics', action='store_true')
    arguments = parser.parse_args()

    config = load_config(arguments.config, arguments.set)
    directory = ROOT / 'out' / 'scratch' / 'sandbox' / arguments.label
    directory.mkdir(parents=True, exist_ok=True)
    bundle = Bundle.open(ROOT / config['base_bundle'])
    world = world_fields(bundle)
    names = arguments.domains.split(',') if arguments.domains else config['default_domains']
    physics_on = config['physics']['enabled'] and not arguments.no_physics
    sections, timings = [], {}
    for domain in make_domains(config, bundle, names):
        started = time.perf_counter()
        fields, timing = procedural(domain, world, config, config['seed'])
        images = []
        if config['render'].get('before_physics', False) or not physics_on:
            before = render(domain, fields['surface'], fields['sea'], None, directory, '-procedural')
            images += [('relief before physics', before['relief'])]
        if physics_on:
            physics_started = time.perf_counter()
            elevation, discharge, report = run_physics(domain, fields, config, directory)
            timing['physics'] = time.perf_counter() - physics_started
            timing['physics_report'] = report
            # The sea keeps its procedural depth: the marine deposition's straight walks are not what is judged.
            elevation = np.where(fields['sea'], fields['surface'], elevation)
            after = render(domain, elevation, fields['sea'], discharge, directory)
            images = [('elevation', after['elevation']), ('relief (hillshade, NW light)', after['relief']),
                      ('discharge (log, stretched from the median to the maximum)', after['discharge'])] + images
        else:
            images = [('elevation (procedural, no physics)', before['elevation'])] + images
        timing['total'] = time.perf_counter() - started
        timings[domain.name] = {key: (round(value, 2) if isinstance(value, float) else value) for key, value in timing.items()}
        title = 'Whole world' if domain.whole else f'Window {domain.name}'
        sections.append((f'{title} ({domain.columns if domain.whole else domain.crop[3] - domain.crop[2]} columns, '
                         f'{domain.fine_km:.2f} km cells, k {domain.k})', images))
        if domain.whole and config['render'].get('stages', True):
            stage_images = render_stages(domain, fields, elevation if physics_on else None, config, directory)
            sections.append(('Whole world, stage by stage (each map the surface so far)', stage_images))
        print(f'{domain.name}: {timings[domain.name]}', flush=True)
    (directory / 'config.json').write_text(json.dumps(config, indent=1), encoding='utf-8')
    write_page(directory, arguments.label, sections, config, timings)
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
