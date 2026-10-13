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

Config groups: domains, seafloor, warp, structure, zones, relief, coast, micro, uplands, hardness, filter,
physics, sea_level, glacial, precipitation, render. sea_level: a history of [Myr, m] points (piecewise
linear, a repeated time a jump) the physics runs through (--sea-history), the default the main erosion at
0 m, a lowstand at -120 m for 3 Myr, then the rise to 0; when enabled the physics' duration is the
history's last time (global and windows alike). At the lowstand the shelf emerges and the rivers and
glaciers cut down to the lowstand; the rise drowns their lower valleys (estuaries, rias, longer fjords).
keep_shelf: the sea shallower than the lowest level keeps the physics' surface instead of its procedural
depth, so the valleys cut across the exposed shelf stay drawn.

  uv run python scripts/sandbox/terrain_sandbox.py --label try1 [--config FILE] [--set group.key=value ...]
      [--domains global,A,B] [--no-physics]
Writes out/scratch/sandbox/<label>/: index.html, the PNGs, config.json (the values used) and timings.
"""
import argparse
import copy
import json
import math
import shutil
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
import overlays  # noqa: E402

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
            'uplift': uplift, 'sediment': bundle.load('geosphere.sediment').astype(float),
            'continental_crust': bundle.load('geosphere.continental_crust').astype(float),
            'crust_class': bundle.load('geosphere.crust_class'), 'crust_age': bundle.load('geosphere.crust_age').astype(float),
            'boundary_type': bundle.load('geosphere.boundary_type')}


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
    rain = np.maximum(sample(world['precipitation'], domain, mx, my), 0.0)
    sea_seed = sample(world['sea'].astype(float), domain, mx, my) > 0.5
    plain_base = base
    ocean_age, spreading = None, None
    if config['seafloor'].get('shelf', {}).get('enabled', False):
        base = shelf_profile(base, np.clip(sample(world['ocean']['oceanic'], domain, mx, my), 0.0, 1.0),
                             config['seafloor']['shelf'])
    if config['seafloor']['enabled']:
        base, ocean_age, spreading = seafloor(domain, world, world['ocean'], base, mx, my, config['seafloor'], seed)
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
    x0, y0 = x, y  # the plain coordinates, before the warp
    x, y = warped(x, y, 0.5 * np.angle(warp_doubled), warp_alignment, config['warp'], seed, domain)
    ratio, long_axis = warp_stretch(domain, x, y)
    stages = {'base': plain_base, 'seafloor': base, 'warp': (x, y), 'stretch': (ratio, long_axis),
              'strike': (theta, alignment), 'ocean_age': ocean_age, 'spreading': spreading}
    # The belt's zones (multiscale plan S5.5): the foreland (moderate strain) is folded and layered, the core
    # (high strain or thick crust) is crystalline: homogeneous rock and a relief stretched less along the
    # strike, so its rivers can branch (dendritic) instead of following bands.
    zones = config['zones']
    crust = sample(world['continental_crust'], domain, mx, my)
    foreland = smoothstep(strain, *zones['foreland_strain']) * (1 - smoothstep(strain, *zones['core_strain']))
    core = np.maximum(smoothstep(strain, *zones['core_strain']),
                      smoothstep(crust, *zones['core_crust_m']) * smoothstep(strain, *zones['foreland_strain']))
    foreland = foreland * (1 - core)
    stages['zones'] = (foreland, core)
    relief = config['relief']
    relief_alignment = alignment * (1 - zones['core_isotropy'] * core)
    noise = oriented_fbm(x, y, theta, relief_alignment, {**relief, **structure}, seed + 20, domain)
    noise = (1 - relief['ridged']) * noise + relief['ridged'] * ridged(noise) / 0.6
    amplitude = relief['base_m'] + relief['elevation_fraction'] * np.maximum(base, 0.0) + \
        relief['strain_m'] * alignment
    # The sea's smaller relief fades in with depth (from +200 m down to -500 m): a switch at the base's zero
    # contour jumped the noise along the smooth old coastline, a crease in the sea floor and the coast.
    wet = 1.0 - smoothstep(base, -500.0, 200.0)
    amplitude = amplitude * (1.0 - (1.0 - relief['sea_fraction']) * wet)
    surface = base + amplitude * noise
    stages['relief'] = surface
    uplands = config['uplands']
    if uplands['lift_m'] > 0:  # continental interiors raised in patches: a non-negative lift, none at the coast
        land_share = sample(ndimage.gaussian_filter((~world['sea']).astype(float), uplands['coast_sigma_cells'],
                                                    mode=('nearest', 'wrap')), domain, mx, my)
        inland = np.clip((land_share - 0.5) * 2.0, 0.0, 1.0)
        large = fbm(x, y, uplands, seed + 50, domain)
        surface = surface + uplands['lift_m'] * inland * np.clip(0.5 + 0.5 * large, 0.0, 1.5) * \
            smoothstep(base, 0.0, 300.0)  # faded in above the coast, no step at the base's zero contour
    stages['uplands'] = surface
    erodibility = hardness(domain, world, config['hardness'], structure, x, y, theta, alignment * foreland, mx, my, seed)
    stages['hardness'] = erodibility
    if config['coast']['amplitude_m'] > 0:
        # coast.warp_fraction: the share of the warp the coast term's noise follows (1: all, v14; less: its
        # features keep the trend but vary more along a strike-parallel shore).
        fraction = config['coast'].get('warp_fraction', 1.0)
        xc, yc = (x, y) if fraction == 1.0 else (x0 + fraction * (x - x0), y0 + fraction * (y - y0))
        surface = surface + coast_detail(domain, base, surface, xc, yc, amplitude, erodibility, config["coast"], seed)
    micro = config.get('micro', {'enabled': False})
    if micro['enabled'] and micro.get('when', 'after') == 'before':
        surface = surface + micro_texture(domain, surface, sea_seed, x, y, x0, y0, ratio, micro, seed)
    stages['coast'] = surface
    stages['precipitation'] = rain
    timings['noise'] = time.perf_counter() - started

    started = time.perf_counter()
    filtered = surface
    spec = config['filter']
    if spec['enabled'] and spec.get('when', 'before') == 'before':
        filtered = surface + erosion_pass(domain, surface, x, y, amplitude, erodibility, spec, seed)
    stages['filter'] = filtered
    timings['filter'] = time.perf_counter() - started

    land_or_low = filtered <= 0.0
    labels, _ = ndimage.label(land_or_low, structure=np.ones((3, 3)))
    sea_labels = np.unique(labels[(labels > 0) & sea_seed])
    sea = np.isin(labels, sea_labels[sea_labels > 0])
    pinned = np.zeros_like(sea)
    if not domain.whole:
        # The outer ring drains the window (sea at or below -1 m); the ring inside it is pinned at the lowest
        # terrain within half a world cell (a valley floor's height), so a river leaving the window keeps the
        # land's height as its base level. With the sea ring alone, every large river crossing the edge was
        # graded to the sea level and cut a canyon down to it (the user's "rivers jumping up", v10 to v13).
        ring = np.zeros_like(sea)
        ring[0, :] = ring[-1, :] = ring[:, 0] = ring[:, -1] = True
        pinned[1:-1, 1:-1] = True
        pinned[2:-2, 2:-2] = False
        pinned &= ~sea
        size = max(3, domain.k // 2) | 1
        floor = ndimage.minimum_filter(filtered, size=size, mode='nearest')
        filtered = np.where(pinned, np.maximum(floor, 1.0), filtered)
        sea |= ring
        filtered = np.where(ring, np.minimum(filtered, -1.0), filtered)
    return {'base': base, 'surface': filtered, 'sea': sea, 'erodibility': erodibility, 'uplift': uplift,
            'strain': strain, 'stages': stages, 'precipitation': rain, 'coords': (x, y), 'amplitude': amplitude,
            'pinned': pinned, 'warped': (x, y), 'plain': (x0, y0), 'ratio': ratio}, timings


def warp_stretch(domain, x, y):
    """How the warp distorts the noise drawn on its coordinates: the singular values s1 >= s2 of its Jacobian
    d(x', y') / d(x, y) (the planar kilometres of the grid against the warped ones). A noise feature of size L
    appears L / s in the map, so the anisotropy is s1 / s2 and features are long along the right singular
    vector of s2 (the map direction in which the warped coordinates change slowest). Returns (s1 / s2, the
    long axis's angle from east toward south, radians)."""
    dx_c, dx_r = np.gradient(x, axis=1) / domain.fine_km, np.gradient(x, axis=0) / domain.fine_km
    dy_c, dy_r = np.gradient(y, axis=1) / domain.fine_km, np.gradient(y, axis=0) / domain.fine_km
    p, q, r = dx_c ** 2 + dy_c ** 2, dx_r ** 2 + dy_r ** 2, dx_c * dx_r + dy_c * dy_r  # J^T J = [[p, r], [r, q]]
    root = np.sqrt(((p - q) / 2.0) ** 2 + r ** 2)
    largest, smallest = (p + q) / 2.0 + root, np.maximum((p + q) / 2.0 - root, 1e-12)
    ratio = np.sqrt(largest / smallest)
    long_axis = 0.5 * np.arctan2(2.0 * r, p - q) + math.pi / 2.0  # perpendicular to the fastest change
    return ratio, long_axis


def micro_texture(domain, surface, sea_seed, x, y, x0, y0, ratio, spec, seed, open_sea=None):
    """A little fine texture on the shoreline where the warp stretches most (option 1 of the user's
    2026-10-12 talk). The coast term's displacement is drawn on the fully warped coordinates, where an
    anisotropy of 4 makes it vary only every 150 km along a strike-parallel shore: it moves a whole arc
    together. Here the surface within reach_km of the current shoreline (the sea connected to the base's
    sea) is resampled at positions moved by up to displacement_km, a fine fBm (wavelength_km) drawn at
    warp_fraction of the warp (features elongated along the trend by about the stretch^warp_fraction, not
    erased), plus a small vertical part (amplitude_m); both weighted by smoothstep(stretch ratio, ratio)."""
    fraction = spec['warp_fraction']
    xm, ym = x0 + fraction * (x - x0), y0 + fraction * (y - y0)
    if open_sea is None:
        labels, _ = ndimage.label(surface <= 0.0, structure=np.ones((3, 3)))
        touching = np.unique(labels[(labels > 0) & sea_seed])
        sea = np.isin(labels, touching[touching > 0])
    else:  # after the physics: the open sea's shore only, not the glacial troughs it flooded inland (G-coast)
        sea = open_sea
    distance = np.where(sea, ndimage.distance_transform_edt(sea), ndimage.distance_transform_edt(~sea))
    near = np.exp(-(distance * domain.fine_km / spec['reach_km']) ** 2)
    weight = smoothstep(ratio, *spec['ratio']) * near
    scale = spec['displacement_km'] / domain.fine_km * weight
    dx = scale * fbm(xm, ym, spec, seed + 75, domain)
    dy = scale * fbm(xm, ym, spec, seed + 76, domain)
    rows, columns = np.indices(surface.shape, dtype=float)
    mode = 'grid-wrap' if domain.whole else 'nearest'
    moved = ndimage.map_coordinates(surface, [rows + dy, columns + dx], order=1, mode=mode)
    vertical = spec['amplitude_m'] * weight * fbm(xm, ym, spec, seed + 77, domain)
    return moved + vertical - surface


def coast_detail(domain, base, surface, x, y, relief_m, erodibility, spec, seed):
    """Relief concentrated around the sea level, so the shoreline is a fractal line rather than a contour of
    smooth fields: a rough fBm (gain near 0.6 keeps the fine octaves strong, a ridged share makes headlands
    and skerries) times exp(-(surface / band_m)^2) on the surface so far. Its amplitude follows the terrain:
    hard rock (erodibility^hardness_exponent, negative: hard coasts are rocky, soft ones smooth), the local
    relief (1 + relief_gain x relief / relief_reference_m: mountains meeting the sea) and cold coasts (1 +
    cold_gain x smoothstep(|latitude|, cold_latitude): the skerry and fjord coasts of former ice). Pure
    function of position, so every zoom draws the same coast down to its grid's resolution."""
    noise = fbm(x, y, spec, seed + 70, domain)
    noise = (1 - spec['ridged']) * noise + spec['ridged'] * ridged(noise) / 0.6
    latitude = np.abs(domain.north_deg - (np.arange(domain.rows) + 0.5) * domain.pixel_deg)[:, None]
    strength = spec['amplitude_m'] * erodibility ** spec['hardness_exponent'] * \
        (1.0 + spec['relief_gain'] * np.clip(relief_m / spec['relief_reference_m'], 0.0, 2.0)) * \
        (1.0 + spec['cold_gain'] * smoothstep(latitude, *spec['cold_latitude']))
    # Near the coast only: the band by height alone also took every low inland valley (G-coast: valley floors
    # pushed below the sea and flooded 100 km inland); the distance to the base's coastline (the zero contour
    # of the smooth base, both sides) fades it beyond reach_km.
    land = base > 0.0
    distance = np.where(land, ndimage.distance_transform_edt(land), ndimage.distance_transform_edt(~land))
    near = np.exp(-(distance * domain.fine_km / spec['reach_km']) ** 2)
    strength = strength * near
    band = np.exp(-(surface / spec['band_m']) ** 2)
    detail = strength * band * noise
    if spec.get('warp_km', 0) > 0:
        # A vertical noise moves a steep shoreline by only its amplitude over the slope; a horizontal fractal
        # displacement of the surface near the coast (twice the band) bends it whatever the slope. The surface
        # is resampled at the displaced positions (bilinear), the displacement scaled as the vertical term.
        warp_spec = {**spec, 'wavelength_km': spec['warp_wavelength_km'], 'octaves': spec['warp_octaves']}
        wide = np.exp(-(surface / (2.0 * spec['band_m'])) ** 2)
        scale = spec['warp_km'] / domain.fine_km * wide * strength / spec['amplitude_m']
        dx = scale * fbm(x, y, warp_spec, seed + 71, domain)
        dy = scale * fbm(x, y, warp_spec, seed + 72, domain)
        rows, columns = np.indices(surface.shape, dtype=float)
        mode = 'grid-wrap' if domain.whole else 'nearest'
        moved = ndimage.map_coordinates(surface + detail, [rows + dy, columns + dx], order=1, mode=mode)
        return moved - surface
    return detail


def hardness_domains(world, spec):
    """ln relative erodibility of the rock domains on the world grid, from the tectonic fields:
    sedimentary basins soft (+sediment_soft x tanh(sediment / sediment_scale_m)); shields hard (crust of
    normal thickness, shield_crust_m, without sediment or shortening: -shield_hard); orogen cores hard
    (thick crust above core_crust_m under shortening: -core_hard); thinned crust soft (below thin_crust_m:
    +thin_soft). Smoothed over one world cell so the domains have no cell edges."""
    sediment = np.tanh(world['sediment'] / spec['sediment_scale_m'])
    crust = world['continental_crust']
    strained = np.clip(world['shortening'] / spec['strain_full'], 0.0, 1.0)
    low, high = spec['shield_crust_m']
    normal = np.clip((crust - low) / (high - low), 0.0, 1.0) * np.clip((high + 10000.0 - crust) / 10000.0, 0.0, 1.0)
    shield = normal * (1 - sediment) * (1 - strained)
    core = np.clip((crust - spec['core_crust_m']) / 10000.0, 0.0, 1.0) * strained
    thin = np.clip((spec['thin_crust_m'] - crust) / 8000.0, 0.0, 1.0) * (crust > 0)
    # Basins the engine leaves without sediment (v9 against GLiM: Earth's soft rock lies in interior basins):
    # low continental interior away from the belts, and the lowland in front of a belt (a foreland basin).
    elevation = world['elevation']
    share = ndimage.gaussian_filter((~world['sea']).astype(float), 2.0, mode=('nearest', 'wrap'))
    inland = np.clip((share - 0.5) * 2.0, 0.0, 1.0)
    interior = np.clip((spec['basin_elevation_m'] - elevation) / spec['basin_elevation_m'], 0.0, 1.0) * inland * \
        (1 - strained)
    near_belt = ndimage.gaussian_filter(strained, spec['foreland_sigma_cells'], mode=('nearest', 'wrap'))
    foreland = np.clip(2.0 * (near_belt - strained), 0.0, 1.0) * np.clip((1500.0 - elevation) / 1000.0, 0.0, 1.0) * \
        (elevation > 0)
    value = spec['sediment_soft'] * sediment - spec['shield_hard'] * shield - spec['core_hard'] * core + \
        spec['thin_soft'] * thin + spec['basin_soft'] * interior + spec['foreland_soft'] * foreland
    return ndimage.gaussian_filter(value, 1.0, mode=('nearest', 'wrap'))


def hardness(domain, world, spec, structure, x, y, theta, banded, mx, my, seed):
    """Relative erodibility K / K0: exp of the rock domains (hardness_domains), a weak large-scale noise, a
    mosaic at the scale of map units (about 60 km) and, where `banded` (the belts' foreland alignment), hard
    and soft bands across the strike (fold-scale layers stretched along it, the material of trellis
    drainage), their strength varied along the belt. Bands finer than a grid shows fade out like every octave."""
    value = sample(hardness_domains(world, spec), domain, mx, my)
    value += spec['noise']['sigma'] * fbm(x, y, spec['noise'], seed + 30, domain)
    value += spec['mosaic']['sigma'] * fbm(x, y, spec['mosaic'], seed + 32, domain)
    bands = spec['bands']
    layered = oriented_fbm(x, y, theta, np.ones_like(banded), {**bands, 'orientations': structure['orientations'],
                                                               'sharpness': structure['sharpness']}, seed + 31, domain)
    contrast = bands['contrast']
    variation = np.clip(0.5 + 0.7 * fbm(x, y, bands['variation'], seed + 33, domain), 0.0, 1.0)
    value += bands['amplitude'] * banded * variation * np.tanh(contrast * layered) / math.tanh(contrast)
    return np.exp(value)


def precipitation(world, spec, radius_m):
    """Annual precipitation (m/yr) on the world grid, a proxy with no climate model: a zonal profile
    (spec['zonal']: [|latitude|, m/yr] pairs, Earth's zonal means roughly) carried by east-west winds (trade
    winds and polar easterlies blow west, the westerlies between westerly_band blow east). Along each row,
    downwind: the air's moisture m recharges over the sea (recharge per cell), rains out a little per 1000 km
    over land (loss_per_1000_km) and much more where the air rises (orographic_loss per km of rise), and the
    rain is zonal x m x (1 + orographic_gain x rise in km), less with altitude (exp(-elevation /
    altitude_scale_m)). So windward coasts and slopes are wet and the lee of a range is dry."""
    elevation = np.maximum(ndimage.gaussian_filter(world['elevation'], 1.0, mode=('nearest', 'wrap')), 0.0)
    sea = world['sea']
    rows, columns = elevation.shape
    latitude = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    zonal = np.interp(np.abs(latitude), *np.array(spec['zonal']).T)
    cell_km = 2.0 * math.pi * radius_m / 1000.0 / columns * np.maximum(np.cos(np.radians(latitude)), 0.05)
    low, high = spec['westerly_band']
    rain = np.zeros_like(elevation)
    for row in range(rows):
        eastward = low <= abs(latitude[row]) < high
        order = np.arange(columns) if eastward else np.arange(columns)[::-1]
        heights, wet = elevation[row, order], sea[row, order]
        moisture = 1.0
        for sweep in range(2):  # the second lap starts from the first lap's moisture: periodic rows
            for i in range(columns):
                rise = max(0.0, heights[i] - heights[i - 1]) / 1000.0
                if wet[i]:
                    moisture += (1.0 - moisture) * spec['recharge']
                    value = zonal[row] * moisture
                else:
                    value = zonal[row] * moisture * (1.0 + spec['orographic_gain'] * rise) * \
                        math.exp(-heights[i] / spec['altitude_scale_m'])
                    loss = spec['loss_per_1000_km'] * cell_km[row] / 1000.0 + spec['orographic_loss'] * rise
                    moisture *= math.exp(-loss)
                if sweep == 1:
                    rain[row, order[i]] = value
    return ndimage.gaussian_filter(rain, (2.5, 1.0), mode=('nearest', 'wrap'))  # across rows too: no row streaks


def smoothstep(values, low, high):
    t = np.clip((values - low) / (high - low), 0.0, 1.0)
    return t * t * (3 - 2 * t)


def ocean_age_depth(age_myr):
    """Sea-floor depth (m, positive down) of oceanic crust from its age: GDH1 (Stein and Stein 1992),
    2600 + 365 sqrt(t) up to 20 Myr, 5651 - 2473 exp(-0.0278 t) after."""
    t = np.maximum(age_myr, 0.0)
    return np.where(t < 20.0, 2600.0 + 365.0 * np.sqrt(t), 5651.0 - 2473.0 * np.exp(-0.0278 * t))


def ocean_world(world, spec, seed, radius_km):
    """World-grid fields of the sea floor: the oceanic share (crust class 1), the spreading axis (the
    direction age increases in, as doubled-angle components weighted by the age gradient, smoothed) and the
    trench indicator (convergent boundary cells beside oceanic crust, smoothed). With half_rate_km_myr, also
    the distance to the spreading axis on a grid ridge_grid_factor times finer than the world's (sea-v4):
    the engine's crust ages about 200 Myr within 450 km of a ridge (plates at 0.7 cm/yr, a tenth of
    Earth's), so its ridges are narrow; the sea floor reads the age as this distance over an Earth-like
    half rate instead, wherever that is younger. The axis is where the age read through the segments
    (segment_shift, so the crest steps at the transforms) is a minimum along the spreading direction (no
    higher than axis_step_cells ahead and behind), on crust younger than young_myr; the distance is planar
    (the sandbox's equirectangular km)."""
    oceanic = (world['crust_class'] == 1).astype(float)
    age = ndimage.gaussian_filter(world['crust_age'], 1.0, mode=('nearest', 'wrap'))
    gy, gx = np.gradient(age)
    strength = np.hypot(gx, gy) * oceanic
    doubled = strength * np.exp(2j * np.arctan2(gy, gx))
    doubled = ndimage.gaussian_filter(doubled.real, 1.5, mode=('nearest', 'wrap')) + \
        1j * ndimage.gaussian_filter(doubled.imag, 1.5, mode=('nearest', 'wrap'))
    near_ocean = ndimage.binary_dilation(oceanic > 0, structure=np.ones((3, 3)))
    trench = ((world['boundary_type'] == 1) & near_ocean).astype(float)
    trench = ndimage.gaussian_filter(trench, 1.0, mode=('nearest', 'wrap'))
    ocean = {'oceanic': oceanic, 'doubled': doubled, 'trench': trench,
             'strength': ndimage.gaussian_filter(strength, 1.5, mode=('nearest', 'wrap'))}
    if spec.get('half_rate_km_myr', 0.0) > 0.0:
        factor = spec['ridge_grid_factor']
        H, W = oceanic.shape
        grid = Domain('ridge-grid', W * factor, H * factor, 0.0, 0.0, factor, W, H, radius_km)
        mx, my = bspline_matrix(grid.x, W, True), bspline_matrix(grid.y, H, False)
        x, y = grid.mesh_km()
        share = np.clip(sample(oceanic, grid, mx, my), 0.0, 1.0)
        axis_doubled = sample(doubled.real, grid, mx, my) + 1j * sample(doubled.imag, grid, mx, my)
        spreading = 0.5 * np.angle(axis_doubled)
        xs, ys, _ = segment_shift(grid, x, y, spreading, share, spec, seed)
        grid_age = np.maximum(sample_at(world['crust_age'], grid, xs, ys), 0.0)
        # A minimum along the spreading direction: a minimum in every direction found only dots on the
        # engine's young lines (sea-v4), whose distances drew round swells instead of ridges.
        step = spec['axis_step_cells'] * grid.fine_km
        dx, dy = step * np.cos(spreading), step * np.sin(spreading)
        ahead = sample_at(world['crust_age'], grid, xs + dx, ys + dy)
        behind = sample_at(world['crust_age'], grid, xs - dx, ys - dy)
        axis = (share > 0.5) & (grid_age < spec['young_myr']) & (grid_age <= ahead) & (grid_age <= behind)
        tiled = np.concatenate([~axis] * 3, axis=1)  # columns wrap
        distance = ndimage.distance_transform_edt(tiled)[:, grid.columns:2 * grid.columns] * grid.fine_km
        ocean['ridge_distance'] = distance if axis.any() else np.full_like(distance, np.inf)
        ocean['ridge_factor'] = factor
    return ocean


def damp_abyss(world, spec):
    """The deep sea off oceanic crust keeps the engine's depths, which the B-spline draws as smooth blotches
    hundreds of km wide (sea-v3); Earth's abyssal plains are flat. Below start_m the depth's departure from
    its sea-only mean over sigma_cells world cells is kept at the share `keep`."""
    elevation, sea = world['elevation'], world['sea']
    depth = -elevation
    mode = ('nearest', 'wrap')
    weight = ndimage.gaussian_filter(sea.astype(float), spec['sigma_cells'], mode=mode)
    mean = ndimage.gaussian_filter(np.where(sea, depth, 0.0), spec['sigma_cells'], mode=mode) / np.maximum(weight, 1e-6)
    deep = smoothstep(depth, spec['start_m'], spec['start_m'] + 1000.0) * sea * (world['crust_class'] != 1)
    damped = mean + spec['keep'] * (depth - mean)
    return elevation - deep * (damped - depth)


def shelf_profile(base, oceanic, spec):
    """A shelf break on the sea off oceanic crust (sea-v4): the base's depth d (m) maps to shelf_m x d /
    break_m down to break_m, then steepens to the base's own depth over slope_m of it, so the shelf stays
    near shelf_m deep and the slope beyond its edge is several times steeper than the base's broad ramp."""
    depth = -base
    shelf, edge, slope = spec['shelf_m'], spec['break_m'], spec['slope_m']
    profile = np.where(depth <= edge, shelf * np.clip(depth / edge, 0.0, 1.0),
                       shelf + (depth - shelf) * smoothstep(depth, edge, edge + slope))
    if spec.get('foot_m'):
        # sea-v7b (the hypsometric curve: twice Earth's area between -2500 and 0 m): the slope runs from the
        # shelf to foot_m over the base's break_m to break_m + slope_m, linearly, and the depths beyond are
        # compressed into foot_m to bottom_m, so the mid-depth sea goes to the abyss as on Earth's curve.
        foot, bottom = spec['foot_m'], spec['bottom_m']
        end = edge + slope
        steep = shelf + (depth - edge) * (foot - shelf) / slope
        deep = foot + (depth - end) * (bottom - foot) / max(bottom - end, 1.0)
        profile = np.where(depth <= edge, profile, np.where(depth <= end, steep, np.maximum(deep, depth * 0 + foot)))
    weight = (base < 0.0) * (1.0 - oceanic)
    return base + weight * (-profile - base)


def sample_at(field, domain, x_km, y_km):
    """A world field at arbitrary planar positions (km, the domain's mesh units): cubic spline, columns
    wrapped, rows clamped."""
    pad = 4
    padded = np.concatenate([field[:, -pad:], field, field[:, :pad]], axis=1)
    columns = x_km / domain.cell_km - 0.5 + pad
    columns = np.mod(columns - pad, domain.world_columns) + pad
    rows = y_km / domain.cell_km - 0.5
    return ndimage.map_coordinates(padded, [rows, columns], order=3, mode='nearest')


def sample_fine(field, factor, domain, x_km, y_km):
    """A field on a grid `factor` times finer than the world's, at planar positions (km): linear (a
    distance's V keeps its crest), columns wrapped, rows clamped."""
    pad = 2
    padded = np.concatenate([field[:, -pad:], field, field[:, :pad]], axis=1)
    cell = domain.cell_km / factor
    columns = np.mod(x_km / cell - 0.5, field.shape[1]) + pad
    rows = y_km / cell - 0.5
    return ndimage.map_coordinates(padded, [rows, columns], order=1, mode='nearest')


def segment_shift(domain, x, y, theta, oceanic, spec, seed):
    """Ridge-transform segments: for a bank of spreading directions phi_b (weighted by cos^sharpness 2(theta -
    phi_b)), the along-ridge coordinate is cut into segments segment_km long, each shifted along the
    spreading direction by up to offset_km. Returns the shifted positions and the fracture-zone trough
    (0 to 1, scaled by its transform's offset). fracture_warp (sea-v6): the segments are cut on positions
    moved by a noise of amplitude_km, so the fracture zones curve and wander instead of running as straight,
    evenly spaced lines (the user: window O's lines looked artificial); fracture_patch: the trough comes and
    goes along its length, as Earth's fracture zones do."""
    shift_x, shift_y, trough, weights = (np.zeros_like(x) for _ in range(4))
    warp = spec.get('fracture_warp')
    if warp:
        xw = x + warp['amplitude_km'] * fbm(x, y, warp, seed + 130, domain)
        yw = y + warp['amplitude_km'] * fbm(x, y, warp, seed + 131, domain)
    else:
        xw, yw = x, y
    count, length = spec['orientations'], spec['segment_km']
    for b in range(count):
        phi = math.pi * b / count
        w = np.maximum(0.0, np.cos(2.0 * (theta - phi))) ** spec['sharpness']
        active = (oceanic > 0.01) & (w > 1e-3)
        if not np.any(active):
            continue
        xa, ya = xw[active], yw[active]
        along_ridge = (-math.sin(phi) * xa + math.cos(phi) * ya) / length
        jitter = spec.get('segment_jitter', 0.0)
        if jitter > 0.0:  # irregular spacing: the along-ridge coordinate stretched by a 1-D noise of itself
            along_ridge = along_ridge + jitter * gradient_noise(along_ridge * 0.6, np.full_like(along_ridge, 7.3 * b),
                                                                seed + 133)
        segment = np.floor(along_ridge).astype(np.int64)
        bank = np.full_like(segment, b)
        offset = spec['offset_km'] * erosion_filter.cell_hash(segment, bank, seed + 120)[0]
        shift_x[active] += w[active] * offset * math.cos(phi)
        shift_y[active] += w[active] * offset * math.sin(phi)
        # A fracture zone's trough scales with its transform's offset: the jump between the two segments it
        # separates (a comb of equal troughs at every boundary read as artificial, sea-v1).
        boundary = np.round(along_ridge).astype(np.int64)
        jump = np.abs(erosion_filter.cell_hash(boundary, bank, seed + 120)[0] -
                      erosion_filter.cell_hash(boundary - 1, bank, seed + 120)[0]) / 2.0
        edge = np.abs(along_ridge - boundary) * length  # km to the nearest segment boundary
        trough[active] += w[active] * jump * np.exp(-(edge / spec['fracture_width_km']) ** 2)
        weights[active] += w[active]
    weights = np.maximum(weights, 1e-12)
    trough = trough / weights
    if spec.get('fracture_patch'):
        trough = trough * np.clip(0.5 + 0.7 * fbm(x, y, spec['fracture_patch'], seed + 132, domain), 0.0, 1.5)
    return x + shift_x / weights, y + shift_y / weights, trough


def seamounts(x, y, spec, seed):
    """Scattered volcanic cones (sea-v6): one candidate per spacing_km lattice cell at a hashed position,
    kept with `probability`; size u in 0..1 sets the radius (radius_km range) and the height (height_m
    range, by u^2: many small, few large); profile (1 - r / R)^profile. The tallest cone wins where two
    overlap."""
    spacing = spec['spacing_km']
    ix, iy = np.floor(x / spacing).astype(np.int64), np.floor(y / spacing).astype(np.int64)
    out = np.zeros_like(x)
    (r0, r1), (h0, h1) = spec['radius_km'], spec['height_m']
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            cx, cy = ix + dx, iy + dy
            u, v = erosion_filter.cell_hash(cx, cy, seed + 140)
            keep, size = erosion_filter.cell_hash(cx, cy, seed + 141)
            size = 0.5 + 0.5 * size
            px, py = (cx + 0.5 + 0.45 * u) * spacing, (cy + 0.5 + 0.45 * v) * spacing
            radius = r0 + (r1 - r0) * size
            cone = (h0 + (h1 - h0) * size ** 2) * np.clip(1.0 - np.hypot(x - px, y - py) / radius, 0.0, 1.0) ** spec['profile']
            out = np.maximum(out, np.where(0.5 + 0.5 * keep < spec['probability'], cone, 0.0))
    return out


def seafloor(domain, world, ocean, base, mx, my, spec, seed):
    """The sea floor on oceanic crust (the user's 2026-10-13 request): depth from the crust age (GDH1),
    with the age read through ridge-transform segments (segment_shift: the ridge steps at each transform and
    the age step carries on across the plate as a fracture zone), an axial valley on the youngest crust,
    fracture-zone troughs, abyssal hills elongated along the ridge (fading with age as sediment buries them)
    and trenches on the ocean side of convergent boundaries. Blended into the base where the crust is
    oceanic and the base is deeper than shallow_m. Returns (the new base, the age read, the spreading axis)."""
    x, y = domain.mesh_km()
    oceanic = np.clip(sample(ocean['oceanic'], domain, mx, my), 0.0, 1.0)
    doubled = sample(ocean['doubled'].real, domain, mx, my) + 1j * sample(ocean['doubled'].imag, domain, mx, my)
    theta = 0.5 * np.angle(doubled)  # the spreading axis
    xs, ys, trough = segment_shift(domain, x, y, theta, oceanic, spec, seed)
    age = np.maximum(sample_at(world['crust_age'], domain, xs, ys), 0.0)
    if 'ridge_distance' in ocean:
        # sea-v4: the age as the distance to the spreading axis over an Earth-like half rate, where younger.
        distance = sample_fine(ocean['ridge_distance'], ocean['ridge_factor'], domain, x, y)
        age = np.minimum(age, distance / spec['half_rate_km_myr'])
    else:
        # The world's age is per 156 km cell, so near a ridge its minimum is a broad smooth trough that never
        # reaches 0: no crest, no axial valley, the transform steps lost in it (sea-v1). Young crust is aged
        # relative to the lowest age within ridge_radius_km, fully at the axis and not at all beyond young_myr.
        size = max(3, int(round(2.0 * spec['ridge_radius_km'] / domain.fine_km)) | 1)
        lowest = ndimage.minimum_filter(age, size=size, mode=('nearest', 'wrap') if domain.whole else 'nearest')
        young = 1.0 - smoothstep(age, 0.5 * spec['young_myr'], spec['young_myr'])
        age = np.maximum(age - lowest * young, 0.0)
    depth = ocean_age_depth(age)
    depth = depth + spec['axial_valley_m'] * np.exp(-(age / spec['axial_age_myr']) ** 2)
    depth = depth + spec['fracture_trough_m'] * trough * np.exp(-age / spec['fracture_age_myr'])
    # The hills are stretched along the ridge only where the spreading direction is coherent (sea-v6: where
    # the bank's orientation turned fast, the stretched noise drew swirls and fingerprints).
    coherence = np.ones_like(x)
    if spec.get('hills_coherence', False):
        coherence = np.clip(np.abs(doubled) / np.maximum(sample(ocean['strength'], domain, mx, my), 1e-9), 0.0, 1.0)
    hills = oriented_fbm(xs, ys, theta, coherence, spec['hills'], seed + 121, domain)
    amplitude = spec['hills']['amplitude_m'] * (spec['hills']['old_fraction'] + (1.0 - spec['hills']['old_fraction'])
                                                * np.exp(-age / spec['hills']['age_scale_myr']))
    if spec['hills'].get('patch'):
        amplitude = amplitude * np.clip(0.7 + 0.5 * fbm(x, y, spec['hills']['patch'], seed + 122, domain), 0.15, 1.6)
    floor = -depth + amplitude * hills
    if spec.get('seamounts', {}).get('enabled', False):
        floor = floor + seamounts(x, y, spec['seamounts'], seed)
    weight = oceanic * smoothstep(-base, spec['shallow_m'], 2.0 * spec['shallow_m']) * spec['blend']
    out = base + weight * (floor - base)
    # The trench: a V along the crest of the smoothed indicator (its value over its local maximum within a
    # world cell), so it is even along the boundary, not the whole smoothed band (sea-v1: beaded troughs
    # hundreds of km wide); trench_crest's lower end sets its width.
    trench = sample(ocean['trench'], domain, mx, my)
    size = max(3, int(round(domain.cell_km / domain.fine_km)) | 1)
    crest = trench / np.maximum(ndimage.maximum_filter(trench, size=size, mode=('nearest', 'wrap') if domain.whole
                                                       else 'nearest'), 1e-6)
    sharp = smoothstep(crest, *spec['trench_crest'])
    out = out - spec['trench_m'] * sharp * sharp * smoothstep(trench, *spec['trench_band']) * (base < 0)
    return out, age, theta


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


def filter_field(domain, fields, ice, spec, seed):
    """The filter's strength as a field (spec['field']): a slow noise (exp of sigma x fBm, regional variety),
    the rock (erodibility^hardness_exponent: soft rock carries denser gullies, badlands), the precipitation
    ((P / reference)^exponent clipped: wetter ground is more dissected) and no gullies under today's ice
    (1 - ice / ice_full_m: a glacier smooths its bed). 1 where every factor is neutral."""
    field_spec = spec['field']
    x, y = fields['coords']
    noise = field_spec['noise']
    field = np.exp(noise['sigma'] * fbm(x, y, noise, seed + 60, domain))
    field = field * fields['erodibility'] ** field_spec['hardness_exponent']
    low, high = field_spec['precipitation_clip']
    rain = np.clip(fields['precipitation'] / field_spec['precipitation_reference_m_yr'], low, high)
    field = field * rain ** field_spec['precipitation_exponent']
    if ice is not None:
        field = field * (1.0 - np.clip(ice / field_spec['ice_full_m'], 0.0, 1.0))
    return field


def erosion_pass(domain, height_m, x, y, amplitude_m, erodibility, spec, seed, field=None):
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
    if field is not None:
        strength = strength * field
    if 'slope_fade' in spec:  # gullies on slopes only: on flats they have no direction and crinkle the plain
        real_slope = np.sqrt(slope_x * slope_x + slope_y * slope_y) / vs
        strength = strength * smoothstep(real_slope, *spec['slope_fade'])
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
    runoff = np.maximum(fields['precipitation'] * config['precipitation']['runoff_fraction'],
                        config['precipitation']['runoff_min_m_yr'])
    for key, values in (('elevation', fields['surface']), ('sea', fields['sea'].astype(np.float32)),
                        ('erodibility', fields['erodibility']), ('uplift', fields['uplift'] * spec['uplift_factor']),
                        ('runoff', runoff), ('pinned', fields['pinned'].astype(np.float32))):
        files[key] = stem.with_name(f'{domain.name}-{key}.f32')
        np.asarray(values, dtype='<f4').tofile(files[key])
    out_elevation = stem.with_name(f'{domain.name}-out-elevation.f32')
    out_discharge = stem.with_name(f'{domain.name}-out-discharge.f32')
    sea_level = config.get('sea_level', {'enabled': False})
    duration = sea_level['history'][-1][0] if sea_level['enabled'] else spec['duration_myr']
    command = [str(PHYSICS), '--columns', str(domain.columns), '--rows', str(domain.rows),
               '--elevation', str(files['elevation']), '--sea', str(files['sea']),
               '--erodibility', str(files['erodibility']), '--uplift', str(files['uplift']),
               '--runoff-file', str(files['runoff']), '--duration', str(duration),
               '--step', str(spec['step_myr']), '--k', str(spec['k']), '--m', str(spec['m']),
               '--deposition', str(spec['deposition']), '--diffusivity', str(spec['diffusivity_km2_myr']),
               '--delta-top', str(spec['delta_top_m']), '--threshold-slope', str(spec['threshold_slope']),
               '--pinned', str(files['pinned']), '--radius-km', str(domain.radius_km),
               '--random-routing', '1' if spec.get('random_routing', False) else '0', '--seed', str(config['seed']),
               '--out-elevation', str(out_elevation), '--out-discharge', str(out_discharge)]
    if not domain.whole:
        command += ['--north', str(domain.north_deg), '--pixel-deg', str(domain.pixel_deg)]
    if sea_level['enabled']:
        command += ['--sea-history', ','.join(f'{time}:{level}' for time, level in sea_level['history'])]
    glacial = config['glacial']
    out_ice = stem.with_name(f'{domain.name}-out-ice.f32')
    ela = None
    if glacial['enabled']:
        ela = equilibrium_line(domain, fields['precipitation'], glacial)
        files['ela'] = stem.with_name(f'{domain.name}-ela.f32')
        np.asarray(ela, dtype='<f4').tofile(files['ela'])
        command += ['--ela', str(files['ela']), '--glacial-ratio', str(glacial['erodibility_ratio']),
                    '--ice-full-range', str(glacial['full_ice_range_m']),
                    '--ice-accumulation', str(glacial['ice_accumulation']), '--out-ice', str(out_ice)]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'terrain_sandbox failed on {domain.name}: {result.stderr.strip()}')
    shape = (domain.rows, domain.columns)
    elevation = np.fromfile(out_elevation, dtype='<f4').reshape(shape).astype(np.float64)
    discharge = np.fromfile(out_discharge, dtype='<f4').reshape(shape).astype(np.float64)
    ice = np.fromfile(out_ice, dtype='<f4').reshape(shape).astype(np.float64) if glacial['enabled'] else None
    for path in list(files.values()) + [out_elevation, out_discharge] + ([out_ice] if glacial['enabled'] else []):
        path.unlink()
    if ice is not None:  # the maps show today's glaciers: the ice where the surface is above today's line
        ice = np.where(elevation + ice > ela + glacial['present_shift_m'], ice, 0.0)
    return elevation, discharge, ice, result.stdout.strip()


def equilibrium_line(domain, precipitation, spec):
    """The glacial maximum's equilibrium-line altitude (m) per cell, from the latitude and the precipitation
    proxy (no climate here): ela_equator_m up to ela_flat_latitude (the tropics' line stays near 4.5 to 5 km),
    falling linearly to the sea level at ela_sea_latitude; lowered by ela_per_doubling_m for each doubling of
    the precipitation over ela_reference_m_yr (a wet coast's glaciers reach lower, Ohmura et al. 1992); never
    below ela_min_m, so polar lowlands carry no ice sheet that would grind them below the sea for the whole
    run (v10 and v11)."""
    latitude = np.abs(domain.north_deg - (np.arange(domain.rows) + 0.5) * domain.pixel_deg)
    flat, at_sea = spec['ela_flat_latitude'], spec['ela_sea_latitude']
    line = spec['ela_equator_m'] * np.clip((at_sea - latitude) / (at_sea - flat), 0.0, 1.0) ** spec['ela_power']
    wet = spec['ela_per_doubling_m'] * np.log2(np.maximum(precipitation, 0.05) / spec['ela_reference_m_yr'])
    return np.maximum(line[:, None] - wet, spec['ela_min_m'])


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
PRECIPITATION_STOPS = [(0.0, (170, 120, 70)), (0.3, (220, 200, 140)), (1.0, (90, 160, 80)), (3.0, (40, 90, 170))]
FIELD_STOPS = [(0.0, (10, 10, 12)), (1.0, (235, 235, 235)), (2.5, (240, 140, 40))]


def render_stages(domain, fields, final_elevation, config, directory):
    """One map per stage of the procedural pass and the physics, each the surface so far."""
    sea = crop(domain, fields['sea'])
    stages = fields['stages']
    images = []
    for index, (key, caption) in enumerate((('base', 'base: B-spline of the tectonic world'),
                                            ('seafloor', 'sea floor: depth from crust age, ridge-transform segments, abyssal hills, trenches'),
                                            ('warp', 'warp: a grid carried by the displacement'),
                                            ('relief', 'relief noise'), ('uplands', 'uplands lift'),
                                            ('hardness', 'rock hardness: ln erodibility over 0.7 (brown hard, blue soft)'),
                                            ('coast', 'coastal detail'),
                                            ('precipitation', 'precipitation proxy, m/yr (brown 0, green 1, blue 3 and more)'),
                                            ('filter', 'erosion filter'), ('physics', 'drainage and erosion')), 1):
        if key == 'warp':
            rgb = relief_rgb(domain, crop(domain, stages['base']), sea)
            spacing = config['render']['warp_grid_km']
            width = 0.75 * domain.fine_km / spacing
            xw, yw = (crop(domain, values) / spacing for values in stages['warp'])
            lines = (np.abs(xw - np.round(xw)) < width) | (np.abs(yw - np.round(yw)) < width)
            rgb[lines] = (255, 255, 255)
        elif key == 'hardness':
            values = np.log(crop(domain, stages['hardness'])) / 0.7  # +-0.7 in ln K spans the palette
            rgb = colourize(values, HARDNESS_STOPS)
            rgb[sea] = (18, 22, 30)  # the rock of the sea floor is not judged
        elif key == 'precipitation':
            rgb = colourize(crop(domain, stages['precipitation']), PRECIPITATION_STOPS)
            rgb[sea] = (18, 22, 30)
        elif key == 'physics':
            if final_elevation is None:
                continue
            rgb = relief_rgb(domain, crop(domain, final_elevation), sea)
        else:
            rgb = relief_rgb(domain, crop(domain, stages[key]), sea)
        name = f'{domain.name}-stage{index}-{key}.png'
        save_map(domain, directory / name, rgb)
        images.append((f'{index}. {caption}', name))
        earth = ROOT / config['render'].get('earth_hardness', '')
        if key == 'hardness' and earth.is_file():  # the Earth beside it, same palette and scale
            shutil.copyfile(earth, directory / 'earth-hardness.png')
            images.append(('5b. Earth: GLiM erodibility (Moosdorf et al. 2018 index; refs/processed/lithology)',
                           'earth-hardness.png'))
    return images


STRETCH_STOPS = [(1.0, (12, 16, 40)), (1.5, (40, 70, 140)), (2.5, (220, 120, 40)), (4.0, (255, 240, 200))]


def nice_spacing_km(width_km, lines=12):
    """A round grid spacing (1, 2, 2.5, 5 x 10^n km) giving about `lines` lines across the width."""
    target = width_km / lines
    power = 10.0 ** math.floor(math.log10(target))
    return min((step * power for step in (1, 2, 2.5, 5, 10)), key=lambda value: abs(value - target))


def glyphs(rgb, angle, length_px, every_px, colour=(255, 255, 255)):
    """Short line marks every `every_px` pixels along `angle` (radians from east toward south), each
    `length_px` long (an array, so a mark's length can show a strength). Returns uint8 RGB."""
    from PIL import Image, ImageDraw
    image = Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8))
    draw = ImageDraw.Draw(image)
    rows, columns = angle.shape
    for r in range(every_px // 2, rows, every_px):
        for c in range(every_px // 2, columns, every_px):
            half = 0.5 * float(length_px[r, c])
            if half < 0.5:
                continue
            dx, dy = half * math.cos(angle[r, c]), half * math.sin(angle[r, c])
            draw.line((c - dx, r - dy, c + dx, r + dy), fill=colour, width=2)
    return np.asarray(image)


def render_diagnostics(domain, fields, config, directory):
    """The fields that give the terrain its direction, for any domain: the warp (a grid carried by the
    displacement, spaced to the map), the warp's stretch (colour: the anisotropy s1 / s2 a noise drawn on the
    warped coordinates takes; marks along the direction its features are long, longer where it is stronger),
    the strike (hue: the strike direction, doubled so opposite directions match; brightness: the alignment
    that mixes in the stretched noise bank and the bands) and the rock hardness."""
    stages = fields['stages']
    sea = crop(domain, fields['sea'])
    images = []
    width_km = (domain.crop[3] - domain.crop[2] if domain.crop else domain.columns) * domain.fine_km
    spacing = config['render']['warp_grid_km'] if domain.whole else nice_spacing_km(width_km)
    rgb = relief_rgb(domain, crop(domain, stages['base']), sea)
    line = 0.75 * domain.fine_km / spacing
    xw, yw = (crop(domain, values) / spacing for values in stages['warp'])
    rgb[(np.abs(xw - np.round(xw)) < line) | (np.abs(yw - np.round(yw)) < line)] = (255, 255, 255)
    images.append((f'warp: a {spacing:g} km grid carried by the displacement, on the base relief', 'warp-grid'))
    maps = [rgb]
    ratio, axis = (crop(domain, values) for values in stages['stretch'])
    every = max(16, min(ratio.shape) // 32)
    strength = np.clip(np.log(ratio) / math.log(4.0), 0.0, 1.0) * every * 0.9
    maps.append(glyphs(colourize(ratio, STRETCH_STOPS), axis, strength, every))
    images.append(('warp stretch: anisotropy s1 / s2 (blue 1.5, orange 2.5, white 4 and more); marks along the '
                   'long axis of the features it draws', 'warp-stretch'))
    theta, alignment = (crop(domain, values) for values in stages['strike'])
    strike = theta + math.pi / 2.0
    hue = (np.mod(2.0 * strike, 2.0 * math.pi) / (2.0 * math.pi))[..., None]
    wheel = 0.5 + 0.5 * np.cos(2.0 * math.pi * (hue - np.array([0.0, 1.0 / 3.0, 2.0 / 3.0])))
    maps.append(glyphs(255.0 * wheel * (0.15 + 0.85 * np.clip(alignment, 0, 1))[..., None], strike,
                       np.clip(alignment, 0, 1) * every * 0.9, every))
    images.append(('strike: hue the strike direction (doubled), brightness the alignment that mixes in the '
                   'stretched noise and the bands; marks along the strike', 'strike'))
    hardness_rgb = colourize(np.log(crop(domain, stages['hardness'])) / 0.7, HARDNESS_STOPS)
    hardness_rgb[sea] = (18, 22, 30)
    maps.append(hardness_rgb)
    images.append(('rock hardness: ln erodibility over 0.7 (brown hard, blue soft)', 'hardness'))
    out = []
    for (caption, key), rgb in zip(images, maps):
        name = f'{domain.name}-diag-{key}.png'
        save_map(domain, directory / name, rgb)
        out.append((caption, name))
    return out


def save_map(domain, path, rgb):
    """A map with its scale references (overlays.py): on a window the faint grid of the world cells it was
    refined from (one line every k pixels; the crop starts on a world-cell edge) and on every map a scale
    bar measured east-west at the map's centre latitude."""
    if not domain.whole:
        rgb = overlays.grid(rgb, domain.k)
        latitude = domain.north_deg - domain.rows * domain.pixel_deg / 2.0
        note = f'(east-west, {abs(latitude):.0f} {"N" if latitude >= 0 else "S"}; grid: world cells)'
        km_per_px = domain.fine_km * math.cos(math.radians(latitude))
    else:
        note, km_per_px = '(east-west, at the equator)', domain.fine_km
    save_rgb(path, overlays.scale_bar(rgb, km_per_px, note))


def with_ice(rgb, ice):
    """Glacier ice thicker than 100 m as a pale blue-white wash on the land colours (thin valley tongues
    left out: they drew every channel as a white thread)."""
    if ice is None:
        return rgb
    alpha = np.clip((ice - 100.0) / 300.0, 0.0, 0.6)[..., None]
    return rgb * (1 - alpha) + np.array([236.0, 242.0, 248.0]) * alpha


def render(domain, elevation, sea, discharge, directory, suffix='', ice=None):
    elevation, sea = crop(domain, elevation), crop(domain, sea)
    ice = None if ice is None else crop(domain, ice)
    base_rgb = with_ice(terrain_base(elevation, sea), ice)
    names = {}
    names['elevation'] = f'{domain.name}{suffix}-elevation.png'
    save_map(domain, directory / names['elevation'], base_rgb)
    names['relief'] = f'{domain.name}{suffix}-relief.png'
    coslat = np.maximum(crop(domain, domain.coslat[:, None] * np.ones((1, domain.columns)))[:, 0], 0.25)
    dy = domain.fine_km * 1000.0
    shade = shade_factor(hillshade(elevation, dy, dy * coslat[:elevation.shape[0]], domain.exaggeration))
    save_map(domain, directory / names['relief'], np.clip(base_rgb * shade[..., None], 0, 255))  # the ice shaded too
    if discharge is not None:
        discharge = crop(domain, discharge)
        logarithm = np.where(~sea & (discharge > 0), np.log10(np.maximum(discharge, 1e-9)), np.nan)
        # The palette's 0..6 stretched over this map's median to maximum, so a window's rivers are readable.
        low, high = np.nanpercentile(logarithm, 50), np.nanmax(logarithm)
        logarithm = (logarithm - low) / max(high - low, 1e-6) * 6.0
        rgb = palette_rgb(logarithm, 'tectonics:discharge_log').astype(float)
        rgb[sea] = (18, 22, 30)
        names['discharge'] = f'{domain.name}{suffix}-discharge.png'
        save_map(domain, directory / names['discharge'], rgb)
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
    radius_m = float(bundle.manifest.get('parameters', {}).get('world.body.radius_km', 6371.0)) * 1000.0
    world['precipitation'] = precipitation(world, config['precipitation'], radius_m)
    if config['seafloor'].get('abyss', {}).get('enabled', False):
        world['elevation'] = damp_abyss(world, config['seafloor']['abyss'])
    world['ocean'] = ocean_world(world, config['seafloor'], config['seed'], radius_m / 1000.0)
    names = arguments.domains.split(',') if arguments.domains else config['default_domains']
    physics_on = config['physics']['enabled'] and not arguments.no_physics
    sections, timings = [], {}
    for domain in make_domains(config, bundle, names):
        started = time.perf_counter()
        fields, timing = procedural(domain, world, config, config['seed'])
        images = []
        if config['render'].get('before_physics', False) or not physics_on:
            before = render(domain, fields['surface'], fields['sea'], None, directory, '-procedural')
            if not physics_on:  # for scripts/sandbox/hypsometry.py and other measurements
                np.save(directory / f'{domain.name}-elevation.npy', crop(domain, fields['surface']).astype(np.float32))
            images += [('relief before physics', before['relief'])]
        if physics_on:
            physics_started = time.perf_counter()
            elevation, discharge, ice, report = run_physics(domain, fields, config, directory)
            timing['physics'] = time.perf_counter() - physics_started
            timing['physics_report'] = report
            # The sea keeps its procedural depth: the marine deposition's straight walks are not what is judged.
            # Land cut below the sea and connected to it (a fjord, a drowned river mouth) becomes sea.
            procedural_sea = fields['sea'].copy()  # the open sea before the physics flooded its troughs
            sea_level = config.get('sea_level', {'enabled': False})
            restored = fields['sea']
            if sea_level['enabled'] and sea_level.get('keep_shelf', False):
                # The shelf the lowstand exposed keeps the physics' surface (the valleys cut across it, drowned
                # by the rise); only the sea deeper than the lowest level gets its procedural depth back. A shelf
                # cell the physics left above the sea (a lowstand deposit) is land, as in the physics' routing.
                restored = fields['sea'] & (fields['surface'] <= min(level for _, level in sea_level['history']))
                if sea_level['keep_shelf'] == 'incision':
                    # Only the cutting stays: the lowstand's rivers filled F-troughs' shallow bay with deposits
                    # that kept it land after the rise (sl-v1, the user: "the bay emptied"). The exposed shelf
                    # keeps the lower of the physics' surface and its procedural depth, so it is sea again.
                    shelf = fields['sea'] & ~restored
                    elevation = np.where(shelf, np.minimum(elevation, fields['surface']), elevation)
            elevation = np.where(restored, fields['surface'], elevation)
            labels, _ = ndimage.label(elevation <= 0.0, structure=np.ones((3, 3)))
            touching = np.unique(labels[(labels > 0) & fields['sea']])
            sea = restored | np.isin(labels, touching[touching > 0])
            fields['sea'] = sea
            coast = config['coast']
            if coast['amplitude_m'] > 0 and coast.get('after_physics', False):
                # The coast term again on the final shoreline: troughs the physics cut below the sea and the sea
                # flooded have the smooth walls of an ice-carved valley (F-troughs' inlet), which the term before
                # the physics never saw. Then the sea is flooded again from the sea it had.
                x, y = fields['coords']
                elevation = elevation + coast_detail(domain, elevation, elevation, x, y, fields['amplitude'],
                                                     fields['erodibility'], coast, config['seed'] + 1)
                labels, _ = ndimage.label(elevation <= 0.0, structure=np.ones((3, 3)))
                touching = np.unique(labels[(labels > 0) & sea])
                sea = np.isin(labels, touching[touching > 0])
                fields['sea'] = sea
            micro = config.get('micro', {'enabled': False})
            if micro['enabled'] and micro.get('when', 'after') == 'after':
                # The shoreline texture on the final coast: drawn before the physics, 20 Myr of erosion and
                # marine filling took its few-km headlands and pockets back to a smooth arc (F-troughs, v16).
                (xw, yw), (xp, yp) = fields['warped'], fields['plain']
                # The open sea: the final sea less every part narrower than about open_width_km (a morphological
                # opening), so the glacial troughs it flooded inland (G-coast) are not shores, while a bay is.
                steps = max(1, int(round(micro['open_width_km'] / 2.0 / domain.fine_km)))
                open_sea = ndimage.binary_opening(sea, structure=np.ones((3, 3)), iterations=steps) & procedural_sea.any()
                elevation = elevation + micro_texture(domain, elevation, sea, xw, yw, xp, yp, fields['ratio'], micro,
                                                      config['seed'], open_sea=open_sea)
                labels, _ = ndimage.label(elevation <= 0.0, structure=np.ones((3, 3)))
                touching = np.unique(labels[(labels > 0) & sea])
                sea = np.isin(labels, touching[touching > 0])
                fields['sea'] = sea
            spec = config['filter']
            extra = []
            np.save(directory / f'{domain.name}-elevation-physics.npy', crop(domain, elevation).astype(np.float32))
            if spec['enabled'] and spec.get('when', 'before') == 'after':  # fine texture on the eroded surface
                x, y = fields['coords']
                field = filter_field(domain, fields, ice, spec, config['seed'])
                detail = erosion_pass(domain, elevation, x, y, fields['amplitude'], fields['erodibility'], spec,
                                      config['seed'], field=field)
                land = ~sea
                if spec.get('mean_preserve', 'none') == 'smooth':
                    # Take out the filter's broad component (its mean over about a world cell, land only), so
                    # the averages over world cells stay close to the physics' (not exact: the pipeline's
                    # projection does that).
                    sigma = domain.k / 2.0
                    weight = ndimage.gaussian_filter(land.astype(float), sigma, mode='nearest')
                    broad = ndimage.gaussian_filter(np.where(land, detail, 0.0), sigma, mode='nearest') / \
                        np.maximum(weight, 1e-6)
                    detail = detail - broad
                elevation = np.where(sea, elevation, np.maximum(elevation + detail, 1.0))
                if not domain.whole:
                    save_map(domain, directory / f'{domain.name}-filter-field.png',
                             colourize(crop(domain, field), FIELD_STOPS))
                    extra = [('filter strength field (black 0, white 1, orange 2 and more)', f'{domain.name}-filter-field.png')]
            np.save(directory / f'{domain.name}-elevation.npy', crop(domain, elevation).astype(np.float32))
            after = render(domain, elevation, sea, discharge, directory, ice=ice)
            images = [('elevation', after['elevation']), ('relief (hillshade, NW light)', after['relief']),
                      ('discharge (log, stretched from the median to the maximum)', after['discharge'])] + extra + images
        else:
            images = [('elevation (procedural, no physics)', before['elevation'])] + images
        timing['total'] = time.perf_counter() - started
        timings[domain.name] = {key: (round(value, 2) if isinstance(value, float) else value) for key, value in timing.items()}
        title = 'Whole world' if domain.whole else f'Window {domain.name}'
        sections.append((f'{title} ({domain.columns if domain.whole else domain.crop[3] - domain.crop[2]} columns, '
                         f'{domain.fine_km:.2f} km cells, k {domain.k})', images))
        if config['render'].get('diagnostics', True):
            sections.append((f'{title}: what gives the terrain its direction', render_diagnostics(domain, fields, config,
                                                                                                    directory)))
        if domain.whole and config['render'].get('stages', True):
            stage_images = render_stages(domain, fields, elevation if physics_on else None, config, directory)
            sections.append(('Whole world, stage by stage (each map the surface so far)', stage_images))
        print(f'{domain.name}: {timings[domain.name]}', flush=True)
    (directory / 'config.json').write_text(json.dumps(config, indent=1), encoding='utf-8')
    write_page(directory, arguments.label, sections, config, timings)
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
