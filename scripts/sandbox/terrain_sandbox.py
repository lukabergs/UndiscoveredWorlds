"""Terrain sandbox: procedural terrain laid on the tectonic world, then the pipeline's drainage and erosion
(terrain_sandbox.exe, src/cli/terrain_sandbox_main.cpp), rendered on the whole world and on windows.
It is for iterating on the look in seconds, not a pipeline stage: nothing is registered and the multiscale
contract (exact block means, pinned rivers, seamless windows) is not enforced. The base is a world bundle
stopped before the climate (config base_bundle, out/scratch/sandbox/base-v9: procedural-world, seed 1,
256 columns, the W3 tectonics baseline at 400 Myr, --steps ..terrain.roughness).

Order per domain (the whole world at a preview resolution, each window at a scale k with a margin):
the base surface (cubic B-spline of the world's surface elevation), a slow position warp, the relief noise
(isotropic, blended into noise stretched along the strike where the crust is shortened), the erosion
filter (erosion_filter.py) on base + relief, a rock-hardness field (relative erodibility, banded along the
strike), then routing and stream power. Every noise is a function of absolute position, so a window is a
zoom of the world preview; octaves finer than a domain's grid fade out instead of aliasing.

Config groups: domains, seafloor, warp, structure, zones, relief, coast, micro, uplands, hardness, platform, filter,
physics, sea_level, history, glacial, precipitation, lakes, rivers, render, nesting.
platform (roadmap T1 v19): a tectonic stability field (0 to 1) and provinces (oceanic, shield, platform, sedimentary
basin, orogen, rift / active margin) from the world's fields (provinces_world), sampled smoothly per domain and saved
as <domain>-stability.npy and <domain>-province.npy; the relief amplitude on land is damped by x 1 - relief_damp x
stability (the filters read that amplitude) and the physics' erodibility by x 1 - (1 - erodibility_factor) x stability.
physics.aggrade_subsidence (--aggrade-subsidence): land cells' negative uplift is filled by sediment instead of
lowering them (the per-cell total basin_fill_m, the report's basin_fill_km3). sea_level: a history of [Myr, m] points
(piecewise linear, a repeated time a jump) the physics runs through (--sea-history), the default the main
erosion at 0 m, a lowstand at -120 m for 3.5 Myr, then the rise to 0; when enabled the physics' duration is
the history's last time (global and windows alike). At the lowstand the shelf emerges and the rivers and
glaciers cut down to the lowstand; the rise drowns their lower valleys (estuaries, rias, longer fjords).
final_glacial_myr: the history's last stretch is a final glacial pulse (--no-infill-after) in which closed
depressions get no deposition, so the ice re-digs the overdeepenings the lowstand's rivers filled (the
linear lakes). keep_shelf: the sea shallower than the lowest level keeps the physics' surface instead of
its procedural depth ("incision": the lower of the two), so the valleys cut across the exposed shelf stay
drawn. lakes: after the physics and the micro texture, every closed depression not joined to the sea is a
lake at its spill level (Priority-Flood); the filter after the physics spares them, and they are found again
on the final surface and drawn flat. rivers: drawn on the elevation and relief maps with a width from the
discharge (hydraulic geometry) times a stated exaggeration.
history (roadmap T1.2b; supersedes sea_level when enabled): phases, each from start_myr to the next one's
start (the last to end.time_myr, the run's duration, where the sea jumps to end.sea_m): the sea level sea_m
(--sea-history, constant within a phase), the ice (true: the glacial ELA, false: none, a number: metres added
to the ELA) with glacial_ratio (--ice-history), the step step_myr (--step, then --step-after where it changes)
and no_infill (from the first such phase to the end: --no-infill-after); keep_shelf as sea_level's.
physics.lake_trap (--lake-trap): every closed depression of the routing keeps the sediment reaching it, filled
from its floor up to at most its spill level, and only the excess passes on (not in no_infill steps).
physics.lake_trap_ice_share (--lake-trap-ice-share): a depression with at least that share of its cells under ice
does not trap (the ice carries its sediment on, so glacial troughs stay deep and drown as fjords); absent: every
lake traps (the trap alone filled G's troughs and removed its fjords, 2026-10-16).
physics.lake_fill_gradient (--lake-fill-gradient, with physics.lake_fill_noise: --lake-fill-noise, the relief's
fill noise): the trap's fill slopes up from the lake's exits at that gradient, so a filled lake floor is no flat
(the ice cut straight grooves along a flat's routing lines). relief.fill_depressions: before the physics, every
closed depression not joined to the sea is filled from its spill, sloped at relief.fill_gradient with
relief.fill_noise (amplitude_m on an fBm) on it (sloped_fill): the procedural pits as sediment-filled basins.
The physics' input surface is saved as <domain>-elevation-procedural.npy.
filter.fill_pits (roadmap T1.6, default off): the closed pits the post-physics filter digs (at the foot of slopes,
where its strength falls with the slope: slope_fade) are filled as the relief's (sloped_fill) from the sea and the physics'
lakes at their spill level. render.save_extra (default off): also saves <domain>-discharge.npy (the final routing's,
m3/s) and <domain>-fill_m.npy (the depth relief.fill_depressions added).
Every physics run also writes per-cell totals (--out-fields), saved cropped as <domain>-<name>.npy (eroded_m,
deposited_m, glacial_eroded_m, ice_time_myr, ice_max_m, landslide_m, last_eroded_m) and drawn as three maps.

Nesting (roadmap T1.8, default off: nesting.enabled false draws every domain as above, byte for byte). One zoom
mechanism for every level below the global domain: a region of a region is a map zoomed in. Levels: L0 the
256-column world bundle (156 km cells), L1 the global domain (2048 columns, 19.5 km), L2 the k 64 windows (2.44 km,
factor 8 from L1), L3 the k 256 windows (0.61 km, factor 4 from L2); L4, 100 m, is nested_window.py. A window domain
names its parent in domains.<name>.parent; with nesting on, the domains run parents first (run_order) and a child
reads its parent's arrays from the same run's directory, or from --parent-run LABEL for a parent not in --domains
(otherwise the parent is run too).
- Rule A (L1, a domain without a parent): today's procedural chain (base, sea floor, warp, relief, uplands, coast,
  fill) and the full history, unchanged.
- Rule B (every nested level, nested_start): start = the Catmull-Rom prolongation of the parent's physics surface
  (<parent>-elevation-physics.npy: before the post-physics filter, which is a render detail and digs pits;
  clamped to the 2 x 2 bracketing parents) + the level's detail (the relief noise and the coast term restricted to
  the level's band of octaves, fbm's band_parent_km: in `complement` mode the share of each octave the parent's
  grid faded out, so the levels' shares sum to the finest grid's noise; `cut` keeps the octaves up to
  band.cut_parent_cells parent cells; on the warped coordinates with the oriented bank, so the directionality is
  the parent's; the ridged share split as mix(whole) - mix(whole - band); the sandbox's amplitude law with the
  stability damp; hardness, uplift, runoff, ELA and provinces evaluated as functions of position by procedural();
  no base, sea floor or uplands: they are in the parent; the relief band and the coast displacement fade to 0 within
  one parent cell of the parent's prolonged shore) + the sloped fill of the new pits (the parent's lakes kept); the
  start's sea only under the parent's sea and its coastal parent cells (children of other parent land cells stay
  land; those below projection.land_floor_m are raised to it; a parent sea cell that is a one-cell drowned valley
  of the coarse D8, nesting.one_cell_inlets, counts as land at 1 m, its cell mean included) -> the projection of the
  whole start onto the parent's cell means (project_and_refill: nesting.project_smooth's smooth pre-step of the
  multiscale plan section 1, projection.smooth_passes passes at projection.relaxation, then the exact remainder in
  proportion to each child's height above its block's lowest sea cell, uniform in a block without one; land never
  below projection.land_floor_m, sea never above projection.sea_ceiling_m, no child rising above its bracketing
  parents' maximum, the rest of a block's departure forced over its children (flips counted) but never more than
  projection.spike_clip_m beyond the bracketing parents' range (clipped blocks counted, not exact); the outer sea
  ring and the pinned ring fixed; then the new closed pits refilled, the parent's lakes and closed depressions
  kept, the exact remainder again and a second refill, the block means checked per block class after each step)
  -> the physics for the tail of the one history from physics.levels[k].from_myr (k 64 from
  20 Myr: the Quaternary only, the Neogene's valleys come from the parent; k 256 from 12.6; the level's sea and ice
  histories start at its own start) with the level's diffusivity, K and landslide threshold (physics.levels[k];
  on unnested windows too with physics.use_levels), the pinned ring at the parent's prolonged heights, the parent's
  rivers entering as runoff on the first movable cell inside the pinned ring (nesting.inflow_raster; the exe has
  no inflow option, so the discharge goes through the runoff raster, and the ELA there is held at the cell's
  start height so the inflow makes no glacier or ablation sink) -> the sea restored, the micro texture with the
  level's band only -> with nesting.project_after_physics (default false since round three of T1.8: with it on,
  t18-nest3 still drew the parent's cells as tiles, block-shaped lakes and terraces; t18-nest3-noafter did not) the
  physics surface projected again onto the parent's cell means by the same sequence, a rise moving only land that
  is no channel, lake or ring cell, a fall lowering the channels too (nested_project_after: the contract's "keep the
  parent's cell means" on the final map) -> saved as
  -elevation-physics.npy (what a child reads), then the filter and the lakes as today. Each nested domain writes
  <domain>-nest.json (departures before and after the projection and after the physics, the inflow and its
  continuity) and <domain>-start-relief.png; with nesting on every domain also saves <domain>-sea.npy and
  <domain>-discharge.npy for its children.
- Why two rules: the world bundle's elevation is the legacy terrain chain's output, which this sandbox exists to
  replace (the uplands and the depression fill change the hypsometry on purpose), so L1 does not keep L0's cell
  means; it is the first level the sandbox owns. After roadmap T2 the split becomes the multiscale plan's own:
  rule A the global stages, rule B the regional levels.
- Open (for the user, 2026-10-16): the global domain's physics (K 3000 on 19.5 km cells, 22.6 Myr) lowers G-coast's
  box from a 2562 m procedural mean to 850 m, and every nested window inherits that denudation.

  uv run python scripts/sandbox/terrain_sandbox.py --label try1 [--config FILE] [--set group.key=value ...]
      [--domains global,A,B] [--no-physics] [--parent-run LABEL]
Writes out/scratch/sandbox/<label>/: index.html, the PNGs, config.json (the values used) and timings.
"""
import argparse
import copy
import heapq
import json
import math
import shutil
import subprocess
import sys
import time
from collections import deque
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
import nesting  # noqa: E402
import overlays  # noqa: E402

DEFAULT_CONFIG = ROOT / 'configs' / 'sandbox' / 'terrain.json'
PHYSICS = ROOT / 'out' / 'build' / 'headless' / 'bin' / 'Release' / 'terrain_sandbox.exe'
NOISE_STD = 0.22  # standard deviation of one octave of the gradient noise below, measured


# ---------------------------------------------------------------- configuration

def load_config(path, overrides):
    """The config with GROUP.KEY=VALUE overrides; a numeric part indexes a list (history.phases.1.glacial_ratio=35)."""
    config = json.loads(Path(path).read_text(encoding='utf-8'))
    for item in overrides:
        key, _, value = item.partition('=')
        node = config
        parts = [int(part) if part.isdigit() else part for part in key.split('.')]
        for part in parts[:-1]:
            node = node[part]
        try:
            node[parts[-1]] = json.loads(value)
        except json.JSONDecodeError:
            node[parts[-1]] = value
    return config


PRESETS = ROOT / 'configs' / 'sandbox' / 'presets.json'


def preset_overrides(names):
    """The --set list of each named preset in configs/sandbox/presets.json, in order."""
    if not names:
        return []
    presets = json.loads(PRESETS.read_text(encoding='utf-8'))['presets']
    overrides = []
    for name in names.split(','):
        if name not in presets:
            raise SystemExit(f'Unknown preset {name}; known: {", ".join(sorted(presets))}')
        overrides += presets[name]['set']
    return overrides


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
        elif nesting_on(config) and spec.get('parent'):
            # A nested level (rule B): the margin in parent cells (nesting.margin_parent_cells by k), so the domain
            # lies on the parent's cell edges; the physics values of its level (physics.levels by k) over the window's.
            c0, r0, c1, r1 = spec['pixels']
            k, parent_k = spec['scale'], domain_k(config, spec['parent'], W)
            if k % parent_k:
                raise ValueError(f'{name}: scale {k} is not a multiple of its parent\'s {parent_k}')
            factor = k // parent_k
            margin = int(config['nesting'].get('margin_parent_cells', {}).get(str(k), spec['margin'] * parent_k))
            columns = ((c1 - c0 + 1) * parent_k + 2 * margin) * factor
            rows = ((r1 - r0 + 1) * parent_k + 2 * margin) * factor
            crop = (margin * factor, rows - margin * factor, margin * factor, columns - margin * factor)
            domains.append(Domain(name, columns, rows, c0 - margin / parent_k, r0 - margin / parent_k, k, W, H,
                                  radius_km, crop=crop, exaggeration=spec.get('exaggeration', 3.0),
                                  physics={**config['physics']['window'], **level_physics(config, k)}))
        else:
            c0, r0, c1, r1 = spec['pixels']
            k, m = spec['scale'], spec['margin']
            columns, rows = (c1 - c0 + 1 + 2 * m) * k, (r1 - r0 + 1 + 2 * m) * k
            crop = (m * k, rows - m * k, m * k, columns - m * k)
            physics = config['physics']['window']
            if config['physics'].get('use_levels', False):  # the level's kappa, K and threshold on an unnested window
                physics = {**physics, **level_physics(config, k)}
            domains.append(Domain(name, columns, rows, c0 - m, r0 - m, k, W, H, radius_km, crop=crop,
                                  exaggeration=spec.get('exaggeration', 3.0), physics=physics))
    return domains


# ---------------------------------------------------------------- nesting (rule B)

def nesting_on(config):
    return bool(config.get('nesting', {}).get('enabled', False))


def domain_k(config, name, world_columns):
    spec = config['domains'][name]
    return spec['columns'] // world_columns if name == 'global' else spec['scale']


def domain_depth(config, name):
    """0 for a domain without a parent (rule A), else its parent's depth + 1."""
    parent = config['domains'][name].get('parent')
    return 0 if not parent else 1 + domain_depth(config, parent)


def level_physics(config, k):
    """The physics values of the level with k cells per world cell (physics.levels[str(k)]): diffusivity, k and the
    landslide threshold go to the domain's physics group (a nested domain's always, an unnested window's with
    physics.use_levels); step_myr and from_myr to a nested domain's history only (level_config)."""
    level = config['physics'].get('levels', {}).get(str(k), {})
    return {key: level[key] for key in ('diffusivity_km2_myr', 'k', 'threshold_slope') if key in level}


def level_config(config, domain):
    """The config a nested level's physics runs with: the tail of the one history from physics.levels[k].from_myr
    (phases shifted to start at 0, the earlier ones dropped) and, with step_myr, no phase step longer than it. The
    sea and ice histories are built from these phases, so the level's ice history starts at its own start (a phase
    cut by from_myr keeps its ice and sea level from the level's time 0)."""
    level = config['physics'].get('levels', {}).get(str(domain.k), {})
    start, step = float(level.get('from_myr', 0.0)), level.get('step_myr')
    if start <= 0.0 and step is None:
        return config
    run = copy.deepcopy(config)
    history = run.get('history', {'enabled': False})
    if history.get('enabled', False):
        end = history['end']['time_myr']
        if start >= end:
            raise ValueError(f'{domain.name}: from_myr {start} is not before the history\'s end {end}')
        starts = [phase['start_myr'] for phase in history['phases']] + [end]
        phases = []
        for phase, until in zip(history['phases'], starts[1:]):
            if until <= start:
                continue
            phases.append({**phase, 'start_myr': max(phase['start_myr'] - start, 0.0)})
        if step is not None:
            phases = [{**phase, 'step_myr': min(phase.get('step_myr') or step, step)} for phase in phases]
        history['phases'] = phases
        history['end'] = {**history['end'], 'time_myr': end - start}
    else:
        window = run['physics']['window']
        window['duration_myr'] = window['duration_myr'] - start
        if step is not None:
            window['step_myr'] = min(window['step_myr'], step)
    return run


def run_order(config, names, parent_run):
    """The domains to run, parents first: a nested domain whose parent is neither requested nor read from
    --parent-run brings its parent (and so on up) into the run."""
    wanted = list(names)
    if nesting_on(config):
        queue = list(names)
        while queue:
            parent = config['domains'][queue.pop()].get('parent')
            if parent and parent not in wanted and parent_run is None:
                wanted.append(parent)
                queue.append(parent)
        wanted.sort(key=lambda name: domain_depth(config, name))  # stable: the requested order within a level
    return wanted


def load_parent(config, directory, name, world_columns):
    """A parent domain's saved arrays (cropped to its kept window) and where they lie: its physics surface
    (-elevation-physics.npy, before the post-physics filter), sea (-sea.npy, saved with nesting; else the cells at
    or below 0 m joined to the crop's border), lakes (-lake-surface.npy, NaN off lakes; none when missing) and
    discharge (-discharge.npy, m3/s; None when missing: no inflow)."""
    elevation = np.load(directory / f'{name}-elevation-physics.npy').astype(np.float64)
    path = directory / f'{name}-sea.npy'
    if path.is_file():
        sea = np.load(path)
    else:
        sea = nesting.flood(elevation <= 0.0, nesting.ring(elevation.shape, 1))
    path = directory / f'{name}-lake-surface.npy'
    lakes = np.load(path) if path.is_file() else np.full(elevation.shape, np.nan)
    path = directory / f'{name}-discharge.npy'
    discharge = np.load(path).astype(np.float64) if path.is_file() else None
    spec = config['domains'][name]
    origin = (0.0, 0.0) if name == 'global' else (float(spec['pixels'][0]), float(spec['pixels'][1]))
    return {'name': name, 'directory': str(directory), 'elevation': elevation, 'sea': sea, 'lakes': lakes,
            'discharge': discharge, 'k': domain_k(config, name, world_columns), 'origin': origin}


def band_spec(config, parent_fine_km):
    band = config['nesting'].get('band', {})
    return {'band_parent_km': parent_fine_km, 'band_mode': band.get('mode', 'complement'),
            'band_cut_cells': band.get('cut_parent_cells', 2.0)}


def nested_start(domain, world, config, parent, seed):
    """Rule B's start surface of a nested level (docstring, "Nesting"). Returns the fields (as procedural's, with
    the start surface, its sea, the pinned ring, the inflow and the ELA at the inflow cells), the timings and the
    nesting diagnostics."""
    timings, info = {}, {}
    pk = parent['k']
    f = domain.k // pk
    prow, pcol = domain.rows // f, domain.columns // f
    x0 = domain.x[0] - 0.5 / domain.k
    y0 = domain.y[0] - 0.5 / domain.k
    pc = int(round((x0 - parent['origin'][0]) * pk))
    pr = int(round((y0 - parent['origin'][1]) * pk))
    P = parent['elevation']
    if not (0 <= pr and pr + prow <= P.shape[0] and 0 <= pc and pc + pcol <= P.shape[1]):
        raise ValueError(f'{domain.name}: parent cells [{pr}, {pr + prow}) x [{pc}, {pc + pcol}) outside '
                         f'{parent["name"]}\'s saved {P.shape}')
    margin = domain.crop[0] // f
    inner = (slice(margin, prow - margin), slice(margin, pcol - margin))
    info.update(parent=parent['name'], parent_directory=parent['directory'], factor=f,
                parent_cells=[pr, pc, prow, pcol], margin_parent_cells=margin)
    spec = config['nesting'].get('projection', {})
    land_floor = spec.get('land_floor_m', 1.0)
    # A parent sea cell that is a one-cell drowned valley of the coarse D8 (nesting.one_cell_inlets: a one-cell-wide
    # chain not touching the open sea) is land at the sea level + 1 m for this level: its classification, its
    # prolongation and its cell mean (else the projection pulls its children back below the sea and the physics
    # floods them again: the staircase inlets of round one).
    stair = nesting.one_cell_inlets(parent['sea'])
    box = (slice(pr, pr + prow), slice(pc, pc + pcol))
    info['parent_one_cell_inlets'] = {'cells': int(stair[box].sum()),
                                      'largest_lift_m': float((1.0 - P[box])[stair[box]].max()) if stair[box].any()
                                      else 0.0}
    P = np.where(stair, 1.0, P)
    parent_sea = parent['sea'] & ~stair

    started = time.perf_counter()
    prolonged = nesting.prolong(P, f, pr, pc, prow, pcol)
    under_sea = nesting.nearest(parent_sea, f, pr, pc, prow, pcol)
    # Children of a parent land cell stay land unless the parent's sea reaches them: sea only under the parent's
    # sea and its 8-neighbouring (coastal) parent cells, flooded from under the parent's sea.
    reached = parent_sea | ndimage.binary_dilation(parent_sea, structure=np.ones((3, 3)))
    allowed = nesting.nearest(reached, f, pr, pc, prow, pcol)
    sea0 = nesting.flood((prolonged <= 0.0) & allowed, under_sea & (prolonged <= 0.0))
    # The level's relief band and the coast displacement fade to 0 within one parent cell of the parent's shore
    # (the prolonged shore), so the detail does not move the coast across parent cells.
    if sea0.any() and not sea0.all():
        shore = np.where(sea0, ndimage.distance_transform_edt(sea0), ndimage.distance_transform_edt(~sea0))
        fade = smoothstep(shore / f, 0.0, 1.0)
    else:
        fade = np.ones(prolonged.shape)
    timings['prolong'] = time.perf_counter() - started

    # The fields that are functions of position (rock, uplift, rain, provinces, warp, the relief's amplitude and
    # noise) as procedural() draws them; no base, sea floor, uplands, coast, filter or fill: the surface is the
    # parent's.
    level = copy.deepcopy(config)
    level['seafloor']['enabled'] = False
    level['seafloor'].setdefault('shelf', {})['enabled'] = False
    level['uplands']['lift_m'] = 0
    level['coast'] = {**level['coast'], 'amplitude_m': 0}
    level['filter'] = {**level['filter'], 'enabled': False}
    level['relief'] = {**level['relief'], 'fill_depressions': False}
    level['micro'] = {**level.get('micro', {}), 'enabled': False}
    fields, timing = procedural(domain, world, level, seed)
    timings.update({key: value for key, value in timing.items() if isinstance(value, float)})

    # The level's detail: the relief noise's band and the coast term's (band_spec), both on the warped coordinates
    # and the oriented bank, so the directionality is the parent's.
    started = time.perf_counter()
    stages = fields['stages']
    x, y = fields['coords']
    theta, alignment = stages['strike']
    _, core = stages['zones']
    relief, structure = config['relief'], config['structure']
    band = band_spec(config, domain.cell_km / pk)
    relief_alignment = alignment * (1 - config['zones']['core_isotropy'] * core)
    whole = stages['relief_noise']
    part = oriented_fbm(x, y, theta, relief_alignment, {**relief, **structure, **band}, seed + 20, domain)
    mix = lambda n: (1 - relief['ridged']) * n + relief['ridged'] * ridged(n) / 0.6  # noqa: E731
    detail = fields['amplitude'] * (mix(whole) - mix(whole - part)) * fade
    surface = prolonged + detail
    coast = config['coast']
    if coast['amplitude_m'] > 0:
        fraction = coast.get('warp_fraction', 1.0)
        (xp, yp) = fields['plain']
        xc, yc = (x, y) if fraction == 1.0 else (xp + fraction * (x - xp), yp + fraction * (y - yp))
        surface = surface + coast_detail(domain, fields['base'], surface, xc, yc, fields['amplitude'],
                                         fields['erodibility'], {**coast, **band}, seed, displacement_fade=fade)
    land0 = ~sea0
    info['detail_rms_m'] = {'relief_land': float(np.sqrt((detail[land0] ** 2).mean())) if land0.any() else 0.0,
                            'total_land': float(np.sqrt(((surface - prolonged)[land0] ** 2).mean())) if land0.any() else 0.0}
    timings['detail'] = time.perf_counter() - started

    # Sea, the outlet ring and the pinned ring (nested_window.py's step 5).
    sea = nesting.flood((surface <= 0.0) & allowed, sea0)
    outer = nesting.ring(surface.shape, 1)
    pinned = nesting.ring(surface.shape, 2) & ~outer & ~sea
    surface = np.where(pinned, prolonged, surface)
    surface = np.where(outer, np.minimum(prolonged, -1.0), surface)
    sea |= outer
    # Land children of a parent land cell below the land floor (the detail's troughs) are raised to it, so the
    # physics' sea flood does not take them; a parent cell below 0 m (an inland basin) keeps its children.
    low = ~sea & ~pinned & (surface < land_floor) & nesting.nearest(P > 0.0, f, pr, pc, prow, pcol)
    surface = np.where(low, land_floor, surface)
    info['start_land_raised_to_floor'] = int(np.count_nonzero(low))

    # The new pits filled (sloped, as relief.fill_depressions), the parent's lakes kept.
    started = time.perf_counter()
    filled = sloped_fill(domain, surface, sea, relief['fill_gradient'], fields['fill_noise'])
    labels, count = ndimage.label(filled - surface > 1e-3, structure=np.ones((3, 3)))
    under_lake = np.isfinite(nesting.nearest(parent['lakes'], f, pr, pc, prow, pcol))
    # The parent's own closed depressions (its Priority-Flood raises them by more than 1 m; the box's edge drains)
    # are kept as its lakes are: a refill would otherwise undo the projection of a parent pit (the global's one-cell
    # pits to -644 m) every time.
    parent_box = P[pr:pr + prow, pc:pc + pcol]
    drains = parent_sea[pr:pr + prow, pc:pc + pcol] | nesting.ring(parent_box.shape, 1)
    hollow = (priority_flood(parent_box, drains, False) - parent_box > 1.0) & ~drains
    under_lake = under_lake | np.repeat(np.repeat(hollow, f, axis=0), f, axis=1)
    info['parent_depression_cells'] = int(hollow.sum())
    kept = np.unique(labels[under_lake & (labels > 0)])
    kept_mask = np.isin(labels, kept[kept > 0])
    unfilled = surface
    surface = np.where(kept_mask, surface, filled)
    surface = np.where(pinned, prolonged, surface)
    fields['fill_stats'] = {**fill_statistics(domain, unfilled, surface, sea), 'depressions': int(count),
                            'kept_with_parent_lake': int(np.count_nonzero(kept))}
    timings['fill'] = time.perf_counter() - started

    # The exact projection onto the parent's cell means (land never below land_floor_m, sea never above
    # sea_ceiling_m; nesting.project_bounded).
    targets = P[pr:pr + prow, pc:pc + pcol]
    movable = ~pinned & ~outer
    before = nesting.departure(surface, targets, movable, f, inner)
    # No child rises above its bracketing parents' maximum (the rest of a block's departure forced over all its
    # children, flipping sea children to land), and none more than spike_clip_m above it (the block's departure
    # clipped there: such a block is not exact; counted).
    bottom, top = nesting.bracket_extremes(P, f, pr, pc, prow, pcol)
    original_sea = sea

    def reflood(values, water):  # forced blocks flip children: the sea flooded again from what stayed at or below 0
        return nesting.flood(values <= 0.0, water & (values <= 0.0)) | outer

    surface, sea, clipped, stats, checks = project_and_refill(
        domain, surface, targets, lambda water: (movable, movable, ~water, water), sea, reflood, (bottom, top), spec, f,
        under_lake, pinned, prolonged, config, fields['fill_noise'], np.zeros_like(sea))
    stats['flipped_to_land'] = int(np.count_nonzero(original_sea & ~sea))
    stats['flipped_to_sea'] = int(np.count_nonzero(~original_sea & sea))
    info['start_block_check_m'] = checks
    everywhere = (slice(0, prow), slice(0, pcol))
    info['start_departure_m'] = {'before_projection': before,
                                 'after_projection_kept': nesting.departure(surface, targets, movable, f, inner),
                                 'after_projection_all': nesting.departure(surface, targets, movable, f, everywhere),
                                 'after_projection_unclipped': nesting.departure(surface, targets, movable & ~clipped, f,
                                                                                 everywhere),
                                 'projection': stats,
                                 'prolongation': nesting.departure(prolonged, targets, movable, f, inner)}

    # The parent's rivers entering at the edge (nesting.inflow_raster), the ELA held at the inflow cells' start height.
    inflow, cells, records = None, None, []
    if parent['discharge'] is not None and config['nesting'].get('inflow', {}).get('enabled', True):
        found = nesting.entries(parent['discharge'], pr, pc, prow, pcol,
                                config['nesting'].get('inflow', {}).get('min_m3_s', 1.0))
        areas = nesting.row_areas_m2(domain.north_deg, domain.pixel_deg, domain.rows, domain.radius_km)
        inflow, cells, records = nesting.inflow_raster(found, pr, pc, f, surface.shape, surface, movable & ~sea, areas)
    info['inflow'] = {'entries': len(records), 'placed': sum(1 for r in records if r['child_cell']),
                      'parent_m3_s_total': float(sum(r['parent_m3_s'] for r in records)),
                      'placed_m3_s_total': float(sum(r['parent_m3_s'] for r in records if r['child_cell'])),
                      'records': records}
    fields.update(surface=surface, sea=sea, pinned=pinned, inflow_m_yr=inflow,
                  ela_fixed=None if cells is None or not cells.any() else (cells, surface))
    stages['filter'] = unfilled  # save_extra's fill_m: what the fill added
    stages['fill'] = surface
    fields['nest'] = {'targets': targets, 'inner': inner, 'factor': f, 'prolonged': prolonged, 'top': top,
                      'bottom': bottom,
                      'outer': outer, 'parent_lakes': under_lake,
                      'parent_discharge': None if parent['discharge'] is None else
                      parent['discharge'][pr:pr + prow, pc:pc + pcol]}
    return fields, timings, info


def project_and_refill(domain, surface, targets, masks, sea, reflood, extremes, spec, f, keep, pinned, pinned_values,
                       config, fill_noise, lake):
    """Rule B's projection sequence (the start's and the one after the physics): nesting.project_smooth (the smooth
    pre-step and the weighted exact remainder), the new closed pits refilled (refill_pits), the weighted exact
    remainder again (nesting.project_weighted; the refill raised the pits' blocks), the pits refilled once more, so
    the lake finder sees no pit the projection made. `masks(sea)` gives (movable_up, movable_down, land, reference);
    `reflood(values, sea)` the sea after a projection. Returns the surface, the sea, the clipped mask, the counts
    and the block-mean check per block class after each step."""
    floor, ceiling = spec.get('land_floor_m', 1.0), spec.get('sea_ceiling_m', -1.0)
    bottom, top = extremes  # the bracketing parents' least and greatest: no child forced 300 m beyond either
    clip, clip_low = top + spec.get('spike_clip_m', 300.0), bottom - spec.get('spike_clip_m', 300.0)
    rise, fall, land, reference = masks(sea)
    projected, stats, clipped = nesting.project_smooth(
        surface, targets, rise, fall, land, reference, f, floor, ceiling, top, clip, spec.get('smooth_passes', 6),
        spec.get('relaxation', 1.4), spec.get('weight_floor_m', 10.0), clip_low=clip_low)
    sea = reflood(projected, sea)
    checks = {'after_projection': nesting.block_check(projected, targets, f, fall & ~sea, sea, lake, clipped)}
    refilled, refill = refill_pits(domain, projected, sea, keep, pinned, pinned_values, config, fill_noise)
    checks['after_refill'] = nesting.block_check(refilled, targets, f, fall & ~sea, sea, lake, clipped)
    rise, fall, land, reference = masks(sea)
    weights = nesting.reference_weights(refilled, reference, rise | fall, f)
    again, stats2, clipped2 = nesting.project_weighted(refilled, targets, rise, fall, land, weights, f, floor, ceiling,
                                                       cap=top, clip=clip, clip_low=clip_low)
    sea = reflood(again, sea)
    clipped = clipped | clipped2
    checks['after_second_remainder'] = nesting.block_check(again, targets, f, fall & ~sea, sea, lake, clipped)
    final, refill2 = refill_pits(domain, again, sea, keep, pinned, pinned_values, config, fill_noise)
    checks['final_after_second_refill'] = nesting.block_check(final, targets, f, fall & ~sea, sea, lake, clipped)
    moved = (final - surface)[rise | fall]
    stats.update(refill=refill, second_remainder={key: stats2[key] for key in ('forced_blocks', 'clipped_blocks',
                                                                               'largest_shift_m')},
                 second_refill=refill2, total_correction_max_m=float(np.abs(moved).max()) if moved.size else 0.0,
                 total_correction_rms_m=float(np.sqrt((moved ** 2).mean())) if moved.size else 0.0)
    return final, sea, clipped, stats, checks


def refill_pits(domain, surface, sea, keep, pinned, pinned_values, config, fill_noise):
    """The closed pits a projection left, filled with the sloped fill (relief.fill_gradient and its noise) from the
    sea, except the depressions overlapping `keep` (the parent's lakes, the physics' lakes); the pinned ring back at
    its values. Returns the surface and counts."""
    filled = sloped_fill(domain, surface, sea, config['relief']['fill_gradient'], fill_noise)
    labels, count = ndimage.label(filled - surface > 1e-3, structure=np.ones((3, 3)))
    kept = np.unique(labels[keep & (labels > 0)])
    out = np.where(np.isin(labels, kept[kept > 0]), surface, filled)
    out = np.where(pinned, pinned_values, out)
    raised = out - surface
    return out, {'pits': int(count), 'kept': int(np.count_nonzero(kept[kept > 0])),
                 'raised_cells': int(np.count_nonzero(raised > 1e-3)), 'raised_max_m': float(raised.max())}


def nested_project_after(domain, fields, elevation, sea, discharge, config):
    """The contract's "keep the parent's cell means" on the final map (nesting.project_after_physics): the child's
    physics surface (sea restored, micro texture on) projected onto the parent's cell means as the start is
    (nesting.project_smooth: the smooth pre-step, then the exact remainder in proportion to the height above each
    block's lowest channel, sea or lake cell; the land floor, the bracketing parents' maximum, the spike clip). A
    rise moves the land that is no channel, lake or ring cell; a fall also lowers the channels (in proportion to
    their height above the block's lowest), so a lowering leaves no river in a slot; the sea and the lakes stay
    fixed. Then the new closed pits are filled (refill_pits; the physics' and the parent's lakes kept). Returns the
    surface and the report (correction max and rms, the block-mean check per block class after the projection and
    after the refill)."""
    nest = fields['nest']
    f, targets = nest['factor'], nest['targets']
    spec = config['nesting'].get('projection', {})
    lakes = config.get('lakes', {'enabled': False})
    lake = np.isfinite(find_lakes(domain, elevation, sea, lakes)) if lakes['enabled'] else np.zeros_like(sea)
    channel = discharge >= river_threshold(domain, config.get('rivers', {'min_discharge_m3_s': 30.0}))
    free = ~sea & ~lake & ~fields['pinned'] & ~nest['outer']
    refilled, _, _, stats, checks = project_and_refill(
        domain, elevation, targets, lambda water: (free & ~channel, free, ~water, channel | water | lake), sea,
        lambda values, water: water, (nest['bottom'], nest['top']), spec, f, lake | nest['parent_lakes'], fields['pinned'], elevation,
        config, fields['fill_noise'], lake)
    report = {key: stats[key] for key in ('correction_max_m', 'correction_rms_m', 'correction_mean_m',
                                          'total_correction_max_m', 'total_correction_rms_m')}
    report.update(channel_cells=int(channel.sum()), lake_cells=int(lake.sum()), projection=stats, block_check_m=checks)
    return refilled, report


def nested_after_physics(domain, fields, elevation, sea, discharge, info):
    """Rule B's checks after the physics: the departure of the parent-cell means (the physics' surface with the
    sea restored, kept window, blocks more than half land), the inflow's continuity (the child's discharge at each
    entry cell and the most in the inside parent cell, against the parent's) and the discharge consistency in the
    kept window (the most child discharge per parent cell against the parent's, where the parent has 100 m3/s)."""
    nest = fields['nest']
    f, inner, targets = nest['factor'], nest['inner'], nest['targets']
    change = (nesting.block_mean(elevation, f) - targets)[inner]
    land = (nesting.block_mean((~sea).astype(float), f) > 0.5)[inner]
    out = {}
    if land.any():
        out['physics_departure_m_land_blocks'] = {
            'mean': float(change[land].mean()), 'p5': float(np.percentile(change[land], 5)),
            'p50': float(np.percentile(change[land], 50)), 'p95': float(np.percentile(change[land], 95)),
            'max_abs': float(np.abs(change[land]).max())}
    rows, columns = discharge.shape
    block_max = discharge.reshape(rows // f, f, columns // f, f).max(axis=(1, 3))
    checks = []
    for record in info['inflow']['records']:
        if not record['child_cell']:
            continue
        i, j = record['child_cell']
        r, c = record['inside'][0] - info['parent_cells'][0], record['inside'][1] - info['parent_cells'][1]
        parent_inside = None if nest['parent_discharge'] is None else float(nest['parent_discharge'][r, c])
        checks.append({'parent_entering_m3_s': record['parent_m3_s'], 'child_at_entry_m3_s': float(discharge[i, j]),
                       'parent_inside_m3_s': parent_inside, 'child_inside_max_m3_s': float(block_max[r, c])})
    if checks:
        entry = np.array([c['child_at_entry_m3_s'] / c['parent_entering_m3_s'] for c in checks])
        big = [c for c in checks if c['parent_entering_m3_s'] >= 100.0]
        out['inflow_continuity'] = {
            'entries': len(checks), 'entry_ratio_p10_p50_p90': [float(np.percentile(entry, q)) for q in (10, 50, 90)],
            'largest': sorted(checks, key=lambda c: -c['parent_entering_m3_s'])[:8],
            'entries_100_m3_s': len(big)}
    if nest['parent_discharge'] is not None:
        parent_q = nest['parent_discharge'][inner]
        child_q = block_max[inner]
        rivers = parent_q >= 100.0
        if rivers.any():
            ratio = child_q[rivers] / parent_q[rivers]
            out['discharge_kept_parent_cells_100_m3_s'] = {
                'cells': int(rivers.sum()), 'child_max_over_parent_p10_p50_p90':
                    [float(np.percentile(ratio, q)) for q in (10, 50, 90)],
                'parent_max_m3_s': float(parent_q.max()), 'child_max_m3_s': float(child_q.max())}
    return out


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
    optional = lambda name, dtype: bundle.load(f'geosphere.{name}').astype(dtype) if bundle.has(f'geosphere.{name}') \
        else np.zeros(elevation.shape, dtype=dtype)  # noqa: E731
    return {'elevation': elevation, 'sea': sea, 'shortening': np.maximum(strain, 0.0), 'direction': direction,
            'uplift': uplift, 'sediment': bundle.load('geosphere.sediment').astype(float),
            'continental_crust': bundle.load('geosphere.continental_crust').astype(float),
            'crust_class': bundle.load('geosphere.crust_class'), 'crust_age': bundle.load('geosphere.crust_age').astype(float),
            'boundary_type': bundle.load('geosphere.boundary_type'),
            'accumulated_strain': optional('accumulated_strain', float),
            'deforming_region_type': optional('deforming_region_type', np.uint8)}


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


def band_weight(wavelength_km, spec):
    """A nested level's share of an octave (nesting, rule B): with spec['band_parent_km'] (the parent level's cell),
    `complement` mode keeps 1 - the parent grid's octave fade (what the parent could not draw, so the levels' shares
    of every octave sum to the finest grid's), `cut` mode keeps the octaves no longer than band_cut_cells parent
    cells. 1 without the key."""
    parent_km = spec.get('band_parent_km')
    if parent_km is None:
        return 1.0
    if spec.get('band_mode', 'complement') == 'cut':
        return 1.0 if wavelength_km <= spec.get('band_cut_cells', 2.0) * parent_km else 0.0
    min_cells = spec.get('min_cells', 2.5)
    return 1.0 - float(np.clip((wavelength_km / parent_km - min_cells) / min_cells, 0.0, 1.0))


def fbm(x, y, spec, seed, domain, stretch=1.0):
    """Fractal sum with unit standard deviation (about); x across and y along a stretch of `stretch`.
    Each octave is rotated by a fixed irrational angle so no lattice direction repeats. With
    spec['band_parent_km'] (a nested level, rule B) only the level's band of octaves is summed (band_weight), with
    the normalisation of the whole sum on this grid, so the band is the part of this grid's noise the parent lacks."""
    total = np.zeros_like(x)
    norm = 0.0
    wavelength, amplitude = spec['wavelength_km'], 1.0
    banded = spec.get('band_parent_km') is not None
    for octave in range(spec['octaves']):
        fade = octave_fade(wavelength, domain, spec.get('min_cells', 2.5))
        if banded and fade > 0.0 and band_weight(wavelength, spec) <= 0.0:
            norm += (amplitude * fade) ** 2
            wavelength /= spec.get('lacunarity', 2.0)
            amplitude *= spec.get('gain', 0.5)
            continue
        if fade > 0.0:
            if stretch == 1.0:  # each octave rotated by the golden angle
                angle = 2.399963 * octave
                ca, sa = math.cos(angle), math.sin(angle)
                u, v = (ca * x - sa * y) / wavelength, (sa * x + ca * y) / wavelength
            else:  # a stretched octave keeps the stretch's axes: shifted, not rotated
                u, v = x / wavelength + 0.37 * octave, y / (wavelength * stretch) + 0.61 * octave
            if banded:
                total += amplitude * fade * band_weight(wavelength, spec) * gradient_noise(u, v, seed * 131 + octave)
            else:
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
    stages['relief_noise'] = noise  # before the ridged share: a nested level splits it by band (nested_start)
    noise = (1 - relief['ridged']) * noise + relief['ridged'] * ridged(noise) / 0.6
    amplitude = relief['base_m'] + relief['elevation_fraction'] * np.maximum(base, 0.0) + \
        relief['strain_m'] * alignment
    # The sea's smaller relief fades in with depth (from +200 m down to -500 m): a switch at the base's zero
    # contour jumped the noise along the smooth old coastline, a crease in the sea floor and the coast.
    wet = 1.0 - smoothstep(base, -500.0, 200.0)
    amplitude = amplitude * (1.0 - (1.0 - relief['sea_fraction']) * wet)
    # Stable platforms keep low relief (platform group, roadmap T1 v19): the relief amplitude damped by the tectonic
    # stability on land (x 1 - relief_damp x stability; the sea floor unchanged). The filters read this amplitude too
    # (strength x amplitude / relief_reference_m), so their gullies are damped alike.
    platform = config.get('platform', {'enabled': False})
    stability, province = None, None
    if platform['enabled']:
        stability, province = provinces_domain(domain, world['provinces'], platform, x, y, mx, my, seed)
        amplitude = amplitude * (1.0 - platform['relief_damp'] * stability * (1.0 - wet))
        stages['provinces'] = (stability, province)
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
    # The physics' erodibility on stable ground (old basement and flat-lying cover): x erodibility_factor at full
    # stability, so the rivers do not dissect the platforms into badlands. The drawn hardness and the coast and filter
    # terms keep the rock's own value.
    physics_erodibility = erodibility if stability is None else \
        erodibility * (1.0 - (1.0 - platform['erodibility_factor']) * stability)
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
    # The procedural pits as sediment-filled basins (relief.fill_depressions, roadmap T1 v18): 23-29 % of the land
    # lay in closed depressions up to 5.8 km deep (the ridged relief noise), and rivers entering one incised toward
    # its floor, below the sea. The fill slopes up from the spill at relief.fill_gradient with relief.fill_noise on
    # it; the same noise goes to the physics' lake fill (physics.lake_fill_gradient).
    relief_spec = config['relief']
    noise_spec = relief_spec.get('fill_noise', {})
    fill_noise = None
    if noise_spec.get('amplitude_m', 0.0) > 0.0:
        fill_noise = noise_spec['amplitude_m'] * fbm(x0, y0, noise_spec, seed + 90, domain)
    fill_stats = {}
    if relief_spec.get('fill_depressions', False):
        started = time.perf_counter()
        unfilled = filtered
        filtered = sloped_fill(domain, filtered, sea, relief_spec['fill_gradient'], fill_noise)
        fill_stats = fill_statistics(domain, unfilled, filtered, sea)
        timings['fill'] = time.perf_counter() - started
    stages['fill'] = filtered
    return {'base': base, 'surface': filtered, 'sea': sea, 'erodibility': erodibility, 'uplift': uplift,
            'physics_erodibility': physics_erodibility, 'stability': stability, 'province': province,
            'strain': strain, 'stages': stages, 'precipitation': rain, 'coords': (x, y), 'amplitude': amplitude,
            'pinned': pinned, 'warped': (x, y), 'plain': (x0, y0), 'ratio': ratio, 'fill_noise': fill_noise,
            'fill_stats': fill_stats}, timings


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


def coast_detail(domain, base, surface, x, y, relief_m, erodibility, spec, seed, displacement_fade=None):
    """Relief concentrated around the sea level, so the shoreline is a fractal line rather than a contour of
    smooth fields: a rough fBm (gain near 0.6 keeps the fine octaves strong, a ridged share makes headlands
    and skerries) times exp(-(surface / band_m)^2) on the surface so far. Its amplitude follows the terrain:
    hard rock (erodibility^hardness_exponent, negative: hard coasts are rocky, soft ones smooth), the local
    relief (1 + relief_gain x relief / relief_reference_m: mountains meeting the sea) and cold coasts (1 +
    cold_gain x smoothstep(|latitude|, cold_latitude): the skerry and fjord coasts of former ice). Pure
    function of position, so every zoom draws the same coast down to its grid's resolution. With
    spec['band_parent_km'] (a nested level, rule B) the vertical term is the difference the level's band makes to the
    ridged mix (of the whole noise on this grid less the band) and the displacement is the band's alone."""
    noise = fbm(x, y, spec, seed + 70, domain)
    if spec.get('band_parent_km') is not None:
        whole = fbm(x, y, {key: value for key, value in spec.items() if key != 'band_parent_km'}, seed + 70, domain)
        mix = lambda n: (1 - spec['ridged']) * n + spec['ridged'] * ridged(n) / 0.6  # noqa: E731
        noise = mix(whole) - mix(whole - noise)
    else:
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
        if displacement_fade is not None:  # a nested level: none near the parent's shore (nested_start)
            scale = scale * displacement_fade
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


# ---------------------------------------------------------------- provinces (craton and platform field)

PROVINCES = ('oceanic', 'shield', 'platform', 'sedimentary basin', 'orogen', 'rift / active margin')
PROVINCE_COLOURS = ((40, 62, 110), (196, 84, 84), (226, 190, 120), (246, 238, 160), (128, 76, 150), (70, 170, 165))


def boundary_distance_km(boundary_type, code, cell_km):
    """Distance (km) from every world cell to the nearest boundary cell of `code`, wrapping east-west (the engine's
    planar cells, cell_km each way); infinite when the world has none."""
    mask = boundary_type == code
    if not mask.any():
        return np.full(mask.shape, np.inf)
    columns = mask.shape[1]
    distance = ndimage.distance_transform_edt(~np.tile(mask, (1, 3)))[:, columns:2 * columns]
    return distance * cell_km


def provinces_world(world, spec, cell_km):
    """The tectonic stability (0 to 1) and the province memberships on the world grid (roadmap T1 v19). Stability is
    a product of factors, each 1 on stable ground: continental crust of normal thickness (smoothstep up through
    crust_m[0:2], down through crust_m[2:4]: oceanic, thinned and orogen-thickened crust are not stable), little
    shortening (the strain, shortening over the cell width, below strain), little recent strain (the engine's
    decaying accumulated_strain), far from active boundaries (convergent, divergent, transform: smoothstep of the
    distance in km through boundary_km[type]), little uplift or subsidence (|rock_uplift_rate| through uplift_m_myr),
    outside the engine's deforming regions (x 1 - deforming_weight) and older crust (1 - age_weight below age_myr[0],
    1 above age_myr[1]); then smoothed over smoothing_cells. The provinces: oceanic (continental crust thinner than
    provinces.ocean_crust_m); stable ground (stability at least provinces.stable) is a sedimentary basin where it
    subsides (rock uplift at most provinces.basin_uplift_m_myr) or carries thick sediment (provinces.basin_sediment_m),
    a platform with sediment cover (provinces.cover_m), else a shield (exposed basement); unstable ground an orogen
    where the compressional factors (strain, recent strain, convergent boundary, thick crust, uplift, a collision
    region) explain more of the instability than the extensional ones (divergent and transform boundaries, thin crust,
    subsidence, a rift region), else a rift or active margin. Returns stability, memberships (one per PROVINCES,
    smoothed like the stability) and the hard codes."""
    smooth = lambda field, sigma: ndimage.gaussian_filter(field, sigma, mode=('nearest', 'wrap'))  # noqa: E731
    crust = world['continental_crust']
    c0, c1, c2, c3 = spec['crust_m']
    thin = 1.0 - smoothstep(crust, c0, c1)
    thick = smoothstep(crust, c2, c3)
    strain = smoothstep(world['shortening'], *spec['strain'])
    recent = smoothstep(world['accumulated_strain'], *spec['accumulated_strain'])
    rate = world['uplift']
    moving = smoothstep(np.abs(rate), *spec['uplift_m_myr'])
    near = {}
    for name, code in (('convergent', 1), ('divergent', 2), ('transform', 3)):
        near[name] = 1.0 - smoothstep(boundary_distance_km(world['boundary_type'], code, cell_km), *spec['boundary_km'][name])
    deforming = world['deforming_region_type']
    a0, a1 = spec['age_myr']
    age = 1.0 - spec['age_weight'] * (1.0 - smoothstep(world['crust_age'], a0, a1))
    raw = (1 - thin) * (1 - thick) * (1 - strain) * (1 - recent) * (1 - moving) * age * \
        (1 - near['convergent']) * (1 - near['divergent']) * (1 - near['transform']) * \
        (1 - spec['deforming_weight'] * (deforming > 0))
    stability = np.clip(smooth(raw, spec['smoothing_cells']), 0.0, 1.0)
    rules = spec['provinces']
    oceanic = crust < rules['ocean_crust_m']
    stable = ~oceanic & (stability >= rules['stable'])
    basin = stable & ((rate <= rules['basin_uplift_m_myr']) | (world['sediment'] >= rules['basin_sediment_m']))
    platform = stable & ~basin & (world['sediment'] >= rules['cover_m'])
    shield = stable & ~basin & ~platform
    compression = np.maximum.reduce([strain, recent, near['convergent'], thick, moving * (rate > 0),
                                     (deforming == 2).astype(float)])
    extension = np.maximum.reduce([near['divergent'], near['transform'], thin, moving * (rate < 0),
                                   (deforming == 1).astype(float)])
    unstable = ~oceanic & ~stable
    orogen = unstable & (compression >= extension)
    rift = unstable & ~orogen
    codes = np.zeros(crust.shape, dtype=np.uint8)
    for code, mask in enumerate((oceanic, shield, platform, basin, orogen, rift)):
        codes[mask] = code
    memberships = [smooth((codes == code).astype(float), spec['smoothing_cells']) for code in range(len(PROVINCES))]
    return {'stability': stability, 'memberships': memberships, 'codes': codes}


def provinces_domain(domain, provinces, spec, x, y, mx, my, seed):
    """The stability (cubic B-spline of the world's, smooth) and the province codes on a domain: each province's
    smoothed membership B-spline sampled, plus provinces.edge_noise x an fBm per province on the (warped) noise
    coordinates so the edges are not smooth contours, then the largest."""
    stability = np.clip(sample(provinces['stability'], domain, mx, my), 0.0, 1.0)
    noise_spec = spec['provinces']['noise']
    amplitude = spec['provinces']['edge_noise']
    scores = [sample(field, domain, mx, my) + (amplitude * fbm(x, y, noise_spec, seed + 95 + code, domain)
                                               if amplitude > 0 else 0.0)
              for code, field in enumerate(provinces['memberships'])]
    return stability, np.argmax(np.stack(scores), axis=0).astype(np.uint8)


def province_shares(codes, land, weights):
    """Area shares (%) of the provinces on land (weights: the cells' areas)."""
    total = float(weights[land].sum())
    return {name: round(100.0 * float(weights[land & (codes == code)].sum()) / max(total, 1e-9), 1)
            for code, name in enumerate(PROVINCES)}


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

NO_ICE_OFFSET_M = 100000.0  # an ELA offset that leaves no ice (a phase with "ice": false)
TOTAL_FIELDS = ('eroded_m', 'deposited_m', 'glacial_eroded_m', 'ice_time_myr', 'ice_max_m', 'landslide_m',
                'last_eroded_m')


def schedule(config):
    """The physics' time schedule, or None for the plain run at 0 m: from the history group when enabled
    (roadmap T1.2b; it supersedes sea_level), else from sea_level when enabled. A dict with sea_history
    ([Myr, m] points), duration, no_infill_after (Myr or None), keep_shelf and, from the history only,
    ice_history ([Myr, ELA offset m, glacial ratio] per phase), step (the first phase's step, Myr) and
    step_after ([Myr, step] where the step changes)."""
    history = config.get('history', {'enabled': False})
    if history.get('enabled', False):
        phases, end = history['phases'], history['end']
        starts = [phase['start_myr'] for phase in phases] + [end['time_myr']]
        if any(b <= a for a, b in zip(starts, starts[1:])):
            raise ValueError('history: phase starts and the end must increase')
        sea, ice, step_after = [], [], []
        step = phases[0].get('step_myr')
        current = step
        ratio = config['glacial']['erodibility_ratio']
        for phase, until in zip(phases, starts[1:]):
            sea += [[phase['start_myr'], phase['sea_m']], [until, phase['sea_m']]]
            value = phase.get('ice', True)
            offset = NO_ICE_OFFSET_M if value is False else 0.0 if value is True else float(value)
            ratio = phase.get('glacial_ratio') or ratio
            ice.append([phase['start_myr'], offset, ratio])
            if phase.get('step_myr') and phase['step_myr'] != current:
                current = phase['step_myr']
                step_after.append([phase['start_myr'], current])
        sea.append([end['time_myr'], end['sea_m']])
        flags = [bool(phase.get('no_infill', False)) for phase in phases]
        no_infill_after = None
        if any(flags):
            first = flags.index(True)
            if not all(flags[first:]):
                raise ValueError('history: no_infill must hold from its first phase to the end')
            no_infill_after = phases[first]['start_myr']
        return {'sea_history': sea, 'duration': end['time_myr'], 'no_infill_after': no_infill_after,
                'keep_shelf': history.get('keep_shelf', 'incision'), 'ice_history': ice, 'step': step,
                'step_after': step_after}
    sea_level = config.get('sea_level', {'enabled': False})
    if sea_level.get('enabled', False):
        duration = sea_level['history'][-1][0]
        final = sea_level.get('final_glacial_myr', 0)
        return {'sea_history': sea_level['history'], 'duration': duration,
                'no_infill_after': duration - final if final > 0 else None,
                'keep_shelf': sea_level.get('keep_shelf', False)}
    return None


def run_physics(domain, fields, config, directory):
    spec = {**config['physics'], **domain.physics}
    stem = directory / f'{domain.name}'
    files = {}
    runoff = np.maximum(fields['precipitation'] * config['precipitation']['runoff_fraction'],
                        config['precipitation']['runoff_min_m_yr'])
    if fields.get('inflow_m_yr') is not None:  # a nested level's parent rivers entering at its edge (rule B)
        runoff = runoff + fields['inflow_m_yr']
    for key, values in (('elevation', fields['surface']), ('sea', fields['sea'].astype(np.float32)),
                        ('erodibility', fields.get('physics_erodibility', fields['erodibility'])),
                        ('uplift', fields['uplift'] * spec['uplift_factor']),
                        ('runoff', runoff), ('pinned', fields['pinned'].astype(np.float32))):
        files[key] = stem.with_name(f'{domain.name}-{key}.f32')
        np.asarray(values, dtype='<f4').tofile(files[key])
    out_elevation = stem.with_name(f'{domain.name}-out-elevation.f32')
    out_discharge = stem.with_name(f'{domain.name}-out-discharge.f32')
    plan = schedule(config)
    duration = plan['duration'] if plan else spec['duration_myr']
    step = (plan or {}).get('step') or spec['step_myr']
    totals_prefix = stem.with_name(f'{domain.name}-totals')
    command = [str(PHYSICS), '--columns', str(domain.columns), '--rows', str(domain.rows),
               '--elevation', str(files['elevation']), '--sea', str(files['sea']),
               '--erodibility', str(files['erodibility']), '--uplift', str(files['uplift']),
               '--runoff-file', str(files['runoff']), '--duration', str(duration),
               '--step', str(step), '--k', str(spec['k']), '--m', str(spec['m']),
               '--deposition', str(spec['deposition']), '--diffusivity', str(spec['diffusivity_km2_myr']),
               '--delta-top', str(spec['delta_top_m']), '--threshold-slope', str(spec['threshold_slope']),
               '--pinned', str(files['pinned']), '--radius-km', str(domain.radius_km),
               '--random-routing', '1' if spec.get('random_routing', False) else '0', '--seed', str(config['seed']),
               '--out-elevation', str(out_elevation), '--out-discharge', str(out_discharge),
               '--out-fields', str(totals_prefix)]
    if not domain.whole:
        command += ['--north', str(domain.north_deg), '--pixel-deg', str(domain.pixel_deg)]
    if spec.get('lake_trap', False):  # lakes keep the sediment reaching them (filled up to their spill)
        command += ['--lake-trap', '1']
        if spec.get('lake_trap_ice_share') is not None:  # lakes under ice stay open (fjords; T1.7)
            command += ['--lake-trap-ice-share', str(spec['lake_trap_ice_share'])]
        if spec.get('lake_fill_gradient', 0.0) > 0.0:  # the fill sloped up from the spill, with the fill noise
            command += ['--lake-fill-gradient', str(spec['lake_fill_gradient'])]
            if spec.get('lake_fill_noise', False) and fields.get('fill_noise') is not None:
                files['fill_noise'] = stem.with_name(f'{domain.name}-fill-noise.f32')
                np.asarray(fields['fill_noise'], dtype='<f4').tofile(files['fill_noise'])
                command += ['--lake-fill-noise', str(files['fill_noise'])]
    names = TOTAL_FIELDS
    if spec.get('aggrade_subsidence', False):  # land cells' subsidence filled by sediment (basin_fill_m)
        command += ['--aggrade-subsidence', '1']
        names = TOTAL_FIELDS + ('basin_fill_m',)
    if plan:
        command += ['--sea-history', ','.join(f'{time}:{level}' for time, level in plan['sea_history'])]
        if plan['no_infill_after'] is not None:  # the final glacial pulse: no infill of closed basins
            command += ['--no-infill-after', str(plan['no_infill_after'])]
        if plan.get('step_after'):
            command += ['--step-after', ','.join(f'{time}:{length}' for time, length in plan['step_after'])]
    glacial = config['glacial']
    out_ice = stem.with_name(f'{domain.name}-out-ice.f32')
    ela = None
    if glacial['enabled']:
        ela = equilibrium_line(domain, fields['precipitation'], glacial)
        if fields.get('ela_fixed') is not None:
            # The inflow cells' line at their start height: the exe's ice source is the runoff times (surface - ELA)
            # / full_ice_range, so a river's whole discharge would otherwise become a glacier or an ablation sink.
            mask, values = fields['ela_fixed']
            ela = np.where(mask, values, ela)
        files['ela'] = stem.with_name(f'{domain.name}-ela.f32')
        np.asarray(ela, dtype='<f4').tofile(files['ela'])
        command += ['--ela', str(files['ela']), '--glacial-ratio', str(glacial['erodibility_ratio']),
                    '--ice-full-range', str(glacial['full_ice_range_m']),
                    '--ice-accumulation', str(glacial['ice_accumulation']), '--out-ice', str(out_ice)]
        if plan and plan.get('ice_history'):
            command += ['--ice-history', ','.join(f'{time}:{offset}:{ratio}' for time, offset, ratio in plan['ice_history'])]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'terrain_sandbox failed on {domain.name}: {result.stderr.strip()}')
    shape = (domain.rows, domain.columns)
    elevation = np.fromfile(out_elevation, dtype='<f4').reshape(shape).astype(np.float64)
    discharge = np.fromfile(out_discharge, dtype='<f4').reshape(shape).astype(np.float64)
    ice = np.fromfile(out_ice, dtype='<f4').reshape(shape).astype(np.float64) if glacial['enabled'] else None
    total_files = [Path(f'{totals_prefix}-{name}.f32') for name in names]
    totals = {name: np.fromfile(path, dtype='<f4').reshape(shape) for name, path in zip(names, total_files)}
    for path in list(files.values()) + [out_elevation, out_discharge] + ([out_ice] if glacial['enabled'] else []) + \
            total_files:
        path.unlink()
    if ice is not None:  # the maps show today's glaciers: the ice where the surface is above today's line
        ice = np.where(elevation + ice > ela + glacial['present_shift_m'], ice, 0.0)
    return elevation, discharge, ice, totals, result.stdout.strip()


def priority_flood(elevation, sea, wrap):
    """The filled surface of a terrain (Priority-Flood, Barnes et al. 2014, 8-connected): every land cell
    raised to the lowest level at which water on it reaches a sea cell, the sea cells the outlets at their own
    elevation. Cells raised (in a depression) go through a plain queue, the others through the heap. A whole
    world wraps east-west through a padding of an eighth of its width on each side (a lake wider than that
    across the date line would be cut there); its polar rows are no outlets."""
    rows, columns = elevation.shape
    pad = columns // 8 if wrap else 0
    z = np.pad(np.asarray(elevation, dtype=float), ((1, 1), (pad + 1, pad + 1)),
               mode='wrap' if wrap else 'edge')
    water = np.pad(sea, ((1, 1), (pad + 1, pad + 1)), mode='wrap' if wrap else 'edge')
    done = water.copy()
    done[0, :] = done[-1, :] = done[:, 0] = done[:, -1] = True  # the frame: no cell, no outlet
    width = z.shape[1]
    offsets = (-width - 1, -width, -width + 1, -1, 1, width - 1, width, width + 1)
    shore = water & ndimage.binary_dilation(~water, structure=np.ones((3, 3)))
    shore[0, :] = shore[-1, :] = shore[:, 0] = shore[:, -1] = False
    height = z.ravel().tolist()
    filled = list(height)
    closed = done.ravel().tolist()
    heap = [(height[i], i) for i in np.flatnonzero(shore.ravel()).tolist()]
    heapq.heapify(heap)
    pit = deque()
    while heap or pit:
        if pit:
            cell = pit.popleft()
            level = filled[cell]
        else:
            level, cell = heapq.heappop(heap)
        for offset in offsets:
            j = cell + offset
            if closed[j]:
                continue
            closed[j] = True
            if height[j] <= level:
                filled[j] = level
                pit.append(j)
            else:
                heapq.heappush(heap, (height[j], j))
    filled = np.asarray(filled).reshape(z.shape)[1:-1, pad + 1:pad + 1 + columns]
    return np.where(sea, elevation, filled)


def sloped_fill(domain, elevation, sea, gradient, noise=None):
    """The procedural surface with its closed depressions filled as sediment-filled basins (relief.fill_depressions,
    roadmap T1 v18): a Priority-Flood from the sea (as priority_flood) in which a cell lower than the cell that
    reaches it is raised to that cell's level plus max(0.1 g d, g d + n_j - n_i) (g the gradient, d the centre
    distance on the sphere, n the noise, m), and goes back into the heap at that level. The fill rises away from
    the spill at about g, with the noise's few metres on it so the routing on the floor follows the noise instead
    of the flood's straight lines; every raised cell has a lower neighbour (its parent), so no closed depression
    is left. The rising fill also buries low ground just above the flat spill level (the fill laps onto the
    basin's margins). Returns the filled surface (the sea and every cell that was not raised unchanged)."""
    rows, columns = elevation.shape
    wrap = domain.whole
    pad = columns // 8 if wrap else 0
    padding = ((1, 1), (pad + 1, pad + 1))
    mode = 'wrap' if wrap else 'edge'
    z = np.pad(np.asarray(elevation, dtype=float), padding, mode=mode)
    water = np.pad(sea, padding, mode=mode)
    n = np.pad(np.zeros(elevation.shape) if noise is None else np.asarray(noise, dtype=float), padding, mode=mode)
    done = water.copy()
    done[0, :] = done[-1, :] = done[:, 0] = done[:, -1] = True
    width = z.shape[1]
    dy = domain.fine_km * 1000.0
    dx = dy * np.pad(domain.coslat, (1, 1), mode='edge')  # per padded row
    offsets = (-width - 1, -width, -width + 1, -1, 1, width - 1, width, width + 1)
    rises = [[gradient * math.hypot(dx[row] if dc else 0.0, dy if dr else 0.0)
              for dr, dc in ((1, 1), (1, 0), (1, 1), (0, 1), (0, 1), (1, 1), (1, 0), (1, 1))] for row in range(z.shape[0])]
    shore = water & ndimage.binary_dilation(~water, structure=np.ones((3, 3)))
    shore[0, :] = shore[-1, :] = shore[:, 0] = shore[:, -1] = False
    height = z.ravel().tolist()
    lift = n.ravel().tolist()
    filled = list(height)
    closed = done.ravel().tolist()
    heap = [(height[i], i) for i in np.flatnonzero(shore.ravel()).tolist()]
    heapq.heapify(heap)
    while heap:
        level, cell = heapq.heappop(heap)
        run = rises[cell // width]
        for k, offset in enumerate(offsets):
            j = cell + offset
            if closed[j]:
                continue
            closed[j] = True
            if height[j] <= level:
                step = run[k] + lift[j] - lift[cell]
                filled[j] = level + max(0.1 * run[k], step)
                heapq.heappush(heap, (filled[j], j))
            else:
                heapq.heappush(heap, (height[j], j))
    filled = np.asarray(filled).reshape(z.shape)[1:-1, pad + 1:pad + 1 + columns]
    return np.where(sea, elevation, filled)


def fill_statistics(domain, before, after, sea):
    """What relief.fill_depressions added on the kept domain: volume (km^3), the raised share of the land (cells
    raised by more than 1 m), mean and largest raise (m), and the land's area-weighted mean before and after."""
    raise_m, land = crop(domain, after - before), ~crop(domain, sea)
    area = domain.fine_km ** 2 * crop(domain, domain.coslat[:, None] * np.ones((1, domain.columns)))
    raised = land & (raise_m > 1.0)
    return {'fill_km3': round(float((raise_m * area)[land].sum()) / 1e3, 1),
            'fill_land_pct': round(100.0 * float(area[raised].sum() / max(area[land].sum(), 1e-9)), 2),
            'fill_mean_m': round(float((raise_m * area)[raised].sum() / max(area[raised].sum(), 1e-9)), 1),
            'fill_max_m': round(float(raise_m.max()), 1),
            'land_mean_before_m': round(float((crop(domain, before) * area)[land].sum() / area[land].sum()), 1),
            'land_mean_after_m': round(float((crop(domain, after) * area)[land].sum() / area[land].sum()), 1)}


def find_lakes(domain, elevation, sea, spec, physics_lakes=None):
    """Lakes at their spill level (roadmap T1.3): the closed depressions the Priority-Flood fills (filled above
    the terrain by more than a millimetre), 8-connected, kept when deeper than min_depth_m somewhere and at
    least min_cells cells. Every lake is wet (no water balance in the sandbox). physics_lakes (a mask), when
    given, keeps only the depressions that overlap it: on the surface after the filter, the lakes the physics
    left, at the spill level the filter leaves them, without the filter's own small pits (thousands of them
    in a k 256 window, t13-control). Returns the lake surface (m; NaN off lakes)."""
    filled = priority_flood(elevation, sea, domain.whole)
    depth = np.where(sea, 0.0, filled - elevation)
    labels, count = ndimage.label(depth > 1e-3, structure=np.ones((3, 3)))
    if count == 0:
        return np.full(elevation.shape, np.nan)
    index = np.arange(1, count + 1)
    deepest = ndimage.maximum(depth, labels, index)
    cells = ndimage.sum(np.ones_like(depth), labels, index)
    keep = np.zeros(count + 1, bool)
    keep[1:] = (deepest > spec['min_depth_m']) & (cells >= spec['min_cells'])
    if physics_lakes is not None:
        overlap = np.zeros(count + 1, bool)
        overlap[np.unique(labels[physics_lakes & (labels > 0)])] = True
        keep &= overlap
    return np.where(keep[labels], filled, np.nan)


def lake_statistics(domain, elevation, surface):
    """Lake count, area (km^2) and the deepest lake (m) in the kept window."""
    surface, elevation = crop(domain, surface), crop(domain, elevation)
    lake = np.isfinite(surface)
    coslat = crop(domain, domain.coslat[:, None] * np.ones((1, domain.columns)))
    area = np.where(lake, domain.fine_km ** 2 * coslat, 0.0)
    labels, count = ndimage.label(lake, structure=np.ones((3, 3)))
    depth = np.where(lake, surface - elevation, 0.0)
    return {'lakes': int(count), 'lake_area_km2': round(float(area.sum()), 1),
            'deepest_lake_m': round(float(depth.max()) if count else 0.0, 1),
            'largest_lake_km2': round(float(ndimage.sum(area, labels, np.arange(1, count + 1)).max()) if count else 0.0, 1)}


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
                                            ('filter', 'erosion filter'),
                                            ('fill', 'closed depressions filled (sloped from the spill)'),
                                            ('physics', 'drainage and erosion')), 1):
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


# Lakes and rivers (roadmap T1.3), one water family distinct from the sea's pale blues and greens: the lakes in
# the lake-surface palette's colour at 0 m (configs/tectonics-palettes.json, lake_surface_m), the rivers a
# little lighter (the river_flow palette's 100 stop), so a river is seen entering its lake.
LAKE_RGB = np.array([30.0, 80.0, 160.0])
RIVER_RGB = np.array([40.0, 90.0, 170.0])


def river_threshold(domain, spec):
    """The least discharge drawn as a river (m3/s): min_discharge_m3_s, one value or one per kind of domain
    ({"global": ..., "window": ...}, as the physics group), so the world map's 20 km pixels are not a mesh."""
    value = spec['min_discharge_m3_s']
    return float(value['global' if domain.whole else 'window'] if isinstance(value, dict) else value)


def water_note(domain, lakes, rivers):
    """The caption's account of the water drawn on the elevation and relief maps."""
    parts = []
    if lakes.get('enabled', False):
        parts.append(f'lakes (dark blue) at their spill level, flat and unshaded (the closed depressions the '
                     f'physics left, deeper than {lakes["min_depth_m"]:g} m, at least {lakes["min_cells"]} cells)')
    if rivers.get('enabled', False):
        parts.append(f'rivers from {river_threshold(domain, rivers):g} m3/s, width {rivers["width_coefficient"]:g} '
                     f'Q^{rivers["width_exponent"]:g} m (hydraulic geometry) exaggerated {rivers["exaggeration"]:g} '
                     f'times, at least one pixel')
    return '; ' + '; '.join(parts) if parts else ''


def river_coverage(domain, discharge, water, spec):
    """Coverage 0..1 of the rivers on the kept window: every land cell with a discharge of at least
    min_discharge_m3_s joined to its downstream neighbour (the neighbour with the most discharge, when it has
    more than the cell) by a segment w = width_coefficient x Q^width_exponent m wide (Leopold and Maddock's
    hydraulic geometry; Q the cell's discharge, m3/s) times exaggeration, at least one pixel; the pixel's size
    is its north-south one. Drawn at `supersample` times the resolution with round joins, then box-averaged
    (antialiased). `water` (sea and lakes) takes no river."""
    from PIL import Image, ImageDraw
    rows, columns = discharge.shape
    q = np.where(water, 0.0, discharge)
    padded = np.pad(discharge, 1, constant_values=-np.inf)
    neighbours = np.stack([padded[1 + dr:1 + dr + rows, 1 + dc:1 + dc + columns]
                           for dr in (-1, 0, 1) for dc in (-1, 0, 1)])
    neighbours[4] = -np.inf  # the cell itself
    best = np.argmax(neighbours, axis=0)
    downstream = np.take_along_axis(neighbours, best[None], 0)[0] > discharge
    river = q >= river_threshold(domain, spec)
    scale = int(spec.get('supersample', 4))
    width_px = np.maximum(spec['width_coefficient'] * np.maximum(q, 0.0) ** spec['width_exponent'] *
                          spec['exaggeration'] / (domain.fine_km * 1000.0), 1.0) * scale
    image = Image.new('L', (columns * scale, rows * scale), 0)
    draw = ImageDraw.Draw(image)
    for r, c in zip(*np.nonzero(river)):
        x0, y0 = (c + 0.5) * scale, (r + 0.5) * scale
        width = float(width_px[r, c])
        half = width / 2.0
        ends = [(x0, y0)]
        if downstream[r, c]:  # a river's last cell (into the sea or a lake) is only its round end
            dr, dc = divmod(int(best[r, c]), 3)
            x1, y1 = (c + dc - 1 + 0.5) * scale, (r + dr - 1 + 0.5) * scale
            draw.line((x0, y0, x1, y1), fill=255, width=max(1, int(round(width))))
            ends.append((x1, y1))
        for x, y in ends:
            draw.ellipse((x - half, y - half, x + half, y + half), fill=255)
    coverage = np.asarray(image.resize((columns, rows), Image.BOX), dtype=float) / 255.0
    return np.where(water, 0.0, coverage)


def render(domain, elevation, sea, discharge, directory, suffix='', ice=None, lake_surface=None, rivers=None):
    elevation, sea = crop(domain, elevation), crop(domain, sea)
    ice = None if ice is None else crop(domain, ice)
    lake = np.zeros_like(sea) if lake_surface is None else np.isfinite(crop(domain, lake_surface))
    if lake.any():  # the relief of a lake is its flat surface
        elevation = np.where(lake, crop(domain, lake_surface), elevation)
    base_rgb = with_ice(terrain_base(elevation, sea), ice)
    coverage = None
    if rivers is not None and discharge is not None:
        coverage = river_coverage(domain, crop(domain, discharge), sea | lake, rivers)[..., None]

    def water(rgb):  # rivers, then the lakes over them
        if coverage is not None:
            rgb = rgb * (1.0 - coverage) + RIVER_RGB * coverage
        rgb[lake] = LAKE_RGB
        return rgb
    names = {}
    names['elevation'] = f'{domain.name}{suffix}-elevation.png'
    save_map(domain, directory / names['elevation'], water(base_rgb.copy()))
    names['relief'] = f'{domain.name}{suffix}-relief.png'
    coslat = np.maximum(crop(domain, domain.coslat[:, None] * np.ones((1, domain.columns)))[:, 0], 0.25)
    dy = domain.fine_km * 1000.0
    shade = shade_factor(hillshade(elevation, dy, dy * coslat[:elevation.shape[0]], domain.exaggeration))
    save_map(domain, directory / names['relief'], np.clip(water(base_rgb * shade[..., None]), 0, 255))  # the ice shaded too
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


DEPTH_STOPS = [(0.0, (12, 12, 18)), (0.25, (80, 30, 110)), (0.5, (190, 55, 80)), (0.75, (245, 150, 45)),
               (1.0, (252, 248, 205))]  # sequential, black to pale yellow, on 0..1 of the map's range
SHARE_STOPS = [(0.0, (12, 12, 18)), (0.5, (50, 110, 185)), (1.0, (235, 245, 255))]  # black to ice white


def render_totals(domain, totals, sea, directory):
    """The physics' per-cell totals on the kept window: the total eroded depth (stream-power erosion with its
    hillslope term, plus the landslides), the deposit thickness (land and marine) and the glacial share of the
    stream-power erosion. Each on a sequential ramp from 0 to the window's 99th percentile (the share 0 to 1);
    the final sea darkened toward blue so the coast reads."""
    sea = crop(domain, sea)
    values = {name: crop(domain, field).astype(np.float64) for name, field in totals.items()}

    def ramp(field, top, stops):
        rgb = colourize(np.clip(field / max(top, 1e-9), 0.0, 1.0), stops).astype(float)
        rgb[sea] = 0.5 * rgb[sea] + 0.5 * np.array([20.0, 40.0, 80.0])
        return rgb

    images = []
    eroded = values['eroded_m'] + values['landslide_m']
    deposited = values['deposited_m']
    stream = values['eroded_m']
    share = np.where(stream > 1.0, values['glacial_eroded_m'] / np.maximum(stream, 1e-9), np.nan)
    ice_time, ice_max = values['ice_time_myr'], values['ice_max_m']
    fill = values.get('basin_fill_m')
    basin = () if fill is None else (
        ('basin-fill', fill, float(max(np.percentile(fill, 99), 1.0)), DEPTH_STOPS,
         'basin fill, m (the subsidence of land cells filled by sediment, physics.aggrade_subsidence)'),)
    for key, field, top, stops, caption in (
            ('eroded', eroded, float(np.percentile(eroded, 99)), DEPTH_STOPS,
             'total eroded depth, m (stream power with its hillslope term, plus landslides)'),
            ('deposited', deposited, float(np.percentile(deposited, 99)), DEPTH_STOPS,
             'deposit thickness, m (land and marine; the final pulse\'s infill of closed basins taken off)'),
            ('glacial-share', share, 1.0, SHARE_STOPS,
             'glacial share of the stream-power erosion (erosion in steps under more than 10 m of ice over all '
             'of it; cells eroded less than 1 m black)')) + basin:
        rgb = ramp(np.nan_to_num(field), top, stops)
        name = f'{domain.name}-total-{key}.png'
        save_map(domain, directory / name, rgb)
        if key == 'glacial-share':
            text = (f'{caption}: black 0 to white 1; ice above 10 m for 0 to {ice_time.max():.2f} Myr, thickest '
                    f'{ice_max.max():.0f} m; final sea darkened')
        else:
            text = (f'{caption}: black 0 to pale yellow {top:.0f} m (99th percentile); range {field.min():.0f} to '
                    f'{field.max():.0f} m; final sea darkened')
        images.append((text, name))
    return images


def render_provinces(domain, fields, directory):
    """The province map (categorical, with its legend in the caption) and the stability (black 0 to white 1) on the
    kept window, the procedural sea darkened; saves <domain>-stability.npy and <domain>-province.npy (cropped). Returns
    the images and the provinces' shares of the procedural land (%, area-weighted)."""
    stability, province = crop(domain, fields['stability']), crop(domain, fields['province'])
    np.save(directory / f'{domain.name}-stability.npy', stability.astype(np.float32))
    np.save(directory / f'{domain.name}-province.npy', province.astype(np.uint8))
    sea = crop(domain, fields['sea'])
    weights = crop(domain, np.repeat(domain.coslat[:, None], domain.columns, axis=1))
    shares = province_shares(province, ~sea, weights)
    rgb = np.asarray(PROVINCE_COLOURS, dtype=float)[province]
    rgb[sea] = 0.45 * rgb[sea] + 0.55 * np.array([15.0, 20.0, 35.0])
    save_map(domain, directory / f'{domain.name}-province.png', rgb)
    legend = ' '.join(f'<span style="background:rgb{colour};display:inline-block;width:12px;height:12px"></span> '
                      f'{name} {shares[name]:.1f} %' for name, colour in zip(PROVINCES, PROVINCE_COLOURS))
    grey = colourize(stability, [(0.0, (10, 10, 14)), (0.5, (120, 120, 110)), (1.0, (250, 248, 235))]).astype(float)
    grey[sea] = 0.45 * grey[sea] + 0.55 * np.array([15.0, 40.0, 80.0])
    save_map(domain, directory / f'{domain.name}-stability.png', grey)
    return [(f'tectonic provinces (sea darkened; shares of the procedural land): {legend}', f'{domain.name}-province.png'),
            ('tectonic stability: black 0 (orogens, active margins) to white 1 (stable craton); it damps the relief '
             'amplitude and the physics\' erodibility; sea darkened', f'{domain.name}-stability.png')], shares


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
    note = (f'Seed 1; base: {config.get("base_bundle")}, procedural-world at 256 columns after the tectonics (base-v9: '
            'the W3 baseline, 400 Myr) and the terrain stages, before the climate and the global erosion. Uniform runoff, no climate. Windows are zooms of the same '
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
    parser.add_argument('--preset', default=None, metavar='NAME[,NAME]',
                        help='named looks from configs/sandbox/presets.json, applied before --set')
    parser.add_argument('--domains', default=None, help="comma-separated, default the config's default_domains")
    parser.add_argument('--no-physics', action='store_true')
    parser.add_argument('--parent-run', default=None, metavar='LABEL',
                        help='with nesting: read the parents not in --domains from out/scratch/sandbox/LABEL '
                             '(default: run them here)')
    arguments = parser.parse_args()

    config = load_config(arguments.config, preset_overrides(arguments.preset) + arguments.set)
    directory = ROOT / 'out' / 'scratch' / 'sandbox' / arguments.label
    directory.mkdir(parents=True, exist_ok=True)
    bundle = Bundle.open(ROOT / config['base_bundle'])
    world = world_fields(bundle)
    radius_m = float(bundle.manifest.get('parameters', {}).get('world.body.radius_km', 6371.0)) * 1000.0
    world['precipitation'] = precipitation(world, config['precipitation'], radius_m)
    if config['seafloor'].get('abyss', {}).get('enabled', False):
        world['elevation'] = damp_abyss(world, config['seafloor']['abyss'])
    world['ocean'] = ocean_world(world, config['seafloor'], config['seed'], radius_m / 1000.0)
    if config.get('platform', {}).get('enabled', False):
        world['provinces'] = provinces_world(world, config['platform'], 2.0 * math.pi * radius_m / 1000.0 / bundle.columns)
    names = arguments.domains.split(',') if arguments.domains else config['default_domains']
    names = run_order(config, names, arguments.parent_run)
    physics_on = config['physics']['enabled'] and not arguments.no_physics
    sections, timings = [], {}
    for domain in make_domains(config, bundle, names):
        started = time.perf_counter()
        parent_name = config['domains'][domain.name].get('parent') if nesting_on(config) else None
        run_config, nest_info, band = config, None, {}
        if parent_name:  # rule B (docstring, "Nesting")
            source = directory if parent_name in names else ROOT / 'out' / 'scratch' / 'sandbox' / arguments.parent_run
            parent = load_parent(config, source, parent_name, bundle.columns)
            fields, timing, nest_info = nested_start(domain, world, config, parent, config['seed'])
            del parent
            run_config = level_config(config, domain)
            band = band_spec(config, domain.fine_km * nest_info['factor'])
        else:
            fields, timing = procedural(domain, world, config, config['seed'])
        timing.update(fields['fill_stats'])
        images = []
        province_images = []
        if fields['province'] is not None:
            province_images, timing['province_land_pct'] = render_provinces(domain, fields, directory)
        if config['render'].get('before_physics', False) or not physics_on:
            before = render(domain, fields['surface'], fields['sea'], None, directory, '-procedural')
            if not physics_on:  # for scripts/sandbox/hypsometry.py and other measurements
                np.save(directory / f'{domain.name}-elevation.npy', crop(domain, fields['surface']).astype(np.float32))
            images += [('relief before physics', before['relief'])]
        if physics_on:
            physics_started = time.perf_counter()
            elevation, discharge, ice, totals, report = run_physics(domain, fields, run_config, directory)
            timing['physics'] = time.perf_counter() - physics_started
            timing['physics_report'] = report
            # The physics' input surface (after the depression fill), for measuring what the physics changed.
            np.save(directory / f'{domain.name}-elevation-procedural.npy', crop(domain, fields['surface']).astype(np.float32))
            for name, values in totals.items():  # per-cell totals for later mineral work
                np.save(directory / f'{domain.name}-{name}.npy', crop(domain, values).astype(np.float32))
            if config['render'].get('save_extra', False) or nesting_on(config):  # for measurements (roadmap
                # T1.6) and a nested child's inflow: the final routing's discharge (m3/s, the rivers drawn)
                np.save(directory / f'{domain.name}-discharge.npy', crop(domain, discharge).astype(np.float32))
            if config['render'].get('save_extra', False):  # the depth relief.fill_depressions added (m)
                np.save(directory / f'{domain.name}-fill_m.npy',
                        crop(domain, fields['surface'] - fields['stages']['filter']).astype(np.float32))
            # The sea keeps its procedural depth: the marine deposition's straight walks are not what is judged.
            # Land cut below the sea and connected to it (a fjord, a drowned river mouth) becomes sea.
            procedural_sea = fields['sea'].copy()  # the open sea before the physics flooded its troughs
            plan = schedule(run_config)
            restored = fields['sea']
            if plan and plan.get('keep_shelf', False):
                # The shelf the lowstand exposed keeps the physics' surface (the valleys cut across it, drowned
                # by the rise); only the sea deeper than the lowest level gets its procedural depth back. A shelf
                # cell the physics left above the sea (a lowstand deposit) is land, as in the physics' routing.
                restored = fields['sea'] & (fields['surface'] <= min(level for _, level in plan['sea_history']))
                if plan['keep_shelf'] == 'incision':
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
            if nest_info is not None:  # rule B's checks after the physics (nested_after_physics)
                nest_info['after_physics'] = nested_after_physics(domain, fields, elevation, sea, discharge, nest_info)
            coast = {**config['coast'], **band}
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
            micro = {**config.get('micro', {'enabled': False}), **band}  # a nested level: its band only
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
            if nest_info is not None and config['nesting'].get('project_after_physics', False):
                started_projection = time.perf_counter()
                elevation, nest_info['after_physics_projection'] = nested_project_after(domain, fields, elevation, sea,
                                                                                        discharge, config)
                timing['projection_after_physics'] = time.perf_counter() - started_projection
            np.save(directory / f'{domain.name}-elevation-physics.npy', crop(domain, elevation).astype(np.float32))
            if nesting_on(config):  # a nested child's parent sea
                np.save(directory / f'{domain.name}-sea.npy', crop(domain, sea))
            lakes = config.get('lakes', {'enabled': False})
            lake = np.zeros_like(sea)
            levels = None
            if lakes['enabled']:  # the closed depressions the physics left, before the filter (which spares them)
                levels = find_lakes(domain, elevation, sea, lakes)
                lake = np.isfinite(levels)
                timing['lakes_before_filter'] = int(ndimage.label(crop(domain, lake), structure=np.ones((3, 3)))[1])
            if spec['enabled'] and spec.get('when', 'before') == 'after':  # fine texture on the eroded surface
                x, y = fields['coords']
                field = filter_field(domain, fields, ice, spec, config['seed'])
                detail = erosion_pass(domain, elevation, x, y, fields['amplitude'], fields['erodibility'], spec,
                                      config['seed'], field=field)
                land = ~sea & ~lake  # lakes are water to the filter, as the sea: no texture on their beds
                if spec.get('mean_preserve', 'none') == 'smooth':
                    # Take out the filter's broad component (its mean over about a world cell, land only), so
                    # the averages over world cells stay close to the physics' (not exact: the pipeline's
                    # projection does that).
                    sigma = domain.k / 2.0
                    weight = ndimage.gaussian_filter(land.astype(float), sigma, mode='nearest')
                    broad = ndimage.gaussian_filter(np.where(land, detail, 0.0), sigma, mode='nearest') / \
                        np.maximum(weight, 1e-6)
                    detail = detail - broad
                elevation = np.where(land, np.maximum(elevation + detail, 1.0), elevation)
                if spec.get('fill_pits', False):
                    # filter.fill_pits (roadmap T1.6): the filter's lowering (its gullies and height_offset, both
                    # scaled by a strength that slope_fade raises with the slope) is deep on a slope and nil on the
                    # flatter floor below it, so the foot of the slope becomes a closed pit (G-coast: 7 % of the land,
                    # up to 1.2 km deep). They are filled as the relief's pits (sloped_fill: relief.fill_gradient with
                    # the fill noise) from the sea and the physics' lakes, which stand at their spill level before
                    # the filter.
                    water = sea | lake
                    outlets = elevation if levels is None else np.where(lake, levels, elevation)
                    filled = sloped_fill(domain, outlets, water, config['relief']['fill_gradient'], fields['fill_noise'])
                    raised = crop(domain, np.where(water, 0.0, filled - outlets))
                    timing['filter_pits_land_pct'] = round(100.0 * float((raised > 1.0).sum()) /
                                                           max(float((~crop(domain, water)).sum()), 1.0), 2)
                    timing['filter_pits_max_m'] = round(float(raised.max()), 1)
                    elevation = np.where(water, elevation, filled)
                if not domain.whole:
                    save_map(domain, directory / f'{domain.name}-filter-field.png',
                             colourize(crop(domain, field), FIELD_STOPS))
                    extra = [('filter strength field (black 0, white 1, orange 2 and more)', f'{domain.name}-filter-field.png')]
            np.save(directory / f'{domain.name}-elevation.npy', crop(domain, elevation).astype(np.float32))
            lake_surface = None
            if lakes['enabled']:  # the physics' lakes at their spill level on the final surface
                lake_surface = find_lakes(domain, elevation, sea, lakes, physics_lakes=lake)
                np.save(directory / f'{domain.name}-lake-surface.npy', crop(domain, lake_surface).astype(np.float32))
                timing.update(lake_statistics(domain, elevation, lake_surface))
            rivers = config.get('rivers', {'enabled': False})
            after = render(domain, elevation, sea, discharge, directory, ice=ice, lake_surface=lake_surface,
                           rivers=rivers if rivers['enabled'] else None)
            water = water_note(domain, lakes, rivers)
            images = [(f'elevation{water}', after['elevation']), (f'relief (hillshade, NW light){water}', after['relief']),
                      ('discharge (log, stretched from the median to the maximum)', after['discharge'])] + extra + \
                render_totals(domain, totals, sea, directory) + province_images + images
        else:
            images = [('elevation (procedural, no physics)', before['elevation'])] + province_images + images
        if nest_info is not None:
            (directory / f'{domain.name}-nest.json').write_text(json.dumps(nest_info, indent=1) + '\n', encoding='utf-8')
            start = nest_info['start_departure_m']
            timing['nest'] = {'parent': nest_info['parent'], 'factor': nest_info['factor'],
                              'start_departure_max_m': start['after_projection_unclipped']['max_abs'],
                              'start_block_check_after_refill': nest_info['start_block_check_m']['after_refill'],
                              'start_raised_to_floor': nest_info['start_land_raised_to_floor'],
                              'after_physics_projection': {
                                  key: value for key, value in nest_info.get('after_physics_projection', {}).items()
                                  if key != 'projection'},
                              'before_projection_max_m': start['before_projection']['max_abs'],
                              'projection': start['projection'],
                              'inflow_entries': nest_info['inflow']['placed'],
                              'inflow_m3_s': round(nest_info['inflow']['placed_m3_s_total'], 1),
                              'after_physics': {key: value for key, value in nest_info.get('after_physics', {}).items()
                                                if key != 'inflow_continuity'}}
            start_images = render(domain, fields['surface'], fields['sea'], None, directory, '-start')
            images.append(('start (rule B): the parent prolonged, the level\'s detail, the new pits filled, projected '
                           'onto the parent\'s cell means', start_images['relief']))
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
