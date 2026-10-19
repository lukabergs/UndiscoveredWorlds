"""A nested sandbox window (terrain roadmap T1.5, the 100 m test): a window finer than its parent window, started
from the parent's eroded surface, with the detail of its own cell size and a short physics. A prototype of one level
of the multiscale hierarchy (docs/architecture/multiscale-plan.md section 1), in the sandbox, not the pipeline.

Steps (the box in the parent's cropped pixels, the margin in parent pixels, the factor an integer):
  1. The parent's eroded surface: by default <parent>/<domain>-elevation-physics.npy (the physics' output with the
     sea restored, the coast term and the micro texture, before the erosion filter); --parent-surface final takes
     <domain>-elevation.npy (after the filter, what the parent draws). Not the final one by default: the filter is a
     render detail that the hierarchy replaces by each level's own detail and physics (multiscale plan section 1, Phase
     F subsumed), and it is not drainage-consistent: on t15-parent's G-coast it leaves closed pits on 12 % of the land
     around the fjord box, up to 834 m deep (p99 243 m; 1.3 % and p99 5 m before it). Its sea: cells at or below 0 m
     8-connected to a cell at or below 0 m on the crop's border (the parent saves no sea mask; inland basins below 0 m
     stay land).
  2. Prolongation by the factor: Catmull-Rom through the parent cell centres (rows and columns separably, the child at
     offset o lying (o + 0.5) / factor - 0.5 parent cells from its parent's centre), each child clamped to the range of
     its 2 x 2 bracketing parents (refinement_levels.hpp's prolong with prolongation_limit 1). Not mean-preserving; its
     parent-mean error is reported.
  3. The child's sea: the parent's sea re-flooded at the child's cells (cells at or below 0 m 8-connected to a child
     cell at or below 0 m under a parent sea cell).
  4. Detail of the child's cell size (nest.detail): fBm (terrain_sandbox.fbm, the sandbox's gradient noise and octave
     fade) from 2 parent cells (the parent grid's shortest wavelength) down to 2 child cells, with the relief's ridged
     share; drawn in local true km (east-west scaled by the cosine of the window's centre latitude) so it is isotropic
     on the ground. Amplitude on land: relief_fraction x the parent's local relief (max - min over
     relief_window_parent_cells^2 parent cells, sea at 0 m, prolonged) + floor_m; 0 on sea cells.
     Its mean over every parent cell is then removed: smooth_passes passes subtracting relaxation x the Catmull-Rom
     prolongation of the block means (over land children, so sea stays 0), then the exact remainder spread evenly over
     the land children: the detail's parent-cell means are 0 to rounding (reported). With normalise the detail is then
     scaled by one factor so its land rms is the amplitude's (the noise's mix and the mean removal leave about 0.3 of
     it), which keeps the parent-cell means at 0. The defaults are Earth's (scripts/sandbox/glo30_residual.py on the
     GLO-30 Sognefjord and Patagonia windows at these pixels, a calibration reference recorded in nested.json's
     detail.source): the 100 m surface minus its prolonged 611 m box means, less its land mean per box, has rms
     0.052 x the 3 x 3 local relief + 11.8 m (R^2 0.93 over 20 bins, both windows), and a spectrum slope (log power
     against log wavelength, 1.2 to 0.25 km) of 3.70, which this detail has at gain 0.52 with the ridged share 0.6.
  5. Sea again on the detailed surface (flooded from step 3's sea); the outer ring is sea (at most -1 m: the outlet) and
     the ring inside it is pinned at the prolonged parent (terrain_sandbox.exe holds pinned land at its start height).
  6. The rock (physics erodibility with the platform factor), uplift, precipitation proxy and lake-fill noise are the
     sandbox's own functions of position, evaluated on the child grid (terrain_sandbox.procedural with the filter and
     the depression fill off; its surface is not used).
  6b. With nest.fill_depressions, the start surface's closed depressions are filled (terrain_sandbox.sloped_fill, the
     relief's fill gradient and fill noise) except those holding a parent lake (the parent's lake surface, upsampled):
     the detail's pits and the parent's own (with --parent-surface final, the filter's: lakes on 10 % of the child's
     land in a first run) would otherwise stay as lakes, as the procedural pits did before relief.fill_depressions.
  6c. With nest.project_means, the whole start surface is projected onto the parent's cell means (the multiscale
     contract, docs/architecture/multiscale-plan.md section 1): the departure of every parent cell's child mean (sea
     children included, at their fixed values) moved over its land children outside the pinned ring, by the detail's
     smooth passes and exact remainder (project_means). It removes the prolongation's and the fill's departures; the
     start's parent-cell means are then the parent's to rounding wherever a parent cell has a movable child (reported
     before and after, and the physics' change after).
  7. terrain_sandbox.exe through terrain_sandbox.run_physics with the parent's physics values, except the scale-dependent
     ones (nest.physics): hillslope diffusivity (Earth's soil creep, not the grid parameter of 2.4 km cells), the
     landslide threshold slope, the step. Default: rivers and hillslopes at sea 0, no ice, for --duration Myr. With
     --glacial: the last glacial maximum (nest.pulse: ice at the glacial ELA, the sea at -120 m, glacial_ratio, for
     lgm_myr) then sea 0 without ice to --duration (the sandbox's history group, keep_shelf "incision").
  8. As the sandbox after its physics: the sea keeps its start depth (below the lowest level), troughs cut below the sea
     and joined to it become sea; no coast term, micro texture or filter (parent-level details); lakes at their spill
     level (Priority-Flood, nest.lakes).
Outputs in out/scratch/sandbox/<label>/ (the sandbox's names): <child>-elevation.npy (final bed), -elevation-procedural.npy
(the physics' input), -sea.npy, -lake-surface.npy, -discharge.npy (m3/s), the totals, -grid.json (the pixel geometry, for glo30_windows.py
--like and window_metrics.py), <child>-relief.png, -elevation.png, -discharge.png, the totals' maps; the parent's box
as <domain>-parent-*: the prolonged surface (the child's start without detail), the parent's sea and lakes (upsampled
by nearest), its rivers routed on the prolonged surface (a zero-duration physics call); nest.json (diagnostics,
timings) and index.html. The child is named <domain>-nest.

  uv run python scripts/sandbox/nested_window.py --parent t15-parent --domain G-coast --box 60,555,324,714
      --factor 6 --margin 16 [--duration 0.5] [--glacial] --label t15-nest-a [--set nest.detail.relief_fraction=0.1]
      [--exe PATH]
"""
import argparse
import copy
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS / 'sandbox'))
import terrain_sandbox as ts  # noqa: E402  (puts scripts/ and scripts/evaluation on the path)
from uwtools.bundle import Bundle  # noqa: E402

SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
DEFAULT_CONFIG = ROOT / 'configs' / 'sandbox' / 'nested.json'
# The grid operations and the projection live in nesting.py, shared with terrain_sandbox.py's nesting (rule B).
from nesting import (EIGHT, block_mean, catmull_rom_matrix, departure, flood, local_relief,  # noqa: E402,F401
                     nearest, project_means, prolong, remove_parent_means, ring)


# ---------------------------------------------------------------- main

def load_config(parent_directory, path, overrides):
    config = json.loads((parent_directory / 'config.json').read_text(encoding='utf-8'))
    config.update(json.loads(Path(path).read_text(encoding='utf-8')))
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


def prepare_world(config):
    """The sandbox's world fields, as terrain_sandbox.main prepares them."""
    bundle = Bundle.open(ROOT / config['base_bundle'])
    world = ts.world_fields(bundle)
    radius_m = float(bundle.manifest.get('parameters', {}).get('world.body.radius_km', 6371.0)) * 1000.0
    world['precipitation'] = ts.precipitation(world, config['precipitation'], radius_m)
    if config['seafloor'].get('abyss', {}).get('enabled', False):
        world['elevation'] = ts.damp_abyss(world, config['seafloor']['abyss'])
    world['ocean'] = ts.ocean_world(world, config['seafloor'], config['seed'], radius_m / 1000.0)
    if config.get('platform', {}).get('enabled', False):
        world['provinces'] = ts.provinces_world(world, config['platform'], 2.0 * math.pi * radius_m / 1000.0 / bundle.columns)
    return bundle, world


def physics_config(config, duration, glacial):
    """The parent's config with the nest's scale-dependent physics and the run's history."""
    nest = config['nest']
    run = copy.deepcopy(config)
    window = run['physics']['window']
    window.update(duration_myr=duration, step_myr=nest['physics']['step_myr'],
                  threshold_slope=nest['physics']['threshold_slope'])
    run['physics']['diffusivity_km2_myr'] = nest['physics']['diffusivity_km2_myr']
    run['sea_level'] = {**run.get('sea_level', {}), 'enabled': False}
    if glacial:
        pulse = nest['pulse']
        run['glacial']['enabled'] = True
        run['history'] = {'enabled': True, 'keep_shelf': 'incision', 'phases': [
            {'name': 'last glacial maximum', 'start_myr': 0.0, 'sea_m': pulse['sea_m'], 'ice': True,
             'glacial_ratio': pulse['glacial_ratio'], 'step_myr': pulse['lgm_step_myr'], 'no_infill': False},
            {'name': 'postglacial: rivers only', 'start_myr': pulse['lgm_myr'], 'sea_m': 0, 'ice': False,
             'glacial_ratio': pulse['glacial_ratio'], 'step_myr': nest['physics']['step_myr'], 'no_infill': False}],
            'end': {'time_myr': duration, 'sea_m': 0}}
    else:
        run['glacial']['enabled'] = False
        run['history'] = {**run.get('history', {}), 'enabled': False}
    return run


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--parent', required=True, help='the parent run label under out/scratch/sandbox')
    parser.add_argument('--domain', required=True)
    parser.add_argument('--box', required=True, help="c0,r0,c1,r1: the parent's cropped pixels, inclusive")
    parser.add_argument('--factor', type=int, default=6)
    parser.add_argument('--margin', type=int, default=16, help='parent pixels around the box')
    parser.add_argument('--duration', type=float, default=0.5, help='Myr')
    parser.add_argument('--glacial', action='store_true', help='a last-glacial-maximum pulse first (nest.pulse)')
    parser.add_argument('--label', required=True)
    parser.add_argument('--config', default=str(DEFAULT_CONFIG))
    parser.add_argument('--set', action='append', default=[], metavar='GROUP.KEY=JSON')
    parser.add_argument('--exe', default=str(ts.PHYSICS), help='terrain_sandbox.exe to run')
    parser.add_argument('--parent-surface', choices=('physics', 'final'), default='physics',
                        help="the parent's surface before the erosion filter (default) or after it")
    arguments = parser.parse_args()
    started_all = time.perf_counter()
    timings = {}

    parent_directory = SANDBOX / arguments.parent
    config = load_config(parent_directory, arguments.config, arguments.set)
    nest = config['nest']
    ts.PHYSICS = Path(arguments.exe)
    directory = SANDBOX / arguments.label
    directory.mkdir(parents=True, exist_ok=True)
    f, margin = arguments.factor, arguments.margin
    c0, r0, c1, r1 = (int(value) for value in arguments.box.split(','))

    started = time.perf_counter()
    bundle, world = prepare_world(config)
    parent_domain = ts.make_domains(config, bundle, [arguments.domain])[0]
    timings['world'] = time.perf_counter() - started
    suffix = '-elevation-physics' if arguments.parent_surface == 'physics' else '-elevation'
    parent = np.load(parent_directory / f'{arguments.domain}{suffix}.npy').astype(np.float64)
    parent_lakes = np.load(parent_directory / f'{arguments.domain}-lake-surface.npy')
    if not (0 <= c0 - margin and c1 + margin < parent.shape[1] and 0 <= r0 - margin and r1 + margin < parent.shape[0]):
        raise SystemExit(f'box with margin outside the parent crop {parent.shape}')
    border = ring(parent.shape, 1)
    parent_sea = flood(parent <= 0.0, border)

    # The child domain: the box plus the margin, on the parent's grid lines.
    pr, pc = r0 - margin, c0 - margin  # the child domain's first parent cell (crop indices)
    prow, pcol = r1 - r0 + 1 + 2 * margin, c1 - c0 + 1 + 2 * margin
    k = parent_domain.k * f
    crop_row0, crop_col0 = parent_domain.crop[0], parent_domain.crop[2]
    x0 = parent_domain.x[crop_col0 + pc] - 0.5 / parent_domain.k
    y0 = parent_domain.y[crop_row0 + pr] - 0.5 / parent_domain.k
    rows, columns = prow * f, pcol * f
    keep = (margin * f, rows - margin * f, margin * f, columns - margin * f)
    child = ts.Domain(f'{arguments.domain}-nest', columns, rows, x0, y0, k, parent_domain.world_columns,
                      parent_domain.world_rows, parent_domain.radius_km, crop=keep, exaggeration=nest['exaggeration'],
                      physics=copy.deepcopy(config['physics']['window']))
    latitude_c = child.north_deg - rows * child.pixel_deg / 2.0
    pixel_ns_km = child.fine_km
    pixel_ew_km = child.fine_km * math.cos(math.radians(latitude_c))

    # 2. Prolongation; 3. the child's sea.
    started = time.perf_counter()
    prolonged = prolong(parent, f, pr, pc, prow, pcol)
    error = block_mean(prolonged, f) - parent[pr:pr + prow, pc:pc + pcol]
    inner = (slice(margin, prow - margin), slice(margin, pcol - margin))
    under_sea = nearest(parent_sea, f, pr, pc, prow, pcol)
    sea0 = flood(prolonged <= 0.0, under_sea & (prolonged <= 0.0))
    timings['prolong'] = time.perf_counter() - started

    # 4. Detail.
    started = time.perf_counter()
    spec = nest['detail']
    relief_parent = local_relief(parent, int(spec['relief_window_parent_cells']))
    relief = np.maximum(prolong(relief_parent, f, pr, pc, prow, pcol), 0.0)
    amplitude = spec['relief_fraction'] * relief + spec['floor_m']
    band = {'wavelength_km': spec['longest_parent_cells'] * parent_domain.fine_km, 'octaves': spec['octaves'],
            'lacunarity': spec['lacunarity'], 'gain': spec['gain'], 'min_cells': spec['min_cells']}
    # At the cells' unit vectors in km (the sandbox's noise on the sphere, sphere-mesh.md section 11).
    noise = ts.fbm(child.points(), band, config['seed'] + spec['seed_offset'], child)
    noise = (1 - spec['ridged']) * noise + spec['ridged'] * ts.ridged(noise) / 0.6
    land0 = ~sea0
    raw = amplitude * noise * land0
    detail = remove_parent_means(raw, land0, f, spec['smooth_passes'], spec['relaxation'])
    detail_scale = 1.0
    if spec.get('normalise', False):  # the detail's land rms made the amplitude's (docstring, step 4)
        detail_scale = float(np.sqrt((amplitude[land0] ** 2).mean() / (detail[land0] ** 2).mean()))
        detail *= detail_scale
    detail_error = block_mean(detail, f)
    surface = prolonged + detail
    timings['detail'] = time.perf_counter() - started

    # 5. Sea, outlet ring, pinned ring.
    sea = flood(surface <= 0.0, sea0)
    outer = ring(surface.shape, 1)
    pinned = ring(surface.shape, 2) & ~outer & ~sea
    surface = np.where(pinned, prolonged, surface)
    surface = np.where(outer, np.minimum(prolonged, -1.0), surface)
    sea |= outer

    # 6. The rock, uplift, precipitation and fill noise on the child grid.
    started = time.perf_counter()
    procedural_config = copy.deepcopy(config)
    procedural_config['filter']['enabled'] = False
    procedural_config['relief']['fill_depressions'] = False
    fields, _ = ts.procedural(child, world, procedural_config, config['seed'])
    timings['rock_and_rain'] = time.perf_counter() - started

    # 6b. The inherited and the detail's pits filled (sloped, as relief.fill_depressions) except the depressions that
    # hold a parent lake (docstring, step 6b).
    fill_stats = {}
    if nest.get('fill_depressions', True):
        started = time.perf_counter()
        filled = ts.sloped_fill(child, surface, sea, config['relief']['fill_gradient'], fields['fill_noise'])
        labels, count = ndimage.label(filled - surface > 1e-3, structure=EIGHT)
        under_lake = np.isfinite(nearest(parent_lakes, f, pr, pc, prow, pcol))
        kept = np.unique(labels[under_lake & (labels > 0)])
        kept_mask = np.isin(labels, kept[kept > 0])
        raise_m = np.where(kept_mask, 0.0, filled - surface)
        surface = np.where(kept_mask, surface, filled)
        surface = np.where(pinned, prolonged, surface)
        inside = ts.crop(child, ~sea)
        fill_stats = {'depressions': int(count), 'kept_with_parent_lake': int(np.count_nonzero(kept)),
                      'filled_land_share': float((ts.crop(child, raise_m) > 1e-3)[inside].mean()),
                      'fill_mean_m': float(ts.crop(child, raise_m)[ts.crop(child, raise_m) > 1e-3].mean()),
                      'fill_max_m': float(raise_m.max())}
        timings['fill'] = time.perf_counter() - started

    # 6c. The whole start projected onto the parent's cell means (docstring, step 6c).
    targets = parent[pr:pr + prow, pc:pc + pcol]
    movable = ~sea & ~pinned
    start_departure = {'before_projection': departure(surface, targets, movable, f, inner)}
    if nest.get('project_means', False):
        surface = project_means(surface, targets, movable, f, spec['smooth_passes'], spec['relaxation'])
        start_departure['after_projection'] = departure(surface, targets, movable, f, inner)
        start_departure['land_at_or_below_0_m'] = int(np.count_nonzero(movable & (surface <= 0.0)))
    physics_fields = {'surface': surface, 'sea': sea, 'erodibility': fields['erodibility'],
                      'physics_erodibility': fields['physics_erodibility'], 'uplift': fields['uplift'],
                      'precipitation': fields['precipitation'], 'pinned': pinned, 'fill_noise': fields['fill_noise']}
    del fields

    # 7. Physics.
    run = physics_config(config, arguments.duration, arguments.glacial)
    child.physics = run['physics']['window']  # run_physics reads the domain's group over the physics group
    started = time.perf_counter()
    elevation, discharge, ice, totals, report = ts.run_physics(child, physics_fields, run, directory)
    timings['physics'] = time.perf_counter() - started

    # 8. After the physics, as the sandbox.
    started = time.perf_counter()
    plan = ts.schedule(run)
    restored = sea
    if plan and plan.get('keep_shelf', False):
        restored = sea & (surface <= min(level for _, level in plan['sea_history']))
        if plan['keep_shelf'] == 'incision':
            shelf = sea & ~restored
            elevation = np.where(shelf, np.minimum(elevation, surface), elevation)
    elevation = np.where(restored, surface, elevation)
    final_sea = restored | flood(elevation <= 0.0, sea)
    lake_surface = ts.find_lakes(child, elevation, final_sea, nest['lakes'])
    lake_stats = ts.lake_statistics(child, elevation, lake_surface)
    timings['lakes'] = time.perf_counter() - started

    started = time.perf_counter()
    crop = lambda values: ts.crop(child, values)  # noqa: E731
    np.save(directory / f'{child.name}-elevation-procedural.npy', crop(surface).astype(np.float32))
    np.save(directory / f'{child.name}-elevation.npy', crop(elevation).astype(np.float32))
    np.save(directory / f'{child.name}-sea.npy', crop(final_sea))
    np.save(directory / f'{child.name}-lake-surface.npy', crop(lake_surface).astype(np.float32))
    np.save(directory / f'{child.name}-discharge.npy', crop(discharge).astype(np.float32))  # m3/s
    for name, values in totals.items():
        np.save(directory / f'{child.name}-{name}.npy', crop(values).astype(np.float32))
    grid = {'rows': keep[1] - keep[0], 'columns': keep[3] - keep[2], 'pixel_ns_km': pixel_ns_km,
            'pixel_ew_km': pixel_ew_km, 'centre_latitude': latitude_c, 'pixel_deg': child.pixel_deg,
            'north_deg': child.north_deg - keep[0] * child.pixel_deg, 'exaggeration': nest['exaggeration'],
            'note': 'the kept window; east-west km at the centre latitude'}
    for name in (child.name, f'{arguments.domain}-parent'):
        (directory / f'{name}-grid.json').write_text(json.dumps(grid, indent=1) + '\n', encoding='utf-8')
    rivers = config.get('rivers', {'enabled': False})
    water = ts.water_note(child, {**nest['lakes'], 'enabled': True}, rivers)
    after = ts.render(child, elevation, final_sea, discharge, directory, ice=ice, lake_surface=lake_surface,
                      rivers=rivers if rivers['enabled'] else None)
    images = [(f'child elevation{water}', after['elevation']), (f'child relief (hillshade, NW light){water}', after['relief']),
              ('child discharge (log, stretched from the median to the maximum)', after['discharge'])]
    images += ts.render_totals(child, totals, final_sea, directory)
    start_names = ts.render(child, surface, sea, None, directory, suffix='-procedural')
    images.append(('child start: the prolonged parent plus the detail, before the physics', start_names['relief']))

    # The parent's box at the child's pixels: the prolonged surface, the parent's sea and lakes by nearest, its rivers
    # routed on the prolonged surface by a zero-duration physics call.
    parent_view = ts.Domain(f'{arguments.domain}-parent', columns, rows, x0, y0, k, child.world_columns,
                            child.world_rows, child.radius_km, crop=keep, exaggeration=nest['exaggeration'],
                            physics=None)
    view_sea = nearest(parent_sea, f, pr, pc, prow, pcol)
    view_lakes = nearest(parent_lakes, f, pr, pc, prow, pcol)
    route = physics_config(config, 0.0, False)
    parent_view.physics = route['physics']['window']
    view_fields = dict(physics_fields, surface=prolonged, sea=view_sea | outer, pinned=np.zeros_like(pinned))
    _, view_discharge, _, _, _ = ts.run_physics(parent_view, view_fields, route, directory)
    view = ts.render(parent_view, prolonged, view_sea, view_discharge, directory, lake_surface=view_lakes,
                     rivers=rivers if rivers['enabled'] else None)
    np.save(directory / f'{parent_view.name}-elevation.npy', crop(prolonged).astype(np.float32))
    np.save(directory / f'{parent_view.name}-sea.npy', crop(view_sea))
    np.save(directory / f'{parent_view.name}-lake-surface.npy', crop(view_lakes).astype(np.float32))
    np.save(directory / f'{parent_view.name}-discharge.npy', crop(view_discharge).astype(np.float32))
    images.append((f'parent {arguments.parent} box, prolonged (Catmull-Rom, bracket-clamped) to the child pixels; the '
                   'parent\'s sea and lakes; rivers routed on the prolonged surface', view['relief']))
    images.append(('parent box elevation', view['elevation']))
    timings['render'] = time.perf_counter() - started

    # Diagnostics: the prolongation's and the detail's parent-mean errors, the physics' change of the parent means.
    change = block_mean(elevation, f) - parent[pr:pr + prow, pc:pc + pcol]
    land_blocks = block_mean((~final_sea).astype(float), f)[inner] > 0.5
    diagnostics = {
        'parent': arguments.parent, 'parent_surface': arguments.parent_surface, 'domain': arguments.domain, 'box': [c0, r0, c1, r1], 'factor': f, 'margin': margin,
        'child_cells': [rows, columns], 'kept': grid, 'glacial': arguments.glacial, 'duration_myr': arguments.duration,
        'prolongation_parent_mean_error_m': {'max_abs': float(np.abs(error[inner]).max()),
                                             'rms': float(np.sqrt((error[inner] ** 2).mean()))},
        'detail_parent_mean_max_abs_m': float(np.abs(detail_error).max()),
        'detail_rms_m_land': float(np.sqrt((detail[land0] ** 2).mean())), 'detail_rms_scale': detail_scale,
        'start_parent_mean_departure_m': start_departure,
        'detail_amplitude_m_land': {'p50': float(np.percentile(amplitude[land0], 50)),
                                    'p90': float(np.percentile(amplitude[land0], 90))},
        'physics_parent_mean_change_m_land': {
            'mean': float(change[inner][land_blocks].mean()), 'p5': float(np.percentile(change[inner][land_blocks], 5)),
            'p95': float(np.percentile(change[inner][land_blocks], 95)),
            'max_abs': float(np.abs(change[inner][land_blocks]).max())},
        'physics_report': report, 'lakes': lake_stats, 'start_fill': fill_stats,
        'physics_command_values': {'k': run['physics']['k'], 'm': run['physics']['m'],
                                   'deposition': run['physics']['deposition'],
                                   'diffusivity_km2_myr': run['physics']['diffusivity_km2_myr'],
                                   'window': run['physics']['window'],
                                   'history': run['history'] if run['history'].get('enabled') else None},
        'nest': nest, 'exe': str(arguments.exe)}
    timings['total'] = time.perf_counter() - started_all
    diagnostics['timings_s'] = {key: round(value, 1) for key, value in timings.items()}
    (directory / 'nest.json').write_text(json.dumps(diagnostics, indent=1) + '\n', encoding='utf-8')
    (directory / 'config.json').write_text(json.dumps(run, indent=1), encoding='utf-8')
    note = (f'Nested window {child.name}: parent {arguments.parent} ({arguments.domain}, k {parent_domain.k}, '
            f'{parent_domain.fine_km:.3f} km) box {arguments.box} refined {f} times (k {k}, {child.fine_km * 1000:.0f} m '
            f'north-south, {pixel_ew_km * 1000:.0f} m east-west at {latitude_c:.2f} N), margin {margin} parent pixels; '
            f'physics {arguments.duration} Myr{" with a glacial pulse" if arguments.glacial else " at sea 0, no ice"}. '
            'Written by scripts/sandbox/nested_window.py.')
    html = ''.join(f'<figure><a href="{path}"><img src="{path}"></a><figcaption>{caption}</figcaption></figure>'
                   for caption, path in images)
    (directory / 'index.html').write_text(ts.PAGE.format(
        label=arguments.label, note=note, sections=f'<h2>{child.name}</h2><div class="row">{html}</div>',
        config=json.dumps(run, indent=1), timings=json.dumps(diagnostics, indent=1)), encoding='utf-8')
    print(json.dumps({key: diagnostics[key] for key in ('prolongation_parent_mean_error_m', 'detail_parent_mean_max_abs_m',
                                                         'detail_rms_m_land', 'detail_rms_scale', 'start_fill',
                                                         'start_parent_mean_departure_m', 'physics_parent_mean_change_m_land',
                                                         'physics_report', 'lakes', 'timings_s')}, indent=1))
    print(directory / 'index.html')


if __name__ == '__main__':
    main()
