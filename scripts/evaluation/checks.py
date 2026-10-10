"""Automatic pass/fail checks on run bundles, configured by evaluation profiles.

Each check returns {'check', 'status' ('pass' | 'fail' | 'skip' | 'reported'), 'summary', 'details'}.
A check with "report": true in the profile is 'reported' where it would fail: shown, not failing; a
check that only reports numbers (resolution-invariance) is always 'reported'.
A check is skipped, not failed, when its input is absent (for example no repeat run) or it does not
apply (a hierarchy check on a single pass; its summary says "not applicable").
"""
import fnmatch
import json
import math

import numpy as np

SPATIAL_LAYOUTS = {'grid', 'seasonal', 'monthly', 'levelseasonal'}


def result(check, status, summary, details=None):
    return {'check': check, 'status': status, 'summary': summary, 'details': details or []}


def select(bundle, patterns):
    """Field ids matching any glob pattern, in manifest order."""
    patterns = [patterns] if isinstance(patterns, str) else patterns
    return [identifier for identifier in bundle.fields
            if any(fnmatch.fnmatchcase(identifier, pattern) for pattern in patterns)]


def check_grid(bundle, spec, context):
    """Cell-centred W x W/2 grid; every spatial field ends in (rows, columns) of the
    surface grid or of a native model grid with the same 2:1 aspect."""
    grid = bundle.manifest['grid']
    problems = []
    if grid.get('registration') != 'cell-centred':
        problems.append(f"registration is {grid.get('registration')!r}")
    if bundle.columns != 2 * bundle.rows:
        problems.append(f'grid is {bundle.columns} x {bundle.rows}, not W x W/2')
    native = {}
    for identifier, info in bundle.fields.items():
        if info.layout not in SPATIAL_LAYOUTS:
            continue
        if len(info.shape) < 2:
            problems.append(f'{identifier}: shape {list(info.shape)} has no grid')
            continue
        rows, columns = info.shape[-2:]
        if (rows, columns) == (bundle.rows, bundle.columns):
            continue
        if columns == 2 * rows:
            native[identifier] = f'{columns} x {rows}'
        else:
            problems.append(f'{identifier}: {columns} x {rows} is not W x W/2')
    details = problems + [f'{identifier}: native grid {size}' for identifier, size in native.items()]
    summary = f'{bundle.columns} x {bundle.rows} cell-centred' if not problems else f'{len(problems)} problem(s)'
    return result('grid', 'fail' if problems else 'pass', summary, details)


def check_required(bundle, spec, context):
    missing = [identifier for identifier in spec['fields'] if identifier not in bundle.fields]
    missing += [identifier for identifier in spec.get('tables', []) if identifier not in bundle.tables]
    count = len(spec['fields']) + len(spec.get('tables', []))
    return result('required', 'fail' if missing else 'pass',
                  f'{count - len(missing)} of {count} present', [f'missing {m}' for m in missing])


def check_finite(bundle, spec, context):
    fields = [identifier for identifier in select(bundle, spec['fields'])
              if bundle.fields[identifier].dtype[1] == 'f']
    bad = []
    for identifier in fields:
        values = bundle.load(identifier)
        count = int(np.size(values) - np.count_nonzero(np.isfinite(values)))
        if count:
            bad.append(f'{identifier}: {count} non-finite value(s)')
    if not fields:
        return result('finite', 'skip', 'no floating-point fields selected')
    return result('finite', 'fail' if bad else 'pass', f'{len(fields) - len(bad)} of {len(fields)} fields finite', bad)


def mask_values(bundle, values, mask):
    """Restrict a spatial field to sea or land cells of the surface grid."""
    if mask is None:
        return values.ravel()
    sea = bundle.sea_fraction(*values.shape[-2:]) > 0.9
    cells = sea if mask == 'sea' else ~sea
    return values[..., cells].ravel()


def check_range(bundle, spec, context):
    identifier = spec['field']
    name = f"range {identifier}" + (f" ({spec['mask']})" if spec.get('mask') else '')
    if identifier not in bundle.fields:
        return result(name, 'skip', 'field absent')
    values = mask_values(bundle, bundle.load(identifier).astype(float), spec.get('mask'))
    values = values[np.isfinite(values)]
    if not values.size:
        return result(name, 'skip', 'no values')
    low, high = spec.get('min', -np.inf), spec.get('max', np.inf)
    # Absolute tolerance for float round-off at a physical bound (e.g. ice fraction 1).
    tolerance = spec.get('tolerance', 0)
    outside = int(np.count_nonzero((values < low - tolerance) | (values > high + tolerance)))
    summary = f'min {values.min():.6g}, max {values.max():.6g}; allowed [{low}, {high}]' + \
        (f' ± {tolerance:g}' if tolerance else '')
    details = [f'{outside} of {values.size} values outside'] if outside else []
    return result(name, 'fail' if outside else 'pass', summary, details)


def check_written(bundle, spec, context):
    """Listed fields were written by a step of this run, not left at their initial values."""
    unwritten = [identifier for identifier in select(bundle, spec['fields'])
                 if not bundle.fields[identifier].written_by]
    selected = select(bundle, spec['fields'])
    return result('written', 'fail' if unwritten else 'pass',
                  f'{len(selected) - len(unwritten)} of {len(selected)} written by a step',
                  [f'{identifier} not written' for identifier in unwritten])


def compare_hashes(bundle, other, identifiers):
    changed = []
    for identifier in identifiers:
        if identifier not in other.fields:
            changed.append(f'{identifier} absent from comparison bundle')
        elif bundle.fields[identifier].hash != other.fields[identifier].hash:
            changed.append(f'{identifier} differs')
    return changed


def check_frozen(bundle, spec, context):
    """Fields of frozen components equal those of the input bundle (imports, earlier stage)."""
    other = context.get('input')
    if other is None:
        return result('frozen', 'skip', 'no input bundle given')
    identifiers = [identifier for identifier in bundle.fields
                   if bundle.fields[identifier].component in spec['components']]
    identifiers = [identifier for identifier in identifiers if identifier not in spec.get('except', [])]
    changed = compare_hashes(bundle, other, identifiers)
    return result('frozen', 'fail' if changed else 'pass',
                  f"{len(identifiers) - len(changed)} of {len(identifiers)} fields of {', '.join(spec['components'])} unchanged",
                  changed)


def compare_bundles(bundle, other):
    """Differences between every field and table of two bundles, including those present in only one.
    No table is excluded: the bundle tables hold no timings or paths (repeat pairs of 2026-10-05 match)."""
    changed = compare_hashes(bundle, other, list(bundle.fields))
    changed += [f'{identifier} absent from the run bundle' for identifier in other.fields if identifier not in bundle.fields]
    for identifier, entry in bundle.tables.items():
        if identifier not in other.tables:
            changed.append(f'{identifier} absent from comparison bundle')
        elif entry['hash'] != other.tables[identifier]['hash']:
            changed.append(f'{identifier} differs')
    changed += [f'{identifier} absent from the run bundle' for identifier in other.tables if identifier not in bundle.tables]
    total = len(set(bundle.fields) | set(other.fields)) + len(set(bundle.tables) | set(other.tables))
    return changed, f'{total - len(changed)} of {total} fields and tables identical'


def check_repeat(bundle, spec, context):
    """A second run with the same seed, recipe and code produced identical fields and tables."""
    other = context.get('repeat')
    if other is None:
        return result('repeat', 'skip', 'no repeat bundle given')
    if other.seed != bundle.seed or other.recipe != bundle.recipe:
        return result('repeat', 'fail', 'repeat run used a different seed or recipe')
    changed, summary = compare_bundles(bundle, other)
    return result('repeat', 'fail' if changed else 'pass', summary, changed)


def check_recipe_equality(bundle, spec, context):
    """A run given every value of the module's recipe file as an option equals the run given none:
    the compiled defaults are the recipe (experiment.py --explicit-recipe makes the second run)."""
    other = context.get('explicit')
    if other is None:
        return result('recipe-equality', 'skip', 'no explicit-recipe bundle given')
    if other.seed != bundle.seed or other.recipe != bundle.recipe:
        return result('recipe-equality', 'fail', 'explicit-recipe run used a different seed or recipe')
    changed, summary = compare_bundles(bundle, other)
    return result('recipe-equality', 'fail' if changed else 'pass', summary, changed)


def json_table(bundle, identifier):
    return json.loads(bundle.table_text(identifier))


def lookup(document, key):
    """A value of a JSON document by dotted key (run_scenario.margin_processes)."""
    for part in key.split('.'):
        document = document[part]
    return document


def check_record(bundle, spec, context):
    """Values of a JSON table: dotted keys equal to given values ("equals") or to each other ("same")."""
    table = spec['table']
    if table not in bundle.tables:
        return result('record', 'skip', f'{table} absent')
    document = json_table(bundle, table)
    problems = []

    def value(key):
        """The value, or None after recording the problem: a missing key or a null value is never as expected."""
        try:
            found = lookup(document, key)
        except (KeyError, TypeError):
            problems.append(f'{key} missing')
            return None
        if found is None:
            problems.append(f'{key} is null')
        return found

    for key, expected in spec.get('equals', {}).items():
        found = value(key)
        if found is not None and found != expected:
            problems.append(f'{key} is {found!r}, expected {expected!r}')
    for left, right in spec.get('same', []):
        a, b = value(left), value(right)
        if a is not None and b is not None and a != b:
            problems.append(f'{left} is {a!r}, {right} is {b!r}')
    count = len(spec.get('equals', {})) + len(spec.get('same', []))
    return result('record', 'fail' if problems else 'pass', f'{max(count - len(problems), 0)} of {count} values of {table} as expected',
                  problems)


def check_budget(bundle, spec, context):
    """The column budget of the tectonic snapshot record closes: |relative mass residual| <= the limit."""
    table = spec.get('table', 'geosphere.tectonic_snapshot')
    if table not in bundle.tables:
        return result('budget', 'skip', f'{table} absent')
    budget = json_table(bundle, table).get('column_budget', {})
    if not budget.get('enabled'):
        return result('budget', 'fail', 'no column budget recorded')
    residual, limit = budget['relative_mass_residual'], spec['max_relative_residual']
    closed = np.isfinite(residual) and abs(residual) <= limit
    return result('budget', 'pass' if closed else 'fail', f'relative mass residual {residual:.3g}; allowed {limit:g}')


def check_minerals(bundle, spec, context):
    """The mineral budget counts as many deposits (inherited and formed) as the deposit table lists."""
    tables = ('geosphere.mineral_budget', 'geosphere.mineral_deposits')
    if any(table not in bundle.tables for table in tables):
        return result('minerals', 'skip', 'mineral tables absent')
    budget = json_table(bundle, tables[0])
    if budget.get('schema') != 'mineral-budget/v1':
        return result('minerals', 'fail', f"budget schema is {budget.get('schema')!r}")
    counted = budget['inherited']['count'] + sum(item['count'] for item in budget['formed'].values())
    listed = len(bundle.table_rows(tables[1]))
    return result('minerals', 'pass' if counted == listed else 'fail',
                  f'budget counts {counted} deposits; the table lists {listed}')


def absent_fields(bundle, identifiers):
    return [identifier for identifier in identifiers if identifier not in bundle.fields]


def whole_metres(values):
    """Round to the nearest whole metre, half away from zero (std::lround, the planet's wholemetres)."""
    return np.where(values >= 0, np.floor(values + 0.5), np.ceil(values - 0.5))


def land_mask(bundle, spec):
    """Land cells: elevation in whole metres above sea level (the rule of planet::basebelowsea)."""
    elevation = bundle.load(spec['elevation']).astype(float)
    return whole_metres(elevation) > float(bundle.scalar(spec['sea_level']))


OCEAN_FRACTION_PARAMETER = 'landscape.routing.ocean_minimum_area_fraction'
OCEAN_FRACTION_DEFAULT = 0.002


def ocean_mask(bundle, spec):
    """The routing's sea (erosion plan H3, landscape::oceanmask): every 8-connected region of
    below-sea cells (longitude wraps, no pole crossing) whose area is at least the manifest's
    landscape.routing.ocean_minimum_area_fraction (default 0.002) of the planet, and always the
    largest region by area (the lower cell index on a tie). Cell areas of the cell-centred grid:
    sin(north edge) - sin(south edge) of the row, the same for every column."""
    below = ~land_mask(bundle, spec)
    rows, columns = below.shape
    fraction = float(bundle.manifest.get('parameters', {}).get(OCEAN_FRACTION_PARAMETER, OCEAN_FRACTION_DEFAULT))
    edges = np.radians(90.0 - 180.0 * np.arange(rows + 1) / rows)
    rowarea = np.sin(edges[:-1]) - np.sin(edges[1:])
    region = np.full(below.shape, -1, dtype=np.int64)
    areas = []
    for start_y, start_x in zip(*np.nonzero(below)):
        if region[start_y, start_x] >= 0:
            continue
        ident = len(areas)
        region[start_y, start_x] = ident
        stack = [(int(start_y), int(start_x))]
        area = 0.0
        while stack:
            y, x = stack.pop()
            area += rowarea[y]
            for dy in (-1, 0, 1):
                ny = y + dy
                if ny < 0 or ny >= rows:
                    continue
                for dx in (-1, 0, 1):
                    nx = (x + dx) % columns
                    if below[ny, nx] and region[ny, nx] < 0:
                        region[ny, nx] = ident
                        stack.append((ny, nx))
        areas.append(area)
    if not areas:
        return np.zeros(below.shape, dtype=bool)
    areas = np.array(areas)
    sea = areas >= fraction * columns * float(rowarea.sum())
    sea[int(np.argmax(areas))] = True  # argmax takes the first of equal areas
    return np.where(region >= 0, sea[np.maximum(region, 0)], False)


def check_volume_budget(bundle, spec, context):
    """A volume total equals the sum of its parts: |total - sum(parts)| <= limit x max(|total|, 1e-12)."""
    missing = absent_fields(bundle, [spec['total']] + list(spec['parts']))
    if missing:
        return result('volume-budget', 'skip', f"{', '.join(missing)} absent")
    total = float(bundle.scalar(spec['total']))
    parts = sum(float(bundle.scalar(identifier)) for identifier in spec['parts'])
    residual = total - parts
    allowed = spec['max_relative_residual'] * max(abs(total), 1e-12)
    closed = bool(np.isfinite(residual)) and abs(residual) <= allowed
    return result('volume-budget', 'pass' if closed else 'fail',
                  f"{spec['total']} {total:.6g}, parts {parts:.6g}; residual {residual:.3g}, allowed {allowed:.3g}")


def check_drained(bundle, spec, context):
    """Every land cell (not the routing's sea, ocean_mask: an enclosed below-sea basin is land)
    has a flow direction or lies in a depression."""
    missing = absent_fields(bundle, [spec['direction'], spec['depression'], spec['elevation'], spec['sea_level']])
    if missing:
        return result('drained', 'skip', f"{', '.join(missing)} absent")
    land = ~ocean_mask(bundle, spec)
    direction = bundle.load(spec['direction'])
    depression = bundle.load(spec['depression'])
    undrained = int(np.count_nonzero(land & (direction == 0) & (depression == 0)))
    return result('drained', 'fail' if undrained else 'pass',
                  f'{undrained} of {int(np.count_nonzero(land))} land cells undrained')


def check_lake_above_terrain(bundle, spec, context):
    """Every lake cell (a nonzero lake surface) lies at or above its terrain: no lake surface below
    its bed. Compared in float32, as the fields are stored."""
    missing = absent_fields(bundle, [spec['lake'], spec['elevation']])
    if missing:
        return result('lake-above-terrain', 'skip', f"{', '.join(missing)} absent")
    surface = bundle.load(spec['lake']).astype(np.float32)
    elevation = bundle.load(spec['elevation']).astype(np.float32)
    lake = surface != 0
    below = lake & (surface < elevation)
    rows, columns = np.nonzero(below)
    depth = (elevation - surface)[below]
    order = np.argsort(-depth, kind='stable')[:10]
    details = [f'row {int(rows[i])}, column {int(columns[i])}: surface {float(surface[rows[i], columns[i]]):.2f}, '
               f'terrain {float(elevation[rows[i], columns[i]]):.2f}' for i in order]
    summary = f'{int(below.sum())} of {int(lake.sum())} lake cells below their terrain'
    if below.any():
        summary += f' (deepest {float(depth.max()):.2f} m)'
    return result('lake-above-terrain', 'fail' if below.any() else 'pass', summary, details)


# Flow direction codes as river_direction: (dx, dy) with east = +x, south = +y (north row first).
DIRECTION_OFFSETS = {1: (0, -1), 2: (1, -1), 3: (1, 0), 4: (1, 1), 5: (0, 1), 6: (-1, 1), 7: (-1, 0), 8: (-1, -1)}


def receiver_values(values, dx, dy):
    """values at the neighbour (x + dx, y + dy) of every cell; longitude wraps."""
    return np.roll(values, shift=(-dy, -dx), axis=(0, 1))


def has_neighbour(shape, dy):
    """False on the polar row that has no neighbour in direction dy."""
    rows = np.ones(shape, dtype=bool)
    if dy < 0:
        rows[0, :] = False
    elif dy > 0:
        rows[-1, :] = False
    return rows


def check_monotone(bundle, spec, context):
    """Discharge does not decrease from a land cell to its land receiver along flow_direction
    (land: not the routing's sea, ocean_mask)."""
    missing = absent_fields(bundle, [spec['direction'], spec['discharge'], spec['elevation'], spec['sea_level']])
    if missing:
        return result('monotone', 'skip', f"{', '.join(missing)} absent")
    land = ~ocean_mask(bundle, spec)
    direction = bundle.load(spec['direction'])
    discharge = bundle.load(spec['discharge']).astype(float)
    tolerance = 1e-6 * np.maximum(1.0, discharge)
    violations = 0
    links = 0
    for code, (dx, dy) in DIRECTION_OFFSETS.items():
        linked = land & (direction == code) & has_neighbour(land.shape, dy) & receiver_values(land, dx, dy)
        links += int(np.count_nonzero(linked))
        violations += int(np.count_nonzero(linked & (receiver_values(discharge, dx, dy) < discharge - tolerance)))
    return result('monotone', 'fail' if violations else 'pass',
                  f'{violations} of {links} land-to-land links lose discharge')


def check_regional_window(bundle, spec, context):
    """A regional bundle (world_pipeline region): its window block matches its grid (the kept
    global cells times the scale, the fine origin the first kept cell times the scale) on a
    W x W/2 world, every grid field has the bundle's grid, and its source block names the
    world bundle's manifest hash and the bundle's seed."""
    window, source = bundle.manifest.get('window'), bundle.manifest.get('source')
    if window is None or source is None:
        return result('regional-window', 'fail', 'no window or source block')
    problems = []
    scale, width, height = window['scale'], window['world_columns'], window['world_rows']
    if width != 2 * height or scale < 1:
        problems.append(f'world grid {width} x {height} at scale {scale}')
    if 'pixel_columns' in window:
        # A resampled bundle (--metres-per-pixel, multiscale plan S2.4): pixels of pixel_cells source-level cells
        # per edge from the kept box's north-west corner, partial edge pixels dropped.
        size = float(window.get('pixel_cells') or 0)
        expected = [math.floor(window[key] * scale / size * (1 + 1e-12)) if size > 0 else None for key in ('columns', 'rows')]
        if [window['pixel_columns'], window['pixel_rows']] != expected:
            problems.append(f"pixels {window['pixel_columns']} x {window['pixel_rows']} are not "
                            f"floor({window['columns']} x {window['rows']} cells at scale {scale} / {size:g} cells per pixel)")
        if (bundle.columns, bundle.rows) != (window['pixel_columns'], window['pixel_rows']):
            problems.append(f"grid {bundle.columns} x {bundle.rows} is not {window['pixel_columns']} x {window['pixel_rows']} pixels")
    elif (bundle.columns, bundle.rows) != (window['columns'] * scale, window['rows'] * scale):
        problems.append(f"grid {bundle.columns} x {bundle.rows} is not {window['columns']} x {window['rows']} cells at scale {scale}")
    if (window['fine_column0'], window['fine_row0']) != (window['column0'] * scale, window['row0'] * scale):
        problems.append('the fine origin is not the first kept cell times the scale')
    if not (0 <= window['column0'] < width and 0 <= window['row0'] and window['row0'] + window['rows'] <= height
            and window['columns'] + 2 * window['halo'] <= width):
        problems.append('the window or its halo lies off the world grid')
    for identifier, info in bundle.fields.items():
        if info.layout == 'grid' and tuple(info.shape) != (bundle.rows, bundle.columns):
            problems.append(f'{identifier}: shape {list(info.shape)}')
    digest = str(source.get('manifest_sha256', ''))
    if len(digest) != 64 or source.get('seed') != bundle.seed:
        problems.append('the source block lacks the manifest hash or names another seed')
    pixels = (f", {window['pixel_columns']} x {window['pixel_rows']} pixels of {float(window['pixel_cells']):g} cells"
              if 'pixel_columns' in window and window.get('pixel_cells') else '')
    summary = (f"{window['columns']} x {window['rows']} global cells at scale {scale} on {width} x {height}{pixels}"
               if not problems else f'{len(problems)} problem(s)')
    return result('regional-window', 'fail' if problems else 'pass', summary, problems)


# A stored error equal to its bound passes (a tie of children rounding the same way, multiscale plan
# section 6, S2.3); the slack only absorbs the decimal round trip of the manifest's numbers.
BOUND_SLACK = 1.0 + 2.0 ** -20


def hierarchy_levels(bundle):
    """The per-level diagnostics of a refinement hierarchy bundle (the manifest's diagnostics.levels,
    by level), or None for a single pass."""
    levels = (bundle.manifest.get('diagnostics') or {}).get('levels')
    return sorted(levels, key=lambda entry: entry['level']) if levels else None


def within_bound(value, bound):
    return bool(np.isfinite(value)) and float(value) <= float(bound) * BOUND_SLACK


def check_parent_mean_exact(bundle, spec, context):
    """Every level of a hierarchy bundle keeps its parents' values (multiscale plan S2.3): the level's
    largest parent-mean error is at most its own meanbound_m (half a float step of the kept box's
    largest |elevation|, read from the manifest, never recomputed), and its global-mean error at most
    the sum of the bounds of levels 1 to it. A single pass is not applicable."""
    name = 'parent-mean-exact'
    levels = hierarchy_levels(bundle)
    if levels is None:
        return result(name, 'skip', 'not applicable: a single pass (no diagnostics.levels)')
    problems, details, total, measured = [], [], 0.0, []
    overall = np.nan
    for entry in levels:
        level = entry.get('level')
        keys = ('maxparentmeanerror_m', 'meanbound_m', 'globalmeanerror_m')
        absent = [key for key in keys if not isinstance(entry.get(key), (int, float)) or isinstance(entry.get(key), bool)]
        if absent:
            problems.append(f"level {level}: {', '.join(absent)} missing or null")
            overall = np.nan
            continue
        error, bound, overall = (float(entry[key]) for key in keys)
        total += bound
        details.append(f'level {level}: parent-mean error {error:.6g} m, bound {bound:.6g} m; '
                       f'global-mean error {overall:.6g} m, bound {total:.6g} m')
        if not within_bound(error, bound):
            problems.append(f'level {level}: parent-mean error {error:.6g} m above its bound {bound:.6g} m')
        if not within_bound(overall, total):
            problems.append(f'level {level}: global-mean error {overall:.6g} m above the bounds\' sum {total:.6g} m')
        measured.append((error, bound))
    worst = max(measured, default=(np.nan, np.nan))
    summary = (f"{len(levels)} level(s): largest parent-mean error {worst[0]:.6g} m (bound {worst[1]:.6g} m); "
               f"final global-mean error {overall:.6g} m (bounds' sum {total:.6g} m)")
    return result(name, 'fail' if problems else 'pass', summary if not problems else f'{len(problems)} problem(s); ' + summary,
                  problems + details)


def box_average(values, factor):
    """Unweighted means of factor x factor blocks, in float64 (rows and columns divisible by factor)."""
    rows, columns = values.shape
    return values.astype(np.float64).reshape(rows // factor, factor, columns // factor, factor).mean(axis=(1, 3))


def shared_offsets(first, count, other_first, other_count, period=None):
    """Offsets (into the first window, into the other) of the global cells both windows keep; columns
    wrap with period (the world's width), rows do not."""
    cells = np.arange(count)
    if period is None:
        other = first + cells - other_first
    else:
        other = (first + cells - other_first) % period
    kept = (other >= 0) & (other < other_count)
    return cells[kept], other[kept]


def check_level_consistency(bundle, spec, context):
    """A coarser hierarchy bundle of the same world (context 'coarser', evaluate --coarser): the run's
    elevation box-averaged to the coarser level (factor^(L_run - L_coarser) cells per side, float64,
    unweighted) equals the coarser elevation on the global cells both keep, within the sum of the run's
    meanbound_m over levels L_coarser + 1 to L_run (the coarser run has no level between; its level
    L_coarser is the run's, bit for bit, multiscale plan S2.3). Without a coarser bundle, or when either
    is a single pass, not applicable."""
    name = 'level-consistency'
    coarser = context.get('coarser')
    if coarser is None:
        return result(name, 'skip', 'not applicable: no coarser bundle given (evaluate --coarser BUNDLE)')
    levels, coarse_levels = hierarchy_levels(bundle), hierarchy_levels(coarser)
    if levels is None or coarse_levels is None:
        which = 'the run' if levels is None else 'the coarser bundle'
        return result(name, 'skip', f'not applicable: {which} is a single pass (no diagnostics.levels)')
    window, coarse = bundle.manifest.get('window') or {}, coarser.manifest.get('window') or {}
    if 'pixel_columns' in window or 'pixel_columns' in coarse:
        which = 'the run' if 'pixel_columns' in window else 'the coarser bundle'
        return result(name, 'skip', f'not applicable: {which} is on a pixel grid (--metres-per-pixel), not on its level\'s cells')
    problems = []
    source, coarse_source = bundle.manifest.get('source') or {}, coarser.manifest.get('source') or {}
    if source.get('manifest_sha256') != coarse_source.get('manifest_sha256') or bundle.seed != coarser.seed:
        problems.append('the bundles come from different world bundles or seeds')
    if (window.get('world_columns'), window.get('world_rows')) != (coarse.get('world_columns'), coarse.get('world_rows')):
        problems.append('the world grids differ')
    if window.get('factor') != coarse.get('factor'):
        problems.append(f"factor {window.get('factor')} against {coarse.get('factor')}")
    if problems:
        return result(name, 'fail', 'not the same world and hierarchy', problems)
    level, coarse_level, factor = window['level'], coarse['level'], window['factor']
    if level <= coarse_level:
        return result(name, 'fail', f'the run (level {level}) is not finer than the coarser bundle (level {coarse_level})')
    ratio = factor ** (level - coarse_level)
    if window['scale'] != coarse['scale'] * ratio:
        return result(name, 'fail', f"scale {window['scale']} is not {coarse['scale']} x {factor}^{level - coarse_level}")
    bounds = {entry['level']: entry.get('meanbound_m') for entry in levels}
    missing = [step for step in range(coarse_level + 1, level + 1) if bounds.get(step) is None]
    if missing:
        return result(name, 'fail', f"no meanbound_m for level(s) {', '.join(map(str, missing))}")
    allowed = sum(float(bounds[step]) for step in range(coarse_level + 1, level + 1))
    coarse_columns, columns = shared_offsets(coarse['column0'], coarse['columns'], window['column0'], window['columns'],
                                             window['world_columns'])
    coarse_rows, rows = shared_offsets(coarse['row0'], coarse['rows'], window['row0'], window['rows'])
    if not columns.size or not rows.size:
        return result(name, 'fail', 'the windows share no global cell')
    k = coarse['scale']
    expand = lambda cells: (cells[:, None] * k + np.arange(k)[None, :]).ravel()
    identifier = spec.get('elevation', 'geosphere.base_elevation')
    fine_values, coarse_values = bundle.load(identifier), coarser.load(identifier)
    shapes = [f"{label}: {identifier} is {list(values.shape)}, not {[w['rows'] * w['scale'], w['columns'] * w['scale']]} "
              f"(the window's cells times the scale)"
              for label, values, w in (('the run', fine_values, window), ('the coarser bundle', coarse_values, coarse))
              if tuple(values.shape) != (w['rows'] * w['scale'], w['columns'] * w['scale'])]
    if shapes:
        return result(name, 'fail', 'grid does not match the window', shapes)
    averaged = box_average(fine_values, ratio)[np.ix_(expand(rows), expand(columns))]
    reference = coarse_values.astype(np.float64)[np.ix_(expand(coarse_rows), expand(coarse_columns))]
    difference = np.abs(averaged - reference)
    at = np.unravel_index(int(np.argmax(difference)), difference.shape)
    largest = float(difference[at])
    row, column = int(expand(coarse_rows)[at[0]]), int(expand(coarse_columns)[at[1]])
    world_row = coarse['row0'] + row // k
    world_column = (coarse['column0'] + column // k) % coarse['world_columns']
    above = int(np.count_nonzero(difference > allowed * BOUND_SLACK))
    summary = (f'largest difference of the box average from the coarser {largest:.6g} m at the coarser bundle\'s row {row}, column {column} '
               f'(global cell {world_column}, {world_row}); allowed {allowed:.6g} m (levels {coarse_level + 1} to {level}); '
               f'{rows.size} x {columns.size} global cells shared, factor {ratio}')
    details = [f'{above} of {difference.size} cells above the allowance'] if above else []
    return result(name, 'pass' if within_bound(largest, allowed) else 'fail', summary, details)


HYPSOMETRY_PERCENTILES = (5, 25, 50, 75, 95)


def landscape_summary(bundle, spec):
    """The land's hypsometry (percentiles of elevation above sea level) and channel share (channel
    cells over land cells; land: not geosphere.sea) of a regional bundle, over its whole grid."""
    sea_level = float(bundle.scalar(spec.get('sea_level', 'geosphere.sea_level')))
    elevation = bundle.load(spec.get('elevation', 'geosphere.base_elevation')).astype(np.float64) - sea_level
    land = bundle.load(spec.get('sea', 'geosphere.sea')) == 0
    channel = bundle.load(spec.get('channel', 'hydrosphere.channel_direction')) != 0
    count = int(np.count_nonzero(land))
    percentiles = np.percentile(elevation[land], HYPSOMETRY_PERCENTILES) if count else [None] * len(HYPSOMETRY_PERCENTILES)
    window = bundle.manifest.get('window') or {}
    return {'scale': window.get('scale'), 'level': window.get('level'), 'cells': int(land.size), 'land_cells': count,
            'hypsometry_m': {f'p{p}': None if value is None else float(value) for p, value in zip(HYPSOMETRY_PERCENTILES, percentiles)},
            'channel_share': float(np.count_nonzero(channel & land)) / count if count else None}


def check_resolution_invariance(bundle, spec, context):
    """Report only (status 'reported', never a failure): the land's hypsometry and channel share of the
    run and, with evaluate --coarser, of the coarser run beside it (multiscale plan section 6: the channel
    share and the valley spacing depend on the level until S3.2)."""
    name = 'resolution-invariance'
    identifiers = [spec.get('elevation', 'geosphere.base_elevation'), spec.get('sea_level', 'geosphere.sea_level'),
                   spec.get('sea', 'geosphere.sea'), spec.get('channel', 'hydrosphere.channel_direction')]
    missing = absent_fields(bundle, identifiers)
    if missing:
        return result(name, 'skip', f"{', '.join(missing)} absent")
    entries = {'run': landscape_summary(bundle, spec)}
    coarser = context.get('coarser')
    if coarser is not None and not absent_fields(coarser, identifiers):
        entries['coarser'] = landscape_summary(coarser, spec)
    label = lambda key, entry: f"{key} (scale {entry['scale']})"
    number = lambda value, text: '–' if value is None else text.format(value)
    details = [f"hypsometry of the land (m above sea level, {'/'.join(f'p{p}' for p in HYPSOMETRY_PERCENTILES)}): " +
               '; '.join(f"{label(key, entry)} " + '/'.join(number(value, '{:.0f}') for value in entry['hypsometry_m'].values())
                         for key, entry in entries.items()),
               'channel share (channel cells / land cells): ' +
               '; '.join(f"{label(key, entry)} {number(entry['channel_share'], '{:.4f}')} of {entry['land_cells']}"
                         for key, entry in entries.items())]
    summary = 'channel share ' + ', '.join(f"{label(key, entry)} {number(entry['channel_share'], '{:.4f}')}"
                                           for key, entry in entries.items()) + \
        '; median land height ' + ', '.join(f"{key} {number(entry['hypsometry_m']['p50'], '{:.0f}')} m" for key, entry in entries.items())
    outcome = result(name, 'reported', summary, details)
    outcome['values'] = entries
    return outcome


CHECKS = {'grid': check_grid, 'required': check_required, 'finite': check_finite, 'range': check_range,
          'written': check_written, 'frozen': check_frozen, 'repeat': check_repeat,
          'recipe-equality': check_recipe_equality, 'record': check_record, 'budget': check_budget,
          'minerals': check_minerals, 'volume-budget': check_volume_budget, 'drained': check_drained,
          'monotone': check_monotone, 'lake-above-terrain': check_lake_above_terrain,
          'regional-window': check_regional_window, 'parent-mean-exact': check_parent_mean_exact,
          'level-consistency': check_level_consistency, 'resolution-invariance': check_resolution_invariance}


def run_checks(bundle, specs, context):
    results = []
    for spec in specs:
        kind = spec['check']
        if kind not in CHECKS:
            raise ValueError(f'Unknown check {kind!r}; known: {", ".join(CHECKS)}')
        outcome = CHECKS[kind](bundle, spec, context)
        if spec.get('report') and outcome['status'] == 'fail':
            outcome['status'] = 'reported'
        results.append(outcome)
    return results
