"""Automatic pass/fail checks on run bundles, configured by evaluation profiles.

Each check returns {'check', 'status' ('pass' | 'fail' | 'skip' | 'reported'), 'summary', 'details'}.
A check with "report": true in the profile is 'reported' where it would fail: shown, not failing.
A check is skipped, not failed, when its input is absent (for example no repeat run).
"""
import fnmatch
import json

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
    """Every land cell has a flow direction or lies in a depression."""
    missing = absent_fields(bundle, [spec['direction'], spec['depression'], spec['elevation'], spec['sea_level']])
    if missing:
        return result('drained', 'skip', f"{', '.join(missing)} absent")
    land = land_mask(bundle, spec)
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
    """Discharge does not decrease from a land cell to its land receiver along flow_direction."""
    missing = absent_fields(bundle, [spec['direction'], spec['discharge'], spec['elevation'], spec['sea_level']])
    if missing:
        return result('monotone', 'skip', f"{', '.join(missing)} absent")
    land = land_mask(bundle, spec)
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


CHECKS = {'grid': check_grid, 'required': check_required, 'finite': check_finite, 'range': check_range,
          'written': check_written, 'frozen': check_frozen, 'repeat': check_repeat,
          'recipe-equality': check_recipe_equality, 'record': check_record, 'budget': check_budget,
          'minerals': check_minerals, 'volume-budget': check_volume_budget, 'drained': check_drained,
          'monotone': check_monotone, 'lake-above-terrain': check_lake_above_terrain}


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
