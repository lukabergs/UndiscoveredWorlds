"""Automatic pass/fail checks on run bundles, configured by evaluation profiles.

Each check returns {'check', 'status' ('pass' | 'fail' | 'skip' | 'reported'), 'summary', 'details'}.
A check with "report": true in the profile is 'reported' where it would fail: shown, not failing.
A check is skipped, not failed, when its input is absent (for example no repeat run).
"""
import fnmatch

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


def check_repeat(bundle, spec, context):
    """A second run with the same seed, recipe and code produced identical fields."""
    other = context.get('repeat')
    if other is None:
        return result('repeat', 'skip', 'no repeat bundle given')
    if other.seed != bundle.seed or other.recipe != bundle.recipe:
        return result('repeat', 'fail', 'repeat run used a different seed or recipe')
    changed = compare_hashes(bundle, other, list(bundle.fields))
    return result('repeat', 'fail' if changed else 'pass',
                  f'{len(bundle.fields) - len(changed)} of {len(bundle.fields)} fields identical', changed)


CHECKS = {'grid': check_grid, 'required': check_required, 'finite': check_finite, 'range': check_range,
          'written': check_written, 'frozen': check_frozen, 'repeat': check_repeat}


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
