"""Compare two run bundles field by field.

    uv run --offline scripts/evaluation/compare_bundles.py BUNDLE_A BUNDLE_B [--fields PREFIX ...] [--json OUT]

A bundle is a directory written by `world_pipeline --bundle` (schema uw-run-bundle/v2, or v1
read as grid kind latlon; docs/architecture/sphere-mesh.md section 7). The native fields are
compared exactly; for differing fields the report gives the number of differing values and the
largest absolute difference. Tables are compared as text. The lat-lon exports and the cell table
of cube bundles are compared too (an entry on one side only is missing). Of the grid block the
kind (v1 is latlon), face_cells and export_columns must agree; the schema, the radius and the
layout strings are not compared, so a v1 and a v2 lat-lon bundle with equal fields are identical.
Exit status: 0 when the grid agrees and every compared field, export, cell array and table is
identical, 1 otherwise.

--legacy-encoding compares a bundle from before the float migration of the climate
fields and base elevation (integers, some scaled) with a later one: a float field is put
back into the former encoding (scaled, rounded, clipped to the integer type) and the
largest difference in physical units is reported as lost_precision.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import SCHEMAS, read_manifest

# Manifest lists compared entry by entry: native fields, lat-lon exports, the cell table, tables.
KINDS = ('fields', 'exports', 'cells', 'tables')
# Grid block members that must agree (the radius and the v1/v2 layout strings are not compared).
GRID_KEYS = ('kind', 'face_cells', 'export_columns')


def grid_value(manifest, key):
    grid = manifest.get('grid', {})
    return grid.get(key, 'latlon') if key == 'kind' else grid.get(key)


def load_manifest(bundle):
    manifest = read_manifest(bundle)
    if manifest.get('schema') not in SCHEMAS:
        raise SystemExit(f'{bundle}: unsupported schema {manifest.get("schema")!r}')
    return manifest


# Factors of the former int16 encodings of climate fields (climatelegacyscale in climate_world.hpp).
LEGACY_SCALES = {'atmosphere.vertical_velocity': 100, 'hydrosphere.ocean_current_east': 100,
                 'hydrosphere.ocean_current_south': 100, 'atmosphere.maritime_influence': 1000,
                 'atmosphere.maritime_thermal_anomaly': 10, 'atmosphere.convergence': 100,
                 'atmosphere.orographic_uplift': 4, 'atmosphere.orographic_subsidence': 4}


def legacy_encoded(values, identifier, integer_type):
    scale = LEGACY_SCALES.get(identifier, 1)
    limits = np.iinfo(integer_type)
    # The product is taken in float32 as the accessors do; half away from zero, as std::lround.
    scaled = np.clip(values.astype(np.float32) * np.float32(scale), limits.min, limits.max).astype(np.float64)
    return (np.sign(scaled) * np.floor(np.abs(scaled) + 0.5)).astype(np.int64)


def selected(identifier, prefixes):
    return not prefixes or any(identifier.startswith(prefix) for prefix in prefixes)


def compare_field(left, right, entry_left, entry_right, legacy=False):
    a = np.load(left / entry_left['file'])
    b = np.load(right / entry_right['file'])
    if legacy and a.shape == b.shape and a.dtype.kind != b.dtype.kind and {a.dtype.kind, b.dtype.kind} == {'i', 'f'}:
        integers, floats = (a, b) if a.dtype.kind == 'i' else (b, a)
        encoded = legacy_encoded(floats, entry_left['id'], integers.dtype)
        differing = int(np.count_nonzero(encoded != integers))
        scale = LEGACY_SCALES.get(entry_left['id'], 1)
        lost = float(np.max(np.abs(floats.astype(np.float64) - integers / scale))) if floats.size else 0.0
        if differing == 0:
            return {'status': 'identical', 'legacy': True, 'lost_precision': lost}
        first = np.unravel_index(int(np.argmax(encoded != integers)), integers.shape)
        return {'status': 'different', 'differing': differing, 'total': int(integers.size), 'legacy': True,
                'max_abs_difference': float(np.max(np.abs(encoded - integers))), 'lost_precision': lost,
                'first_index': [int(i) for i in first]}
    if a.dtype != b.dtype or a.shape != b.shape:
        return {'status': 'shape', 'left': f'{a.dtype}{list(a.shape)}', 'right': f'{b.dtype}{list(b.shape)}'}
    if a.dtype.kind == 'f':
        same = (a == b) | (np.isnan(a) & np.isnan(b))
    else:
        same = a == b
    differing = int(np.size(same) - np.count_nonzero(same))
    if differing == 0:
        return {'status': 'identical'}
    diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
    diff = np.where(same, 0.0, diff)
    first = np.unravel_index(int(np.argmax(~same)), a.shape) if a.ndim else ()
    return {'status': 'different', 'differing': differing, 'total': int(np.size(a)),
            'max_abs_difference': float(np.nanmax(diff)), 'first_index': [int(i) for i in first]}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    parser.add_argument('--fields', nargs='*', default=[], help='Only compare ids starting with these prefixes')
    parser.add_argument('--json', type=Path, help='Write the full comparison as JSON')
    parser.add_argument('--legacy-encoding', action='store_true',
                        help='Compare float fields with their former int16 encoding in the other bundle')
    args = parser.parse_args()

    left, right = load_manifest(args.left), load_manifest(args.right)
    result = {'left': str(args.left), 'right': str(args.right), 'fields': {}, 'tables': {}}
    # The grid: kind (v1 is latlon), face cells and export columns must agree.
    result['grid'] = {}
    for key in GRID_KEYS:
        values = [grid_value(manifest, key) for manifest in (left, right)]
        if values[0] != values[1]:
            result['grid'][key] = {'status': 'different', 'left': values[0], 'right': values[1]}
    for kind in KINDS:
        if kind in ('exports', 'cells'):
            # Only cube bundles have them; an entry on one side only is missing on the other.
            if kind not in left and kind not in right:
                continue
            result[kind] = {}
        lefts = {entry['id']: entry for entry in left.get(kind, [])}
        rights = {entry['id']: entry for entry in right.get(kind, [])}
        for identifier in sorted(set(lefts) | set(rights)):
            if kind != 'cells' and not selected(identifier, args.fields):
                continue
            if identifier not in lefts or identifier not in rights:
                result[kind][identifier] = {'status': 'missing', 'side': 'left' if identifier not in lefts else 'right'}
            elif kind in ('fields', 'exports', 'cells'):
                result[kind][identifier] = compare_field(args.left, args.right, lefts[identifier], rights[identifier],
                                                         args.legacy_encoding and kind == 'fields')
            else:
                same = (args.left / lefts[identifier]['file']).read_bytes() == (args.right / rights[identifier]['file']).read_bytes()
                result[kind][identifier] = {'status': 'identical' if same else 'different'}

    problems = 0
    for kind in KINDS:
        prefix = {'exports': 'latlon:', 'cells': 'grid:'}.get(kind, '')
        for identifier, outcome in result.get(kind, {}).items():
            if outcome['status'] == 'identical':
                continue
            problems += 1
            details = ', '.join(f'{key}={value}' for key, value in outcome.items() if key != 'status')
            print(f'{outcome["status"]:>9}  {prefix}{identifier}  {details}')
    compared = sum(len(result.get(kind, {})) for kind in KINDS)
    extra = ', exports and cell arrays' if 'exports' in result or 'cells' in result else ''
    for key, outcome in result['grid'].items():
        problems += 1
        print(f'different  grid.{key}  left={outcome["left"]}, right={outcome["right"]}')
    print(f'{compared - problems + len(result["grid"])} of {compared} fields{extra} and tables identical')
    legacy = {identifier: outcome for identifier, outcome in result['fields'].items() if outcome.get('legacy')}
    if legacy:
        refined = sum(1 for identifier, outcome in legacy.items()
                      if outcome['lost_precision'] > 1e-3 / LEGACY_SCALES.get(identifier, 1))
        print(f'{len(legacy)} fields compared in their former integer encoding; {refined} now hold finer values')
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
