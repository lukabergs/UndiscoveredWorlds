"""Compare two regional bundles on their overlap.

    uv run --offline scripts/evaluation/regional_overlap.py A B [--json]

A regional bundle is a run bundle (uw-run-bundle/v1 or v2, grid kind latlon) written by `world_pipeline region`
whose manifest has a `window` block: the scale k, the world grid W x H and the absolute
fine indices (fine_column0, fine_row0) of the bundle's first cell on the planet-wide
kW x kH grid (docs/architecture/erosion-hydrology-plan.md, Phase E). The two bundles are
aligned through those origins (columns wrap at kW, so a window across the date line
aligns as any other), and every field the two hold is compared exactly on the cells they
share (scalars as one value). Two cube-window bundles (grid kind cube-window,
docs/architecture/sphere-mesh.md section 7) of one level are aligned by their cell_ids, the
level keys of their cells, whatever faces the windows cover. Printed per field: the cells
compared, the differing cells and the largest absolute difference; --json prints the result
as JSON instead.

Exit status: 0 when every common field is identical on the overlap, 1 on any difference,
2 when the bundles cannot be compared (not regional bundles, another scale, world grid,
halo, source manifest hash or manifest parameters, or no overlap).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import SCHEMAS, grid_kind, read_manifest


class Refused(Exception):
    """The two bundles cannot be compared."""


def window_of(path):
    manifest = read_manifest(path)
    if manifest.get('schema') not in SCHEMAS:
        raise Refused(f'{path}: unsupported schema {manifest.get("schema")!r}')
    if grid_kind(manifest) not in ('latlon', 'cube-window'):
        raise Refused(f'{path}: grid kind {grid_kind(manifest)!r}; only lat-lon and cube-window bundles are windows')
    window = manifest.get('window')
    if not isinstance(window, dict):
        raise Refused(f'{path}: no window block; not a regional bundle')
    return manifest, window


def axis_overlap(origin_a, length_a, origin_b, length_b, period=None):
    """Indices into A and into B of the fine cells both cover along one axis; with a period
    (the columns) absolute indices wrap."""
    absolute = origin_a + np.arange(length_a)
    if period is not None:
        absolute %= period
        position = (absolute - origin_b) % period
    else:
        position = absolute - origin_b
    inside = (position >= 0) & (position < length_b)
    return np.arange(length_a)[inside], position[inside]


def overlap(path_a, path_b):
    """{'columns', 'rows', 'fields': {id: {...}}} of two regional bundles; raises Refused."""
    manifest_a, window_a = window_of(path_a)
    manifest_b, window_b = window_of(path_b)
    for path, window in ((path_a, window_a), (path_b, window_b)):
        if 'pixel_columns' in window:
            raise Refused(f'{path} is resampled to {window.get("metres_per_pixel")} m pixels (--metres-per-pixel): its pixels '
                          'are anchored at its own window, so two windows are not comparable cell by cell; compare the '
                          'level bundles they were resampled from')
    for key in ('scale', 'world_columns', 'world_rows', 'halo'):
        if window_a.get(key) != window_b.get(key):
            raise Refused(f'the bundles differ in {key}: {window_a.get(key)} and {window_b.get(key)}')
    parameters_a, parameters_b = manifest_a.get('parameters', {}), manifest_b.get('parameters', {})
    if parameters_a != parameters_b:
        names = sorted(name for name in set(parameters_a) | set(parameters_b) if parameters_a.get(name) != parameters_b.get(name))
        raise Refused(f"the bundles were made with other parameters: {', '.join(names[:5])}"
                      + (f' and {len(names) - 5} more' if len(names) > 5 else ''))
    hash_a = manifest_a.get('source', {}).get('manifest_sha256')
    hash_b = manifest_b.get('source', {}).get('manifest_sha256')
    if not hash_a or hash_a != hash_b:
        raise Refused(f'the bundles come from different world bundles (source manifest_sha256 {hash_a} and {hash_b})')
    scale = window_a['scale']
    grid_a, grid_b = manifest_a['grid'], manifest_b['grid']
    kind = grid_kind(manifest_a)
    if grid_kind(manifest_b) != kind:
        raise Refused(f'the bundles are of grid kinds {kind!r} and {grid_kind(manifest_b)!r}')
    if kind == 'cube-window':
        return cube_overlap(path_a, path_b, manifest_a, manifest_b, scale)
    fine_columns = scale * window_a['world_columns']
    columns_a, columns_b = axis_overlap(window_a['fine_column0'], grid_a['columns'], window_b['fine_column0'],
                                        grid_b['columns'], fine_columns)
    rows_a, rows_b = axis_overlap(window_a['fine_row0'], grid_a['rows'], window_b['fine_row0'], grid_b['rows'])
    if not columns_a.size or not rows_a.size:
        raise Refused('the windows do not overlap')
    entries_a = {entry['id']: entry for entry in manifest_a['fields']}
    entries_b = {entry['id']: entry for entry in manifest_b['fields']}
    fields = {}
    for identifier in sorted(set(entries_a) & set(entries_b)):
        a = np.load(Path(path_a) / entries_a[identifier]['file'])
        b = np.load(Path(path_b) / entries_b[identifier]['file'])
        if a.dtype != b.dtype:
            fields[identifier] = {'status': 'dtype', 'a': str(a.dtype), 'b': str(b.dtype)}
            continue
        if a.ndim == 0 and b.ndim == 0:
            left, right = a.reshape(1), b.reshape(1)
        elif a.shape == (grid_a['rows'], grid_a['columns']) and b.shape == (grid_b['rows'], grid_b['columns']):
            left, right = a[np.ix_(rows_a, columns_a)], b[np.ix_(rows_b, columns_b)]
        else:
            fields[identifier] = {'status': 'shape', 'a': list(a.shape), 'b': list(b.shape)}
            continue
        fields[identifier] = compared(left, right)
    return {'a': str(path_a), 'b': str(path_b), 'scale': scale,
            'columns': int(columns_a.size), 'rows': int(rows_a.size),
            'only_a': sorted(set(entries_a) - set(entries_b)), 'only_b': sorted(set(entries_b) - set(entries_a)),
            'fields': fields}


def compared(left, right):
    """The outcome of two aligned arrays: cells, differing cells, largest absolute difference."""
    same = left == right
    if left.dtype.kind == 'f':
        same |= np.isnan(left) & np.isnan(right)
    differing = int(same.size - np.count_nonzero(same))
    difference = np.where(same, 0.0, np.abs(left.astype(np.float64) - right.astype(np.float64)))
    return {'status': 'identical' if differing == 0 else 'different', 'cells': int(same.size),
            'differing': differing, 'max_abs_difference': float(np.nanmax(difference)) if differing else 0.0}


def cube_overlap(path_a, path_b, manifest_a, manifest_b, scale):
    """The overlap of two cube-window bundles: the cells whose level keys (grid/cell_ids.npy) both list."""
    grid_a, grid_b = manifest_a['grid'], manifest_b['grid']
    for key in ('face_cells', 'level', 'base_face_cells', 'orientation'):
        if grid_a.get(key) != grid_b.get(key):
            raise Refused(f'the cube windows differ in {key}: {grid_a.get(key)} and {grid_b.get(key)}')

    def cell_ids(path, manifest):
        entries = {entry['id']: entry for entry in manifest.get('cells', [])}
        if 'cell_ids' not in entries:
            raise Refused(f'{path}: a cube-window bundle without grid/cell_ids.npy')
        return np.load(Path(path) / entries['cell_ids']['file'])

    ids_a, ids_b = cell_ids(path_a, manifest_a), cell_ids(path_b, manifest_b)
    shared, index_a, index_b = np.intersect1d(ids_a, ids_b, assume_unique=True, return_indices=True)
    if not shared.size:
        raise Refused('the windows do not overlap')
    count_a, count_b = grid_a['cell_count'], grid_b['cell_count']
    entries_a = {entry['id']: entry for entry in manifest_a['fields']}
    entries_b = {entry['id']: entry for entry in manifest_b['fields']}
    fields = {}
    for identifier in sorted(set(entries_a) & set(entries_b)):
        a = np.load(Path(path_a) / entries_a[identifier]['file'])
        b = np.load(Path(path_b) / entries_b[identifier]['file'])
        if a.dtype != b.dtype:
            fields[identifier] = {'status': 'dtype', 'a': str(a.dtype), 'b': str(b.dtype)}
            continue
        if a.ndim == 0 and b.ndim == 0:
            left, right = a.reshape(1), b.reshape(1)
        elif a.ndim >= 1 and a.shape[-1] == count_a and b.shape[-1] == count_b and a.shape[:-1] == b.shape[:-1]:
            left, right = a[..., index_a], b[..., index_b]  # [layer,] n: the cells on the last axis
        else:
            fields[identifier] = {'status': 'shape', 'a': list(a.shape), 'b': list(b.shape)}
            continue
        fields[identifier] = compared(left, right)
    return {'a': str(path_a), 'b': str(path_b), 'scale': scale, 'kind': 'cube-window', 'cells': int(shared.size),
            'only_a': sorted(set(entries_a) - set(entries_b)), 'only_b': sorted(set(entries_b) - set(entries_a)),
            'fields': fields}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('a', type=Path)
    parser.add_argument('b', type=Path)
    parser.add_argument('--json', action='store_true', help='Print the comparison as JSON')
    args = parser.parse_args(argv)
    try:
        result = overlap(args.a, args.b)
    except Refused as error:
        print(f'Refused: {error}', file=sys.stderr)
        return 2
    problems = [identifier for identifier, outcome in result['fields'].items() if outcome['status'] != 'identical']
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        if result.get('kind') == 'cube-window':
            print(f"Overlap: {result['cells']} level cells at scale {result['scale']}")
        else:
            print(f"Overlap: {result['columns']} x {result['rows']} fine cells at scale {result['scale']}")
        for identifier, outcome in result['fields'].items():
            if outcome['status'] in ('identical', 'different'):
                print(f"{outcome['status']:>9}  {identifier}  cells={outcome['cells']} differing={outcome['differing']} "
                      f"max_abs_difference={outcome['max_abs_difference']:.6g}")
            else:
                print(f"{outcome['status']:>9}  {identifier}  a={outcome['a']} b={outcome['b']}")
        for side in ('a', 'b'):
            if result[f'only_{side}']:
                print(f"only in {side.upper()}: {', '.join(result[f'only_{side}'])}")
        print(f'{len(result["fields"]) - len(problems)} of {len(result["fields"])} common fields identical on the overlap')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
