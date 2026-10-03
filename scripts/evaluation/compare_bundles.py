# /// script
# dependencies = ["numpy"]
# ///
"""Compare two run bundles field by field.

    uv run --offline scripts/evaluation/compare_bundles.py BUNDLE_A BUNDLE_B [--fields PREFIX ...] [--json OUT]

A bundle is a directory written by `world_pipeline --bundle` (schema uw-run-bundle/v1).
Fields are compared exactly; for differing fields the report gives the number of
differing values and the largest absolute difference. Tables are compared as text.
Exit status: 0 when every compared field and table is identical, 1 otherwise.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def load_manifest(bundle):
    manifest = json.loads((bundle / 'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema') != 'uw-run-bundle/v1':
        raise SystemExit(f'{bundle}: unsupported schema {manifest.get("schema")!r}')
    return manifest


def selected(identifier, prefixes):
    return not prefixes or any(identifier.startswith(prefix) for prefix in prefixes)


def compare_field(left, right, entry_left, entry_right):
    a = np.load(left / entry_left['file'])
    b = np.load(right / entry_right['file'])
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
    args = parser.parse_args()

    left, right = load_manifest(args.left), load_manifest(args.right)
    result = {'left': str(args.left), 'right': str(args.right), 'fields': {}, 'tables': {}}
    for kind in ('fields', 'tables'):
        lefts = {entry['id']: entry for entry in left[kind]}
        rights = {entry['id']: entry for entry in right[kind]}
        for identifier in sorted(set(lefts) | set(rights)):
            if not selected(identifier, args.fields):
                continue
            if identifier not in lefts or identifier not in rights:
                result[kind][identifier] = {'status': 'missing', 'side': 'left' if identifier not in lefts else 'right'}
            elif kind == 'fields':
                result[kind][identifier] = compare_field(args.left, args.right, lefts[identifier], rights[identifier])
            else:
                same = (args.left / lefts[identifier]['file']).read_bytes() == (args.right / rights[identifier]['file']).read_bytes()
                result[kind][identifier] = {'status': 'identical' if same else 'different'}

    problems = 0
    for kind in ('fields', 'tables'):
        for identifier, outcome in result[kind].items():
            if outcome['status'] == 'identical':
                continue
            problems += 1
            details = ', '.join(f'{key}={value}' for key, value in outcome.items() if key != 'status')
            print(f'{outcome["status"]:>9}  {identifier}  {details}')
    compared = len(result['fields']) + len(result['tables'])
    print(f'{compared - problems} of {compared} fields and tables identical')
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
