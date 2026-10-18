"""Face-edge seam statistics of a cube bundle's native fields (terrain roadmap T8 "Acceptance", cube imprint;
docs/architecture/sphere-mesh.md sections 2, 3 and 7).

For each field, the gradient of the native (6, N, N) field in face coordinates through the bundle's cell table
(uwtools.cell_grid.face_gradients: central differences over the edge-slot neighbours, great-circle distances
between centres, per km), split into the cells within `--band` cells of a face edge (default 2: the two outer
rows and columns of every face) and the face interiors. Per part: cells, median and mean |grad| (field units
per km), and the aspect harmonics in face coordinates on cells steeper than `--min-slope`: A8 = |mean e^(8 i
theta)| and A4 = mean cos(4 theta). ratio_median and ratio_mean are edge over interior (the acceptance: within
10 % of 1), a4_difference edge minus interior. SPEC is a field id or a sum of ids joined by '+', as stored;
--above V keeps the cells whose summed value is above V (e.g. sea level for land only).

Noise floor of the 10 % gate (measured 2026-10-19, out/scratch/sphere/W2-EV/seams-spread.json): on a seamless
isotropic random field (tests/scripts/test_cube_evaluation.py, wavenumber 20) ratio_median ranged 0.94 to 1.10
over seeds 1-5 at N 64 and 0.94 to 1.11 at N 128 (ratio_mean 0.98 to 1.07): the band samples about an eighth
of the sphere, so a field's own variability moves the ratio by up to 11 %, independent of N. A field with the
same |grad| everywhere gives 1.000. A ratio outside 0.9-1.1 on one realisation is therefore not yet a seam;
compare seeds, or the same field's ratio against a seamless reference of matching spectrum.

  uv run python scripts/evaluation/cube_seams.py BUNDLE --field SPEC[,SPEC] [--band 2] [--min-slope 0]
      [--above V] [--json OUT]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import cell_grid  # noqa: E402
from uwtools.bundle import Bundle  # noqa: E402
from uwtools.grid import EARTH_RADIUS_KM  # noqa: E402


def native_field(bundle, spec):
    total = None
    for identifier in spec.split('+'):
        if not bundle.has(identifier):
            raise SystemExit(f'{identifier}: not a field of {bundle.path}')
        info = bundle.fields[identifier]
        if info.layout != 'grid':
            raise SystemExit(f'{identifier}: layout {info.layout}, not a grid field')
        values = bundle.load_native(identifier).astype(np.float64)
        total = values if total is None else total + values
    return total


def part(magnitude, theta, mask, min_slope):
    steep = mask & (magnitude > min_slope)
    a8, a4 = cell_grid.aspect_harmonics(theta[steep])
    values = magnitude[mask]
    return {'cells': int(mask.sum()), 'aspect_cells': int(steep.sum()),
            'median_grad': float(np.median(values)) if values.size else None,
            'mean_grad': float(np.mean(values)) if values.size else None,
            'a8': a8, 'a4': a4}


def ratio(edge, interior):
    return edge / interior if edge is not None and interior else None


def seams(bundle, spec, band=2, min_slope=0.0, above=None):
    """The seam record of one field spec of a cube bundle."""
    if bundle.grid_kind != 'cube':
        raise SystemExit(f'{bundle.path}: grid kind {bundle.grid_kind}; seam statistics need a cube bundle')
    n = bundle.face_cells
    z = native_field(bundle, spec)
    if z.shape != (6, n, n):
        raise SystemExit(f'{spec}: native shape {list(z.shape)} is not (6, {n}, {n})')
    g_u, g_v, valid = cell_grid.face_gradients(bundle.cells(), z, bundle.radius_km or EARTH_RADIUS_KM)
    magnitude = np.where(valid, np.hypot(g_u, g_v), 0.0)
    theta = np.arctan2(np.where(valid, g_v, 0.0), np.where(valid, g_u, 0.0))
    keep = valid & (z.reshape(-1) > above if above is not None else True)
    near = cell_grid.face_edge_distance(n).reshape(-1) < band
    edge, interior = part(magnitude, theta, keep & near, min_slope), part(magnitude, theta, keep & ~near, min_slope)
    return {'edge': edge, 'interior': interior,
            'ratio_median': ratio(edge['median_grad'], interior['median_grad']),
            'ratio_mean': ratio(edge['mean_grad'], interior['mean_grad']),
            'a4_difference': (edge['a4'] - interior['a4']) if edge['a4'] is not None and interior['a4'] is not None
            else None}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('bundle', type=Path)
    parser.add_argument('--field', required=True, help='field ids (or sums joined by +), comma-separated')
    parser.add_argument('--band', type=int, default=2, help='cells from a face edge counted as the edge band')
    parser.add_argument('--min-slope', type=float, default=0.0, help='aspect only above this |grad| (units per km)')
    parser.add_argument('--above', type=float, help='keep cells whose value is above this')
    parser.add_argument('--json', type=Path)
    args = parser.parse_args(argv)
    bundle = Bundle.open(args.bundle)
    result = {'bundle': args.bundle.as_posix(),
              'grid': {'kind': bundle.grid_kind, 'face_cells': bundle.face_cells, 'radius_km': bundle.radius_km},
              'band_cells': args.band, 'min_slope': args.min_slope, 'above': args.above, 'fields': {}}
    for spec in args.field.split(','):
        record = seams(bundle, spec, args.band, args.min_slope, args.above)
        result['fields'][spec] = record
        print(f"{spec}: median |grad| edge {record['edge']['median_grad']:.4g}, interior "
              f"{record['interior']['median_grad']:.4g}, ratio {record['ratio_median']:.3f}; A4 edge "
              f"{record['edge']['a4']}, interior {record['interior']['a4']}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1), encoding='utf-8')
    return result


if __name__ == '__main__':
    main()
