"""Face-edge seam statistics of a cube bundle's native fields (terrain roadmap T8 "Acceptance", cube imprint;
docs/architecture/sphere-mesh.md sections 2, 3 and 7).

For each field, the gradient of the native (6, N, N) field in face coordinates through the bundle's cell table
(uwtools.cell_grid.face_gradients: central differences over the edge-slot neighbours, great-circle distances
between centres, per km), split into the cells within `--band` cells of a face edge (default 2: the two outer
rows and columns of every face) and the face interiors. Per part: cells, median and mean |grad| (field units
per km), and the aspect harmonics in face coordinates on cells steeper than `--min-slope`: A8 = |mean e^(8 i
theta)| and A4 = mean cos(4 theta). ratio_median and ratio_mean are edge over interior, a4_difference edge
minus interior. SPEC is a field id or a sum of ids joined by '+', as stored; --above V keeps the cells whose
summed value is above V (e.g. sea level for land only).

Seam gate by rank (decision of T8 W3-seam, replacing the earlier gate of a ratio within 10 % of 1): the band
samples about an eighth of the sphere, so a field's own variability moves ratio_median by up to 11 % on a
seamless field (measured 2026-10-19, out/scratch/sphere/W2-EV/seams-spread.json), independent of N. The same
band statistics are therefore computed for `--rotations` randomly rotated virtual cubes of the same N (default
200, seed `--seed`, deterministic; a cell is in the virtual band when its centre lies within `--band` cells of
an edge of the rotated cube; the gradient, its aspect angle in the true face coordinates and the kept cells are
unchanged). The rotated bands sample the same geography with no relation to the grid's faces, so they give the
spread of a seamless field of this geography. Each statistic (ratio_median, a4_difference) records under
`null` the 5/50/95 % of the rotations and the true value's percentile rank (the fraction of rotations below
it); it is flagged as a seam only when the true value lies above the null's 95th percentile. `seam` is true
when either statistic is flagged.

  uv run python scripts/evaluation/cube_seams.py BUNDLE --field SPEC[,SPEC] [--band 2] [--min-slope 0]
      [--above V] [--rotations 200] [--seed 1] [--json OUT]
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


def random_rotation(rng):
    """A uniformly random rotation matrix (from a normalised Gaussian quaternion)."""
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def virtual_edge_distance(points, n):
    """Distance in whole cells of each unit vector to the nearest edge of its face of an axis-aligned
    equiangular cube of N n."""
    a = np.abs(points)
    rows = np.arange(len(points))
    k = np.argmax(a, axis=1)
    major = a[rows, k]
    others = np.stack([a[rows, (k + 1) % 3], a[rows, (k + 2) % 3]], axis=1)
    angle = np.arctan(others / major[:, None]).max(axis=1)
    return np.floor((np.pi / 4 - angle) / (np.pi / 2 / n))


def band_statistics(magnitude, cos4, keep, steep, near):
    """(ratio_median, a4_difference) of the band `near` against the rest, over the kept cells."""
    edge, interior = magnitude[keep & near], magnitude[keep & ~near]
    middle = np.median(interior) if interior.size else 0.0
    ratio_median = float(np.median(edge) / middle) if edge.size and middle else None
    edge_a4, interior_a4 = cos4[steep & near], cos4[steep & ~near]
    a4 = float(edge_a4.mean() - interior_a4.mean()) if edge_a4.size and interior_a4.size else None
    return ratio_median, a4


def rank_record(true_value, null):
    """The true value's percentile rank in the null, the null's 5/50/95 % and the seam flag (above the 95th)."""
    null = np.asarray([v for v in null if v is not None], dtype=np.float64)
    if true_value is None or not null.size:
        return {'rank': None, 'null_p5_p50_p95': None, 'seam': False}
    p5, p50, p95 = (float(v) for v in np.percentile(null, [5, 50, 95]))
    return {'rank': float((null < true_value).mean()), 'null_p5_p50_p95': [p5, p50, p95],
            'seam': bool(true_value > p95)}


def seams(bundle, spec, band=2, min_slope=0.0, above=None, rotations=200, seed=1):
    """The seam record of one field spec of a cube bundle, with its rank against rotated virtual cubes."""
    if bundle.grid_kind != 'cube':
        raise SystemExit(f'{bundle.path}: grid kind {bundle.grid_kind}; seam statistics need a cube bundle')
    n = bundle.face_cells
    z = native_field(bundle, spec)
    if z.shape != (6, n, n):
        raise SystemExit(f'{spec}: native shape {list(z.shape)} is not (6, {n}, {n})')
    cells = bundle.cells()
    g_u, g_v, valid = cell_grid.face_gradients(cells, z, bundle.radius_km or EARTH_RADIUS_KM)
    magnitude = np.where(valid, np.hypot(g_u, g_v), 0.0)
    theta = np.arctan2(np.where(valid, g_v, 0.0), np.where(valid, g_u, 0.0))
    keep = valid & (z.reshape(-1) > above if above is not None else True)
    near = cell_grid.face_edge_distance(n).reshape(-1) < band
    edge, interior = part(magnitude, theta, keep & near, min_slope), part(magnitude, theta, keep & ~near, min_slope)
    record = {'edge': edge, 'interior': interior,
              'ratio_median': ratio(edge['median_grad'], interior['median_grad']),
              'ratio_mean': ratio(edge['mean_grad'], interior['mean_grad']),
              'a4_difference': (edge['a4'] - interior['a4']) if edge['a4'] is not None and interior['a4'] is not None
              else None}
    steep = keep & (magnitude > min_slope)
    cos4 = np.cos(4.0 * theta)
    centre = np.asarray(cells['centre'], dtype=np.float64).reshape(-1, 3)
    centre = centre / np.linalg.norm(centre, axis=1, keepdims=True)
    rng = np.random.default_rng(seed)
    null = [band_statistics(magnitude, cos4, keep, steep,
                            virtual_edge_distance(centre @ random_rotation(rng).T, n) < band)
            for _ in range(rotations)]
    record['null'] = {'rotations': rotations, 'seed': seed,
                      'ratio_median': rank_record(record['ratio_median'], [v[0] for v in null]),
                      'a4_difference': rank_record(record['a4_difference'], [v[1] for v in null])}
    record['seam'] = record['null']['ratio_median']['seam'] or record['null']['a4_difference']['seam']
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('bundle', type=Path)
    parser.add_argument('--field', required=True, help='field ids (or sums joined by +), comma-separated')
    parser.add_argument('--band', type=int, default=2, help='cells from a face edge counted as the edge band')
    parser.add_argument('--min-slope', type=float, default=0.0, help='aspect only above this |grad| (units per km)')
    parser.add_argument('--above', type=float, help='keep cells whose value is above this')
    parser.add_argument('--rotations', type=int, default=200, help='rotated virtual cubes of the null')
    parser.add_argument('--seed', type=int, default=1, help='seed of the rotations')
    parser.add_argument('--json', type=Path)
    args = parser.parse_args(argv)
    bundle = Bundle.open(args.bundle)
    result = {'bundle': args.bundle.as_posix(),
              'grid': {'kind': bundle.grid_kind, 'face_cells': bundle.face_cells, 'radius_km': bundle.radius_km},
              'band_cells': args.band, 'min_slope': args.min_slope, 'above': args.above,
              'rotations': args.rotations, 'seed': args.seed, 'fields': {}}
    for spec in args.field.split(','):
        record = seams(bundle, spec, args.band, args.min_slope, args.above, args.rotations, args.seed)
        result['fields'][spec] = record
        print(f"{spec}: median |grad| edge {record['edge']['median_grad']:.4g}, interior "
              f"{record['interior']['median_grad']:.4g}, ratio {record['ratio_median']:.3f}; A4 edge "
              f"{record['edge']['a4']}, interior {record['interior']['a4']}")
        for key in ('ratio_median', 'a4_difference'):
            gate = record['null'][key]
            null = gate['null_p5_p50_p95']
            print(f"  {key} {record[key]}: rank {gate['rank']}, null 5/50/95 % "
                  f"{[round(v, 4) for v in null] if null else None}, seam {gate['seam']}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1), encoding='utf-8')
    return result


if __name__ == '__main__':
    main()
