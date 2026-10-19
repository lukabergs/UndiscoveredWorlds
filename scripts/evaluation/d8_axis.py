"""Grid-axis preference of a bundle's D8 directions (terrain roadmap T8 W3.4, the face-axis measure of W2-H's
review; docs/architecture/sphere-mesh.md section 2, slots and bearings).

For every cell whose direction code k is 1..8 with a neighbour in slot k, phi = bearing(c, k) - bearing(c, 1):
the angle of the flow from the cell's own slot-1 axis (on the cube the face's v axis, also across face edges;
on lat-lon the row axis, north, so phi is the compass bearing). The measure is the mean of cos 4 phi: +1 when
every flow follows a grid axis (slots 1, 3, 5, 7), -1 when every flow is diagonal, 0 without a preference for
either. It is given per band of |latitude| 0-30, 30-60, 60-90, per cube face and over all counted cells, with
the fraction of flows along an axis, against random directions: on the same cells a code drawn uniformly from
the cell's existing slots, `--draws` times (seed `--seed`, deterministic), whose exact mean (`random_mean`, the
grid's own geometry) and 5/50/95 % over the draws are recorded with the measured value's rank among the draws.

Bearings as the C++ grid's `bearing_rad` (section 2): on the cube the initial great-circle bearing to the
neighbour's centre from the cell table; on lat-lon the planar atan2(dx, dy) of the slot, dx the cell width of
the row (the mean of the two rows for a diagonal) and dy the row height.

  uv run python scripts/evaluation/d8_axis.py BUNDLE [BUNDLE ...] [--field hydrosphere.river_direction]
      [--draws 200] [--seed 1] [--json OUT]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import Bundle  # noqa: E402

# Slots 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW, 7 W, 8 NW as (east, north) steps.
STEPS = ((0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1))
BANDS = ((0.0, 30.0), (30.0, 60.0), (60.0, 90.0))


def cube_bearings(cells):
    """(n, 8) compass bearings in radians of the eight slots at each cell (NaN where the slot has no
    neighbour), and the cells' latitudes."""
    centre = np.asarray(cells['centre'], dtype=np.float64)
    centre = centre / np.linalg.norm(centre, axis=1, keepdims=True)
    neighbours = np.asarray(cells['neighbours'])
    east = np.cross(np.array([0.0, 0.0, 1.0]), centre)
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    north = np.cross(centre, east)
    bearings = np.full(neighbours.shape, np.nan)
    for slot in range(8):
        present = neighbours[:, slot] >= 0
        target = centre[neighbours[present, slot]]
        c = centre[present]
        # The initial great-circle direction to the target is the tangent component of target - c's chord.
        t = target - np.einsum('ij,ij->i', target, c)[:, None] * c
        bearings[present, slot] = np.arctan2(np.einsum('ij,ij->i', t, east[present]),
                                             np.einsum('ij,ij->i', t, north[present]))
    return bearings, np.asarray(cells['latitude_deg'], dtype=np.float64)


def latlon_bearings(rows, columns, latitudes):
    """(rows * columns, 8) planar slot bearings of a lat-lon grid (NaN beyond a pole) and the latitudes."""
    width = np.cos(np.deg2rad(latitudes)) * 2.0 * np.pi / columns
    height = np.pi / rows
    row_bearings = np.full((rows, 8), np.nan)
    for slot, (de, dn) in enumerate(STEPS):
        for r in range(rows):
            other = r - dn
            if not 0 <= other < rows:
                continue
            dx = width[r] if dn == 0 else 0.5 * (width[r] + width[other])
            row_bearings[r, slot] = np.arctan2(de * dx, dn * height)
    bearings = np.repeat(row_bearings, columns, axis=0)
    return bearings, np.repeat(latitudes, columns)


def groups(latitude, face):
    """Named boolean masks over the counted cells: all, the latitude bands, the faces (cube)."""
    result = {'all': np.ones(latitude.size, dtype=bool)}
    absolute = np.abs(latitude)
    for low, high in BANDS:
        result[f'|lat| {low:g}-{high:g}'] = (absolute >= low) & ((absolute < high) if high < 90.0 else True)
    if face is not None:
        for f in range(6):
            result[f'face {f}'] = face == f
    return result


def measure(bundle, field='hydrosphere.river_direction', draws=200, seed=1):
    """The record of one bundle: per group the measured mean cos 4 phi, the axis fraction and the random null."""
    if not bundle.has(field):
        raise SystemExit(f'{field}: not a field of {bundle.path}')
    codes = np.asarray(bundle.load_native(field)).reshape(-1).astype(np.int64)
    kind = bundle.grid_kind
    if kind == 'cube':
        bearings, latitude = cube_bearings(bundle.cells())
        n = bundle.face_cells
        face = np.arange(codes.size) // (n * n)
    elif kind == 'latlon':
        bearings, latitude = latlon_bearings(bundle.rows, bundle.columns, np.asarray(bundle.latitudes()))
        face = None
    else:
        raise SystemExit(f'{bundle.path}: grid kind {kind}; the measure needs a cube or lat-lon world bundle')
    if bearings.shape[0] != codes.size:
        raise SystemExit(f'{field}: {codes.size} values for {bearings.shape[0]} cells')
    relative = bearings - bearings[:, :1]  # phi of every slot from the cell's slot 1
    cos4 = np.cos(4.0 * relative)  # NaN where the slot has no neighbour
    ids = np.arange(codes.size)
    coded = (codes >= 1) & (codes <= 8)
    value = np.full(codes.size, np.nan)
    value[coded] = cos4[ids[coded], codes[coded] - 1]
    counted = coded & np.isfinite(value)
    value, cos4, latitude = value[counted], cos4[counted], latitude[counted]
    axis = (codes[counted] % 2) == 1
    face = face[counted] if face is not None else None
    present = np.isfinite(cos4)
    exact = np.where(present, cos4, 0.0).sum(axis=1) / present.sum(axis=1)  # the expectation of a uniform code
    # Random codes: the k-th existing slot of each cell for k uniform in 0 .. (existing slots - 1).
    rng = np.random.default_rng(seed)
    order = np.argsort(~present, axis=1, kind='stable')  # existing slots first, in slot order
    available = present.sum(axis=1)
    rows = np.arange(value.size)
    random_values = np.empty((draws, value.size))
    for d in range(draws):
        pick = order[rows, (rng.random(value.size) * available).astype(np.int64)]
        random_values[d] = cos4[rows, pick]
    record = {'grid': {'kind': kind, 'face_cells': bundle.face_cells, 'columns': bundle.columns,
                       'rows': bundle.rows},
              'field': field, 'cells': int(counted.sum()), 'groups': {}}
    for name, mask in groups(latitude, face).items():
        if not mask.any():
            record['groups'][name] = {'cells': 0}
            continue
        mean = float(value[mask].mean())
        null = random_values[:, mask].mean(axis=1)
        p5, p50, p95 = (float(v) for v in np.percentile(null, [5, 50, 95]))
        record['groups'][name] = {'cells': int(mask.sum()), 'mean_cos4': mean,
                                  'axis_fraction': float(axis[mask].mean()),
                                  'random_mean': float(exact[mask].mean()),
                                  'random_p5_p50_p95': [p5, p50, p95],
                                  'rank': float((null < mean).mean())}
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('bundles', type=Path, nargs='+')
    parser.add_argument('--field', default='hydrosphere.river_direction', help='a D8 direction-code field')
    parser.add_argument('--draws', type=int, default=200, help='random direction fields of the null')
    parser.add_argument('--seed', type=int, default=1, help='seed of the random directions')
    parser.add_argument('--json', type=Path)
    args = parser.parse_args(argv)
    result = {'field': args.field, 'draws': args.draws, 'seed': args.seed, 'bundles': {}}
    for path in args.bundles:
        bundle = Bundle.open(path)
        record = measure(bundle, args.field, args.draws, args.seed)
        record['seed_of_world'] = bundle.seed
        result['bundles'][path.as_posix()] = record
        print(f"{path.as_posix()} ({record['grid']['kind']}, {record['cells']} cells)")
        for name, group in record['groups'].items():
            if not group['cells']:
                continue
            p5, _, p95 = group['random_p5_p50_p95']
            print(f"  {name}: {group['cells']} cells, mean cos 4 phi {group['mean_cos4']:.3f} (axis "
                  f"{group['axis_fraction']:.3f}); random {group['random_mean']:.3f} [{p5:.3f}, {p95:.3f}]")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1), encoding='utf-8')
    return result


if __name__ == '__main__':
    main()
