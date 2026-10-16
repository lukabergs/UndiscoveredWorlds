"""Item 3 (T1.6): straight river segments. Rivers as terrain_sandbox.river_coverage draws them: land cells (not sea,
not lake) with discharge >= rivers.min_discharge_m3_s (global / window), joined to the neighbour with the most
discharge when it has more. A straight run: consecutive river cells with the same step direction (one of 8); each
cell's run length is the run it lies in (cells). Reported: river cells, the share in runs of at least 5 and 8
cells, the mean run, on cells the procedural fill raised (fill_m > 1 m: the filled basins) and the rest.

  uv run python scripts/sandbox/straight_rivers.py OUT.json LABEL1,LABEL2 DOMAIN1,DOMAIN2
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path('C:/dev/undiscovered-worlds')
sys.path.insert(0, str(ROOT / 'scripts' / 'sandbox'))
import terrain_metrics as tm  # noqa: E402

OFFSETS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def runs(discharge, water, threshold, wrap):
    rows, columns = discharge.shape
    q = np.where(water, 0.0, discharge)
    river = q >= threshold
    pad = np.pad(discharge, 1, mode='constant', constant_values=-np.inf)
    neighbours = np.stack([pad[1 + dr:1 + dr + rows, 1 + dc:1 + dc + columns] for dr, dc in OFFSETS])
    best = np.argmax(neighbours, axis=0)
    has = np.take_along_axis(neighbours, best[None], 0)[0] > discharge
    direction = np.where(river & has, best, -1)
    r, c = np.indices((rows, columns))
    dr = np.array([o[0] for o in OFFSETS])[np.maximum(direction, 0)]
    dc = np.array([o[1] for o in OFFSETS])[np.maximum(direction, 0)]
    nr, nc = r + dr, c + dc
    nc = np.mod(nc, columns) if wrap else np.clip(nc, 0, columns - 1)
    nr = np.clip(nr, 0, rows - 1)
    same_next = (direction >= 0) & (direction[nr, nc] == direction)
    ahead = np.ones((rows, columns), dtype=np.int32)
    for _ in range(400):
        new = np.where(same_next, 1 + ahead[nr, nc], 1)
        if np.array_equal(new, ahead):
            break
        ahead = new
    # behind: the upstream cell with the same direction is unique (the cell one step back along it)
    pr, pc = r - dr, c - dc
    pc = np.mod(pc, columns) if wrap else np.clip(pc, 0, columns - 1)
    pr = np.clip(pr, 0, rows - 1)
    same_prev = (direction >= 0) & (direction[pr, pc] == direction) & ((pr != r) | (pc != c))
    behind = np.ones((rows, columns), dtype=np.int32)
    for _ in range(400):
        new = np.where(same_prev, 1 + behind[pr, pc], 1)
        if np.array_equal(new, behind):
            break
        behind = new
    length = ahead + behind - 1
    # Chord over path along 8 downstream steps (cell units, diagonals sqrt 2), on cells with 8 river steps ahead.
    steps = 8
    cr, cc = r.copy(), c.copy()
    path = np.zeros((rows, columns))
    full = river & (direction >= 0)
    for _ in range(steps):
        d = direction[cr, cc]
        ok = d >= 0
        full &= ok
        sr = np.array([o[0] for o in OFFSETS])[np.maximum(d, 0)]
        sc = np.array([o[1] for o in OFFSETS])[np.maximum(d, 0)]
        path += np.where(ok, np.hypot(sr, sc), 0.0)
        cr, cc = np.clip(cr + np.where(ok, sr, 0), 0, rows - 1), np.clip(cc + np.where(ok, sc, 0), 0, columns - 1)
    chord = np.hypot(cr - r, cc - c)
    runs.straightness = np.where(full, chord / np.maximum(path, 1e-9), np.nan)
    return river & (direction >= 0), length


def summary(mask, length):
    n = int(mask.sum())
    if n == 0:
        return {'cells': 0}
    v = length[mask]
    s = runs.straightness[mask]
    s = s[np.isfinite(s)]
    return {'cells': n, 'run5_pct': round(100 * float((v >= 5).mean()), 1), 'run8_pct': round(100 * float((v >= 8).mean()), 1),
            'mean_run': round(float(v.mean()), 2), 'chord8_median': round(float(np.median(s)), 3) if s.size else None,
            'chord8_ge_0.98_pct': round(100 * float((s >= 0.98).mean()), 1) if s.size else None}


def main():
    out, labels, domains = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
    result = {}
    for label in labels:
        d = tm.SANDBOX / label
        config = json.loads((d / 'config.json').read_text())
        for domain in domains:
            if not (d / f'{domain}-discharge.npy').is_file():
                continue
            q = np.load(d / f'{domain}-discharge.npy').astype(float)
            fill = np.load(d / f'{domain}-fill_m.npy').astype(float)
            phys = np.load(d / f'{domain}-elevation-physics.npy').astype(float)
            pre = np.load(d / f'{domain}-elevation-procedural.npy').astype(float)
            lake = np.isfinite(np.load(d / f'{domain}-lake-surface.npy'))
            sea = tm.sea_from(phys, seed=tm.sea_from(pre))
            spec = config['rivers']['min_discharge_m3_s']
            threshold = spec['global' if domain == 'global' else 'window'] if isinstance(spec, dict) else spec
            river, length = runs(q, sea | lake, threshold, False)  # as the renderer: no wrap
            filled = fill > 1.0
            result[f'{label}/{domain}'] = {'all': summary(river, length), 'on fill': summary(river & filled, length),
                                           'off fill': summary(river & ~filled, length),
                                           'fill_land_pct': round(100 * float(filled[~sea].mean()), 1)}
            print(label, domain, result[f'{label}/{domain}'], flush=True)
    Path(out).write_text(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
