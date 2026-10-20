"""Item 3 (T1.6): straight river segments. Rivers as terrain_sandbox.river_coverage draws them: land cells (not sea,
not lake) with discharge >= rivers.min_discharge_m3_s (global / window), joined to the neighbour with the most
discharge when it has more. A straight run: consecutive river cells with the same step direction (one of 8); each
cell's run length is the run it lies in (cells). Reported: river cells, the share in runs of at least 5 and 8
cells, the mean run, on cells the procedural fill raised (fill_m > 1 m: the filled basins) and the rest.

Heading harmonics (item 3 step 0): on river cells with 16 river steps downstream, the chord to the cell 16 steps
on, its heading theta in true km (east-west scaled by the cosine of the chord's mean latitude; in grid cells
for '_grid'); A8 = |mean e^(8 i theta)| (0 isotropic, 1 every chord on one of the 8 grid directions), A4 = mean
cos(4 theta) (+ axes over diagonals, - diagonals over axes), a 16-bin heading histogram (% of chords, bin k
centred on k x 22.5 degrees counter-clockwise from east).

Routing classes (item 3 step 0) of the physics surface (-elevation-physics.npy), as flow_graph.cpp routes it:
the surface Priority-Flood filled (terrain_sandbox.priority_flood, the sea the outlets, the global wrapping);
'depression' a cell the fill raised (> 0 m; 'depression_1m' > 1 m), 'flat' a cell not raised with neither a
lower neighbour on the filled surface nor a sea neighbour, the rest 'descent'. On descent cells the steepest
D-infinity facet and its fraction (the share of randomiseddirection, flow_graph.cpp, offsets in km);
'axis_locked' a fraction within 0.05 of 0 or 1 (the draw keeps one neighbour). Shares in % of the cells
(river cells: all, on fill, off fill; and land cells). A domain without -discharge.npy (render.save_extra off)
gets the land cells' classes only.

Calibration (--calibrate): synthetic straight and gently curved lines (heading change uniform within +- 60
degrees over the line) at uniformly random angles (seed 1), each sampled every 0.05 cells into the cells it
crosses, the discharge the order along it, drawn by the same rule (runs()), the river the chain drawn from the
line's first cell (the crossed cells it skips dropped); reported: A8 / A4, the histogram,
the straight-run shares and A8 / A4 of the ideal chords (the continuous lines' own headings: the sampling floor).

A cube domain (its <domain>-*-native.npy arrays; master plan step 5.1) is measured on its own cells through its cell
table (<label>/cells-<N> for the global, <label>/window-<domain> for a window, whose kept cells alone are counted):
the same river rule over the eight slots; a straight run continues where the next slot is the step's straight
continuation (the opposite of its reverse slot, as CubeGrid.walks); the chord heading in the start cell's grid frame,
'a8'/'a4' from the true chord (the great-circle bearing to the end cell minus the start cell's slot-1 bearing, as
d8_axis.py), '_grid' from the lattice steps carried into the start cell's frame across face edges, '_compass' the
true chord's compass heading; routing classes with Priority-Flood on the slots (a window's outer ring, the cells
missing a cardinal slot, drains as in the exe), facets between consecutive slots on the table's tangent offsets;
'd8_axis' the mean cos 4 phi of the river steps (phi from the cell's slot-1 bearing; d8_axis.py's measure) with the
uniform-slot expectation, by |latitude| band; 'land d8_axis' the same on every land cell with a higher-discharge
neighbour (the receiver inferred from the discharge: approximate off the rivers).

  uv run python scripts/sandbox/straight_rivers.py OUT.json LABEL1,LABEL2 DOMAIN1,DOMAIN2
  uv run python scripts/sandbox/straight_rivers.py --calibrate OUT.json [LINES] [LENGTH]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path('C:/dev/undiscovered-worlds')
sys.path.insert(0, str(ROOT / 'scripts' / 'sandbox'))
import terrain_metrics as tm  # noqa: E402

OFFSETS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
CHORD_STEPS = 16
# Circular order (N, NE, E, SE, S, SW, W, NW): consecutive pairs are the D-infinity facets.
RING = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
LOCK = 0.05


def runs(discharge, water, threshold, wrap):
    rows, columns = discharge.shape
    q = np.where(water, 0.0, discharge)
    river = q >= threshold
    pad = np.pad(discharge, 1, mode='constant', constant_values=-np.inf)
    neighbours = np.stack([pad[1 + dr:1 + dr + rows, 1 + dc:1 + dc + columns] for dr, dc in OFFSETS])
    best = np.argmax(neighbours, axis=0)
    has = np.take_along_axis(neighbours, best[None], 0)[0] > discharge
    direction = np.where(river & has, best, -1)
    runs.direction = direction
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


def chord_headings(direction, lat, dy_km, steps=CHORD_STEPS):
    """Headings (radians, counter-clockwise from east) of the chord to the cell `steps` river steps downstream:
    (true km, grid cells); NaN where the path ends sooner. No wrap (as the renderer)."""
    rows, columns = direction.shape
    r, c = np.indices((rows, columns))
    cr, cc = r.copy(), c.copy()
    full = direction >= 0
    sr_of = np.array([o[0] for o in OFFSETS])
    sc_of = np.array([o[1] for o in OFFSETS])
    for _ in range(steps):
        d = direction[cr, cc]
        ok = d >= 0
        full &= ok
        cr = np.clip(cr + np.where(ok, sr_of[np.maximum(d, 0)], 0), 0, rows - 1)
        cc = np.clip(cc + np.where(ok, sc_of[np.maximum(d, 0)], 0), 0, columns - 1)
    lat = np.asarray(lat, dtype=float)
    mean_lat = 0.5 * (lat[r] + lat[cr])
    north, east = -(cr - r).astype(float), (cc - c).astype(float)
    km = np.arctan2(north * dy_km, east * dy_km * np.cos(np.radians(mean_lat)))
    grid = np.arctan2(north, east)
    return np.where(full, km, np.nan), np.where(full, grid, np.nan)


def harmonics(theta):
    theta = theta[np.isfinite(theta)]
    if theta.size == 0:
        return {'chords': 0}
    a8 = float(np.abs(np.mean(np.exp(8j * theta))))
    a4 = float(np.mean(np.cos(4.0 * theta)))
    bins = np.mod(np.floor((np.degrees(theta) + 11.25) / 22.5).astype(int), 16)
    histogram = np.bincount(bins, minlength=16) / theta.size
    return {'chords': int(theta.size), 'a8': round(a8, 4), 'a4': round(a4, 4),
            'histogram_pct': [round(100 * float(v), 2) for v in histogram]}


def heading_summary(mask, theta_km, theta_grid):
    out = harmonics(theta_km[mask])
    grid = harmonics(theta_grid[mask])
    if grid.get('chords'):
        out['a8_grid'], out['a4_grid'] = grid['a8'], grid['a4']
    return out


def neighbour(values, dr, dc, wrap):
    """values shifted so that [r, c] holds values[r + dr, c + dc]; and whether that neighbour exists."""
    rows, columns = values.shape
    r, c = np.indices((rows, columns))
    nr, nc = r + dr, c + dc
    valid = (nr >= 0) & (nr < rows)
    if wrap:
        nc = np.mod(nc, columns)
    else:
        valid &= (nc >= 0) & (nc < columns)
    return values[np.clip(nr, 0, rows - 1), np.clip(nc, 0, columns - 1)], valid


def routing_classes(physics, sea, lat, dy_km, wrap):
    """Masks 'depression', 'depression_1m', 'flat', 'descent' and the facet fraction (NaN off descent) of the
    physics surface (module docstring)."""
    import terrain_sandbox as ts
    filled = np.asarray(ts.priority_flood(physics, sea, wrap), dtype=float)
    raise_m = np.where(sea, 0.0, filled - physics)
    height = np.where(sea, 0.0, filled)  # a sea neighbour at sea level, as randomiseddirection
    dx = (dy_km * np.cos(np.radians(np.asarray(lat, dtype=float))))[:, None] * np.ones((1, physics.shape[1]))
    lower = np.zeros(physics.shape, dtype=bool)
    for dr, dc in RING:
        h, valid = neighbour(height, dr, dc, wrap)
        s, _ = neighbour(sea, dr, dc, wrap)
        lower |= valid & (s | (h < filled))
    best = np.zeros(physics.shape)
    fraction = np.full(physics.shape, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        for k in range(8):
            (dr1, dc1), (dr2, dc2) = RING[k], RING[(k + 1) % 8]
            h1, v1 = neighbour(height, dr1, dc1, wrap)
            h2, v2 = neighbour(height, dr2, dc2, wrap)
            drop1, drop2 = filled - h1, filled - h2
            p1e, p1n, p2e, p2n = dc1 * dx, -dr1 * dy_km, dc2 * dx, -dr2 * dy_km
            det = p1e * p2n - p1n * p2e
            ve = (drop1 * p2n - drop2 * p1n) / det
            vn = (p1e * drop2 - p2e * drop1) / det
            cross1 = p1e * vn - p1n * ve
            cross2 = ve * p2n - vn * p2e
            inside = ((ve != 0) | (vn != 0)) & (cross1 * det >= 0) & (cross2 * det >= 0)
            slope1, slope2 = drop1 / np.hypot(p1e, p1n), drop2 / np.hypot(p2e, p2n)
            slope = np.where(inside, np.hypot(ve, vn), np.maximum(slope1, slope2))
            share = np.where(inside, np.clip(np.abs(cross1) / (np.abs(cross1) + np.abs(cross2)), 0.0, 1.0),
                             np.where(slope2 > slope1, 1.0, 0.0))
            better = v1 & v2 & (slope > best)
            best = np.where(better, slope, best)
            fraction = np.where(better, share, fraction)
    depression = ~sea & (raise_m > 0.0)
    flat = ~sea & ~depression & ~lower
    descent = ~sea & ~depression & ~flat
    fraction = np.where(descent, fraction, np.nan)
    return {'depression': depression, 'depression_1m': ~sea & (raise_m > 1.0), 'flat': flat, 'descent': descent,
            'fraction': fraction}


def routing_summary(mask, classes):
    n = int(mask.sum())
    if n == 0:
        return {'cells': 0}
    f = classes['fraction'][mask]
    locked = np.isfinite(f) & ((f <= LOCK) | (f >= 1.0 - LOCK))
    defined = int(np.isfinite(f).sum())

    def pct(m):
        return round(100 * float(m[mask].mean()), 1)
    return {'cells': n, 'depression_pct': pct(classes['depression']), 'depression_1m_pct': pct(classes['depression_1m']),
            'flat_pct': pct(classes['flat']), 'depression_or_flat_pct': pct(classes['depression'] | classes['flat']),
            'descent_pct': pct(classes['descent']), 'axis_locked_pct': round(100 * float(locked.sum()) / n, 1),
            'axis_locked_pct_of_descent': round(100 * float(locked.sum()) / defined, 1) if defined else None,
            'fraction_p50_distance_to_lock': round(float(np.median(np.minimum(f[np.isfinite(f)], 1 - f[np.isfinite(f)]))), 3)
            if defined else None}


def cube_grid(d, domain, config):
    """(CubeGrid, kept mask) of a cube domain's saved arrays in label directory d."""
    import cube_domain
    sys.path.insert(0, str(ROOT / 'scripts'))
    from uwtools.bundle import Bundle
    path = Path(config['base_bundle'])
    base = Bundle.open(path if path.is_absolute() else ROOT / path)
    radius = float(base.radius_km or 6371.0)
    if domain == 'global':
        n = np.load(d / 'global-elevation-physics-native.npy', mmap_mode='r').shape[-1]
        grid = cube_domain.CubeGrid.from_directory(d / f'cells-{n}', n, radius)
        return grid, np.ones(grid.n, dtype=bool)
    target = d / f'window-{domain}'
    keys = np.load(target / 'keys.npy')
    grid = cube_domain.CubeGrid.from_directory(target, base.face_cells * config['domains'][domain]['scale'], radius, keys)
    return grid, np.load(target / 'kept.npy').astype(bool)


def cube_sea(grid, physics, procedural, edge):
    """terrain_metrics.sea_from on the slots: cells <= 0 m connected to the procedural sea (its components of at
    least SEA_MIN_CELLS or touching the domain's edge)."""
    labels, count = grid.label(procedural <= 0.0)
    sizes = np.bincount(labels, minlength=count + 1)
    keep = (sizes >= tm.SEA_MIN_CELLS) | np.isin(np.arange(count + 1), labels[edge])
    keep[0] = False
    seed = keep[labels]
    labels, count = grid.label(physics <= 0.0)
    hit = np.zeros(count + 1, dtype=bool)
    hit[labels[seed]] = True
    hit[0] = False
    return hit[labels]


def cube_continuation(grid, direction):
    """The slot (0-based) at the receiver that continues each step straight (the opposite of its reverse slot)."""
    s = np.maximum(direction, 0)
    reverse = grid.reverse[np.arange(grid.n), s]
    return (reverse - 1 + 4) % 8


def cube_runs(grid, discharge, water, threshold):
    """runs() on the slots: (river mask, run length); sets cube_runs.direction (0-based slot, -1 none),
    .receiver and runs.straightness (chord / path over 8 steps, km)."""
    n = grid.n
    ids = np.arange(n)
    q = np.where(water, 0.0, discharge)
    river = q >= threshold
    around = np.where(grid.valid, discharge[grid.safe], -np.inf)
    best = np.argmax(around, axis=1)
    has = around[ids, best] > discharge
    direction = np.where(river & has, best, -1)
    receiver = np.where(direction >= 0, grid.safe[ids, np.maximum(direction, 0)], ids)
    cont = cube_continuation(grid, direction)
    same_next = (direction >= 0) & (direction[receiver] == cont)
    ahead = np.ones(n, dtype=np.int32)
    for _ in range(400):
        new = np.where(same_next, 1 + ahead[receiver], 1)
        if np.array_equal(new, ahead):
            break
        ahead = new
    predecessor = np.full(n, -1)
    predecessor[receiver[same_next]] = ids[same_next]
    behind = np.ones(n, dtype=np.int32)
    for _ in range(400):
        new = np.where(predecessor >= 0, 1 + behind[np.maximum(predecessor, 0)], 1)
        if np.array_equal(new, behind):
            break
        behind = new
    length = ahead + behind - 1
    cur = ids.copy()
    path = np.zeros(n)
    full = river & (direction >= 0)
    for _ in range(8):
        dd = direction[cur]
        ok = dd >= 0
        full &= ok
        path += np.where(ok, grid.slot_km[cur, np.maximum(dd, 0)], 0.0)
        cur = np.where(ok, receiver[cur], cur)
    chord = grid.radius_km * np.arccos(np.clip(np.einsum('nd,nd->n', grid.centre[cur], grid.centre), -1.0, 1.0))
    runs.straightness = np.where(full, chord / np.maximum(path, 1e-9), np.nan)
    cube_runs.direction, cube_runs.receiver = direction, receiver
    return river & (direction >= 0), length


def cube_bearings(grid):
    """(n, 8) compass bearings of the slots (d8_axis.cube_bearings on the domain's table)."""
    sys.path.insert(0, str(ROOT / 'scripts' / 'evaluation'))
    import d8_axis
    bearings, _ = d8_axis.cube_bearings({'centre': grid.centre, 'neighbours': grid.neighbours,
                                         'latitude_deg': np.degrees(grid.latitude)})
    return bearings


def cube_chord_headings(grid, direction, receiver, bearings, steps=CHORD_STEPS):
    """Headings (radians, counter-clockwise from the start cell's slot 3 axis, so slot 1 is pi / 2) of the chord to
    the cell `steps` river steps downstream: (true chord, lattice steps, compass: counter-clockwise from east)."""
    n = grid.n
    ids = np.arange(n)
    cur = ids.copy()
    rotation = np.zeros(n, dtype=np.int64)
    east, north = np.zeros(n), np.zeros(n)
    full = direction >= 0
    step_east = np.array([0, 1, 1, 1, 0, -1, -1, -1], dtype=float)
    step_north = np.array([1, 1, 0, -1, -1, -1, 0, 1], dtype=float)
    for _ in range(steps):
        dd = direction[cur]
        ok = dd >= 0
        full &= ok
        local = np.maximum(dd, 0)
        start_slot = (local - rotation) % 8
        east += np.where(ok, step_east[start_slot], 0.0)
        north += np.where(ok, step_north[start_slot], 0.0)
        reverse = grid.reverse[cur, local]
        cont = (reverse - 1 + 4) % 8
        rotation = np.where(ok, (rotation + cont - local) % 8, rotation)
        cur = np.where(ok, receiver[cur], cur)
    c, e = grid.centre, grid.centre[cur]
    t = e - np.einsum('nd,nd->n', e, c)[:, None] * c
    bearing = np.arctan2(np.einsum('nd,dn->n', t, grid.east), np.einsum('nd,dn->n', t, grid.north))
    true = 0.5 * np.pi - (bearing - bearings[:, 0])
    compass = 0.5 * np.pi - bearing
    lattice = np.arctan2(north, east)
    return (np.where(full, true, np.nan), np.where(full, lattice, np.nan), np.where(full, compass, np.nan))


def cube_routing_classes(grid, physics, sea, outlet):
    """routing_classes on the slots: the filled surface from the outlets (sea and a window's edge), facets between
    consecutive slots (slot order N, NE, E, ... is the ring) on the table's tangent offsets (km)."""
    filled = np.asarray(grid.priority_flood(physics, outlet), dtype=float)
    raise_m = np.where(outlet, 0.0, filled - physics)
    height = np.where(outlet, 0.0, filled)
    h = np.where(grid.valid, height[grid.safe], np.inf)
    o = grid.valid & outlet[grid.safe]
    lower = (o | (h < filled[:, None])).any(1)
    best = np.zeros(grid.n)
    fraction = np.full(grid.n, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        for k in range(8):
            k2 = (k + 1) % 8
            drop1, drop2 = filled - h[:, k], filled - h[:, k2]
            p1e, p1n = grid.offset_east[:, k], grid.offset_north[:, k]
            p2e, p2n = grid.offset_east[:, k2], grid.offset_north[:, k2]
            det = p1e * p2n - p1n * p2e
            ve = (drop1 * p2n - drop2 * p1n) / det
            vn = (p1e * drop2 - p2e * drop1) / det
            cross1 = p1e * vn - p1n * ve
            cross2 = ve * p2n - vn * p2e
            inside = ((ve != 0) | (vn != 0)) & (cross1 * det >= 0) & (cross2 * det >= 0)
            slope1, slope2 = drop1 / grid.slot_km[:, k], drop2 / grid.slot_km[:, k2]
            slope = np.where(inside, np.hypot(ve, vn), np.maximum(slope1, slope2))
            share = np.where(inside, np.clip(np.abs(cross1) / (np.abs(cross1) + np.abs(cross2)), 0.0, 1.0),
                             np.where(slope2 > slope1, 1.0, 0.0))
            better = grid.valid[:, k] & grid.valid[:, k2] & (slope > best)
            best = np.where(better, slope, best)
            fraction = np.where(better, share, fraction)
    depression = ~outlet & (raise_m > 0.0)
    flat = ~outlet & ~depression & ~lower
    descent = ~outlet & ~depression & ~flat
    return {'depression': depression, 'depression_1m': ~outlet & (raise_m > 1.0), 'flat': flat, 'descent': descent,
            'fraction': np.where(descent, fraction, np.nan)}


def cube_axis(mask, direction, bearings, latitude_deg):
    """d8_axis.py's mean cos 4 phi of the steps in `mask`, with the uniform-slot expectation, by |latitude| band."""
    ids = np.flatnonzero(mask & (direction >= 0))
    relative = bearings[ids] - bearings[ids, :1]
    cos4 = np.cos(4.0 * relative)
    value = cos4[np.arange(ids.size), direction[ids]]
    present = np.isfinite(cos4)
    exact = np.where(present, cos4, 0.0).sum(1) / present.sum(1)
    axis = (direction[ids] % 2) == 0
    absolute = np.abs(latitude_deg[ids])
    out = {}
    for name, sub in (('all', np.ones(ids.size, dtype=bool)), ('|lat| 0-30', absolute < 30),
                      ('|lat| 30-60', (absolute >= 30) & (absolute < 60)), ('|lat| 60-90', absolute >= 60)):
        if sub.any():
            out[name] = {'cells': int(sub.sum()), 'mean_cos4': round(float(value[sub].mean()), 4),
                         'axis_fraction': round(float(axis[sub].mean()), 4),
                         'random_mean': round(float(exact[sub].mean()), 4)}
        else:
            out[name] = {'cells': 0}
    return out


def measure_cube(d, label, domain, config):
    """The record of a cube domain (module docstring)."""
    grid, kept = cube_grid(d, domain, config)

    def load(name):
        return np.load(d / f'{domain}-{name}-native.npy').astype(float).reshape(-1)
    phys, pre = load('elevation-physics'), load('elevation-procedural')
    lake = np.isfinite(load('lake-surface'))
    edge = ~grid.valid[:, list((0, 2, 4, 6))].all(1)
    sea = cube_sea(grid, phys, pre, edge)
    classes = cube_routing_classes(grid, phys, sea, sea | edge)
    land = ~sea & ~lake & kept
    record = {'grid': {'kind': 'cube' if grid.whole else 'cube-window', 'face_cells': grid.N, 'cells': grid.n,
                       'kept_cells': int(kept.sum())}}
    if not (d / f'{domain}-discharge-native.npy').is_file():
        record.update({'note': 'no -discharge-native.npy (render.save_extra off): land cells only',
                       'land routing': routing_summary(land, classes)})
        return record
    q = load('discharge')
    fill = load('fill_m')
    spec = config['rivers']['min_discharge_m3_s']
    threshold = spec['global' if domain == 'global' else 'window'] if isinstance(spec, dict) else spec
    river, length = cube_runs(grid, q, sea | lake, threshold)
    river &= kept
    direction, receiver = cube_runs.direction, cube_runs.receiver
    bearings = cube_bearings(grid)
    theta_km, theta_grid, theta_compass = cube_chord_headings(grid, direction, receiver, bearings)
    fields = {'theta_km': theta_km, 'theta_grid': theta_grid, 'classes': classes}
    filled = fill > 1.0
    record.update({'all': summary(river, length, fields), 'on fill': summary(river & filled, length, fields),
                   'off fill': summary(river & ~filled, length, fields),
                   'fill_land_pct': round(100 * float(filled[~sea & kept].mean()), 1),
                   'land routing': routing_summary(land, classes),
                   'heading16_compass': harmonics(theta_compass[river]),
                   'd8_axis': cube_axis(river, direction, bearings, np.degrees(grid.latitude))})
    # Every land cell's step to its highest-discharge neighbour (when higher): the routing tree, inferred.
    around = np.where(grid.valid, q[grid.safe], -np.inf)
    best = np.argmax(around, axis=1)
    higher = around[np.arange(grid.n), best] > q
    record['land d8_axis'] = cube_axis(land & higher, np.where(higher, best, -1), bearings, np.degrees(grid.latitude))
    return record


def summary(mask, length, fields=None):
    n = int(mask.sum())
    if n == 0:
        return {'cells': 0}
    v = length[mask]
    s = runs.straightness[mask]
    s = s[np.isfinite(s)]
    out = {'cells': n, 'run5_pct': round(100 * float((v >= 5).mean()), 1), 'run8_pct': round(100 * float((v >= 8).mean()), 1),
           'mean_run': round(float(v.mean()), 2), 'chord8_median': round(float(np.median(s)), 3) if s.size else None,
           'chord8_ge_0.98_pct': round(100 * float((s >= 0.98).mean()), 1) if s.size else None}
    if fields is not None:
        out['heading16'] = heading_summary(mask, fields['theta_km'], fields['theta_grid'])
        out['routing'] = routing_summary(mask, fields['classes'])
    return out


def arc_cells(rng, length, curved):
    """The cells (row, column) an arc of `length` cells crosses, in order, sampled every 0.05 cells; and the
    ideal chord headings (radians, x east, y south as rows) of its cells' positions `CHORD_STEPS` steps on."""
    theta0 = rng.uniform(0.0, 2.0 * math.pi)
    turn = rng.uniform(-math.pi / 3.0, math.pi / 3.0) if curved else 0.0
    s = np.arange(0.0, length, 0.05)
    heading = theta0 + turn * s / length
    if curved and turn != 0.0:
        kappa = turn / length
        x = (np.sin(heading) - math.sin(theta0)) / kappa
        y = -(np.cos(heading) - math.cos(theta0)) / kappa
    else:
        x, y = s * math.cos(theta0), s * math.sin(theta0)
    x, y = x + rng.uniform(0, 1), y + rng.uniform(0, 1)
    cx, cy = np.floor(x).astype(int), np.floor(y).astype(int)
    keep = np.ones(cx.size, dtype=bool)
    keep[1:] = (cx[1:] != cx[:-1]) | (cy[1:] != cy[:-1])
    return cy[keep], cx[keep], x, y


def calibrate(out, lines=3000, length=128):
    rng = np.random.default_rng(1)
    result = {'seed': 1, 'lines_per_set': lines, 'line_length_cells': length, 'chord_steps': CHORD_STEPS,
              'sampling_cells': 0.05, 'curved_turn_deg': 'uniform -60..60 over the line'}
    for name, curved in (('straight', False), ('curved', True)):
        thetas, ideals, lengths = [], [], []
        for _ in range(lines):
            rr, cc, x, y = arc_cells(rng, length, curved)
            r0, c0 = rr.min() - 2, cc.min() - 2
            rows, columns = rr.max() - r0 + 3, cc.max() - c0 + 3
            q = np.zeros((rows, columns))
            order = np.arange(1, rr.size + 1, dtype=float)
            q[rr - r0, cc - c0] = order  # a cell crossed twice keeps its later order
            runs(q, np.zeros(q.shape, dtype=bool), 0.5, False)
            # The river: the chain drawn from the line's first cell (a D8 tree's one cell per step), not the
            # crossed cells the drawing skips; on it the directions are unchanged.
            chain = np.zeros(q.shape, dtype=bool)
            r, c = rr[0] - r0, cc[0] - c0
            while True:
                chain[r, c] = True
                d = runs.direction[r, c]
                if d < 0:
                    break
                r, c = r + OFFSETS[d][0], c + OFFSETS[d][1]
            river, run_length = runs(np.where(chain, q, 0.0), np.zeros(q.shape, dtype=bool), 0.5, False)
            theta_km, _ = chord_headings(runs.direction, np.zeros(rows), 1.0)
            # North up: rows run south, so the drawn heading is atan2(-drow, dcol); the ideal line likewise.
            thetas.append(theta_km[river])
            lengths.append(run_length[river])
            ideal = np.arctan2(-(y[400:] - y[:-400]), x[400:] - x[:-400])  # 20-cell chords of the continuous line
            ideals.append(ideal[::20])
        theta = np.concatenate(thetas)
        v = np.concatenate(lengths)
        entry = harmonics(theta)
        ideal = harmonics(np.concatenate(ideals))
        entry.update({'ideal_a8': ideal['a8'], 'ideal_a4': ideal['a4'], 'cells': int(v.size),
                      'run5_pct': round(100 * float((v >= 5).mean()), 1), 'run8_pct': round(100 * float((v >= 8).mean()), 1),
                      'mean_run': round(float(v.mean()), 2)})
        result[name] = entry
        print(name, {k: entry[k] for k in ('chords', 'a8', 'a4', 'ideal_a8', 'ideal_a4', 'run5_pct', 'run8_pct')},
              flush=True)
    both = result['straight'], result['curved']
    result['threshold_note'] = ('acceptance "within 2x isotropic": A8 <= 2 x max(a8 of the two sets); |A4| <= 2 x max(|a4|)'
                                ' of the two sets')
    result['isotropic_a8_max'] = max(b['a8'] for b in both)
    result['isotropic_abs_a4_max'] = max(abs(b['a4']) for b in both)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(result, indent=1))


def main():
    if sys.argv[1] == '--calibrate':
        calibrate(sys.argv[2], *[int(a) for a in sys.argv[3:5]])
        return
    out, labels, domains = sys.argv[1], sys.argv[2].split(','), sys.argv[3].split(',')
    result = {}
    for label in labels:
        d = tm.SANDBOX / label
        config = json.loads((d / 'config.json').read_text())
        for domain in domains:
            if (d / f'{domain}-elevation-physics-native.npy').is_file():  # a cube domain: on its own cells
                result[f'{label}/{domain}'] = measure_cube(d, label, domain, config)
                print(label, domain, result[f'{label}/{domain}'], flush=True)
                continue
            if not (d / f'{domain}-elevation-physics.npy').is_file():
                continue
            phys = np.load(d / f'{domain}-elevation-physics.npy').astype(float)
            pre = np.load(d / f'{domain}-elevation-procedural.npy').astype(float)
            lake = np.isfinite(np.load(d / f'{domain}-lake-surface.npy'))
            sea = tm.sea_from(phys, seed=tm.sea_from(pre))
            _, lat, dy = tm.geometry(config, domain, phys.shape)
            classes = routing_classes(phys, sea, lat, dy, domain == 'global')
            land = ~sea & ~lake
            if not (d / f'{domain}-discharge.npy').is_file():
                result[f'{label}/{domain}'] = {'note': 'no -discharge.npy (render.save_extra off): land cells only',
                                               'land routing': routing_summary(land, classes)}
                print(label, domain, result[f'{label}/{domain}'], flush=True)
                continue
            q = np.load(d / f'{domain}-discharge.npy').astype(float)
            fill = np.load(d / f'{domain}-fill_m.npy').astype(float)
            spec = config['rivers']['min_discharge_m3_s']
            threshold = spec['global' if domain == 'global' else 'window'] if isinstance(spec, dict) else spec
            river, length = runs(q, sea | lake, threshold, False)  # as the renderer: no wrap
            theta_km, theta_grid = chord_headings(runs.direction, lat, dy)
            fields = {'theta_km': theta_km, 'theta_grid': theta_grid, 'classes': classes}
            filled = fill > 1.0
            result[f'{label}/{domain}'] = {'all': summary(river, length, fields), 'on fill': summary(river & filled, length, fields),
                                           'off fill': summary(river & ~filled, length, fields),
                                           'fill_land_pct': round(100 * float(filled[~sea].mean()), 1),
                                           'land routing': routing_summary(land, classes)}
            print(label, domain, result[f'{label}/{domain}'], flush=True)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
