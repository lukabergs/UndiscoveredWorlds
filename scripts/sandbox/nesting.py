"""The nesting mechanism's grid operations, shared by terrain_sandbox.py (rule B, every level below the global
domain) and nested_window.py (the 100 m test): Catmull-Rom prolongation of a parent array, block means, the sea
flood, the window rings, the projections onto the parent's cell means and the inflow of a parent's rivers at a
child domain's edge. Pure numpy / scipy, no sandbox state (multiscale plan section 2, "The exact projection").

Index convention: a child domain covers the parent cells [row0, row0 + rows) x [col0, col0 + columns) of a parent
array; child cell (i, j) lies in parent cell (row0 + i // factor, col0 + j // factor).
"""
import math

import numpy as np
from scipy import ndimage

EIGHT = np.ones((3, 3), bool)
SECONDS_PER_YEAR = 365.25 * 86400.0  # terrain_sandbox_main.cpp's secondsperyear


# ---------------------------------------------------------------- grids

def catmull_rom_matrix(count, factor, start, size):
    """Rows of Catmull-Rom weights from `size` parent samples to count x factor children covering the parents
    start .. start + count - 1 (indices clamped to the array), and the two bracketing parent indices per child."""
    children = np.arange(count * factor)
    p = start + (children + 0.5) / factor - 0.5
    i0 = np.floor(p).astype(np.int64)
    t = p - i0
    weights = [(-t ** 3 + 2 * t ** 2 - t) / 2, (3 * t ** 3 - 5 * t ** 2 + 2) / 2, (-3 * t ** 3 + 4 * t ** 2 + t) / 2,
               (t ** 3 - t ** 2) / 2]
    matrix = np.zeros((children.size, size))
    for offset, w in zip(range(-1, 3), weights):
        np.add.at(matrix, (children, np.clip(i0 + offset, 0, size - 1)), w)
    return matrix, np.clip(i0, 0, size - 1), np.clip(i0 + 1, 0, size - 1)


def prolong(parent, factor, row0, col0, rows, columns, limit=True):
    """The parent array's cells [row0, row0 + rows) x [col0, col0 + columns) prolonged by `factor` (Catmull-Rom,
    clamped to the 2 x 2 bracketing parents with `limit`)."""
    my, r0, r1 = catmull_rom_matrix(rows, factor, row0, parent.shape[0])
    mx, c0, c1 = catmull_rom_matrix(columns, factor, col0, parent.shape[1])
    values = my @ np.asarray(parent, dtype=np.float64) @ mx.T
    if limit:
        corners = [parent[np.ix_(r, c)] for r in (r0, r1) for c in (c0, c1)]
        values = np.clip(values, np.minimum.reduce(corners), np.maximum.reduce(corners))
    return values


def nearest(parent, factor, row0, col0, rows, columns):
    return np.repeat(np.repeat(parent[row0:row0 + rows, col0:col0 + columns], factor, axis=0), factor, axis=1)


def block_mean(values, factor):
    rows, columns = values.shape
    return values.reshape(rows // factor, factor, columns // factor, factor).mean(axis=(1, 3))


def flood(low, seeds):
    """The cells of `low` 8-connected to a seed cell of `low`."""
    labels, _ = ndimage.label(low, structure=EIGHT)
    touching = np.unique(labels[seeds & (labels > 0)])
    return np.isin(labels, touching[touching > 0])


def ring(shape, depth):
    mask = np.zeros(shape, bool)
    mask[:depth, :] = mask[-depth:, :] = mask[:, :depth] = mask[:, -depth:] = True
    return mask


# ---------------------------------------------------------------- projections onto the parent's cell means

def project_means(values, targets, movable, factor, passes, relaxation):
    """Values with their mean over every parent cell moved to `targets` exactly, moving the `movable` children only:
    `passes` smooth passes (relaxation x the Catmull-Rom prolongation of the block departures over the movable share),
    then the remainder spread evenly over each block's movable children. Blocks without one keep their departure.
    (nested_window.py's projection; unbounded: a few land children beside deep sea carry a block's whole departure.)"""
    share = block_mean(movable.astype(float), factor)
    rows, columns = share.shape
    with np.errstate(invalid='ignore', divide='ignore'):
        for _ in range(passes):
            per_cell = np.where(share > 0, (block_mean(values, factor) - targets) / share, 0.0)
            values = values - relaxation * movable * prolong(per_cell, factor, 0, 0, rows, columns, limit=False)
        per_cell = np.where(share > 0, (block_mean(values, factor) - targets) / share, 0.0)
    return values - movable * nearest(per_cell, factor, 0, 0, rows, columns)


def remove_parent_means(detail, land, factor, passes, relaxation):
    """Detail with its mean over every parent cell removed exactly, moving land children only (project_means to 0).
    Blocks without land carry no detail."""
    return project_means(np.where(land, detail, 0.0), 0.0, land, factor, passes, relaxation)


def to_blocks(values, factor):
    """(rows, columns) -> (blocks, factor^2), blocks row-major over the parent cells."""
    rows, columns = values.shape
    return values.reshape(rows // factor, factor, columns // factor, factor).transpose(0, 2, 1, 3).reshape(
        (rows // factor) * (columns // factor), factor * factor)


def from_blocks(blocks, factor, shape):
    rows, columns = shape
    return blocks.reshape(rows // factor, columns // factor, factor, factor).transpose(0, 2, 1, 3).reshape(shape)


def bracket_extremes(parent, factor, row0, col0, rows, columns):
    """The least and the greatest of the 2 x 2 bracketing parents of every child (the range prolong clamps to)."""
    _, r0, r1 = catmull_rom_matrix(rows, factor, row0, parent.shape[0])
    _, c0, c1 = catmull_rom_matrix(columns, factor, col0, parent.shape[1])
    corners = [parent[np.ix_(r, c)] for r in (r0, r1) for c in (c0, c1)]
    return np.minimum.reduce(corners), np.maximum.reduce(corners)


def project_bounded(values, targets, movable, land, factor, land_floor_m, sea_ceiling_m, cap=None, clip=None):
    """Values with their mean over every parent cell set to `targets` exactly (multiscale plan section 3, S2.3's
    "bounds, not weights"): the same shift on every movable child of a parent cell (land and sea alike), each child
    capped by its headroom in the shift's direction: a land child never goes below land_floor_m, a sea child never
    above sea_ceiling_m, and with `cap` (per child, e.g. the bracketing parents' maximum) no child rises above it
    (without: upward on land and downward at sea are unbounded). The shift is found by water-filling inside the
    block (children with less headroom stop at it, the rest share the remainder). When the caps exhaust, the
    remainder is spread evenly over the movable children ignoring the caps (counted as forced; it may flip a child's
    class); with `clip` (per child) an upward remainder stops where a child would rise above its clip, and the block
    keeps the rest of its departure (counted as clipped: those blocks are not exact). Blocks without a movable child
    keep their departure. Returns the values and a dict of counts."""
    shape = values.shape
    target = np.asarray(targets, dtype=np.float64).ravel()
    v = to_blocks(np.asarray(values, dtype=np.float64), factor)
    move = to_blocks(movable, factor)
    is_land = to_blocks(land, factor)
    n = factor * factor
    departure = target - v.mean(axis=1)  # the block mean must rise by this
    need = np.abs(departure) * n  # the sum of the shifts' magnitudes
    up = departure > 0
    big = 1e12
    ceiling = np.where(is_land, big, sea_ceiling_m)
    if cap is not None:
        ceiling = np.minimum(ceiling, to_blocks(np.asarray(cap, dtype=np.float64), factor))
    headroom_up = np.maximum(ceiling - v, 0.0)
    headroom_down = np.where(is_land, np.maximum(v - land_floor_m, 0.0), big)
    headroom = np.where(up[:, None], headroom_up, headroom_down)
    headroom = np.where(move, headroom, 0.0)
    order = np.sort(headroom, axis=1)
    below = np.zeros_like(order)  # the sum of the smaller headrooms (exclusive: no 1e12 cap rounds it)
    below[:, 1:] = np.cumsum(order[:, :-1], axis=1)
    reach = below + order * (n - np.arange(n))[None, :]  # the total moved if the shift were order[:, j]
    feasible = reach >= need[:, None]
    first = np.where(feasible.any(axis=1), feasible.argmax(axis=1), n)
    rows_index = np.arange(v.shape[0])
    safe = np.minimum(first, n - 1)
    shift = np.where(first < n, (need - below[rows_index, safe]) / np.maximum(n - safe, 1), big)
    moved = np.minimum(headroom, shift[:, None])
    count = move.sum(axis=1)
    forced = (first >= n) & (count > 0) & (need > 0)
    remainder = np.where(forced, need - headroom.sum(axis=1), 0.0) / np.maximum(count, 1)
    clipped = np.zeros_like(forced)
    if clip is not None:
        room = np.where(move, to_blocks(np.asarray(clip, dtype=np.float64), factor) - (v + moved), big)
        allowed = np.maximum(room.min(axis=1), 0.0)
        clipped = forced & up & (remainder > allowed)
        remainder = np.where(clipped, allowed, remainder)
    moved = moved + np.where(move, remainder[:, None], 0.0)
    signed = np.where(up[:, None], moved, -moved)
    signed = np.where(move & (count > 0)[:, None], signed, 0.0)
    out = from_blocks(v + signed, factor, shape)
    residual = np.abs(target - (v + signed).mean(axis=1))
    stats = {'blocks': int(v.shape[0]), 'blocks_without_movable_child': int(np.count_nonzero(count == 0)),
             'forced_blocks': int(np.count_nonzero(forced)), 'clipped_blocks': int(np.count_nonzero(clipped)),
             'clipped_residual_max_m': float(residual[clipped].max()) if clipped.any() else 0.0,
             'sea_children_moved': int(np.count_nonzero(move & ~is_land & (np.abs(signed) > 1e-9))),
             'largest_shift_m': float(np.abs(signed).max()) if signed.size else 0.0,
             'shift_p99_m': float(np.percentile(np.abs(signed[move]), 99)) if move.any() else 0.0}
    return out, stats, from_blocks(np.repeat(clipped[:, None], n, axis=1), factor, shape)


def class_bounds(values, land, land_floor_m, sea_ceiling_m, cap):
    """Per child, the lowest and highest value a projection may move it to: land never below land_floor_m, sea never
    above sea_ceiling_m, nothing above `cap` (None: no cap); a child already outside keeps its own value as the bound
    (a projection never pushes it back)."""
    big = 1e12
    lower = np.where(land, land_floor_m, -big)
    upper = np.where(land, big, sea_ceiling_m)
    if cap is not None:
        upper = np.minimum(upper, cap)
    return np.minimum(values, lower), np.maximum(values, upper)


def reference_weights(values, reference, movable, factor):
    """The exact remainder's weights: each child's height above its block's lowest reference cell (channel, sea or
    lake; never below 0), 1 in a block without a reference cell or whose movable children all sit at it."""
    v = to_blocks(values, factor)
    ref = to_blocks(reference, factor)
    move = to_blocks(movable, factor)
    lowest = np.where(ref, v, np.inf).min(axis=1)
    w = np.where(np.isfinite(lowest)[:, None], np.maximum(v - np.where(np.isfinite(lowest), lowest, 0.0)[:, None], 0.0),
                 1.0)
    w = np.where(move, w, 0.0)
    flat = w.sum(axis=1) <= 0.0
    w = np.where(flat[:, None] & move, 1.0, w)
    return from_blocks(w, factor, values.shape)


def project_weighted(values, targets, movable_up, movable_down, land, weights, factor, land_floor_m, sea_ceiling_m,
                     cap=None, clip=None, clip_low=None):
    """Values with their mean over every parent cell set to `targets` exactly, the correction proportional to
    `weights` (multiscale plan section 2: F = F* + lambda w / mean(w)): a block that must rise moves its movable_up
    children, one that must fall its movable_down children; each child stops at its bound (class_bounds: the land
    floor, the sea ceiling, `cap`), the others share the rest (weighted water-filling). When the bounds exhaust, the
    rest is spread over the moving children in proportion to their weights ignoring the bounds (forced; uniform
    when the weights are all 0), but an upward rest stops where a child would pass `clip` and a downward one where
    it would pass below `clip_low` (the block keeps that departure: clipped). Blocks without a moving child keep their departure. Returns the values, counts and the
    clipped blocks' mask."""
    shape = values.shape
    target = np.asarray(targets, dtype=np.float64).ravel()
    values = np.asarray(values, dtype=np.float64)
    lower, upper = class_bounds(values, land, land_floor_m, sea_ceiling_m, cap)
    v = to_blocks(values, factor)
    n = factor * factor
    departure = target - v.mean(axis=1)
    need = np.abs(departure) * n
    up = departure > 0
    move = np.where(up[:, None], to_blocks(movable_up, factor), to_blocks(movable_down, factor))
    headroom = np.where(up[:, None], to_blocks(upper, factor) - v, v - to_blocks(lower, factor))
    headroom = np.where(move, np.maximum(headroom, 0.0), 0.0)
    w = np.where(move, to_blocks(np.asarray(weights, dtype=np.float64), factor), 0.0)
    flat = (w.sum(axis=1) <= 0.0)[:, None] & move
    w = np.where(flat, 1.0, w)
    big_t = 1e30
    t = np.where(w > 0, headroom / np.where(w > 0, w, 1.0), big_t)  # the multiplier at which a child stops
    order = np.argsort(t, axis=1)
    t_s = np.take_along_axis(t, order, axis=1)
    h_s = np.take_along_axis(headroom, order, axis=1)
    w_s = np.take_along_axis(w, order, axis=1)
    below = np.zeros_like(h_s)
    below[:, 1:] = np.cumsum(np.where(w_s[:, :-1] > 0, h_s[:, :-1], 0.0), axis=1)
    above = np.cumsum(w_s[:, ::-1], axis=1)[:, ::-1]
    reach = np.where(t_s < big_t, below + t_s * above, np.inf)
    feasible = (reach >= need[:, None]) & (above > 0)
    first = np.where(feasible.any(axis=1), feasible.argmax(axis=1), n)
    rows_index = np.arange(v.shape[0])
    safe = np.minimum(first, n - 1)
    lam = np.where(first < n, (need - below[rows_index, safe]) / np.maximum(above[rows_index, safe], 1e-300), 0.0)
    moved = np.where(first[:, None] < n, np.minimum(headroom, lam[:, None] * w), headroom * (w > 0))
    count = move.sum(axis=1)
    forced = (first >= n) & (count > 0) & (need > 0)
    rest = np.where(forced, need - moved.sum(axis=1), 0.0)
    total_w = w.sum(axis=1)
    lam_f = rest / np.maximum(total_w, 1e-300)
    clipped = np.zeros_like(forced)
    if clip is not None:
        room = np.where(move & (w > 0), (to_blocks(np.asarray(clip, dtype=np.float64), factor) - (v + moved)) /
                        np.where(w > 0, w, 1.0), np.inf)
        allowed = np.maximum(room.min(axis=1), 0.0)
        clipped = forced & up & (lam_f > allowed)
        lam_f = np.where(clipped, allowed, lam_f)
    if clip_low is not None:  # a downward rest stops where a child would pass below clip_low (a hole)
        room = np.where(move & (w > 0), ((v - moved) - to_blocks(np.asarray(clip_low, dtype=np.float64), factor)) /
                        np.where(w > 0, w, 1.0), np.inf)
        allowed = np.maximum(room.min(axis=1), 0.0)
        low = forced & ~up & (lam_f > allowed)
        lam_f = np.where(low, allowed, lam_f)
        clipped = clipped | low
    moved = moved + np.where(forced[:, None], lam_f[:, None] * w, 0.0)
    signed = np.where(up[:, None], moved, -moved)
    signed = np.where(move, signed, 0.0)
    out = from_blocks(v + signed, factor, shape)
    residual = np.abs(target - (v + signed).mean(axis=1))
    stats = {'blocks': int(v.shape[0]), 'blocks_without_moving_child': int(np.count_nonzero(count == 0)),
             'forced_blocks': int(np.count_nonzero(forced)), 'clipped_blocks': int(np.count_nonzero(clipped)),
             'clipped_residual_max_m': float(residual[clipped].max()) if clipped.any() else 0.0,
             'largest_shift_m': float(np.abs(signed).max()) if signed.size else 0.0}
    return out, stats, from_blocks(np.repeat(clipped[:, None], n, axis=1), factor, shape)


def smooth_prestep(values, targets, movable_up, movable_down, land, factor, land_floor_m, sea_ceiling_m, cap, passes,
                   relaxation, weight_floor_m):
    """The multiscale plan section 1 pre-step: `passes` passes adding relaxation x G / M x w, G the Catmull-Rom
    prolongation of the block departures (each over its block's moving share), M the prolongation of the blocks'
    mean weight over their moving children, w each child's weight: its height above the lowest value within one
    parent cell around it plus weight_floor_m (continuous, no block edges). Both prolongations are clamped to their
    bracketing parents, so the ratio stays bounded (dividing before prolonging made a flat block's large ratio land
    on a steep neighbour's large weights: runaway shifts of 15 km). Where the field raises, the movable_up children
    move; where it lowers, the movable_down ones; every child stays within its class bounds. The per-parent
    multipliers so vary smoothly; the exact remainder follows (project_weighted)."""
    values = np.asarray(values, dtype=np.float64)
    lower, upper = class_bounds(values, land, land_floor_m, sea_ceiling_m, cap)
    rows, columns = targets.shape
    for _ in range(passes):
        local = ndimage.minimum_filter(values, size=2 * factor + 1, mode='nearest')
        w = values - local + weight_floor_m
        gap = targets - block_mean(values, factor)
        share = np.where(gap > 0, block_mean(movable_up.astype(float), factor),
                         block_mean(movable_down.astype(float), factor))
        mean_w = np.where(gap > 0, block_mean(w * movable_up, factor), block_mean(w * movable_down, factor))
        has = (share > 0) & (mean_w > 0)
        per_share = np.where(has, gap / np.where(has, share, 1.0), 0.0)
        weight_mean = np.where(has, mean_w / np.where(has, share, 1.0), 0.0)  # mean weight of a moving child
        G = prolong(per_share, factor, 0, 0, rows, columns, limit=True)
        M = prolong(weight_mean, factor, 0, 0, rows, columns, limit=True)
        field = np.where(M > 0, G / np.where(M > 0, M, 1.0), 0.0)
        moving = np.where(field > 0, movable_up, movable_down)
        values = np.clip(values + relaxation * field * w * moving, lower, upper)
    return values


def project_smooth(values, targets, movable_up, movable_down, land, reference, factor, land_floor_m, sea_ceiling_m,
                   cap, clip, passes, relaxation, weight_floor_m, clip_low=None):
    """smooth_prestep, then the exact remainder weighted by the height above each block's lowest reference cell
    (reference_weights, project_weighted). Returns the values, counts and the clipped blocks' mask."""
    smoothed = smooth_prestep(values, targets, movable_up, movable_down, land, factor, land_floor_m, sea_ceiling_m,
                              cap, passes, relaxation, weight_floor_m)
    pre_gap = np.abs(block_mean(smoothed, factor) - targets)
    weights = reference_weights(smoothed, reference, movable_up | movable_down, factor)
    out, stats, clipped = project_weighted(smoothed, targets, movable_up, movable_down, land, weights, factor,
                                           land_floor_m, sea_ceiling_m, cap=cap, clip=clip, clip_low=clip_low)
    moved = (out - values)[movable_up | movable_down]
    stats.update(prestep_departure_rms_m=float(np.sqrt((pre_gap ** 2).mean())),
                 prestep_departure_max_m=float(pre_gap.max()),
                 correction_max_m=float(np.abs(moved).max()) if moved.size else 0.0,
                 correction_rms_m=float(np.sqrt((moved ** 2).mean())) if moved.size else 0.0,
                 correction_mean_m=float(moved.mean()) if moved.size else 0.0)
    return out, stats, clipped


def block_check(values, targets, factor, moving_land, sea, lake, clipped, inner=None):
    """The parent-cell means' departure per block class: clipped, movable land (a moving land child), all sea, all
    lake, other (no moving child: channels, rings, mixed water). Count, max |.| and rms per class."""
    gap = np.abs(block_mean(values, factor) - targets)
    is_clipped = block_mean(clipped.astype(float), factor) > 0
    has_land = (block_mean(moving_land.astype(float), factor) > 0) & ~is_clipped
    all_sea = (block_mean(sea.astype(float), factor) >= 1.0) & ~has_land & ~is_clipped
    all_lake = (block_mean(lake.astype(float), factor) >= 1.0) & ~has_land & ~is_clipped & ~all_sea
    other = ~(is_clipped | has_land | all_sea | all_lake)
    region = np.ones(gap.shape, bool)
    if inner is not None:
        region = np.zeros(gap.shape, bool)
        region[inner] = True
    out = {}
    for name, mask in (('movable_land', has_land), ('all_sea', all_sea), ('all_lake', all_lake), ('clipped', is_clipped),
                       ('other', other)):
        mask = mask & region
        values_m = gap[mask]
        out[name] = {'blocks': int(mask.sum()), 'max_abs_m': float(values_m.max()) if values_m.size else 0.0,
                     'rms_m': float(np.sqrt((values_m ** 2).mean())) if values_m.size else 0.0}
    return out


def one_cell_inlets(sea):
    """Parent sea cells that are one-cell drowned valleys of the coarse D8: at most 3 sea cells among their 8
    neighbours (a one-cell-wide chain, staircase corners included) and not touching the open sea (the sea that
    survives a 3 x 3 morphological opening, dilated by one cell)."""
    sea = np.asarray(sea, bool)
    count = ndimage.convolve(sea.astype(np.int32), np.ones((3, 3), np.int32), mode='constant') - sea
    open_sea = ndimage.binary_opening(sea, structure=np.ones((3, 3)))
    near_open = ndimage.binary_dilation(open_sea, structure=np.ones((3, 3)))
    return sea & (count <= 3) & ~near_open


def departure(values, parent, movable, factor, inner):
    """max |.| and rms of the parent-cell means' departure from the parent, over the inner blocks with a movable child."""
    gap = (block_mean(values, factor) - parent)[inner][block_mean(movable.astype(float), factor)[inner] > 0]
    return {'max_abs': float(np.abs(gap).max()), 'rms': float(np.sqrt((gap ** 2).mean()))}


def local_relief(parent, window):
    """max - min of the parent surface (sea at 0 m) over window x window parent cells."""
    surface = np.maximum(parent, 0.0)
    return ndimage.maximum_filter(surface, size=window, mode='nearest') - \
        ndimage.minimum_filter(surface, size=window, mode='nearest')


# ---------------------------------------------------------------- inflow at a child domain's edge

def row_areas_m2(north_deg, pixel_deg, rows, radius_km):
    """Cell area per row of a window (m^2), as terrain_sandbox_main.cpp's windowmetric computes it."""
    radius_m = radius_km * 1000.0
    pixel = math.radians(pixel_deg)
    north = math.radians(north_deg) - np.arange(rows) * pixel
    return radius_m * radius_m * pixel * (np.sin(north) - np.sin(north - pixel))


def entries(discharge, row0, col0, rows, columns, minimum_m3_s):
    """The parent's rivers entering a child domain: parent cells in the ring just outside the domain's parent cells
    [row0, row0 + rows) x [col0, col0 + columns) whose downstream neighbour (the 8-neighbour with the most discharge,
    when it has more than the cell: the rule river_coverage draws with) lies inside, with at least minimum_m3_s.
    Returns a list of (outside (r, c), inside (r, c), discharge m3/s) in the parent array's indices."""
    found = []
    height, width = discharge.shape
    for r in range(row0 - 1, row0 + rows + 1):
        for c in range(col0 - 1, col0 + columns + 1):
            inside_r, inside_c = row0 <= r < row0 + rows, col0 <= c < col0 + columns
            if (inside_r and inside_c) or not (0 <= r < height and 0 <= c < width):
                continue
            q = float(discharge[r, c])
            if not q >= minimum_m3_s:
                continue
            best, target = q, None
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    rr, cc = r + dr, c + dc
                    if (dr or dc) and 0 <= rr < height and 0 <= cc < width and discharge[rr, cc] > best:
                        best, target = float(discharge[rr, cc]), (rr, cc)
            if target and row0 <= target[0] < row0 + rows and col0 <= target[1] < col0 + columns:
                found.append(((r, c), target, q))
    return found


def inflow_raster(found, row0, col0, factor, shape, surface, candidate, areas_m2):
    """Each entry's discharge as runoff (m/yr) on one child cell: the lowest `candidate` cell (land, movable) of the
    inside parent cell at child index 2 from the domain edge on the side the river comes from (the first cell inside
    the outer sea ring and the pinned ring). Returns the runoff raster (m/yr), the cells' mask and a record per entry
    (child cell or None when the inside parent cell has no candidate there: a river entering at the sea)."""
    rows, columns = shape
    inflow = np.zeros(shape)
    records = []
    for outside, inside, q in found:
        br, bc = (inside[0] - row0) * factor, (inside[1] - col0) * factor
        cells = []
        if outside[0] < row0:
            cells += [(2, j) for j in range(bc, bc + factor)]
        if outside[0] >= row0 + rows // factor:
            cells += [(rows - 3, j) for j in range(bc, bc + factor)]
        if outside[1] < col0:
            cells += [(i, 2) for i in range(br, br + factor)]
        if outside[1] >= col0 + columns // factor:
            cells += [(i, columns - 3) for i in range(br, br + factor)]
        cells = [(i, j) for i, j in cells if candidate[i, j]]
        cell = min(cells, key=lambda ij: surface[ij]) if cells else None
        if cell is not None:
            inflow[cell] += q * SECONDS_PER_YEAR / areas_m2[cell[0]]
        records.append({'outside': list(outside), 'inside': list(inside), 'parent_m3_s': q,
                        'child_cell': list(cell) if cell else None})
    return inflow, inflow > 0, records
