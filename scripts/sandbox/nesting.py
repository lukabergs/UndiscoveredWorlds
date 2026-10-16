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


def project_bounded(values, targets, movable, land, factor, land_floor_m, sea_ceiling_m):
    """Values with their mean over every parent cell set to `targets` exactly (multiscale plan section 3, S2.3's
    "bounds, not weights"): the same shift on every movable child of a parent cell (land and sea alike), each child
    capped by its headroom in the shift's direction: a land child never goes below land_floor_m, a sea child never
    above sea_ceiling_m (upward on land and downward at sea are unbounded), so no child changes class. The shift is
    found by water-filling inside the block (children with less headroom stop at it, the rest share the remainder).
    When the caps exhaust (rare: a block whose land would have to go below the sea), the remainder is spread evenly
    over the movable children ignoring the caps (counted as forced). Blocks without a movable child keep their
    departure. Returns the values and a dict of counts."""
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
    headroom_up = np.where(is_land, big, np.maximum(sea_ceiling_m - v, 0.0))
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
    moved = moved + np.where(move, remainder[:, None], 0.0)
    signed = np.where(up[:, None], moved, -moved)
    signed = np.where(move & (count > 0)[:, None], signed, 0.0)
    out = from_blocks(v + signed, factor, shape)
    stats = {'blocks': int(v.shape[0]), 'blocks_without_movable_child': int(np.count_nonzero(count == 0)),
             'forced_blocks': int(np.count_nonzero(forced)),
             'sea_children_moved': int(np.count_nonzero(move & ~is_land & (np.abs(signed) > 1e-9))),
             'largest_shift_m': float(np.abs(signed).max()) if signed.size else 0.0,
             'shift_p99_m': float(np.percentile(np.abs(signed[move]), 99)) if move.any() else 0.0}
    return out, stats


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
