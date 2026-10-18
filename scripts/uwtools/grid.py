"""Cell-centred W x W/2 grids: centre latitudes and longitudes, area weights.

Rows run north to south, columns east from -180 degrees (docs/architecture/world-model.md).
Two weight forms are in use and are not numerically identical: the exact band area
sin(north) - sin(south) and the cosine of the centre latitude. The band edges are computed
in radians or in degrees by different callers; the two differ in the last bits.
remap_area_mean is the validation summary's conservative reduction (also used for region masks).
"""
import math

import numpy as np

EARTH_RADIUS_KM = 6371.0


def cell_km(scenario, columns, radius_km=None):
    """East-west size in km of one of `columns` cells at the equator. The torus engine's own metric,
    run_scenario.world_width_km / columns, while its scenario records it (so the lat-lon numbers are
    unchanged); without it the sphere's 2 pi R / columns, R the grid's radius (6371 km when absent)."""
    width_km = (scenario or {}).get('world_width_km')
    if width_km is not None:
        return float(width_km) / columns
    return 2.0 * math.pi * float(radius_km or EARTH_RADIUS_KM) / columns


def latitudes(rows):
    """Cell-centre latitudes in degrees, north row first."""
    return 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows


def longitudes(columns):
    return -180.0 + (np.arange(columns) + 0.5) * 360.0 / columns


def cosine_latitude(rows):
    """Cosine of each row's centre latitude."""
    return np.cos(np.deg2rad(latitudes(rows)))


def band_area(rows, degree_edges=False):
    """sin(north edge) - sin(south edge) of each row: its area on the unit sphere per radian of longitude.

    Edges come from linspace in radians, or with degree_edges from linspace in degrees converted to radians.
    """
    if degree_edges:
        edges = np.sin(np.deg2rad(np.linspace(90.0, -90.0, rows + 1)))
    else:
        edges = np.sin(np.linspace(np.pi / 2, -np.pi / 2, rows + 1))
    return edges[:-1] - edges[1:]


def area_weights(rows, columns):
    """Band areas broadcast to rows x columns (read-only view, not normalised)."""
    return np.broadcast_to(band_area(rows)[:, None], (rows, columns))


def remap_area_mean(a, ny, nx):
    """Area mean of the finite source values over each target cell of a cell-centred ny x nx grid.

    Latitude overlaps are exact band areas with degree edges, longitude overlaps periodic
    fractions; target cells without finite source values are NaN. a is 2-D or a stack of 2-D fields.
    """
    sy, sx = a.shape[-2:]
    sb, tb = np.sin(np.deg2rad(np.linspace(90, -90, sy + 1))), np.sin(np.deg2rad(np.linspace(90, -90, ny + 1)))
    lat = np.maximum(0, np.minimum(tb[:-1, None], sb[None, :-1]) - np.maximum(tb[1:, None], sb[None, 1:]))
    tx, ss = np.linspace(0, 1, nx + 1), np.linspace(0, 1, sx + 1)
    lon = np.maximum(0, np.minimum(tx[1:, None], ss[None, 1:]) - np.maximum(tx[:-1, None], ss[None, :-1]))
    def one(v):
        finite = np.isfinite(v)
        if sx % nx == 0:
            k = sx // nx
            num = np.where(finite, v, 0).reshape(sy, nx, k).sum(2)
            den = finite.reshape(sy, nx, k).sum(2)
            num, den = lat @ num, lat @ den
        else:
            num, den = lat @ np.where(finite, v, 0) @ lon.T, lat @ finite.astype(float) @ lon.T
        return np.divide(num, den, out=np.full((ny, nx), np.nan), where=den > 0)
    if a.ndim == 2:
        return one(a)
    return np.stack([one(v) for v in a])
