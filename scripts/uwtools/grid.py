"""Cell-centred W x W/2 grids: centre latitudes and longitudes, area weights.

Rows run north to south, columns east from -180 degrees (docs/architecture/world-model.md).
Two weight forms are in use and are not numerically identical: the exact band area
sin(north) - sin(south) and the cosine of the centre latitude. The band edges are computed
in radians or in degrees by different callers; the two differ in the last bits.
"""
import numpy as np


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
