"""Copernicus DEM GLO-30 (30 m) tiles: names, paths and a mosaicking box reader.

    refs/source/copernicus-glo30/<tile>.tif                one degree tiles (Cloud Optimized GeoTIFF)
    refs/source/copernicus-glo30/<tile>.tif.receipt.json   download receipts
    refs/source/copernicus-glo30/<tile>.tif.absent.json    tiles the provider does not ship (all sea)

A tile is named by its south-west corner: Copernicus_DSM_COG_10_N61_00_E006_00_DEM covers
61..62 N, 6..7 E. Rows are 1 arc-second; the east-west spacing widens with latitude
(3600 columns up to 50 degrees, 2400 to 60, 1800 to 70, 1200 to 80, 720 to 85, 360 to 90),
and each tile is read by its own transform. Pixels are point-registered on whole arc-second
lines (the transform's origin is half a pixel west and north of the corner), so neighbouring
tiles abut without overlap. Elevations are metres on EGM2008; the sea is 0 inside a tile.

Copernicus DEM GLO-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018
provided under COPERNICUS by the European Union and ESA; all rights reserved.
"""
import contextlib
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TILE_ROOT = ROOT / 'refs/source/copernicus-glo30'
URL = 'https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif'
SOURCE = ('Copernicus DEM GLO-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH '
          '2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved')
# Areas fetched for the 100 m nested test (terrain roadmap T1.5): west, south, east, north.
AREAS = {
    'sognefjord': (4, 60, 9, 62.5),
    'patagonia': (-76, -52, -72, -47),
}
# Box edges within this fraction of a pixel of a pixel edge snap to it.
SNAP = 1e-6
# Wider than half the widest pixel (10 arc-seconds): how far read_box looks for neighbours.
MARGIN_DEG = 10.0 / 3600.0


def tile_name(lat, lon):
    """Name of the tile holding (lat, lon), from the floor of both."""
    south, west = int(math.floor(lat)), int(math.floor(lon))
    if not (-90 <= south < 90 and -180 <= west < 180):
        raise ValueError(f'Outside the globe: {lat}, {lon}')
    return (f"Copernicus_DSM_COG_10_{'N' if south >= 0 else 'S'}{abs(south):02d}_00_"
            f"{'E' if west >= 0 else 'W'}{abs(west):03d}_00_DEM")


def tile_url(name):
    return URL.format(name=name)


def tile_path(lat, lon, root=None):
    return Path(root or TILE_ROOT) / f'{tile_name(lat, lon)}.tif'


def absent_path(path):
    path = Path(path)
    return path.with_name(path.name + '.absent.json')


def check_box(west, south, east, north):
    if not (-180 <= west < east <= 180):
        raise ValueError(f'Longitudes must satisfy -180 <= west < east <= 180: {west}, {east}')
    if not -90 <= south < north <= 90:
        raise ValueError(f'Latitudes must satisfy -90 <= south < north <= 90: {south}, {north}')


def tiles_for_box(west, south, east, north):
    """(lat, lon) south-west corners of the tiles whose interior meets the box, north row first,
    west to east."""
    check_box(west, south, east, north)
    return [(lat, lon)
            for lat in range(math.ceil(north) - 1, math.floor(south) - 1, -1)
            for lon in range(math.floor(west), math.ceil(east))]


def read_box(west, south, east, north, root=None):
    """Elevation of the box mosaicked from the tiles.

    Returns (elevation, extent): float32 rows north first, nodata and tiles recorded as absent
    as NaN. Rows are the tiles' own 1 arc-second rows; columns use the finest east-west spacing
    among the tiles met, and a tile with a wider spacing (a more poleward band) is sampled by
    nearest pixel. The box snaps outwards to whole pixels; extent {west, south, east, north,
    pixel_deg_x, pixel_deg_y} gives the returned pixels' outer edges. A tile meeting the box
    that is neither downloaded nor recorded absent raises FileNotFoundError; a neighbour that
    only lends the snapped edge row or column is read when present and is NaN otherwise.
    """
    import rasterio
    from rasterio.windows import Window

    check_box(west, south, east, north)
    required = set(tiles_for_box(west, south, east, north))
    # Snapping outwards can take a pixel row or column from a neighbouring tile (a box edge on
    # a whole degree is a pixel centre); such neighbours are read when present.
    margin = MARGIN_DEG
    paths, missing = [], []
    for lat, lon in tiles_for_box(max(west - margin, -180), max(south - margin, -90),
                                  min(east + margin, 180), min(north + margin, 90)):
        path = tile_path(lat, lon, root)
        if path.is_file():
            paths.append(path)
        elif (lat, lon) in required and not absent_path(path).is_file():
            missing.append(str(path))
    if missing:
        raise FileNotFoundError('GLO-30 tiles not downloaded (scripts/refs/fetch_glo30.py): '
                                + ', '.join(missing))
    with contextlib.ExitStack() as stack:
        datasets = [stack.enter_context(rasterio.open(path)) for path in paths]
        for dataset in datasets:
            t = dataset.transform
            if t.b or t.d or t.a <= 0 or t.e >= 0:
                raise ValueError(f'{dataset.name}: not a north-up grid')
        dy = 1.0 / 3600.0
        dx = min((d.transform.a for d in datasets), default=dy)
        # Pixel edges: whole-degree corners shifted half a pixel west and north.
        phase_x, phase_y = -dx / 2, dy / 2
        i0 = math.floor((west - phase_x) / dx + SNAP)
        i1 = math.ceil((east - phase_x) / dx - SNAP)
        # Row edges are the latitudes phase_y + k * dy; output row r has its north edge at
        # north_edge - r * dy.
        top = math.ceil((north - phase_y) / dy - SNAP)
        bottom = math.floor((south - phase_y) / dy + SNAP)
        rows = top - bottom
        elevation = np.full((rows, i1 - i0), np.nan, dtype=np.float32)
        north_edge = phase_y + top * dy
        centre_lon = phase_x + (np.arange(i0, i1) + 0.5) * dx
        for dataset in datasets:
            t = dataset.transform
            if not math.isclose(-t.e, dy, rel_tol=1e-9):
                raise ValueError(f'{dataset.name}: rows are not 1 arc-second')
            # Output rows inside the tile: their centres lie between the tile's edges.
            r0 = max(0, int(round((north_edge - t.f) / dy)))
            r1 = min(rows, int(round((north_edge - (t.f - dataset.height * dy)) / dy)))
            if r0 >= r1:
                continue
            columns = np.floor((centre_lon - t.c) / t.a).astype(np.int64)
            inside = (columns >= 0) & (columns < dataset.width)
            if not inside.any():
                continue
            c_lo, c_hi = int(columns[inside].min()), int(columns[inside].max()) + 1
            row_lo = int(round((t.f - north_edge) / dy)) + r0
            window = Window(c_lo, row_lo, c_hi - c_lo, r1 - r0)
            data = dataset.read(1, window=window).astype(np.float32)
            if dataset.nodata is not None:
                data[data == dataset.nodata] = np.nan
            elevation[r0:r1, inside] = data[:, columns[inside] - c_lo]
    extent = dict(west=round(phase_x + i0 * dx, 10), south=round(phase_y + bottom * dy, 10),
                  east=round(phase_x + i1 * dx, 10), north=round(north_edge, 10),
                  pixel_deg_x=dx, pixel_deg_y=dy)
    return elevation, extent
