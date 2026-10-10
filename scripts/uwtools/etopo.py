"""ETOPO 2022 surface elevation: the 60-arc-second globe and the 15-arc-second tiles.

    refs/source/ETOPO/ETOPO_2022_v1_60s_N90W180_surface.tif     the globe, 21600 x 10800
    refs/source/ETOPO/15s/ETOPO_2022_v1_15s_<tile>_surface.tif  15 degree tiles
    refs/source/ETOPO/15s/<file>.receipt.json                   download receipts

A tile is named by its north and west edges (N45E000 covers 30..45 N, 0..15 E). A point
belongs to the tile whose west edge is floor(lon / 15) * 15 and whose north edge is
ceil(lat / 15) * 15; longitude 180 is longitude -180 and the south pole belongs to the
S75 row. Boxes are [west, south, east, north] in degrees; west > east crosses the date
line. Elevations are metres on EGM2008, sea below 0, nodata -99999 read as NaN. Both
sources are read on their own pixel grid (pixel-is-area); nothing is resampled here.
"""
import contextlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ETOPO = ROOT / 'refs/source/ETOPO'
GLOBE_60 = ETOPO / 'ETOPO_2022_v1_60s_N90W180_surface.tif'
TILE_ROOT = ETOPO / '15s'
CATALOGUE = ROOT / 'configs/evaluation/earth-windows.json'
CATALOGUE_SCHEMA = 'uw-earth-windows/v1'
TILE_DEG = 15
TILE_URL = ('https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/'
            '15s_surface_elev_gtif/ETOPO_2022_v1_15s_{name}_surface.tif')
NODATA = -99999.0
SOURCES = ('etopo60', 'etopo15')
# Box edges within this fraction of a pixel of a grid line snap to it.
SNAP = 1e-6


def tile_name(lon, lat):
    """Name of the 15 degree tile holding (lon, lat), e.g. N45E000."""
    if not -90 <= lat <= 90:
        raise ValueError(f'Latitude out of range: {lat}')
    north = max(int(math.ceil(lat / TILE_DEG)) * TILE_DEG, -90 + TILE_DEG)
    west = int(math.floor(_wrap(lon) / TILE_DEG)) * TILE_DEG
    return _name(west, north)


def _name(west, north):
    return (f"{'N' if north >= 0 else 'S'}{abs(north):02d}"
            f"{'E' if west >= 0 else 'W'}{abs(west):03d}")


def _wrap(lon):
    """Longitude in [-180, 180)."""
    return (lon + 180.0) % 360.0 - 180.0


def tile_bounds(name):
    """(west, south, east, north) of a tile name."""
    if (len(name) != 7 or name[0] not in 'NS' or name[3] not in 'EW'
            or not name[1:3].isdigit() or not name[4:7].isdigit()):
        raise ValueError(f'Not a tile name: {name}')
    north = int(name[1:3]) * (1 if name[0] == 'N' else -1)
    west = int(name[4:7]) * (1 if name[3] == 'E' else -1)
    if north % TILE_DEG or west % TILE_DEG or not (-75 <= north <= 90 and -180 <= west <= 165):
        raise ValueError(f'Not a tile name: {name}')
    return west, north - TILE_DEG, west + TILE_DEG, north


def check_box(west, south, east, north):
    if not (-180 <= west <= 180 and -180 <= east <= 180):
        raise ValueError(f'Longitude out of range: {west}, {east}')
    if not -90 <= south < north <= 90:
        raise ValueError(f'Latitudes must satisfy -90 <= south < north <= 90: {south}, {north}')
    if west == east:
        raise ValueError('A box needs west != east (use -180, 180 for the whole circle)')


def _spans(west, east):
    return [(west, east)] if west < east else [(west, 180), (-180, east)]


def tiles_for_box(west, south, east, north):
    """Tiles whose interior meets the box, north row first, west to east (across the date line
    the part west of it first)."""
    check_box(west, south, east, north)
    names = []
    for top in range(90, -90, -TILE_DEG):
        if not (top - TILE_DEG < north and top > south):
            continue
        for a, b in _spans(west, east):
            for left in range(-180, 180, TILE_DEG):
                if left < b and left + TILE_DEG > a:
                    name = _name(left, top)
                    if name not in names:
                        names.append(name)
    return names


def tile_url(name):
    tile_bounds(name)
    return TILE_URL.format(name=name)


def tile_path(name, root=None):
    tile_bounds(name)
    return Path(root or TILE_ROOT) / f'ETOPO_2022_v1_15s_{name}_surface.tif'


def load_catalogue(path=None):
    data = json.loads(Path(path or CATALOGUE).read_text(encoding='utf-8'))
    if data.get('schema') != CATALOGUE_SCHEMA:
        raise ValueError(f'Catalogue schema is not {CATALOGUE_SCHEMA}')
    return data['windows']


def box_contains(box, lon, lat):
    west, south, east, north = box
    inside_lon = west <= lon <= east if west < east else (lon >= west or lon <= east)
    return inside_lon and south <= lat <= north


def _grid(dataset):
    """Pixel size and the grid's phase: the grid line nearest -180 and the one nearest 90."""
    a, b, c, d, e, f = tuple(dataset.transform)[:6]
    if b or d or not math.isclose(-e, a, rel_tol=1e-9) or a <= 0:
        raise ValueError(f'{dataset.name}: not a north-up grid with square pixels')
    phase_x = c - a * round((c + 180.0) / a)
    phase_y = f - a * round((f - 90.0) / a)
    return a, phase_x, phase_y


def read_box(west, south, east, north, source):
    """Elevation of the box from 'etopo60' (the globe) or 'etopo15' (the tiles).

    Returns (elevation, extent): float32 rows north first, nodata as NaN, on the source's own
    pixels, the box snapped outwards to them; extent {west, south, east, north, pixel_deg}
    of the returned pixels' outer edges (west > east when they cross the date line).
    """
    import rasterio
    from rasterio.windows import Window

    check_box(west, south, east, north)
    if source == 'etopo60':
        paths = [GLOBE_60]
    elif source == 'etopo15':
        paths = [tile_path(name) for name in tiles_for_box(west, south, east, north)]
    else:
        raise ValueError(f'Unknown source {source!r}; expected one of {SOURCES}')
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError('ETOPO data not downloaded (scripts/refs/fetch_etopo_tiles.py): '
                                + ', '.join(missing))
    with contextlib.ExitStack() as stack:
        datasets = [stack.enter_context(rasterio.open(path)) for path in paths]
        pixel, phase_x, phase_y = _grid(datasets[0])
        for dataset in datasets[1:]:
            other = _grid(dataset)
            if not all(math.isclose(u, v, rel_tol=0, abs_tol=pixel * SNAP)
                       for u, v in zip(other, (pixel, phase_x, phase_y))):
                raise ValueError(f'{dataset.name}: pixel grid differs from {datasets[0].name}')
        circle = int(round(360.0 / pixel))
        east_unwrapped = east if east > west else east + 360.0
        i0 = math.floor((west - phase_x) / pixel + SNAP)
        i1 = math.ceil((east_unwrapped - phase_x) / pixel - SNAP)
        j0 = math.floor((phase_y - north) / pixel + SNAP)
        j1 = math.ceil((phase_y - south) / pixel - SNAP)
        placed = []
        for dataset in datasets:
            col0 = int(round((dataset.transform.c - phase_x) / pixel))
            row0 = int(round((phase_y - dataset.transform.f) / pixel))
            placed.append((dataset, col0, row0))
        j0 = max(j0, min(row0 for _, _, row0 in placed))
        j1 = min(j1, max(row0 + d.height for d, _, row0 in placed))
        elevation = np.full((j1 - j0, i1 - i0), np.nan, dtype=np.float32)
        for dataset, col0, row0 in placed:
            r0, r1 = max(j0, row0), min(j1, row0 + dataset.height)
            if r0 >= r1:
                continue
            for shift in (0, circle):
                c0 = max(i0, col0 + shift)
                c1 = min(i1, col0 + shift + dataset.width)
                if c0 >= c1:
                    continue
                window = Window(c0 - col0 - shift, r0 - row0, c1 - c0, r1 - r0)
                data = dataset.read(1, window=window).astype(np.float32)
                nodata = dataset.nodata if dataset.nodata is not None else NODATA
                data[(data == nodata) | (data <= NODATA)] = np.nan
                elevation[r0 - j0:r1 - j0, c0 - i0:c1 - i0] = data
    extent_west = _wrap(phase_x + i0 * pixel)
    extent_east = phase_x + i1 * pixel
    if extent_east > 180.0 + pixel * SNAP:
        extent_east -= 360.0
    extent = dict(west=round(extent_west, 10), south=round(phase_y - j1 * pixel, 10),
                  east=round(extent_east, 10), north=round(phase_y - j0 * pixel, 10),
                  pixel_deg=pixel)
    return elevation, extent
