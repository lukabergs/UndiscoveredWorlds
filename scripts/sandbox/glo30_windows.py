"""Copernicus GLO-30 reference windows at a sandbox window's pixels (terrain roadmap T1.5): the Earth drawn at the
same pixel size and extent in kilometres as a nested sandbox window, so a 100 m level can be judged beside it.

The target grid copies a window's geometry (its <domain>-grid.json, written by nested_window.py, or explicit
numbers): rows x columns pixels of pixel_ns_km north-south and pixel_ew_km east-west, centred on the given
longitude and latitude (a pixel corner). The pixels are square in neither degrees nor km: north-south they are
pixel_ns_km / R radians; east-west pixel_ew_km / (R cos(centre latitude)) radians, so every pixel has the window's
km size at the centre latitude (the sandbox window's pixels are square in degrees: their east-west km size is the
centre latitude's, as here). R is the Earth's 6371 km.

Elevation: GLO-30 (uwtools.glo30.read_box: 1 arc-second rows, the tiles' own east-west spacing, 1 to 2 arc-seconds
here) box-averaged to the target pixels, area-weighted in degrees (fractional source pixels, the cumulative-sum
method of evaluation/earth_windows.py, with separate east-west and north-south source spacings); absent tiles (all
sea) and nodata are NaN and left out of the means. GLO-30 has no bathymetry (the sea is 0 m, about 130 shore pixels
below 0 are artefacts), so the sea comes from ETOPO 2022 at 15 arc-seconds, as earth_atlas.py does: ETOPO interpolated
bilinearly at the pixel centres (its 460 m pixels are coarser than the target's: a box mean would copy its blocks), its
ocean the pixels below 0 m whose 60 arc-second cell is in the world-ocean mask
(earth_windows.global_ocean) plus the below-0 m pixels 4-connected to them; then the sea is every pixel whose GLO-30
mean is at or below 0 m (or NaN) 4-connected to that ocean, so the fjords are sea up to their heads. Sea pixels take
ETOPO's depth (at most -1 m); land keeps GLO-30. No lakes, rivers or ice are drawn (GLO-30's lakes are flat land).

Writes refs/processed/earth-atlas/glo30-<name>-<km>km.npy (elevation, float32, rows north first),
glo30-<name>-<km>km-sea.npy (bool), glo30-<name>-<km>km.json (receipt: box, grid, method, source tiles) and
glo30-<name>-<km>km-{elevation,relief}.png (the sandbox's colours and hillshade, the window's exaggeration, a scale
bar), so compare_runs.py --earth glo30-<name>-<km>km shows them. <km> is the north-south extent, rounded.

  uv run python scripts/sandbox/glo30_windows.py --name sognefjord --centre 7.0,61.2
      --like out/scratch/sandbox/t15-nest-a/G-coast-nest-grid.json
  uv run python scripts/sandbox/glo30_windows.py --name x --centre LON,LAT --rows 960 --columns 1590
      --pixel-km 0.1018,0.0637 [--exaggeration 2]
"""
import argparse
import datetime
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'evaluation'))
sys.path.insert(0, str(SCRIPTS / 'sandbox'))
from uwtools import etopo, glo30  # noqa: E402
from uwtools.map_rendering import save_rgb  # noqa: E402
from maps import hillshade, shade_factor, terrain_base  # noqa: E402
import earth_windows as ew  # noqa: E402
import overlays  # noqa: E402

ATLAS = ROOT / 'refs' / 'processed' / 'earth-atlas'
RADIUS_M = ew.EARTH_RADIUS_M


def target_edges(longitude, latitude, rows, columns, pixel_ns_km, pixel_ew_km):
    """(longitude edges west to east, latitude edges north to south, dx_deg, dy_deg) of the target grid."""
    dy = math.degrees(pixel_ns_km * 1000.0 / RADIUS_M)
    dx = math.degrees(pixel_ew_km * 1000.0 / (RADIUS_M * math.cos(math.radians(latitude))))
    lon = longitude + (np.arange(columns + 1) - columns / 2.0) * dx
    lat = latitude + (rows / 2.0 - np.arange(rows + 1)) * dy
    return lon, lat, dx, dy


def box_average(values, west, north, dx_source, dy_source, lon_edges, lat_edges):
    """Area-weighted (in degrees) means of a source raster with east-west spacing dx_source and north-south
    dy_source (outer edges west, north) over target pixels between the edges; NaN source pixels left out."""
    values = np.asarray(values, dtype=float)
    rows, columns = values.shape
    x = np.clip((lon_edges - west) / dx_source, 0.0, columns)
    y = np.clip((north - lat_edges) / dy_source, 0.0, rows)
    finite = np.isfinite(values)
    sums = ew._integrated(ew._integrated(np.where(finite, values, 0.0), x, 1), y, 0)
    weights = ew._integrated(ew._integrated(finite.astype(float), x, 1), y, 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(weights > 1e-6, sums / np.where(weights > 0, weights, 1.0), np.nan)


def bilinear(values, extent, lon_edges, lat_edges):
    """A source coarser than the target (ETOPO's 15 arc-seconds under 100 m pixels) interpolated bilinearly at the
    pixel centres: a box mean would copy its blocks (stair steps in the shaded sea)."""
    step = float(extent['pixel_deg'])
    centres_lon = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    centres_lat = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    column = (centres_lon - float(extent['west'])) / step - 0.5
    row = (float(extent['north']) - centres_lat) / step - 0.5
    rows, columns = np.meshgrid(row, column, indexing='ij')
    filled = np.where(np.isfinite(values), values, 0.0)
    return ndimage.map_coordinates(filled, [rows, columns], order=1, mode='nearest')


def etopo_ocean(etopo_mean, lon_edges, lat_edges):
    """ETOPO's ocean on the target pixels: below 0 m and in the world-ocean mask's 60 arc-second cell, plus the
    below-0 m pixels 4-connected to those (earth_atlas.window's rule for the 15 arc-second windows)."""
    mask, extent = ew.global_ocean()
    step = float(extent['pixel_deg'])
    centres_lon = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    centres_lat = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    columns = np.floor(((centres_lon - float(extent['west'])) % 360.0) / step).astype(int) % mask.shape[1]
    rows = np.clip(np.floor((float(extent['north']) - centres_lat) / step).astype(int), 0, mask.shape[0] - 1)
    below = etopo_mean < 0.0
    sea = below & (np.asarray(mask[np.ix_(rows, columns)]) != 0)
    labels, _ = ndimage.label(below)
    return np.isin(labels, np.unique(labels[sea & (labels > 0)])) & (labels > 0)


def reference_window(longitude, latitude, rows, columns, pixel_ns_km, pixel_ew_km):
    """(elevation m, sea mask, receipt dict) of the GLO-30 window on the target grid."""
    lon, lat, dx, dy = target_edges(longitude, latitude, rows, columns, pixel_ns_km, pixel_ew_km)
    pad = 4.0 / 3600.0  # a few source pixels beyond the edges, so edge pixels are whole
    west, east, south, north = lon[0] - pad, lon[-1] + pad, lat[-1] - pad, lat[0] + pad
    source, extent = glo30.read_box(west, south, east, north)
    elevation = box_average(source, extent['west'], extent['north'], extent['pixel_deg_x'], extent['pixel_deg_y'],
                            lon, lat)
    del source
    pad15 = 30.0 / 3600.0
    array, extent15 = etopo.read_box(west - pad15, south - pad15, east + pad15, north + pad15, 'etopo15')
    etopo_mean = bilinear(array, extent15, lon, lat)
    del array
    ocean = etopo_ocean(etopo_mean, lon, lat)
    candidate = ~np.isfinite(elevation) | (elevation <= 0.0)
    labels, _ = ndimage.label(candidate)
    sea = np.isin(labels, np.unique(labels[ocean & candidate & (labels > 0)])) & (labels > 0)
    shore_artefacts = int(np.count_nonzero(np.isfinite(elevation) & (elevation < 0.0) & ~sea))
    out = np.where(sea, np.minimum(np.nan_to_num(etopo_mean, nan=-1.0), -1.0), elevation)
    out = np.where(np.isfinite(out), out, 0.0)  # an absent pixel off the sea (none expected) as 0 m land
    glo_files = sorted({glo30.tile_path(lat_, lon_).name for lat_, lon_ in glo30.tiles_for_box(west, south, east, north)})
    receipt = {
        'schema': 'uw-glo30-window/v1',
        'centre': [longitude, latitude], 'rows': rows, 'columns': columns,
        'pixel_ns_km': pixel_ns_km, 'pixel_ew_km': pixel_ew_km,
        'pixel_deg_x': dx, 'pixel_deg_y': dy,
        'box': {'west': float(lon[0]), 'south': float(lat[-1]), 'east': float(lon[-1]), 'north': float(lat[0])},
        'extent_km': [round(rows * pixel_ns_km, 2), round(columns * pixel_ew_km, 2)],
        'method': ('GLO-30 box-averaged to the pixels, area-weighted in degrees (fractional source pixels, NaN left '
                   'out); sea = pixels at or below 0 m (or NaN) 4-connected to ETOPO 2022 15" ocean (below 0 m in '
                   'the 60" world-ocean mask, plus below-0 m pixels 4-connected to it), with ETOPO 15" bilinear '
                   'depth (at most -1 m)'),
        'glo30_source': {'tiles': glo_files, 'root': 'refs/source/copernicus-glo30',
                         'pixel_deg_x': extent['pixel_deg_x'], 'pixel_deg_y': extent['pixel_deg_y'],
                         'licence': glo30.SOURCE},
        'etopo_source': {'tiles': [p.name for p in ew.etopo_files('etopo15', {'west': west, 'south': south,
                                                                               'east': east, 'north': north})]},
        'sea_share': round(float(sea.mean()), 4),
        'below_zero_land_pixels': shore_artefacts,
        'code': {'module': 'scripts/sandbox/glo30_windows.py', 'revision': ew.code_revision()},
        'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
    }
    return out.astype(np.float32), sea, receipt


def save_maps(directory, stem, elevation, sea, pixel_ns_km, pixel_ew_km, exaggeration, note, lake_surface=None):
    """<stem>-elevation.png and <stem>-relief.png with a scale bar (km per pixel east-west)."""
    shown = elevation if lake_surface is None else np.where(np.isfinite(lake_surface), lake_surface, elevation)
    base = terrain_base(shown, sea)
    relief = np.clip(base * shade_factor(hillshade(shown, pixel_ns_km * 1000.0, pixel_ew_km * 1000.0,
                                                   exaggeration))[..., None], 0, 255)
    if lake_surface is not None:
        lake = np.isfinite(lake_surface)
        base[lake] = (30.0, 80.0, 160.0)
        relief[lake] = (30.0, 80.0, 160.0)
    for kind, rgb in (('elevation', base), ('relief', relief)):
        save_rgb(directory / f'{stem}-{kind}.png', overlays.scale_bar(rgb, pixel_ew_km, note))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', required=True)
    parser.add_argument('--centre', required=True, help='LON,LAT (degrees, east and north positive)')
    parser.add_argument('--like', help="a nested window's <domain>-grid.json (rows, columns, pixel km)")
    parser.add_argument('--rows', type=int)
    parser.add_argument('--columns', type=int)
    parser.add_argument('--pixel-km', help='NS,EW km per pixel')
    parser.add_argument('--exaggeration', type=float, default=None)
    parser.add_argument('--out', default=str(ATLAS))
    arguments = parser.parse_args()
    longitude, latitude = (float(value) for value in arguments.centre.split(','))
    exaggeration = 2.0
    if arguments.like:
        grid = json.loads(Path(arguments.like).read_text(encoding='utf-8'))
        rows, columns = grid['rows'], grid['columns']
        pixel_ns_km, pixel_ew_km = grid['pixel_ns_km'], grid['pixel_ew_km']
        exaggeration = grid.get('exaggeration', exaggeration)
    else:
        rows, columns = arguments.rows, arguments.columns
        pixel_ns_km, pixel_ew_km = (float(value) for value in arguments.pixel_km.split(','))
    if arguments.exaggeration is not None:
        exaggeration = arguments.exaggeration
    started = time.perf_counter()
    elevation, sea, receipt = reference_window(longitude, latitude, rows, columns, pixel_ns_km, pixel_ew_km)
    receipt['exaggeration'] = exaggeration
    receipt['like'] = arguments.like
    stem = f'glo30-{arguments.name}-{round(rows * pixel_ns_km):d}km'
    out = Path(arguments.out)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / f'{stem}.npy', elevation)
    np.save(out / f'{stem}-sea.npy', sea)
    hemisphere = 'N' if latitude >= 0 else 'S'
    save_maps(out, stem, elevation, sea, pixel_ns_km, pixel_ew_km, exaggeration,
              f'(GLO-30 {arguments.name}, {abs(latitude):.1f} {hemisphere}; {pixel_ns_km * 1000:.0f} x '
              f'{pixel_ew_km * 1000:.0f} m pixels)')
    receipt['seconds'] = round(time.perf_counter() - started, 1)
    (out / f'{stem}.json').write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    print(f'{stem}: sea {receipt["sea_share"]:.3f}, {receipt["seconds"]} s -> {out / stem}.npy')


if __name__ == '__main__':
    main()
