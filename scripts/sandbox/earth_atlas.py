"""Earth reference atlas for the terrain sandbox: the Earth drawn once at the sandbox's scales, so its maps
can be compared with what the Earth shows without rendering the Earth on every run.

Three groups, each item an elevation map and a relief map (the sandbox's colours and light: terrain_base on
tectonics:elevation, times shade_factor(hillshade) of scripts/evaluation/maps.py, east-west pixel size with
cos(latitude) clamped to at least 0.25, no outline):
  1. the whole Earth at 128 to 2048 columns (W x W/2, cell-centred, rows north first, columns from 180 W),
     box-averaged from the 60 arc-second globe; vertical exaggeration 8;
  2. regional windows at the sandbox's window zoom (180 / 128 / 64 degrees per pixel, about 2.44 km),
     1280 x 1024 pixels, from the 60 arc-second globe; exaggeration 3;
  3. two closer windows at k 256 (about 610 m), 1024 x 1024 pixels, from the 15 arc-second tiles
     (scripts/refs/fetch_etopo_tiles.py); exaggeration 2.
Elevations are read and box-averaged (area-weighted in degrees) with scripts/evaluation/earth_windows.py; a
window's grid is anchored at its centre (a pixel corner). Sea is the world ocean (earth_windows.global_ocean,
its cache under refs/processed/earth-windows/ocean-60s/); below-0 m basins off it are land. On the 15
arc-second windows, whose pixels are finer than the mask's cells, the below-0 m pixels 4-connected to a pixel
the mask calls ocean are sea too, so the fjords are sea.

  uv run python scripts/sandbox/earth_atlas.py [--out refs/processed/earth-atlas]
Writes the PNGs, manifest.json and index.html (local-only under refs/).
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
from uwtools import etopo  # noqa: E402
from uwtools.map_rendering import save_rgb  # noqa: E402
from maps import hillshade, shade_factor, terrain_base  # noqa: E402
import earth_windows as ew  # noqa: E402
import overlays  # noqa: E402  (scripts/sandbox, the script's own directory)

RADIUS_M = ew.EARTH_RADIUS_M
GLOBE_COLUMNS = (128, 256, 512, 1024, 2048)
GLOBE_EXAGGERATION = 8.0
BAND_DEG = 22.5  # the globe is read in latitude bands of this height (a whole number of pixels at every width)

WINDOW_DEG = 180.0 / 128 / 64  # a 256-column world at k 64
WINDOW_PIXELS = (1280, 1024)
WINDOW_EXAGGERATION = 3.0
WINDOWS = [  # (name, title, longitude, latitude)
    ('andes', 'Central Andes and the Bolivian orocline', -67.0, -18.0),
    ('himalaya', 'Himalaya, Tibetan Plateau and Hengduan Shan', 92.0, 31.0),
    ('norway', 'Norwegian fjords', 12.0, 64.0),
    ('patagonia', 'Patagonian fjords and icefields', -72.0, -49.0),
    ('zagros', 'Zagros fold belt', 52.0, 30.0),
    ('appalachians', 'Appalachian valley and ridge', -78.0, 39.0),
    ('alps', 'The Alps', 10.0, 46.0),
    ('basin-range', 'Basin and Range with the Colorado Plateau', -113.0, 38.0),
    ('east-africa', 'East African Rift', 36.0, -2.0),
]

CLOSE_DEG = 180.0 / 128 / 256  # k 256
CLOSE_PIXELS = (1024, 1024)
CLOSE_EXAGGERATION = 2.0
CLOSE = [
    ('sognefjord', 'Norwegian fjords: Sognefjord', 6.5, 61.5),
    ('everest', 'Himalaya at Everest', 87.0, 28.0),
]

SOURCE_TEXT = {'etopo60': 'ETOPO 2022 60 arc-second globe (about 1.85 km)',
               'etopo15': 'ETOPO 2022 15 arc-second tiles (about 460 m)'}


def relative(path):
    path = Path(path)
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def relief_rgb(elevation, sea, north, pixel_deg, exaggeration):
    """Elevation colours times the hillshade, as the sandbox's relief_rgb (cos latitude at least 0.25)."""
    rows = elevation.shape[0]
    dy = RADIUS_M * math.radians(pixel_deg)
    latitude = north - (np.arange(rows) + 0.5) * pixel_deg
    dx = dy * np.maximum(np.cos(np.radians(latitude)), 0.25)
    shade = shade_factor(hillshade(elevation, dy, dx, exaggeration))
    return np.clip(terrain_base(elevation, sea) * shade[..., None], 0, 255)


def to_uint8(rgb):
    return np.clip(np.rint(rgb), 0, 255).astype(np.uint8)


WORLD_CELL_DEG = 360.0 / 256.0  # the sandbox's world grid: the cells its windows are refined from


def annotate(rgb, north, pixel_deg, west=None, centre_latitude=0.0):
    """The sandbox's scale references (overlays.py): on a window (west given) the faint grid of the
    256-column world's cells, at their true longitudes and latitudes; on every image a scale bar measured
    east-west at the centre latitude."""
    if west is not None:
        step = WORLD_CELL_DEG / pixel_deg
        rgb = overlays.grid(rgb, step, ((-180.0 - west) / pixel_deg) % step, ((north - 90.0) / pixel_deg) % step)
        note = (f'(east-west, {abs(centre_latitude):.0f} {"N" if centre_latitude >= 0 else "S"}; '
                'grid: 256-column world cells)')
    else:
        note = '(east-west, at the equator)'
    km_per_px = RADIUS_M * math.radians(pixel_deg) * math.cos(math.radians(centre_latitude)) / 1000.0
    return overlays.scale_bar(rgb, km_per_px, note)


def save_pair(out, stem, elevation, sea, north, pixel_deg, exaggeration, west=None, centre_latitude=0.0):
    save_rgb(out / f'{stem}-elevation.png',
             annotate(to_uint8(terrain_base(elevation, sea)), north, pixel_deg, west, centre_latitude))
    save_rgb(out / f'{stem}-relief.png',
             annotate(to_uint8(relief_rgb(elevation, sea, north, pixel_deg, exaggeration)), north, pixel_deg, west,
                      centre_latitude))
    return f'{stem}-elevation.png', f'{stem}-relief.png'


def extent_text(extent):
    def lon(value):
        return f"{abs(value):.2f}°{'E' if value >= 0 else 'W'}"

    def lat(value):
        return f"{abs(value):.2f}°{'N' if value >= 0 else 'S'}"
    return f"{lon(extent['west'])} to {lon(extent['east'])}, {lat(extent['south'])} to {lat(extent['north'])}"


def centre_text(longitude, latitude):
    return f"{abs(latitude):g}°{'N' if latitude >= 0 else 'S'} {abs(longitude):g}°{'E' if longitude >= 0 else 'W'}"


# ---------------------------------------------------------------- the whole Earth

def globe(out, ocean, entries):
    sums = {columns: [] for columns in GLOBE_COLUMNS}
    bands = int(round(180.0 / BAND_DEG))
    for band in range(bands):
        north = 90.0 - band * BAND_DEG
        array, extent = etopo.read_box(-180.0, north - BAND_DEG, 180.0, north, 'etopo60')
        for columns in GLOBE_COLUMNS:
            pixel = 360.0 / columns
            band_rows = int(round(BAND_DEG / pixel))
            longitudes = -180.0 + np.arange(columns + 1) * pixel
            latitudes = north - np.arange(band_rows + 1) * pixel
            sums[columns].append(ew.box_average(array, extent, longitudes, latitudes).astype(np.float32))
        del array
        print(f'globe band {band + 1}/{bands}', flush=True)
    for columns in GLOBE_COLUMNS:
        elevation = np.concatenate(sums[columns], axis=0)
        pixel = 360.0 / columns
        record = {'centre': [-180.0, 90.0], 'pixel_deg': pixel, 'columns_index': [0, columns],
                  'rows_index': [0, columns // 2]}
        sea = ew.window_sea(elevation, record, ocean)
        files = save_pair(out, f'globe-{columns:04d}', elevation, sea, 90.0, pixel, GLOBE_EXAGGERATION)
        km = RADIUS_M * math.radians(pixel) / 1000.0
        caption = (f'The whole Earth, {columns} x {columns // 2} pixels; {pixel:g}° per pixel '
                   f'({km:,.1f} km north-south); box-averaged from the {SOURCE_TEXT["etopo60"]}')
        add_entries(entries, 'globe', files, caption, f'globe-{columns}', [0.0, 0.0],
                    {'west': -180.0, 'south': -90.0, 'east': 180.0, 'north': 90.0}, pixel, columns, columns // 2,
                    GLOBE_EXAGGERATION, [etopo.GLOBE_60])
        print(f'globe {columns} written', flush=True)


# ---------------------------------------------------------------- windows

def window(out, ocean, entries, group, name, title, longitude, latitude, pixel_deg, pixels, exaggeration, source):
    entry = {'name': name, 'centre': [longitude, latitude], 'box': [-180.0, -90.0, 180.0, 90.0]}
    request = {'columns': pixels[0], 'rows': pixels[1], 'metres_per_pixel': RADIUS_M * math.radians(pixel_deg)}
    grid = ew.window_grid(entry, request, clip=False)
    extent = grid['extent']
    sources = ew.etopo_files(source, extent)
    missing = [relative(path) for path in sources if not Path(path).is_file()]
    if missing:
        print(f'{name}: skipped, missing {", ".join(missing)}', flush=True)
        return f'{title}: skipped (missing {", ".join(missing)})'
    array, source_extent = etopo.read_box(extent['west'], extent['south'], extent['east'], extent['north'], source)
    elevation = ew.box_average(array, source_extent, *ew.edges(grid)).astype(np.float32)
    sea = ew.window_sea(elevation, grid, ocean)
    if source == 'etopo15':  # pixels finer than the mask: below-0 m pixels joined to the ocean are sea
        labels, _ = ndimage.label(elevation < 0.0)
        sea = np.isin(labels, np.unique(labels[sea & (labels > 0)])) & (labels > 0)
    files = save_pair(out, f'{group}-{name}', elevation, sea, extent['north'], grid['pixel_deg'], exaggeration,
                      west=extent['west'], centre_latitude=latitude)
    km = request['metres_per_pixel'] / 1000.0
    caption = (f'{title}; centre {centre_text(longitude, latitude)}; {extent_text(extent)} '
               f'({pixels[0] * pixel_deg:g}° x {pixels[1] * pixel_deg:g}°, {pixels[0]} x {pixels[1]} pixels); '
               f'{pixel_deg:g}° per pixel ({km:,.2f} km north-south); box-averaged from the {SOURCE_TEXT[source]}')
    add_entries(entries, group, files, caption, name, [longitude, latitude], extent, grid['pixel_deg'],
                grid['pixels'][0], grid['pixels'][1], exaggeration, sources)
    print(f'{group} {name} written', flush=True)
    return None


def add_entries(entries, group, files, caption, name, centre, extent, pixel_deg, columns, rows, exaggeration, sources):
    for file, style in zip(files, ('elevation', 'relief')):
        entries.append({'file': file, 'group': group, 'kind': f'{group}-{style}', 'name': name, 'caption': caption,
                        'centre': centre, 'extent': extent, 'degrees_per_pixel': pixel_deg,
                        'km_per_pixel_north_south': RADIUS_M * math.radians(pixel_deg) / 1000.0,
                        'columns': int(columns), 'rows': int(rows),
                        'vertical_exaggeration': exaggeration if style == 'relief' else None,
                        'source_files': [relative(path) for path in sources]})


# ---------------------------------------------------------------- page

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Earth atlas</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}h2{{margin:24px 0 8px}}
.row{{display:flex;gap:8px;flex-wrap:wrap}}figure{{margin:0 0 12px;flex:1 1 100%;min-width:280px}}
img{{width:100%;image-rendering:pixelated;cursor:zoom-in}}figcaption{{color:#999;margin:4px 0}}
pre{{background:#0d0f12;padding:8px;overflow:auto;font-size:12px}}</style></head><body>
<h1>Earth atlas for the terrain sandbox</h1><p>{note}</p>{sections}</body></html>"""

SECTIONS = [('globe', 'The whole Earth (128 to 2048 columns; relief exaggeration 8)'),
            ('window', 'Regional windows at the sandbox window zoom (k 64 of a 256-column world, about 2.44 km; exaggeration 3)'),
            ('close', 'Closer windows (k 256, about 610 m; exaggeration 2)')]


def write_page(out, entries, skipped, created):
    html = []
    for group, heading in SECTIONS:
        figures = ''.join(f'<figure><a href="{entry["file"]}"><img src="{entry["file"]}" loading="lazy"></a>'
                          f'<figcaption>{entry["kind"].split("-")[-1]}: {entry["caption"]}</figcaption></figure>'
                          for entry in entries if entry['group'] == group)
        if figures:
            html.append(f'<h2>{heading}</h2><div class="row">{figures}</div>')
    if skipped:
        html.append('<h2>Skipped</h2><pre>' + '\n'.join(skipped) + '</pre>')
    note = ('ETOPO 2022 surface elevation (ice surface, lake beds), area-weighted box means in degrees. Colours: '
            'tectonics:elevation (terrain_base); relief: times the hillshade (light from 315°, 45° high) as the '
            "sandbox's maps. Sea is the world ocean (below 0 m and connected to it on the 60 arc-second globe; on "
            'the 15 arc-second windows also the below-0 m pixels joined to it); basins below 0 m off it are land. '
            f'No rivers, lakes or ice drawn. Written {created} by scripts/sandbox/earth_atlas.py.')
    (out / 'index.html').write_text(PAGE.format(note=note, sections=''.join(html)), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--out', default=str(ROOT / 'refs' / 'processed' / 'earth-atlas'))
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    ocean = ew.global_ocean()
    entries, skipped = [], []
    globe(out, ocean, entries)
    for name, title, longitude, latitude in WINDOWS:
        window(out, ocean, entries, 'window', name, title, longitude, latitude, WINDOW_DEG, WINDOW_PIXELS,
               WINDOW_EXAGGERATION, 'etopo60')
    for name, title, longitude, latitude in CLOSE:
        message = window(out, ocean, entries, 'close', name, title, longitude, latitude, CLOSE_DEG, CLOSE_PIXELS,
                         CLOSE_EXAGGERATION, 'etopo15')
        if message:
            skipped.append(message)
    created = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    seconds = time.perf_counter() - started
    manifest = {'schema': 'uw-earth-atlas/v1', 'created': created, 'seconds': round(seconds, 1),
                'code': {'module': 'scripts/sandbox/earth_atlas.py', 'revision': ew.code_revision()},
                'earth_radius_m': RADIUS_M, 'ocean_mask': 'refs/processed/earth-windows/ocean-60s/ocean.npy',
                'skipped': skipped, 'images': entries}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
    write_page(out, entries, skipped, created)
    print(f'{len(entries)} images in {relative(out)} in {seconds:.1f} s', flush=True)


if __name__ == '__main__':
    main()
