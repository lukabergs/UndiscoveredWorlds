"""Terrain metrics of terrain sandbox runs against Earth: one table per run label and domain, measured the
same way every time (the definitions of the T1 diagnoses of docs/architecture/terrain-roadmap.md: v17, lakes
and pits, v18, v19, v20).

  uv run python scripts/sandbox/terrain_metrics.py --name NAME --labels v18,v19
      [--domains global,B,F,G] [--earth-windows norway,patagonia,alps | none]

A label (an id, <id>-<tag>, an earlier label or an absolute path: naming.find) reads the main tree's
out/scratch/sandbox/<name>/: config.json (the domain geometry), <domain>-elevation.npy (the
final surface), -elevation-physics.npy (the physics' surface before the lakes and rivers are drawn; the final
surface when missing), -elevation-procedural.npy (the physics input; without it the drowned land is not
measured and the sea is not seeded), -lake-surface.npy (NaN off the lakes), the totals -glacial_eroded_m,
-ice_time_myr, -eroded_m and -province.npy when present. A missing domain is skipped. Writes summary.json and
index.html to the main tree's out/evaluation/<id>-metrics-<NAME>/ (scripts/uwtools/naming.py, kind page: the tag
metrics-<NAME> follows the tag rule); earlier pages of that tag are then marked superseded (naming.finish).

Geometry: the base world is 256 x 128 cells; the global domain has k = columns / 256 cells per world cell, a
window the config's inclusive world pixels [c0, r0, c1, r1] at its `scale` k. Cells are dy x dy cos(latitude)
km with dy = 6371 km x (180 / 128 / k) degrees in radians; every share and median is weighted by that area.

Masks:
- sea: cells at or below 0 m, 8-connected components; seeded by the physics input's sea (the components of
  the physics surface that touch it); without the physics input (and for the input itself) the components
  of at least 2000 cells or touching the domain edge. Sea is taken from the physics surface.
- land: not sea and final elevation above 0 m. Lakes are not sea; "not sea" is the denominator of the lake
  share and of the coastline.

Metrics (per label and domain; Earth where defined):
- land %: land area / domain area; land median: area-weighted median of the final elevation on land
  (interpolated quantile).
- land bands: area shares of the whole domain surface (not of the land) with elevation in [0, 200),
  [200, 500), [500, 1000), [1000, 2000), [2000, 3000), [3000, 9000) m (spherical band areas on the global
  grid, the cell areas in a window); no sea mask.
- lakes: lake cells are finite lake-surface cells; share = lake area / not-sea area; count = 8-connected
  lake components; the largest by area with its deepest point (surface minus final elevation).
- drowned land: cells sea after the physics that are not sea in the physics input, % of the input's land
  (its not-sea area); glacially cut share: of the drowned area, glacial_eroded > 50 m or ice_time > 0;
  median elevation of the drowned cells on the physics surface (unweighted, metres, negative).
- coastline: sea / not-sea 4-neighbour edges in km (east-west neighbours share an edge of dy km, north-south
  ones dy cos(latitude of the edge); no wrap across the date line) per 1000 km^2 of not-sea area. Earth
  3.7 on the 2048 grid.
- plains: local relief = max - min over 3 x 3 cells of the elevation with sea set to 0 m (columns wrap on
  the global domain, rows and window edges repeat the nearest cell); plain = land with relief under 100 m;
  plain share of land; the area-weighted median plain patch (8-connected; the patch size at which half the
  plain area lies in patches up to that size) in M km^2; the land's relief median; per elevation band
  (lo, hi] of the land: land share, relief median, plain share.
- ice and erosion: iced land = input land (not-sea of the physics input, else land) with ice_time > 0, %
  of the input land; mean glacial_eroded on it; mean eroded_m on the input land.
- provinces (when -province.npy exists): per code (0 oceanic, 1 shield, 2 platform, 3 sedimentary basin,
  4 orogen, 5 rift / active margin) the share of the land, its plain share and relief median.

Earth: the global domain against refs/processed/heightmap/2048.tif (ETOPO 2022; the path of
uwtools.tectonics_earth_reference; land is above 0 m, sea the large or edge-touching components at or below
0 m; no lakes, drowned land or ice). The cached hypsography reference (GSHHG mask, which counts lakes and land
below the sea) is recorded beside it. Windows against the atlas windows of scripts/sandbox/earth_atlas.py
(1280 x 1024 pixels at 180 / 128 / 64 degrees, box-averaged from the ETOPO 60 arc-second globe through
scripts/evaluation/earth_windows.py; sea is the world ocean, below-0 m basins off it are land).
"""
import argparse
import datetime
import html
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'evaluation'))
from uwtools import grid, naming  # noqa: E402
from uwtools.tectonics_earth_reference import load_earth_reference, reference_inputs  # noqa: E402

WORLD_COLUMNS, WORLD_ROWS, RADIUS_KM = 256, 128, 6371.0
BANDS = (0, 200, 500, 1000, 2000, 3000, 9000)
ELEVATION_BANDS = tuple(zip(BANDS[:-1], BANDS[1:]))
PLAIN_M = 100.0
GLACIAL_M = 50.0
SEA_MIN_CELLS = 2000
EIGHT = np.ones((3, 3))
PROVINCES = ('oceanic', 'shield', 'platform', 'sedimentary basin', 'orogen', 'rift / active margin')
DEFAULT_EARTH_WINDOWS = 'norway,patagonia,alps'


def sea_from(elevation, seed=None, min_cells=SEA_MIN_CELLS):
    """Cells <= 0 m (8-connected) in a component touching seed, else in a component that is large or touches
    the domain edge."""
    labels, count = ndimage.label(elevation <= 0.0, structure=EIGHT)
    if seed is not None:
        keep = np.unique(labels[(labels > 0) & seed])
        return np.isin(labels, keep[keep > 0])
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    edge = np.unique(np.concatenate([labels[0], labels[-1], labels[:, 0], labels[:, -1]]))
    keep = set(np.nonzero(sizes >= min_cells)[0].tolist()) | set(edge.tolist())
    keep.discard(0)
    return np.isin(labels, list(keep))


def geometry(config, domain, shape):
    """(cell area km^2, centre latitudes, dy km) of a domain's saved array."""
    rows, columns = shape
    if domain == 'global':
        k = columns // WORLD_COLUMNS
        north = 90.0
    else:
        c0, r0, c1, r1 = config['domains'][domain]['pixels']
        k = config['domains'][domain]['scale']
        if (((c1 - c0) % WORLD_COLUMNS + 1) * k, (r1 - r0 + 1) * k) != (columns, rows):
            raise SystemExit(f'{domain}: array {shape} does not match the config pixels {[c0, r0, c1, r1]} at k {k}')
        north = 90.0 - r0 * 180.0 / WORLD_ROWS
    pixel_deg = 180.0 / WORLD_ROWS / k
    lat = north - (np.arange(rows) + 0.5) * pixel_deg
    dy_km = RADIUS_KM * math.radians(pixel_deg)
    area = dy_km * dy_km * np.cos(np.radians(lat))[:, None] * np.ones((1, columns))
    return area, lat, dy_km


def wquant(values, weights, qs):
    order = np.argsort(values)
    v, w = values[order], weights[order]
    c = np.cumsum(w) / w.sum()
    return [float(np.interp(q, c, v)) for q in qs]


def local_relief(z, size, wrap):
    mode = ('nearest', 'wrap') if wrap else 'nearest'
    return ndimage.maximum_filter(z, size=size, mode=mode) - ndimage.minimum_filter(z, size=size, mode=mode)


def coastline_km(sea, lat, dy_km):
    east_west = float((sea[:, :-1] != sea[:, 1:]).sum()) * dy_km
    lat_edge = 0.5 * (lat[:-1] + lat[1:])
    north_south = float(((sea[:-1] != sea[1:]).sum(axis=1) * dy_km * np.cos(np.radians(lat_edge))).sum())
    return east_west + north_south


def share(numerator, denominator):
    return float(100.0 * numerator / denominator) if denominator > 0 else None


def land_metrics(z, land, sea, area, lat, dy_km, wrap, band_weights):
    """The metrics defined for worlds and Earth alike."""
    out = {'land_pct': share(area[land].sum(), area.sum()),
           'land_median_m': wquant(z[land], area[land], (0.5,))[0] if land.any() else None}
    histogram, _ = np.histogram(z, bins=BANDS, weights=band_weights)
    out['band_pct_of_surface'] = {f'{lo}-{hi}': float(100 * v / band_weights.sum())
                                  for (lo, hi), v in zip(ELEVATION_BANDS, histogram)}
    not_sea = float(area[~sea].sum())
    out['coast_km_per_1000km2'] = 1000.0 * coastline_km(sea, lat, dy_km) / not_sea if not_sea > 0 else None
    relief = local_relief(np.where(land, z, 0.0), 3, wrap)
    plain = land & (relief < PLAIN_M)
    labels, count = ndimage.label(plain, structure=EIGHT)
    out['plain_pct_of_land'] = share(area[plain].sum(), area[land].sum())
    out['plain_patch_median_Mkm2'] = None
    if count:
        sizes = np.sort(ndimage.sum(area, labels, index=np.arange(1, count + 1)))
        cumulative = np.cumsum(sizes) / sizes.sum()
        out['plain_patch_median_Mkm2'] = float(sizes[np.searchsorted(cumulative, 0.5)]) / 1e6
    out['relief_median_m'] = wquant(relief[land], area[land], (0.5,))[0] if land.any() else None
    by_band = {}
    for lo, hi in ELEVATION_BANDS:
        m = land & (z > lo) & (z <= hi)
        if area[m].sum() == 0:
            continue
        by_band[f'{lo}-{hi}'] = {'land_pct': share(area[m].sum(), area[land].sum()),
                                 'relief_median_m': wquant(relief[m], area[m], (0.5,))[0],
                                 'plain_pct': share(area[m & plain].sum(), area[m].sum())}
    out['relief_by_band'] = by_band
    return out, relief


def optional(directory, domain, field):
    path = directory / f'{domain}-{field}.npy'
    return np.load(path).astype(np.float64) if path.is_file() else None


def run_directory(label):
    """The main tree's sandbox run a label names (naming.find: an id, <id>-<tag>, an earlier label or an absolute
    path); SystemExit when there is none."""
    try:
        return naming.find('sandbox', label)
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(str(error))


def run_metrics(label, domain):
    directory = run_directory(label)
    config = json.loads((directory / 'config.json').read_text(encoding='utf-8'))
    final = np.load(directory / f'{domain}-elevation.npy').astype(np.float64)
    notes = []
    physics = optional(directory, domain, 'elevation-physics')
    if physics is None:
        physics = final
        notes.append('no elevation-physics: sea from the final surface')
    pre = optional(directory, domain, 'elevation-procedural')
    area, lat, dy_km = geometry(config, domain, final.shape)
    wrap = domain == 'global'
    if pre is not None:
        pre_sea = sea_from(pre)
        sea = sea_from(physics, seed=pre_sea)
        pre_land = ~pre_sea
    else:
        pre_sea = None
        sea = sea_from(physics)
        pre_land = ~sea
        notes.append('no elevation-procedural (the physics input): drowned land not measured, sea not seeded, '
                     'input land = final not-sea')
    land = ~sea & (final > 0)
    band_weights = grid.area_weights(*final.shape) if wrap else area
    out, relief = land_metrics(final, land, sea, area, lat, dy_km, wrap, band_weights)
    out['shape'] = list(final.shape)

    not_sea = float(area[~sea].sum())
    surface = optional(directory, domain, 'lake-surface')
    out['lakes'] = None
    if surface is not None:
        lake = np.isfinite(surface)
        labels, count = ndimage.label(lake, structure=EIGHT)
        lakes = {'pct_of_land': share(area[lake].sum(), not_sea), 'count': int(count),
                 'largest_km2': None, 'largest_max_depth_m': None}
        if count:
            index = np.arange(1, count + 1)
            depth = np.where(lake, surface - final, 0.0)
            areas = ndimage.sum(area, labels, index)
            deepest = ndimage.maximum(depth, labels, index)
            j = int(np.argmax(areas))
            lakes['largest_km2'] = float(areas[j])
            lakes['largest_max_depth_m'] = float(deepest[j])
        out['lakes'] = lakes
    else:
        notes.append('no lake-surface')

    glacial = optional(directory, domain, 'glacial_eroded_m')
    ice = optional(directory, domain, 'ice_time_myr')
    out['drowned'] = None
    if pre_sea is not None:
        drowned = sea & pre_land
        w = area[drowned]
        entry = {'pct_of_input_land': share(w.sum(), area[pre_land].sum()), 'glacial_pct': None,
                 'median_elevation_m': float(np.median(physics[drowned])) if drowned.any() else None}
        if glacial is not None and ice is not None:
            cut = (glacial > GLACIAL_M) | (ice > 0)
            entry['glacial_pct'] = float(100 * area[drowned & cut].sum() / max(w.sum(), 1e-9))
        out['drowned'] = entry

    out['ice'] = None
    if glacial is not None and ice is not None:
        iced = pre_land & (ice > 0)
        out['ice'] = {'iced_pct_of_input_land': share(area[iced].sum(), area[pre_land].sum()),
                      'glacial_mean_m_iced': float((glacial * area)[iced].sum() / max(area[iced].sum(), 1e-9))}
    eroded = optional(directory, domain, 'eroded_m')
    out['eroded_mean_m_input_land'] = (float((eroded * area)[pre_land].sum() / max(area[pre_land].sum(), 1e-9))
                                       if eroded is not None else None)

    out['provinces'] = None
    path = directory / f'{domain}-province.npy'
    if path.is_file():
        province = np.load(path)
        per = {}
        for code, name in enumerate(PROVINCES):
            mask = land & (province == code)
            if area[mask].sum() <= 0:
                continue
            per[name] = {'land_pct': share(area[mask].sum(), area[land].sum()),
                         'plain_pct': share(area[mask & (relief < PLAIN_M)].sum(), area[mask].sum()),
                         'relief_median_m': wquant(relief[mask], area[mask], (0.5,))[0]}
        out['provinces'] = per
    out['notes'] = notes
    out['inputs'] = {'seed': config.get('seed'), 'base_bundle': config.get('base_bundle'),
                     'domain': config['domains'].get(domain)}
    return out


def earth_global():
    from PIL import Image
    path = reference_inputs(ROOT, 2048)[0]
    with Image.open(path) as image:
        z = np.array(image).astype(np.float64)
    rows, columns = z.shape
    dy_km = RADIUS_KM * math.radians(180.0 / rows)
    lat = grid.latitudes(rows)
    area = dy_km * dy_km * np.cos(np.radians(lat))[:, None] * np.ones((1, columns))
    sea = sea_from(z)
    out, _ = land_metrics(z, ~sea & (z > 0), sea, area, lat, dy_km, True, grid.area_weights(rows, columns))
    out['source'] = path.relative_to(ROOT).as_posix()
    out['shape'] = [rows, columns]
    return out


def earth_window(name):
    import earth_atlas as atlas
    import earth_windows as ew
    from uwtools import etopo
    matches = [item for item in atlas.WINDOWS if item[0] == name]
    if not matches:
        raise SystemExit(f'unknown Earth window {name}; earth_atlas.WINDOWS: {[item[0] for item in atlas.WINDOWS]}')
    _, title, longitude, latitude = matches[0]
    entry = {'name': name, 'centre': [longitude, latitude], 'box': [-180.0, -90.0, 180.0, 90.0]}
    request = {'columns': atlas.WINDOW_PIXELS[0], 'rows': atlas.WINDOW_PIXELS[1],
               'metres_per_pixel': atlas.RADIUS_M * math.radians(atlas.WINDOW_DEG)}
    window = ew.window_grid(entry, request, clip=False)
    extent = window['extent']
    array, source_extent = etopo.read_box(extent['west'], extent['south'], extent['east'], extent['north'], 'etopo60')
    z = ew.box_average(array, source_extent, *ew.edges(window)).astype(np.float64)
    sea = ew.window_sea(z, window, ew.global_ocean())
    lat = extent['north'] - (np.arange(z.shape[0]) + 0.5) * window['pixel_deg']
    dy_km = RADIUS_KM * math.radians(window['pixel_deg'])
    area = dy_km * dy_km * np.cos(np.radians(lat))[:, None] * np.ones((1, z.shape[1]))
    out, _ = land_metrics(z, ~sea & (z > 0), sea, area, lat, dy_km, False, area)
    out.update({'title': title, 'centre': [longitude, latitude], 'extent': extent, 'shape': list(z.shape),
                'source': 'etopo60'})
    return out


def code_revision():
    try:
        head = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', 'scripts/sandbox/terrain_metrics.py'], cwd=ROOT,
                               capture_output=True, text=True, check=True).stdout.strip()
        return head + ('+dirty' if dirty else '')
    except (OSError, subprocess.CalledProcessError):
        return None


PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>Terrain metrics {name}</title>
<style>body{{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}}
table{{border-collapse:collapse;margin:12px 0}}td,th{{border-bottom:1px solid #333;padding:4px 10px;text-align:right}}
th{{vertical-align:bottom}}td:first-child,th:first-child,td:nth-child(2),th:nth-child(2){{text-align:left}}
tr.earth td{{color:#9fc4f0}}p.note{{color:#aaa;max-width:70em}}</style></head><body><h1>Terrain metrics: {name}</h1>
<p class="note">Definitions: the docstring of scripts/sandbox/terrain_metrics.py. Shares weighted by spherical
cell area; land is not-sea above 0 m; sea is taken from the physics surface, seeded by the physics input's sea.
Earth: {earth}. Code {revision}; {created}; {seconds:.0f} s.</p>{tables}{notes}</body></html>"""


def number(value, digits=1):
    if value is None:
        return ''
    if digits == 0:
        return f'{value:,.0f}'
    return f'{value:.{digits}f}'


def get(record, *keys):
    for key in keys:
        if record is None:
            return None
        record = record.get(key) if isinstance(record, dict) else None
    return record


def table(title, headers, rows):
    head = ''.join(f'<th>{html.escape(h)}</th>' for h in headers)
    body = ''.join(f'<tr class="{css}">' + ''.join(f'<td>{html.escape(str(c))}</td>' for c in cells) + '</tr>'
                   for css, cells in rows)
    return f'<h2>{html.escape(title)}</h2><table><tr>{head}</tr>{body}</table>'


def page_rows(summary):
    """(css class, world, domain, record) in display order: per domain the labels, Earth after the global
    domain; the Earth windows last."""
    rows = []
    for domain in summary['domains']:
        for label in summary['labels']:
            record = summary['runs'][label].get(domain)
            if record is not None:
                rows.append(('', label, domain, record))
        if domain == 'global' and summary['earth'].get('global'):
            rows.append(('earth', 'Earth (ETOPO 2048)', 'global', summary['earth']['global']))
    for name, record in summary['earth'].get('windows', {}).items():
        rows.append(('earth', f'Earth: {name}', 'window', record))
    return rows


def render(summary):
    rows = page_rows(summary)
    bands = [f'{lo}-{hi}' for lo, hi in ELEVATION_BANDS]
    band_titles = ['0-200', '200-500', '500-1000', '1000-2000', '2000-3000', '>3000']
    tables = [
        table('Land and height (bands: % of the domain surface)',
              ['world', 'domain', 'land %', 'land median m'] + band_titles,
              [(css, [world, domain, number(r['land_pct']), number(r['land_median_m'], 0)]
                + [number(r['band_pct_of_surface'][b]) for b in bands]) for css, world, domain, r in rows]),
        table('Water and coast',
              ['world', 'domain', 'lakes % of land', 'lakes', 'largest lake km2', 'its max depth m',
               'drowned % of input land', 'drowned glacially cut %', 'drowned median elevation m',
               'coast km per 1000 km2'],
              [(css, [world, domain, number(get(r, 'lakes', 'pct_of_land'), 2), number(get(r, 'lakes', 'count'), 0),
                      number(get(r, 'lakes', 'largest_km2'), 0), number(get(r, 'lakes', 'largest_max_depth_m'), 0),
                      number(get(r, 'drowned', 'pct_of_input_land'), 2), number(get(r, 'drowned', 'glacial_pct')),
                      number(get(r, 'drowned', 'median_elevation_m'), 0), number(r['coast_km_per_1000km2'], 2)])
               for css, world, domain, r in rows]),
        table('Roughness (relief over 3 x 3 cells; plain: under 100 m; relief median m per elevation band of the land)',
              ['world', 'domain', 'plain % of land', 'plain patch median M km2', 'relief median m'] + band_titles,
              [(css, [world, domain, number(r['plain_pct_of_land']), number(r['plain_patch_median_Mkm2'], 2),
                      number(r['relief_median_m'], 0)]
                + [number(get(r, 'relief_by_band', b, 'relief_median_m'), 0) for b in bands])
               for css, world, domain, r in rows]),
        table('Ice, erosion, provinces (province: % of the land)',
              ['world', 'domain', 'iced % of input land', 'glacial erosion on iced land m', 'eroded on input land m']
              + list(PROVINCES),
              [(css, [world, domain, number(get(r, 'ice', 'iced_pct_of_input_land')),
                      number(get(r, 'ice', 'glacial_mean_m_iced'), 0), number(r.get('eroded_mean_m_input_land'), 0)]
                + [number(get(r, 'provinces', p, 'land_pct')) for p in PROVINCES])
               for css, world, domain, r in rows]),
    ]
    notes = [f'<li>{html.escape(world)} {html.escape(domain)}: {html.escape(note)}</li>'
             for _, world, domain, r in rows for note in r.get('notes', [])]
    reference = summary['earth'].get('reference')
    earth = 'refs/processed/heightmap/2048.tif with land above 0 m'
    if reference:
        earth += (f"; the cached hypsography reference (GSHHG mask) has land {100 * reference['land_fraction']:.1f} % "
                  f"and a land median of {reference['land_median_m']:.0f} m")
    return PAGE.format(name=html.escape(summary['name']), earth=html.escape(earth),
                       revision=html.escape(str(summary['code_revision'])), created=summary['created'],
                       seconds=summary['seconds'], tables=''.join(tables),
                       notes=f'<h2>Notes</h2><ul>{"".join(notes)}</ul>' if notes else '')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', required=True)
    parser.add_argument('--labels', required=True)
    parser.add_argument('--domains', default='global,B,F,G')
    parser.add_argument('--earth-windows', default=DEFAULT_EARTH_WINDOWS,
                        help="earth_atlas.WINDOWS names, comma-separated; 'none' for no Earth windows")
    arguments = parser.parse_args()
    try:
        tag = naming.check_tag(f'metrics-{arguments.name}')
    except ValueError as error:
        parser.error(str(error))
    start = time.perf_counter()
    labels = [item for item in arguments.labels.split(',') if item]
    domains = [item for item in arguments.domains.split(',') if item]
    windows = [] if arguments.earth_windows in ('', 'none') else arguments.earth_windows.split(',')

    summary = {'name': arguments.name, 'labels': labels, 'domains': domains, 'runs': {}, 'earth': {},
               'definitions': 'scripts/sandbox/terrain_metrics.py docstring', 'code_revision': code_revision(),
               'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    for label in labels:
        if not (run_directory(label) / 'config.json').is_file():
            raise SystemExit(f'no run {run_directory(label)}')
        summary['runs'][label] = {}
        for domain in domains:
            if not (run_directory(label) / f'{domain}-elevation.npy').is_file():
                continue
            summary['runs'][label][domain] = run_metrics(label, domain)
            print(f'{label} {domain}', flush=True)
    if 'global' in domains:
        summary['earth']['global'] = earth_global()
        metrics = load_earth_reference(ROOT, 2048)['metrics']
        summary['earth']['reference'] = {'land_fraction': metrics['land_fraction'],
                                         'land_median_m': metrics['land_elevation_m']['median']}
        print('earth global', flush=True)
    if windows:
        summary['earth']['windows'] = {}
        for name in windows:
            summary['earth']['windows'][name] = earth_window(name)
            print(f'earth {name}', flush=True)
    summary['seconds'] = time.perf_counter() - start

    named = naming.new_directory('page', tag)
    directory = named.path
    (directory / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n', encoding='utf-8')
    (directory / 'index.html').write_text(render(summary), encoding='utf-8')
    print(directory / 'index.html')
    print(f'{summary["seconds"]:.1f} s')
    naming.finish(named, True)  # earlier pages of the tag marked superseded


if __name__ == '__main__':
    main()
