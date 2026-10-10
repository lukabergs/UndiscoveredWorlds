"""Earth reference windows beside regional maps (multiscale plan S1.2; section 3 "Evaluation").

A regional map (a bundle with a window block) is shown beside a window of the Earth from the
catalogue configs/evaluation/earth-windows.json (schema uw-earth-windows/v1): the entry named by
--earth NAME, otherwise the entry whose metres_per_pixel range holds the map's and whose centre
latitude is nearest the window's centre latitude (none: no Earth column, and the page says so).

The window has the map's pixel count and north-south metres per pixel, its pixels square in
degrees as the regional grid's cells are, centred on the entry's centre on the Earth's radius
(6,371 km). From the 60-arc-second globe it has the map's full pixel count (cut only at the poles);
from the 15-arc-second tiles it is cut to the entry's box, the area whose tiles are downloaded
(multiscale plan S1.2, the orchestrator's decision). An entry named by --earth whose
metres_per_pixel range does not hold the map's is drawn and its caption says so. Its grid is anchored at the
centre (a pixel corner), so windows of one entry and pixel size nest. The elevation is read with
uwtools.etopo.read_box from ETOPO 2022 at 60 arc-seconds when the pixel is at least 1,850 m, at
15 arc-seconds otherwise, and box-averaged (area-weighted) to the pixels. Processed windows are
cached under refs/processed/earth-windows/<name>/<metres>m/ (elevation.npy and manifest.json:
source, extents, pixel size, catalogue entry, code version) and reused when the cached window
covers the request on the same grid, with the same CACHE_VERSION and source files (sizes and modification
times). Sea is the world ocean: a window pixel below 0 m whose 60-arc-second globe cell is in the
ocean mask (global_ocean: the globe's cells below 0 m 4-connected, across the date line, to a seed in
the open Pacific at 0 N 160 W; cached as refs/processed/earth-windows/ocean-60s/ocean.npy under the
same staleness rule). Basins below 0 m off the ocean (the Caspian, the Dead Sea, Death Valley, the
Qattara depression, Lake Superior's bed) are drawn as land; there are no lakes, rivers or ice.
Each window is drawn in the style of the profile's regional maps it stands beside: `terrain`
(regional_terrain's sea and land on tectonics:elevation and the landmass outline) and `relief`
(the same times the hillshade of maps.regional_relief, the map's vertical exaggeration).
Map entries have kind "earth" (never "reference", the climate validation's).
"""
import datetime
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import map_rendering, palettes
from maps import RELIEF_EXAGGERATION, hillshade, legend, pixel_sizes_m, shade_factor, shore, terrain_base

ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = 'configs/evaluation/earth-windows.json'
CATALOGUE_SCHEMA = 'uw-earth-windows/v1'
CACHE = ROOT / 'refs/processed/earth-windows'
WINDOW_SCHEMA = 'uw-earth-window/v1'
# Bump when the processing of a cached window (the grid, the box average) changes: older caches are recomputed.
CACHE_VERSION = 1
EARTH_RADIUS_M = 6_371_000.0
# The 60-arc-second globe is used from this pixel size up (its own pixel is about 1,850 m north-south).
COARSE_SOURCE_FROM_M = 1850.0
SOURCE_LABELS = {'etopo60': '60 arc-second', 'etopo15': '15 arc-second'}
# Composite of a profile map -> style of the Earth window beside it.
STYLES = {'regional_terrain': 'terrain', 'regional_relief': 'relief'}
NEAREST = "the entry whose metres_per_pixel range holds the map's and whose centre latitude is nearest the window's"
EPSILON = 1e-9


def map_request(manifest):
    """The regional map's grid from its bundle manifest: columns, rows, pixel size in degrees (north-south),
    north-south metres per pixel on the planet's radius (parameters world.body.radius_km, 6371 km when
    absent), the window's centre latitude and extent. A resampled bundle (window pixel_columns, multiscale
    plan S2.4) has pixels of metres_per_pixel from the kept box's north-west corner with the partial edge
    pixels dropped: its extent and centre are the pixels', not the box's."""
    window, grid = manifest['window'], manifest['grid']
    columns, rows = int(grid['columns']), int(grid['rows'])
    north, south = float(window['north']), float(window['south'])
    radius_m = 1000.0 * float(manifest.get('parameters', {}).get('world.body.radius_km', EARTH_RADIUS_M / 1000.0))
    extent = {key: float(window[key]) for key in ('west', 'south', 'east', 'north')}
    if 'pixel_columns' in window:
        pixel_deg = math.degrees(float(window['metres_per_pixel']) / radius_m)
        south = north - rows * pixel_deg
        east = extent['west'] + columns * pixel_deg
        extent.update(south=south, east=east - 360.0 if east > 180.0 else east)  # west > east across the date line
    else:
        pixel_deg = (north - south) / rows
    return {'columns': columns, 'rows': rows, 'pixel_deg': pixel_deg, 'radius_m': radius_m,
            'metres_per_pixel': radius_m * math.radians(pixel_deg), 'latitude': 0.5 * (north + south),
            'extent': extent}


def source_for(metres_per_pixel):
    """'etopo60' when a pixel is at least 1,850 m, else 'etopo15'."""
    return 'etopo60' if metres_per_pixel >= COARSE_SOURCE_FROM_M else 'etopo15'


def load_catalogue(path=None):
    path = ROOT / CATALOGUE if path is None else Path(path)
    catalogue = json.loads(path.read_text(encoding='utf-8'))
    if catalogue.get('schema') != CATALOGUE_SCHEMA:
        raise ValueError(f'{path}: expected schema {CATALOGUE_SCHEMA}')
    return catalogue


def choose_entry(catalogue, metres_per_pixel, latitude, name=None):
    """The catalogue entry named `name` (ValueError when there is none), otherwise the entry whose
    metres_per_pixel [min, max] holds the map's and whose centre latitude is nearest `latitude` (the
    first in catalogue order on a tie), or None."""
    windows = catalogue['windows']
    if name is not None:
        for entry in windows:
            if entry['name'] == name:
                return entry
        raise ValueError(f"No Earth window {name!r}; known: {', '.join(entry['name'] for entry in windows)}")
    candidates = [entry for entry in windows
                  if float(entry['metres_per_pixel'][0]) <= metres_per_pixel <= float(entry['metres_per_pixel'][1])]
    return min(candidates, key=lambda entry: abs(float(entry['centre'][1]) - latitude), default=None)


def normalised(longitude):
    return (longitude + 180.0) % 360.0 - 180.0


def window_grid(entry, request, clip=None):
    """The Earth window of a map request: pixel size (degrees, from the map's metres per pixel on the
    Earth's radius), the column and row index ranges [i0, i1) and [j0, j1) of pixels counted from the
    entry's centre (column i spans centre longitude + [i, i + 1] pixels, row j spans centre latitude -
    [j + 1, j] pixels), the map's pixel count cut at the poles and, with clip, to the pixels wholly
    inside the entry's box, and the extent in degrees (west > east across the date line). clip defaults
    to the source being 'etopo15' (only the tiles meeting the box are downloaded)."""
    pixel_deg = math.degrees(request['metres_per_pixel'] / EARTH_RADIUS_M)
    longitude, latitude = (float(value) for value in entry['centre'])
    west, south, east, north = (float(value) for value in entry['box'])
    west_unwrapped = longitude - (longitude - west) % 360.0
    width = (east - west) % 360.0 or 360.0
    if longitude - west_unwrapped > width + EPSILON or not south <= latitude <= north:
        raise ValueError(f"Earth window {entry['name']}: its centre is outside its box")
    if clip is None:
        clip = source_for(request['metres_per_pixel']) == 'etopo15'
    if not clip:
        west_unwrapped, width, south, north = longitude - 180.0, 360.0, -90.0, 90.0
    columns, rows = request['columns'], request['rows']
    i0 = max(-(columns // 2), math.ceil((west_unwrapped - longitude) / pixel_deg - EPSILON))
    i1 = min(columns - columns // 2, math.floor((west_unwrapped + width - longitude) / pixel_deg + EPSILON))
    j0 = max(-(rows // 2), math.ceil((latitude - north) / pixel_deg - EPSILON))
    j1 = min(rows - rows // 2, math.floor((latitude - south) / pixel_deg + EPSILON))
    if i1 <= i0 or j1 <= j0:
        raise ValueError(f"Earth window {entry['name']}: no pixel of {request['metres_per_pixel']:.0f} m fits its box")
    east_edge = normalised(longitude + i1 * pixel_deg)
    extent = {'west': normalised(longitude + i0 * pixel_deg), 'south': latitude - j1 * pixel_deg,
              'east': 180.0 if east_edge == -180.0 else east_edge, 'north': latitude - j0 * pixel_deg}
    return {'pixel_deg': pixel_deg, 'centre': [longitude, latitude], 'columns_index': [i0, i1], 'rows_index': [j0, j1],
            'pixels': [i1 - i0, j1 - j0], 'extent': extent}


def edges(grid):
    """(longitude edges west to east, unwrapped from the centre; latitude edges north to south) of a grid."""
    (i0, i1), (j0, j1) = grid['columns_index'], grid['rows_index']
    longitude, latitude = grid['centre']
    return (longitude + np.arange(i0, i1 + 1) * grid['pixel_deg'], latitude - np.arange(j0, j1 + 1) * grid['pixel_deg'])


def _integrated(values, positions, axis):
    """Sums of values (piecewise constant per cell along axis) between consecutive fractional cell
    positions, from the cumulative sum interpolated at each position."""
    count = values.shape[axis]
    cumulative = np.concatenate([np.zeros_like(np.take(values, [0], axis=axis)), np.cumsum(values, axis=axis)], axis=axis)
    index = np.minimum(np.floor(positions).astype(int), count - 1)
    fraction = positions - index
    if axis == 1:
        at = cumulative[:, index] + fraction[None, :] * values[:, index]
    else:
        at = cumulative[index, :] + fraction[:, None] * values[index, :]
    return np.diff(at, axis=axis)


def box_average(values, extent, longitude_edges, latitude_edges):
    """Area-weighted mean of a source raster (north row first; extent {west, north, pixel_deg} of its
    outer edges, west > east across the date line) over target pixels between the given edges (longitudes
    west to east, unwrapped; latitudes north to south). NaN source pixels are left out; a target pixel
    with no finite source pixel is NaN. Weights are in degrees, as the sources' pixels are."""
    values = np.asarray(values, dtype=float)
    rows, columns = values.shape
    step = float(extent['pixel_deg'])
    longitude_edges, latitude_edges = np.asarray(longitude_edges, float), np.asarray(latitude_edges, float)
    west = float(extent['west'])
    # Unwrap the source's west edge to the window's: the nearest turn for a source narrower than the globe
    # (its edge may lie a pixel east of the window's), the turn at or west of the window's for the globe.
    turns = (longitude_edges[0] - west) / 360.0
    west += 360.0 * (math.floor(turns + EPSILON) if columns * step >= 360.0 - EPSILON else round(turns))
    x = np.clip((longitude_edges - west) / step, 0.0, columns)
    y = np.clip((float(extent['north']) - latitude_edges) / step, 0.0, rows)
    finite = np.isfinite(values)
    sums = _integrated(_integrated(np.where(finite, values, 0.0), x, 1), y, 0)
    weights = _integrated(_integrated(finite.astype(float), x, 1), y, 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(weights > 1e-6, sums / np.where(weights > 0, weights, 1.0), np.nan)


def code_revision():
    completed = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True)
    return completed.stdout.strip() if completed.returncode == 0 else None


def cache_directory(name, metres_per_pixel, root=None):
    return (CACHE if root is None else Path(root)) / name / f'{round(metres_per_pixel)}m'


def etopo_reader():
    """uwtools.etopo.read_box (S1.1), imported when a window is first read."""
    from uwtools import etopo
    return etopo.read_box


def etopo_files(source, extent):
    """The ETOPO files a window reads (uwtools.etopo): the 60-arc-second globe, or the 15-arc-second tiles
    meeting the extent."""
    from uwtools import etopo
    if source == 'etopo60':
        return [etopo.GLOBE_60]
    return [etopo.tile_path(name) for name in etopo.tiles_for_box(extent['west'], extent['south'], extent['east'], extent['north'])]


def file_signature(paths):
    """[{path, bytes, mtime_ns}] of the source files (None for a missing one): a cache is stale when they change."""
    signature = []
    for path in paths:
        path = Path(path)
        stat = path.stat() if path.is_file() else None
        signature.append({'path': path.as_posix(), 'bytes': stat.st_size if stat else None,
                          'mtime_ns': stat.st_mtime_ns if stat else None})
    return signature


def cached_window(directory, source, grid, sources):
    """(elevation, manifest) of a cached window of this CACHE_VERSION that covers the grid on the same anchor
    and pixel size, read from source files of the same signature (sources), else None."""
    path = directory / 'manifest.json'
    if not path.is_file() or not (directory / 'elevation.npy').is_file():
        return None
    manifest = json.loads(path.read_text(encoding='utf-8'))
    (i0, i1), (j0, j1) = grid['columns_index'], grid['rows_index']
    (c0, c1), (r0, r1) = manifest.get('columns_index', (0, 0)), manifest.get('rows_index', (0, 0))
    same = (manifest.get('schema') == WINDOW_SCHEMA and manifest.get('cache_version') == CACHE_VERSION
            and manifest.get('source') == source and manifest.get('source_files') == sources
            and manifest.get('centre') == grid['centre']
            and abs(float(manifest.get('pixel_deg', 0.0)) - grid['pixel_deg']) <= 1e-12 * grid['pixel_deg'])
    if not (same and c0 <= i0 and i1 <= c1 and r0 <= j0 and j1 <= r1):
        return None
    elevation = np.load(directory / 'elevation.npy')
    return elevation[j0 - r0:j1 - r0, i0 - c0:i1 - c0], manifest


def read_window(entry, request, reader=None, cache_root=None, files=None):
    """(elevation float32 north row first, record) of an entry's window for a map request, from the cache
    when it covers the request, else read with reader (default uwtools.etopo.read_box: (west, south, east,
    north, source) -> (array, extent {west, south, east, north, pixel_deg})), box-averaged and cached.
    files: (source, extent) -> the source files read (default etopo_files), whose sizes and modification
    times the cache records."""
    metres = request['metres_per_pixel']
    source = source_for(metres)
    grid = window_grid(entry, request)
    directory = cache_directory(entry['name'], metres, cache_root)
    sources = file_signature((files or etopo_files)(source, grid['extent']))
    record = {'schema': WINDOW_SCHEMA, 'cache_version': CACHE_VERSION, 'entry': entry, 'source': source,
              'source_label': SOURCE_LABELS[source], 'source_files': sources, 'metres_per_pixel': metres,
              'earth_radius_m': EARTH_RADIUS_M, **grid, 'map_pixels': [request['columns'], request['rows']]}
    found = cached_window(directory, source, grid, sources)
    if found is not None:
        elevation, manifest = found
        return elevation, dict(record, cache=directory.as_posix(), reused=True, source_extent=manifest.get('source_extent'),
                               code=manifest.get('code'))
    extent = grid['extent']
    array, source_extent = (reader or etopo_reader())(extent['west'], extent['south'], extent['east'], extent['north'], source)
    longitudes, latitudes = edges(grid)
    elevation = box_average(array, source_extent, longitudes, latitudes).astype(np.float32)
    manifest = {key: value for key, value in record.items() if key != 'map_pixels'}
    manifest.update({'source_extent': {key: float(value) for key, value in dict(source_extent).items()},
                     'method': 'area-weighted mean of the source pixels over each window pixel (degrees); NaN source pixels left out',
                     'code': {'module': 'scripts/evaluation/earth_windows.py', 'schema': WINDOW_SCHEMA, 'revision': code_revision()},
                     'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')})
    directory.mkdir(parents=True, exist_ok=True)
    np.save(directory / 'elevation.npy', elevation)
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    return elevation, dict(record, cache=directory.as_posix(), reused=False, source_extent=manifest['source_extent'],
                           code=manifest['code'])


OCEAN_SEED = (-160.0, 0.0)  # (longitude, latitude) in the open Pacific


def ocean_label(below, extent, seed=OCEAN_SEED):
    """uint8 mask of the cells of a global grid's `below` (rows north first, columns east from extent's
    west) 4-connected to the seed cell, the west and east columns joined across the date line."""
    from scipy import ndimage
    labels, _ = ndimage.label(below)  # the default structure is 4-connected
    rows, columns = below.shape
    step = float(extent['pixel_deg'])
    row = min(max(int(math.floor((float(extent['north']) - seed[1]) / step)), 0), rows - 1)
    column = int(math.floor(((seed[0] - float(extent['west'])) % 360.0) / step)) % columns
    if not below[row, column]:
        raise ValueError(f'The ocean seed {seed} is not below 0 m on the globe')
    joined = {}
    for west, east in zip(labels[:, 0], labels[:, -1]):
        if west and east and west != east:
            joined.setdefault(int(west), set()).add(int(east))
            joined.setdefault(int(east), set()).add(int(west))
    ocean, frontier = {int(labels[row, column])}, [int(labels[row, column])]
    while frontier:
        for other in joined.get(frontier.pop(), ()):
            if other not in ocean:
                ocean.add(other)
                frontier.append(other)
    return np.isin(labels, sorted(ocean)).astype(np.uint8)


def read_globe_below():
    """(cells below 0 m of the 60-arc-second globe, extent {west, north, pixel_deg}), read in strips."""
    import rasterio
    from rasterio.windows import Window
    from uwtools import etopo
    with rasterio.open(etopo.GLOBE_60) as dataset:
        below = np.zeros((dataset.height, dataset.width), dtype=bool)
        for start in range(0, dataset.height, 1024):
            height = min(1024, dataset.height - start)
            data = dataset.read(1, window=Window(0, start, dataset.width, height))
            valid = data > etopo.NODATA if dataset.nodata is None else (data != dataset.nodata) & (data > etopo.NODATA)
            below[start:start + height] = valid & (data < 0)
        transform = dataset.transform
    return below, {'west': float(transform.c), 'north': float(transform.f), 'pixel_deg': float(transform.a)}


def global_ocean(cache_root=None, globe=None, files=None):
    """(ocean mask uint8 of the 60-arc-second globe, its extent {west, north, pixel_deg}), from the cache
    <cache root>/ocean-60s/ (ocean.npy, manifest.json) when its CACHE_VERSION and source file signature
    match, else computed with ocean_label from globe() (default read_globe_below) and cached."""
    directory = (CACHE if cache_root is None else Path(cache_root)) / 'ocean-60s'
    sources = file_signature(files() if files else etopo_files('etopo60', None))
    path = directory / 'manifest.json'
    if path.is_file() and (directory / 'ocean.npy').is_file():
        manifest = json.loads(path.read_text(encoding='utf-8'))
        if manifest.get('cache_version') == CACHE_VERSION and manifest.get('source_files') == sources:
            return np.load(directory / 'ocean.npy', mmap_mode='r'), manifest['extent']
    below, extent = (globe or read_globe_below)()
    mask = ocean_label(below, extent)
    directory.mkdir(parents=True, exist_ok=True)
    np.save(directory / 'ocean.npy', mask)
    manifest = {'schema': 'uw-earth-ocean/v1', 'cache_version': CACHE_VERSION, 'source_files': sources,
                'extent': extent, 'shape': list(mask.shape), 'seed': list(OCEAN_SEED),
                'method': 'cells below 0 m 4-connected to the seed, the date line joined (scipy.ndimage.label)',
                'code': {'module': 'scripts/evaluation/earth_windows.py', 'revision': code_revision()},
                'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return mask, extent


def window_sea(elevation, record, ocean):
    """Pixels below 0 m whose centre lies in a cell of the ocean mask ((mask, extent) of global_ocean);
    every pixel below 0 m when ocean is None."""
    below = np.asarray(elevation, dtype=float) < 0.0
    if ocean is None:
        return below
    mask, extent = ocean
    longitudes, latitudes = edges(record)
    step = float(extent['pixel_deg'])
    centres = 0.5 * (longitudes[:-1] + longitudes[1:])
    columns = np.floor(((centres - float(extent['west'])) % 360.0) / step).astype(int) % mask.shape[1]
    rows = np.clip(np.floor((float(extent['north']) - 0.5 * (latitudes[:-1] + latitudes[1:])) / step).astype(int),
                   0, mask.shape[0] - 1)
    return below & (np.asarray(mask[np.ix_(rows, columns)]) != 0)


def render_window(elevation, record, style, vertical_exaggeration=RELIEF_EXAGGERATION, ocean=None):
    """(RGB, legend note) of an Earth window: sea (window_sea: the world ocean of `ocean`) and land (lake
    beds included) on tectonics:elevation as regional_terrain draws them, times the hillshade of
    maps.regional_relief in the relief style, and the landmass outline in black."""
    elevation = np.asarray(elevation, dtype=float)
    sea = window_sea(elevation, record, ocean)
    rgb = terrain_base(elevation, sea)
    note = ('Earth (ETOPO 2022 surface) on the elevation scale; sea = the world ocean (below 0 m and connected to it '
            'on the 60 arc-second globe); lake beds as land; no rivers or ice; black: the landmass outline')
    if style == 'relief':
        dy, dx = pixel_sizes_m(record['extent']['north'], record['pixel_deg'], elevation.shape[0], EARTH_RADIUS_M)
        rgb = rgb * shade_factor(hillshade(elevation, dy, dx, vertical_exaggeration))[..., None]
    elif style != 'terrain':
        raise ValueError(f'Unknown Earth window style {style!r}; known: terrain, relief')
    rgb[shore(~sea)] = 0.0
    return np.clip(np.rint(rgb), 0, 255).astype(np.uint8), note


def caption(record):
    """The Earth column's caption; an entry drawn outside its metres_per_pixel range (by --earth) says so."""
    low, high = (float(value) for value in record['entry']['metres_per_pixel'])
    metres = record['metres_per_pixel']
    text = f"Earth: {record['entry']['title']} (ETOPO 2022 {record['source_label']}, {metres:,.0f} m/px)"
    if not low <= metres <= high:
        text += f' (outside the range where its features show: {low:,.0f}–{high:,.0f} m/px)'
    return text


def earth_specs(profile):
    """The profile maps an Earth window stands beside: composites of STYLES, restricted to the names of
    the profile's earth.maps when it lists them."""
    names = profile['earth'].get('maps')
    return [spec for spec in profile.get('maps', []) if spec.get('composite') in STYLES and (names is None or spec['name'] in names)]


def render_maps(bundle, profile, directory, name=None, reader=None, catalogue=None, cache_root=None, tag='', files=None,
                ocean=None):
    """(map entries of kind "earth", the page's earth block) of a regional bundle for a profile with an
    "earth" block ({"catalogue": path, "maps": [names]}). The PNGs go to directory/maps/<map>-earth<tag>.png,
    enlarged by the same whole factor as the map beside them. The block names the entry, how it was chosen,
    the source, the pixel size and count, the extent, the cache and the features to look for, or says why
    there is no Earth window."""
    config = profile.get('earth')
    if not config:
        return [], None
    if bundle.manifest.get('window') is None:
        return [], {'entry': None, 'note': 'No Earth window: the bundle is not a regional map (no window block).'}
    request = map_request(bundle.manifest)
    metres = request['metres_per_pixel']
    if catalogue is None:
        path = ROOT / config.get('catalogue', CATALOGUE)
        if not path.is_file():
            return [], {'entry': None, 'note': f'No Earth window: the catalogue {path.relative_to(ROOT).as_posix()} does not exist.'}
        catalogue = load_catalogue(path)
    entry = choose_entry(catalogue, metres, request['latitude'], name)
    if entry is None:
        return [], {'entry': None, 'note': f'No Earth window: no catalogue entry holds {metres:,.0f} m per pixel.'}
    try:
        elevation, record = read_window(entry, request, reader, cache_root, files)
        ocean = global_ocean(cache_root) if ocean is None else ocean
    except (ImportError, OSError, ValueError) as error:
        return [], {'entry': entry['name'], 'note': f"No Earth window: {entry['name']} could not be read ({type(error).__name__}: {error})."}
    title = caption(record)
    extent = tuple(record['extent'][key] for key in ('west', 'south', 'east', 'north'))
    factor = -(-1024 // request['columns'])  # the map's enlargement factor (map_rendering.frame's)
    written = []
    for spec in earth_specs(profile):
        style = STYLES[spec['composite']]
        rgb, note = render_window(elevation, record, style, spec.get('options', {}).get('vertical_exaggeration', RELIEF_EXAGGERATION),
                                  ocean)
        relative = Path('maps') / f"{spec['name']}-earth{tag}.png"
        (Path(directory) / relative).parent.mkdir(parents=True, exist_ok=True)
        image = map_rendering.frame(rgb, title, legend(spec['palette'], note), width=factor * record['pixels'][0], extent=extent)
        image.save(Path(directory) / relative)
        # The framed map beside it is wider by the columns the window lacks (the margins are the same):
        # the page draws the window at this fraction of the map's width, so both share one scale.
        fraction = image.width / (image.width + factor * (request['columns'] - record['pixels'][0]))
        written.append({'map': spec['name'], 'kind': 'earth', 'season': None, 'season_label': None, 'title': title,
                        'path': relative.as_posix(), 'palette': spec['palette'], 'earth': entry['name'], 'style': style,
                        'pixels': record['pixels'], 'display_fraction': fraction})
    if written:
        palettes.write_sidecar(Path(directory) / 'maps')
    cache = Path(record['cache'])
    info = {'entry': entry['name'], 'title': entry['title'], 'caption': title, 'chosen_by': '--earth' if name else NEAREST,
            'source': record['source'], 'source_label': record['source_label'], 'metres_per_pixel': metres,
            'pixels': record['pixels'], 'map_pixels': record['map_pixels'], 'extent': record['extent'],
            'cache': cache.relative_to(ROOT).as_posix() if cache.is_relative_to(ROOT) else cache.as_posix(),
            'reused': record['reused'], 'features': entry.get('features', []),
            'note': f"Earth window {entry['name']} ({entry['title']}), chosen by {'--earth' if name else NEAREST}."}
    return written, info
