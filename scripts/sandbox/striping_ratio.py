"""Striping evidence (not a fix): on the global domain's land, the median |dz/dx| over the median |dz/dy| (true km,
dx = dy cos latitude, central differences) by latitude band; Earth's ETOPO 2048 about 0.9 to 1.2 (roadmap v17).

Slope-aspect harmonics (item 3 step 0), by the same latitude bands on land: the aspect theta of the central-
difference gradient in true km on cells steeper than 0.1 m/km, A8 = |mean e^(8 i theta)| and A4 = mean
cos(4 theta) (both unchanged by a reflection across an axis, so the gradient's sign convention does not matter);
for a label's final surface (-elevation.npy, land as above), its procedural input (-elevation-procedural.npy,
land: its own sea (terrain_metrics.sea_from) and > 0 m) and its physics surface; with --earth, ETOPO at the
global's resolution (refs/processed/heightmap/<columns>.tif, terrain_metrics' loader, read only), sea by
terrain_metrics.sea_from. --json OUT writes the records.

A label is a directory under out/scratch/sandbox or an absolute sandbox directory. --bundle DIR --field SPEC[,SPEC]
measures bundle grid fields (T8.0): SPEC is a field id or a sum of ids joined by '+', each term in metres above sea
level (units 'm above datum' minus the bundle's geosphere.sea_level; geosphere.tectonic_height, sea level at 1,
times the metres per engine unit of the snapshot's heightmap_stats; other units as they are), land as for the
procedural input. --earth-columns C1,C2 measures ETOPO at those columns: refs/processed/heightmap/<C>.tif, or
where none exists the area mean of the 2048 raster (uwtools.grid.remap_area_mean); written as earth-<C>.

  uv run python scripts/sandbox/striping_ratio.py LABEL1,LABEL2 [--earth] [--earth-columns C1,C2]
      [--bundle DIR --field SPEC[,SPEC]] [--json OUT]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path('C:/dev/undiscovered-worlds')
sys.path.insert(0, str(ROOT / 'scripts' / 'sandbox'))
import terrain_metrics as tm  # noqa: E402

BANDS = ((0, 30), (30, 50), (50, 60), (60, 70), (70, 80))
MIN_SLOPE_M_KM = 0.1


def gradients(z, lat, dy):
    """(d z / d east, d z / d north) in m/km, central differences, east-west wrapping."""
    gn = -np.gradient(z, axis=0) / dy
    ge = ((np.roll(z, -1, 1) - np.roll(z, 1, 1)) / 2.0) / (dy * np.cos(np.radians(lat))[:, None])
    return ge, gn


def measure(z, land, lat, dy):
    ge, gn = gradients(z, lat, dy)
    gx, gy = np.abs(ge), np.abs(gn)
    slope = np.hypot(ge, gn)
    theta = np.arctan2(gn, ge)
    out = {}
    for lo, hi in BANDS:
        band = land & (np.abs(lat)[:, None] >= lo) & (np.abs(lat)[:, None] < hi)
        if band.sum() <= 100:
            continue
        steep = band & (slope > MIN_SLOPE_M_KM)
        t = theta[steep]
        out[f'{lo}-{hi}'] = {'cells': int(band.sum()), 'aspect_cells': int(steep.sum()),
                             'ratio': round(float(np.median(gx[band]) / max(np.median(gy[band]), 1e-9)), 3),
                             'a8': round(float(np.abs(np.mean(np.exp(8j * t)))), 4) if t.size else None,
                             'a4': round(float(np.mean(np.cos(4.0 * t))), 4) if t.size else None}
    steep = land & (np.abs(lat)[:, None] < 80) & (slope > MIN_SLOPE_M_KM)
    t = theta[steep]
    out['all-0-80'] = {'aspect_cells': int(steep.sum()), 'a8': round(float(np.abs(np.mean(np.exp(8j * t)))), 4),
                       'a4': round(float(np.mean(np.cos(4.0 * t))), 4)}
    return out


def earth(columns):
    from PIL import Image
    from uwtools.tectonics_earth_reference import reference_inputs
    path = reference_inputs(ROOT, columns)[0]
    derived = None
    if not path.is_file():
        from uwtools import grid
        path = reference_inputs(ROOT, 2048)[0]
        derived = f'area mean of the 2048 raster to {columns} x {columns // 2} (uwtools.grid.remap_area_mean)'
    with Image.open(path) as image:
        z = np.array(image).astype(np.float64)
    if derived:
        z = grid.remap_area_mean(z, columns // 2, columns)
    rows = z.shape[0]
    lat = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    dy = tm.RADIUS_KM * math.radians(180.0 / rows)
    land = ~tm.sea_from(z) & (z > 0)
    record = {'source': path.relative_to(ROOT).as_posix(), 'shape': list(z.shape)}
    if derived:
        record['derived'] = derived
    record['bands'] = measure(z, land, lat, dy)
    return record


def bundle_term(directory, manifest, field_id):
    """(field in metres above sea level, conversion text) of one bundle grid field."""
    entries = {f['id']: f for f in manifest['fields']}
    if field_id not in entries:
        raise SystemExit(f'{field_id}: not a field of {directory}')
    entry = entries[field_id]
    if entry['layout'] != 'grid':
        raise SystemExit(f'{field_id}: layout {entry["layout"]}, not a grid field')
    a = np.load(directory / entry['file']).astype(np.float64)
    if field_id == 'geosphere.tectonic_height':
        stats = json.loads((directory / 'tables/geosphere/tectonic_snapshot.json').read_text())['heightmap_stats']
        k = (stats['metric_max'] - stats['metric_min']) / (stats['internal_max'] - stats['internal_min'])
        return (a - 1.0) * k, f'({field_id} - 1) * {k:.1f} m per engine unit (snapshot heightmap_stats)'
    if entry['units'] == 'm above datum':
        sea = float(np.load(directory / entries['geosphere.sea_level']['file']))
        return a - sea, f'{field_id} - geosphere.sea_level ({sea:.0f} m)'
    return a, f'{field_id} ({entry["units"] or "no units"}, as stored)'


def bundle_field(directory, spec):
    manifest = json.loads((directory / 'manifest.json').read_text())
    terms = [bundle_term(directory, manifest, field_id) for field_id in spec.split('+')]
    z = sum(t[0] for t in terms)
    rows = z.shape[0]
    lat = 90.0 - (np.arange(rows) + 0.5) * 180.0 / rows
    dy = tm.RADIUS_KM * math.radians(180.0 / rows)
    land = ~tm.sea_from(z) & (z > 0)
    return {'conversion': ' + '.join(t[1] for t in terms), 'shape': list(z.shape),
            'land_share': round(float(land.mean()), 4), 'bands': measure(z, land, lat, dy)}


def option(args, name):
    if name not in args:
        return None
    i = args.index(name)
    value = args[i + 1]
    del args[i:i + 2]
    return value


def main():
    args = sys.argv[1:]
    out = option(args, '--json')
    bundle = option(args, '--bundle')
    fields = option(args, '--field')
    earth_columns = option(args, '--earth-columns')
    with_earth = '--earth' in args
    if with_earth:
        args.remove('--earth')
    result = {'min_slope_m_km': MIN_SLOPE_M_KM}
    columns = None
    if bundle:
        directory = Path(bundle)
        result['bundle'] = {'path': directory.as_posix(), 'fields': {}}
        for spec in (fields or 'geosphere.tectonic_height').split(','):
            record = bundle_field(directory, spec)
            result['bundle']['fields'][spec] = record
            print(spec, '; '.join(f'{band}: ' + (f'ratio {v["ratio"]:.2f} ' if 'ratio' in v else '')
                                  + f'A8 {v["a8"]} A4 {v["a4"]}' for band, v in record['bands'].items()))
    for label in (args[0].split(',') if args else []):
        d = Path(label) if Path(label).is_absolute() else tm.SANDBOX / label
        config = json.loads((d / 'config.json').read_text())
        z = np.load(d / 'global-elevation.npy').astype(float)
        phys = np.load(d / 'global-elevation-physics.npy').astype(float)
        pre = np.load(d / 'global-elevation-procedural.npy').astype(float)
        columns = z.shape[1]
        area, lat, dy = tm.geometry(config, 'global', z.shape)
        land = ~tm.sea_from(phys, seed=tm.sea_from(pre)) & (z > 0)
        record = {'final': measure(z, land, lat, dy),
                  'procedural': measure(pre, ~tm.sea_from(pre) & (pre > 0), lat, dy),
                  'physics': measure(phys, land, lat, dy)}
        result[d.name] = record
        print(d.name,'; '.join(f'{band}: {v["ratio"]:.2f}' for band, v in record['final'].items() if 'ratio' in v))
        for surface, bands in record.items():
            print(' ', surface, '; '.join(f'{band}: A8 {v["a8"]} A4 {v["a4"]}' for band, v in bands.items()))
    if with_earth:
        result['earth'] = earth(columns or 2048)
        print('earth', '; '.join(f'{band}: A8 {v["a8"]} A4 {v["a4"]}' + (f' ratio {v["ratio"]:.2f}' if 'ratio' in v else '')
                                 for band, v in result['earth']['bands'].items()))
    for c in (earth_columns.split(',') if earth_columns else []):
        record = earth(int(c))
        result[f'earth-{c}'] = record
        print(f'earth-{c}', '; '.join(f'{band}: A8 {v["a8"]} A4 {v["a4"]}' + (f' ratio {v["ratio"]:.2f}' if 'ratio' in v else '')
                                      for band, v in record['bands'].items()))
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
