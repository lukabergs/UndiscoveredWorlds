"""Simple terrain metrics of windows on a common pixel grid (terrain roadmap T1.5): sandbox windows (nested_window.py's
outputs) and Earth reference windows (glo30_windows.py's), so a level can be compared with the Earth at the same pixel
size in km.

Per window, on its land (not sea; lakes at their surface, as drawn): the slope (|gradient| by central differences, the
window's north-south and east-west pixel sizes in metres) at the 50th and 90th percentiles; the local relief (max - min
of the surface, sea at 0 m, in a moving box of about 1 km and 5 km on a side, the pixel counts from the pixel sizes) at
the median and 90th percentile of the land cells; the land's hypsometry (percentiles 10, 25, 50, 75, 90 and the
maximum; shares of the land in 0-200, 200-500, 500-1000, 1000-1500, above 1500 m); the sea share of the window.

Items: LABEL/DOMAIN for out/scratch/sandbox/LABEL/DOMAIN-{elevation,sea,lake-surface}.npy and DOMAIN-grid.json, or
atlas:STEM for refs/processed/earth-atlas/STEM.npy, STEM-sea.npy and STEM.json.

  uv run python scripts/sandbox/window_metrics.py --items t15-nest-a/G-coast-nest,atlas:glo30-sognefjord-98km
      --out out/scratch/sandbox/compare-t15/metrics.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
SANDBOX = ROOT / 'out' / 'scratch' / 'sandbox'
ATLAS = ROOT / 'refs' / 'processed' / 'earth-atlas'
BANDS = ((0, 200), (200, 500), (500, 1000), (1000, 1500), (1500, None))


def load(item):
    if item.startswith('atlas:'):
        stem = item.split(':', 1)[1]
        elevation = np.load(ATLAS / f'{stem}.npy').astype(np.float64)
        sea = np.load(ATLAS / f'{stem}-sea.npy')
        grid = json.loads((ATLAS / f'{stem}.json').read_text(encoding='utf-8'))
        lake = None
    else:
        label, domain = item.split('/', 1)
        stem = SANDBOX / label / domain
        elevation = np.load(f'{stem}-elevation.npy').astype(np.float64)
        sea = np.load(f'{stem}-sea.npy')
        grid = json.loads(Path(f'{stem}-grid.json').read_text(encoding='utf-8'))
        lake_path = Path(f'{stem}-lake-surface.npy')
        lake = np.load(lake_path) if lake_path.is_file() else None
    if lake is not None:
        elevation = np.where(np.isfinite(lake), lake, elevation)
    return elevation, sea.astype(bool), float(grid['pixel_ns_km']), float(grid['pixel_ew_km'])


def relief(surface, rows, columns):
    size = (max(1, rows), max(1, columns))
    return ndimage.maximum_filter(surface, size=size, mode='nearest') - ndimage.minimum_filter(surface, size=size, mode='nearest')


def metrics(elevation, sea, ns_km, ew_km):
    land = ~sea
    surface = np.where(sea, 0.0, elevation)
    gy, gx = np.gradient(surface)
    slope = np.hypot(gx / (ew_km * 1000.0), gy / (ns_km * 1000.0))
    out = {'sea_share': round(float(sea.mean()), 4),
           'slope_p50': round(float(np.percentile(slope[land], 50)), 4),
           'slope_p90': round(float(np.percentile(slope[land], 90)), 4)}
    for km in (1, 5):
        values = relief(surface, int(round(km / ns_km)), int(round(km / ew_km)))[land]
        out[f'relief_{km}km_p50_m'] = round(float(np.percentile(values, 50)), 1)
        out[f'relief_{km}km_p90_m'] = round(float(np.percentile(values, 90)), 1)
    heights = elevation[land]
    for q in (10, 25, 50, 75, 90):
        out[f'land_p{q}_m'] = round(float(np.percentile(heights, q)), 1)
    out['land_max_m'] = round(float(heights.max()), 1)
    for low, high in BANDS:
        inside = (heights >= low) & ((heights < high) if high is not None else True)
        out[f'land_{low}_{high or "up"}_pct'] = round(100.0 * float(inside.mean()), 1)
    out['pixel_km'] = [round(ns_km, 4), round(ew_km, 4)]
    out['pixels'] = list(elevation.shape)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--items', required=True)
    parser.add_argument('--out', required=True)
    arguments = parser.parse_args()
    table = {}
    for item in arguments.items.split(','):
        table[item] = metrics(*load(item))
        print(item, table[item], flush=True)
    out = Path(arguments.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    document = {'schema': 'uw-window-metrics/v1', 'module': 'scripts/sandbox/window_metrics.py',
                'method': __doc__.split('\n\n')[1].replace('\n', ' '), 'items': table}
    out.write_text(json.dumps(document, indent=1) + '\n', encoding='utf-8')
    print(out)


if __name__ == '__main__':
    main()
