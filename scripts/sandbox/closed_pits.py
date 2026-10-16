"""Item 1 (T1.6): closed depressions on a sandbox domain's final surface against its physics surface.

  uv run python scripts/sandbox/closed_pits.py LABEL DOMAIN [LABEL2 ...]

Sea: <=0 m components of the physics surface touching the procedural input's sea (terrain_metrics.sea_from).
Pits: Priority-Flood with the sea and (variant 'edge') the crop's border cells as outlets; a pit cell is land
with fill depth > 1 m. Lakes: finite -lake-surface cells. Reports area shares (of land), depth p50/p99/max,
lake overlap, and the filter's change (final - physics) on pit cells.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path('C:/dev/undiscovered-worlds')
sys.path.insert(0, str(ROOT / 'scripts' / 'sandbox'))
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'scripts' / 'evaluation'))
import terrain_metrics as tm  # noqa: E402
import terrain_sandbox as ts  # noqa: E402


def flood(z, sea, edge_outlet):
    outlet = sea.copy()
    if edge_outlet:
        outlet[0, :] = outlet[-1, :] = outlet[:, 0] = outlet[:, -1] = True
    filled = ts.priority_flood(z, outlet, False)
    return np.where(outlet, 0.0, filled - z)


def stats(depth, land, area, lake=None):
    pit = land & (depth > 1.0)
    out = {'pit_pct_of_land': round(100 * area[pit].sum() / area[land].sum(), 2),
           'depth_p50_m': round(float(np.percentile(depth[pit], 50)), 1) if pit.any() else 0,
           'depth_p99_m': round(float(np.percentile(depth[pit], 99)), 1) if pit.any() else 0,
           'depth_max_m': round(float(depth[pit].max()), 1) if pit.any() else 0,
           'pits': int(ndimage.label(pit, structure=np.ones((3, 3)))[1])}
    if lake is not None:
        out['pit_in_lake_pct'] = round(100 * area[pit & lake].sum() / max(area[pit].sum(), 1e-9), 1)
        out['pit_not_lake_pct_of_land'] = round(100 * area[pit & ~lake].sum() / area[land].sum(), 2)
    return out


def main():
    domain = sys.argv[2]
    for label in [sys.argv[1]] + sys.argv[3:]:
        d = tm.SANDBOX / label
        config = json.loads((d / 'config.json').read_text())
        final = np.load(d / f'{domain}-elevation.npy').astype(float)
        physics = np.load(d / f'{domain}-elevation-physics.npy').astype(float)
        pre = np.load(d / f'{domain}-elevation-procedural.npy').astype(float)
        surface = np.load(d / f'{domain}-lake-surface.npy')
        lake = np.isfinite(surface)
        area, lat, dy = tm.geometry(config, domain, final.shape)
        sea = tm.sea_from(physics, seed=tm.sea_from(pre))
        land = ~sea & (final > 0)
        result = {'label': label, 'domain': domain, 'lake_pct_of_land': round(100 * area[lake].sum() / area[~sea].sum(), 2),
                  'lakes': int(ndimage.label(lake, structure=np.ones((3, 3)))[1])}
        for edge in (False, True):
            key = 'edge-outlet' if edge else 'sea-only'
            result[key] = {'physics': stats(flood(physics, sea, edge), land, area),
                           'final': stats(flood(final, sea, edge), land, area, lake)}
        depth = flood(final, sea, True)
        pit = land & (depth > 1.0) & ~lake
        change = final - physics
        result['filter_change_m'] = {'land_p1': round(float(np.percentile(change[land], 1)), 1),
                                     'land_p99': round(float(np.percentile(change[land], 99)), 1),
                                     'on_new_pits_p50': round(float(np.percentile(change[pit], 50)), 1) if pit.any() else 0,
                                     'on_new_pits_p5': round(float(np.percentile(change[pit], 5)), 1) if pit.any() else 0}
        # Where: distance to the sea (km) and the physics slope at the pit cells.
        dist = ndimage.distance_transform_edt(~sea) * dy
        gy, gx = np.gradient(physics, dy * 1000.0)
        slope = np.hypot(gx, gy)
        if pit.any():
            result['new_pits_where'] = {'dist_sea_km_p50': round(float(np.median(dist[pit])), 1),
                                        'physics_elev_p50_m': round(float(np.median(physics[pit])), 0),
                                        'physics_slope_p50': round(float(np.median(slope[pit])), 3),
                                        'land_slope_p50': round(float(np.median(slope[land])), 3)}
        print(json.dumps(result, indent=1))


if __name__ == '__main__':
    main()
