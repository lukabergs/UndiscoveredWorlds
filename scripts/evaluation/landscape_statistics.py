"""Reported numbers of the landscape modules (erosion plan, Phase C, step C0).

Two validation kinds that measure a bundle without an observed reference map, so that every
step of Phase C reports the same numbers from the evaluation instead of from ad hoc scripts:

- drainage_pattern: the routing's direction statistics (the share of land draining due east or
  west, the share continuing in its receiver's direction), lakes, closed lakes, endorheic land,
  deltas, wetlands and the land-mean denudation rate. Context, not a target: the numbers are
  compared with the baseline run.
- ocean_depths: the area-weighted share of the ocean by depth band against Earth's (ETOPO 2022,
  the hypsometry reference's inputs), with the mean and the median depth.
"""
import sys
from pathlib import Path

import numpy as np
from PIL import Image

from checks import DIRECTION_OFFSETS, has_neighbour, receiver_values, whole_metres

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from uwtools import grid  # noqa: E402
from uwtools.tectonics_earth_reference import reference_inputs  # noqa: E402

OCEAN_BANDS_KM = (0.0, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, np.inf)
DRAINAGE_FIELDS = ('hydrosphere.flow_direction', 'geosphere.base_elevation', 'geosphere.sea_level')


def band_label(lower, upper):
    return f'{lower:g}-{upper:g} km' if np.isfinite(upper) else f'over {lower:g} km'


def depth_shares(depth_m, weights):
    """Percent of the weighted area per band of OCEAN_BANDS_KM (flat keys, as the report tables want),
    with the mean and median depth (m)."""
    depth_m, weights = np.asarray(depth_m, dtype=float).ravel(), np.asarray(weights, dtype=float).ravel()
    total = weights.sum()
    shares = {}
    for lower, upper in zip(OCEAN_BANDS_KM[:-1], OCEAN_BANDS_KM[1:]):
        inside = (depth_m >= 1000 * lower) & (depth_m < 1000 * upper)
        shares[f'ocean {band_label(lower, upper)} percent'] = float(100 * weights[inside].sum() / total) if total else None
    order = np.argsort(depth_m, kind='stable')
    cumulative = np.cumsum(weights[order])
    median = float(depth_m[order][np.searchsorted(cumulative, cumulative[-1] / 2)]) if total else None
    return {**shares, 'mean_depth_m': float(np.average(depth_m, weights=weights)) if total else None,
            'median_depth_m': median}


def earth_ocean_depths(width):
    """Earth's ocean depth bands from the hypsometry reference's heightmap and non-ocean mask."""
    heightmap, mask, *_ = reference_inputs(ROOT, width)
    with Image.open(heightmap) as image:
        elevation = np.array(image, dtype=float)
    with Image.open(mask) as image:
        ocean = np.array(image) == 0
    weights = grid.area_weights(*elevation.shape)
    return depth_shares(-elevation[ocean], weights[ocean])


def measure_ocean_depths(bundle):
    from validation import surface_elevation
    sea = bundle.sea_mask()
    weights = bundle.area_weights()
    return depth_shares(-surface_elevation(bundle)[sea], weights[sea])


def compare_ocean_depths(bundle, spec, baseline=None):
    width = int(spec.get('reference_width', 2048))
    earth = earth_ocean_depths(width)
    result = {'reference': f'Earth ocean depth by band, ETOPO 2022 at {width} x {width // 2} (the hypsometry inputs)',
              'earth': earth}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is not None:
            measured = measure_ocean_depths(source)
            result[key] = {'metrics': measured, 'earth_comparison': {
                f'{name} minus Earth': (value - earth[name] if value is not None and earth[name] is not None else None)
                for name, value in measured.items()}}
    return result


def land_mask(bundle):
    """Land: base elevation in whole metres above sea level (planet::basebelowsea)."""
    return whole_metres(bundle.load('geosphere.base_elevation').astype(float)) > float(bundle.scalar('geosphere.sea_level'))


def measure_drainage_pattern(bundle):
    land = land_mask(bundle)
    direction = bundle.load('hydrosphere.flow_direction').astype(int)
    shape = direction.shape
    weights = bundle.area_weights()
    land_cells = int(np.count_nonzero(land))
    metrics = {'land_cells': land_cells}
    if land_cells:
        east_west = land & ((direction == 3) | (direction == 7))
        metrics['share_east_west'] = float(np.count_nonzero(east_west) / land_cells)
        # A land cell continues straight when its receiver is a land cell draining in the same direction.
        land_receiver = np.zeros(shape, dtype=bool)
        straight = np.zeros(shape, dtype=bool)
        for code, (dx, dy) in DIRECTION_OFFSETS.items():
            mine = land & (direction == code) & has_neighbour(shape, dy)
            receiver_land = receiver_values(land, dx, dy)
            receiver_direction = receiver_values(direction, dx, dy)
            land_receiver |= mine & receiver_land
            straight |= mine & receiver_land & (receiver_direction == code)
        receivers = int(np.count_nonzero(land_receiver))
        metrics['cells_with_land_receiver'] = receivers
        metrics['share_straight'] = float(np.count_nonzero(straight) / receivers) if receivers else None
    if bundle.has('hydrosphere.lake_surface'):
        metrics['lake_cells'] = int(np.count_nonzero(bundle.load('hydrosphere.lake_surface') != 0))
    if bundle.has('hydrosphere.special_features'):
        special = bundle.load('hydrosphere.special_features')
        metrics['salt_lake_cells'] = int(np.count_nonzero(special == 100))
        metrics['salt_pan_cells'] = int(np.count_nonzero(special == 110))
        metrics['wetland_cells'] = int(np.count_nonzero((special >= 130) & (special <= 132)))
    if bundle.has('hydrosphere.basin_class'):
        metrics['endorheic_cells'] = int(np.count_nonzero(land & (bundle.load('hydrosphere.basin_class') == 2)))
    if bundle.has('hydrosphere.delta_direction'):
        metrics['delta_cells'] = int(np.count_nonzero(bundle.load('hydrosphere.delta_direction') != 0))
    if bundle.has('geosphere.erosion_rate') and land_cells:
        rate = bundle.load('geosphere.erosion_rate').astype(float)
        metrics['land_mean_erosion_rate_m_myr'] = float(np.mean(rate[land]))
        metrics['land_area_mean_erosion_rate_m_myr'] = float(np.average(rate[land], weights=weights[land]))
    return metrics


def compare_drainage_pattern(bundle, spec, baseline=None):
    result = {'reference': 'none: reported numbers of the routing and the hydrology, compared with the baseline run'}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is not None:
            result[key] = {'earth_comparison': {}, 'metrics': measure_drainage_pattern(source)}
    return result


def main():
    """Print the measurements of bundles (checkpoints included) as one JSON document per bundle."""
    import argparse
    import json
    from uwtools.bundle import Bundle
    parser = argparse.ArgumentParser(description='Reported landscape numbers of run bundles: '
                                     'uv run --offline scripts/evaluation/landscape_statistics.py '
                                     '{ocean-depths,drainage-pattern} BUNDLE [BUNDLE ...]')
    parser.add_argument('measure', choices=('ocean-depths', 'drainage-pattern'))
    parser.add_argument('bundles', nargs='+', type=Path)
    args = parser.parse_args()
    measure = measure_ocean_depths if args.measure == 'ocean-depths' else measure_drainage_pattern
    for path in args.bundles:
        print(json.dumps({'bundle': path.as_posix(), 'metrics': measure(Bundle.open(path))}))


if __name__ == '__main__':
    main()
