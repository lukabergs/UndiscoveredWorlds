"""Reported numbers of the landscape modules (erosion plan, Phase C, step C0).

Two validation kinds that measure a bundle without an observed reference map, so that every
step of Phase C reports the same numbers from the evaluation instead of from ad hoc scripts:

- drainage_pattern: the routing's direction statistics (the share of land draining due east or
  west, overall and in the latitude bands 0-30, 30-60 and 60-90 degrees; the share continuing in
  its receiver's direction; the land-mean angle between the receiver displacement and the downhill
  direction of the filled surface, step C3), lakes, closed lakes, endorheic land,
  deltas, wetlands and the land-mean denudation rate. Context, not a target: the numbers are
  compared with the baseline run.
- ocean_depths: the area-weighted share of the ocean by depth band against Earth's (ETOPO 2022,
  the hypsometry reference's inputs), with the mean and the median depth.
- glaciation: the glacial fields of the erosion run (erosion plan Phase G): the glaciated and the
  present-day ice-covered share of the land, the mean glaciation, glaciated sea cells (sea with
  glacial denudation) and their depths, glaciated cells below the sea level, glaciated and
  ice-covered lakes, the ELA above the sea level per 15-degree latitude band and the
  glacial share of the denudation. Reported numbers, no threshold.
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
from uwtools.paths import DATA_ROOT  # noqa: E402

OCEAN_BANDS_KM = (0.0, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, np.inf)
DRAINAGE_FIELDS = ('hydrosphere.flow_direction', 'geosphere.base_elevation', 'geosphere.sea_level')
GLACIATION_FIELDS = ('cryosphere.glaciation', 'cryosphere.ice_thickness', 'cryosphere.equilibrium_line_altitude',
                     'geosphere.glacial_denudation', 'geosphere.base_elevation', 'geosphere.sea_level')
ELA_BAND_DEGREES = 15


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
    heightmap, mask, *_ = reference_inputs(DATA_ROOT, width)
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


def direction_deviation(land, direction, filled, sea_level):
    """Land-mean absolute angle (degrees) between each land cell's receiver displacement and the
    local downhill direction of the filled surface (sea cells at the sea level), by central
    differences (one-sided on the polar rows), on the sphere metric of the routing: east-west
    spacing cos(latitude) per row (grid.latitudes), a diagonal at the mean spacing of its two rows,
    north-south spacing pi / rows; outlets point towards their sea cell. Cells without a direction
    or without a gradient are skipped; None when no cell remains."""
    rows, columns = direction.shape
    dx = np.cos(np.deg2rad(grid.latitudes(rows))) * 2.0 * np.pi / columns
    dy = np.pi / rows
    surface = np.where(land, filled, sea_level)
    gradient_east = (np.roll(surface, -1, axis=1) - np.roll(surface, 1, axis=1)) / (2.0 * dx[:, None])
    north = np.vstack([surface[:1], surface[:-1]])
    south = np.vstack([surface[1:], surface[-1:]])
    spacing = np.full(rows, 2.0 * dy)
    spacing[0] = spacing[-1] = dy
    gradient_north = (north - south) / spacing[:, None]
    east = np.zeros(direction.shape)
    northward = np.zeros(direction.shape)
    row_index = np.arange(rows)[:, None] * np.ones(direction.shape, dtype=int)
    for code, (ox, oy) in DIRECTION_OFFSETS.items():
        mine = direction == code
        if oy == 0:
            spacing_east = dx[row_index]
        else:
            spacing_east = 0.5 * (dx[row_index] + dx[np.clip(row_index + oy, 0, rows - 1)])
        east = np.where(mine, ox * spacing_east, east)
        northward = np.where(mine, -oy * dy, northward)
    downhill_east, downhill_north = -gradient_east, -gradient_north
    valid = land & (direction >= 1) & (direction <= 8) & ((downhill_east != 0) | (downhill_north != 0))
    if not np.any(valid):
        return None
    dot = east * downhill_east + northward * downhill_north
    cross = east * downhill_north - northward * downhill_east
    return float(np.mean(np.degrees(np.abs(np.arctan2(cross, dot)))[valid]))


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
        latitude = np.abs(grid.latitudes(shape[0]))[:, None] * np.ones(shape)
        for lower, upper in ((0, 30), (30, 60), (60, 90)):
            band = land & (latitude >= lower) & ((latitude < upper) if upper < 90 else (latitude <= upper))
            cells = int(np.count_nonzero(band))
            metrics[f'share_east_west_{lower}_{upper}'] = (float(np.count_nonzero(band & east_west) / cells)
                                                           if cells else None)
        if bundle.has('hydrosphere.filled_elevation'):
            metrics['direction_deviation_deg'] = direction_deviation(
                land, direction, bundle.load('hydrosphere.filled_elevation').astype(float),
                float(bundle.scalar('geosphere.sea_level')))
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


def latitude_name(degrees):
    return f'{abs(degrees)}{"S" if degrees < 0 else "N"}'


def measure_glaciation(bundle):
    """The glacial fields over the land of land_mask: shares by area (bundle.area_weights), glaciated
    sea cells (not land, no lake surface, glacial denudation > 0; their area from planetbody.radius,
    None without it; the minimum and median depth below the sea level and the count deeper than
    100 m), glaciated cells whose bed ends at or below the sea level (land_mask false), the largest
    glacial denudation, lake cells (a nonzero lake surface, the lake checks' rule) with glaciation
    and with present-day ice,
    the area-weighted mean ELA above the sea level per 15-degree band (None where a band has no
    land, and everywhere when the ELA is zero on every cell: the glacial term was off) and the land-mean glacial denudation with its share of the land-mean denudation."""
    land = land_mask(bundle)
    weights = bundle.area_weights()
    glaciation = bundle.load('cryosphere.glaciation').astype(float)
    ice = bundle.load('cryosphere.ice_thickness').astype(float)
    glacial = bundle.load('geosphere.glacial_denudation').astype(float)
    sea_level = float(bundle.scalar('geosphere.sea_level'))
    lake = (bundle.load('hydrosphere.lake_surface') != 0) if bundle.has('hydrosphere.lake_surface') \
        else np.zeros(land.shape, dtype=bool)
    land_area = weights[land].sum()
    metrics = {'land_cells': int(np.count_nonzero(land))}

    def land_share(mask):
        return float(weights[land & mask].sum() / land_area) if land_area else None

    metrics['glaciated_land_share'] = land_share(glaciation > 0)
    metrics['land_mean_glaciation'] = float(np.average(glaciation[land], weights=weights[land])) if land_area else None
    metrics['ice_covered_land_share'] = land_share(ice > 0)
    metrics['glaciated_land_cells'] = int(np.count_nonzero(land & (glaciation > 0)))
    metrics['ice_covered_land_cells'] = int(np.count_nonzero(land & (ice > 0)))
    metrics['maximum_ice_thickness_m'] = float(ice.max()) if ice.size else None
    # Glaciated sea: sea at the end with glacial denudation. Fjords cut from land and seabed
    # eroded under grounded ice are not told apart: the bundle does not hold the start's sea.
    sea = ~land & ~lake & (glacial > 0)
    depth = sea_level - bundle.load('geosphere.base_elevation').astype(float)[sea]
    metrics['glaciated_sea_cells'] = int(np.count_nonzero(sea))
    if bundle.has('planetbody.radius'):
        radius_km = float(bundle.scalar('planetbody.radius')) / 1000.0
        metrics['glaciated_sea_area_km2'] = float(weights[sea].sum() * 4.0 * np.pi * radius_km ** 2)
    else:
        metrics['glaciated_sea_area_km2'] = None
    metrics['glaciated_sea_min_depth_m'] = float(depth.min()) if depth.size else None
    metrics['glaciated_sea_median_depth_m'] = float(np.median(depth)) if depth.size else None
    metrics['glaciated_sea_cells_deeper_than_100_m'] = int(np.count_nonzero(depth > 100.0))
    metrics['glaciated_cells_below_sea_level'] = int(np.count_nonzero(~land & (glaciation > 0)))
    metrics['maximum_glacial_denudation_m'] = float(glacial.max()) if glacial.size else None
    metrics['glaciated_lake_cells'] = int(np.count_nonzero(lake & (glaciation > 0)))
    metrics['ice_on_lake_cells'] = int(np.count_nonzero(lake & (ice > 0)))
    ela = bundle.load('cryosphere.equilibrium_line_altitude').astype(float)
    computed = bool(np.any(ela != 0))  # zero everywhere: the glacial term was off, no ELA
    ela = ela - sea_level
    latitude = grid.latitudes(land.shape[0])[:, None] * np.ones(land.shape)
    for lower in range(-90, 90, ELA_BAND_DEGREES):
        upper = lower + ELA_BAND_DEGREES
        band = land & (latitude >= lower) & ((latitude < upper) if upper < 90 else (latitude <= upper))
        key = f'ela_above_sea_m_{latitude_name(lower)}_{latitude_name(upper)}'
        metrics[key] = float(np.average(ela[band], weights=weights[band])) if computed and np.any(band) else None
    if land_area:
        glacial_mean = float(np.average(glacial[land], weights=weights[land]))
        metrics['land_mean_glacial_denudation_m'] = glacial_mean
        if bundle.has('geosphere.denudation'):
            denudation_mean = float(np.average(bundle.load('geosphere.denudation').astype(float)[land],
                                               weights=weights[land]))
            metrics['land_mean_denudation_m'] = denudation_mean
            metrics['glacial_denudation_share'] = glacial_mean / denudation_mean if denudation_mean > 0 else None
    return metrics


def compare_glaciation(bundle, spec, baseline=None):
    result = {'reference': 'none: reported numbers of the glacial fields (erosion plan Phase G), compared with the '
                           'baseline run'}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is not None:
            result[key] = {'earth_comparison': {}, 'metrics': measure_glaciation(source)}
    return result


def main():
    """Print the measurements of bundles (checkpoints included) as one JSON document per bundle."""
    import argparse
    import json
    from uwtools.bundle import Bundle
    parser = argparse.ArgumentParser(description='Reported landscape numbers of run bundles: '
                                     'uv run --offline scripts/evaluation/landscape_statistics.py '
                                     '{ocean-depths,drainage-pattern,glaciation} BUNDLE [BUNDLE ...]')
    parser.add_argument('measure', choices=('ocean-depths', 'drainage-pattern', 'glaciation'))
    parser.add_argument('bundles', nargs='+', type=Path)
    args = parser.parse_args()
    measure = {'ocean-depths': measure_ocean_depths, 'drainage-pattern': measure_drainage_pattern,
               'glaciation': measure_glaciation}[args.measure]
    for path in args.bundles:
        print(json.dumps({'bundle': path.as_posix(), 'metrics': measure(Bundle.open(path))}))


if __name__ == '__main__':
    main()
