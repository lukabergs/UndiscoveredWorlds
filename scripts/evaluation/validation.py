"""Validation: comparisons of bundle fields with observed Earth references.

Map comparisons use the per-field evaluation of score.py, so each field's loss is
defined as in global-balanced-v1. Distribution comparisons (hypsometry) use the tectonics Earth
reference measurements; river_mouths compares the routing's largest outlets with the authored
table configs/river-mouths.json. References are never inputs to a simulation.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

from checks import DIRECTION_OFFSETS, has_neighbour, receiver_values, whole_metres
from references import MONTH_DAYS
from score import DEFINITION as PALETTE, VERSION as SCORE_VERSION, WEIGHTS, evaluate, magnitude, weights

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from uwtools import grid, regions  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402

SELECTOR = re.compile(r'^([a-z]+\.[a-z0-9_]+)(?:\[(\d+|@[A-Za-z0-9_.-]+)\])?$')


def parse_selector(selector):
    """'component.name', 'component.name[i]' (i indexes the first axis, e.g. the level)
    or 'component.name[@level]' (an atmosphere level by name, from the atmosphere.levels table)."""
    match = SELECTOR.match(selector)
    if not match:
        raise ValueError(f'Bad field selector {selector!r}')
    index = match.group(2)
    if index is None:
        return match.group(1), None
    return match.group(1), index[1:] if index.startswith('@') else int(index)


def level_index(bundle, name):
    for row in bundle.table_rows('atmosphere.levels'):
        if row['name'] == name:
            return int(row['index'])
    raise KeyError(f'No atmosphere level {name!r} in {bundle.path}')


def available_selector(bundle, selector):
    identifier, index = parse_selector(selector)
    if identifier not in bundle.fields:
        return False
    if isinstance(index, str):
        return 'atmosphere.levels' in bundle.tables and any(row['name'] == index for row in bundle.table_rows('atmosphere.levels'))
    return True


def load_selected(bundle, selector):
    identifier, index = parse_selector(selector)
    values = bundle.load(identifier).astype(float)
    if index is None:
        return values
    return values[level_index(bundle, index) if isinstance(index, str) else index]


def surface_elevation(bundle):
    """Surface elevation relative to sea level (m), as planet::map minus the sea level."""
    base = bundle.load('geosphere.base_elevation').astype(float)
    sea_level = float(bundle.scalar('geosphere.sea_level'))
    volcano = np.abs(bundle.load('geosphere.volcano').astype(float))
    land = base + np.maximum(volcano, bundle.load('geosphere.mountain_height')) + \
        bundle.load('geosphere.extra_elevation') + bundle.load('geosphere.crater_rims') + \
        bundle.load('geosphere.crater_centres')
    sea = base + bundle.load('geosphere.ocean_ridge_height') + volcano
    return np.where(base <= sea_level, sea, land) - sea_level


def season_field(values, spec, season):
    """One season of a seasonal (4, rows, columns) field, or the mean over seasons.

    Monthly (12, rows, columns) fields follow the reference's season rule, so model and
    reference cover the same months.
    """
    if values.shape[0] == 12:
        rule = spec['season']
        if rule == 'month':
            return values[season * 3]
        months = np.arange(12) if rule == 'annual' else np.arange(season * 3, season * 3 + 3) % 12
        return np.average(values[months], axis=0, weights=MONTH_DAYS[months])
    if spec.get('reduce') == 'season-mean':
        return values.mean(axis=0)
    return values[season]


def convert(values, spec, season, cell_weights):
    scale = spec.get('scale', 1)
    if scale == 'per-day-quarter':
        # Seasonal amounts are the mean monthly amount of the season's quarter; the
        # reference quarter is the day-weighted mean rate over the same three months.
        months = np.arange(season * 3, season * 3 + 3) % 12
        values = values * 3 / MONTH_DAYS[months].sum()
    else:
        values = values * float(scale)
    for step in spec.get('transform', []):
        if step == 'remove-area-mean':
            finite = np.isfinite(values)
            values = values - np.sum(np.where(finite, values, 0) * cell_weights) / np.sum(cell_weights * finite)
        else:
            raise ValueError(f'Unknown transform {step!r}')
    return values


def model_field(bundle, spec, season):
    parts = []
    for selector in spec['model']:
        values = season_field(load_selected(bundle, selector), spec, season)
        parts.append(convert(values, spec, season, weights(values.shape)))
    return np.array(parts) if spec['kind'] == 'vector' else parts[0]


def field_mask(bundle, spec, season, shape):
    mask = np.ones(shape, bool)
    rule = spec.get('mask')
    if rule == 'sea':
        # Wet where land covers less than a tenth of the cell, as the climate report's mask.
        mask &= bundle.sea_fraction(*shape) > 0.9
    elif isinstance(rule, dict):
        values = season_field(load_selected(bundle, rule['field']), spec, season)
        mask &= values >= rule.get('min', 0.5)
    elif rule is not None:
        raise ValueError(f"{spec['id']}: unknown mask {rule!r}")
    return mask


def reference_field(spec, references, season, width):
    parts = [references.season(name, spec['season'], season, width) for name in spec['reference']]
    return np.array(parts) if spec['kind'] == 'vector' else parts[0]


def seasons_for(spec):
    return [0] if spec['season'] == 'annual' or spec.get('reduce') == 'season-mean' else list(range(4))


def available(bundle, spec):
    return all(available_selector(bundle, selector) for selector in spec['model'])


class RegionalTable:
    """Area-weighted bias and RMSE of model minus reference per region of configs/regions.json.

    Reported, never part of a score. Project boxes select cells by centre; AR6 regions weight
    cells by their overlap fraction with the rasterised mask on the field's grid.
    """

    def __init__(self):
        self.config = regions.load()
        self.labels = None
        self.grids = {}
        self.fields = {}

    def weights(self, shape):
        if shape not in self.grids:
            if self.labels is None:
                self.labels = regions.ar6_labels(self.config)
            ar6 = regions.ar6_fractions(*shape, config=self.config, labels=self.labels)
            self.grids[shape] = {'project': regions.box_masks(*shape, config=self.config),
                                 'ar6': {item['acronym']: ar6[item['id']] for item in regions.ar6(self.config)}}
        return self.grids[shape]

    def add(self, field, season, values, reference, mask):
        cell_weights = weights(reference.shape)
        entry = {}
        for name, members in self.weights(reference.shape).items():
            rows = {}
            for region, fraction in members.items():
                inside = mask & (fraction > 0)
                if not inside.any():
                    continue
                w = (cell_weights * fraction)[inside]
                row = {'cells': int(inside.sum())}
                for key, value in values.items():
                    difference = (value - reference)[inside]
                    row[key] = {'bias': float(np.average(difference, weights=w)),
                                'rmse': float(np.sqrt(np.average(difference * difference, weights=w)))}
                rows[region] = row
            entry[name] = rows
        self.fields.setdefault(field, {})[str(season)] = entry

    def result(self):
        sets = self.config['sets']
        mask = regions.ar6_mask_path(self.config)
        return {'note': 'Reported, never part of a score.',
                'method': ('Per scalar entry and season: area-weighted mean (bias) and root mean square (rmse) of model '
                           'minus reference over the cells of the entry\'s common mask in each region. Project boxes take '
                           'the cells whose centre lies in the box; AR6 regions weight cells by their overlap fraction with '
                           'the rasterised mask. Regions without valid cells are omitted.'),
                'config': regions.PATH.relative_to(ROOT).as_posix(), 'config_sha256': sha256(regions.PATH),
                'ar6_mask': mask.relative_to(ROOT).as_posix(), 'ar6_mask_sha256': sha256(mask),
                'sets': {'project': {box['id']: box['name'] for box in sets['project']['regions']},
                         'ar6': {item['acronym']: item['name'] for item in sets['ar6']['regions']}},
                'fields': self.fields}


def compare_maps(bundle, spec, references, baseline=None, regional=None):
    """Per-season evaluation rows for the run and, with a common mask, for the baseline.

    With a RegionalTable, scalar entries also add their regional errors on the same mask."""
    # A comparison may borrow another field's scale (e.g. a second 850 hPa wind).
    scale = spec.get('metric_scale', spec['id'])
    if scale not in PALETTE['scales']:
        raise ValueError(f"{spec['id']}: map comparisons need a scale in configs/climate-palettes.json")
    runs = {'run': bundle} if baseline is None else {'run': bundle, 'baseline': baseline}
    seasons = {}
    for season in seasons_for(spec):
        values = {key: model_field(source, spec, season) for key, source in runs.items()}
        shape = magnitude(values['run']).shape
        if baseline is not None and magnitude(values['baseline']).shape != shape:
            raise ValueError(f"{spec['id']}: run and baseline grids differ")
        reference = reference_field(spec, references, season, shape[-1])
        mask = np.isfinite(magnitude(reference))
        for key, source in runs.items():
            mask &= field_mask(source, spec, season, shape) & np.isfinite(magnitude(values[key]))
        seasons[str(season)] = {key: evaluate(value, reference, mask, scale) for key, value in values.items()}
        if regional is not None and spec['kind'] == 'scalar':
            regional.add(spec['id'], season, values, reference, mask)
    return seasons


def compare_hypsometry(bundle, spec, baseline=None):
    from uwtools.tectonics_earth_reference import compare_to_earth, load_earth_reference, measure_earth
    reference = load_earth_reference(ROOT, spec.get('reference_width', 2048))
    result = {'reference': reference['title'], 'reference_cache': reference['cache']}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is None:
            continue
        elevation = surface_elevation(source)
        metrics = measure_earth(elevation, ~source.sea_mask())
        result[key] = {'metrics': {name: metrics[name] for name in
                                   ('land_fraction', 'elevation_m', 'land_elevation_m', 'ocean_depth_m',
                                    'area_above_9000_m')},
                       'earth_comparison': compare_to_earth(metrics, reference)}
    return result


RIVER_FIELDS = ('hydrosphere.flow_direction', 'hydrosphere.discharge', 'geosphere.base_elevation',
                'geosphere.sea_level', 'planetbody.radius')


def river_outlets(bundle):
    """(longitude, latitude, discharge) of the land cells whose flow direction points to a sea cell,
    largest discharge first (ties in row-major order). Sea is the base elevation in whole metres at
    or below sea level; longitude wraps and a polar row has no neighbour beyond the pole."""
    sea = whole_metres(bundle.load('geosphere.base_elevation').astype(float)) <= float(bundle.scalar('geosphere.sea_level'))
    direction = bundle.load('hydrosphere.flow_direction')
    outlet = np.zeros(sea.shape, dtype=bool)
    for code, (dx, dy) in DIRECTION_OFFSETS.items():
        outlet |= ~sea & (direction == code) & has_neighbour(sea.shape, dy) & receiver_values(sea, dx, dy)
    y, x = np.nonzero(outlet)
    discharge = bundle.load('hydrosphere.discharge').astype(float)[y, x]
    order = np.argsort(-discharge, kind='stable')
    rows, columns = sea.shape
    return grid.longitudes(columns)[x[order]], grid.latitudes(rows)[y[order]], discharge[order]


def great_circle_km(lon1, lat1, lon2, lat2, radius_m):
    """Haversine distance (km) on a sphere of the given radius; arguments in degrees, broadcast."""
    lon1, lat1, lon2, lat2 = (np.deg2rad(np.asarray(value, dtype=float)) for value in (lon1, lat1, lon2, lat2))
    h = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * radius_m * np.arcsin(np.sqrt(np.clip(h, 0, 1))) / 1000


def measure_river_mouths(bundle, rivers, top, match_km):
    lon, lat, discharge = river_outlets(bundle)
    count = len(discharge)
    lon, lat, discharge = lon[:top], lat[:top], discharge[:top]
    radius = float(bundle.scalar('planetbody.radius'))
    listed, comparison = [], {}
    for river in rivers:
        distances = great_circle_km(river['lon'], river['lat'], lon, lat, radius)
        nearest = int(np.argmin(distances)) if len(distances) else None
        distance = float(distances[nearest]) if nearest is not None else None
        matched = distance is not None and distance <= match_km
        entry = {'name': river['name'], 'lon': float(river['lon']), 'lat': float(river['lat']),
                 'distance_km': distance, 'matched': bool(matched),
                 'outlet_discharge_m3_s': float(discharge[nearest]) if matched else None,
                 'observed_discharge_m3_s': river['discharge_m3_s']}
        if nearest is not None:
            entry['outlet'] = {'lon': float(lon[nearest]), 'lat': float(lat[nearest])}
        listed.append(entry)
        comparison[f"{river['name']}: distance_km"] = distance
        comparison[f"{river['name']}: matched outlet discharge_m3_s"] = entry['outlet_discharge_m3_s']
    distances = [entry['distance_km'] for entry in listed if entry['distance_km'] is not None]
    metrics = {'matched': sum(entry['matched'] for entry in listed), 'table_rivers': len(listed),
               'median_distance_km': float(np.median(distances)) if distances else None,
               'largest_outlet_discharge_m3_s': float(discharge[0]) if count else None,
               'outlets': count, 'top': int(top), 'match_km': float(match_km)}
    return {'metrics': metrics, 'earth_comparison': comparison, 'rivers': listed,
            'top_outlets': [{'lon': float(a), 'lat': float(b), 'discharge_m3_s': float(q)}
                            for a, b, q in zip(lon, lat, discharge)]}


def compare_river_mouths(bundle, spec, baseline=None):
    """The `top` largest-discharge outlets of the routing against an authored table of Earth's river mouths.

    A table river is matched when one of the top outlets lies within match_km of its mouth;
    the median distance is over the table rivers, each to its nearest top outlet."""
    path = ROOT / spec['table']
    table = json.loads(path.read_text(encoding='utf-8'))
    top, match_km = int(spec.get('top', 20)), float(spec.get('match_km', 500))
    result = {'reference': f"{path.name}: {len(table['rivers'])} authored river mouths, top {top} outlets, "
                           f"match within {match_km:g} km",
              'table': path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix(),
              'table_sha256': sha256(path), 'note': table.get('description', '')}
    for key, source in (('run', bundle), ('baseline', baseline)):
        if source is not None:
            result[key] = measure_river_mouths(source, table['rivers'], top, match_km)
    return result


def field_losses(comparison, key):
    """Season-mean loss of one map comparison, as rank() in the climate report."""
    entries = [season[key]['loss'] for season in comparison.values()]
    return float(np.mean(entries))


def score(comparisons, key):
    """global-balanced-v1 when every weighted field was compared; otherwise the same
    formula over the compared fields, renormalized and reported with its coverage."""
    weighted = dict(WEIGHTS)
    losses = {field: field_losses(comparisons[field], key) for field in weighted if field in comparisons}
    covered = sum(weighted[field] for field in losses)
    missing = sorted(set(weighted) - set(losses))
    total = 100 * sum(weighted[field] * loss for field, loss in losses.items())
    entry = {'version': SCORE_VERSION, 'inputs': 'bundle fields', 'lower_is_better': True, 'field_losses': losses,
             'field_weights': {field: weighted[field] for field in losses}, 'covered_weight': covered,
             'missing_fields': missing}
    if missing:
        entry['partial_score'] = total / covered if covered else None
        entry['note'] = ('Not a global-balanced-v1 score: the bundle lacks ' + ', '.join(missing) +
                         '. The partial score renormalizes the same formula over the compared fields and is '
                         'comparable only between runs with the same coverage.')
    else:
        entry['score'] = total
        entry['note'] = ('Formula and weights of global-balanced-v1, on bundle fields at their native grids. '
                         'Comparable between bundle evaluations; registered climate rankings extract their inputs '
                         'from the run diagnostics on other grids, so the values are not interchangeable.')
    return entry


def run_validation(bundle, specs, references, baseline=None):
    comparisons, skipped, distributions = {}, [], {}
    regional = None
    for spec in specs:
        if spec['kind'] == 'hypsometry':
            distributions[spec['id']] = compare_hypsometry(bundle, spec, baseline)
        elif spec['kind'] in ('drainage_pattern', 'ocean_depths', 'glaciation'):
            import landscape_statistics
            needed = {'drainage_pattern': landscape_statistics.DRAINAGE_FIELDS,
                      'glaciation': landscape_statistics.GLACIATION_FIELDS}.get(spec['kind'], ())
            # A baseline from before the glacial fields (erosion plan Phase G) had no glacial
            # term: the run's numbers are still reported, without a baseline side.
            glacialbaseline = baseline
            if spec['kind'] == 'glaciation' and baseline is not None and \
                    any(field not in baseline.fields for field in needed):
                glacialbaseline = None
            missing = [field for field in needed if field not in bundle.fields or
                       (baseline is not None and spec['kind'] != 'glaciation' and field not in baseline.fields)]
            if missing:
                skipped.append({'id': spec['id'], 'reason': 'model fields absent: ' + ', '.join(missing)})
            elif spec['kind'] == 'drainage_pattern':
                distributions[spec['id']] = landscape_statistics.compare_drainage_pattern(bundle, spec, baseline)
            elif spec['kind'] == 'glaciation':
                distributions[spec['id']] = landscape_statistics.compare_glaciation(bundle, spec, glacialbaseline)
            else:
                distributions[spec['id']] = landscape_statistics.compare_ocean_depths(bundle, spec, baseline)
        elif spec['kind'] in ('upstream_area', 'basins', 'internal_drainage', 'discharge'):
            import hydrology_reference
            reason = hydrology_reference.skip_reason(spec['kind'], bundle, spec, baseline)
            if reason:
                skipped.append({'id': spec['id'], 'reason': reason})
            else:
                distributions[spec['id']] = hydrology_reference.compare(spec['kind'], bundle, spec, baseline)
        elif spec['kind'] == 'river_mouths':
            missing = [field for field in RIVER_FIELDS if field not in bundle.fields or
                       (baseline is not None and field not in baseline.fields)]
            if missing:
                skipped.append({'id': spec['id'], 'reason': 'model fields absent: ' + ', '.join(missing)})
            else:
                distributions[spec['id']] = compare_river_mouths(bundle, spec, baseline)
        elif not available(bundle, spec) or (baseline is not None and not available(baseline, spec)):
            skipped.append({'id': spec['id'], 'reason': 'model fields absent: ' + ', '.join(spec['model'])})
        else:
            if regional is None and spec['kind'] == 'scalar':
                regional = RegionalTable()
            comparisons[spec['id']] = compare_maps(bundle, spec, references, baseline, regional)
    scores = {}
    if comparisons:
        scores['run'] = score(comparisons, 'run')
        if baseline is not None:
            scores['baseline'] = score(comparisons, 'baseline')
    result = {'maps': comparisons, 'distributions': distributions, 'skipped': skipped, 'scores': scores}
    if regional is not None and regional.fields:
        result['regions'] = regional.result()
    return result
