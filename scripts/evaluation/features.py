"""Feature diagnostics and the footprint of a run against its parent. Reported, never part of a score.

Features are read from the fields of the profile's validation entries, with the entry's
own reduction, for the run, the reference and the parent:

- ITCZ latitude: latitude of the largest zonal-mean precipitation within 30 degrees of the
  equator, per season (entry "rain").
- 500 hPa jet: latitude and speed of the largest zonal-mean eastward wind between 15 and 75
  degrees in each hemisphere, per season (entry "500_wind").
- Monsoon reversal: share of the area between 35 S and 35 N where the surface wind turns by
  120 degrees or more between January and July and blows at 1 m/s or more in both (the
  direction criterion of Ramage's monsoon definition; entry "surface_wind").
- Koeppen groups: share of the land area in each major group A to E (model only: no
  processed Koeppen reference exists, so there is no agreement measure).

The footprint says where a run differs from its parent: per validated field and season,
the share of the common area whose values differ by more than one hundredth of the field's
display range, and the largest difference.
"""
import numpy as np

from uwtools import grid, palettes
from uwtools.bundle import SEASONS
from validation import available, field_mask, magnitude, model_field, reference_field, seasons_for

KOPPEN_GROUPS = {'A': range(1, 5), 'B': range(5, 9), 'C': range(9, 18), 'D': range(18, 30), 'E': range(30, 32)}
METHOD = __doc__.split('\n\n', 1)[1]


def zonal_mean(values):
    """Mean along each row over its finite cells (NaN for a row without any)."""
    finite = np.isfinite(values)
    count = finite.sum(axis=1)
    return np.where(count > 0, np.where(finite, values, 0).sum(axis=1) / np.maximum(count, 1), np.nan)


def itcz_latitude(rain):
    latitudes = grid.latitudes(rain.shape[0])
    profile = np.where(np.abs(latitudes) <= 30, zonal_mean(rain), np.nan)
    return None if not np.isfinite(profile).any() else float(latitudes[int(np.nanargmax(profile))])


def jet(east, hemisphere):
    """(latitude, speed) of the largest zonal-mean eastward wind of a hemisphere (+1 north, -1 south)."""
    latitudes = grid.latitudes(east.shape[0])
    band = (latitudes * hemisphere >= 15) & (latitudes * hemisphere <= 75)
    profile = np.where(band, zonal_mean(east), np.nan)
    if not np.isfinite(profile).any():
        return None, None
    row = int(np.nanargmax(profile))
    return float(latitudes[row]), float(profile[row])


def monsoon_reversal(january, july):
    """Percent of the area between 35 S and 35 N whose wind turns by at least 120 degrees between the two."""
    speeds = np.hypot(*january), np.hypot(*july)
    cosine = (january[0] * july[0] + january[1] * july[1]) / np.maximum(speeds[0] * speeds[1], 1e-30)
    latitudes = grid.latitudes(january.shape[-2])[:, None]
    tropics = (np.abs(latitudes) <= 35) & np.isfinite(cosine)
    weights = grid.area_weights(*january.shape[-2:])
    reversed_ = tropics & (cosine <= np.cos(np.deg2rad(120))) & (speeds[0] >= 1) & (speeds[1] >= 1)
    total = float(weights[tropics].sum())
    return None if total == 0 else 100 * float(weights[reversed_].sum()) / total


def koppen_groups(bundle):
    if 'atmosphere.climate_class' not in bundle.fields:
        return None
    classes = bundle.load('atmosphere.climate_class')
    weights = grid.area_weights(*classes.shape)
    land = classes > 0
    total = float(weights[land].sum())
    if total == 0:
        return None
    return {group: 100 * float(weights[np.isin(classes, list(codes))].sum()) / total for group, codes in KOPPEN_GROUPS.items()}


def features(bundle, specs, references, baseline=None):
    """Rows {feature, season, unit, run, reference, baseline} of the features the bundle's fields allow."""
    entries = {spec['id']: spec for spec in specs if spec.get('kind') in ('scalar', 'vector')}
    sources = {'run': bundle} if baseline is None else {'run': bundle, 'baseline': baseline}
    rows = []

    def usable(name):
        return name in entries and all(available(source, entries[name]) for source in sources.values())

    def add(feature, season, unit, values):
        if any(value is not None for value in values.values()):
            rows.append({'feature': feature, 'season': season, 'unit': unit, **{key: values.get(key) for key in ('run', 'reference', 'baseline')}})

    def fields(name, season):
        spec = entries[name]
        result = {key: model_field(source, spec, season) for key, source in sources.items()}
        result['reference'] = reference_field(spec, references, season, magnitude(result['run']).shape[-1])
        return result

    if usable('rain'):
        for season in seasons_for(entries['rain']):
            label = SEASONS[season] if entries['rain']['season'] != 'quarter' else ('Jan–Mar', 'Apr–Jun', 'Jul–Sep', 'Oct–Dec')[season]
            add('ITCZ latitude', label, 'degrees north', {key: itcz_latitude(value) for key, value in fields('rain', season).items()})
    if usable('500_wind'):
        for season in seasons_for(entries['500_wind']):
            values = fields('500_wind', season)
            for hemisphere, name in ((1, 'northern'), (-1, 'southern')):
                found = {key: jet(value[0], hemisphere) for key, value in values.items()}
                add(f'500 hPa jet latitude, {name}', SEASONS[season], 'degrees north', {key: value[0] for key, value in found.items()})
                add(f'500 hPa jet speed, {name}', SEASONS[season], 'm/s', {key: value[1] for key, value in found.items()})
    if usable('surface_wind') and len(seasons_for(entries['surface_wind'])) == 4:
        january, july = fields('surface_wind', 0), fields('surface_wind', 2)
        add('Monsoon reversal area', 'January to July', '% of 35S–35N', {key: monsoon_reversal(january[key], july[key]) for key in january})
    groups = {key: koppen_groups(source) for key, source in sources.items()}
    if groups['run'] is not None:
        for group in KOPPEN_GROUPS:
            add(f'Köppen group {group}', 'annual', '% of land', {key: value[group] for key, value in groups.items() if value is not None})
    return {'note': 'Reported, never part of a score.', 'method': METHOD, 'rows': rows}


def footprint(bundle, baseline, specs):
    """Per validated field and season: how much of the common area differs between run and parent."""
    scales = palettes.climate()['scales']
    rows = []
    for spec in specs:
        if spec.get('kind') not in ('scalar', 'vector') or not (available(bundle, spec) and available(baseline, spec)):
            continue
        low, high = scales[spec.get('metric_scale', spec['id'])]
        threshold = (high - low) / 100
        for season in seasons_for(spec):
            run, parent = magnitude(model_field(bundle, spec, season)), magnitude(model_field(baseline, spec, season))
            if run.shape != parent.shape:
                continue
            mask = np.isfinite(run) & np.isfinite(parent) & field_mask(bundle, spec, season, run.shape) & field_mask(baseline, spec, season, run.shape)
            if not mask.any():
                continue
            weights = grid.area_weights(*run.shape)
            difference = np.abs(run - parent)
            changed = mask & (difference > threshold)
            row = int(np.argmax(np.where(mask, difference, -1)) // run.shape[1])
            rows.append({'field': spec['id'], 'season': 'annual' if len(seasons_for(spec)) == 1 else SEASONS[season],
                         'threshold': threshold, 'changed_area_percent': 100 * float(weights[changed].sum()) / float(weights[mask].sum()),
                         'largest_difference': float(difference[mask].max()),
                         'largest_difference_latitude': float(grid.latitudes(run.shape[0])[row])})
    return {'note': 'Where the run differs from its parent; vector fields by speed.',
            'method': 'Share of the common area whose values differ by more than one hundredth of the display range of the field '
                      '(configs/climate-palettes.json), the largest difference and its latitude.',
            'rows': rows}
