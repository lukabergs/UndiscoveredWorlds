"""Adapt archived surface-feedback experiments without inventing native crust state."""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import palettes, runs
from uwtools.hashing import sha256
from uwtools.paths import ROOT
from tectonics_report_data import (load_bundle, stats, local_relief, coast_distance,
                                   coastal_terrain, terrain_detail, island_sizes,
                                   hypsography_shape)


def climate_field_specs():
    palette = palettes.climate()
    result = {}
    for key, title, units, scale, colours in [
        ('rain', 'Annual rainfall', 'mm/year', 'annual_rain', 'positive'),
        ('temperature', 'Annual air temperature', '°C', 'air_temperature', 'signed'),
        ('wind', 'Annual surface wind speed', 'm/s; mean seasonal speed', 'surface_wind', 'positive'),
        ('snow', 'Snow coverage', 'annual mean fraction', None, 'positive'),
        ('ice', 'Permanent land ice', 'fraction', None, 'positive'),
    ]:
        low, high = palette['scales'][scale] if scale else (0, 1)
        result[key] = dict(title=title, units=units, stops=[
            [low + t * (high-low), rgb] for t, rgb in palette[colours]])
    return result


def surface_metrics(native, parent_fields, z):
    """Equal-cell diagnostics for comparison to the parent torus, not sphere areas."""
    if not np.isfinite(z).all():
        raise ValueError('Nonfinite feedback elevation')
    land = z >= 0
    fields = dict(elevation=z, local_relief=local_relief(z), plate_id=parent_fields['plate_id'])
    curve = dict(fraction_below=np.linspace(0, 1, 201).tolist(),
                 elevation_m=np.quantile(z, np.linspace(0, 1, 201)).tolist())
    return dict(
        checks={'nonfinite_elevation': 0},
        weighting='equal cells for parent comparison; erosion itself uses spherical areas',
        land_fraction=float(land.mean()), elevation_m=stats(z),
        land_elevation_m=stats(z[land]), ocean_depth_m=stats(-z[~land]),
        hypsography=curve, hypsography_shape=hypsography_shape(curve),
        terrain_detail=terrain_detail(fields), coastal_terrain=coastal_terrain(native, fields),
        islands=island_sizes(z), area_above_9000_m=float((z > 9000).mean()),
        area_below_minus_11000_m=float((z < -11000).mean()),
        neighbour_jump_m=stats(np.concatenate([abs(z-np.roll(z, 1, axis)).ravel() for axis in (0, 1)])),
        fraction_below={str(level): float((z < level).mean()) for level in
                        (-6000, -4000, -2000, -1000, -258, -200, -160, 0, 1000, 2000, 4000)},
        shoreline_within_1m=float((abs(z) <= 1).mean()),
        elevation_area_bands=dict(
            shelf_0_200_m=float(((z < 0) & (z >= -200)).mean()),
            upper_ocean_200_2000_m=float(((z < -200) & (z >= -2000)).mean()),
            deep_ocean_below_6000_m=float((z < -6000).mean()),
            lowland_0_1000_m=float(((z >= 0) & (z < 1000)).mean()),
            highland_above_2000_m=float((z >= 2000).mean())))


def feedback_snapshots(record):
    archive = ROOT / f'runs/reports/tectonics/{record["id"]}'
    field_path = archive / 'comparison-fields.npz'
    if sha256(field_path) != record['fields_sha256']:
        raise ValueError('Feedback comparison field hash mismatch')
    with np.load(field_path) as source:
        endpoints = {key: source[key] for key in ('initial', 'final')}
    parent = runs.tectonics(record['parent'])
    final = next(s for s in parent['snapshots'] if s['label'] == 'final')
    native, parent_fields, _ = load_bundle(final['path'], record['parent_final_hashes'])
    if not np.array_equal(endpoints['initial'], parent_fields['elevation']):
        raise ValueError('Feedback initial terrain does not match archived parent')
    for index in range(record['cycles'] + 1):
        label = 'initial' if index == 0 else 'final' if index == record['cycles'] else f'surface-{index}'
        climate_path = archive / f'climate-{index}/annual-fields.npz'
        with np.load(climate_path) as source:
            fields = {key: source[key] for key in climate_field_specs()}
        if label in endpoints:
            z = endpoints[label]
        else:
            path = archive / f'erosion-{index-1}/output/elevation.f64'
            z = np.fromfile(path, dtype='<f8').reshape(endpoints['initial'].shape)
        if any(value.shape != z.shape or not np.isfinite(value).all() for value in fields.values()):
            raise ValueError('Invalid climate grid')
        fields.update(elevation=z, local_relief=local_relief(z),
                      coast_distance_km=coast_distance(z >= 0, native.get('run_scenario', {}).get('world_width_km', 40000) / z.shape[1]),
                      _boundary_overlay=parent_fields['boundary_type'])
        metrics = surface_metrics(native, parent_fields, z)
        metadata = dict(label=label,
                        time_myr=record['tectonic_endpoint_time_myr'] + index * record['surface_interval_myr'],
                        update_count=None, height=z.shape[0], width=z.shape[1],
                        run_scenario={'surface_model': 'climate-feedback; frozen plates'},
                        phase=f'{index * record["surface_interval_myr"]} Myr surface erosion after run {record["parent"]}',
                        boundary_note=f'Frozen run {record["parent"]} boundary overlay; no plate motion in this phase.')
        hashes = dict(comparison_npz=record['fields_sha256'], climate_npz=sha256(climate_path))
        yield metadata, fields, hashes, metrics


def feedback_controls(record):
    """Return measured control surfaces and differences for a separate visual panel."""
    archive = ROOT / f'runs/reports/tectonics/{record["id"]}'
    paths = {'comparison': archive / 'comparison-fields.npz',
             'matched': archive / 'matched-comparison-fields.npz'}
    if sha256(paths['comparison']) != record['fields_sha256']:
        raise ValueError('Feedback comparison hash mismatch')
    if sha256(paths['matched']) != record['matched_control']['fields_sha256']:
        raise ValueError('Matched-control hash mismatch')
    with np.load(paths['comparison']) as source:
        data = dict(source)
    with np.load(paths['matched']) as source:
        data.update(matched_control=source['matched_control'], spatial_effect=source['spatial_effect'])
    review = json.loads((archive / 'review.json').read_text())
    return data, review
