"""Plate measures of a world bundle (terrain roadmap T8.0): plate count, the true spherical area of each plate,
the share of the largest, the count above a share (5 %), and area-weighted RMS speeds.

Areas: each cell of the cell-centred W x W/2 grid has its exact spherical area R^2 (2 pi / W)(sin north - sin south)
(uwtools.grid.area_weights), R the bundle's planetbody.radius. Plates are the geosphere.plate_id grid.

Speeds: the bundle has no plate-velocity grid; the per-plate record is the tectonic snapshot's plates list
(tables/geosphere/tectonic_snapshot.json: linear_velocity in world cells per update, x east, y south;
angular_velocity). A cell is run_scenario.world_width_km / width km and an update delta_time_myr, so the
engine's speed is |v| * that / dt in km/Myr (= mm/yr; motion_diagnostics.max_speed_cm_yr is the largest).
Two RMS speeds, weighted by true cell area: 'engine' in the engine's flat metric (a column and a row the same
length), 'sphere' the same cell displacement on the sphere the map is drawn on (east component times
cos(latitude): a column is 2 pi R cos(lat) / W, a row pi R / H). Rotation is not added: angular_velocity is
reported, and a bundle with any non-zero value is flagged. geosphere.deformation_velocity_x/_y (the deforming
part, same units) is reported separately when present.

--earth FILE adds the summary of an Earth reference (configs/evaluation/plate-earth-reference.json).

  uv run python scripts/uwtools/plate_metrics.py BUNDLE [--threshold 0.05] [--earth FILE] [--json OUT]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import grid  # noqa: E402


def field(directory, entries, field_id):
    return np.load(directory / entries[field_id]['file'])


def area_weighted_rms(speed, area):
    return float(np.sqrt(np.sum(area * speed * speed) / np.sum(area)))


def measure(directory, threshold):
    manifest = json.loads((directory / 'manifest.json').read_text())
    entries = {f['id']: f for f in manifest['fields']}
    plate = field(directory, entries, 'geosphere.plate_id')
    rows, columns = plate.shape
    radius_km = float(field(directory, entries, 'planetbody.radius')) / 1000.0
    area = np.asarray(grid.area_weights(rows, columns)) * (2.0 * math.pi / columns) * radius_km ** 2
    total = float(area.sum())
    ids = np.unique(plate)
    plates = []
    for i in ids:
        mask = plate == i
        plates.append({'plate_id': int(i), 'cells': int(mask.sum()), 'area_km2': round(float(area[mask].sum()), 0),
                       'share': round(float(area[mask].sum()) / total, 4),
                       'cell_share': round(float(mask.mean()), 4)})
    shares = sorted((p['share'] for p in plates), reverse=True)
    record = {'bundle': directory.as_posix(), 'grid': [columns, rows], 'radius_km': radius_km,
              'surface_km2': round(total, 0), 'plate_count': len(plates), 'largest_share': shares[0],
              'threshold': threshold, 'plates_above_threshold': sum(s > threshold for s in shares),
              'shares_descending': shares}
    lat = np.radians(grid.latitudes(rows))[:, None]
    record['fields'] = {
        'plate': sorted(k for k in entries if 'plate' in k),
        'velocity': sorted(k for k in entries if 'velocity' in k and k.startswith('geosphere.')),
        'note': 'no plate-velocity grid in the bundle; per-plate velocities from the tectonic snapshot table'}
    snapshot_path = directory / 'tables/geosphere/tectonic_snapshot.json'
    if snapshot_path.is_file():
        snapshot = json.loads(snapshot_path.read_text())
        scenario = snapshot['run_scenario']
        km_per_cell = scenario['world_width_km'] / snapshot['width']
        mm_yr = km_per_cell / snapshot['delta_time_myr']  # one cell per update, in km/Myr = mm/yr
        velocity = {p['plate_id']: p for p in snapshot['plates']}
        missing = [int(i) for i in ids if int(i) not in velocity]
        vx = np.zeros(plate.shape)
        vy = np.zeros(plate.shape)
        for i, p in velocity.items():
            mask = plate == i
            vx[mask], vy[mask] = p['linear_velocity']
        engine = np.hypot(vx, vy) * mm_yr / 10.0
        sphere = np.hypot(vx * (2 * math.pi * radius_km * np.cos(lat) / columns),
                          vy * math.pi * radius_km / rows) / snapshot['delta_time_myr'] / 10.0
        for p in plates:
            v = velocity.get(p['plate_id'])
            if v:
                p['speed_cm_yr_engine'] = round(float(np.hypot(*v['linear_velocity']) * mm_yr / 10.0), 3)
                p['angular_velocity'] = v['angular_velocity']
        angular = max(abs(p['angular_velocity']) for p in snapshot['plates'])
        record['speed'] = {
            'source': 'tectonic_snapshot.json plates[].linear_velocity (world cells per update)',
            'km_per_cell': km_per_cell, 'delta_time_myr': snapshot['delta_time_myr'],
            'time_myr': snapshot['time_myr'],
            'rms_cm_yr_engine': round(area_weighted_rms(engine, area), 3),
            'rms_cm_yr_sphere': round(area_weighted_rms(sphere, area), 3),
            'mean_cm_yr_engine': round(float(np.sum(area * engine) / total), 3),
            'max_cm_yr_engine': round(float(engine.max()), 3),
            'diagnostics_max_speed_cm_yr': snapshot.get('motion_diagnostics', {}).get('max_speed_cm_yr'),
            'diagnostics_mean_speed_cm_yr': snapshot.get('motion_diagnostics', {}).get('mean_speed_cm_yr'),
            'max_abs_angular_velocity': angular,
            'rotation_included': False,
            'rotation_flag': 'non-zero angular velocity: the RMS lacks the rotational part' if angular else None,
            'plates_without_velocity': missing}
        if 'geosphere.deformation_velocity_x' in entries:
            dx = field(directory, entries, 'geosphere.deformation_velocity_x').astype(np.float64)
            dy = field(directory, entries, 'geosphere.deformation_velocity_y').astype(np.float64)
            record['speed']['deformation_rms_cm_yr_engine'] = round(
                area_weighted_rms(np.hypot(dx, dy) * mm_yr / 10.0, area), 3)
    record['plates'] = plates
    return record


def earth_summary(path):
    reference = json.loads(Path(path).read_text())
    return {'file': Path(path).as_posix(), 'status': reference.get('status'), **reference.get('summary', {})}


def main():
    args = sys.argv[1:]
    def option(name, default=None):
        if name not in args:
            return default
        i = args.index(name)
        value = args[i + 1]
        del args[i:i + 2]
        return value
    out = option('--json')
    threshold = float(option('--threshold', '0.05'))
    earth = option('--earth')
    if len(args) != 1:
        raise SystemExit(__doc__)
    record = measure(Path(args[0]), threshold)
    if earth:
        record['earth'] = earth_summary(earth)
    s = record.get('speed', {})
    print(f"{record['bundle']}: {record['plate_count']} plates, largest {record['largest_share']:.3f}, "
          f"{record['plates_above_threshold']} above {threshold:.0%}; RMS {s.get('rms_cm_yr_engine')} cm/yr engine, "
          f"{s.get('rms_cm_yr_sphere')} cm/yr sphere; velocity fields {record['fields']['velocity']}")
    if out:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text(json.dumps(record, indent=1))


if __name__ == '__main__':
    main()
