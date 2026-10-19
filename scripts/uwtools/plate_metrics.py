"""Plate measures of a world bundle (terrain roadmap T8.0): plate count, the true spherical area of each plate,
the share of the largest, the count above a share (5 %), and area-weighted RMS speeds.

Bundles v1 and v2 (uwtools.bundle). Areas: on a lat-lon bundle each cell of the cell-centred W x W/2 grid has its
exact spherical area R^2 (2 pi / W)(sin north - sin south) (uwtools.grid.area_weights), R the bundle's
planetbody.radius; on a cube bundle the cell table's area_m2 and latitude_deg (true areas, written by C++).
Plates are the native geosphere.plate_id field ((rows, columns) or (6, N, N)).

Speeds: the bundle has no plate-velocity grid; the per-plate record is the tectonic snapshot's plates list
(tables/geosphere/tectonic_snapshot.json). Contract schema 8 (the cube engine) gives each plate's rotation vector
rotation_rad_s (rad/s, the planet's axes: +X at 0 N 0 E, +Y at 0 N 90 E, +Z north), its centroid
centroid_lon_lat_deg and the centroid's velocity_east_north_cm_yr: the speed of every cell is then |omega x R p|
at its own centre p (rotation included), and rms_cm_yr_sphere is its area-weighted RMS; per plate the centroid
speed and the speed's area RMS are reported. Records before schema 8 (the torus engine) give linear_velocity in
engine cells per update, x east, y south, and angular_velocity, read as follows. An engine cell is uwtools.grid.cell_km: run_scenario.world_width_km / width km while the
scenario records it (the torus engine), else 2 pi R / width; an update is delta_time_myr, so the
engine's speed is |v| * that / dt in km/Myr (= mm/yr; motion_diagnostics.max_speed_cm_yr is the largest).
Two RMS speeds, weighted by true cell area: 'engine' in the engine's flat metric (a column and a row the same
length), 'sphere' the same cell displacement on the sphere the map is drawn on (east component times
cos(latitude): a column is 2 pi R cos(lat) / W, a row pi R / H; W and H the plate grid's on lat-lon, the
snapshot's engine width and height on a cube, lat each cell's); rms_cm_yr_sphere is the baseline the sphere
engine's area-weighted RMS speed is compared with. per_plate_mean_cm_yr_engine is the snapshot's
motion_diagnostics.mean_speed_cm_yr: the unweighted mean over plates. Rotation is not added: angular_velocity is
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
from uwtools.bundle import Bundle  # noqa: E402


def area_weighted_rms(speed, area):
    return float(np.sqrt(np.sum(area * speed * speed) / np.sum(area)))


SECONDS_PER_MYR = 3.15576e13  # tectonics::dynamics::kSecondsPerMyr
CM_YR_PER_M_S = 100.0 * SECONDS_PER_MYR / 1.0e6


def cell_unit_vectors(bundle, shape):
    """Each cell centre of the bundle's native grid as a unit vector, shape + (3,)."""
    if bundle.grid_kind == 'cube':
        cells = bundle.cells()
        lat = np.radians(np.asarray(cells['latitude_deg'], dtype=np.float64)).reshape(shape)
        lon = np.radians(np.asarray(cells['longitude_deg'], dtype=np.float64)).reshape(shape)
    else:
        rows, columns = shape
        lat, lon = np.meshgrid(np.radians(grid.latitudes(rows)), np.radians(grid.longitudes(columns)), indexing='ij')
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def rotation_speeds(bundle, snapshot, plate, area, radius_km, plates, ids):
    """Contract schema 8: the speed |omega x R p| of every cell from its plate's rotation, in cm/yr."""
    omega = np.zeros(plate.shape + (3,))
    known = {i: p for i, p in enumerate(snapshot['plates'])}
    for i, p in known.items():
        omega[plate == i] = p['rotation_rad_s']
    p = cell_unit_vectors(bundle, plate.shape)
    speed = np.linalg.norm(np.cross(omega, radius_km * 1000.0 * p), axis=-1) * CM_YR_PER_M_S
    for entry in plates:
        v = known.get(entry['plate_id'])
        if v:
            mask = plate == entry['plate_id']
            entry['centroid_lon_lat_deg'] = v['centroid_lon_lat_deg']
            entry['speed_cm_yr_centroid'] = round(float(np.hypot(*v['velocity_east_north_cm_yr'])), 3)
            entry['speed_cm_yr_rms'] = round(area_weighted_rms(speed[mask], area[mask]), 3)
    result = {
        'source': 'tectonic_snapshot.json plates[].rotation_rad_s (contract schema 8): |omega x R p| per cell',
        'time_myr': snapshot['time_myr'], 'delta_time_myr': snapshot['delta_time_myr'],
        'rms_cm_yr_sphere': round(area_weighted_rms(speed, area), 3),
        'mean_cm_yr_sphere': round(float(np.sum(area * speed) / np.sum(area)), 3),
        'max_cm_yr_sphere': round(float(speed.max()), 3),
        'diagnostics_max_speed_cm_yr': snapshot.get('motion_diagnostics', {}).get('max_speed_cm_yr'),
        'per_plate_mean_cm_yr_engine': snapshot.get('motion_diagnostics', {}).get('mean_speed_cm_yr'),
        'rotation_included': True,
        'plates_without_velocity': [int(i) for i in ids if int(i) not in known]}
    if 'geosphere.deformation_velocity_x' in bundle.fields:
        # (east, north) in world cells per update: on lat-lon the equatorial cell, on a cube the mean side.
        cells = plate.size
        cell_km = (2.0 * math.pi * radius_km / plate.shape[-1] if bundle.grid_kind != 'cube'
                   else math.sqrt(4.0 * math.pi * radius_km ** 2 / cells))
        dx = bundle.load_native('geosphere.deformation_velocity_x').astype(np.float64)
        dy = bundle.load_native('geosphere.deformation_velocity_y').astype(np.float64)
        result['deformation_rms_cm_yr'] = round(
            area_weighted_rms(np.hypot(dx, dy) * cell_km / snapshot['delta_time_myr'] / 10.0, area), 3)
    return result


def measure(directory, threshold):
    bundle = Bundle.open(directory)
    entries = bundle.fields
    plate = bundle.load_native('geosphere.plate_id')
    radius_km = float(bundle.load_native('planetbody.radius')) / 1000.0
    cube = bundle.grid_kind == 'cube'
    if cube:
        cells = bundle.cells()
        area = np.asarray(cells['area_m2'], dtype=np.float64).reshape(plate.shape) / 1.0e6
        lat = np.radians(np.asarray(cells['latitude_deg'], dtype=np.float64)).reshape(plate.shape)
        grid_record = {'kind': 'cube', 'face_cells': bundle.face_cells, 'export': [bundle.columns, bundle.rows]}
    else:
        rows, columns = plate.shape
        area = np.asarray(grid.area_weights(rows, columns)) * (2.0 * math.pi / columns) * radius_km ** 2
        lat = np.radians(grid.latitudes(rows))[:, None]
        grid_record = [columns, rows]
    total = float(area.sum())
    ids = np.unique(plate)
    plates = []
    for i in ids:
        mask = plate == i
        plates.append({'plate_id': int(i), 'cells': int(mask.sum()), 'area_km2': round(float(area[mask].sum()), 0),
                       'share': round(float(area[mask].sum()) / total, 4),
                       'cell_share': round(float(mask.mean()), 4)})
    shares = sorted((p['share'] for p in plates), reverse=True)
    record = {'bundle': directory.as_posix(), 'grid': grid_record, 'radius_km': radius_km,
              'surface_km2': round(total, 0), 'plate_count': len(plates), 'largest_share': shares[0],
              'threshold': threshold, 'plates_above_threshold': sum(s > threshold for s in shares),
              'shares_descending': shares}
    record['fields'] = {
        'plate': sorted(k for k in entries if 'plate' in k),
        'velocity': sorted(k for k in entries if 'velocity' in k and k.startswith('geosphere.')),
        'note': 'no plate-velocity grid in the bundle; per-plate velocities from the tectonic snapshot table'}
    snapshot_path = directory / 'tables/geosphere/tectonic_snapshot.json'
    snapshot = json.loads(snapshot_path.read_text()) if snapshot_path.is_file() else None
    if snapshot is not None and snapshot['plates'] and 'rotation_rad_s' in snapshot['plates'][0]:
        record['speed'] = rotation_speeds(bundle, snapshot, plate, area, radius_km, plates, ids)
    elif snapshot is not None:
        scenario = snapshot['run_scenario']
        km_per_cell = grid.cell_km(scenario, snapshot['width'], radius_km, bundle.grid_kind)
        mm_yr = km_per_cell / snapshot['delta_time_myr']  # one cell per update, in km/Myr = mm/yr
        if cube:
            columns, rows = snapshot['width'], snapshot['height']  # the engine's cells
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
            'per_plate_mean_cm_yr_engine': snapshot.get('motion_diagnostics', {}).get('mean_speed_cm_yr'),
            'max_abs_angular_velocity': angular,
            'rotation_included': False,
            'rotation_flag': 'non-zero angular velocity: the RMS lacks the rotational part' if angular else None,
            'plates_without_velocity': missing}
        if 'geosphere.deformation_velocity_x' in entries:
            dx = bundle.load_native('geosphere.deformation_velocity_x').astype(np.float64)
            dy = bundle.load_native('geosphere.deformation_velocity_y').astype(np.float64)
            record['speed']['deformation_rms_cm_yr_engine'] = round(
                area_weighted_rms(np.hypot(dx, dy) * mm_yr / 10.0, area), 3)
    record['plates'] = plates
    return record


def earth_summary(path):
    reference = json.loads(Path(path).read_text())
    areas = sum(reference['areas_sr_pb2002'].values())
    assert abs(areas - 4.0 * math.pi) < 1e-4, f'{path}: plate areas sum to {areas:.5f} sr, not 4 pi'
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
