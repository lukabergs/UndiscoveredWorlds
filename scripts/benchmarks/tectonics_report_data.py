"""Read native tectonics fields and measure them independently of preview colours."""
import json
from collections import deque
from pathlib import Path

import numpy as np

from tectonics_experiment import digest

DTYPES = {
    "heightmap": "<f4", "plate_id": "<u4", "crust_age_steps": "<u4",
    "crust_age_myr": "<f4", "crust_thickness": "<f4", "crust_class": "u1",
    "uplift_tendency": "<f4", "subsidence_tendency": "<f4", "accumulated_strain": "<f4",
    "boundary_type": "u1", "boundary_distance": "<u2", "boundary_segment_id": "<u4",
    "nearest_boundary_id": "<u4", "deforming_region_id": "<u4", "deforming_region_type": "u1",
    "deformation_rate": "<f4", "deformation_velocity_x": "<f4", "deformation_velocity_y": "<f4",
    "convergence_score": "u1", "divergence_score": "u1", "shear_score": "u1", "geologic_regime": "u1",
}
OPTIONAL_DTYPES = {"rift_damage": "<f4", "continental_crust_m": "<f4", "oceanic_crust_m": "<f4", "sediment_m": "<f4"}


def elevation_m(height, sea_level):
    """Native physical conversion, relative to sea datum; NO upper image-format clamp.

    Matches lithosphere::internalToPhysicalMeters, before uint16 rounding.
    Sub-oceanic-base values map to the native physical floor and are counted separately.
    """
    height = np.asarray(height, dtype=np.float64)
    base = float(np.float32(.1))
    ceiling = float(np.nextafter(np.float32(1), np.float32(.1)))
    ocean = np.maximum(0, (height - base) / (ceiling - base)) * max(0, sea_level - 1)
    absolute = np.where(height < 1, ocean, sea_level + (height - 1) * (65535 - sea_level))
    return absolute - sea_level


def load_bundle(path, expected_hashes=None):
    path = Path(path)
    manifest = json.loads((path / "manifest.json").read_text())
    if manifest.get("schema") != "tectonic-snapshot-bundle/v1" or manifest.get("contract_schema_version") != 6:
        raise ValueError("Unsupported native snapshot schema")
    shape = (int(manifest["height"]), int(manifest["width"]))
    if min(shape) <= 0:
        raise ValueError("Empty native grid")
    fields, hashes = {}, {}
    for name, filename in manifest["files"].items():
        source = path / filename
        if not source.resolve().is_relative_to(path.resolve()):
            raise ValueError("Field path escapes bundle")
        hashes[name] = digest(source)
        if expected_hashes is not None and expected_hashes.get(name) != hashes[name]:
            raise ValueError(f"Archived field hash mismatch: {name}")
        if name not in DTYPES and name not in OPTIONAL_DTYPES:
            continue
        dtype = np.dtype((DTYPES | OPTIONAL_DTYPES)[name])
        if source.stat().st_size != int(np.prod(shape)) * dtype.itemsize:
            raise ValueError(f"Wrong field size: {name}")
        fields[name] = np.fromfile(source, dtype=dtype).reshape(shape)
        if not np.isfinite(fields[name]).all():
            raise ValueError(f"Nonfinite field: {name}")
    if set(DTYPES) - fields.keys():
        raise ValueError("Incomplete native field stack")
    if manifest.get("column_budget", {}).get("enabled") and {"continental_crust_m", "oceanic_crust_m", "sediment_m"} - fields.keys():
        raise ValueError("Incomplete physical crust fields")
    fields["elevation"] = elevation_m(fields["heightmap"], manifest["sea_level_m"])
    fields["local_relief"] = local_relief(fields["elevation"])
    cell_km = manifest.get("run_scenario", {}).get("world_width_km", 40000) / shape[1]
    fields["coast_distance_km"] = coast_distance(fields["elevation"] >= 0, cell_km)
    identities = manifest.get("topology_diagnostics", {}).get("stable_ids")
    if identities is not None:
        if (len(identities) != len(manifest["plates"]) or len(set(identities)) != len(identities)
                or any(type(i) is not int or i < 0 for i in identities)
                or fields["plate_id"].max() >= len(identities)):
            raise ValueError("Invalid persistent plate identities")
        fields["plate_identity"] = np.asarray(identities, dtype=np.uint32)[fields["plate_id"]]
    else:
        fields["plate_identity"] = fields["plate_id"].copy()  # Historical runs only have indices.
    highlands = fields["elevation"] >= 3000
    fields["mountain_zones"] = np.where(highlands,
        np.where(near_convergence(fields["boundary_type"]), 2, 1), 0).astype(np.uint8)
    return manifest, fields, hashes


def stats(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return None
    q = np.quantile(values, [0, .01, .5, .95, .99, 1])
    return dict(zip(("min", "p01", "median", "p95", "p99", "max"), q.tolist()), mean=float(values.mean()))


def local_relief(elevation):
    """Signed departure from a 5 x 5 cell mean, for diagnosis only."""
    return elevation - sum(np.roll(elevation, (dy, dx), axis=(0, 1))
                           for dy in range(-2, 3) for dx in range(-2, 3)) / 25


def coast_distance(land, cell_km):
    """Four-neighbour inland distance on the torus, from the land/water face.

    -1 denotes water or undefined distance on an entirely dry world.
    This is a grid distance, not a subcell coastline reconstruction.
    """
    land = np.asarray(land, dtype=bool)
    h, w = land.shape
    coastal = land & (np.roll(~land, 1, 0) | np.roll(~land, -1, 0)
                      | np.roll(~land, 1, 1) | np.roll(~land, -1, 1))
    steps = np.full(land.shape, -1, dtype=np.int32)
    steps[coastal] = 0
    queue = deque(zip(*np.nonzero(coastal)))
    while queue:
        y, x = queue.popleft()
        for ny, nx in (((y-1) % h, x), ((y+1) % h, x), (y, (x-1) % w), (y, (x+1) % w)):
            if land[ny, nx] and steps[ny, nx] < 0:
                steps[ny, nx] = steps[y, x] + 1
                queue.append((ny, nx))
    return np.where(steps >= 0, (steps + .5) * cell_km, -1.0)


def coastal_terrain(manifest, fields):
    z = fields['elevation']
    cell_km = manifest.get('run_scenario', {}).get('world_width_km', 40000) / z.shape[1]
    distance = fields.get('coast_distance_km')
    if distance is None:
        distance = coast_distance(z >= 0, cell_km)
    coastal = (distance >= 0) & (distance < cell_km)
    high = z >= 3000
    jumps = []
    for axis in (0, 1):
        other = np.roll(z, 1, axis)
        jumps.extend(abs(z - other)[(z >= 0) != (other >= 0)].tolist())
    profile = []
    for low in range(0, 1600, 200):
        sample = z[(distance >= low) & (distance < low + 200)]
        summary = stats(sample)
        profile.append({'distance_km': low + 100, 'cells': int(sample.size),
                        'mean_m': summary['mean'] if summary else None,
                        'p95_m': summary['p95'] if summary else None})
    active = coastal & near_convergence(fields['boundary_type']) if 'boundary_type' in fields else np.zeros(z.shape, dtype=bool)
    return {'coastal_land_m': stats(z[coastal]), 'coastal_face_jump_m': stats(jumps),
            'active_coastal_land_m': stats(z[active]), 'active_coastal_cells': int(active.sum()),
            'coastal_highland_cells': int((coastal & high).sum()),
            'highland_fraction_on_coast': float((coastal & high).sum() / max(1, high.sum())),
            'highland_mean_inland_km': float(distance[high & (distance >= 0)].mean())
                if (high & (distance >= 0)).any() else None,
            'profile': profile, 'cell_width_km': cell_km}


def plate_events(manifest):
    """Validate ancestry independently of raster connectivity; disappearance is allowed."""
    topology = manifest.get("topology_diagnostics", {})
    events = topology.get("events", [])
    if "stable_ids" not in topology:
        return {"available": False, "splits": 0, "fusions": 0, "invalid": 0, "events": []}
    initial = manifest.get("run_scenario", {}).get("plate_count", len(topology["stable_ids"]))
    active, known = set(range(initial)), set(range(initial))
    invalid, last = 0, -1.0
    for event in events:
        a, b, result = (event.get(k) for k in ("parent_a", "parent_b", "result"))
        time, kind = event.get("time_myr", -1), event.get("kind")
        valid = (kind in ("split", "fusion") and type(result) is int and result >= 0
                 and a in active and b in active and result not in known
                 and np.isfinite(time) and last <= time <= manifest["time_myr"]
                 and event.get("cells", 0) > 0
                 and ((kind == "split" and a == b) or (kind == "fusion" and a != b)))
        invalid += int(not valid)
        if not valid:
            continue
        if kind == "fusion":
            active.remove(a)
            active.remove(b)
        active.add(result)
        known.add(result)
        last = time
    identities = topology["stable_ids"]
    invalid += int(len(set(identities)) != len(identities) or not set(identities) <= active)
    return {"available": True, "splits": sum(e.get("kind") == "split" for e in events),
            "fusions": sum(e.get("kind") == "fusion" for e in events), "invalid": invalid, "events": events}


def plate_topology(owners):
    """Connected material-owner regions on the native torus, not image outlines."""
    owners = np.asarray(owners)
    height, width = owners.shape
    visited = np.zeros(owners.shape, dtype=bool)
    components = {}
    for y in range(height):
        for x in range(width):
            if visited[y, x]:
                continue
            owner = int(owners[y, x])
            visited[y, x] = True
            queue, size = deque([(y, x)]), 0
            while queue:
                cy, cx = queue.popleft()
                size += 1
                for ny, nx in (((cy - 1) % height, cx), ((cy + 1) % height, cx),
                               (cy, (cx - 1) % width), (cy, (cx + 1) % width)):
                    if not visited[ny, nx] and owners[ny, nx] == owner:
                        visited[ny, nx] = True
                        queue.append((ny, nx))
            components.setdefault(owner, []).append(size)
    east, south = np.roll(owners, -1, 1), np.roll(owners, -1, 0)
    diagonal = np.roll(east, -1, 0)
    return {"components_per_plate": {str(p): len(sizes) for p, sizes in components.items()},
            "detached_components": sum(len(sizes) - 1 for sizes in components.values()),
            "detached_cells": sum(sum(sizes) - max(sizes) for sizes in components.values()),
            "boundary_faces": int((owners != east).sum() + (owners != south).sum()),
            "checkerboard_vertices": int(((owners == diagonal) & (east == south) & (owners != east)).sum())}


def terrain_detail(fields):
    z, owners = fields["elevation"], fields["plate_id"]
    residual = local_relief(z)
    interior = z > 0
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            interior &= (np.roll(owners, (dy, dx), axis=(0, 1)) == owners)
            interior &= np.roll(z > 0, (dy, dx), axis=(0, 1))
    ocean = z < 0
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            ocean &= (np.roll(owners, (dy, dx), axis=(0, 1)) == owners)
            ocean &= np.roll(z < 0, (dy, dx), axis=(0, 1))
    return {"local_relief_rms_m": float(np.sqrt(np.mean(residual**2))),
            "interior_land_cells": int(interior.sum()),
            "interior_land_relief_rms_m": float(np.sqrt(np.mean(residual[interior]**2))) if interior.any() else None,
            "interior_ocean_cells": int(ocean.sum()),
            "interior_ocean_relief_rms_m": float(np.sqrt(np.mean(residual[ocean]**2))) if ocean.any() else None}


def island_sizes(elevation):
    """Resolved four-connected landmasses on both periodic seams; not plate fragments."""
    land = elevation >= 0
    visited = np.zeros(land.shape, dtype=bool)
    h, w = land.shape
    sizes = []
    for y, x in zip(*np.nonzero(land)):
        if visited[y, x]:
            continue
        visited[y, x] = True
        queue, size = deque([(y, x)]), 0
        while queue:
            cy, cx = queue.popleft()
            size += 1
            for ny, nx in (((cy-1) % h, cx), ((cy+1) % h, cx), (cy, (cx-1) % w), (cy, (cx+1) % w)):
                if land[ny, nx] and not visited[ny, nx]:
                    visited[ny, nx] = True
                    queue.append((ny, nx))
        sizes.append(size)
    return {"landmass_count": len(sizes), "single_cell_islands": sum(s == 1 for s in sizes),
            "islands_2_to_64_cells": sum(2 <= s <= 64 for s in sizes),
            "island_cells_2_to_64": sum(s for s in sizes if 2 <= s <= 64),
            "largest_landmasses_cells": sorted(sizes, reverse=True)[:10]}


def hypsography_shape(curve):
    """Measure jumps in the common 0.5% quantile grid, excluding extreme tails."""
    z = np.asarray(curve["elevation_m"])
    q = np.asarray(curve["fraction_below"])
    steps = np.diff(z)
    selected = np.flatnonzero((q[:-1] >= .05) & (q[1:] <= .95))
    largest = int(selected[np.argmax(steps[selected])])
    return {"max_central_quantile_step_m": float(steps[largest]),
            "max_step_fraction_below": [float(q[largest]), float(q[largest + 1])],
            "elevation_span_35_40_percent_m": float(np.interp(.4, q, z) - np.interp(.35, q, z))}


def largest_connected_patch(mask):
    """Four-connected native cells, with both periodic seams included."""
    remaining = np.asarray(mask, dtype=bool).copy()
    height, width = remaining.shape
    largest = 0
    for y, x in zip(*np.nonzero(remaining)):
        if not remaining[y, x]:
            continue
        remaining[y, x] = False
        queue, size = deque([(y, x)]), 0
        while queue:
            cy, cx = queue.popleft()
            size += 1
            for ny, nx in (((cy - 1) % height, cx), ((cy + 1) % height, cx),
                           (cy, (cx - 1) % width), (cy, (cx + 1) % width)):
                if remaining[ny, nx]:
                    remaining[ny, nx] = False
                    queue.append((ny, nx))
        largest = max(largest, size)
    return largest


def near_convergence(boundary_type):
    convergence = boundary_type == 1
    near = np.zeros(convergence.shape, dtype=bool)
    for dy in range(-3, 4):
        for dx in range(-3, 4):
            if dx * dx + dy * dy <= 9:
                near |= np.roll(convergence, (dy, dx), axis=(0, 1))
    return near


def mountain_activity(manifest, fields):
    z = fields["elevation"]
    convergence = fields["boundary_type"] == 1
    near = near_convergence(fields["boundary_type"])
    highland = z >= 3000
    speeds = [float(np.hypot(*p["linear_velocity"])) for p in manifest["plates"] if "linear_velocity" in p]
    return {"area_above_3000_m": float(highland.mean()),
            "area_above_5000_m": float((z >= 5000).mean()),
            "area_above_8000_m": float((z >= 8000).mean()),
            "largest_highland_patch_cells": largest_connected_patch(highland),
            "highland_fraction_near_convergence": float(near[highland].mean()) if highland.any() else None,
            "land_near_convergence_m": stats(z[near & (z > 0)]),
            "convergent_boundary_cells": int(convergence.sum()),
            "mean_plate_speed": float(np.mean(speeds)) if speeds else None,
            "moving_plates": sum(speed > .01 for speed in speeds) if speeds else None,
            "plate_count": len(manifest["plates"])}


def diagnostics(manifest, fields):
    z, raw = fields["elevation"], fields["heightmap"]
    land = raw >= 1
    checks = {}
    for name in ("uplift_tendency", "subsidence_tendency", "accumulated_strain", "deformation_rate"):
        checks[name + "_outside_0_1"] = int(((fields[name] < 0) | (fields[name] > 1)).sum())
    for name, upper in (("crust_class", 3), ("boundary_type", 4), ("geologic_regime", 7),
                        ("deforming_region_type", 2), ("convergence_score", 100),
                        ("divergence_score", 100), ("shear_score", 100)):
        checks[name + "_invalid"] = int((fields[name] > upper).sum())
    checks["invalid_plate_owner"] = int((fields["plate_id"] >= len(manifest["plates"])).sum())
    checks["negative_crust_age"] = int((fields["crust_age_myr"] < 0).sum())
    age_limit = manifest.get("iteration_count", 0) * manifest.get("delta_time_myr", 0)
    checks["crust_age_exceeds_elapsed_and_prehistory"] = int(
        (fields["crust_age_myr"] > age_limit + max(.001, age_limit * 1e-6)).sum()) if age_limit > 0 else 0
    checks["negative_crust_thickness"] = int((fields["crust_thickness"] < 0).sum())
    if "rift_damage" in fields:
        checks["rift_damage_outside_0_1"] = int(((fields["rift_damage"] < 0) | (fields["rift_damage"] > 1)).sum())
    checks["oceanic_continental_collision"] = int(((fields["crust_class"] == 1) & (fields["geologic_regime"] == 2)).sum())
    motion = manifest.get("motion_diagnostics", {})
    budget = manifest.get("column_budget", {})
    checks["column_mass_budget_invalid"] = 0
    if budget.get("enabled"):
        required = ("initial_mass_kg", "mass_kg", "generated_ocean_kg", "subducted_ocean_kg",
                    "eroded_kg", "deposited_kg", "relative_mass_residual")
        valid = all(np.isfinite(budget.get(k, np.nan)) for k in required)
        if valid:
            initial = budget["initial_mass_kg"]
            valid = (initial > 0 and all(budget[k] >= 0 for k in required[:-1])
                     and abs(budget["relative_mass_residual"]) < 1e-9
                     and abs(budget["mass_kg"] - initial - budget["generated_ocean_kg"]
                             + budget["subducted_ocean_kg"] + budget.get("subducted_sediment_kg", 0)) < initial * 1e-9
                     and abs(budget["eroded_kg"] - budget["deposited_kg"]) < initial * 1e-9)
            for subset, total in (("generated_arc_kg", "generated_ocean_kg"), ("river_eroded_kg", "eroded_kg")):
                if subset in budget:
                    valid = valid and np.isfinite(budget[subset]) and 0 <= budget[subset] <= budget[total] + initial*1e-12
            for transfer in ("redistributed_continent_kg", "redistributed_sediment_kg", "noise_redistributed_kg"):
                if transfer in budget:
                    valid = valid and np.isfinite(budget[transfer]) and budget[transfer] >= 0
            for field in ("continental_crust_m", "oceanic_crust_m", "sediment_m"):
                valid = valid and field in fields and bool((fields[field] >= 0).all())
            if valid:
                cell_m = manifest["run_scenario"]["world_width_km"] * 1000 / manifest["width"]
                raster_mass = sum(float(fields[k].astype(float).sum()) * rho for k, rho in
                    (("continental_crust_m", 2800), ("oceanic_crust_m", 2900), ("sediment_m", 2700))) * cell_m**2
                valid = abs(raster_mass - budget["mass_kg"]) < initial * 1e-7
        checks["column_mass_budget_invalid"] = int(not valid)
    checks["force_balance_invalid"] = 0
    if motion.get("enabled"):
        required = ("relative_force_residual", "driving_power_w", "dissipation_w",
                    "mean_speed_cm_yr", "max_speed_cm_yr")
        valid = all(np.isfinite(motion.get(key, np.nan)) for key in required)
        if valid:
            power, loss = motion["driving_power_w"], motion["dissipation_w"]
            valid = (0 <= motion["relative_force_residual"] < 1e-8 and loss >= 0
                     and abs(power - loss) <= max(1, abs(power)) * 1e-8
                     and 0 <= motion["mean_speed_cm_yr"] <= motion["max_speed_cm_yr"] + 1e-12)
            if manifest.get('run_scenario', {}).get('consistent_motion'):
                valid = valid and 1 <= motion.get('nonlinear_iterations', 0) <= 96
                for key in ('velocity_consistency_cm_yr', 'max_velocity_consistency_cm_yr'):
                    value = motion.get(key, np.nan)
                    valid = valid and np.isfinite(value) and 0 <= value < 1e-4
                for subset, total in (('reversing_face_steps', 'persistent_face_steps'),
                                      ('strong_reversing_face_steps', 'strong_face_steps')):
                    valid = valid and 0 <= motion.get(subset, -1) <= motion.get(total, -1)
        checks["force_balance_invalid"] = int(not valid)
    topology = plate_topology(fields["plate_id"])
    checks["disconnected_plate_regions"] = topology["detached_components"]
    events = plate_events(manifest)
    checks["plate_ancestry_invalid"] = events["invalid"]
    invariants_pass = not any(checks.values())
    # These are codec representability checks, separate from physical plausibility.
    checks["above_uint16_codec"] = int((raw > 2).sum())
    checks["below_native_physical_floor"] = int((raw < np.float32(.1)).sum())
    dx = abs(z - np.roll(z, 1, axis=1))
    dy = abs(z - np.roll(z, 1, axis=0))
    local_mean = sum(np.roll(z, shift, axis=axis) for axis in (0, 1) for shift in (-1, 1)) / 4
    interior_land = land.copy()
    for axis in (0, 1):
        for shift in (-1, 1):
            interior_land &= np.roll(land, shift, axis=axis)
    land_jumps = np.concatenate((dx[land & np.roll(land, 1, axis=1)],
                                 dy[land & np.roll(land, 1, axis=0)]))
    ocean_jumps = np.concatenate((dx[~land & ~np.roll(land, 1, axis=1)],
                                  dy[~land & ~np.roll(land, 1, axis=0)]))
    quantiles = np.linspace(0, 1, 201)
    curve = {"fraction_below": quantiles.tolist(), "elevation_m": np.quantile(z, quantiles).tolist()}
    return {
        "checks": checks, "invariants_pass": invariants_pass,
        "codec_range_preserved": not (checks["above_uint16_codec"] or checks["below_native_physical_floor"]),
        "weighting": "equal-area native torus cells; not spherical latitude weights",
        "land_fraction": float(land.mean()), "elevation_m": stats(z), "land_elevation_m": stats(z[land]),
        "ocean_depth_m": stats(-z[~land]), "raw_height": stats(raw),
        "elevation_area_bands": {
            "shelf_0_200_m": float(((z<0)&(z>=-200)).mean()),
            "upper_ocean_200_2000_m": float(((z<-200)&(z>=-2000)).mean()),
            "deep_ocean_below_6000_m": float((z<-6000).mean()),
            "lowland_0_1000_m": float(((z>=0)&(z<1000)).mean()),
            "highland_above_2000_m": float((z>=2000).mean())},
        "ocean_at_convergent_boundary_m": stats(-z[(~land)&(fields['boundary_type']==1)]),
        "ocean_at_divergent_boundary_m": stats(-z[(~land)&(fields['boundary_type']==2)]),
        "crust_age_myr": stats(fields["crust_age_myr"]),
        "seafloor_age_myr": stats(fields["crust_age_myr"][(fields["crust_class"] == 1) & ~land]),
        "crust_column_sum_native": float(fields["crust_thickness"].astype(float).sum()),
        "cells_above_9000_m": int((z > 9000).sum()),
        "area_above_9000_m": float((z > 9000).mean()),
        "area_below_minus_11000_m": float((z < -11000).mean()),
        "cells_at_uint16_ceiling": int((raw >= 2).sum()),
        "interior_land_peaks_over_neighbour_mean_2000_m": int((interior_land & (z - local_mean > 2000)).sum()),
        "neighbour_jump_m": stats(np.concatenate((dx.ravel(), dy.ravel()))),
        "land_neighbour_jump_m": stats(land_jumps), "ocean_neighbour_jump_m": stats(ocean_jumps),
        "x_seam_jump_m": stats(dx[:, 0]), "y_seam_jump_m": stats(dy[0]),
        "hypsography": curve, "hypsography_shape": hypsography_shape(curve),
        "mountain_activity": mountain_activity(manifest, fields),
        "motion_budget": manifest.get("motion_diagnostics", {}),
        "column_budget": budget,
        "plate_topology": topology,
        "plate_events": events,
        "terrain_detail": terrain_detail(fields),
        "islands": island_sizes(fields["elevation"]),
        "coastal_terrain": coastal_terrain(manifest, fields),
        "topology_repairs": manifest.get("topology_diagnostics", {}),
        "evolution_mode": "continuous" if manifest.get("run_scenario", {}).get("evolution_duration_myr", 0) > 0 else "cycles",
        "plate_cells": np.bincount(fields["plate_id"].ravel(), minlength=len(manifest["plates"])).tolist(),
        "boundary_cells": np.bincount(fields["boundary_type"].ravel(), minlength=5).tolist(),
    }


def compare_fields(first, second):
    shared = first.keys() & second.keys()
    if (first.keys() ^ second.keys()) - OPTIONAL_DTYPES.keys() or any(first[k].shape != second[k].shape for k in shared):
        raise ValueError("Comparison requires matching native field stacks and grids")
    result = {}
    for name in sorted(shared):
        a, b = first[name], second[name]
        item = {"changed_cells": int(np.count_nonzero(a != b))}
        if name not in {"plate_id", "crust_class", "geologic_regime", "boundary_type",
                        "mountain_zones",
                        "boundary_segment_id", "nearest_boundary_id", "deforming_region_id", "deforming_region_type"}:
            delta = b.astype(float) - a.astype(float)
            item.update(mean_delta=float(delta.mean()), mae=float(abs(delta).mean()),
                        rmse=float(np.sqrt(np.mean(delta * delta))), max_abs=float(abs(delta).max()))
        result[name] = item
    return result
