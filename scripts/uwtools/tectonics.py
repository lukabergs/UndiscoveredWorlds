"""Tectonic snapshots: run bundles of the tectonics recipe and the former snapshot bundles.

The reports keep the field names of the former snapshot files (BUNDLE_FIELDS maps them to
catalogue ids). A run's snapshots are its checkpoint bundles and its final bundle.
"""
import json
from pathlib import Path

from .bundle import SCHEMAS as RUN_BUNDLE_SCHEMAS, read_manifest
from .hashing import sha256

# Catalogue fields of a run bundle (uw-run-bundle/v1 or v2) under the names the former
# snapshot files had; the reports keep using those names.
BUNDLE_FIELDS = {
    "heightmap": "geosphere.tectonic_height", "plate_id": "geosphere.plate_id",
    "crust_age_steps": "geosphere.crust_age_steps", "crust_age_myr": "geosphere.crust_age",
    "crust_thickness": "geosphere.crust_thickness", "crust_class": "geosphere.crust_class",
    "uplift_tendency": "geosphere.uplift_tendency", "subsidence_tendency": "geosphere.subsidence_tendency",
    "rock_uplift_rate": "geosphere.rock_uplift_rate",
    "shortening": "geosphere.shortening", "shortening_direction": "geosphere.shortening_direction",
    "accumulated_strain": "geosphere.accumulated_strain", "boundary_type": "geosphere.boundary_type",
    "boundary_distance": "geosphere.boundary_distance", "boundary_segment_id": "geosphere.boundary_segment",
    "nearest_boundary_id": "geosphere.nearest_boundary", "deforming_region_id": "geosphere.deforming_region",
    "deforming_region_type": "geosphere.deforming_region_type", "deformation_rate": "geosphere.deformation_rate",
    "deformation_velocity_x": "geosphere.deformation_velocity_x",
    "deformation_velocity_y": "geosphere.deformation_velocity_y",
    "convergence_score": "geosphere.tectonic_convergence", "divergence_score": "geosphere.tectonic_divergence",
    "shear_score": "geosphere.tectonic_shear", "geologic_regime": "geosphere.geologic_regime",
    "rift_damage": "geosphere.rift_damage", "plume_influence": "geosphere.plume_influence",
    "continental_crust_m": "geosphere.continental_crust", "oceanic_crust_m": "geosphere.oceanic_crust",
    "sediment_m": "geosphere.sediment", "continental_area": "geosphere.continental_area",
}
SNAPSHOT_TABLE = "geosphere.tectonic_snapshot"


def is_run_bundle(bundle):
    return read_manifest(bundle).get("schema") in RUN_BUNDLE_SCHEMAS


def bundle_files(bundle):
    """Field name -> file of a snapshot, for a run bundle or a former snapshot bundle."""
    manifest = read_manifest(bundle)
    if manifest.get("schema") not in RUN_BUNDLE_SCHEMAS:
        return dict(manifest["files"])
    entries = {entry["id"]: entry["file"] for entry in manifest["fields"]}
    return {name: entries[identifier] for name, identifier in BUNDLE_FIELDS.items() if identifier in entries}


def snapshot_manifest(bundle):
    """The engine's snapshot record: label, time, scenario, diagnostics, plates."""
    bundle = Path(bundle)
    manifest = read_manifest(bundle)
    if manifest.get("schema") not in RUN_BUNDLE_SCHEMAS:
        return manifest
    table = next(entry for entry in manifest["tables"] if entry["id"] == SNAPSHOT_TABLE)
    return json.loads((bundle / table["file"]).read_text())


def field_hashes(bundle):
    return {name: sha256(Path(bundle) / filename) for name, filename in bundle_files(bundle).items()}


def run_snapshots(directory, checkpoints="checkpoints", final="bundle"):
    """The bundles of one run in time order: label, update and iteration count, time, directory (relative)."""
    directory = Path(directory)
    found = []
    folders = sorted((directory / checkpoints).iterdir()) if (directory / checkpoints).is_dir() else []
    for folder in folders + [directory / final]:
        if not (folder / "manifest.json").is_file():
            continue
        record = snapshot_manifest(folder)
        found.append({"update_count": record["update_count"], "iteration_count": record["iteration_count"],
                      "time_myr": record["time_myr"], "directory": folder.relative_to(directory).as_posix()})
    found.sort(key=lambda item: (item["update_count"], item["directory"] == final))
    snapshots = []
    for item in found:
        # The end of the run is both the last periodic checkpoint and the final bundle: keep the final one.
        if snapshots and snapshots[-1]["update_count"] == item["update_count"]:
            snapshots.pop()
        label = "final" if item["directory"] == final else "initial" if item["update_count"] == 0 \
            else f"update{item['update_count']:06d}"
        snapshots.append(dict(label=label, **item))
    return snapshots
