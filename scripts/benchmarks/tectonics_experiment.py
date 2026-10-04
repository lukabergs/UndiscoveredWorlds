# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Run one numbered tectonics experiment; preserve its recipe, binary and native fields."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths
from run_storage import POWERSHELL, ensure_archive_directory
import run_log

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


# Catalogue fields of a run bundle (uw-run-bundle/v1) under the names the former
# snapshot files had; the reports keep using those names.
BUNDLE_FIELDS = {
    "heightmap": "geosphere.tectonic_height", "plate_id": "geosphere.plate_id",
    "crust_age_steps": "geosphere.crust_age_steps", "crust_age_myr": "geosphere.crust_age",
    "crust_thickness": "geosphere.crust_thickness", "crust_class": "geosphere.crust_class",
    "uplift_tendency": "geosphere.uplift_tendency", "subsidence_tendency": "geosphere.subsidence_tendency",
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
    "sediment_m": "geosphere.sediment",
}
SNAPSHOT_TABLE = "geosphere.tectonic_snapshot"
RUN_BUNDLE_SCHEMA = "uw-run-bundle/v1"
# world_pipeline needs these beside it when it runs from a run's archive.
RUNTIME_LIBRARIES = ("libpng16.dll", "zlib1.dll")


def command(executable, recipe, output):
    """world_pipeline run --recipe tectonics: an initial bundle, one every checkpoint interval and the final one."""
    if recipe["height"] * 2 != recipe["width"]:
        raise ValueError("The world grid is W x W/2")
    scenario = recipe["scenario"]
    interval = recipe["checkpoint_every_updates"] * float(scenario.get("delta-time-myr", 1.0))
    args = [str(executable), "run", "--recipe", "tectonics", "--seed", str(recipe["seed"]),
            "--resolution", str(recipe["width"])]
    for key, value in scenario.items():
        args += ["--" + key, str(value)]
    args += ["--checkpoint-at", "geology.geosphere_initialization", "--checkpoint-every", f"{interval:g}",
             "--checkpoint-dir", str(Path(output) / "snapshots"), "--bundle", str(Path(output) / "final")]
    return args


def legacy_command(executable, recipe):
    """tectonic_pipeline batch-export, for the archived executable of a run made before world_pipeline took over."""
    args = [str(executable), "batch-export"]
    for key in ("seed", "width", "height"):
        args += ["--" + key, str(recipe[key])]
    for key, value in recipe["scenario"].items():
        args += ["--" + key, str(value)]
    args += ["--checkpoint-every-updates", str(recipe["checkpoint_every_updates"])]
    return args


def is_run_bundle(bundle):
    return json.loads((Path(bundle) / "manifest.json").read_text()).get("schema") == RUN_BUNDLE_SCHEMA


def bundle_files(bundle):
    """Field name -> file of a snapshot, for a run bundle or a former snapshot bundle."""
    manifest = json.loads((Path(bundle) / "manifest.json").read_text())
    if manifest.get("schema") != RUN_BUNDLE_SCHEMA:
        return dict(manifest["files"])
    entries = {entry["id"]: entry["file"] for entry in manifest["fields"]}
    return {name: entries[identifier] for name, identifier in BUNDLE_FIELDS.items() if identifier in entries}


def snapshot_manifest(bundle):
    """The engine's snapshot record: label, time, scenario, diagnostics, plates."""
    bundle = Path(bundle)
    manifest = json.loads((bundle / "manifest.json").read_text())
    if manifest.get("schema") != RUN_BUNDLE_SCHEMA:
        return manifest
    table = next(entry for entry in manifest["tables"] if entry["id"] == SNAPSHOT_TABLE)
    return json.loads((bundle / table["file"]).read_text())


def field_hashes(bundle):
    return {name: digest(Path(bundle) / filename) for name, filename in bundle_files(bundle).items()}


def run_snapshots(output):
    """The bundles of one run in time order, described as the former batch manifest did."""
    output = Path(output)
    found = []
    for directory in sorted((output / "snapshots").iterdir()) + [output / "final"]:
        if not (directory / "manifest.json").is_file():
            continue
        record = snapshot_manifest(directory)
        found.append({"update_count": record["update_count"], "iteration_count": record["iteration_count"],
                      "time_myr": record["time_myr"], "directory": directory.relative_to(output).as_posix()})
    found.sort(key=lambda item: (item["update_count"], item["directory"] == "final"))
    snapshots = []
    for item in found:
        # The end of the run is both the last periodic checkpoint and the final bundle: keep the final one.
        if snapshots and snapshots[-1]["update_count"] == item["update_count"]:
            snapshots.pop()
        label = "final" if item["directory"] == "final" else "initial" if item["update_count"] == 0 \
            else f"update{item['update_count']:06d}"
        snapshots.append(dict(label=label, **item))
    return snapshots


def sample_run(run_id, interval):
    """Add denser observations by replaying the archived binary; never alter dynamics."""
    if interval < 1:
        raise ValueError("Checkpoint interval must be positive")
    manifest_path = ROOT / f"runs/manifests/tectonics/{run_id}.json"
    record = json.loads(manifest_path.read_text())
    archive = ROOT / f"runs/reports/tectonics/{run_id}"
    legacy = (archive / "tectonic_pipeline.exe").is_file()
    executable = archive / ("tectonic_pipeline.exe" if legacy else "world_pipeline.exe")
    if digest(executable) != record["executable_sha256"]:
        raise ValueError("Archived executable hash mismatch")
    root = ROOT / f"runs/checkpoints/tectonics/{run_id}"
    output = root / f"sampling-{interval}"
    if output.exists():
        raise ValueError("Sampling output already exists")
    recipe = dict(record["recipe"], checkpoint_every_updates=interval)
    args = legacy_command(executable, recipe) + ["--output", str(output)] if legacy else command(executable, recipe, output)
    with (archive / f"sampling-{interval}.log").open("x", encoding="utf-8") as log:
        subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    batch = json.loads((output / "batch_manifest.json").read_text()) if legacy else {"snapshots": run_snapshots(output)}
    for item in batch["snapshots"]:
        hashes = field_hashes(output / item["directory"])
        if item["label"] in ("initial", "final"):
            if hashes != record["field_sha256"][item["label"]]:
                raise ValueError("Sampling replay mismatch; original run retained")
    existing = {item["update_count"] for item in record["snapshots"]}
    for item in batch["snapshots"]:
        if item["update_count"] not in existing:
            item["directory"] = f"sampling-{interval}/" + item["directory"]
            record["snapshots"].append(item)
            record["field_sha256"][item["label"]] = field_hashes(root / item["directory"])
    record["snapshots"].sort(key=lambda item: item["update_count"])
    record.setdefault("sampling_replays", []).append({"command": args, "endpoints_exact_match": True})
    write_json(manifest_path, record)
    print(f"Added verified {interval}-update sampling to tectonics run {run_id}")


def run(args):
    recipe = json.loads(args.recipe.read_text())
    executable = args.executable.absolute()
    if not executable.is_file():
        raise ValueError(f"Build world_pipeline first: {executable}")
    if recipe["seed"] != 1 and not args.allow_seed_comparison:
        raise ValueError("Use --allow-seed-comparison for an explicitly requested seed experiment")
    registry_path = ROOT / "runs/registry/tectonics.json"
    registry = json.loads(registry_path.read_text()) if registry_path.exists() else {
        "schema": "tectonics-run-registry/v1", "seed": 1, "baseline": None, "runs": []}
    # A seed sensitivity experiment must not replace the accepted baseline.
    run_id = max((item["id"] for item in registry["runs"]), default=0) + 1
    archive = ensure_archive_directory(ROOT / "runs/reports/tectonics") / str(run_id)
    # Route new module roots through the same storage policy as climate.
    if (ROOT / "runs/storage.json").exists():
        for kind in ("maps", "fields", "checkpoints", "logs", "diagnostics"):
            parent = ROOT / "runs" / kind
            if parent.is_junction():
                (parent / "tectonics").mkdir(exist_ok=True)
                continue
            subprocess.run([POWERSHELL, "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                            str(ROOT / "scripts/benchmarks/configure-run-storage.ps1"),
                            "-RelativePaths", f"runs/{kind}/tectonics"], cwd=ROOT, check=True)
    checkpoint_dir = ROOT / f"runs/checkpoints/tectonics/{run_id}"
    if archive.exists() or checkpoint_dir.exists():
        raise ValueError("Run artifacts already exist; refusing to overwrite evidence")
    archive.mkdir(parents=True)
    log_path = ROOT / f"runs/logs/tectonics/{run_id}.log"
    manifest_path = ROOT / f"runs/manifests/tectonics/{run_id}.json"
    record = {"id": run_id, "status": "running", "recipe": recipe,
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "git_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
              "executable_sha256": digest(executable),
              "geometry": {"registration": "cell-centred", "layout": "row-major",
                           "topology": "planar torus", "wrap_x": True, "wrap_y": True,
                           "statistics_weighting": "equal native cell areas"},
              "stopping_rule": "continuous evolution until recipe duration; record actual time and updates",
              "determinism_scope": "same source, recipe, executable, platform and toolchain",
              "native_scenario_metadata_incomplete": True}
    registry["runs"].append({"id": run_id, "status": "running", "manifest": manifest_path.relative_to(ROOT).as_posix()})
    write_json(registry_path, registry)
    write_json(archive / "recipe.json", recipe)
    shutil.copy2(executable, archive / executable.name)
    for library in RUNTIME_LIBRARIES:
        if (executable.parent / library).is_file():
            shutil.copy2(executable.parent / library, archive / library)
    # Include the actual working files, including uncommitted task changes.
    inputs = [ROOT / "CMakeLists.txt", ROOT / "CMakePresets.json"]
    for folder in ("src/simulations/geology", "src/validation/tectonics", "scripts/benchmarks",
                   "tests/simulations/geology", "configs"):
        inputs += [p for p in (ROOT / folder).rglob("*")
                   if p.is_file() and p.suffix in {".cpp", ".hpp", ".h", ".py", ".json", ".txt", ".html", ".js", ".cmake"}]
    with zipfile.ZipFile(archive / "source.zip", "w", zipfile.ZIP_DEFLATED) as package:
        for path in sorted(set(inputs)):
            package.write(path, path.relative_to(ROOT).as_posix())
    record["source_zip_sha256"] = digest(archive / "source.zip")
    build_dir = executable.parent.parent
    cache = build_dir / "CMakeCache.txt"
    if cache.exists():
        shutil.copy2(cache, archive / "CMakeCache.txt")
    # Always execute the archived binary, so later builds cannot change this run
    # or its replay after provenance has already been recorded.
    archived_executable = archive / executable.name
    native_command = command(archived_executable, recipe, checkpoint_dir)
    record["command"] = native_command
    write_json(manifest_path, record)
    start = time.monotonic()
    try:
        with log_path.open("x", encoding="utf-8") as log:
            log.write(run_log.header("tectonics", archived_executable, native_command))
            log.flush()
            subprocess.run(native_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
        record["elapsed_seconds"] = time.monotonic() - start
        record["snapshots"] = run_snapshots(checkpoint_dir)
        record["field_sha256"] = {item["label"]: field_hashes(checkpoint_dir / item["directory"])
                                  for item in record["snapshots"]}
        if args.verify_replay:
            final = next(item for item in record["snapshots"] if item["label"] == "final")
            replay = archive / "replay-final"
            # Replay is verification of this numbered run, not a second experiment: the final bundle only.
            replay_command = command(archived_executable, recipe, replay)
            first = replay_command.index("--checkpoint-at")
            replay_command = replay_command[:first] + ["--bundle", str(replay)]
            with (archive / "replay.log").open("x", encoding="utf-8") as log:
                subprocess.run(replay_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
            expected = record["field_sha256"][final["label"]]
            actual = field_hashes(replay)
            replay_manifest = snapshot_manifest(replay)
            original = snapshot_manifest(checkpoint_dir / final["directory"])
            metadata_keys = ("run_scenario", "time_myr", "update_count", "iteration_count")
            if "motion_diagnostics" in original:
                metadata_keys += ("motion_diagnostics",)
            if "column_budget" in original:
                metadata_keys += ("column_budget",)
            if "topology_diagnostics" in original:
                metadata_keys += ("topology_diagnostics",)
            exact = expected == actual and all(original[k] == replay_manifest[k] for k in metadata_keys)
            record["replay"] = {"exact_match": exact, "field_sha256": actual, "command": replay_command}
            if not exact:
                raise RuntimeError("Deterministic replay mismatch; evidence retained")
        else:
            record["replay"] = {"status": "skipped", "reason": "Replay verification was not requested"}
        record["status"] = "complete"
    except Exception as error:
        record["status"] = "failed"
        record["error"] = str(error)
        raise
    finally:
        write_json(manifest_path, record)
        registry["runs"][-1]["status"] = record["status"]
        if record["status"] == "complete" and registry["baseline"] is None:
            registry["baseline"] = run_id
        write_json(registry_path, registry)
    print(f"Tectonics run {run_id}: {manifest_path}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recipe", type=Path, default=ROOT / "configs/tectonics-baseline.json")
    parser.add_argument("--executable", type=Path, default=paths.executable("world_pipeline"))
    parser.add_argument("--verify-replay", action=argparse.BooleanOptionalAction, default=False,
                        help="Opt in to a second simulation for exact replay verification (default: off)")
    parser.add_argument("--allow-seed-comparison", action="store_true")
    parser.add_argument("--sample-run", type=int, help="Resample an archived run using its original executable")
    parser.add_argument("--checkpoint-every", type=int, default=5)
    args = parser.parse_args()
    if args.sample_run is not None:
        sample_run(args.sample_run, args.checkpoint_every)
    else:
        run(args)


if __name__ == "__main__":
    main()
