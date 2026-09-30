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
import time
import zipfile

from run_storage import ensure_archive_directory

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def command(executable, recipe):
    args = [str(executable), "batch-export"]
    for key in ("seed", "width", "height"):
        args += ["--" + key, str(recipe[key])]
    for key, value in recipe["scenario"].items():
        args += ["--" + key, str(value)]
    args += ["--checkpoint-every-updates", str(recipe["checkpoint_every_updates"])]
    return args


def field_hashes(bundle):
    manifest = json.loads((bundle / "manifest.json").read_text())
    return {name: digest(bundle / filename) for name, filename in manifest["files"].items()}


def sample_run(run_id, interval):
    """Add denser observations by replaying the archived binary; never alter dynamics."""
    if interval < 1:
        raise ValueError("Checkpoint interval must be positive")
    manifest_path = ROOT / f"runs/manifests/tectonics/{run_id}.json"
    record = json.loads(manifest_path.read_text())
    archive = ROOT / f"runs/reports/tectonics/{run_id}"
    executable = archive / "tectonic_pipeline.exe"
    if digest(executable) != record["executable_sha256"]:
        raise ValueError("Archived executable hash mismatch")
    root = ROOT / f"runs/checkpoints/tectonics/{run_id}"
    output = root / f"sampling-{interval}"
    if output.exists():
        raise ValueError("Sampling output already exists")
    recipe = dict(record["recipe"], checkpoint_every_updates=interval)
    args = command(executable, recipe) + ["--output", str(output)]
    with (archive / f"sampling-{interval}.log").open("x", encoding="utf-8") as log:
        subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    batch = json.loads((output / "batch_manifest.json").read_text())
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
        raise ValueError(f"Build tectonic_pipeline first: {executable}")
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
            subprocess.run(["pwsh", "-NoProfile", "-File",
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
    # Include the actual working files, including uncommitted task changes.
    inputs = [ROOT / "CMakeLists.txt", ROOT / "CMakePresets.json"]
    for folder in ("src/simulations/geology", "tools/tectonics", "scripts/benchmarks",
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
    native_command = command(archived_executable, recipe) + ["--output", str(checkpoint_dir)]
    record["command"] = native_command
    write_json(manifest_path, record)
    start = time.monotonic()
    try:
        with log_path.open("x", encoding="utf-8") as log:
            log.write(subprocess.list2cmdline(native_command) + "\n")
            log.flush()
            subprocess.run(native_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
        record["elapsed_seconds"] = time.monotonic() - start
        batch = json.loads((checkpoint_dir / "batch_manifest.json").read_text())
        record["snapshots"] = batch["snapshots"]
        record["field_sha256"] = {item["label"]: field_hashes(checkpoint_dir / item["directory"])
                                  for item in record["snapshots"]}
        if args.verify_replay:
            final = next(item for item in record["snapshots"] if item["label"] == "final")
            replay = archive / "replay-final"
            replay_command = command(archived_executable, recipe)
            # Replay is verification of this numbered run, not a second experiment.
            replay_command[1] = "inspect"
            replay_command = replay_command[:-2] + ["--final", "--bundle-output", str(replay)]
            with (archive / "replay.log").open("x", encoding="utf-8") as log:
                subprocess.run(replay_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
            expected = record["field_sha256"][final["label"]]
            actual = field_hashes(replay)
            replay_manifest = json.loads((replay / "manifest.json").read_text())
            original = json.loads((checkpoint_dir / final["directory"] / "manifest.json").read_text())
            metadata_keys = ("run_scenario", "time_myr", "update_count", "iteration_count", "cycle_count")
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
    parser.add_argument("--executable", type=Path, default=ROOT / "out/build/headless-modules/bin/tectonic_pipeline.exe")
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
