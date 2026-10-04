"""Migrate existing climate maps to one-based seasonal paths, without overwrites.

Dry run by default. --apply journals and stages file renames so January/April
folders cannot overwrite each other. Re-running resumes an interrupted move.
"""

import argparse
import json
from pathlib import Path

from reference_map_layout import product_path

ROOT = Path(__file__).resolve().parents[2]
JOURNAL = ROOT / "runs/manifests/climate-map-layout.json"


def legacy_run_directory(parts):
    if parts == ["koppen"]:
        return "koppen"
    if parts == ["rain"]:
        return "rain/annual/s"
    if parts == ["air_temperature"]:
        return "air_temp/annual/s"
    first = parts[0]
    if first == "wind":
        style = {"part": "particles", "err": "error", "cons": "consistency",
                 "u": "east", "v": "north"}.get(parts[1], parts[1])
        tail = "/".join(parts[2:]).replace("_", "/").split("/")
        return f"wind/{style}/{int(tail[0])+1}/{tail[1]}"
    tail = parts[-1].split("_")
    season = int(tail[0]) + 1
    if first in ("s_div", "u_div"):
        return f"wind/divergence/{season}/{first[0]}"
    if first == "moisture":
        quantity = {"conv": "convergence", "column_water": "water", "flux": "flux"}[parts[1]]
        return f"moisture/{quantity}/{season}/{tail[1] if len(tail)>1 else 'column'}"
    return {"sst": f"sea_temp/{season}/s", "ocean": f"ocean/current/{season}/s",
            "ascent": f"wind/ascent/{season}/u", "heating": f"energy/heating/{season}/column"}[first]


def checked(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT) or path == ROOT:
        raise ValueError(f"Path escaped workspace: {path}")
    return path


def plan():
    manifest_path = ROOT / "refs/processed/metadata/reference-maps/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    moves = {}
    for resolution in manifest["outputs"]:
        for product in resolution["products"]:
            for record in (*product["maps"], product["csv"]):
                old = record["file"]
                if not old.startswith("climate/"):
                    continue
                kind = {".png": "maps", ".tif": "fields", ".csv": "csv"}[Path(old).suffix]
                new = product_path(product["name"], resolution["width"], kind).as_posix()
                if old != new:
                    moves["refs/processed/" + old] = "refs/processed/" + new
    for kind in ("maps", "fields"):
        base = ROOT / "runs" / kind / "climate"
        for source in base.rglob("*"):
            if not source.is_file():
                continue
            parts = source.relative_to(base).parts
            if source.suffix == ".txt":
                destination = ROOT / "runs/work/legacy-map-notes" / kind / Path(*parts)
            elif "unregistered" in parts:
                destination = ROOT / "runs/work/legacy-maps" / kind / Path(*parts)
            else:
                if not source.stem.isdigit():
                    raise ValueError(f"Unknown run identity: {source}")
                directory = legacy_run_directory(list(parts[:-1]))
                destination = base / directory / source.name
            if source != destination:
                moves[source.relative_to(ROOT).as_posix()] = destination.relative_to(ROOT).as_posix()
    if len(set(moves.values())) != len(moves):
        raise ValueError("Two source files map to the same destination")
    for source, destination in moves.items():
        if not checked(source).is_file():
            raise ValueError(f"Missing source: {source}")
        if checked(destination).exists() and destination not in moves:
            raise ValueError(f"Refusing to overwrite {destination}")
    return moves


def apply(journal):
    moves = journal["moves"]
    staging = ROOT / "runs/work/map-layout-staging"
    # Stage every source before writing destinations: 0 -> 1 -> 2 is a cycle.
    if journal["phase"] == "staging":
        for index, (source, destination) in enumerate(moves.items()):
            src, temp = checked(source), checked(staging / str(index))
            if temp.exists():
                continue
            temp.parent.mkdir(parents=True, exist_ok=True)
            src.rename(temp)
        journal["phase"] = "destinations"
        JOURNAL.write_text(json.dumps(journal, indent=2) + "\n", encoding="utf-8")
    for index, destination in enumerate(moves.values()):
        temp, dst = checked(staging / str(index)), checked(destination)
        if not temp.exists() and dst.is_file():
            continue
        if dst.exists():
            raise ValueError(f"Refusing to overwrite {dst}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        temp.rename(dst)
    manifest_path = ROOT / "refs/processed/metadata/reference-maps/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    def update(value):
        if isinstance(value, dict):
            for key, child in value.items():
                if key == "file" and "refs/processed/" + str(child) in moves:
                    value[key] = moves["refs/processed/" + child].removeprefix("refs/processed/")
                else:
                    update(child)
        elif isinstance(value, list):
            for child in value:
                update(child)

    update(manifest)
    manifest["map_layout_version"] = 1
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    # Remove only now-empty, explicitly checked directories under migrated roots.
    for root in (ROOT / "runs/maps/climate", ROOT / "runs/fields/climate", staging):
        for directory in sorted(root.rglob("*"), key=lambda p: len(p.parts), reverse=True):
            if directory.is_dir() and not any(directory.iterdir()):
                checked(directory).rmdir()
    journal["phase"] = "complete"
    JOURNAL.write_text(json.dumps(journal, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    if JOURNAL.exists():
        journal = json.loads(JOURNAL.read_text(encoding="utf-8"))
        if journal["phase"] == "complete":
            print("Climate map layout already migrated.")
            raise SystemExit(0)
    else:
        journal = {"version": 1, "phase": "staging", "moves": plan()}
    print(f"{len(journal['moves'])} file renames; destination collisions checked.")
    if args.apply:
        JOURNAL.parent.mkdir(parents=True, exist_ok=True)
        JOURNAL.write_text(json.dumps(journal, indent=2) + "\n", encoding="utf-8")
        apply(journal)
        print("Migrated maps and manifest; raster contents unchanged.")
