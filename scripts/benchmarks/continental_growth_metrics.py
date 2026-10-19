"""Continental area, thickness and freeboard of tectonics bundles over time (terrain roadmap T4.12).

Reads a bundle and its `<bundle>-checkpoints` directory and prints one row per checkpoint.
Every measure is named; all are unweighted shares of the engine's planar grid cells:

- col_area: cells whose column `continental_area` is at least 0.5 (needs the field, T4.12 round two)
- col_area_sum: mean column `continental_area`
- col_thick_km: sum `continental_crust` / sum `continental_area` (continental volume per continental area);
  col_thick_sed_km adds the sediment on the continental area (sediment taken as spread evenly over the cell);
  surface_continent_km is the surface mean
- c20_area: cells with `continental_crust` (continent per cell area) of at least 20 km; c20_thick_km their mean
- ghost: cells with `continental_area` >= 0.5 and under 1 km of continent
- land: cells above the engine's sea level; bands of land in m
- cont_land_median: median height of land on continental columns (col_area, or c20 without the field)
- interior_median: median height of land cells more than --interior-cells cells (chessboard distance; a coastal
  cell is at 1) from the sea
- shelf_band, slope_band: cells from -200 to 0 m and from -1000 to -200 m; shelf_thick_km, slope_thick_km their
  continental volume per continental area (their mean continent without the area field); shelf_sed_m, slope_sed_m
  their mean sediment
- f, residuals, shortfall: the snapshot's column budget
- collided_area_rate, stretched_area_rate, arc_area_rate, capped_area_rate: the budget's cumulative continental
  area terms (capped: oceanized minus continentalized) between this row and the previous one in time (the first
  from t = 0), as shares of the surface per Myr; only when the snapshot has them

Usage: uv run python scripts/benchmarks/continental_growth_metrics.py BUNDLE [BUNDLE ...] [--json OUT]
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import ndimage

BANDS = [(0, 200), (200, 500), (500, 1000), (1000, 2000), (2000, None)]
SEA_BANDS = {"shelf": (-200, 0), "slope": (-1000, -200)}
# Cumulative budget area terms, m²: rate name -> (added keys, subtracted keys). Signs: collided and capped
# positive = continental area lost, stretched and arc positive = continental area gained.
AREA_TERMS = {
    "collided_area_rate": (["collided_continental_area_m2"], []),
    "stretched_area_rate": (["stretched_continental_area_m2"], []),
    "arc_area_rate": (["arc_continental_area_m2"], []),
    "capped_area_rate": (["oceanized_continental_area_m2"], ["continentalized_oceanic_area_m2"]),
}


def load(directory: Path, name: str):
    path = directory / "fields" / "geosphere" / f"{name}.npy"
    return np.load(path) if path.exists() else None


def measure(directory: Path, interior_cells: int) -> dict:
    manifest = directory / "manifest.json"
    if manifest.exists() and json.loads(manifest.read_text(encoding="utf-8")).get("grid", {}).get("kind", "latlon") != "latlon":
        raise SystemExit(f"{directory}: a cube bundle; these torus measures (tiling, the torus surface) are ported "
                         "with the engine (E5/E6)")
    continent = load(directory, "continental_crust").astype(np.float64)
    area = load(directory, "continental_area")
    elevation = load(directory, "base_elevation").astype(np.float64) - float(load(directory, "sea_level"))
    row: dict = {}
    c20 = continent >= 20000
    row["c20_area"] = float(c20.mean())
    row["c20_thick_km"] = float(continent[c20].mean() / 1000) if c20.any() else 0.0
    if area is not None:
        area = area.astype(np.float64)
        continental = area >= 0.5
        row["col_area"] = float(continental.mean())
        row["col_area_sum"] = float(area.mean())
        row["col_thick_km"] = float(continent.sum() / max(area.sum(), 1e-30) / 1000)
        # With the sediment on the continental area (CRUST1.0's 36.1 km includes 1.6 km of sediment).
        sediment = load(directory, "sediment").astype(np.float64)
        row["col_thick_sed_km"] = float((continent + sediment * area).sum() / max(area.sum(), 1e-30) / 1000)
        row["surface_continent_km"] = float(continent.mean() / 1000)
        row["ghost"] = float((continental & (continent < 1000)).mean())
        thickness = np.where(area > 0, continent / np.maximum(area, 1e-30), 0)
        row["over50_share"] = float((thickness[continental] > 50000).mean()) if continental.any() else 0.0
        row["over70_share"] = float((thickness[continental] > 70000).mean()) if continental.any() else 0.0
    else:
        continental = c20
    land = elevation > 0
    row["land"] = float(land.mean())
    for low, high in BANDS:
        band = land & (elevation >= low) & ((elevation < high) if high is not None else True)
        row[f"land_{low}_{high or 'up'}"] = float(band.mean())
    sediment_m = load(directory, "sediment")
    for name, (low, high) in SEA_BANDS.items():
        band = (elevation >= low) & (elevation < high)
        row[f"{name}_band"] = float(band.mean())
        if not band.any():
            row[f"{name}_thick_km"] = row[f"{name}_sed_m"] = float("nan")
            continue
        if area is not None:
            row[f"{name}_thick_km"] = float(continent[band].sum() / max(area[band].sum(), 1e-30) / 1000)
        else:
            row[f"{name}_thick_km"] = float(continent[band].mean() / 1000)
        row[f"{name}_sed_m"] = float(sediment_m[band].mean()) if sediment_m is not None else float("nan")
    cont_land = land & continental
    row["cont_land_median"] = float(np.median(elevation[cont_land])) if cont_land.any() else float("nan")
    # Distance to the sea on the torus: tile 3 x 3, take the centre.
    sea = ~land
    tiled = np.tile(sea, (3, 3))
    distance = ndimage.distance_transform_cdt(~tiled, metric="chessboard")
    rows, columns = sea.shape
    distance = distance[rows:2 * rows, columns:2 * columns]
    interior = land & (distance > interior_cells)
    row["interior_share"] = float(interior.mean())
    row["interior_median"] = float(np.median(elevation[interior])) if interior.any() else float("nan")
    snapshot = directory / "tables" / "geosphere" / "tectonic_snapshot.json"
    if snapshot.exists():
        data = json.loads(snapshot.read_text())
        budget = data.get("column_budget", {})
        row["f"] = budget.get("arc_continent_fraction")
        row["mass_residual"] = budget.get("relative_mass_residual")
        row["felsic_residual"] = budget.get("felsic_relative_residual")
        row["shortfall_kg"] = budget.get("continental_growth_shortfall_kg")
        row["felsic_kg"] = budget.get("felsic_kg")
        row["mean_speed_cm_yr"] = data.get("motion_diagnostics", {}).get("mean_speed_cm_yr")
        row["plates"] = len(data.get("plates", []))
        # Cumulative area terms and what their rates need (removed again by series).
        # The engine's surface: the sphere of run_scenario.radius_km (the cube engine, contract
        # schema 8), or the torus of world_width_km in records before it.
        scenario = data.get("run_scenario", {})
        width_m = float(scenario.get("world_width_km", 0)) * 1000
        cells_wide, cells_high = data.get("width", 0), data.get("height", 0)
        if width_m > 0 and cells_wide:
            row["_surface_m2"] = width_m * width_m * cells_high / cells_wide
        elif scenario.get("radius_km"):
            row["_surface_m2"] = 4.0 * math.pi * (float(scenario["radius_km"]) * 1000) ** 2
        else:
            row["_surface_m2"] = None
        row["_engine_time_myr"] = data.get("time_myr")
        for name, (added, subtracted) in AREA_TERMS.items():
            if all(key in budget for key in added + subtracted):
                row[f"_{name}"] = sum(budget[k] for k in added) - sum(budget[k] for k in subtracted)
    return row


def add_area_rates(rows: list[dict]) -> None:
    """Area terms as shares of the surface per Myr since the previous row in time (t = 0 for the first)."""
    for index, row in enumerate(rows):
        time, surface = row.get("_engine_time_myr"), row.get("_surface_m2")
        for name in AREA_TERMS:
            value = row.get(f"_{name}")
            if value is None or time is None or not surface:
                continue
            earlier = [r for r in rows[:index] if r.get("_engine_time_myr") is not None
                       and r["_engine_time_myr"] < time and r.get(f"_{name}") is not None]
            start_time, start_value = (earlier[-1]["_engine_time_myr"], earlier[-1][f"_{name}"]) if earlier else (0.0, 0.0)
            if time > start_time:
                row[name] = (value - start_value) / surface / (time - start_time)
    for row in rows:
        for key in [k for k in row if k.startswith("_")]:
            del row[key]


def series(bundle: Path, interior_cells: int) -> list[dict]:
    rows = []
    checkpoints = bundle.parent / f"{bundle.name}-checkpoints"
    if checkpoints.is_dir():
        for directory in sorted(checkpoints.iterdir()):
            match = re.search(r"-t(\d+(?:\.\d+)?)$", directory.name)
            if match and (directory / "fields").is_dir():
                rows.append({"time_myr": float(match.group(1)), **measure(directory, interior_cells)})
    final = {"time_myr": "final", **measure(bundle, interior_cells)}
    rows.append(final)
    add_area_rates(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("bundles", nargs="+", type=Path)
    parser.add_argument("--interior-cells", type=int, default=1)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--columns", default="time_myr,col_area,c20_area,col_thick_km,c20_thick_km,ghost,land,"
                        "cont_land_median,interior_median,over50_share,f")
    args = parser.parse_args()
    result = {}
    columns = args.columns.split(",")
    for bundle in args.bundles:
        rows = series(bundle, args.interior_cells)
        result[str(bundle)] = rows
        print(f"\n{bundle}")
        print(" ".join(f"{c:>14}" for c in columns))
        for row in rows:
            cells = []
            for c in columns:
                value = row.get(c)
                cells.append(f"{value:>14.4g}" if isinstance(value, float) else f"{str(value):>14}")
            print(" ".join(cells))
    if args.json:
        args.json.write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
