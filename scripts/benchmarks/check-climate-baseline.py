# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Check a completed climate run and score seasonal ocean fields without rerunning it."""

import argparse
import json
import re
import struct
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/benchmarks"))
from climate_map_geometry import expected_map_size, load_wind_geometry
sys.path.insert(0, str(ROOT / "scripts/refs"))
from cell_grid import conservative_remap, latitude_edges


def score(simulated, reference, mask):
    weights = np.broadcast_to(-np.diff(np.sin(np.deg2rad(latitude_edges(simulated.shape[0]))))[:, None], simulated.shape)
    valid = mask & np.isfinite(reference) & np.isfinite(simulated)
    if not valid.any():
        raise ValueError("No overlapping ocean reference cells")
    a, b, w = simulated[valid], reference[valid], weights[valid]
    w = w / w.sum()
    am, bm = np.sum(w * a), np.sum(w * b)
    denominator = np.sqrt(np.sum(w * (a - am)**2) * np.sum(w * (b - bm)**2))
    return dict(cells=int(valid.sum()), simulated_mean=float(am), reference_mean=float(bm),
                bias=float(am - bm), rmse=float(np.sqrt(np.sum(w * (a - b)**2))),
                correlation=float(np.sum(w * (a - am) * (b - bm)) / denominator) if denominator else None)


def reference_quarter(name, season, width):
    path = ROOT / "refs/processed/climate" / f"{name}_monthly.uwclim"
    with path.open("rb") as stream:
        header = stream.read(40)
    version, columns, rows, layers = struct.unpack("<IIII", header[8:24])
    if header[:8] != b"UWCLIM1\0" or version != 1 or layers != 12 or columns != rows * 2:
        raise ValueError(f"Invalid monthly reference header: {path}")
    if path.stat().st_size != 40 + columns * rows * layers * 4:
        raise ValueError(f"Invalid monthly reference payload: {path}")
    values = np.memmap(path, dtype="<f4", mode="r", offset=40, shape=(layers, rows, columns))
    days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    months = range(season * 3, season * 3 + 3)
    fields = np.stack([conservative_remap(values[m], width) for m in months])
    # Keep incomplete reference coverage missing rather than extending it.
    return np.average(fields, axis=0, weights=[days[m] for m in months])


def tropical_structure(directory):
    data = np.genfromtxt(directory / "maps/weather_statistics.csv", names=True, delimiter=",")
    columns = len(np.unique(data["longitude"]))
    rows = columns // 2
    latitude = 90 - (np.arange(rows) + 0.5) * 180 / rows
    selected = np.abs(latitude) < 20
    weights = np.broadcast_to(np.cos(np.deg2rad(latitude[selected]))[:, None], (selected.sum(), columns))
    result = {}
    for season in (0, 2):
        state = data[(data["season"] == season) & (data["layer"] == 0)]
        simulated = [state["mean_u_mps"].reshape(rows, columns), -state["mean_v_south_mps"].reshape(rows, columns)]
        reference = [reference_quarter(name, season, columns) for name in ("era5_u10m", "era5_v10m")]
        entry = {}
        for label, fields in (("simulation", simulated), ("ERA5", reference)):
            u, v = [field[selected] for field in fields]
            full = np.sum(weights * (u*u + v*v))
            regional = np.sum(weights * ((u-u.mean(axis=1, keepdims=True))**2 + (v-v.mean(axis=1, keepdims=True))**2))
            entry[label] = dict(vector_rms_mps=float(np.sqrt(full / weights.sum())),
                regional_vector_rms_mps=float(np.sqrt(regional / weights.sum())),
                zonal_fraction_of_squared_wind=float(1-regional/full) if full else 0.0)
        error = sum((a[selected]-b[selected])**2 for a, b in zip(simulated, reference))
        entry["vector_rmse_mps"] = float(np.sqrt(np.sum(weights*error) / weights.sum()))
        result[str(season+1)] = entry
    return result


def compare_reference_skill(current, baseline):
    """Independent targets: improvements in wind cannot cancel worse rainfall."""
    comparisons = {}
    def add(label, a, b, higher_is_better):
        if a is None or b is None or not np.isfinite([a, b]).all():
            raise ValueError(f"Missing reference metric: {label}")
        delta = a - b
        comparisons[label] = dict(current=a, baseline=b, change=delta,
            higher_is_better=higher_is_better,
            regressed=bool(delta < -1e-6 if higher_is_better else delta > 1e-6))

    for key in ("area_weighted_land_precipitation_correlation", "era5_column_water_correlation",
                "era5_surface_eastward_wind_correlation", "era5_surface_northward_wind_correlation",
                "era5_global_pressure_correlation"):
        add(key, current["historical_metrics"][str(current["run_id"])].get(key),
            baseline["historical_metrics"][str(baseline["run_id"])].get(key), True)
    for season in ("1", "3"):
        add(f"tropical_wind_vector_rmse_mps_quarter_{season}",
            current["tropical_structure"][season]["vector_rmse_mps"],
            baseline["tropical_structure"][season]["vector_rmse_mps"], False)
    for key, value in current["ocean_reference_scores"].items():
        add(f"{key}_rmse", value["rmse"], baseline["ocean_reference_scores"][key]["rmse"], False)
    return dict(baseline_run_id=baseline["run_id"],
        status="regressed" if any(v["regressed"] for v in comparisons.values()) else "no_detected_regression",
        criteria="Each listed reference metric must not deteriorate beyond 1e-6 numerical tolerance. This is a regression gate, not proof of convergence or Earth accuracy.",
        metrics=comparisons)


def check(run_id, compare_run_id=None):
    records = json.loads((ROOT / "runs/registry/climate.json").read_text())["runs"]
    record = next(r for r in records if r["id"] == run_id)
    width = int(record["horizontal_resolution"])
    directory = ROOT / f"runs/diagnostics/climate/{run_id}"
    log = (ROOT / f"runs/logs/climate/{run_id}.log").read_text(encoding="utf-8", errors="replace")
    maps = {}
    wind_geometry = load_wind_geometry(directory, width)
    for path in (ROOT / "runs/maps/climate").rglob(f"{run_id}.png"):
        with Image.open(path) as image:
            image.load()
            if image.size != expected_map_size(path.relative_to(ROOT / "runs/maps/climate"), width, wind_geometry):
                raise ValueError(f"Wrong exported dimensions: {path}")
            maps[path.relative_to(ROOT).as_posix()] = list(image.size)
    if not maps:
        raise ValueError("No map exports")
    with Image.open(ROOT / f"runs/fields/climate/rain/annual/s/{run_id}.tif") as image:
        rain = np.asarray(image)
    if rain.shape != (width // 2, width) or not np.isfinite(rain).all() or (rain < 0).any():
        raise ValueError("Invalid rainfall raster")
    numeric = {}
    for name in ("climate_process_fields", "weather_statistics", "ocean_fields"):
        data = np.genfromtxt(directory / f"maps/{name}.csv", names=True, delimiter=",")
        if not all(np.isfinite(data[field]).all() for field in data.dtype.names):
            raise ValueError(f"Non-finite {name}")
        numeric[name] = len(data)
    mixed_layer = "tropical closure mixed-layer" in record.get("information", "")
    if mixed_layer:
        data = np.genfromtxt(directory / "maps/tropical_boundary_layer.csv", names=True, delimiter=",")
        if data.size == 0 or not all(np.isfinite(data[field]).all() for field in data.dtype.names):
            raise ValueError("Missing or non-finite tropical boundary-layer fields")
        if set(data["season"]) != {0, 1, 2, 3} or (data["solver_residual"] > 1.01e-4).any():
            raise ValueError("Incomplete or unaccepted inversion solve")
        diagnoses = [line for line in log.splitlines() if line.startswith("Tropical boundary layer season=")]
        if not diagnoses or any("accepted=1" not in line for line in diagnoses):
            raise ValueError("A tropical boundary-layer update was rejected")
        if max(float(re.search(r"momentum_residual_mps2=([0-9eE.+-]+)", line)[1]) for line in diagnoses) > 1e-7:
            raise ValueError("Tropical momentum balance failed")
        numeric["tropical_boundary_layer"] = len(data)
    years = [dict(re.findall(r"(\w+)=([0-9eE.+-]+)", line)) for line in log.splitlines()
             if line.startswith("Climate hydrology spinup cycle=")]
    profile = {key: sum(float(year[key]) for year in years) for key in
               ("elapsed_seconds", "transport_seconds", "weather_seconds", "transport_calls", "transport_cell_substeps")}
    profile["years"] = len(years)
    profile["maximum_substeps"] = max(int(year["transport_max_substeps"]) for year in years)
    water_residual = max(abs(float(value)) for value in re.findall(r"(?<!\w)relative_residual=([0-9eE.+-]+)", log))
    if water_residual > 2e-6:
        raise ValueError(f"Water budget does not close: {water_residual}")

    ocean = np.genfromtxt(directory / "maps/ocean_fields.csv", names=True, delimiter=",")
    columns = len(np.unique(ocean["longitude"]))
    ocean_scores = {}
    if "mean_sst_c" in ocean.dtype.names:
        for season in (0, 2):
            state = ocean[ocean["season"] == season]
            with Image.open(ROOT / f"runs/maps/climate/sea_temp/{season + 1}/s/{run_id}.png") as image:
                # The native SST renderer masks land black; 0 C water is not black.
                ocean_mask = np.any(np.asarray(image.convert("RGB")) != 0, axis=-1)
            mask = conservative_remap(ocean_mask.astype(np.float32), columns) > 0.5
            for column, bundle, sign in (("mean_sst_c", "oisst_sst", 1),
                                         ("east_current_mps", "glorys_uo", 1),
                                         ("south_current_mps", "glorys_vo", -1)):
                simulated = sign * state[column].reshape(columns // 2, columns)
                reference = reference_quarter(bundle, season, columns)
                label = column.replace("south_current", "north_current")
                ocean_scores[f"{label}_quarter_{season + 1}"] = score(simulated, reference, mask)
    result = dict(run_id=run_id, maps=maps, numeric_csv_rows=numeric, profile=profile,
                  maximum_water_relative_residual=water_residual,
                  rainfall_range_mm=[float(rain.min()), float(rain.max())],
                  final_coupling=next(line for line in reversed(log.splitlines()) if line.startswith("Climate coupling iteration=")),
                  ocean_reference_scores=ocean_scores,
                  tropical_structure=tropical_structure(directory),
                  tropical_comparison="20 S to 20 N, both land and ocean, area weighted; day-weighted Jan-Mar/Jul-Sep ERA5 means on the native weather-statistics grid.",
                  ocean_comparison="Native ocean grid; day-weighted Jan-Mar/Jul-Sep reference means; majority ocean mask. GLORYS shallowest currents versus model mixed-layer currents; OISST versus mean liquid SST.",
                  historical_metrics={str(r["id"]): r.get("metrics", {}) for r in records
                      if r["id"] in (142, 143, 144, 145, compare_run_id, run_id)},
                  historical_comparison={
                      "role": "Informational only; reference_regression gates only the requested previous run.",
                      "limitations": "Runs 142–145 have identical recorded metrics. Their 512×257 exports, seeds, reference inputs and workloads differ from current runs. A metric leader is not a universally best simulation; zero spatial scores with no compared cells are unavailable.",
                      "runs": {str(r["id"]): {key: r.get(key) for key in
                          ("horizontal_resolution", "datetime", "information", "diagnostics_directory")}
                          for r in records if r["id"] in (142, 143, 144, 145, compare_run_id, run_id)}})
    if compare_run_id is not None:
        baseline = json.loads((ROOT / f"runs/reports/{compare_run_id}-checks.json").read_text())
        if baseline["run_id"] != compare_run_id:
            raise ValueError("Baseline check file has a different run ID")
        result["reference_regression"] = compare_reference_skill(result, baseline)
    destination = ROOT / f"runs/reports/{run_id}-checks.json"
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Run {run_id}: {len(maps)} correctly sized maps; finite fields; water residual {water_residual:.3g}")
    print(json.dumps(dict(profile=profile, ocean_reference_scores=ocean_scores), indent=2))
    print(destination)
    if compare_run_id is not None:
        comparison = result["reference_regression"]
        failed = [key for key, value in comparison["metrics"].items() if value["regressed"]]
        print(f"Reference gate versus run {compare_run_id}: {comparison['status']}; regressed targets: {', '.join(failed) or 'none'}")
        return 2 if failed else 0
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--compare-run-id", type=int,
                        help="Compare saved reference checks; exit 2 if any target regresses (does not rerun the simulation).")
    args = parser.parse_args()
    sys.exit(check(args.run_id, args.compare_run_id))
