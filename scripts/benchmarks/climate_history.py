"""Recorded historical skill with explicit availability and original grids."""

import csv
import math
import struct
from pathlib import Path

HISTORICAL_NOTE = ("Historical context, not a controlled regression gate: seeds, grids, reference inputs "
                   "and model workloads differ. Original PNG grids are retained; no spatial rescaling or "
                   "independent spatial Koeppen validation is implied.")
DRY_KEY = "area_weighted_land_below_1mm_fraction"
DRY_LABEL = "Land below 12 mm/year (%) [stored monthly mean <1 mm]"
METRICS = [
    ("area_weighted_land_precipitation_correlation", "Land rainfall correlation", True,
     "annual_imerg_precipitation_comparison", None, "land"),
    ("area_weighted_land_precipitation_correlation_5degree", "Land rainfall correlation, 5 degree", True,
     "annual_imerg_precipitation_comparison", None, "land"),
    ("land_annual_temperature_rmse_c", "Land annual temperature RMSE (C)", False,
     "monthly_climate_reference_comparison", "temperature_c", None),
]
for key, label, variable, higher in (
    ("era5_column_water_correlation", "Column water correlation", "column_water_kg_m2", True),
    ("era5_surface_eastward_wind_correlation", "Surface eastward wind correlation", "surface_eastward_wind_m_s", True),
    ("era5_surface_northward_wind_correlation", "Surface northward wind correlation", "surface_northward_wind_m_s", True),
    ("era5_transport_eastward_wind_correlation", "Upper eastward wind correlation", "upper_500hpa_eastward_wind_m_s", True),
    ("era5_transport_northward_wind_correlation", "Upper northward wind correlation", "upper_500hpa_northward_wind_m_s", True),
    ("era5_global_pressure_correlation", "Pressure correlation", "sea_level_pressure_anomaly_hpa", True),
    ("era5_global_pressure_rmse_hpa", "Pressure RMSE (hPa)", "sea_level_pressure_anomaly_hpa", False),
    ("era5_vertical_ascent_correlation", "Ascent correlation", "midlevel_ascent_hpa_day", True),
    ("era5_vertical_ascent_rmse_hpa_day", "Ascent RMSE (hPa/day)", "midlevel_ascent_hpa_day", False),
    ("era5_land_precipitation_correlation", "ERA5 land rainfall correlation", "precipitation_mm_month", True),
):
    METRICS.append((key, label, higher, "monthly_physical_reference_comparison", variable,
                    "land" if key == "era5_land_precipitation_correlation" else "global"))


def png_dimensions(path):
    with Path(path).open("rb") as stream:
        header = stream.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Invalid PNG header: {path}")
    return list(struct.unpack(">II", header[16:24]))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def diagnostic_rows(root, run_id, name):
    for path in (root / f"runs/diagnostics/climate/{run_id}/{name}.csv",
                 root / f"runs/diagnostics/{name}/{run_id}.csv"):
        if path.exists():
            with path.open(newline="", encoding="utf-8-sig") as stream:
                return list(csv.DictReader(stream))
    return []


def recorded_metrics(root, record):
    """Require a positive comparison count; absent-reference zeroes are not skill."""
    diagnostics = {name: diagnostic_rows(root, record["id"], name) for name in {m[3] for m in METRICS}}
    values = {}
    for key, _, _, name, variable, scope in METRICS:
        row = next((r for r in diagnostics[name]
                    if (variable is None or r.get("variable") == variable)
                    and (scope is None or r.get("scope") == scope)
                    and r.get("period", "annual_mean") == "annual_mean"), None)
        value = None
        if row and int(row.get("compared_cells", 0)) > 0:
            value = (float(row["area_weighted_rmse"]) if key == "land_annual_temperature_rmse_c"
                     else record.get("metrics", {}).get(key))
        values[key] = value if finite(value) else None
    dry = record.get("metrics", {}).get(DRY_KEY)
    values[DRY_KEY] = dry if finite(dry) else None
    return values


def run_summary(root, record, metrics):
    path = root / f"runs/maps/climate/koppen/{record['id']}.png"
    return dict(run_id=record["id"], original_png_dimensions=png_dimensions(path) if path.exists() else None,
                recorded_horizontal_resolution=record.get("horizontal_resolution"), metrics=metrics,
                information=record.get("information", ""))


def historical_comparison(records, run_id, historical_run_id, root):
    selected = {r["id"]: r for r in records}
    for required in (run_id, historical_run_id):
        if required not in selected:
            raise ValueError(f"Historical comparison run {required} is absent from the registry")
    values = {i: recorded_metrics(root, record) for i, record in selected.items()}
    current, anchor = values[run_id], values[historical_run_id]
    comparisons, leaders = {}, {}
    for key, label, higher, *_ in METRICS:
        a, b = current[key], anchor[key]
        comparisons[key] = dict(label=label, current=a, historical=b,
            change=a-b if a is not None and b is not None else None, higher_is_better=higher)
        available = [(i, fields[key]) for i, fields in values.items() if fields[key] is not None]
        if available:
            best = (max if higher else min)(v for _, v in available)
            leaders[key] = dict(label=label, value=best, higher_is_better=higher,
                run_ids=[i for i, v in available if abs(v-best) <= 1e-10],
                eligible_runs=len(available))
    cohort = [run_summary(root, selected[i], values[i]) for i in range(142, 146) if i in selected]
    return dict(status="historical_context_only", note=HISTORICAL_NOTE,
        anchor=run_summary(root, selected[historical_run_id], anchor),
        current=run_summary(root, selected[run_id], current), metrics=comparisons,
        runs_142_145=cohort, recorded_metric_leaders=leaders,
        leader_limitations="Independent recorded-metric leaders, not an overall best model. "
            "Only finite scores with positive saved reference-comparison counts are ranked; "
            "different resolutions, references and workloads remain non-comparable. Dry fraction is descriptive, not ranked.",
        missing_metrics="Unavailable values remain null; no matching reference cells is not zero error.",
        dry_fraction_label=DRY_LABEL,
        recovered_checkpoint=(dict(commit="10b4aba9d380e379ab99407a73f89eac3992e94c",
            scope="Retained source and registry through run 145; the 142-145 registry metrics match.",
            provenance="Recoverable source checkpoint only. Saved run manifests lack source/binary hashes; "
                       "the exact clean build used by each executable is unverified.")
            if historical_run_id in range(142, 146) else None))
