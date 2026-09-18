"""Reference product coverage and seasonal reductions for benchmark comparisons."""

import importlib.util
import json
from pathlib import Path

import numpy as np

from cell_grid import conservative_remap
from reference_map_rendering import (
    divergence, lic_luminance, lic_preview, particle_intensity, particle_preview,
    scalar_rgb, vector_preview,
)

SEASONS = ((0, "jan"), (3, "apr"), (6, "jul"), (9, "oct"))
MONTH_DAYS = np.array((31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))


def quarter_mean(months, month):
    """Require all three months: a partial quarter is not a quarter mean."""
    return np.average(months[month:month+3], axis=0, weights=MONTH_DAYS[month:month+3])


def load_diagnostics(directory, load_bundle):
    variables = ("u10m", "v10m", "u500", "v500", "u850", "v850", "tcwv",
                 "slp_anom", "w500_ascent", "pr", "viwve", "viwvn")
    fields = {}
    sources = {}
    for variable in variables:
        path = directory / f"era5_{variable}_monthly.uwclim"
        if path.exists():
            values = load_bundle(path, variable)
            if values.shape[0] != 12:
                raise ValueError(f"Monthly diagnostic must have twelve layers: {path}")
            fields[variable] = values
            sources[variable] = path.as_posix()
    for east, north in (("u10m", "v10m"), ("u500", "v500"), ("u850", "v850"), ("viwve", "viwvn")):
        if east in fields and north in fields and fields[east].shape != fields[north].shape:
            raise ValueError(f"Vector component grids differ: {east}, {north}")
    return fields, sources


def paired_remap(east, north, width):
    valid = np.isfinite(east) & np.isfinite(north)
    return (conservative_remap(np.where(valid, east, np.nan), width),
            conservative_remap(np.where(valid, north, np.nan), width))


def export_additional(directory, width, land, load_bundle, emit, families=None):
    """Consume ingestion receipts, one bundle at a time to bound memory."""
    written = set()
    for metadata in sorted(directory.glob("*-additional-preparation.json")):
        receipt = json.loads(metadata.read_text(encoding="utf-8"))
        if families and receipt.get("family") not in families:
            continue
        for record in receipt["products"]:
            if record["variable"] in ("viwve", "viwvn"):
                continue  # Existing vector exporter handles these together.
            name = record["bundle"].removesuffix("_monthly.uwclim")
            if Path(record["bundle"]).name != record["bundle"]:
                raise ValueError(f"Bundle must be a filename: {metadata}")
            bundle = directory / record["bundle"]
            values = load_bundle(bundle, record["variable"])
            if values.shape != (12, receipt["height"], receipt["width"]):
                raise ValueError(f"Bundle dimensions do not match ingestion receipt: {bundle}")
            for month, label in SEASONS:
                reduced = conservative_remap(quarter_mean(values, month), width)
                if record["ocean_only"]:
                    reduced[land] = np.nan
                if record.get("land_only", False):
                    reduced[~land] = np.nan
                emit(f"{label}_{name}", reduced,
                     scalar_rgb(reduced, record["limit"], record["palette"]), record["units"],
                     description=record.get("description", name), source_depth_m=record.get("source_depth_m"),
                     sources=record["sources"], period=record["period"] + f"; quarter starting {label}, day-weighted")
            written.add(name)
    if {"glorys_uo", "glorys_vo"}.issubset(written):
        east = load_bundle(directory / "glorys_uo_monthly.uwclim", "uo")
        north = load_bundle(directory / "glorys_vo_monthly.uwclim", "vo")
        if east.shape != north.shape:
            raise ValueError("Ocean-current component grids differ")
        for month, label in SEASONS:
            u, v = paired_remap(quarter_mean(east, month), quarter_mean(north, month), width)
            u[land] = v[land] = np.nan
            speed = np.hypot(u, v)
            emit(f"{label}_ocean_current_speed", speed,
                 vector_preview(scalar_rgb(speed, 1), u, v), "m/s",
                 sources=[str(directory / f"glorys_{key}_monthly.uwclim") for key in ("uo", "vo")],
                 period=f"three months starting {label}, day-weighted",
                 comparison="Speed of mean GLORYS current vector at shallowest model depth, not mean instantaneous speed")
        written.add("ocean_currents")
    return written


def export_diagnostics(fields, sources, width, land, emit, additional=()):
    """emit(name, values, image, units, **metadata) writes each product family."""
    unavailable = []
    for layer, east_key, north_key, seed in (
        ("surface", "u10m", "v10m", 0x1185B11),
        ("upper", "u500", "v500", 0x1185B91),
        ("850hpa", "u850", "v850", 0x1185C11),
    ):
        if east_key not in fields or north_key not in fields:
            unavailable.append({"product": f"{layer}_wind", "reason": f"Requires {east_key} and {north_key} monthly bundles"})
            continue
        for season, (month, label) in enumerate(SEASONS):
            east, north = paired_remap(fields[east_key][month], fields[north_key][month], width)
            speed = np.hypot(east, north)
            info = dict(sources=[sources[east_key], sources[north_key]], period=f"{label} monthly climatology",
                        comparison="Speed of monthly mean vector; not mean instantaneous speed")
            prefix = f"{label}_{layer}_wind"
            limit = {"surface": 25, "850hpa": 30, "upper": 40}[layer]
            for component, values in (("east", east), ("north", north)):
                emit(f"{prefix}_{component}", values, scalar_rgb(values, limit, "divergence"), "m/s", **info)
            emit(f"{prefix}_speed", speed, scalar_rgb(speed, limit), "m/s", **info)
            # Unlike the interactive benchmark renderer, reference LIC is also
            # produced at 2048 columns. Process one field at a time to bound RAM.
            luminance = lic_luminance(east, north, seed + season * 17)
            emit(f"{prefix}_lic", luminance, lic_preview(east, north, land, luminance, limit),
                 "dimensionless luminance [0,1]", **info,
                 renderer="NumPy RK2 LIC v1; same benchmark settings, independent MT19937 noise sequence",
                 numeric_meaning="LIC texture before colour/arrows; use wind speed/east/north TIFFs for physical comparisons")
            intensity = particle_intensity(east, north, seed + season * 31)
            emit(f"{prefix}_particles", intensity, particle_preview(east, north, land, intensity, limit),
                 "dimensionless trail intensity [0,1]", **info,
                 renderer="NumPy RK2 particle v1; independent noise sequence and rasterization",
                 numeric_meaning="Trail intensity before colour/arrows; not speed or particle concentration")
            # Divergence is a quarter process diagnostic, separate from the
            # representative-month flow visualization.
            qe, qn = paired_remap(quarter_mean(fields[east_key], month),
                                  quarter_mean(fields[north_key], month), width)
            div = divergence(qe, qn) * 86400.0
            emit(f"{label}_{layer}_divergence", div, scalar_rgb(div, 1, "divergence"), "day^-1",
                 sources=info["sources"], period=f"three months starting {label}, day-weighted",
                 comparison="Finite-volume divergence of monthly climatological winds on target grid")

    for variable, name, unit, limit, palette in (
        ("tcwv", "column_water", "kg/m^2", 65, "speed"),
        ("slp_anom", "pressure_anomaly", "hPa", 30, "divergence"),
        ("w500_ascent", "ascent", "hPa/day", 100, "convergence"),
        ("pr", "era5_precipitation", "mm/day", 20, "speed"),
    ):
        if variable not in fields:
            unavailable.append({"product": name, "reason": f"Requires {variable} monthly bundle"})
            continue
        for month, label in SEASONS:
            values = fields[variable]
            if variable == "pr":
                values = values / MONTH_DAYS[:, None, None]
            reduced = conservative_remap(quarter_mean(values, month), width)
            emit(f"{label}_{name}", reduced, scalar_rgb(reduced, limit, palette), unit,
                 sources=[sources[variable]], period=f"three months starting {label}, day-weighted")

    if "viwve" in fields and "viwvn" in fields:
        for month, label in SEASONS:
            east, north = paired_remap(quarter_mean(fields["viwve"], month),
                                      quarter_mean(fields["viwvn"], month), width)
            magnitude = np.hypot(east, north)
            info = dict(sources=[sources["viwve"], sources["viwvn"]],
                        period=f"three months starting {label}, day-weighted",
                        comparison="Vertically integrated water-vapour flux; not column water times mean wind")
            for component, values in (("east", east), ("north", north)):
                emit(f"{label}_column_moisture_flux_{component}", values,
                     scalar_rgb(values, 800, "divergence"), "kg m^-1 s^-1", **info)
            emit(f"{label}_column_moisture_flux", magnitude,
                 vector_preview(scalar_rgb(magnitude, 800), east, north), "kg m^-1 s^-1", **info)
            convergence = -divergence(east, north) * 86400.0
            emit(f"{label}_moisture_flux_convergence", convergence,
                 scalar_rgb(convergence, 20, "convergence"), "mm/day", **info,
                 numerical_method="Negative finite-volume divergence on target grid; benchmark uses time-integrated native face fluxes")
    else:
        unavailable.append({"product": "column_moisture_flux, moisture_flux_convergence",
                            "reason": "Requires era5_viwve_monthly.uwclim and era5_viwvn_monthly.uwclim, east/north integrated vapour flux in kg m^-1 s^-1; tcwv and winds alone are insufficient"})
    for product, reason in (
        ("boundary_moisture_flux, free_moisture_flux", "No layer-resolved humidity flux climatology; a column flux cannot determine the two layers"),
        ("surface_wind_consistency, upper_wind_consistency", "Monthly mean components cannot recover climate-sampled mean speed or directional consistency"),
        ("column_heating", "No matching radiative, sensible and phase-change heating reference"),
        ("surface_wind_vector_error, upper_wind_vector_error", "Run-dependent difference; compare benchmark vectors against exported reference east/north TIFFs"),
    ):
        unavailable.append({"product": product, "reason": reason})
    if not {"oisst_sst", "era5_sst"}.intersection(additional):
        unavailable.append({"product": "sea_surface_temperature", "reason": "No ingested SST reference"})
    if "ocean_currents" not in additional:
        unavailable.append({"product": "ocean_currents", "reason": "No ocean-current reference bundles in the retained inputs"})
    return unavailable


def export_koppen(temperature, precipitation, width, land, emit):
    # Reuse the observed control's classifier and category colours.
    path = Path(__file__).resolve().parents[1] / "benchmarks/benchmark-observed-koppen.py"
    spec = importlib.util.spec_from_file_location("observed_koppen", path)
    classifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(classifier)
    tavg = np.stack([conservative_remap(month, width) for month in temperature])
    rain = np.stack([conservative_remap(month, width) for month in precipitation])
    classes = classifier.classify_monthly(tavg, rain)
    classes[~land] = 0
    colours = np.vstack((np.zeros((1, 3), dtype=np.uint8), classifier.CLIMATE_COLOURS))
    values = classes.astype(np.float32)
    values[land & (classes == 0)] = np.nan
    emit("koppen", values, colours[classes], "category ID",
         comparison="Derived observed 12-month classification, not the independently drawn atlas",
         categories={str(i): code for i, code in enumerate(("ocean", *classifier.CLIMATE_CODES))})
