"""Reference product coverage and seasonal reductions for benchmark comparisons."""

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools.cell_grid import conservative_remap
from uwtools import palettes, references
from reference_map_rendering import (
    divergence, display_rgb, lic_luminance, lic_preview, particle_intensity, particle_preview,
    vector_preview,
)
from prepare_additional_references import display as catalogued_display

SEASONS = ((0, "jan"), (3, "apr"), (6, "jul"), (9, "oct"))
MONTH_DAYS = np.array((31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

# Display of each product: (palette, limit, factor) for reference_map_rendering.display_rgb.
# A quantity with a generated map uses its project palette ("climate:<scale>", range from
# configs/climate-palettes.json); signed components use the signed ramp over the same
# scale's maximum. Ascent (w500) and wind divergence are not the generated moisture
# convergence: they keep fixed limits on the signed ramp.
WIND_SCALES = {"surface": "surface_wind", "850hpa": "850_wind", "upper": "500_wind"}
SCALARS = (
    ("tcwv", "column_water", "kg/m^2", ("climate:water", None, 1)),
    ("slp_anom", "pressure_anomaly", "hPa", ("climate:pressure", None, 1)),
    ("w500_ascent", "ascent", "hPa/day", ("convergence", 100, 1)),
    ("pr", "era5_precipitation", "mm/day", ("climate:rain", None, 1)),
)
DIVERGENCE = ("divergence", 1, 1)


def component_display(scale):
    return ("divergence", palettes.climate_scale(scale)[1], 1)


def product_display(key):
    """Display of an exported diagnostic product, by name without its season prefix."""
    if "_wind_" in key:
        layer, style = key.split("_wind_", 1)
        if style in ("east", "north"):
            return component_display(WIND_SCALES[layer])
        return ("climate:" + WIND_SCALES[layer], None, 1)  # speed; LIC and particle backgrounds
    if key.endswith("_divergence"):
        return DIVERGENCE
    displays = {name: display for _, name, _, display in SCALARS}
    displays.update(column_moisture_flux=("climate:moisture_flux", None, 1),
                    column_moisture_flux_east=component_display("moisture_flux"),
                    column_moisture_flux_north=component_display("moisture_flux"),
                    moisture_flux_convergence=("climate:convergence", None, 1),
                    ocean_current_speed=("climate:current", None, 1))
    return displays[key]


def additional_display(dataset, record):
    """Display of an ingested field: the current catalogue's, else what its receipt recorded."""
    return (catalogued_display(dataset, record["variable"])
            or (record["palette"], record["limit"], record.get("display_factor", 1)))


def recorded(display):
    """Product metadata naming the display its PNG was drawn with; the map index reads it."""
    palette, limit, factor = display
    return dict(display=dict(palette=palette, limit=limit, factor=factor))


def quarter_mean(months, month):
    """Require all three months: a partial quarter is not a quarter mean."""
    return np.average(months[month:month+3], axis=0, weights=MONTH_DAYS[month:month+3])


def load_diagnostics(directory, load_bundle):
    """ERA5 diagnostic fields under the processed reference root `directory`."""
    variables = ("u10m", "v10m", "u500", "v500", "u850", "v850", "tcwv",
                 "slp_anom", "w500_ascent", "pr", "viwve", "viwvn")
    fields = {}
    sources = {}
    for variable in variables:
        path = references.path(f"era5_{variable}", root=directory)
        if path.exists():
            values = load_bundle(path, variable)
            if values.shape[0] != 12:
                raise ValueError(f"Monthly diagnostic must have twelve layers: {path}")
            fields[variable] = values
            sources[variable] = references.display(path)
    for east, north in (("u10m", "v10m"), ("u500", "v500"), ("u850", "v850"), ("viwve", "viwvn")):
        if east in fields and north in fields and fields[east].shape != fields[north].shape:
            raise ValueError(f"Vector component grids differ: {east}, {north}")
    return fields, sources


def paired_remap(east, north, width):
    valid = np.isfinite(east) & np.isfinite(north)
    return (conservative_remap(np.where(valid, east, np.nan), width),
            conservative_remap(np.where(valid, north, np.nan), width))


def export_additional(directory, width, land, load_bundle, emit, families=None):
    """Consume ingestion receipts (in their dataset directories under the processed
    reference root `directory`), one field at a time to bound memory."""
    written = set()
    for metadata in sorted(directory.glob("*/*/*-additional-preparation.json")):
        dataset = metadata.parent.parent.name
        receipt = json.loads(metadata.read_text(encoding="utf-8"))
        if families and receipt.get("family") not in families:
            continue
        for record in receipt["products"]:
            if record["variable"] in ("viwve", "viwvn"):
                continue  # Existing vector exporter handles these together.
            name = f"{dataset}_{record['variable']}"
            bundle = references.path(name, width=receipt["width"], root=directory)
            values = load_bundle(bundle, record["variable"])
            if values.shape != (12, receipt["height"], receipt["width"]):
                raise ValueError(f"Bundle dimensions do not match ingestion receipt: {bundle}")
            for month, label in SEASONS:
                reduced = conservative_remap(quarter_mean(values, month), width)
                if record["ocean_only"]:
                    reduced[land] = np.nan
                if record.get("land_only", False):
                    reduced[~land] = np.nan
                display = additional_display(dataset, record)
                emit(f"{label}_{name}", reduced, display_rgb(reduced, *display), record["units"], **recorded(display),
                     description=record.get("description", name), source_depth_m=record.get("source_depth_m"),
                     sources=record["sources"], period=record["period"] + f"; quarter starting {label}, day-weighted")
            written.add(name)
    if {"glorys_uo", "glorys_vo"}.issubset(written):
        east = load_bundle(references.path("glorys_uo", root=directory), "uo")
        north = load_bundle(references.path("glorys_vo", root=directory), "vo")
        if east.shape != north.shape:
            raise ValueError("Ocean-current component grids differ")
        for month, label in SEASONS:
            u, v = paired_remap(quarter_mean(east, month), quarter_mean(north, month), width)
            u[land] = v[land] = np.nan
            speed = np.hypot(u, v)
            emit(f"{label}_ocean_current_speed", speed,
                 vector_preview(display_rgb(speed, *product_display("ocean_current_speed")), u, v), "m/s",
                 **recorded(product_display("ocean_current_speed")),
                 sources=[str(references.path(f"glorys_{key}", root=directory)) for key in ("uo", "vo")],
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
            limit = palettes.climate_scale(WIND_SCALES[layer])[1]
            for component, values in (("east", east), ("north", north)):
                display = product_display(f"{layer}_wind_{component}")
                emit(f"{prefix}_{component}", values, display_rgb(values, *display), "m/s", **info, **recorded(display))
            display = product_display(f"{layer}_wind_speed")
            emit(f"{prefix}_speed", speed, display_rgb(speed, *display), "m/s", **info, **recorded(display))
            # Unlike the interactive benchmark renderer, reference LIC is also
            # produced at 2048 columns. Process one field at a time to bound RAM.
            luminance = lic_luminance(east, north, seed + season * 17)
            emit(f"{prefix}_lic", luminance, lic_preview(east, north, land, luminance, limit),
                 "dimensionless luminance [0,1]", **info, **recorded(display),
                 renderer="NumPy RK2 LIC v1; same benchmark settings, independent MT19937 noise sequence",
                 numeric_meaning="LIC texture before colour/arrows; use wind speed/east/north TIFFs for physical comparisons")
            intensity = particle_intensity(east, north, seed + season * 31)
            emit(f"{prefix}_particles", intensity, particle_preview(east, north, land, intensity, limit),
                 "dimensionless trail intensity [0,1]", **info, **recorded(display),
                 renderer="NumPy RK2 particle v1; independent noise sequence and rasterization",
                 numeric_meaning="Trail intensity before colour/arrows; not speed or particle concentration")
            # Divergence is a quarter process diagnostic, separate from the
            # representative-month flow visualization.
            qe, qn = paired_remap(quarter_mean(fields[east_key], month),
                                  quarter_mean(fields[north_key], month), width)
            div = divergence(qe, qn) * 86400.0
            emit(f"{label}_{layer}_divergence", div, display_rgb(div, *DIVERGENCE), "day^-1", **recorded(DIVERGENCE),
                 sources=info["sources"], period=f"three months starting {label}, day-weighted",
                 comparison="Finite-volume divergence of monthly climatological winds on target grid")

    for variable, name, unit, display in SCALARS:
        if variable not in fields:
            unavailable.append({"product": name, "reason": f"Requires {variable} monthly bundle"})
            continue
        for month, label in SEASONS:
            values = fields[variable]
            if variable == "pr":
                values = values / MONTH_DAYS[:, None, None]
            reduced = conservative_remap(quarter_mean(values, month), width)
            emit(f"{label}_{name}", reduced, display_rgb(reduced, *display), unit, **recorded(display),
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
                display = product_display(f"column_moisture_flux_{component}")
                emit(f"{label}_column_moisture_flux_{component}", values,
                     display_rgb(values, *display), "kg m^-1 s^-1", **info, **recorded(display))
            display = product_display("column_moisture_flux")
            emit(f"{label}_column_moisture_flux", magnitude,
                 vector_preview(display_rgb(magnitude, *display), east, north),
                 "kg m^-1 s^-1", **info, **recorded(display))
            convergence = -divergence(east, north) * 86400.0
            emit(f"{label}_moisture_flux_convergence", convergence,
                 display_rgb(convergence, *product_display("moisture_flux_convergence")), "mm/day", **info,
                 **recorded(product_display("moisture_flux_convergence")),
                 numerical_method="Negative finite-volume divergence on target grid; benchmark uses time-integrated native face fluxes")
    else:
        unavailable.append({"product": "column_moisture_flux, moisture_flux_convergence",
                            "reason": "Requires era5.viwve_monthly and era5.viwvn_monthly, east/north integrated vapour flux in kg m^-1 s^-1; tcwv and winds alone are insufficient"})
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
    path = Path(__file__).resolve().parents[1] / "benchmarks/benchmark_observed_koppen.py"
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
