#!/usr/bin/env python3
"""Create Earth benchmark import maps at arbitrary even horizontal resolutions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from PIL import Image

from cell_grid import conservative_remap, grid_dimensions, latitude_centres
from reference_climate_maps import export_additional, export_diagnostics, export_koppen, load_diagnostics
from reference_map_rendering import write_float32_geotiff
from reference_map_layout import product_path
from climate_palettes import VERSION as PALETTE_VERSION

DEFAULT_WIDTHS = (128, 256, 512, 1024, 2048)
UWCLIM_MAGIC = b"UWCLIM1\0"


def load_single_channel(path: Path, dtype: np.dtype) -> np.ndarray:
    with Image.open(path) as image:
        values = np.asarray(image, dtype=dtype)
    if values.ndim != 2:
        raise ValueError(f"Expected a single-channel raster: {path}")
    return values


def prepare_topography(
    elevation: np.ndarray,
    land_mask: np.ndarray,
    width: int,
) -> tuple[np.ndarray, np.ndarray]:
    if elevation.shape != land_mask.shape:
        raise ValueError(
            f"Elevation and land-mask dimensions differ: "
            f"{elevation.shape} != {land_mask.shape}"
        )
    if elevation.shape[1] != elevation.shape[0] * 2:
        raise ValueError(
            "Expected a 2:1 cell-centred Earth source without stored pole rows"
        )

    _, height = grid_dimensions(width)
    land = land_mask.astype(np.float32)
    ocean = 1.0 - land
    land_weight = conservative_remap(land, width, height)
    ocean_weight = conservative_remap(ocean, width, height)
    land_sum = conservative_remap(np.where(land_mask, elevation, 0.0), width, height)
    ocean_sum = conservative_remap(np.where(land_mask, 0.0, elevation), width, height)
    selected_land = land_weight >= ocean_weight

    land_metres = np.zeros(selected_land.shape, dtype=np.uint16)
    sea_metres = np.zeros(selected_land.shape, dtype=np.uint16)
    land_values = land_sum / np.maximum(land_weight, 1.0e-6)
    ocean_values = ocean_sum / np.maximum(ocean_weight, 1.0e-6)
    land_metres[selected_land] = np.clip(
        np.rint(np.maximum(land_values[selected_land], 1.0)), 1, 65535
    ).astype(np.uint16)
    sea_metres[~selected_land] = np.clip(
        np.rint(np.maximum(-ocean_values[~selected_land], 1.0)), 1, 65535
    ).astype(np.uint16)
    return land_metres, sea_metres


def write_tiff16(path: Path, values: np.ndarray) -> None:
    """Write the small uncompressed uint16 TIFF subset used by the importer."""
    height, width = values.shape
    samples = values.astype("<u2", copy=False).tobytes(order="C")
    entries = (
        (256, 4, width),
        (257, 4, height),
        (258, 3, 16),
        (259, 3, 1),
        (262, 3, 1),
        (273, 4, 8 + 2 + 11 * 12 + 4),
        (277, 3, 1),
        (278, 4, height),
        (279, 4, len(samples)),
        (284, 3, 1),
        (339, 3, 1),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as output:
        output.write(b"II")
        output.write(struct.pack("<HIH", 42, 8, len(entries)))
        for tag, kind, value in entries:
            output.write(struct.pack("<HHI", tag, kind, 1))
            if kind == 3:
                output.write(struct.pack("<H", value) + b"\0\0")
            else:
                output.write(struct.pack("<I", value))
        output.write(struct.pack("<I", 0))
        output.write(samples)


def file_record(path: Path, root: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {
        "file": path.relative_to(root).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def write_reference_product(root, width, name, values, image, units, category="climate", **metadata):
    _, height = grid_dimensions(width)
    if values.shape != (height, width):
        raise ValueError(f"Invalid product dimensions: {name}")
    png = root / product_path(name, width, "maps", category)
    tiff = root / product_path(name, width, "fields", category)
    csv_path = root / product_path(name, width, "csv", category)
    png.parent.mkdir(parents=True, exist_ok=True)
    if not isinstance(image, Image.Image):
        image = Image.fromarray(image)
    image.save(png)
    write_float32_geotiff(tiff, values)
    write_grid_csv(csv_path, values)
    return dict(name=name, units=units, palette_version=PALETTE_VERSION,
                maps=[file_record(png, root), file_record(tiff, root)],
                csv=file_record(csv_path, root),
                valid_cells=int(np.count_nonzero(np.isfinite(values))), **metadata)


def refresh_additional_maps(args, widths):
    """Add newly ingested fields without rerendering unchanged wind textures."""
    manifest_path = args.output_root / "metadata/reference-maps/manifest.json"
    receipt = json.loads(manifest_path.read_text(encoding="utf-8"))
    fields, sources = {}, {}
    for variable in ("viwve", "viwvn"):
        path = args.diagnostics / f"era5_{variable}_monthly.uwclim"
        if path.exists():
            fields[variable] = load_uwclim(path, variable)
            if fields[variable].shape[0] != 12:
                raise ValueError(f"Expected twelve months: {path}")
            sources[variable] = path.as_posix()
    if len(fields) == 2 and fields["viwve"].shape != fields["viwvn"].shape:
        raise ValueError("Moisture-flux component grids differ")
    for width in widths:
        record = next((r for r in receipt["outputs"] if r["width"] == width), None)
        if record is None:
            raise ValueError(f"Generate the base reference maps for width {width} first")
        height = width // 2
        mask_path = args.output_root / f"land_ocean/maps/earth_land_ocean_{width}x{height}.png"
        mask = load_single_channel(mask_path, np.uint8) > 0
        if mask.shape != (height, width):
            raise ValueError(f"Invalid mask dimensions: {mask_path}")
        products = []

        def emit(*values, **metadata):
            products.append(write_reference_product(args.output_root, width, *values, **metadata))

        additional = export_additional(args.diagnostics, width, mask, load_uwclim, emit,
                                       families=args.additional_family)
        if len(fields) == 2:
            export_diagnostics(fields, sources, width, mask, emit, additional)
        replacements = {p["name"]: p for p in products}
        retained = {p["name"]: p for p in record["products"]}
        retained.update(replacements)
        record["products"] = list(retained.values())
        available = set()
        if len(fields) == 2:
            available.add("column_moisture_flux, moisture_flux_convergence")
        if "ocean_currents" in additional:
            available.add("ocean_currents")
        if {"oisst_sst", "era5_sst"}.intersection(additional):
            available.add("sea_surface_temperature")
            # Migrate the earlier combined availability entry.
            for missing in record["unavailable"]:
                if missing["product"] == "sea_surface_temperature, ocean_currents":
                    missing.update(product="ocean_currents", reason="No ocean-current reference bundles in the retained inputs")
        record["unavailable"] = [p for p in record["unavailable"] if p["product"] not in available]
        print(f"Refreshed {len(products)} additional products at {width}x{height}", flush=True)
    receipt.setdefault("diagnostic_sources", {}).update(sources)
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    temporary.replace(manifest_path)


def load_precipitation(path: Path) -> np.ndarray:
    table = np.genfromtxt(
        path,
        delimiter=",",
        skip_header=1,
        dtype=np.float32,
    )
    if table.ndim != 2 or table.shape[1] < 4:
        raise ValueError(f"Invalid precipitation reference CSV: {path}")
    values = table[:, 2:]
    rows, columns = values.shape
    _, expected_rows = grid_dimensions(columns)
    if rows != expected_rows:
        raise ValueError(
            f"Expected a cell-centred {columns}x{expected_rows} reference: {path}"
        )
    if not np.allclose(table[:, 0], np.arange(rows), rtol=0.0, atol=1.0e-6):
        raise ValueError(f"Invalid row indices in precipitation reference: {path}")
    if not np.allclose(
        table[:, 1], latitude_centres(rows), rtol=0.0, atol=5.0e-6
    ):
        raise ValueError(f"Invalid latitude centres in precipitation reference: {path}")
    return values


def load_uwclim(path: Path, expected_variable: str) -> np.ndarray:
    with path.open("rb") as source:
        if source.read(8) != UWCLIM_MAGIC:
            raise ValueError(f"Invalid UWCLIM magic: {path}")
        version, width, height, layers = struct.unpack("<IIII", source.read(16))
        variable = source.read(16).rstrip(b"\0").decode("ascii")
        if version != 1 or variable != expected_variable:
            raise ValueError(
                f"Unexpected UWCLIM header in {path}: version={version}, variable={variable}"
            )
        grid_dimensions(width, height)
        values = np.frombuffer(source.read(), dtype="<f4")
    expected = width * height * layers
    if values.size != expected:
        raise ValueError(f"Invalid UWCLIM payload size in {path}: {values.size} != {expected}")
    return values.reshape((layers, height, width))


def prepare_precipitation(values: np.ndarray, width: int) -> np.ndarray:
    return conservative_remap(values, width)


def write_grid_csv(path: Path, values: np.ndarray) -> None:
    height, width = values.shape
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["y", "latitude", *(f"x{x}" for x in range(width))])
        for y, (latitude, row) in enumerate(zip(latitude_centres(height), values)):
            writer.writerow(
                [
                    y,
                    f"{latitude:.6f}",
                    *("nan" if not np.isfinite(value) else f"{value:.6f}" for value in row),
                ]
            )


def _palette(values: np.ndarray, anchors: tuple[tuple[float, tuple[int, int, int]], ...]) -> np.ndarray:
    levels = np.asarray([anchor[0] for anchor in anchors], dtype=np.float32)
    colours = np.asarray([anchor[1] for anchor in anchors], dtype=np.float32)
    rgb = np.empty((*values.shape, 3), dtype=np.uint8)
    for channel in range(3):
        rgb[..., channel] = np.clip(
            np.interp(np.where(np.isfinite(values), values, 0), levels, colours[:, channel]), 0.0, 255.0
        ).astype(np.uint8)
    rgb[~np.isfinite(values)] = 0
    return rgb


def visual_palette(quantity: str):
    from climate_palettes import anchors
    temperature = anchors(-54, 40, True)
    precipitation = anchors(0, 6000)
    land = (
        (0.0, (0, 0, 0)),
        (1.0, (66, 92, 45)),
        (1500.0, (153, 137, 94)),
        (4000.0, (193, 177, 150)),
        (6500.0, (245, 244, 239)),
    )
    sea = (
        (0.0, (0, 0, 0)),
        (1.0, (86, 155, 177)),
        (3000.0, (22, 72, 125)),
        (11000.0, (5, 28, 60)),
    )
    return {
        "temperature": temperature,
        "precipitation": precipitation,
        "land": land,
        "sea": sea,
    }[quantity]


def write_visual_map(path: Path, values: np.ndarray, quantity: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(_palette(values.astype(np.float32), visual_palette(quantity)), mode="RGB").save(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-additional", action="store_true",
                        help="Refresh only newly ingested diagnostics and moisture flux in an existing map manifest.")
    parser.add_argument(
        "--width",
        "--widths",
        nargs="+",
        type=int,
        default=DEFAULT_WIDTHS,
        help="Even benchmark widths; height is always width/2.",
    )
    parser.add_argument(
        "--elevation",
        type=Path,
        default=Path("refs/processed/heightmap/2048.tif"),
        help="Signed float Earth elevation/bathymetry master.",
    )
    parser.add_argument(
        "--land-mask",
        type=Path,
        default=Path("refs/processed/land_ocean/2048.png"),
        help="Non-zero for non-ocean cells; dimensions must match --elevation.",
    )
    parser.add_argument(
        "--precipitation",
        type=Path,
        default=Path("refs/processed/climate/imerg_prec_annual.uwclim"),
        help="Cell-centred annual IMERG precipitation UWCLIM bundle.",
    )
    parser.add_argument(
        "--temperature",
        type=Path,
        default=Path("refs/processed/climate/worldclim_tavg_monthly.uwclim"),
        help="Cell-centred monthly WorldClim temperature UWCLIM bundle.",
    )
    parser.add_argument(
        "--diagnostics",
        type=Path,
        default=Path("refs/processed/climate"),
        help="Directory of monthly ERA5 UWCLIM bundles; missing products are reported in the manifest.",
    )
    parser.add_argument(
        "--monthly-precipitation",
        type=Path,
        default=Path("refs/processed/climate/imerg_prec_monthly.uwclim"),
        help="Monthly precipitation (prec, mm/month) for observed Koppen classification.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("refs/processed"),
    )
    parser.add_argument("--additional-family", nargs="+", choices=("era5", "land", "oisst", "marine", "ceres"),
                        help="With --refresh-additional, export only these ingestion families.")
    args = parser.parse_args()
    if args.additional_family and not args.refresh_additional:
        parser.error("--additional-family requires --refresh-additional")
    return args


def main() -> None:
    args = parse_args()
    widths = list(dict.fromkeys(args.width))
    invalid = [width for width in widths if width < 2 or width % 2]
    if invalid:
        raise ValueError(f"Widths must be even and at least 2: {invalid}")
    if args.refresh_additional:
        refresh_additional_maps(args, widths)
        return

    elevation = load_single_channel(args.elevation, np.float32)
    land_mask = load_single_channel(args.land_mask, np.uint8) > 0
    precipitation_layers = load_uwclim(args.precipitation, "prec_annual")
    if precipitation_layers.shape[0] != 1:
        raise ValueError("Annual precipitation bundle must contain exactly one layer")
    precipitation = precipitation_layers[0]
    temperature_layers = load_uwclim(args.temperature, "tavg")
    if temperature_layers.shape[0] != 12:
        raise ValueError("Monthly temperature bundle must contain exactly twelve layers")
    sampled_temperature = temperature_layers[[0, 3, 6, 9]]
    valid_temperature = np.isfinite(sampled_temperature)
    temperature_count = valid_temperature.sum(axis=0)
    temperature = np.divide(
        np.where(valid_temperature, sampled_temperature, 0).sum(axis=0),
        temperature_count,
        out=np.full(sampled_temperature.shape[1:], np.nan, dtype=np.float32),
        where=temperature_count > 0,
    )
    diagnostics, diagnostic_sources = load_diagnostics(args.diagnostics, load_uwclim)
    monthly_precipitation = None
    if args.monthly_precipitation.exists():
        monthly_precipitation = load_uwclim(args.monthly_precipitation, "prec")
        if monthly_precipitation.shape[0] != 12:
            raise ValueError("Monthly precipitation bundle must contain exactly twelve layers")
    args.output_root.mkdir(parents=True, exist_ok=True)
    records = []

    for width in widths:
        _, height = grid_dimensions(width)
        land, sea = prepare_topography(elevation, land_mask, width)
        precipitation_reference = prepare_precipitation(precipitation, width)
        temperature_reference = conservative_remap(temperature, width)
        land_map_directory = args.output_root / "land_ocean" / "maps"
        land_csv_directory = args.output_root / "land_ocean" / "csv"
        climate_map_directory = args.output_root / "climate" / "maps"
        climate_csv_directory = args.output_root / "climate" / "csv"
        land_tiff = land_map_directory / f"earth_land_{width}x{height}.tif"
        sea_tiff = land_map_directory / f"earth_sea_{width}x{height}.tif"
        land_png = land_map_directory / f"earth_land_{width}x{height}.png"
        sea_png = land_map_directory / f"earth_sea_{width}x{height}.png"
        precipitation_png = args.output_root / product_path("precipitation", width, "maps")
        temperature_png = args.output_root / product_path("temperature", width, "maps")
        land_csv = land_csv_directory / f"earth_land_{width}x{height}.csv"
        sea_csv = land_csv_directory / f"earth_sea_{width}x{height}.csv"
        precipitation_csv = args.output_root / product_path("precipitation", width, "csv")
        temperature_csv = args.output_root / product_path("temperature", width, "csv")
        write_tiff16(land_tiff, land)
        write_tiff16(sea_tiff, sea)
        write_visual_map(land_png, land, "land")
        write_visual_map(sea_png, sea, "sea")
        write_visual_map(precipitation_png, precipitation_reference, "precipitation")
        write_visual_map(temperature_png, temperature_reference, "temperature")
        products = []

        def emit(name, values, image, units, category="climate", **metadata):
            products.append(write_reference_product(args.output_root, width, name, values, image,
                                                     units, category=category, **metadata))

        mask = land > 0
        emit("land_ocean", mask.astype(np.float32), mask.astype(np.uint8) * 255,
             "0=ocean, 1=continental (including inland lakes)", category="land_ocean")
        # Existing uint16 import filenames remain compatible with map_imports.
        # Every preview also receives a separate float32 scientific field.
        for name, values, png, unit, category in (
            ("land", land, land_png, "m (rounded import elevation)", "land_ocean"),
            ("sea", sea, sea_png, "m (rounded positive import depth)", "land_ocean"),
            ("precipitation", precipitation_reference, precipitation_png, "mm/year", "climate"),
            ("temperature", temperature_reference, temperature_png, "degrees C", "climate"),
        ):
            with Image.open(png) as preview:
                source = args.precipitation if name == "precipitation" else args.temperature if name == "temperature" else args.elevation
                emit(name, values, preview.copy(), unit, category=category, sources=[source.as_posix()])
        if monthly_precipitation is not None:
            export_koppen(temperature_layers, monthly_precipitation, width, mask, emit)
            products[-1]["sources"] = [args.temperature.as_posix(), args.monthly_precipitation.as_posix()]
        additional = export_additional(args.diagnostics, width, mask, load_uwclim, emit)
        unavailable = export_diagnostics(diagnostics, diagnostic_sources, width, mask, emit, additional)
        if monthly_precipitation is None:
            unavailable.append({"product": "koppen", "reason": f"Missing monthly precipitation: {args.monthly_precipitation}"})

        overlap = int(np.count_nonzero((land > 0) & (sea > 0)))
        uncovered = int(np.count_nonzero((land == 0) & (sea == 0)))
        if overlap or uncovered:
            raise RuntimeError(
                f"Invalid {width}x{height} land/sea coverage: "
                f"overlap={overlap}, uncovered={uncovered}"
            )
        records.append(
            {
                "width": width,
                "height": height,
                "products": products,
                "unavailable": unavailable,
                "latitude_first_degrees": float(latitude_centres(height)[0]),
                "latitude_last_degrees": float(latitude_centres(height)[-1]),
                "land_maps": [
                    file_record(land_tiff, args.output_root),
                    file_record(land_png, args.output_root),
                ],
                "land_csv": file_record(land_csv, args.output_root),
                "sea_maps": [
                    file_record(sea_tiff, args.output_root),
                    file_record(sea_png, args.output_root),
                ],
                "sea_csv": file_record(sea_csv, args.output_root),
                "precipitation_map": file_record(precipitation_png, args.output_root),
                "precipitation_csv": file_record(precipitation_csv, args.output_root),
                "temperature_map": file_record(temperature_png, args.output_root),
                "temperature_csv": file_record(temperature_csv, args.output_root),
                "land_cells": int(np.count_nonzero(land)),
                "sea_cells": int(np.count_nonzero(sea)),
                "precipitation_valid_cells": int(
                    np.count_nonzero(np.isfinite(precipitation_reference))
                ),
                "temperature_valid_cells": int(
                    np.count_nonzero(np.isfinite(temperature_reference))
                ),
            }
        )
        print(
            f"Created {width}x{height}: land={records[-1]['land_cells']} "
            f"sea={records[-1]['sea_cells']} "
            f"precipitation={records[-1]['precipitation_valid_cells']} "
            f"temperature={records[-1]['temperature_valid_cells']}",
            flush=True,
        )
        print(f"  {len(products)} PNG/float32 TIFF/CSV products; {len(unavailable)} unavailable groups (see manifest)", flush=True)

    receipt = {
        "map_layout_version": 1,
        "elevation_source": args.elevation.as_posix(),
        "land_mask_source": args.land_mask.as_posix(),
        "precipitation_source": args.precipitation.as_posix(),
        "temperature_source": args.temperature.as_posix(),
        "diagnostic_sources": diagnostic_sources,
        "grid": "cell-centred global equirectangular W x W/2",
        "latitude_formula": "90 - (y + 0.5) * 180 / height",
        "transformations": {
            "topography": (
                "spherical-area conservative averages computed separately over source "
                "land and ocean coverage; majority coverage selects exactly one class"
            ),
            "precipitation": "masked spherical-area conservative averaging; units remain mm/year",
            "temperature": (
                "Jan/Apr/Jul/Oct WorldClim means followed by masked spherical-area "
                "conservative averaging; units remain degrees C"
            ),
            "visual_palettes": (
                "temperature and precipitation use the simulator benchmark anchors; "
                "land elevation and sea depth use fixed metre palettes"
            ),
            "numeric_fields": "Single-channel IEEE float32 GeoTIFF, EPSG:4326, pixel-is-area, nodata=-9999.9; no display clamping",
            "vectors": "Eastward and northward components in physical units; northward is negated only for screen coordinates",
            "flow_previews": "Headless NumPy LIC/particles following benchmark palettes and integration settings; not pixel-identical to SFML; TIFFs store dimensionless texture intensity",
        },
        "outputs": records,
    }
    receipt_directory = args.output_root / "metadata" / "reference-maps"
    receipt_directory.mkdir(parents=True, exist_ok=True)
    manifest_path = receipt_directory / "manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        retained = [r for r in previous["outputs"] if r["width"] not in widths]
        for resolution in retained:
            for product in resolution["products"]:
                if product["name"] in ("precipitation", "temperature"):
                    product.setdefault("sources", [previous[product["name"] + "_source"]])
        receipt["outputs"] = sorted([*retained, *records], key=lambda r: r["width"])
    manifest_path.write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    for missing in records[0]["unavailable"]:
        print(f"Unavailable {missing['product']}: {missing['reason']}", flush=True)


if __name__ == "__main__":
    main()
