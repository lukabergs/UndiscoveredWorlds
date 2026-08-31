#!/usr/bin/env python3
"""Align ERA5 diagnostics and IPCC AR6 regions to the Earth benchmark grid."""

from __future__ import annotations

import argparse
import calendar
import json
import struct
from pathlib import Path

import numpy as np
import xarray as xr
from PIL import Image, ImageDraw


MAGIC = b"UWCLIM1\0"
ERA5_VARIABLES = {
    "slp": "slp_anom",
    "tcwv": "tcwv",
    "u10m": "u10m",
    "v10m": "v10m",
    "u500": "u500",
    "v500": "v500",
    "u850": "u850",
    "v850": "v850",
    "w500": "w500_ascent",
    "pr": "pr",
}


def resize(values: np.ndarray, width: int, height: int) -> np.ndarray:
    return np.asarray(
        Image.fromarray(values.astype(np.float32), mode="F").resize(
            (width, height), Image.Resampling.BILINEAR
        ),
        dtype=np.float32,
    )


def transform_era5(variable: str, values: np.ndarray, month: int) -> np.ndarray:
    if variable == "slp":
        return values / 100.0
    if variable == "w500":
        # ERA5 omega is positive downward in Pa/s. The simulator stores ascent-positive hPa/day.
        return values * -864.0
    if variable == "pr":
        average_days = sum(calendar.monthrange(year, month)[1] for year in range(2001, 2021)) / 20.0
        return values * (average_days * 86400.0)
    return values


def remove_area_weighted_mean(values: np.ndarray) -> np.ndarray:
    latitudes = np.linspace(90.0, -90.0, values.shape[0], dtype=np.float64)
    weights = np.maximum(0.0, np.cos(np.deg2rad(latitudes)))[:, None]
    finite = np.isfinite(values)
    total_weight = np.sum(weights * finite)
    if total_weight > 0.0:
        values = values - np.sum(np.where(finite, values, 0.0) * weights) / total_weight
    return values


def write_monthly_grid(path: Path, variable: str, months: list[np.ndarray]) -> None:
    height, width = months[0].shape
    with path.open("wb") as output:
        output.write(MAGIC)
        output.write(struct.pack("<IIII", 1, width, height, len(months)))
        output.write(variable.encode("ascii").ljust(16, b"\0"))
        for values in months:
            output.write(values.astype("<f4", copy=False).tobytes(order="C"))


def prepare_era5(source: Path, output: Path, width: int, height: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for source_variable, output_variable in ERA5_VARIABLES.items():
        matches = sorted(source.glob(f"era5_{source_variable}_monthly_climatology_*.nc"))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one ERA5 {source_variable} climatology in {source}, found {matches}")

        dataset = xr.open_dataset(matches[0])
        data = dataset[source_variable]
        months: list[np.ndarray] = []
        statistics: list[dict[str, float | int]] = []

        for month in range(1, 13):
            values = np.asarray(data.sel(month=month), dtype=np.float32)
            values = transform_era5(source_variable, values, month)
            values = resize(values, width, height)
            if source_variable == "slp":
                values = remove_area_weighted_mean(values)
            months.append(values)
            statistics.append(
                {
                    "month": month,
                    "mean": float(np.nanmean(values)),
                    "minimum": float(np.nanmin(values)),
                    "maximum": float(np.nanmax(values)),
                }
            )

        dataset.close()
        output_path = output / f"era5_{output_variable}_monthly.uwclim"
        write_monthly_grid(output_path, output_variable, months)
        records.append(
            {
                "source_variable": source_variable,
                "variable": output_variable,
                "source": matches[0].as_posix(),
                "bundle": output_path.name,
                "months": statistics,
            }
        )
    return records


def polygon_mask(
    rings: list[list[list[float]]], width: int, height: int
) -> Image.Image:
    mask = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(mask)

    def pixels(ring: list[list[float]]) -> list[tuple[float, float]]:
        return [
            (
                (longitude + 180.0) * (width - 1) / 360.0,
                (90.0 - latitude) * (height - 1) / 180.0,
            )
            for longitude, latitude in ring
        ]

    if rings:
        draw.polygon(pixels(rings[0]), fill=255)
    for hole in rings[1:]:
        draw.polygon(pixels(hole), fill=0)
    return mask


def prepare_regions(source: Path, output: Path, width: int, height: int) -> dict[str, object]:
    collection = json.loads(source.read_text(encoding="utf-8"))
    region_image = Image.new("L", (width, height), 0)
    rows = ["id\tacronym\tname\tcontinent\ttype"]

    for feature in collection["features"]:
        properties = feature["properties"]
        region_id = int(properties["id"])
        geometry = feature["geometry"]
        polygons = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]

        for rings in polygons:
            mask = polygon_mask(rings, width, height)
            region_image.paste(region_id + 1, mask=mask)

        rows.append(
            "\t".join(
                (
                    str(region_id),
                    properties["Acronym"],
                    properties["Name"],
                    properties["Continent"],
                    properties["Type"],
                )
            )
        )

    mask_path = output / "ipcc_ar6_regions.png"
    labels_path = output / "ipcc_ar6_regions.tsv"
    region_image.save(mask_path)
    labels_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return {
        "source": source.as_posix(),
        "mask": mask_path.name,
        "labels": labels_path.name,
        "region_count": len(collection["features"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--era5",
        type=Path,
        default=Path("extra/reference/era5-planette-2001-2020"),
    )
    parser.add_argument(
        "--regions",
        type=Path,
        default=Path("extra/reference/ipcc-ar6-wgi-regions-v4/IPCC-WGI-reference-regions-v4.geojson"),
    )
    parser.add_argument("--output", type=Path, default=Path("extra/reference"))
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--height", type=int, default=1025)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    receipt = {
        "target_width": args.width,
        "target_height": args.height,
        "alignment": "global_equirectangular_north_to_south",
        "era5": prepare_era5(args.era5, args.output, args.width, args.height),
        "ipcc_ar6_regions": prepare_regions(
            args.regions, args.output, args.width, args.height
        ),
    }
    (args.output / "diagnostic-reference-preparation.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
