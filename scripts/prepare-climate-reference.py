#!/usr/bin/env python3
"""Convert WorldClim monthly GeoTIFFs into simulator-aligned validation grids."""

from __future__ import annotations

import argparse
import csv
import json
import math
import struct
from pathlib import Path

import numpy as np
from PIL import Image


MAGIC = b"UWCLIM1\0"
NODATA_LIMIT = -1.0e30


def resize_masked(source: Path, width: int, height: int) -> np.ndarray:
    with Image.open(source) as image:
        values = np.asarray(image, dtype=np.float32)
        nodata_text = image.tag_v2.get(42113)

    valid = np.isfinite(values) & (values > NODATA_LIMIT)

    if nodata_text is not None:
        nodata = float(nodata_text)
        valid &= ~np.isclose(values, nodata, rtol=0.0, atol=1.0e-5)
    numerator = np.where(valid, values, 0.0).astype(np.float32)
    weight = valid.astype(np.float32)
    target_size = (width, height)
    numerator_resized = np.asarray(
        Image.fromarray(numerator, mode="F").resize(target_size, Image.Resampling.BILINEAR),
        dtype=np.float32,
    )
    weight_resized = np.asarray(
        Image.fromarray(weight, mode="F").resize(target_size, Image.Resampling.BILINEAR),
        dtype=np.float32,
    )

    result = np.full((height, width), np.nan, dtype=np.float32)
    usable = weight_resized > 0.05
    result[usable] = numerator_resized[usable] / weight_resized[usable]
    return result


def write_bundle(
    grid_directory: Path,
    output_directory: Path,
    variable: str,
    width: int,
    height: int,
) -> tuple[Path, np.ndarray, list[dict[str, float]]]:
    output_path = output_directory / f"worldclim_{variable}_monthly.uwclim"
    climatology_sum = np.zeros((height, width), dtype=np.float64)
    valid_count = np.zeros((height, width), dtype=np.uint8)
    month_stats: list[dict[str, float]] = []

    with output_path.open("wb") as output:
        output.write(MAGIC)
        output.write(struct.pack("<IIII", 1, width, height, 12))
        output.write(variable.encode("ascii").ljust(16, b"\0"))

        for month in range(1, 13):
            source = grid_directory / f"wc2.1_10m_{variable}_{month:02d}.tif"
            if not source.exists():
                raise FileNotFoundError(source)

            values = resize_masked(source, width, height)
            finite = np.isfinite(values)
            climatology_sum[finite] += values[finite]
            valid_count[finite] += 1
            output.write(values.astype("<f4", copy=False).tobytes(order="C"))

            month_stats.append(
                {
                    "month": month,
                    "valid_cells": int(finite.sum()),
                    "mean": float(np.nanmean(values)),
                    "minimum": float(np.nanmin(values)),
                    "maximum": float(np.nanmax(values)),
                }
            )

    climatology = np.full((height, width), np.nan, dtype=np.float32)
    valid = valid_count > 0
    climatology[valid] = (climatology_sum[valid] / valid_count[valid]).astype(np.float32)
    return output_path, climatology, month_stats


def write_precipitation_csv(path: Path, precipitation: np.ndarray) -> None:
    height, width = precipitation.shape
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["y", "latitude", *(f"x{x}" for x in range(width))])

        for y in range(height):
            latitude = 90.0 - 180.0 * y / max(1, height - 1)
            row = ["nan" if not math.isfinite(float(value)) else f"{float(value):.4f}" for value in precipitation[y]]
            writer.writerow([y, f"{latitude:.6f}", *row])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("extra/reference/worldclim-2.1-10m/grids"),
        help="Directory containing extracted wc2.1_10m_*.tif files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("extra/reference"),
        help="Directory for simulator-aligned validation artifacts.",
    )
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--height", type=int, default=1025)
    parser.add_argument(
        "--variables",
        nargs="+",
        choices=("tavg", "prec", "wind", "srad", "vapr"),
        default=("tavg", "prec"),
        help="WorldClim variables to align; defaults to temperature and precipitation.",
    )
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    metadata_path = args.output / "worldclim-preparation.json"
    all_stats: dict[str, object] = {
        "dataset": "worldclim-2.1-1970-2000-10m",
        "target_width": args.width,
        "target_height": args.height,
        "alignment": "global_equirectangular_north_to_south",
        "variables": {},
    }

    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if (
            existing.get("dataset") == all_stats["dataset"]
            and existing.get("target_width") == args.width
            and existing.get("target_height") == args.height
        ):
            all_stats["variables"].update(existing.get("variables", {}))

    for variable in args.variables:
        bundle_path, climatology, stats = write_bundle(
            args.input, args.output, variable, args.width, args.height
        )
        all_stats["variables"][variable] = {
            "bundle": bundle_path.name,
            "months": stats,
        }

        if variable == "prec":
            write_precipitation_csv(
                args.output / "earth_precipitation_grid.csv", climatology
            )

    with metadata_path.open("w", encoding="utf-8") as output:
        json.dump(all_stats, output, indent=2, allow_nan=False)
        output.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
