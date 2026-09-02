#!/usr/bin/env python3
"""Convert NASA IMERG precipitation into simulator-aligned validation grids."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from PIL import Image, TiffImagePlugin


MAGIC = b"UWCLIM1\0"
CLIMATOLOGICAL_MONTH_DAYS = (31.0, 28.25, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0)


def resize_masked(source: Path, width: int, height: int) -> np.ndarray:
    with Image.open(source) as image:
        values = np.asarray(image, dtype=np.float32)

    valid = np.isfinite(values) & (values >= 0.0)
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
    usable = weight_resized > 0.5
    result[usable] = numerator_resized[usable] / weight_resized[usable]
    return result


def write_bundle(path: Path, values: np.ndarray) -> None:
    height, width = values.shape

    with path.open("wb") as output:
        output.write(MAGIC)
        output.write(struct.pack("<IIII", 1, width, height, 1))
        output.write(b"prec_annual".ljust(16, b"\0"))
        output.write(values.astype("<f4", copy=False).tobytes(order="C"))


def write_monthly_bundle(
    input_directory: Path,
    path: Path,
    width: int,
    height: int,
) -> tuple[list[dict[str, float | int]], np.ndarray]:
    month_stats: list[dict[str, float | int]] = []
    annual_total = np.zeros((height, width), dtype=np.float32)
    annual_valid = np.ones((height, width), dtype=bool)

    with path.open("wb") as output:
        output.write(MAGIC)
        output.write(struct.pack("<IIII", 1, width, height, 12))
        output.write(b"prec".ljust(16, b"\0"))

        for month, days in enumerate(CLIMATOLOGICAL_MONTH_DAYS, start=1):
            source = input_directory / f"IMERG-Final.CLIM.2001-2022.{month:02d}.V07B.tif"

            if not source.exists():
                raise FileNotFoundError(source)

            daily_rate = resize_masked(source, width, height)
            monthly_total = daily_rate * np.float32(days)
            finite = np.isfinite(monthly_total)
            annual_total[finite] += monthly_total[finite]
            annual_valid &= finite
            output.write(monthly_total.astype("<f4", copy=False).tobytes(order="C"))
            month_stats.append(
                {
                    "month": month,
                    "days": days,
                    "valid_cells": int(finite.sum()),
                    "mean_daily_rate_mm": float(np.nanmean(daily_rate)),
                    "minimum_daily_rate_mm": float(np.nanmin(daily_rate)),
                    "maximum_daily_rate_mm": float(np.nanmax(daily_rate)),
                    "mean_monthly_total_mm": float(np.nanmean(monthly_total)),
                    "maximum_monthly_total_mm": float(np.nanmax(monthly_total)),
                }
            )

    annual_total[~annual_valid] = np.nan
    return month_stats, annual_total


def write_float32_geotiff(path: Path, values: np.ndarray) -> None:
    height, width = values.shape
    stored = np.where(np.isfinite(values), values, -9999.9).astype(np.float32)
    geotags = TiffImagePlugin.ImageFileDirectory_v2()
    geotags[33550] = (360.0 / width, 180.0 / height, 0.0)
    geotags[33922] = (0.0, 0.0, 0.0, -180.0, 90.0, 0.0)
    geotags[34735] = (
        1, 1, 0, 3,
        1024, 0, 1, 2,
        1025, 0, 1, 1,
        2048, 0, 1, 4326,
    )
    geotags[42113] = "-9999.9"
    Image.fromarray(stored, mode="F").save(
        path, format="TIFF", compression="raw", tiffinfo=geotags
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "extra/reference/climate/sources/imerg-v07b-annual/grids/"
            "IMERG-Final.CLIM.200006-202305.V07B.tif"
        ),
    )
    parser.add_argument(
        "--monthly-input",
        type=Path,
        help=(
            "Directory containing the twelve 2001-2022 monthly climatology GeoTIFFs. "
            "When supplied, writes a twelve-month bundle instead of the annual bundle."
        ),
    )
    parser.add_argument("--output", type=Path, default=Path("extra/reference/climate/processed"))
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--height", type=int, default=1025)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    if args.monthly_input is not None:
        output_path = args.output / "imerg_prec_monthly.uwclim"
        month_stats, annual_total = write_monthly_bundle(
            args.monthly_input, output_path, args.width, args.height
        )
        annual_bundle_path = args.output / "imerg_prec_annual_2001_2022.uwclim"
        annual_tiff_path = args.output / "imerg_precipitation_mm_year.tif"
        write_bundle(annual_bundle_path, annual_total)
        write_float32_geotiff(annual_tiff_path, annual_total)
        metadata = {
            "dataset": "imerg-final-v07b-monthly-climatology-2001-2022",
            "source": "https://gpm.nasa.gov/data/imerg/precipitation-climatology",
            "period": "2001-2022",
            "source_units": "mm/day",
            "bundle_units": "mm/month",
            "calendar": "365.25-day climatological year",
            "missing_value": -9999.9,
            "target_width": args.width,
            "target_height": args.height,
            "alignment": "global_equirectangular_north_to_south",
            "bundle": output_path.name,
            "annual_bundle": annual_bundle_path.name,
            "annual_geotiff": annual_tiff_path.name,
            "annual_geotiff_sample": "single-band IEEE float32",
            "annual_geotiff_units": "mm/year",
            "annual_geotiff_nodata": -9999.9,
            "months": month_stats,
        }
        metadata_path = args.output / "imerg-monthly-preparation.json"
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        return 0

    values = resize_masked(args.input, args.width, args.height)
    output_path = args.output / "imerg_prec_annual.uwclim"
    annual_tiff_path = args.output / "imerg_precipitation_mm_year.tif"
    write_bundle(output_path, values)
    write_float32_geotiff(annual_tiff_path, values)

    finite = np.isfinite(values)
    metadata = {
        "dataset": "imerg-final-v07b-grand-average-2000-2023",
        "source": "https://gpm.nasa.gov/data/imerg/precipitation-climatology",
        "period": "2000-06/2023-05",
        "units": "mm/year",
        "target_width": args.width,
        "target_height": args.height,
        "alignment": "global_equirectangular_north_to_south",
        "bundle": output_path.name,
        "geotiff": annual_tiff_path.name,
        "valid_cells": int(finite.sum()),
        "mean": float(np.nanmean(values)),
        "minimum": float(np.nanmin(values)),
        "maximum": float(np.nanmax(values)),
    }
    metadata_path = args.output / "imerg-preparation.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
