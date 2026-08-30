#!/usr/bin/env python3
"""Convert NASA IMERG annual precipitation into a simulator-aligned validation grid."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from PIL import Image


MAGIC = b"UWCLIM1\0"


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "extra/reference/imerg-v07b-annual/grids/"
            "IMERG-Final.CLIM.200006-202305.V07B.tif"
        ),
    )
    parser.add_argument("--output", type=Path, default=Path("extra/reference"))
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--height", type=int, default=1025)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    values = resize_masked(args.input, args.width, args.height)
    output_path = args.output / "imerg_prec_annual.uwclim"
    write_bundle(output_path, values)

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
