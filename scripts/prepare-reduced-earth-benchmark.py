#!/usr/bin/env python3
"""Create a smaller, mutually exclusive Earth land/sea map pair."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from PIL import Image


def load_height(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        values = np.asarray(image, dtype=np.uint16)

    if values.ndim != 2:
        raise ValueError(f"Expected a single-channel height map: {path}")
    return values


def box_average(values: np.ndarray, target_size: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    present = values > 0
    weighted = np.asarray(
        Image.fromarray(values.astype(np.float32)).resize(target_size, Image.Resampling.BOX),
        dtype=np.float32,
    )
    coverage = np.asarray(
        Image.fromarray(present.astype(np.float32)).resize(target_size, Image.Resampling.BOX),
        dtype=np.float32,
    )
    return weighted, coverage


def reduce_pair(land: np.ndarray, sea: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray]:
    if land.shape != sea.shape:
        raise ValueError(f"Land and sea dimensions differ: {land.shape} != {sea.shape}")
    if np.any((land > 0) & (sea > 0)):
        raise ValueError("Source land and sea maps overlap")

    height, width = land.shape
    target_size = (math.ceil(width / factor), math.ceil(height / factor))
    land_weighted, land_coverage = box_average(land, target_size)
    sea_weighted, sea_coverage = box_average(sea, target_size)

    land_selected = (land_coverage >= sea_coverage) & (land_coverage > 0.0)
    sea_selected = ~land_selected & (sea_coverage > 0.0)
    reduced_land = np.zeros((target_size[1], target_size[0]), dtype=np.uint16)
    reduced_sea = np.zeros_like(reduced_land)

    reduced_land[land_selected] = np.clip(
        np.rint(land_weighted[land_selected] / land_coverage[land_selected]), 1, 65535
    ).astype(np.uint16)
    reduced_sea[sea_selected] = np.clip(
        np.rint(sea_weighted[sea_selected] / sea_coverage[sea_selected]), 1, 65535
    ).astype(np.uint16)
    return reduced_land, reduced_sea


def save_height(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(values).save(path)


def reduce_reference(path: Path, source_shape: tuple[int, int], target_shape: tuple[int, int]) -> np.ndarray:
    source_height, source_width = source_shape
    values = np.genfromtxt(
        path,
        delimiter=",",
        skip_header=1,
        usecols=range(2, source_width + 2),
        dtype=np.float32,
    )
    if values.shape != source_shape:
        raise ValueError(f"Reference dimensions differ: {values.shape} != {source_shape}")

    valid = np.isfinite(values)
    target_size = (target_shape[1], target_shape[0])
    numerator = np.asarray(
        Image.fromarray(np.where(valid, values, 0.0).astype(np.float32)).resize(
            target_size, Image.Resampling.BOX
        ),
        dtype=np.float32,
    )
    weight = np.asarray(
        Image.fromarray(valid.astype(np.float32)).resize(target_size, Image.Resampling.BOX),
        dtype=np.float32,
    )
    result = np.full(target_shape, np.nan, dtype=np.float32)
    usable = weight > 0.01
    result[usable] = numerator[usable] / weight[usable]
    return result


def save_reference(path: Path, values: np.ndarray) -> None:
    height, width = values.shape
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["y", "latitude", *(f"x{x}" for x in range(width))])

        for y, row in enumerate(values):
            latitude = 90.0 - 180.0 * y / max(1, height - 1)
            writer.writerow(
                [
                    y,
                    f"{latitude:.6f}",
                    *("nan" if not np.isfinite(value) else f"{value:.6f}" for value in row),
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor", type=int, default=4)
    parser.add_argument("--land", type=Path, default=Path("extra/img/earth/in/earth_land_l_3.png"))
    parser.add_argument("--sea", type=Path, default=Path("extra/img/earth/in/earth_sea_l_1.png"))
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("extra/reference/earth_precipitation_grid.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("out/reduced-earth-benchmark"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.factor < 2:
        raise ValueError("--factor must be at least 2")

    source_land = load_height(args.land)
    source_sea = load_height(args.sea)
    reduced_land, reduced_sea = reduce_pair(source_land, source_sea, args.factor)
    width = reduced_land.shape[1]
    height = reduced_land.shape[0]
    land_path = args.output_dir / f"earth_land_{width}x{height}.png"
    sea_path = args.output_dir / f"earth_sea_{width}x{height}.png"
    save_height(land_path, reduced_land)
    save_height(sea_path, reduced_sea)
    reference_path = args.output_dir / f"earth_precipitation_{width}x{height}.csv"
    reduced_reference = reduce_reference(
        args.reference, source_land.shape, reduced_land.shape
    )
    save_reference(reference_path, reduced_reference)

    overlap = int(np.count_nonzero((reduced_land > 0) & (reduced_sea > 0)))
    uncovered = int(np.count_nonzero((reduced_land == 0) & (reduced_sea == 0)))
    print(f"Created {width}x{height} Earth benchmark maps at {args.output_dir}")
    print(
        f"land_cells={np.count_nonzero(reduced_land)} "
        f"sea_cells={np.count_nonzero(reduced_sea)} overlap={overlap} uncovered={uncovered}"
    )
    print(f"reference_valid_cells={np.count_nonzero(np.isfinite(reduced_reference))}")


if __name__ == "__main__":
    main()
