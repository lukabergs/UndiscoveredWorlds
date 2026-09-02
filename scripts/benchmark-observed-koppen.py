#!/usr/bin/env python3
"""Classify observed monthly climate grids and compare them with the Earth atlas."""

from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

import numpy as np
from PIL import Image


MAGIC = b"UWCLIM1\0"
CLIMATE_CODES = (
    "Af", "Am", "Aw", "As", "BWh", "BWk", "BSh", "BSk",
    "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc",
    "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd",
    "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF",
)
CLIMATE_COLOURS = np.asarray(
    (
        (0, 0, 254), (1, 119, 255), (70, 169, 250), (70, 169, 250),
        (249, 15, 0), (251, 150, 149), (245, 163, 1), (254, 219, 99),
        (255, 255, 0), (198, 199, 1), (184, 184, 114), (138, 255, 162),
        (86, 199, 112), (30, 150, 66), (192, 254, 109), (76, 255, 93),
        (19, 203, 74), (255, 8, 245), (204, 3, 192), (154, 51, 144),
        (153, 100, 146), (172, 178, 249), (91, 121, 213), (78, 83, 175),
        (54, 3, 130), (0, 255, 245), (32, 200, 250), (0, 126, 125),
        (0, 69, 92), (178, 178, 178), (104, 104, 104),
    ),
    dtype=np.uint8,
)
COLOUR_DISTANCE_LIMIT_SQUARED = 65 * 65


def load_bundle(path: Path, expected_variable: str) -> np.ndarray:
    with path.open("rb") as source:
        if source.read(8) != MAGIC:
            raise ValueError(f"Invalid climate bundle magic: {path}")
        version, width, height, months = struct.unpack("<IIII", source.read(16))
        variable = source.read(16).split(b"\0", 1)[0].decode("ascii")
        values = np.fromfile(source, dtype="<f4")

    if version != 1 or variable != expected_variable or months != 12:
        raise ValueError(
            f"Unexpected bundle header in {path}: version={version}, "
            f"variable={variable!r}, months={months}"
        )
    expected_values = width * height * months
    if values.size != expected_values:
        raise ValueError(
            f"Truncated climate bundle {path}: expected {expected_values} values, "
            f"found {values.size}"
        )
    return values.reshape((months, height, width))


def set_where(target: np.ndarray, mask: np.ndarray, value: int) -> None:
    target[mask & (target == 0)] = value


def classify_monthly(temperature: np.ndarray, precipitation: np.ndarray) -> np.ndarray:
    """Apply the standard 12-month Köppen-Geiger decision rules used by the control."""
    _, height, width = temperature.shape
    valid = np.all(np.isfinite(temperature), axis=0)
    valid &= np.all(np.isfinite(precipitation), axis=0)
    valid &= np.all(precipitation >= 0.0, axis=0)

    mean_temperature = np.mean(temperature, axis=0)
    minimum_temperature = np.min(temperature, axis=0)
    maximum_temperature = np.max(temperature, axis=0)
    months_above_ten = np.count_nonzero(temperature >= 10.0, axis=0)
    annual_precipitation = np.sum(precipitation, axis=0)
    minimum_precipitation = np.min(precipitation, axis=0)

    northern = np.arange(height)[:, None] <= (height - 1) / 2.0
    northern = np.broadcast_to(northern, (height, width))
    north_summer = precipitation[3:9]
    north_winter = np.concatenate((precipitation[9:12], precipitation[0:3]), axis=0)
    south_summer = north_winter
    south_winter = north_summer
    driest_summer = np.where(northern, np.min(north_summer, axis=0), np.min(south_summer, axis=0))
    wettest_summer = np.where(northern, np.max(north_summer, axis=0), np.max(south_summer, axis=0))
    driest_winter = np.where(northern, np.min(north_winter, axis=0), np.min(south_winter, axis=0))
    wettest_winter = np.where(northern, np.max(north_winter, axis=0), np.max(south_winter, axis=0))
    summer_precipitation = np.where(
        northern, np.sum(north_summer, axis=0), np.sum(south_summer, axis=0)
    )
    summer_fraction = np.divide(
        summer_precipitation,
        annual_precipitation,
        out=np.zeros_like(annual_precipitation),
        where=annual_precipitation > 0.0,
    )

    dry_threshold = mean_temperature * 20.0
    dry_threshold += np.where(summer_fraction >= 0.7, 280.0, np.where(summer_fraction >= 0.3, 140.0, 0.0))

    group = np.zeros((height, width), dtype=np.uint8)
    group[valid & (maximum_temperature < 0.0)] = 5
    group[valid & (maximum_temperature >= 0.0) & (maximum_temperature <= 10.0)] = 4
    set_where(group, valid & (annual_precipitation < dry_threshold * 0.5), 2)
    set_where(group, valid & (annual_precipitation <= dry_threshold), 3)
    set_where(group, valid & (minimum_temperature >= 18.0), 1)
    set_where(group, valid & (minimum_temperature > -3.0), 6)
    set_where(group, valid, 7)

    climate = np.zeros((height, width), dtype=np.uint8)

    tropical = group == 1
    set_where(climate, tropical & (minimum_precipitation >= 60.0), 1)
    set_where(
        climate,
        tropical & (minimum_precipitation >= (100.0 - annual_precipitation / 25.0)),
        2,
    )
    set_where(climate, tropical & (driest_summer < driest_winter), 4)
    set_where(climate, tropical, 3)

    climate[(group == 2) & (mean_temperature >= 18.0)] = 5
    climate[(group == 2) & (mean_temperature < 18.0)] = 6
    climate[(group == 3) & (mean_temperature >= 18.0)] = 7
    climate[(group == 3) & (mean_temperature < 18.0)] = 8

    for temperate_group, base, moisture_offsets in (
        (6, 9, (0, 3, 6)),
        (7, 18, (0, 4, 8)),
    ):
        cells = group == temperate_group
        dry_summer = cells & (driest_summer < wettest_winter / 3.0) & (driest_summer < 40.0)
        dry_winter = cells & ~dry_summer & (driest_winter < wettest_summer / 10.0)
        fully_humid = cells & ~dry_summer & ~dry_winter
        for moisture, moisture_offset in zip(
            (dry_summer, dry_winter, fully_humid), moisture_offsets
        ):
            climate[moisture & (maximum_temperature >= 22.0)] = base + moisture_offset
            climate[moisture & (maximum_temperature < 22.0) & (months_above_ten >= 4)] = base + moisture_offset + 1
            climate[
                moisture
                & (maximum_temperature < 22.0)
                & (months_above_ten < 4)
                & (minimum_temperature <= -38.0)
            ] = base + moisture_offset + 3
            climate[
                moisture
                & (maximum_temperature < 22.0)
                & (months_above_ten < 4)
                & (minimum_temperature > -38.0)
            ] = base + moisture_offset + 2

    climate[group == 4] = 30
    climate[group == 5] = 31
    return climate


def decode_reference(reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pixels = reference.astype(np.int16)
    colours = CLIMATE_COLOURS.astype(np.int16)
    best_class = np.zeros(reference.shape[:2], dtype=np.uint8)
    best_distance = np.full(reference.shape[:2], np.iinfo(np.int32).max, dtype=np.int32)

    for index, colour in enumerate(colours, start=1):
        difference = pixels.astype(np.int32) - colour.astype(np.int32)
        distance = np.sum(difference * difference, axis=2, dtype=np.int32)
        closer = distance < best_distance
        best_distance[closer] = distance[closer]
        best_class[closer] = index
    return best_class, best_distance <= COLOUR_DISTANCE_LIMIT_SQUARED


def comparable(classes: np.ndarray) -> np.ndarray:
    result = classes.copy()
    result[result == 4] = 3
    return result


def major_group(classes: np.ndarray) -> np.ndarray:
    result = np.full(classes.shape, -1, dtype=np.int8)
    result[(classes >= 1) & (classes <= 4)] = 0
    result[(classes >= 5) & (classes <= 8)] = 1
    result[(classes >= 9) & (classes <= 17)] = 2
    result[(classes >= 18) & (classes <= 29)] = 3
    result[(classes >= 30) & (classes <= 31)] = 4
    return result


def compare(
    simulated: np.ndarray,
    reference: np.ndarray,
    reference_valid: np.ndarray,
    land: np.ndarray,
) -> dict[str, object]:
    simulated = comparable(simulated)
    reference = comparable(reference)
    compared = land & reference_valid & (simulated > 0)
    latitude = np.linspace(90.0, -90.0, simulated.shape[0])
    weights = np.cos(np.deg2rad(latitude))[:, None]
    cell_weights = np.broadcast_to(weights, simulated.shape)
    total_weight = float(np.sum(cell_weights[compared]))
    exact = simulated == reference
    group_exact = major_group(simulated) == major_group(reference)

    confusion = np.zeros((31, 31), dtype=np.float64)
    np.add.at(
        confusion,
        (simulated[compared].astype(int) - 1, reference[compared].astype(int) - 1),
        cell_weights[compared],
    )
    observed = float(np.trace(confusion) / total_weight) if total_weight else math.nan
    simulated_totals = np.sum(confusion, axis=1)
    reference_totals = np.sum(confusion, axis=0)
    expected = float(np.sum(simulated_totals * reference_totals) / (total_weight * total_weight)) if total_weight else math.nan
    kappa = (observed - expected) / (1.0 - expected) if expected < 1.0 else math.nan

    counts = {
        code: int(np.count_nonzero(land & (simulated == index)))
        for index, code in enumerate(CLIMATE_CODES, start=1)
    }
    return {
        "land_cells": int(np.count_nonzero(land)),
        "classified_land_cells": int(np.count_nonzero(land & (simulated > 0))),
        "compared_land_cells": int(np.count_nonzero(compared)),
        "classified_land_fraction": float(np.count_nonzero(land & (simulated > 0)) / np.count_nonzero(land)),
        "area_weighted_exact_accuracy": float(np.sum(cell_weights[compared & exact]) / total_weight),
        "area_weighted_major_group_accuracy": float(np.sum(cell_weights[compared & group_exact]) / total_weight),
        "area_weighted_spatial_kappa": kappa,
        "class_counts": counts,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--temperature",
        type=Path,
        default=Path("extra/reference/climate/processed/worldclim_tavg_monthly.uwclim"),
    )
    parser.add_argument(
        "--precipitation",
        type=Path,
        default=Path("extra/reference/climate/processed/imerg_prec_monthly.uwclim"),
    )
    parser.add_argument(
        "--land-mask",
        type=Path,
        default=Path("extra/reference/earth/base-maps/earth_land_l_3.png"),
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("extra/reference/climate/previews/0.png"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("extra/validation/observed-climate"),
    )
    parser.add_argument("--name", default="worldclim-temperature_imerg-precipitation")
    args = parser.parse_args()

    temperature = load_bundle(args.temperature, "tavg")
    precipitation = load_bundle(args.precipitation, "prec")
    if temperature.shape != precipitation.shape:
        raise ValueError(
            f"Temperature shape {temperature.shape} does not match precipitation shape {precipitation.shape}"
        )

    land = np.asarray(Image.open(args.land_mask)) > 0
    reference_image = np.asarray(Image.open(args.reference).convert("RGB"))
    if land.shape != temperature.shape[1:] or reference_image.shape[:2] != land.shape:
        raise ValueError("Land, climate, and reference grids are not aligned")

    climate = classify_monthly(temperature, precipitation)
    climate[~land] = 0
    reference, reference_valid = decode_reference(reference_image)
    metrics = compare(climate, reference, reference_valid, land)
    metrics.update(
        {
            "name": args.name,
            "temperature_bundle": str(args.temperature),
            "precipitation_bundle": str(args.precipitation),
            "land_mask": str(args.land_mask),
            "reference": str(args.reference),
            "rules": {
                "months": 12,
                "cold_temperate_boundary_c": -3.0,
                "hot_summer_boundary_c": 22.0,
                "warm_month_boundary_c": 10.0,
                "dry_summer_ratio": 3.0,
                "dry_summer_limit_mm": 40.0,
                "dry_winter_ratio": 10.0,
            },
        }
    )

    output_image = np.zeros((*land.shape, 3), dtype=np.uint8)
    classified = climate > 0
    output_image[classified] = CLIMATE_COLOURS[climate[classified] - 1]
    args.output.mkdir(parents=True, exist_ok=True)
    image_path = args.output / f"{args.name}.png"
    metrics_path = args.output / f"{args.name}.json"
    Image.fromarray(output_image, mode="RGB").save(image_path)
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
