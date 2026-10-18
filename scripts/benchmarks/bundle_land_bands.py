"""Land share, shelf-to-land bands and land median of run bundles (spherical weights), against Earth's bands.

The bundle's surface is `geosphere.base_elevation` minus `geosphere.sea_level` (cell-centred W x W/2 grid, the
lat-lon export on a cube bundle,
north row first, cells weighted by the cosine of latitude). Earth: ETOPO 2048 bands of the terrain roadmap
(sea -3000 m to 0: 6.0, 3.1, 3.2, 5.0; land: 8.0, 8.2, 5.4, 4.4, 1.7, 1.5 % of the surface; land median 420 m).

Usage: uv run python scripts/benchmarks/bundle_land_bands.py DIR [DIR ...]   (a bundle or checkpoint directory)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import Bundle  # noqa: E402

BANDS = [(-3000, -2000), (-2000, -1000), (-1000, -200), (-200, 0),
         (0, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, None)]
EARTH = [6.0, 3.1, 3.2, 5.0, 8.0, 8.2, 5.4, 4.4, 1.7, 1.5]


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, cumulative[-1] / 2)])


def measure(directory: Path) -> dict:
    bundle = Bundle.open(directory)  # v1 or v2; a cube bundle's lat-lon export
    elevation = bundle.load("geosphere.base_elevation").astype(np.float64) - float(bundle.load("geosphere.sea_level"))
    rows = elevation.shape[0]
    latitude = (0.5 - (np.arange(rows) + 0.5) / rows) * np.pi
    weights = np.repeat(np.cos(latitude)[:, None], elevation.shape[1], axis=1)
    weights /= weights.sum()
    land = elevation > 0
    row = {"land_pct": 100 * float(weights[land].sum())}
    for low, high in BANDS:
        band = (elevation > low if low == 0 else elevation >= low) & ((elevation < high) if high is not None else True)
        band &= land if low >= 0 else ~land
        row[f"{low}-{'up' if high is None else high}"] = 100 * float(weights[band].sum())
    row["land_median_m"] = weighted_median(elevation[land], weights[land]) if land.any() else float("nan")
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    names = ["land_pct"] + [f"{low}-{'up' if high is None else high}" for low, high in BANDS] + ["land_median_m"]
    print(f"{'run':<60}" + "".join(f"{n:>11}" for n in names))
    print(f"{'Earth (ETOPO 2048)':<60}{29.2:>11.1f}" + "".join(f"{v:>11.1f}" for v in EARTH) + f"{420:>11.0f}")
    for directory in args.directories:
        row = measure(directory)
        print(f"{str(directory)[-60:]:<60}" + "".join(f"{row[n]:>11.1f}" for n in names))


if __name__ == "__main__":
    main()
