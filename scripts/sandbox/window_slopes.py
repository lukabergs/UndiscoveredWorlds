"""Slope and relief percentiles of sandbox windows at their own cells, beside Earth box-averaged to the same cells.

Items: LABEL/DOMAIN (a sandbox window: `out/scratch/sandbox/LABEL/DOMAIN-elevation.npy`, the grid from its
config.json as fjord_metrics.domain_geometry reads it) or atlas:STEM@F (refs/processed/earth-atlas/STEM.npy with
STEM-sea.npy, box-averaged by the integer factor F to the window's cells: 6 takes the 102 x 64 m GLO-30 windows to
G-coast's 611 x 382 m). Land = cells above 0 m that are not sea (sandbox: terrain_metrics.sea_from; atlas: the
sea mask). Slopes by central differences in true km (east-west scaled by cos latitude); relief = max - min in
square boxes of about 1 km and 5 km (north-south pixels) over land; land height percentiles.

Usage: uv run python scripts/sandbox/window_slopes.py --items t15-parent/G-coast,atlas:glo30-sognefjord-98km@6 [--json OUT]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "sandbox"))
import terrain_metrics as tm  # noqa: E402
from fjord_metrics import domain_geometry  # noqa: E402


def box_average(values: np.ndarray, factor: int) -> np.ndarray:
    rows, columns = (values.shape[0] // factor) * factor, (values.shape[1] // factor) * factor
    return values[:rows, :columns].reshape(rows // factor, factor, columns // factor, factor).mean(axis=(1, 3))


def load(item: str) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Elevation (m), sea mask, north-south pixel km, east-west pixel km."""
    if item.startswith("atlas:"):
        stem, _, factor = item[6:].partition("@")
        factor = int(factor or 1)
        directory = ROOT / "refs/processed/earth-atlas"
        elevation = np.load(directory / f"{stem}.npy").astype(np.float64)
        sea = np.load(directory / f"{stem}-sea.npy").astype(bool)
        receipt = json.loads((directory / f"{stem}.json").read_text(encoding="utf-8"))
        ns = float(receipt.get("pixel_km_ns") or receipt.get("pixel_ns_km") or receipt["pixel_km"][0])
        ew = float(receipt.get("pixel_km_ew") or receipt.get("pixel_ew_km") or receipt["pixel_km"][1])
        if factor > 1:
            elevation = box_average(elevation, factor)
            sea = box_average(sea.astype(np.float64), factor) >= 0.5
            ns, ew = ns * factor, ew * factor
        return elevation, sea, ns, ew
    label, domain = item.split("/")
    elevation = np.load(ROOT / "out/scratch/sandbox" / label / f"{domain}-elevation.npy").astype(np.float64)
    sea = tm.sea_from(elevation)
    ns, latitude = domain_geometry(label, domain)
    return elevation, sea, ns, ns * math.cos(math.radians(latitude))


def measure(item: str) -> dict:
    elevation, sea, ns, ew = load(item)
    land = ~sea & (elevation > 0)
    dzdy, dzdx = np.gradient(elevation, ns * 1000, ew * 1000)
    slope = np.hypot(dzdx, dzdy)
    row = {"pixel_km": f"{ns:.3f}x{ew:.3f}", "sea_pct": 100 * float(sea.mean())}
    for p in (50, 90, 99):
        row[f"slope_p{p}"] = float(np.percentile(slope[land], p)) if land.any() else float("nan")
    for km in (1, 5):
        size = max(3, int(round(km / ns)) | 1)
        relief = ndimage.maximum_filter(elevation, size) - ndimage.minimum_filter(elevation, size)
        row[f"relief_{km}km_p50"] = float(np.percentile(relief[land], 50)) if land.any() else float("nan")
        row[f"relief_{km}km_p90"] = float(np.percentile(relief[land], 90)) if land.any() else float("nan")
    for p in (10, 50, 90):
        row[f"land_p{p}"] = float(np.percentile(elevation[land], p)) if land.any() else float("nan")
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--items", required=True)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows = {item: measure(item) for item in args.items.split(",")}
    names = list(next(iter(rows.values())).keys())
    print(f"{'item':<34}" + "".join(f"{n:>15}" for n in names))
    for item, row in rows.items():
        cells = [f"{v:>15.3f}" if isinstance(v, float) else f"{str(v):>15}" for v in row.values()]
        print(f"{item:<34}" + "".join(cells))
    if args.json:
        args.json.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
