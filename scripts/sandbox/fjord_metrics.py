"""Fjord measures of a sandbox window: how much of the sea is narrow inlets, how far they reach inland.

Sea = cells at or below 0 m connected to the window border (terrain_metrics.sea_from). A morphological opening
with an ellipse of --radius-km (5 km by default; pixels are square in degrees, so the ellipse is wider east-west
by 1 / cos(latitude)) removes everything narrower than twice the radius: what the opening removes is "narrow
sea" (fjords, rias, channels). Reported per run and domain:

- sea_pct: sea share of the window; narrow_pct: narrow sea as a share of the window; narrow_of_sea: its share of the sea
- inlets: connected components of narrow sea larger than 4 km^2; reach_km: the greatest distance (km, chessboard,
  north-south pixel size) of narrow sea from the wide sea, p90 over inlets; longest_km: the single longest
- coast_km_per_1000km2: coastline length per sea area (fjord coasts are high)
- glacial_m_iced: mean glacial erosion (m) on cells iced at any time, when `<domain>-glacial_eroded_m.npy` exists

Usage: uv run python scripts/sandbox/fjord_metrics.py LABEL1,LABEL2 [--domain G-coast] [--radius-km 5] [--json OUT]
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

EARTH_RADIUS_KM = 6371.0


def domain_geometry(label: str, domain: str) -> tuple[float, float]:
    """North-south pixel size (km) and the window's centre latitude (degrees) from the run's config."""
    config = json.loads((ROOT / "out/scratch/sandbox" / label / "config.json").read_text(encoding="utf-8"))
    spec = config["domains"][domain]
    bundle_columns = 256  # the base world's grid (configs/sandbox/terrain.json base_bundle)
    if "columns" in spec:
        pixel_km = 2 * math.pi * EARTH_RADIUS_KM / spec["columns"]
        return pixel_km, 0.0
    c0, r0, c1, r1 = spec["pixels"]
    k = spec["scale"]
    pixel_km = 2 * math.pi * EARTH_RADIUS_KM / bundle_columns / k
    rows = bundle_columns // 2
    latitude = 90 - (r0 + r1 + 1) / 2 / rows * 180
    return pixel_km, latitude


def ellipse(ry: int, rx: int) -> np.ndarray:
    y, x = np.ogrid[-ry:ry + 1, -rx:rx + 1]
    return (y / max(ry, 1)) ** 2 + (x / max(rx, 1)) ** 2 <= 1


def measure(label: str, domain: str, radius_km: float) -> dict:
    directory = ROOT / "out/scratch/sandbox" / label
    elevation = np.load(directory / f"{domain}-elevation.npy")
    pixel_km, latitude = domain_geometry(label, domain)
    ew_km = pixel_km * math.cos(math.radians(latitude))
    sea = tm.sea_from(elevation)
    ry = max(1, round(radius_km / pixel_km))
    rx = max(1, round(radius_km / ew_km))
    wide = ndimage.binary_opening(sea, structure=ellipse(ry, rx))
    narrow = sea & ~wide
    labels, count = ndimage.label(narrow, structure=np.ones((3, 3)))
    pixel_area = pixel_km * ew_km
    sizes = ndimage.sum(np.ones_like(labels), labels, index=np.arange(1, count + 1)) * pixel_area
    keep = np.flatnonzero(sizes >= 4) + 1
    # Distance of narrow sea from the wide sea, in north-south pixels (chessboard).
    distance = ndimage.distance_transform_cdt(~wide, metric="chessboard") * pixel_km
    reaches = [float(distance[labels == i].max()) for i in keep] if len(keep) else []
    coast = (sea != ndimage.binary_erosion(sea)).sum() * pixel_km
    row = {
        "sea_pct": 100 * float(sea.mean()),
        "narrow_pct": 100 * float(narrow.mean()),
        "narrow_of_sea": float(narrow.sum() / max(sea.sum(), 1)),
        "inlets": int(len(keep)),
        "reach_km_p90": float(np.percentile(reaches, 90)) if reaches else 0.0,
        "longest_km": max(reaches) if reaches else 0.0,
        "coast_km_per_1000km2": float(coast / max(sea.sum() * pixel_area, 1) * 1000),
    }
    glacial = directory / f"{domain}-glacial_eroded_m.npy"
    ice = directory / f"{domain}-ice_time_myr.npy"
    if glacial.exists() and ice.exists():
        iced = np.load(ice) > 0
        row["iced_pct_land"] = 100 * float((iced & ~sea).sum() / max((~sea).sum(), 1))
        row["glacial_m_iced"] = float(np.load(glacial)[iced].mean()) if iced.any() else 0.0
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("labels")
    parser.add_argument("--domain", default="G-coast")
    parser.add_argument("--radius-km", type=float, default=5.0)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows = {}
    for label in args.labels.split(","):
        rows[label] = measure(label, args.domain, args.radius_km)
    names = list(next(iter(rows.values())).keys())
    print(f"{'run':<18}" + "".join(f"{n:>20}" for n in names))
    for label, row in rows.items():
        print(f"{label:<18}" + "".join(f"{row.get(n, float('nan')):>20.3f}" for n in names))
    if args.json:
        args.json.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
