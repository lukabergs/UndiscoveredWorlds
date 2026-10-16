"""Resolution agreement of two tectonics runs (terrain roadmap T4.12): global statistics per checkpoint.

Compares two bundles (typically the same settings at two grid resolutions), each with its `<bundle>-checkpoints`
directory (`--checkpoint-every 50 --checkpoint-at 1`). Checkpoints are aligned by time; the final bundles are
compared as one more row. For every metric and checkpoint it prints both values, their difference, the tolerance
and pass or fail.

Agreement means that the global statistics below stay within their tolerances at every checkpoint from --from-myr
on: continental area and thickness, land share and hypsometry, plate motion and count, the share of overthickened
crust and the rates of the continental area terms. It never compares maps: two resolutions place continents,
boundaries and mountains differently, and that is not a disagreement. The t = 1 row (the seeding) is shown apart
and never judged.

Metrics: col_area, col_thick_km, surface_continent_km, land (engine cells, unweighted), over50_share,
mean_speed_cm_yr, plates and the area rates from continental_growth_metrics.series; land_median_m and the ten bands
(percent of the surface, spherical weights) from bundle_land_bands.measure.

Usage: uv run python scripts/benchmarks/resolution_check.py A B [--from-myr 50] [--json OUT]
Exit status 1 when a judged checkpoint fails.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from bundle_land_bands import BANDS as LAND_BANDS  # noqa: E402
from bundle_land_bands import measure  # noqa: E402
from continental_growth_metrics import series  # noqa: E402

# Tolerances: metric -> ("abs", largest difference) or ("rel", largest difference over the larger magnitude).
# Shares are fractions of the cells (0..1); bands are percentage points of the surface.
SHARE_TOLERANCE = 0.02  # col_area, land: 2 % of the cells.
THICKNESS_TOLERANCE_KM = 1.5  # col_thick_km: continental volume per continental area.
SURFACE_CONTINENT_TOLERANCE_KM = 0.5  # surface_continent_km: continental volume per surface area.
LAND_MEDIAN_TOLERANCE_M = 100.0  # Spherical land median.
BAND_TOLERANCE_POINTS = 1.5  # Each of bundle_land_bands' ten bands.
SPEED_TOLERANCE = 0.30  # mean_speed_cm_yr, relative.
PLATES_TOLERANCE = 2  # Plate count.
OVER50_TOLERANCE = 0.03  # Share of continental cells thicker than 50 km.
AREA_RATE_TOLERANCE = 0.30  # The four continental area rates, relative (when both snapshots have them).

TOLERANCES: dict[str, tuple[str, float]] = {
    "col_area": ("abs", SHARE_TOLERANCE),
    "col_thick_km": ("abs", THICKNESS_TOLERANCE_KM),
    "surface_continent_km": ("abs", SURFACE_CONTINENT_TOLERANCE_KM),
    "land": ("abs", SHARE_TOLERANCE),
    "land_median_m": ("abs", LAND_MEDIAN_TOLERANCE_M),
    **{f"{low}-{'up' if high is None else high}": ("abs", BAND_TOLERANCE_POINTS) for low, high in LAND_BANDS},
    "mean_speed_cm_yr": ("rel", SPEED_TOLERANCE),
    "plates": ("abs", PLATES_TOLERANCE),
    "over50_share": ("abs", OVER50_TOLERANCE),
    "collided_area_rate": ("rel", AREA_RATE_TOLERANCE),
    "stretched_area_rate": ("rel", AREA_RATE_TOLERANCE),
    "arc_area_rate": ("rel", AREA_RATE_TOLERANCE),
    "capped_area_rate": ("rel", AREA_RATE_TOLERANCE),
}
OPTIONAL = {"collided_area_rate", "stretched_area_rate", "arc_area_rate", "capped_area_rate"}
SEEDING_MYR = 1.0


def directories(bundle: Path) -> dict:
    """Checkpoint directories by time (as continental_growth_metrics.series finds them), and the final bundle."""
    found: dict = {}
    checkpoints = bundle.parent / f"{bundle.name}-checkpoints"
    if checkpoints.is_dir():
        for directory in sorted(checkpoints.iterdir()):
            match = re.search(r"-t(\d+(?:\.\d+)?)$", directory.name)
            if match and (directory / "fields").is_dir():
                found[float(match.group(1))] = directory
    found["final"] = bundle
    return found


def rows_by_time(bundle: Path) -> dict:
    paths = directories(bundle)
    rows = {row["time_myr"]: row for row in series(bundle, 1)}
    for time, path in paths.items():
        rows.setdefault(time, {}).update(measure(path))
    return rows


def final_time(bundle: Path):
    """The engine time of the final bundle's snapshot, Myr (None without one)."""
    snapshot = bundle / "tables" / "geosphere" / "tectonic_snapshot.json"
    return json.loads(snapshot.read_text()).get("time_myr") if snapshot.exists() else None


def number(value):
    if value is None:
        return None
    return float(value)


def judge(metric: str, a, b) -> dict:
    kind, tolerance = TOLERANCES[metric]
    a, b = number(a), number(b)
    entry = {"metric": metric, "a": a, "b": b, "kind": kind, "tolerance": tolerance}
    if a is None or b is None:
        entry.update(difference=None, status="n/a")
        return entry
    if math.isnan(a) or math.isnan(b):
        entry.update(difference=None, status="n/a" if math.isnan(a) and math.isnan(b) else "fail")
        return entry
    difference = b - a
    limit = tolerance * max(abs(a), abs(b)) if kind == "rel" else tolerance
    entry.update(difference=difference, status="pass" if abs(difference) <= limit + 1e-12 else "fail")
    return entry


def time_label(time) -> str:
    return "final" if time == "final" else f"{time:g}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("a", type=Path)
    parser.add_argument("b", type=Path)
    parser.add_argument("--from-myr", type=float, default=50.0)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows_a, rows_b = rows_by_time(args.a), rows_by_time(args.b)
    numeric = sorted(t for t in rows_a if t != "final" and t in rows_b)
    final_a, final_b = final_time(args.a), final_time(args.b)
    times = numeric + (["final"] if final_a == final_b else [])
    if final_a != final_b:
        print(f"final bundles at different times ({final_a} and {final_b} Myr): not compared")
    only = sorted(set(t for t in list(rows_a) + list(rows_b) if t != "final") - set(numeric))
    print(f"A: {args.a}\nB: {args.b}")
    if only:
        print(f"checkpoints in one run only (not compared): {', '.join(time_label(t) for t in only)}")
    results, failures = [], []
    header = f"{'t_myr':>7} {'metric':<22}{'A':>14}{'B':>14}{'B-A':>14}{'tolerance':>14}  status"
    section = None
    for time in times:
        seeding = time == SEEDING_MYR
        judged = not seeding and (time == "final" or time >= args.from_myr)
        title = "seeding (t = 1, not judged)" if seeding else f"checkpoints (judged from t = {args.from_myr:g})"
        if title != section:
            section = title
            print(f"\n{title}\n{header}")
        for metric in TOLERANCES:
            if metric in OPTIONAL and (rows_a[time].get(metric) is None or rows_b[time].get(metric) is None):
                continue
            entry = judge(metric, rows_a[time].get(metric), rows_b[time].get(metric))
            entry.update(time_myr=time, judged=judged)
            results.append(entry)
            if judged and entry["status"] == "fail":
                failures.append(entry)
            status = entry["status"] if judged or entry["status"] == "n/a" else f"({entry['status']})"
            tolerance = f"{'+-' if entry['kind'] == 'abs' else '+-rel '}{entry['tolerance']:g}"
            cells = [f"{v:>14.5g}" if v is not None else f"{'-':>14}" for v in (entry["a"], entry["b"], entry["difference"])]
            print(f"{time_label(time):>7} {metric:<22}{''.join(cells)}{tolerance:>14}  {status}")
    judged_count = sum(1 for r in results if r["judged"] and r["status"] != "n/a")
    print(f"\n{len(failures)} of {judged_count} judged comparisons fail (from t = {args.from_myr:g}, final included)")
    for entry in failures:
        print(f"  t = {time_label(entry['time_myr'])}: {entry['metric']} A {entry['a']:.5g} B {entry['b']:.5g}")
    if args.json:
        args.json.write_text(json.dumps({"a": str(args.a), "b": str(args.b), "from_myr": args.from_myr,
                                         "tolerances": {k: list(v) for k, v in TOLERANCES.items()},
                                         "results": results, "failures": len(failures)}, indent=1))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
