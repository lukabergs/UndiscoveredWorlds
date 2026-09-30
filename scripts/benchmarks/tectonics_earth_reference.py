"""Distribution-only Earth context; never supplies forcing to the simulator."""
import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

from tectonics_experiment import digest
from tectonics_report_data import hypsography_shape


QUANTILES = np.linspace(0, 1, 201)
CACHE_SCHEMA = "tectonics-earth-reference/v1"


def reference_inputs(root, width):
    return [root / f"refs/processed/heightmap/{width}.tif",
            root / f"refs/processed/land_ocean/maps/earth_land_ocean_{width}x{width // 2}.png",
            root / "refs/README.md",
            root / f"refs/processed/metadata/{width}/verification.json",
            root / f"refs/processed/metadata/{width}/manifest.json"]


def input_hashes(root, width):
    return {p.relative_to(root).as_posix(): digest(p) if p.is_file() else None
            for p in reference_inputs(root, width)}


def weighted_distribution(values, weights):
    values, weights = np.asarray(values, dtype=float), np.asarray(weights, dtype=float)
    if values.size == 0 or values.shape != weights.shape:
        raise ValueError("Distribution needs nonempty matching values and weights")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("Distribution needs finite values and positive area weights")
    order = np.argsort(values.ravel(), kind="stable")
    sorted_values, sorted_weights = values.ravel()[order], weights.ravel()[order]
    cumulative = np.cumsum(sorted_weights)
    cumulative /= cumulative[-1]

    def quantile(q):
        return sorted_values[np.minimum(np.searchsorted(cumulative, q), len(order) - 1)]

    statistics = dict(zip(("min", "p01", "median", "p95", "p99", "max"),
                          quantile([0, .01, .5, .95, .99, 1]).tolist()),
                      mean=float(np.average(values, weights=weights)))
    return statistics, {"fraction_below": QUANTILES.tolist(),
                        "elevation_m": quantile(QUANTILES).tolist()}


def measure_earth(elevation, non_ocean):
    elevation, non_ocean = np.asarray(elevation), np.asarray(non_ocean, dtype=bool)
    if elevation.ndim != 2 or elevation.shape != non_ocean.shape or elevation.shape[1] != 2 * elevation.shape[0]:
        raise ValueError("Earth reference must share a cell-centred global 2:1 grid")
    height = elevation.shape[0]
    # Exact spherical row areas; longitude width and R^2 cancel when normalized.
    edges = np.linspace(np.pi / 2, -np.pi / 2, height + 1)
    area = np.broadcast_to((np.sin(edges[:-1]) - np.sin(edges[1:]))[:, None], elevation.shape)
    metrics, curve = weighted_distribution(elevation, area)
    return {"elevation_m": metrics, "hypsography": curve,
            "land_fraction": float(np.sum(area[non_ocean]) / area.sum()),
            "land_elevation_m": weighted_distribution(elevation[non_ocean], area[non_ocean])[0],
            "ocean_depth_m": weighted_distribution(-elevation[~non_ocean], area[~non_ocean])[0],
            "area_above_9000_m": float(np.sum(area[elevation > 9000]) / area.sum())}


def measure_earth_reference(root, width):
    root = Path(root)
    heightmap, mask, *_ = reference_inputs(root, width)
    with Image.open(heightmap) as image:
        if image.mode != "F" or image.size != (width, width // 2):
            raise ValueError("Earth heightmap must be the signed float32 metre reference")
        # Prepared GeoTIFF: global north-up cell edges, not inclusive pole pixels.
        if not np.allclose(image.tag_v2.get(33550, ()), (360 / width, 360 / width, 0)):
            raise ValueError("Unexpected Earth pixel scale")
        if not np.allclose(image.tag_v2.get(33922, ()), (0, 0, 0, -180, 90, 0)):
            raise ValueError("Unexpected Earth grid origin")
        elevation = np.array(image)
    with Image.open(mask) as image:
        if image.mode != "L" or image.size != (width, width // 2):
            raise ValueError("Earth non-ocean mask must match the heightmap")
        mask_values = np.array(image)
    if not np.isin(mask_values, [0, 255]).all():
        raise ValueError("Earth non-ocean mask must be binary")
    metrics = measure_earth(elevation, mask_values == 255)
    return {"title": f"Earth · ETOPO 2022 · {width} × {width // 2}", "metrics": metrics,
            "shape": list(elevation.shape), "input_sha256": input_hashes(root, width),
            "datum": "EGM2008; signed metres, ice surface", "weighting": "spherical cell area",
            "mask": "GSHHG non-ocean mask; includes lakes and below-sea-level land",
            "notes": "Distribution reference only, not a geographical reconstruction. Earth uses spherical area; native tectonics uses equal torus cells. Reference averaging smooths peaks and trenches; extrema depend on resolution. The local processing manifest is unavailable; current input and verification hashes are recorded."}


def load_earth_reference(root, width=2048, *, refresh=False):
    """Measure once, then read cached statistics; hashes detect changed inputs.

    Hash verification reads bytes but does not decode or sort the source grids.
    Changes to the measurement definition require a new CACHE_SCHEMA version.
    """
    root = Path(root)
    cache = root / f"refs/processed/metadata/{width}/tectonics-hypsography.json"
    if cache.is_file() and not refresh:
        document = json.loads(cache.read_text(encoding="utf-8"))
        reference = document.get("reference", {})
        if (document.get("schema") != CACHE_SCHEMA or reference.get("shape") != [width // 2, width]
                or reference.get("input_sha256") != input_hashes(root, width)):
            raise ValueError("Stale Earth hypsography cache; run tectonics_earth_reference.py "
                             f"--width {width} --refresh")
    else:
        reference = measure_earth_reference(root, width)
        cache.parent.mkdir(parents=True, exist_ok=True)
        temporary = cache.with_suffix(".json.tmp")
        temporary.write_text(json.dumps({"schema": CACHE_SCHEMA, "reference": reference},
                                       indent=2, allow_nan=False) + "\n", encoding="utf-8")
        temporary.replace(cache)
    # Cheap derived curve checks use only the cached 201 points.
    reference["metrics"]["hypsography_shape"] = hypsography_shape(reference["metrics"]["hypsography"])
    return {**reference, "cache": {"path": cache.relative_to(root).as_posix(),
                                   "sha256": digest(cache), "schema": CACHE_SCHEMA}}


def compare_to_earth(metrics, reference):
    earth = reference["metrics"]
    delta = np.asarray(metrics["hypsography"]["elevation_m"]) - np.asarray(earth["hypsography"]["elevation_m"])
    return {"land_fraction_delta_percentage_points": 100 * (metrics["land_fraction"] - earth["land_fraction"]),
            "hypsography_mean_absolute_quantile_gap_m": float(np.trapezoid(abs(delta), QUANTILES)),
            "mean_elevation_delta_m": metrics["elevation_m"]["mean"] - earth["elevation_m"]["mean"],
            "mean_ocean_depth_delta_m": (metrics["ocean_depth_m"]["mean"] - earth["ocean_depth_m"]["mean"]
                                         if metrics["ocean_depth_m"] else None)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare the reusable Earth hypsography reference")
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--refresh", action="store_true", help="Explicitly replace cached measurements")
    args = parser.parse_args()
    reference = load_earth_reference(Path(__file__).resolve().parents[2], args.width, refresh=args.refresh)
    print(json.dumps(reference["cache"], indent=2))
