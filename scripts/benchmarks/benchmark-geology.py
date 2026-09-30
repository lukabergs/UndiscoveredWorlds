# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.2", "rasterio==1.5.1"]
# ///
"""Compare native tectonic bundles with physical reference distributions (no downloads)."""
import argparse
import hashlib
import html
import json
from pathlib import Path

import numpy as np
import rasterio

REPO = Path(__file__).resolve().parents[2]
REFERENCES = {
    "tectonics/seafloor_age": ("Ma", "Oceanic crust age", "comparable"),
    "tectonics/seafloor_full_rate": ("mm yr-1", "Spreading rate at crust formation", "Native velocity is cells/update; no physical length scale"),
    "crust/solid_crust_thickness": ("km", "Solid crust thickness", "Native crust column is not thickness in kilometres"),
    "geology/marine_sediment_thickness": ("m", "Marine sediment thickness", "No exported sediment column in metres"),
    "geothermal/heat_flow": ("mW m-2", "Observed heat flow", "No exported thermal model; observations have uneven coverage"),
}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def area_weights(height, width):
    # Exact relative cell areas on a regular global latitude-longitude raster.
    edges = np.linspace(np.pi / 2, -np.pi / 2, height + 1)
    rows = np.sin(edges[:-1]) - np.sin(edges[1:])
    return np.broadcast_to(rows[:, None], (height, width))


def distribution(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    if values.size == 0 or np.any(~np.isfinite(values)) or np.any(weights <= 0):
        raise ValueError("Distribution requires finite samples with positive weights")
    order = np.argsort(values, kind="stable")
    values, weights = values[order], weights[order]
    cumulative = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
    quantiles = np.interp([0.05, 0.5, 0.95], cumulative, values)
    return {"cells": int(values.size), "mean": float(np.average(values, weights=weights)),
            "min": float(values[0]), "p05": float(quantiles[0]),
            "median": float(quantiles[1]), "p95": float(quantiles[2]), "max": float(values[-1])}


def age_distribution_distance(first, first_weights, second, second_weights):
    # Exact 1-D Wasserstein distance: integral of absolute weighted CDF difference.
    def cdf_at(samples, weights, edges):
        order = np.argsort(samples, kind="stable")
        samples, weights = samples[order], weights[order]
        cumulative = np.concatenate(([0.0], np.cumsum(weights) / weights.sum()))
        return cumulative[np.searchsorted(samples, edges, side="right")]
    edges = np.unique(np.concatenate((first, second)))
    if edges.size < 2:
        return 0.0
    delta = abs(cdf_at(first, first_weights, edges) - cdf_at(second, second_weights, edges))
    return float(np.sum(delta[:-1] * np.diff(edges)))


def read_reference(root, key, width, units, products):
    if products[key]["units"] != units:
        raise ValueError(f"Unexpected units for {key}: {products[key]['units']}")
    path = root / "fields" / key / f"{width}.tif"
    with rasterio.open(path) as dataset:
        if (dataset.crs != rasterio.crs.CRS.from_epsg(4326) or
                dataset.width != width or dataset.height * 2 != width or
                not np.allclose(tuple(dataset.bounds), (-180, -90, 180, 90)) or
                dataset.transform.e >= 0 or dataset.transform.a <= 0):
            raise ValueError(f"Expected a global north-up 2:1 EPSG:4326 grid: {path}")
        field = dataset.read(1, masked=True)
    valid = ~np.ma.getmaskarray(field) & np.isfinite(field.data)
    weights = area_weights(*field.shape)
    return field.data, valid, weights, {
        "path": str(path.resolve()), "sha256": digest(path), "units": units,
        "dataset": products[key]["dataset"],
        "valid_area_fraction": float(weights[valid].sum() / weights.sum()),
        "statistics_over_valid_cells": distribution(field.data[valid], weights[valid]),
    }


def read_bundle(path):
    manifest_path = path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "tectonic-snapshot-bundle/v1" or manifest.get("contract_schema_version") != 6:
        raise ValueError(f"Unsupported native bundle contract: {path}")
    width, height = int(manifest["width"]), int(manifest["height"])
    if width != height * 2 or height < 1:
        raise ValueError("Reference comparison requires a 2:1 planetary interpretation")
    fields = {}
    hashes = {}
    for name, dtype in {"heightmap": "<f4", "crust_age_myr": "<f4", "crust_class": "u1",
                        "geologic_regime": "u1", "convergence_score": "u1",
                        "divergence_score": "u1", "shear_score": "u1"}.items():
        field_path = (path / manifest["files"][name]).resolve()
        if not field_path.is_relative_to(path.resolve()):
            raise ValueError("Bundle field lies outside the bundle")
        field = np.fromfile(field_path, dtype=dtype)
        if field.size != width * height or not np.isfinite(field).all():
            raise ValueError(f"Invalid shape or nonfinite values: {field_path}")
        fields[name] = field.reshape(height, width)
        hashes[name] = digest(field_path)
    if (np.any(fields["crust_age_myr"] < 0) or np.any(fields["crust_class"] > 3) or
            np.any(fields["geologic_regime"] > 7) or
            any(np.any(fields[name] > 100) for name in ("convergence_score", "divergence_score", "shear_score"))):
        raise ValueError(f"Native field values violate the export contract: {path}")
    return manifest, fields, hashes


def benchmark(bundles, reference_root, width):
    catalog = json.loads((reference_root / "manifest.json").read_text(encoding="utf-8"))
    products = {product["key"]: product for product in catalog["products"]}
    report = {"schema": "geology-reference-benchmark/v1", "reference_width": width,
              "interpretation": "Area-weighted distribution comparison, not geographic reconstruction or a calibrated pass/fail score.",
              "limitations": ["Native geometry is a planar torus; latitude-area weights apply only to its planetary map interpretation.",
                              "Initial ages and motion retain resolution/timestep-dependent procedural assumptions.",
                              "Sparse reference statistics describe observed cells, not an unbiased global average.",
                              "Mineral reserves remain potential scores; lithology/composition/ore quantity have no compatible native exports."],
              "references": {}, "runs": []}
    loaded = {}
    for key, (units, title, status) in REFERENCES.items():
        data, valid, weights, summary = read_reference(reference_root, key, width, units, products)
        summary.update(title=title, comparison_status=status)
        report["references"][key] = summary
        loaded[key] = data, valid, weights
    reference_age, valid, reference_weights = loaded["tectonics/seafloor_age"]
    reference_samples, reference_sample_weights = reference_age[valid], reference_weights[valid]
    for path in bundles:
        manifest, fields, hashes = read_bundle(path)
        weights = area_weights(*fields["heightmap"].shape)
        oceanic = fields["crust_class"] == 1
        seafloor = oceanic & (fields["heightmap"] < 1.0)
        invalid = oceanic & (fields["geologic_regime"] == 2)
        run = {"bundle": str(path.resolve()), "manifest_sha256": digest(path / "manifest.json"),
               "field_sha256": hashes, "seed": manifest["seed"], "width": manifest["width"],
               "height": manifest["height"], "time_myr": manifest["time_myr"],
               "scenario": manifest["run_scenario"],
               "oceanic_crust_area_fraction": float(weights[oceanic].sum() / weights.sum()),
               "oceanic_cells_with_continental_collision": int(invalid.sum()),
               "invariants_pass": not bool(invalid.any())}
        if seafloor.any():
            age = fields["crust_age_myr"][seafloor]
            sample_weights = weights[seafloor]
            run["seafloor_age_ma"] = distribution(age, sample_weights)
            run["age_wasserstein_distance_ma"] = age_distribution_distance(
                age, sample_weights, reference_samples, reference_sample_weights)
        else:
            run["seafloor_age_ma"] = None
            run["age_comparison_status"] = "No exposed oceanic crust; age comparison unavailable"
        report["runs"].append(run)
    return report


def render_html(report):
    esc = html.escape
    rows = []
    for run in report["runs"]:
        age = run["seafloor_age_ma"]
        rows.append("<tr>" + "".join(f"<td>{esc(str(value))}</td>" for value in (
            run["seed"], f"{run['width']} × {run['height']}", run["time_myr"],
            f"{age['mean']:.2f}" if age else "Unavailable",
            f"{run['age_wasserstein_distance_ma']:.2f}" if age else "Unavailable",
            run["oceanic_cells_with_continental_collision"])) + "</tr>")
    reference_rows = []
    for ref in report["references"].values():
        stats = ref["statistics_over_valid_cells"]
        reference_rows.append("<tr>" + "".join(f"<td>{esc(str(value))}</td>" for value in (
            ref["title"], ref["units"], f"{100 * ref['valid_area_fraction']:.1f}%",
            f"{stats['mean']:.2f}", ref["comparison_status"])) + "</tr>")
    return """<!doctype html><html lang="en"><meta charset="utf-8"><title>Geology reference benchmark</title>
<style>body{font:16px system-ui;max-width:1150px;margin:40px auto;padding:0 24px;color:#203040;background:#f8fafc}
table{border-collapse:collapse;width:100%;margin:24px 0}th,td{padding:10px;text-align:left;border-bottom:1px solid #ccd5df}th{background:#e8eef5}li{margin:8px 0}</style>
<h1>Geology reference benchmark</h1><p>""" + esc(report["interpretation"]) + """</p>
<table><tr><th>Seed</th><th>Grid</th><th>Time (Myr)</th><th>Mean seafloor age (Ma)</th><th>Age distribution distance (Ma)</th><th>Oceanic collision errors</th></tr>""" + "".join(rows) + """</table>
<p>Age distance is the area-weighted Wasserstein distance; lower means closer distributions, without implying equivalent tectonic histories. There is no fitted acceptance threshold.</p>
<table><tr><th>Reference</th><th>Units</th><th>Valid area</th><th>Mean of valid cells</th><th>Comparison status</th></tr>""" + "".join(reference_rows) + "</table><ul>" + "".join(
        f"<li>{esc(note)}</li>" for note in report["limitations"]) + "</ul></html>"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, nargs="+", required=True)
    parser.add_argument("--references", type=Path, default=REPO / "refs/processed/physical")
    parser.add_argument("--width", type=int, choices=(128, 256, 512, 1024, 2048), default=512)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = benchmark(args.bundle, args.references, args.width)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    args.output.with_suffix(".html").write_text(render_html(report), encoding="utf-8")
    print(f"Wrote {args.output} and {args.output.with_suffix('.html')}")
    return int(any(not run["invariants_pass"] for run in report["runs"]))


if __name__ == "__main__":
    raise SystemExit(main())
