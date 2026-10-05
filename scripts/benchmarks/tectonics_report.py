"""Build a portable HTML comparison from numbered, verified native tectonics runs."""
import argparse
import base64
import colorsys
import json
import os
from pathlib import Path
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import palettes, runs
from uwtools.hashing import sha256, sha256_bytes
from uwtools.paths import ROOT
import map_rendering
from tectonics_report_data import load_bundle, diagnostics, compare_fields, write_json
from tectonics_earth_reference import load_earth_reference, compare_to_earth
from tectonics_feedback_report import feedback_snapshots, feedback_controls, climate_field_specs

PALETTE_PATH = palettes.TECTONICS_PATH


def outputs(run_id):
    """Maps, fields, animations and metrics derived from one run: never inside runs/archive."""
    return runs.output_directory("tectonics", run_id)


def relative(path):
    path = Path(path)
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def report_template():
    return Path(__file__).with_name("tectonics_report_template.html").read_text(encoding="utf-8").replace(
        "__MAP_ANIMATION_JS__", Path(__file__).with_name("map_animation.js").read_text(encoding="utf-8")).replace(
        "__MAP_COMPARISON_JS__", Path(__file__).with_name("map_comparison.js").read_text(encoding="utf-8")).replace(
        "__MAP_SELECTION_JS__", Path(__file__).with_name("map_report_selection.js").read_text(encoding="utf-8"))


def stops_for(spec, definition):
    if "stops" in spec:
        return spec["stops"]
    low, high = spec["range"]
    return [(low + t * (high - low), rgb) for t, rgb in definition["sequential"]]


def image_uri(path):
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def external_asset(data, suffix, output):
    """Immutable shared assets keep report size linear in recorded snapshots."""
    directory = Path(output).parent / "assets"
    directory.mkdir(parents=True, exist_ok=True)
    checksum = sha256_bytes(data)
    name = checksum + suffix
    target = directory / name
    if target.exists():
        if sha256(target) != checksum:
            raise ValueError(f"Report asset changed: {target}")
    else:
        target.write_bytes(data)
    return "assets/" + name


def externalize_report_assets(path):
    """Keep image bytes immutable but load them on demand instead of parsing all APNGs."""
    path = Path(path)
    original = path.read_text(encoding="utf-8")
    payload, _ = json.JSONDecoder().raw_decode(original, original.index("const DATA=") + len("const DATA="))
    asset_directory = path.parent / "assets"
    asset_directory.mkdir(exist_ok=True)
    assets = set()

    def convert(value):
        if isinstance(value, dict):
            return {key: convert(item) for key, item in value.items()}
        if isinstance(value, list):
            return [convert(item) for item in value]
        if isinstance(value, str) and value.startswith("data:image/png;base64,"):
            pixels = base64.b64decode(value.split(",", 1)[1], validate=True)
            name = sha256_bytes(pixels) + ".png"
            target = asset_directory / name
            if target.exists():
                if sha256(target) != name[:-4]:
                    raise ValueError(f"Report image asset changed: {target}")
            else:
                target.write_bytes(pixels)
            assets.add(name)
            return "assets/" + name
        if isinstance(value, str) and value.startswith("assets/") and value.endswith(".png"):
            name = value.removeprefix("assets/")
            if len(name) != 68 or any(c not in "0123456789abcdef" for c in name[:-4]):
                raise ValueError("Invalid content-addressed report image")
            if sha256(asset_directory / name) != name[:-4]:
                raise ValueError("Report image checksum mismatch")
            assets.add(name)
        return value

    payload = convert(payload)
    template_path = Path(__file__).with_name("tectonics_report_template.html")
    rendered = report_template().replace(
        "__REPORT_DATA__", json.dumps(payload, allow_nan=False).replace("<", "\\u003c"))
    path.write_text(rendered, encoding="utf-8")
    metadata_path = path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    metadata["asset_conversion"] = {"mode": "external-content-addressed-png", "directory": "assets",
        "images": sorted(assets), "externalizer_sha256": sha256(Path(__file__)),
        "template_sha256": sha256(template_path), "animation_player_sha256": sha256(Path(__file__).with_name("map_animation.js")),
        "html_sha256": sha256(path)}
    write_json(metadata_path, metadata)
    print(f"External assets: {path.name}, {len(assets)} images, HTML {len(original):,} -> {len(rendered):,} characters")


def run_animations(run_id, states, fieldspec, image_link=image_uri):
    selected = [state for state in states if state["run"] == run_id]
    times = [state["time"] for state in selected]
    # 120 ms/Myr, plus a brief final pause. Zero-time administrative steps
    # remain visible for 120 ms; these are observations, not invented in-betweens.
    durations = [max(120, round((b - a) * 120)) for a, b in zip(times, times[1:])] + [1200]
    directory = outputs(run_id)
    metadata_path = directory / "animations.json"
    previous = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    result = {}
    for field, spec in fieldspec.items():
        if any(field in state.get("missing_fields", []) for state in selected):
            continue
        paths = [directory / f"maps/{field}/{s['label']}.png" for s in selected]
        # Finalization can share model time and pixels with a preceding
        # observation, even at the same update count. Checkpoint labels prevent
        # APNG coalescing those separately recorded observations.
        labels = [f"Run {run_id} | {s['time']:g} Myr | {s['label']} | {spec['title']} ({spec['units']})" for s in selected]
        output = directory / f"maps/{field}/animation.png"
        signature = sha256_bytes(json.dumps({"frames": [sha256(path) for path in paths],
            "labels": labels, "durations": durations,
            "renderer_sha256": sha256(Path(map_rendering.__file__))}, sort_keys=True).encode())
        metadata = previous.get(field, {})
        if (metadata.get("render_signature") != signature or not output.exists()
                or sha256(output) != metadata.get("sha256")):
            metadata = map_rendering.save_animation(output, paths, labels, durations)
        metadata.update(path=relative(output), sha256=sha256(output),
                        times_myr=times, durations_ms=durations, render_signature=signature)
        result[field] = dict(metadata, image=image_link(output))
    write_json(metadata_path,
               {field: {k: v for k, v in value.items() if k != "image"} for field, value in result.items()})
    return result


def write_run_note(run_id, record, states):
    """Create a sealed run's review stub once; preserve subsequent human findings/decisions.

    A run of the current layout has its run page and ledger draft instead."""
    path = ROOT / f"runs/reports/tectonics/{run_id}.md"
    if path.exists() or runs.directory("tectonics", run_id) is not None:
        return
    final = next(s for s in states if s["run"] == run_id and s["label"] == "final")
    m = final["metrics"]
    peak = m["land_elevation_m"]["max"] if m["land_elevation_m"] else "No land"
    text = f"""# Tectonics run {run_id}

[Visual comparison]({run_id}.html)

## Purpose

Inspect the configured procedural tectonics scenario against the established fixed-seed series. Visual acceptance is pending human review.

## Recorded observations

- Seed: {record['recipe']['seed']}; grid: {final['shape'][1]} × {final['shape'][0]}.
- Final state: {final['time']} Myr, {final['updates']} update calls.
- Land fraction: {m['land_fraction']:.3%}; maximum land elevation: {peak} m.
- Native field invariants: {m['invariants_pass']}; codec range preserved: {m['codec_range_preserved']}.
- Exact final replay: {record.get('replay', {}).get('exact_match', 'not checked')}.
- Source, executable, recipe and field hashes: `runs/archive/manifests/tectonics/{run_id}.json`.
- Numerical metrics: `{relative(outputs(run_id) / 'metrics.json')}`.
- Raw execution output: `runs/archive/logs/tectonics/{run_id}.log`.

## Interpretation and decision

Record visual findings and the next controlled experiment here. Numerical correctness does not establish physical validity. Use the explicit column budget, including sources and sinks, when available; crust-thickness proxy sums are not a mass budget.
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def findings(run_id, output):
    """Link from the report to the run's findings: its run page, or a sealed run's note."""
    run = runs.directory("tectonics", run_id)
    note = run / "evaluation/index.html" if run else ROOT / f"runs/reports/tectonics/{run_id}.md"
    return Path(os.path.relpath(note, Path(output).absolute().parent)).as_posix()


def report_snapshots(run_id, record):
    if record.get('kind') == 'tectonics-climate-feedback/v1':
        yield from feedback_snapshots(record)
        return
    for checkpoint in record['snapshots']:
        label = checkpoint['label']
        manifest, fields, hashes = load_bundle(checkpoint['path'], record.get('field_sha256', {}).get(label))
        metrics = diagnostics(manifest, fields)
        metrics['fraction_below'] = {str(level): float((fields['elevation'] < level).mean())
            for level in (-6000, -4000, -2000, -1000, -258, -200, -160, 0, 1000, 2000, 4000)}
        metrics['shoreline_within_1m'] = float((abs(fields['elevation']) <= 1).mean())
        yield dict(manifest, label=label), fields, hashes, metrics


def build_report(run_ids, output, earth_reference=None, maps_only=False, external_assets=False,
                 baseline_run=None, candidate_run=None):
    if any(run is not None and run not in run_ids for run in (baseline_run, candidate_run)):
        raise ValueError('Preferred A/B runs must be included in the report')
    if earth_reference is not None:
        # Derive band estimates from the immutable cached 0.5%-quantile curve.
        # Do not remeasure or silently rewrite the Earth reference cache.
        curve = earth_reference['metrics']['hypsography']
        cdf = lambda z: float(np.interp(z, curve['elevation_m'], curve['fraction_below']))
        cache_path = ROOT / earth_reference['cache']['path']
        thresholds_path = cache_path.with_name('tectonics-shelf-thresholds.json')
        thresholds = None
        if thresholds_path.exists():
            thresholds = json.loads(thresholds_path.read_text())
            input_path = f"refs/processed/heightmap/{earth_reference['shape'][1]}.tif"
            if thresholds['input_sha256'] != earth_reference['input_sha256'][input_path]:
                raise ValueError('Earth shelf thresholds refer to a different heightmap')
            cdf = lambda z: thresholds['fraction_below'][str(z)]
        earth_reference = dict(earth_reference, metrics=dict(earth_reference['metrics'], elevation_area_bands={
            'shelf_0_200_m': cdf(0)-cdf(-200),
            'upper_ocean_200_2000_m': cdf(-200)-cdf(-2000),
            'deep_ocean_below_6000_m': cdf(-6000),
            'lowland_0_1000_m': cdf(1000)-cdf(0),
            'highland_above_2000_m': 1-cdf(2000)}))
        if thresholds:
            earth_reference['metrics']['fraction_below'] = thresholds['fraction_below']
            earth_reference['threshold_cache'] = dict(path=thresholds_path.relative_to(ROOT).as_posix(), sha256=sha256(thresholds_path))
    definition = palettes.tectonics()
    fieldspec = definition["fields"]
    records = [runs.tectonics(run_id) for run_id in run_ids]
    if not records or any(r["status"] != "complete" for r in records):
        raise ValueError("Report requires complete runs")
    if any(r.get('kind') == 'tectonics-climate-feedback/v1' for r in records):
        if not external_assets and not maps_only:
            raise ValueError('Feedback comparisons require --external-assets for heterogeneous fields')
        fieldspec.update(climate_field_specs())
    seeds = {str(r["id"]): r["recipe"]["seed"] for r in records}
    states, arrays, animations, feedback = [], [], {}, {}
    image_link = (lambda path: external_asset(path.read_bytes(), ".png", output)) if external_assets else image_uri
    # Stable ID colours across every checkpoint/run, not a random per-image palette.
    id_colours = {i: [round(c * 255) for c in colorsys.hsv_to_rgb((i * .61803398875) % 1, .52, .92)]
                  for i in range(256)}
    for run_id, record in zip(run_ids, records):
        print(f'Rendering run {run_id}…', flush=True)
        for manifest, fields, hashes, metrics in report_snapshots(run_id, record):
            label = manifest['label']
            if earth_reference is not None:
                metrics["earth_comparison"] = compare_to_earth(metrics, earth_reference)
            state = {"run": run_id, "label": label, "time": manifest["time_myr"],
                     "seed": record["recipe"]["seed"], "runtime_seconds": record.get("elapsed_seconds"),
                     "surface_model": manifest.get("run_scenario", {}).get("surface_model", "unrecorded"),
                     "updates": manifest["update_count"],
                     "shape": [manifest["height"], manifest["width"]], "metrics": metrics,
                     "phase": manifest.get('phase'), "boundary_note": manifest.get('boundary_note'),
                     "maps": {}, "numeric": {}, "display_saturation": {}, "missing_fields": []}
            for field, spec in fieldspec.items():
                if field not in fields:
                    state["missing_fields"].append(field)
                    continue
                data = fields[field]
                map_path = outputs(run_id) / f"maps/{field}/{label}.png"
                numeric = outputs(run_id) / f"fields/{field}/{label}.npy"
                numeric.parent.mkdir(parents=True, exist_ok=True)
                np.save(numeric, data, allow_pickle=False)
                if external_assets and not maps_only:
                    state["numeric"][field] = {
                        "url": external_asset(np.asarray(data, dtype="<f8").tobytes(), ".f64", output),
                        "cells": int(data.size), "categorical": bool(spec.get("ids") or "colours" in spec)}
                if spec.get("ids"):
                    for identity in np.unique(data):
                        if int(identity) not in id_colours:
                            id_colours[int(identity)] = [round(c * 255) for c in
                                colorsys.hsv_to_rgb((int(identity) * .61803398875) % 1, .52, .92)]
                    rgb = map_rendering.categories(data, id_colours)
                elif "colours" in spec:
                    rgb = map_rendering.categories(data, dict(enumerate(spec["colours"])))
                else:
                    stops = stops_for(spec, definition)
                    rgb = map_rendering.colourize(data, stops, definition["missing"])
                    state["display_saturation"][field] = {"below": int((data < stops[0][0]).sum()),
                                                          "above": int((data > stops[-1][0]).sum())}
                map_rendering.save_rgb(map_path, rgb)
                state["maps"][field] = image_link(map_path)
            boundary = fields.get('_boundary_overlay', fields.get('boundary_type'))
            rgba = np.zeros((*boundary.shape, 4), dtype=np.uint8)
            rgba[..., :3] = map_rendering.categories(boundary,
                                                    dict(enumerate(fieldspec["boundary_type"]["colours"])))
            rgba[..., 3] = np.where(boundary > 0, 230, 0)
            overlay = outputs(run_id) / f"maps/boundary_overlay/{label}.png"
            overlay.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(rgba).save(overlay)
            state["overlay"] = image_link(overlay)
            coastline = outputs(run_id) / f"maps/coastline_overlay/{label}.png"
            coastline.parent.mkdir(parents=True, exist_ok=True)
            Image.fromarray(map_rendering.outline(fields["elevation"] >= 0)).save(coastline)
            state["coastline"] = image_link(coastline)
            state["hashes"] = hashes
            states.append(state)
            if not external_assets and not maps_only:
                arrays.append(fields)
        write_json(outputs(run_id) / "metrics.json",
                   {state["label"]: state["metrics"] for state in states if state["run"] == run_id})
        animations[str(run_id)] = run_animations(run_id, states, fieldspec, image_link)
        if record.get('kind') == 'tectonics-climate-feedback/v1':
            control_fields, review = feedback_controls(record)
            controls = {}
            difference_stops = [[-1500, [49, 87, 154]], [0, [229, 229, 219]], [1500, [170, 66, 40]]]
            for key, title in [('initial', 'Inherited run 62 terrain'), ('final', 'Spatial climate feedback'),
                               ('uniform_control', 'Uniform forcing control'), ('matched_control', 'Matched mean forcing control'),
                               ('total_change', 'Spatial feedback minus initial'),
                               ('climate_effect', 'Spatial feedback minus uniform control'),
                               ('spatial_effect', 'Spatial feedback minus matched control')]:
                data = control_fields[key]
                is_difference = key not in review['distributions']
                stops = difference_stops if is_difference else fieldspec['elevation']['stops']
                map_path = outputs(run_id) / f'maps/feedback/{key}.png'
                map_rendering.save_rgb(map_path, map_rendering.colourize(data, stops, definition['missing']))
                controls[key] = dict(title=title, image=image_link(map_path), stops=stops,
                    units='m difference' if is_difference else 'm relative to sea level',
                    below=int((data < stops[0][0]).sum()), above=int((data > stops[-1][0]).sum()),
                    numeric=external_asset(np.asarray(data, dtype='<f8').tobytes(), '.f64', output) if external_assets else None)
            feedback[str(run_id)] = dict(parent=record['parent'], controls=controls,
                duration_myr=record['surface_experiment_duration_myr'], review=review,
                assessment=record.get('physical_assessment'), geometry=record['geometry'], limitations=record['limitations'],
                matched_control_seconds=record.get('matched_control', {}).get('elapsed_seconds'))
        else:
            write_run_note(run_id, record, states)
    if maps_only:
        print(f"Exported maps, numerical fields and APNGs for runs {run_ids}; no HTML generated")
        return states
    # Preserve numeric comparisons, independent of palette saturation and PNG precision.
    comparisons = {}
    for i, first in enumerate(arrays):
        for j, second in enumerate(arrays):
            if states[i]["seed"] == states[j]["seed"] and first["heightmap"].shape == second["heightmap"].shape:
                comparisons[f"{i}:{j}"] = compare_fields(first, second)
    probe = None
    probe_path = runs.derived("tectonics", run_ids[-1], "surface-process-probe.json")
    if probe_path:
        probe = {"data": json.loads(probe_path.read_text()), "sha256": sha256(probe_path),
                 "path": relative(probe_path)}
    regional_reviews = {}
    for run_id in run_ids:
        regional_path = runs.derived("tectonics", run_id, "coastal-window.png")
        if regional_path:
            regional_reviews[str(run_id)] = image_link(regional_path)
    payload = {"states": states, "fields": fieldspec, "sequential": definition["sequential"],
               "regional_reviews": regional_reviews,
               "preferred_runs": [baseline_run or run_ids[0], candidate_run or run_ids[-1]]
                   if baseline_run is not None or candidate_run is not None else None,
               "comparisons": comparisons, "seed": records[0]["recipe"]["seed"], "seeds": seeds,
               "same_seed": len(set(seeds.values())) == 1,
               "run_ids": run_ids, "animations": animations, "earth_reference": earth_reference,
               "feedback": feedback,
               "run_summaries": {str(run_id): {'name': record['recipe'].get('name'),
                    'seconds': record.get('elapsed_seconds'), 'kind': record.get('kind', 'native tectonics')}
                    for run_id, record in zip(run_ids, records)},
               "recipes": {str(run_id): record["recipe"] for run_id, record in zip(run_ids, records)},
               "findings_url": Path(output).with_suffix(".md").name if Path(output).with_suffix(".md").exists() else findings(run_ids[-1], output),
               "surface_probe": probe,
               "replay_exact": all(r.get("replay", {}).get("exact_match", False) for r in records)}
    template = report_template()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(template.replace("__REPORT_DATA__", json.dumps(payload, allow_nan=False).replace("<", "\\u003c")), encoding="utf-8")
    write_json(output.with_suffix(".json"), {
        "schema": "tectonics-visual-report/v1", "runs": run_ids,
        "renderer": map_rendering.VERSION, "palette_sha256": sha256(PALETTE_PATH),
        "source_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(__file__).with_name("map_rendering.py"),
                          Path(__file__).with_name("tectonics_report_data.py"), Path(__file__).with_name("tectonics_earth_reference.py"),
                          Path(__file__).with_name("map_animation.js"), Path(__file__).with_name("map_comparison.js"),
                          Path(__file__).with_name("map_report_selection.js"),
                          Path(__file__).with_name("tectonics_feedback_report.py"),
                          Path(__file__).with_name("tectonics_coastal_review.py"),
                          Path(__file__).with_name("tectonics_boundary_review.py"),
                          Path(__file__).with_name("tectonics_report_template.html"))},
        "states": [{k: v for k, v in state.items() if k not in ("maps", "overlay", "coastline")} for state in states],
        "comparisons": comparisons, "earth_reference": earth_reference, "surface_probe": probe,
        "feedback": feedback,
    })
    print(f"Wrote {output}")
    return states


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--runs", type=int, nargs="+")
    action.add_argument("--externalize", type=Path, nargs="+", help="Convert existing reports without rerunning/re-rendering simulations")
    parser.add_argument("--external-assets", action="store_true", help="Load PNG/APNG assets on demand; keep the sibling assets directory")
    parser.add_argument("--maps-only", action="store_true", help="Export numerical fields, maps and APNGs without an HTML comparison")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--baseline-run", type=int, help="Initially selected A run; does not alter the registry baseline")
    parser.add_argument("--candidate-run", type=int, help="Initially selected B run")
    parser.add_argument("--earth-width", type=int, default=2048,
                        help="Prepared refs heightmap width; 0 explicitly omits Earth")
    args = parser.parse_args()
    if args.externalize:
        for path in args.externalize:
            externalize_report_assets(path)
        return
    if len(set(args.runs)) != len(args.runs):
        parser.error("Run IDs must be unique")
    output = args.output or ROOT / "runs/reports/tectonics" / ("-".join(map(str, args.runs)) + ".html")
    reference = load_earth_reference(ROOT, args.earth_width) if args.earth_width else None
    build_report(args.runs, output, reference, args.maps_only, args.external_assets,
                 args.baseline_run, args.candidate_run)
    if args.external_assets and not args.maps_only:
        externalize_report_assets(output)


if __name__ == "__main__":
    main()
