"""Create an offline side-by-side map gallery for one completed benchmark."""

import argparse
import csv
import html
import json
import math
import os
import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/benchmarks"))
from climate_map_geometry import expected_map_size, load_wind_geometry


def recorded_metric(record, key):
    """Missing support and legacy zero-error placeholders are not perfect scores."""
    value = record.get("metrics", {}).get(key)
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    if "spatial_" in key and record.get("metrics", {}).get("spatial_compared_cells", 0) <= 0:
        return None
    if ("error" in key or "rmse" in key) and value <= 0:
        return None
    if "correlation" in key and not -1 <= value <= 1:
        return None
    return value


def annual_temperature_rmse(run_id):
    path = ROOT / f"runs/diagnostics/climate/{run_id}/monthly_climate_reference_comparison.csv"
    if not path.exists():
        path = ROOT / f"runs/diagnostics/monthly_climate_reference_comparison/{run_id}.csv"
    if path.exists():
        with path.open() as stream:
            row = next((r for r in csv.DictReader(stream) if r["variable"] == "temperature_c" and r["period"] == "annual_mean"), None)
        if row and int(row.get("compared_cells", 0)) > 0:
            value = float(row["area_weighted_rmse"])
            if math.isfinite(value) and value > 0:
                return value
    return None


def record_leader(records, key, higher_is_better, scale=1):
    if higher_is_better is None:
        return "No monotonic target"
    values = [(r["id"], annual_temperature_rmse(r["id"]) if key == "annual_temperature_rmse" else recorded_metric(r, key)) for r in records]
    values = [(i, v) for i, v in values if v is not None]
    if not values:
        return "unavailable"
    best = (max if higher_is_better else min)(v for _, v in values)
    widths = {r["id"]: r.get("horizontal_resolution") for r in records}
    leaders = ", ".join(f"{i} ({widths[i]}-wide output)" if widths[i] else f"{i} (output width unavailable)"
                        for i, v in values if abs(v - best) <= 1e-9)
    return f"{best * scale:.3f} · runs {leaders}"


def baseline_comparison(records, run_id, compare_run_id=145, historical_run_id=145):
    ids = list(dict.fromkeys((142, 143, 144, 145, historical_run_id, compare_run_id, run_id)))
    selected = {r["id"]: r for r in records if r["id"] in ids}
    metrics = [("Köppen distribution weighted relative error ↓", "weighted_relative_error", 1, False),
               ("Land rainfall correlation ↑", "area_weighted_land_precipitation_correlation", 1, True),
               ("Land below 1 mm/month mean (12 mm/year), %", "area_weighted_land_below_1mm_fraction", 100, None),
               ("Surface eastward wind correlation ↑", "era5_surface_eastward_wind_correlation", 1, True),
               ("Surface northward wind correlation ↑", "era5_surface_northward_wind_correlation", 1, True),
               ("Column water correlation ↑", "era5_column_water_correlation", 1, True),
               ("Pressure correlation ↑", "era5_global_pressure_correlation", 1, True),
               ("Land annual temperature RMSE (°C) ↓", "annual_temperature_rmse", 1, False)]
    rows = []
    for label, key, scale, higher_is_better in metrics:
        values = [annual_temperature_rmse(i) if key == "annual_temperature_rmse" else recorded_metric(selected.get(i, {}), key) for i in ids]
        rows.append([label] + [f"{v * scale:.3f}" if v is not None else "unavailable" for v in values] +
                    [record_leader(records, key, higher_is_better, scale)])
    heading = ["Recorded metric"] + [f"Run {i}" for i in ids] + ["Best recorded per metric (all available runs)"]
    markdown = "\n".join("| " + " | ".join(row) + " |" for row in [heading, ["---"] * len(heading), *rows])
    table = "<table>" + "".join("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in row) + "</tr>" for row in [heading, *rows]) + "</table>"
    return table, markdown


def map_panel(path, output, label, title):
    if not path.exists():
        return f'<figure class="missing">{html.escape(label)}: map unavailable.</figure>'
    with path.open("rb") as stream:
        header = stream.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Invalid PNG: {path}")
    width, height = struct.unpack(">II", header[16:24])
    link = relative(path, output)
    return f'<figure><figcaption>{html.escape(label)} · {width}×{height} original pixels</figcaption><a href="{link}"><img loading="lazy" width="{width}" height="{height}" src="{link}" alt="{html.escape(label)}: {title}"></a></figure>'


def relative(path, directory):
    return Path(os.path.relpath(path, directory)).as_posix()


def write_report(run_id, elapsed_seconds=None, compare_run_id=145, historical_run_id=145):
    runs = json.loads((ROOT / "runs/registry/climate.json").read_text(encoding="utf-8"))["runs"]
    record = next(r for r in runs if r.get("id", r.get("run_id")) == run_id)
    width = int(record["horizontal_resolution"])
    output = ROOT / "runs/reports"
    output.mkdir(parents=True, exist_ok=True)
    logpath = ROOT / f"runs/logs/climate/{run_id}.log"
    log = logpath.read_text(encoding="utf-8", errors="replace") if logpath.exists() else ""
    coupling = [line for line in log.splitlines() if line.startswith("Climate coupling iteration=")]
    years = [dict(re.findall(r"(\w+)=([0-9eE.+-]+)", line)) for line in log.splitlines()
             if line.startswith("Climate hydrology spinup cycle=")]
    hydrology_seconds = sum(float(year.get("elapsed_seconds", 0)) for year in years)
    timing = f"Hydrology: {len(years)} simulated years, {hydrology_seconds:.1f} seconds."
    if years and all("transport_seconds" in year for year in years):
        transport_seconds = sum(float(year["transport_seconds"]) for year in years)
        weather_seconds = sum(float(year["weather_seconds"]) for year in years)
        calls = sum(int(year["transport_calls"]) for year in years)
        timing += f" Within those years: transport {transport_seconds:.1f} s ({calls:,} calls), weather {weather_seconds:.1f} s."
    grids = next((line for line in log.splitlines() if line.startswith("Climate hydrology grid output=")), "Grid metadata unavailable.")
    status = coupling[-1] if coupling else "Coupling status unavailable: inspect the run diagnostics."
    baseline_table, baseline_markdown = baseline_comparison(runs, run_id, compare_run_id, historical_run_id)
    if coupling:
        converged = re.search(r"\bconverged=(\d+)", status)
        status_label = "Converged" if converged and converged[1] == "1" else "Not converged"
    else:
        status_label = "Convergence unverified"
    mapsroot = ROOT / "runs/maps/climate"
    refroot = ROOT / "refs/processed/climate/maps/climate"
    wind_geometry = load_wind_geometry(ROOT / f"runs/diagnostics/climate/{run_id}", width)
    dimensions = set()
    cards, links = [], []
    matches = 0
    for path in sorted(mapsroot.rglob(f"{run_id}.png")):
        with path.open("rb") as stream:
            header = stream.read(24)
        map_width, map_height = expected_map_size(path.relative_to(mapsroot), width, wind_geometry)
        if header[:8] != b"\x89PNG\r\n\x1a\n" or struct.unpack(">II", header[16:24]) != (map_width, map_height):
            raise ValueError(f"Unexpected map dimensions: {path}")
        kind = path.parent.relative_to(mapsroot).as_posix()
        dimensions.add((map_width, map_height))
        reference = refroot / kind / f"{map_width}.png"
        runlink = relative(path, output)
        season = next((p for p in kind.split("/") if p in ("1", "2", "3", "4")), "annual")
        title = html.escape(kind)
        runimage = f'<figure><figcaption>Run {run_id} · {map_width}×{map_height}</figcaption><a href="{runlink}"><img loading="lazy" width="{map_width}" height="{map_height}" src="{runlink}" alt="Run {run_id}: {title}"></a></figure>'
        if reference.exists():
            reflink = relative(reference, output)
            referenceimage = f'<figure><figcaption>Reference · {map_width}×{map_height}</figcaption><a href="{reflink}"><img loading="lazy" width="{map_width}" height="{map_height}" src="{reflink}" alt="Reference: {title}"></a></figure>'
            matches += 1
        else:
            referenceimage = '<figure class="missing">No corresponding standalone reference. See the reference guide for unavailable quantities.</figure>'
        comparisons = dict.fromkeys((compare_run_id, historical_run_id))
        historyimage = "".join(map_panel(path.parent / f"{i}.png", output,
            f"{'Previous comparison' if i == compare_run_id else 'Historical'} run {i}", title)
            for i in comparisons if i != run_id)
        cards.append(f'<section data-season="{season}"><h2>{title}</h2><div class="pair">{runimage}{referenceimage}{historyimage}</div></section>')
        links.append(f"- [{kind}]({runlink})")
    if not cards:
        print(f"Run {run_id} has no PNG exports; report includes diagnostics only.")
    dimensions_label = ", ".join(f"{w}×{h}" for w, h in sorted(dimensions)) or "no images"
    geometry_label = f"world {width}×{width//2}"
    wind_palette_note = ""
    if wind_geometry is not None:
        geometry_label += f" · wind images {wind_geometry['columns']}×{wind_geometry['rows']}"
        wind_palette_note = ("Native wind speed scales: surface 0–25 m/s, upper 0–40 m/s. "
                             "Prepared reference speed maps use 0–25 m/s in both layers; upper speed colours are not directly matched. "
                             "Native fields exclude experimental W2/Q3 diagnostic refinement.")
    elapsed = f"{elapsed_seconds/60:.1f} minutes" if elapsed_seconds is not None else "not recorded"
    refindex = relative(ROOT / "refs/processed/climate/maps/index.html", output)
    guide = relative(ROOT / "runs/maps/climate/README.md", output)
    diagnostics = relative(ROOT / f"runs/diagnostics/climate/{run_id}/run_manifest.txt", output)
    information = str(record.get("information", record.get("INFORMATION", "")))
    reference_gate = ""
    checks_path = output / f"{run_id}-checks.json"
    if checks_path.exists():
        comparison = json.loads(checks_path.read_text()).get("reference_regression")
        if comparison and comparison["baseline_run_id"] == compare_run_id:
            failed = [key for key, value in comparison["metrics"].items() if value["regressed"]]
            reference_gate = (f"Reference regression gate versus run {compare_run_id}: " +
                (f"FAILED ({len(failed)} targets deteriorated)." if failed else "No deterioration in the checked targets.") +
                " Numerical validity does not imply accurate climate. See the numerical/reference checks for each target.")
    checks_html = (f'<p><b>{html.escape(reference_gate)}</b> <a href="{run_id}-checks.json">Numerical/reference checks</a></p>'
                   if checks_path.exists() else "")
    checks_markdown = (f"{reference_gate or 'Reference regression gate not run.'} [Numerical/reference checks]({run_id}-checks.json)."
                       if checks_path.exists() else "Numerical/reference checks have not been run.")
    ocean_path = ROOT / f"runs/diagnostics/climate/{run_id}/maps/ocean_fields.csv"
    with ocean_path.open() if ocean_path.exists() else open(os.devnull) as stream:
        mean_sst = "mean_sst_c" in stream.readline()
    ocean_note = ("SST displays the quarter mean used by climate coupling. The companion checker compares day-weighted reference quarters."
                  if mean_sst else "Historical SST displays an end-of-quarter state, which differs from the mean used by climate coupling.")
    spatial_note = ("Native Koeppen spatial scoring compared no cells; the gallery uses the prepared reference maps."
                    if record.get("metrics", {}).get("spatial_compared_cells") == 0 else "")
    historical_note = ("Runs 142–145 have identical recorded metric dictionaries and represent one historical result. "
        "Their maps retain the original 512×257 grid; current maps use cell-centred 2:1 grids. "
        "Run 145 used seed 1 and 128×65 hydrology, with three saved spin-up years and no repeated outer coupling loop; "
        "its convergence test accepted 10.8% snow-mass drift. The current run's grids and total workload are reported above. "
        "Seeds, forcing inputs, grid registration and workloads differ; historical comparisons and record leaders are informational, not regression gates. "
        "Leaders are separate targets across all available records, not a single best simulation. Missing comparisons and unsupported zero-error placeholders are excluded. "
        "Commit 10b4aba contains climate source and data/climate/benchmark_runs.json through run 145; an exact build-to-run mapping is not recorded.")
    document = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Climate run {run_id}</title>
<style>body{{font:16px/1.5 system-ui;background:#f6f8fb;color:#202630;max-width:1120px;margin:2rem auto;padding:0 20px}}h2{{font-size:18px}}section{{margin:28px 0;background:white;padding:12px;border:1px solid #dce1e8}}.pair{{display:grid;grid-template-columns:1fr 1fr;gap:16px}}figure{{margin:0;min-width:0}}img{{width:100%;height:auto;background:black}}.missing{{display:grid;place-content:center;padding:20px;background:#f1f3f7;color:#586376}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf1f6;padding:12px}}.filters{{position:sticky;top:0;background:#eaf0f8;padding:12px}}input,select{{font:inherit;padding:6px;margin-right:12px}}[hidden]{{display:none!important}}a{{color:#126497}}</style>
<style>.pair{{grid-template-columns:repeat(auto-fit,minmax(240px,1fr))}}.table-scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;background:white}}td{{padding:6px 10px;border:1px solid #dce1e8}}tr:first-child{{font-weight:bold}}</style>
<h1>Climate run {run_id} · {geometry_label}</h1><p>{html.escape(information)}</p>
<p><b>{status_label}.</b> Elapsed: {elapsed}. {len(cards)} PNG maps; {matches} matching reference previews. PNG dimensions verified.</p>
{checks_html}
<p>{html.escape(timing)}</p><p>{html.escape(grids)}</p>
<p>{spatial_note}</p>
<p><a href="{guide}">Benchmark commands</a> · <a href="{refindex}">Reference index and palette legends</a> · <a href="climate-history-140-160.md">Earlier simulator history</a> · <a href="{diagnostics}">Run manifest</a> · <a href="{relative(logpath, output)}">Full log</a></p>
<details><summary>Final coupling diagnostics</summary><pre>{html.escape(status)}</pre></details>
<h2>Previous and historical comparisons</h2><div class="table-scroll">{baseline_table}</div><p>{html.escape(historical_note)}</p>
<p>{ocean_note} Ocean currents and SST are separate validation targets; the minimal maps below include both.</p>
<p>1 = January, 2 = April, 3 = July, 4 = October. Wind visualizations are representative months; process maps average their starting quarters.
Reference and simulator LIC/particle renderers have different noise and rasterization. Compare patterns and physical fields, not exact texture pixels.
The retained reference Koeppen classifier uses a −3 °C C/D boundary; the simulator uses 0 °C. SST compares a simulated mixed layer with NOAA OISST.</p>
<p>SST palettes differ: native run −35 °C is blue (#19aff5), 0 °C near-black (#121419), +35 °C orange (#f57823), with land black. Consult the reference index for the OISST palette; equal colours are not equal temperatures across these SST panels.</p>
<p>{wind_palette_note}</p>
<div class="filters"><label>Find <input id="search" placeholder="wind/lic, rain, sea_temp…"></label><label>Month <select id="season"><option value="">All</option><option value="annual">Annual</option><option value="1">1 · January</option><option value="2">2 · April</option><option value="3">3 · July</option><option value="4">4 · October</option></select></label></div>
{"".join(cards)}<script>function filter(){{const q=document.querySelector('#search').value.toLowerCase(),s=document.querySelector('#season').value;document.querySelectorAll('section').forEach(r=>r.hidden=!r.querySelector('h2').textContent.toLowerCase().includes(q)||(s&&r.dataset.season!==s));}}document.querySelectorAll('input,select').forEach(e=>e.addEventListener('input',filter));</script></html>'''
    (output / f"{run_id}.html").write_text(document, encoding="utf-8")
    summary = f"""# Climate benchmark {run_id}

[Comparison gallery: {geometry_label}]({run_id}.html) · [commands]({guide}) · [reference palettes]({refindex}) · [earlier simulator history](climate-history-140-160.md)

- Configuration: {information}
- Elapsed: {elapsed}.
- {timing}
- Grid: {grids}
- Status: **{status_label}**.
- {checks_markdown}
- Artifacts: {len(cards)} PNG maps, verified product dimensions {dimensions_label}; {matches} corresponding reference previews.
{('- ' + wind_palette_note) if wind_palette_note else ''}
- [Log]({relative(logpath, output)}); [run manifest]({diagnostics}).

{baseline_markdown}

{historical_note} {ocean_note}

```text
{status}
```

References retain their original pixels and source conventions. The observed Koeppen C/D threshold is −3 °C, versus 0 °C in the current simulator. LIC/particles have independent rendering noise. Unconverged output is a diagnostic result, not an equilibrated climate.

{spatial_note}

## Maps

""" + "\n".join(links) + "\n"
    (output / f"{run_id}.md").write_text(summary, encoding="utf-8")
    print(f"Run {run_id}: {len(cards)} maps, {matches} reference pairs, {status_label}; {output / f'{run_id}.html'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--elapsed-seconds", type=float)
    parser.add_argument("--legacy-layout", action="store_true", help="Explicit historical report reproduction only; new runs use the shared batch report.")
    parser.add_argument("--compare-run-id", type=int, default=145)
    parser.add_argument("--historical-run-id", type=int, default=145,
                        help="Informational historical map panel, independent of the previous-run regression gate.")
    args = parser.parse_args()
    if args.legacy_layout:
        write_report(args.run_id, args.elapsed_seconds, args.compare_run_id, args.historical_run_id)
    else:
        from climate_run_record import write_run_record
        write_run_record(args.run_id, args.elapsed_seconds)
