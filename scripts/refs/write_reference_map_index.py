"""Build a searchable, offline map index from the export manifest and palettes."""

import argparse
import html
import importlib.util
import json
import os
import re
from pathlib import Path

import numpy as np

from reference_map_rendering import display_stops
from reference_map_layout import product_directory
from reference_climate_maps import additional_display, product_display
from climate_palettes import VERSION


def module_at(path):
    spec = importlib.util.spec_from_file_location(path.stem.replace("-", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


UNKNOWN = "Unknown palette: this product does not record the palette its PNG was drawn with; the legend is the renderer's palette and may not match the image."


def hexadecimal(stops):
    return [(value, "#" + "".join(f"{int(c):02x}" for c in colour)) for value, colour in stops]


def legacy_colours(levels, limit, palette):
    """Anchors of the renderer palettes before climate-colour-v2, for products still rendered with them."""
    stops=((0,(5,12,24)),(.55*limit,(18,125,175)),(limit,(250,235,125))) if palette=='speed' else (
        ((-limit,(25,175,245)),(0,(18,20,25)),(limit,(245,120,35))) if palette=='divergence' else
        ((-limit,(225,75,45)),(0,(18,20,25)),(limit,(45,220,185))))
    rgb=np.rint(np.stack([np.interp(levels,[v for v,_ in stops],[c[i] for _,c in stops]) for i in range(3)],axis=-1)).astype(int)
    return hexadecimal(zip(levels, rgb))


def describe(product, additional, generator):
    name = product["name"]
    key = name[4:] if name[:3] in ("jan", "apr", "jul", "oct") else name
    notes = product.get("numeric_meaning", product.get("comparison", ""))
    legacy=product.get('palette_version')!=VERSION
    # The display a product's PNG was drawn with, recorded when it was written; products
    # written before that fall back to the renderer palettes, flagged as unknown.
    shown = product.get("display")
    if shown is None and key != "koppen":
        notes = (notes + " " if notes else "") + UNKNOWN
    if key in ("temperature", "precipitation"):
        old={'temperature':[(-54,'#e1e1e1'),(-40,'#b900ff'),(-25,'#2d00ff'),(-10,'#0069ff'),(0,'#00dce6'),(10,'#00cd46'),(20,'#ffe600'),(30,'#dc2d00')],
             'precipitation':[(0,'#040a23'),(250,'#0c235f'),(500,'#0f50a0'),(1000,'#0a96cd'),(1500,'#0faf5a'),(2500,'#87cd23'),(3500,'#f5e119'),(4500,'#fa7314'),(6000,'#dc1914')]}
        stops = (hexadecimal(display_stops(**shown)) if shown else
                 old[key] if legacy else hexadecimal(generator.visual_palette(key)))
        description = {"temperature": "Four representative months of WorldClim near-surface air temperature, averaged over land",
                       "precipitation": "IMERG annual precipitation, land and ocean"}[key]
        return description, stops, product["units"], notes
    if key == "koppen":
        return "Observed monthly-temperature/rainfall derived Koeppen-Geiger climate class", [], "category", notes
    if key in additional:
        record = additional[key]
        dataset = key[:-len(record["variable"]) - 1]
        limit, palette = record["limit"], record["palette"]
        display = additional_display(dataset, record)
        description = record.get("description") or key.replace("_", " ")
        if key.startswith("ceres_"):
            words = {"toa": "top of atmosphere", "sfc": "surface", "sw": "shortwave",
                     "lw": "longwave", "all": "all sky", "cre": "cloud radiative effect",
                     "tot": "total", "temp": "temperature", "tau": "optical depth"}
            description = "CERES " + " ".join(words.get(w, w) for w in key.split("_")[1:])
        if key == "oisst_sst":
            description = "NOAA OISST sea-surface temperature (primary SST reference)"
        if key == "era5_vimd":
            description = "ERA5 archived integrated moisture divergence, sign-reversed to convergence"
        if key == "era5_mer":
            description = "ERA5 evaporation, converted to upward-positive water flux"
        if key.startswith("era5land_swvl"):
            depth = {"1": "0-7", "2": "7-28", "3": "28-100", "4": "100-289"}[key[-1]]
            description += f" ({depth} cm below surface)"
    elif "_wind_" in key:
        layer, style = key.split("_wind_", 1)
        description = f"{layer} wind {style}; " + {"surface": "10 m", "upper": "500 hPa", "850hpa": "850 hPa"}[layer]
        limit, palette = 25, "divergence" if style in ("east", "north") else "speed"
        if style in ("lic", "particles"):
            notes += " Background colour is speed tinted by land/ocean and texture; arrows show direction. The PNG has no exact scalar colour lookup."
    elif key.endswith("_divergence"):
        description, limit, palette = key.replace("_", " "), 1, "divergence"
    else:
        description, limit, palette = {
            "column_water": ("Total atmospheric column water vapour", 60, "speed"),
            "pressure_anomaly": ("Sea-level pressure anomaly", 30, "divergence"),
            "ascent": ("500 hPa vertical motion, positive upward", 100, "convergence"),
            "era5_precipitation": ("ERA5 daily precipitation rate", 20, "speed"),
            "column_moisture_flux": ("Magnitude of vertically integrated vapour transport; arrows show direction", 500, "speed"),
            "column_moisture_flux_east": ("Eastward vertically integrated vapour transport", 500, "divergence"),
            "column_moisture_flux_north": ("Northward vertically integrated vapour transport", 500, "divergence"),
            "moisture_flux_convergence": ("Vapour convergence computed on this output grid, positive moisture supply", 20, "convergence"),
            "ocean_current_speed": ("Speed of the mean GLORYS near-surface current vector; arrows show direction", 2, "speed"),
        }[key]
    unit = "m/s background" if key.endswith(("_lic", "_particles")) else product["units"]
    if shown:
        return description, hexadecimal(display_stops(**shown)), unit, notes
    if not legacy or ":" in palette:
        # Current products: the anchors of the palette the renderer used (uwtools.palettes).
        display = display if key in additional else product_display(key)
        return description, hexadecimal(display_stops(*display)), unit, notes
    levels = (0, limit * .55, limit) if palette == "speed" else (-limit, 0, limit)
    return description, legacy_colours(levels, limit, palette), unit, notes


def legend(stops, unit):
    if not stops:
        return '<a href="#koppen">Category legend below</a>'
    lo, hi = stops[0][0], stops[-1][0]
    gradient = ",".join(f"{c} {(v-lo)/(hi-lo)*100:.5f}%" for v, c in stops)
    labels = " · ".join(f"{v:g}: {c}" for v, c in stops)
    return (f'<div class="bar" style="background:linear-gradient(to right,{gradient})"></div>'
            f'<small>{html.escape(unit)}<br>{labels}</small>')


def write_index(root):
    directory = root / "climate/maps"
    manifest = json.loads((root / "metadata/reference-maps/manifest.json").read_text(encoding="utf-8"))
    generator = module_at(Path(__file__).with_name("prepare_reduced_earth_benchmark.py"))
    classifier = module_at(Path(__file__).parents[1] / "benchmarks/benchmark_observed_koppen.py")
    additional = {}
    for receipt in root.glob("*/*/*-additional-preparation.json"):
        for record in json.loads(receipt.read_text(encoding="utf-8"))["products"]:
            additional[f"{receipt.parent.parent.name}_{record['variable']}"] = record
    by_name = {}
    for resolution in manifest["outputs"]:
        for product in resolution["products"]:
            if not product["maps"][0]["file"].startswith("climate/"):
                continue
            by_name.setdefault(product["name"], []).append((resolution["width"], product))
    rows = []
    catalogue = []
    for name, variants in sorted(by_name.items(), key=lambda item: str(product_directory(item[0]))):
        product = variants[0][1]
        description, stops, unit, notes = describe(product, additional, generator)
        catalogue.append(dict(name=name, product=product, description=description, unit=unit,
                              folder=product_directory(name).relative_to("climate").as_posix(),
                              widths=sorted(width for width, _ in variants)))
        links = []
        for width, variant in sorted(variants):
            files = (*variant["maps"], variant["csv"])
            anchors = []
            for record in files:
                path = root / record["file"]
                if not path.is_file():
                    raise ValueError(f"Index target missing: {path}")
                relative = Path(os.path.relpath(path, directory)).as_posix()
                anchors.append(f'<a href="{html.escape(relative)}">{path.suffix[1:].upper()}</a>')
            links.append(f'<span class="resolution" data-width="{width}">{width}: {" / ".join(anchors)}</span>')
        period = product.get("period", "annual / representative-month average" if name != "koppen" else "12-month climatology")
        source_names = sorted({Path(p).name for p in product.get("sources", [])})
        source_names = source_names or [Path(manifest["temperature_source" if name == "temperature" else "precipitation_source"]).name]
        sources = "<br>".join(html.escape(s) for s in source_names)
        if len(source_names) > 3:
            sources = (f'<details><summary>{len(source_names)} source files: '
                       f'{html.escape(source_names[0])} …</summary>{sources}</details>')
        season = {"jan": "1", "apr": "2", "jul": "3", "oct": "4"}.get(name[:3], "annual")
        folder = product_directory(name).relative_to("climate").as_posix()
        rows.append(f'<tr data-season="{season}"><td><b>{html.escape(folder)}</b><p>{html.escape(description)}</p>'
                    f'<small>{html.escape(period)}<br>{sources}</small><p>{html.escape(notes)}</p></td>'
                    f'<td>{legend(stops, unit)}</td><td>{"".join(links)}<small>Numeric units: {html.escape(product["units"])}</small></td></tr>')
    categories = ['<span class="category"><i style="background:#000000"></i>0 ocean</span>']
    for index, (code, rgb) in enumerate(zip(classifier.CLIMATE_CODES, classifier.CLIMATE_COLOURS), 1):
        colour = "#" + "".join(f"{int(c):02x}" for c in rgb)
        categories.append(f'<span class="category"><i style="background:{colour}"></i>{index} {code} {colour}</span>')
    missing = "".join(f'<li>{html.escape(p["product"])}: {html.escape(p["reason"])}</li>'
                      for p in manifest["outputs"][0]["unavailable"])
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>Climate reference map index</title>
<style>body{font:16px/1.5 system-ui;max-width:1500px;margin:2rem auto;padding:0 24px;color:#202632;background:#fafbfc}
h1{margin-bottom:0}p{max-width:1050px}a{color:#176298}table{border-collapse:collapse;width:100%;background:white}
td,th{border-bottom:1px solid #d9dfe8;padding:16px;text-align:left;vertical-align:top}td:first-child{width:48%}
td:nth-child(2){width:30%;min-width:230px}td p{margin:8px 0}small{color:#526075;overflow-wrap:anywhere}
.bar{height:22px;border:1px solid #909aaa;border-radius:3px;margin:6px 0}.resolution{display:block;white-space:nowrap}
.category{display:inline-block;margin:6px 16px 6px 0}.category i{display:inline-block;width:20px;height:20px;border:1px solid #777;vertical-align:middle;margin-right:6px}
.filters{position:sticky;top:0;background:#edf2f8;padding:14px;margin:24px 0}input,select{padding:8px;font:inherit;margin-right:14px}input{width:35%}
[hidden]{display:none!important}</style>
<h1>Climate reference maps</h1>
<p><a href="standard-v2/index.html">Current standard palette set, standalone ocean LIC/particles and common legends</a>. This catalogue also retains historical renderings; use the standard set for current batch comparisons.</p>
<p>Search the quantity or provider, select a month and width, then open a PNG preview or an unclamped float32 GeoTIFF / CSV.
All maps are cell-centred W × W/2, north up, longitude −180° to 180°. <a href="README.md">Path and comparison guide</a>.</p>
<p><b>1 = January, 2 = April, 3 = July, 4 = October.</b> Wind LIC, particles, speed and components use the named month.
Other seasonal products average the three months beginning there, weighted by days. The period in each row is authoritative.</p>
<p>Palette numbers below are fixed display limits, not observed minima/maxima. Values outside them are clipped in PNGs.
Black is missing or masked; TIFF nodata is −9999.9 and CSV missing values are NaN.
LIC and particle TIFFs contain texture intensity [0,1], not wind speed; use the speed/east/north fields for calculations.</p>
<p>These are observations, reanalyses and derived diagnostics with different averaging periods, depths and masks.
The primary SST is NOAA OISST; ERA5 and GLORYS alternatives have provider subfolders. Cloud radiative effect (CRE) is all-sky minus clear-sky net radiative flux.
Net radiation and ERA5 heat fluxes are downward positive; evaporation is upward positive. SW = shortwave, LW = longwave, TOA = top of atmosphere.</p>
<div class="filters"><label>Find <input id="search" placeholder="wind, soil, CERES, shortwave…"></label><label>Month <select id="season"><option value="">All</option><option value="annual">Annual</option><option value="1">1 · January</option><option value="2">2 · April</option><option value="3">3 · July</option><option value="4">4 · October</option></select></label><label>Width <select id="width"><option value="">All</option>WIDTHS</select></label> <span id="count"></span></div>
<table><thead><tr><th>Map, meaning and source</th><th>Palette anchors: value → RGB hex</th><th>Files</th></tr></thead><tbody>ROWS</tbody></table>
<h2 id="koppen">Koeppen category palette</h2>CATEGORIES
<p>A = tropical; B = arid (W desert, S steppe; h hot, k cold); C = temperate; D = continental; E = polar (T tundra, F ice cap).
For A/C/D, f = no dry season, s = dry summer, w = dry winter, m = monsoon. Temperature suffixes: a hot summer, b warm summer, c cool summer, d extremely cold winter.</p>
<h2>Unavailable comparison products</h2><ul>MISSING</ul>
<p>Generated by scripts/refs/write_reference_map_index.py from the retained export manifest and renderer palettes. No external scripts, fonts or network requests.</p>
<script>function filter(){let n=0;const q=document.querySelector('#search').value.toLowerCase(),s=document.querySelector('#season').value,w=document.querySelector('#width').value;
document.querySelectorAll('tbody tr').forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q)||(s&&r.dataset.season!==s);if(!r.hidden)n++});
document.querySelectorAll('.resolution').forEach(r=>r.hidden=w&&r.dataset.width!==w);document.querySelector('#count').textContent=n+' map types';}
document.querySelectorAll('input,select').forEach(e=>e.addEventListener('input',filter));filter();</script></html>'''
    widths = sorted(r["width"] for r in manifest["outputs"])
    document = document.replace("WIDTHS", "".join(f'<option value="{w}" {"selected" if w == 512 else ""}>{w}</option>' for w in widths))
    document = document.replace("ROWS", "\n".join(rows)).replace("CATEGORIES", "".join(categories)).replace("MISSING", missing)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "index.html").write_text(document, encoding="utf-8")
    write_guide(root, directory, manifest, catalogue, widths)
    print(f"Indexed {len(rows)} climate map types at {len(widths)} widths; all file links exist; README.md written.")


REPOSITORY = Path(__file__).resolve().parents[2]
MONTHS = re.compile(r"\b(jan|apr|jul|oct)\b")

# Preparation commands in build order: (processed dataset, command, default refs/source input).
PREPARATION = [
    ("heightmap, land_ocean", "uv run scripts/refs/prepare_earth_maps.py --width 2048", "refs/source/ETOPO, refs/source/GSHHG"),
    ("worldclim", "uv run --offline python scripts/refs/prepare_climate_reference.py --width 2048", "refs/source/worldclim-2.1-10m"),
    ("imerg", "uv run --offline python scripts/refs/prepare_imerg_reference.py --width 2048", "refs/source/imerg-v07b-annual"),
    ("imerg", "uv run --offline python scripts/refs/prepare_imerg_reference.py --width 2048 --monthly-input refs/source/IMERGmonthlyClimatology2001to2022GeoTIFF", "refs/source/IMERGmonthlyClimatology2001to2022GeoTIFF"),
    ("era5", "uv run --offline python scripts/refs/prepare_diagnostic_references.py --width 2048", "refs/source/era5-planette-2001-2020"),
    ("era5", "uv run scripts/refs/prepare_additional_references.py --family era5 --input refs/source/era5-cds-2001-2020 --width 2048", "refs/source/era5-cds-2001-2020"),
    ("era5land", "uv run scripts/refs/prepare_additional_references.py --family land --input refs/source/era5-land-2001-2020 --width 2048", "refs/source/era5-land-2001-2020"),
    ("oisst", "uv run scripts/refs/prepare_additional_references.py --family oisst --input refs/source/noaa-oisst-v2.1-2001-2020 --width 2048", "refs/source/noaa-oisst-v2.1-2001-2020"),
    ("glorys", "uv run scripts/refs/prepare_additional_references.py --family marine --input refs/source/glorys12-2001-2020 --width 2048", "refs/source/glorys12-2001-2020"),
    ("ceres", "uv run scripts/refs/prepare_additional_references.py --family ceres --input refs/source/ceres-ebaf-edition4.2.1 --width 2048", "refs/source/ceres-ebaf-edition4.2.1"),
]
RENDERING = [
    "uv run --offline python scripts/refs/prepare_reduced_earth_benchmark.py --width 128 256 512 1024 2048",
    "uv run --offline scripts/refs/render_standard_climate_maps.py",
    "uv run --offline scripts/refs/render_standard_flow_references.py",
    "uv run --offline scripts/refs/write_reference_map_index.py",
]


def shown(path):
    """A path as written from the repository root (junctions are not followed)."""
    return Path(os.path.relpath(os.path.abspath(path), REPOSITORY)).as_posix()


def cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def source_datasets(root):
    """{processed dataset: refs/source folders or URLs its preparation receipts name}."""
    found = {}
    receipts = [(r.parent.parent.name, json.loads(r.read_text(encoding="utf-8")))
                for r in sorted(root.glob("*/*/*-preparation.json"))]
    diagnostic = root / "climate/diagnostic-reference-preparation.json"
    if diagnostic.is_file():
        receipts.append(("era5", {"era5": json.loads(diagnostic.read_text(encoding="utf-8")).get("era5")}))
    for dataset, receipt in receipts:
        text = json.dumps(receipt).replace("\\\\", "/")
        folders = set(re.findall(r"refs/source/([^/\"]+)", text))
        if not folders and isinstance(receipt.get("source"), str):
            folders = {receipt["source"]}
        found.setdefault(dataset, set()).update(folders)
    return found


def write_guide(root, directory, manifest, catalogue, widths):
    """refs/processed/climate/maps/README.md, the path and comparison guide linked from index.html."""
    origins = source_datasets(root)
    folder_dataset = {folder: dataset for dataset, folders in origins.items() for folder in folders}
    palettes_file = REPOSITORY / "configs/climate-palettes.json"
    definition = json.loads(palettes_file.read_text(encoding="utf-8"))
    profile_file = REPOSITORY / "configs/evaluation/climate.json"
    profile = json.loads(profile_file.read_text(encoding="utf-8"))
    lines = ["# Path and comparison guide", "",
             f"Generated by `scripts/refs/write_reference_map_index.py` with [index.html](index.html) from "
             f"`{shown(root / 'metadata/reference-maps/manifest.json')}`, the preparation receipts, "
             f"`configs/climate-palettes.json` and `configs/evaluation/climate.json`. Do not edit: rerun the script.", ""]

    # Families: one row per map type, the season folder written <s>.
    groups = {}
    for entry in catalogue:
        key = entry["name"][4:] if entry["name"][:3] in ("jan", "apr", "jul", "oct") else entry["name"]
        groups.setdefault(key, []).append(entry)
    families = {}
    for key, entries in groups.items():
        parts = [entry["folder"].split("/") for entry in entries]
        pattern = "/".join(p if all(len(q) == len(parts[0]) and q[i] == p for q in parts) else "<s>"
                           for i, p in enumerate(parts[0]))
        first = entries[0]
        product = first["product"]
        datasets = set()
        for source in product.get("sources", []):
            source = source.replace("\\", "/")
            if source.startswith("refs/source/"):
                folder = source.split("/")[2]
                datasets.add(folder_dataset.get(folder, folder))
            elif source.startswith("refs/processed/"):
                datasets.add(source.split("/")[2])
        if not datasets:
            datasets = {Path(manifest["temperature_source" if key == "temperature" else "precipitation_source"]).parts[2]}
        display = product.get("display")
        if display is None:
            palette = "Koeppen categories" if key == "koppen" else "not recorded"
        elif ":" in display["palette"]:
            palette = display["palette"]
        else:
            palette = f"{display['palette']} (limit {display['limit']}, {product.get('palette_version', 'unversioned')})"
        seasons = "1-4" if len(entries) > 1 else "annual"
        period = MONTHS.sub("the named month", product.get("period", "annual / representative-month average" if key != "koppen" else "12-month climatology"))
        families.setdefault(pattern.split("/")[0], []).append(
            (pattern, first["description"], first["unit"], period, ", ".join(sorted(datasets)), palette, seasons))
    lines += ["## What this folder holds", "",
              f"{len(catalogue)} map products ({len(groups)} map types) in {len(families)} families at widths "
              f"{', '.join(map(str, widths))}, plus the standard 2048-wide set in `standard-v2/`. "
              "Seasons are folder `1` January, `2` April, `3` July, `4` October; winds, pressure and LIC use the named month, "
              "other seasonal products the day-weighted quarter starting there. The searchable index adds palette anchors and links.", ""]
    for family, rows in sorted(families.items()):
        lines += [f"### {family} ({len(rows)})", "",
                  "| Folder | Meaning | Units | Period | Dataset | Palette | Seasons |", "| --- | --- | --- | --- | --- | --- | --- |"]
        lines += ["| " + " | ".join(cell(v) for v in (f"`{row[0]}`",) + row[1:]) + " |" for row in sorted(rows)]
        lines.append("")
    standard = directory / "standard-v2/manifest.json"
    if standard.is_file():
        record = json.loads(standard.read_text(encoding="utf-8"))
        flow = directory / "standard-v2/wind-flow-manifest.json"
        flows = json.loads(flow.read_text(encoding="utf-8")) if flow.is_file() else {}
        lines += ["### standard-v2", "",
                  f"[standard-v2/index.html](standard-v2/index.html): {len(record['reference_maps'])} images "
                  f"({record.get('palette_version')}) plus {len(flows.get('products', {}))} wind LIC/particle images "
                  f"(seed {flows.get('seed')}). {record.get('method', '')}", ""]
        lines += [f"- `{field}`: {state}" for field, state in record.get("coverage", {}).items() if state != "Available"]
        lines.append("")
    unavailable = manifest["outputs"][0]["unavailable"]
    if unavailable:
        lines += ["### Unavailable", ""] + [f"- {p['product']}: {p['reason']}" for p in unavailable] + [""]

    first = catalogue[0]
    width = first["widths"][-1]
    output = next(o for o in manifest["outputs"] if o["width"] == width)
    example = next(p for p in output["products"] if p["name"] == first["name"])
    lines += ["## Path layout", "",
              f"- Grid: {manifest['grid']}; latitude = {manifest['latitude_formula']}; longitude = -180 + (x + 0.5) * 360 / W; "
              "north up, -180 to 180. Map layout version " + str(manifest.get("map_layout_version")) + ".",
              f"- Widths: {', '.join(f'{w} x {w // 2}' for w in widths)}. File names are the width only.",
              "- `maps/climate/<quantity>/<style or season>/<season>/<layer>[/<provider>]/<W>.png` (this folder), with the same "
              "folders under `../fields/` (float32 GeoTIFF) and `../csv/`. Layers: `s` surface, `u` 500 hPa, `850` 850 hPa, "
              "`toa`, `column`; a provider folder (`era5`, `glorys`, `ceres`, ...) marks an alternative to the primary product.",
              f"- Numeric fields: {manifest['transformations'].get('numeric_fields', '')}. CSV missing values are NaN; black PNG pixels are missing.",
              f"- Example ({first['name']}, {width}): " + ", ".join(f"`{shown(root / f['file'])}`" for f in (*example["maps"], example["csv"])),
              f"- Transformations: " + "; ".join(f"{k}: {v}" for k, v in manifest["transformations"].items() if k != "numeric_fields") + ".", ""]

    sidecar = directory / "palette.json"
    recorded = json.loads(sidecar.read_text(encoding="utf-8"))["palettes"] if sidecar.is_file() else {}
    current = sum(1 for e in catalogue if e["product"].get("palette_version") == definition["version"])
    lines += ["## Colour scales", "",
              f"Runs and references use the project palettes (`scripts/uwtools/palettes.py`). `configs/climate-palettes.json` "
              f"is `{definition['version']}`; {current} of {len(catalogue)} products here were drawn with it"
              + (" (the table names the palette of the others)" if current < len(catalogue) else "")
              + "; a palette without `climate:` is a renderer palette with the fixed display limit shown. `palette.json` beside the maps records the palette files in force when "
              "they were rendered: " + ("; ".join(f"{k} {v['version']} ({v['file']}, sha256 {v['sha256'][:12]})" for k, v in recorded.items()) or "none written yet") + ".",
              "Positive scales run blue to red; signed scales run white-magenta-blue-cyan-green-yellow-red with zero cyan; "
              "missing is black; endpoints saturate. Model-minus-reference maps use the `difference` scale.", "",
              "| Scale | Range | Kind | Difference half range |", "| --- | --- | --- | --- |"]
    halves = definition.get("difference", {}).get("half_ranges", {})
    lines += [f"| `{name}` | {low} to {high} | {'signed' if low < 0 else 'positive'} | {halves.get(name, '')} |"
              for name, (low, high) in definition["scales"].items()]
    lines.append("")

    maps_by_validation = {m.get("validation"): m.get("palette", "") for m in profile.get("maps", [])}
    lines += ["## Comparing a run with these references", "",
              "Scores do not read these images: `scripts/uw.py climate` evaluates a run bundle against the prepared fields "
              "`refs/processed/<dataset>/<width>/fields/<variable>_<layout>.npy` (named `<dataset>_<variable>`, read through "
              "`scripts/uwtools/references.py`) on common masks with spherical area weights, scored by global-balanced-v1. "
              "The run page draws model, reference and model minus reference in the palettes above, so a run map is compared "
              "with the PNG of the same quantity and season here.", "", "```powershell",
              "uv run --offline scripts/uw.py climate experiment --resolution 128 --climate-resolution 64",
              "uv run --offline scripts/uw.py climate evaluate --bundle runs/climate/<id>/bundle --baseline runs/climate/<baseline>/bundle --output out/evaluation/climate/<id>",
              "uv run --offline scripts/uw.py climate report --evaluation runs/climate/<baseline>/evaluation runs/climate/<id>/evaluation --output runs/reports/<name>",
              "```", "",
              f"Validation entries of `configs/evaluation/climate.json` ({len(profile.get('validation', []))}):", "",
              "| Entry | Kind | Season | Reference fields | Review palette |", "| --- | --- | --- | --- | --- |"]
    lines += [f"| `{v['id']}` | {v['kind']} | {v.get('season', '')} | {', '.join(v['reference'])} | {maps_by_validation.get(v['id'], '')} |"
              for v in profile.get("validation", [])]
    lines += ["", "Commands and the score are described in `docs/reference-guides/climate-reporting.md` and `benchmark-commands.md`.", ""]

    lines += ["## Rebuilding", "", "Prepared fields first (each writes `refs/processed/<dataset>/<width>/`), then the maps, then this guide:", "",
              "```powershell"] + [command for _, command, _ in PREPARATION] + RENDERING + ["```", "",
              "## Sources", "",
              "| Processed dataset | Preparation input | `refs/source` folders named by its receipts | Map types |", "| --- | --- | --- | --- |"]
    counts = {}
    for rows in families.values():
        for row in rows:
            for dataset in row[4].split(", "):
                counts[dataset] = counts.get(dataset, 0) + 1
    inputs = {}
    for datasets, _, source in PREPARATION:
        for dataset in datasets.split(", "):
            inputs.setdefault(dataset, []).extend(source.split(", "))
    for dataset in sorted(set(inputs) | set(origins) | set(counts)):
        named = ", ".join(f"`{f}`" if not f.startswith("http") else f for f in sorted(origins.get(dataset, ()))) or "not recorded"
        lines.append(f"| `{dataset}` | {', '.join(f'`{s}`' for s in inputs.get(dataset, [])) or '-'} | {named} | {counts.get(dataset, 0)} |")
    lines.append("")
    (directory / "README.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("refs/processed"))
    write_index(parser.parse_args().root)
