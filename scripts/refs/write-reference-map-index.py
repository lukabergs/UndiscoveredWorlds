# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Build a searchable, offline map index from the export manifest and palettes."""

import argparse
import html
import importlib.util
import json
import os
from pathlib import Path

import numpy as np

from reference_map_rendering import scalar_rgb
from reference_map_layout import product_directory
from climate_palettes import anchors, VERSION


def module_at(path):
    spec = importlib.util.spec_from_file_location(path.stem.replace("-", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def colours(levels, limit, palette, legacy=False):
    if legacy:
        stops=((0,(5,12,24)),(.55*limit,(18,125,175)),(limit,(250,235,125))) if palette=='speed' else (
            ((-limit,(25,175,245)),(0,(18,20,25)),(limit,(245,120,35))) if palette=='divergence' else
            ((-limit,(225,75,45)),(0,(18,20,25)),(limit,(45,220,185))))
        rgb=np.rint(np.stack([np.interp(levels,[v for v,_ in stops],[c[i] for _,c in stops]) for i in range(3)],axis=-1)).astype(int)
    else:
        stops=anchors(0 if palette=='speed' else -limit,limit,palette!='speed')
        return [(v,'#'+''.join(f'{c:02x}' for c in rgb)) for v,rgb in stops]
    return [(level, "#" + "".join(f"{int(c):02x}" for c in colour))
            for level, colour in zip(levels, rgb)]


def describe(product, additional, generator):
    name = product["name"]
    key = name[4:] if name[:3] in ("jan", "apr", "jul", "oct") else name
    notes = product.get("numeric_meaning", product.get("comparison", ""))
    legacy=product.get('palette_version')!=VERSION
    if key in ("temperature", "precipitation"):
        old={'temperature':[(-54,'#e1e1e1'),(-40,'#b900ff'),(-25,'#2d00ff'),(-10,'#0069ff'),(0,'#00dce6'),(10,'#00cd46'),(20,'#ffe600'),(30,'#dc2d00')],
             'precipitation':[(0,'#040a23'),(250,'#0c235f'),(500,'#0f50a0'),(1000,'#0a96cd'),(1500,'#0faf5a'),(2500,'#87cd23'),(3500,'#f5e119'),(4500,'#fa7314'),(6000,'#dc1914')]}
        stops = old[key] if legacy else [(v, "#" + "".join(f"{c:02x}" for c in rgb)) for v, rgb in generator.visual_palette(key)]
        description = {"temperature": "Four representative months of WorldClim near-surface air temperature, averaged over land",
                       "precipitation": "IMERG annual precipitation, land and ocean"}[key]
        return description, stops, product["units"], notes
    if key == "koppen":
        return "Observed monthly-temperature/rainfall derived Koeppen-Geiger climate class", [], "category", notes
    if key in additional:
        record = additional[key]
        limit, palette = record["limit"], record["palette"]
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
        if not legacy:limit={'surface':25,'850hpa':30,'upper':40}[layer]
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
    if not legacy:
        limit={'column_water':65,'column_moisture_flux':800,'column_moisture_flux_east':800,'column_moisture_flux_north':800,'ocean_current_speed':1}.get(key,limit)
    levels = (0, limit * .55, limit) if palette == "speed" else (-limit, 0, limit)
    unit = "m/s background" if key.endswith(("_lic", "_particles")) else product["units"]
    return description, colours(levels, limit, palette,legacy), unit, notes


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
    generator = module_at(Path(__file__).with_name("prepare-reduced-earth-benchmark.py"))
    classifier = module_at(Path(__file__).parents[1] / "benchmarks/benchmark-observed-koppen.py")
    additional = {}
    for receipt in (root / "climate").glob("*-additional-preparation.json"):
        for record in json.loads(receipt.read_text(encoding="utf-8"))["products"]:
            additional[record["bundle"].removesuffix("_monthly.uwclim")] = record
    by_name = {}
    for resolution in manifest["outputs"]:
        for product in resolution["products"]:
            if not product["maps"][0]["file"].startswith("climate/"):
                continue
            by_name.setdefault(product["name"], []).append((resolution["width"], product))
    rows = []
    for name, variants in sorted(by_name.items(), key=lambda item: str(product_directory(item[0]))):
        product = variants[0][1]
        description, stops, unit, notes = describe(product, additional, generator)
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
<p>Generated by scripts/refs/write-reference-map-index.py from the retained export manifest and renderer palettes. No external scripts, fonts or network requests.</p>
<script>function filter(){let n=0;const q=document.querySelector('#search').value.toLowerCase(),s=document.querySelector('#season').value,w=document.querySelector('#width').value;
document.querySelectorAll('tbody tr').forEach(r=>{r.hidden=!r.textContent.toLowerCase().includes(q)||(s&&r.dataset.season!==s);if(!r.hidden)n++});
document.querySelectorAll('.resolution').forEach(r=>r.hidden=w&&r.dataset.width!==w);document.querySelector('#count').textContent=n+' map types';}
document.querySelectorAll('input,select').forEach(e=>e.addEventListener('input',filter));filter();</script></html>'''
    widths = sorted(r["width"] for r in manifest["outputs"])
    document = document.replace("WIDTHS", "".join(f'<option value="{w}" {"selected" if w == 512 else ""}>{w}</option>' for w in widths))
    document = document.replace("ROWS", "\n".join(rows)).replace("CATEGORIES", "".join(categories)).replace("MISSING", missing)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "index.html").write_text(document, encoding="utf-8")
    print(f"Indexed {len(rows)} climate map types at {len(widths)} widths; all file links exist.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("refs/processed"))
    write_index(parser.parse_args().root)
