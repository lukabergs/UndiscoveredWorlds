"""Earth rock-erodibility reference from GLiM, for comparison with the terrain sandbox's hardness map.

Reads the GLiM v1.0 level-1 lithology grid at 0.5 degree (Hartmann and Moosdorf 2012,
doi:10.1594/PANGAEA.788537) and the erodibility index of Moosdorf, Cohen and von Hagke 2018
(Applied Geography 101, 36-44, doi:10.1016/j.apgeog.2018.10.010). Writes to
refs/processed/lithology/:

  glim_classes_<W>.npy      uint8 GLiM class codes 1-16, 255 where GLiM has no data
  glim_erodibility_<W>.npy  float32 index / land-area-weighted median index (median land = 1);
                            NaN for no data, water bodies, ice and glaciers, and GLiM's own "nd"
  glim_erodibility_<W>.png  check image of ln(relative erodibility) / 0.7, brown hard, blue soft
  glim_summary.json         class table, GEroID cross-tabulation, median, area percentiles
  README.md                 sources, class table, processing, grid convention

Grid: cell-centred W x W/2 equirectangular, rows from the north, columns from 180 W eastward.
The source (0.5 degree) is coarser than the 2048 grid, so the class of each cell is the source
cell that contains its centre (no cell centre falls on a source edge at W = 2048).
"""

import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.hashing import sha256

ROOT = Path(__file__).resolve().parents[2]
GLIM_ZIP = ROOT / "refs/source/glim-2012/hartmann-moosdorf_2012.zip"
GEROID_TIF = ROOT / "refs/source/geroid-2018/Global_BQART_Lithology_GLiMv1.1_06min_Static.tif"
OUTPUT = ROOT / "refs/processed/lithology"
GLIM_MEMBER = "glim_wgs84_0point5deg.txt.asc"
NO_DATA = 255
LOG_SCALE = 0.7
STOPS = ((-1.0, (196, 120, 60)), (0.0, (232, 226, 205)), (1.0, (60, 120, 170)))
NO_DATA_RGB = (18, 22, 30)

PAPER = "Moosdorf et al. 2018"
SECTION = f"{PAPER}, section 3.1 and abstract"
TABLE = f"{PAPER}, Table 1, GLOBAL AVERAGE (=EI) row"
COMBINED_PA = f"{PAPER}, section 2: acid and intermediate plutonics combined (Table 1 column pa+)"
COMBINED_VB = f"{PAPER}, section 2: basic and intermediate volcanics combined (Table 1 column vb+)"
GEROID = "not in the paper; mode of the authors' GEroID grid (Mendeley Data, doi:10.17632/y8g64m2v52.1) over this class"
EXCLUDED = "no erodibility (water, ice or GLiM no data)"

# code: (abbreviation, name, erodibility index or None, source of the index)
CLASSES = {
    1: ("su", "Unconsolidated sediments", 3.2, f"{SECTION}; {TABLE}"),
    2: ("vb", "Basic volcanic rocks", 1.4, f"{SECTION}; {TABLE} (vb+)"),
    3: ("ss", "Siliciclastic sedimentary rocks", 1.5,
        f"{SECTION} ('siliciclastic rocks of all grain sizes'); Table 1 gives 1.4 for ss (mixed or unknown grain size) "
        "and 1.5 for sssh and ssss, and GEroID mostly carries 1.4 over this class"),
    4: ("pb", "Basic plutonic rocks", 1.5, f"{SECTION}; {TABLE}"),
    5: ("sm", "Mixed sedimentary rocks", 1.2, f"{SECTION}; {TABLE}"),
    6: ("sc", "Carbonate sedimentary rocks", 1.0, f"{SECTION}; {TABLE}"),
    7: ("va", "Acid volcanic rocks", 1.1, f"{SECTION}; {TABLE}"),
    8: ("mt", "Metamorphics", 1.0, f"{SECTION}; {TABLE}"),
    9: ("pa", "Acid plutonic rocks", 1.0, f"{SECTION} (the normalisation class); {TABLE} (pa+)"),
    10: ("vi", "Intermediate volcanic rocks", 1.4, f"{COMBINED_VB}; {TABLE}"),
    11: ("wb", "Water bodies", None, EXCLUDED),
    12: ("py", "Pyroclastics", 1.4, GEROID),
    13: ("pi", "Intermediate plutonic rocks", 1.0, f"{COMBINED_PA}; {TABLE}"),
    14: ("ev", "Evaporites", 1.0, GEROID),
    15: ("nd", "No data", None, EXCLUDED),
    16: ("ig", "Ice and glaciers", None, EXCLUDED),
}

RECEIPTS = {
    GLIM_ZIP: "43b4ce3276b155d804db8ff9fb227d620b4c35015a4cf564eac4d06d2b69d88e",
    GEROID_TIF: "0e409866294d4d66ca7cba45f8532b01681fbd09c59927f7423ccda5b59ed2a3",
}


def check_sources():
    for path, digest in RECEIPTS.items():
        if sha256(path) != digest:
            raise SystemExit(f"{path}: sha256 differs from its receipt")


def read_glim(path):
    """Return the 0.5 degree class grid, rows from the north, columns from 180 W, -1 for no data."""
    with zipfile.ZipFile(path) as archive:
        lines = archive.read(GLIM_MEMBER).decode("ascii").splitlines()
    header = {line.split()[0].lower(): float(line.split()[1]) for line in lines[:6]}
    expected = {"ncols": 720, "nrows": 360, "xllcorner": -180, "yllcorner": -90, "cellsize": 0.5}
    if any(header[key] != value for key, value in expected.items()):
        raise SystemExit(f"unexpected GLiM grid header {header}")
    rows = [np.array(line.split(), dtype=np.int32) for line in lines[6:] if line.strip()]
    grid = np.vstack(rows)
    grid[grid == int(header["nodata_value"])] = -1
    if grid.shape != (360, 720) or not set(np.unique(grid)) <= {-1, *CLASSES}:
        raise SystemExit("unexpected GLiM grid contents")
    return grid


def cell_centres(width):
    height = width // 2
    lon = -180 + (np.arange(width) + 0.5) * 360 / width
    lat = 90 - (np.arange(height) + 0.5) * 180 / height
    return lon, lat


def regrid(source, width):
    """Class of the 0.5 degree source cell containing each target cell centre."""
    lon, lat = cell_centres(width)
    columns = np.floor((lon + 180) / 0.5).astype(int)
    rows = np.floor((90 - lat) / 0.5).astype(int)
    grid = source[np.ix_(rows, columns)]
    return np.where(grid < 0, NO_DATA, grid).astype(np.uint8)


def cell_areas(width):
    """Relative areas of the cells of one column (rows from the north)."""
    height = width // 2
    edges = np.radians(90 - np.arange(height + 1) * 180 / height)
    return np.sin(edges[:-1]) - np.sin(edges[1:])


def weighted_quantile(values, weights, q):
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cumulative, q * cumulative[-1])])


def geroid_crosstab(glim):
    """Share of GEroID values over each GLiM class, sampled at the GLiM 0.5 degree cell centres."""
    import rasterio

    with rasterio.open(GEROID_TIF) as source:
        data = source.read(1)
        nodata = source.nodata
        lon, lat = cell_centres(720)
        lons, lats = np.meshgrid(lon, lat)
        rows, columns = rasterio.transform.rowcol(source.transform, lons.ravel(), lats.ravel())
    rows = np.clip(np.asarray(rows), 0, data.shape[0] - 1)
    columns = np.clip(np.asarray(columns), 0, data.shape[1] - 1)
    sampled = np.round(data[rows, columns].reshape(glim.shape), 2)
    table = {}
    for code in CLASSES:
        values = sampled[(glim == code) & (sampled != nodata)]
        unique, counts = np.unique(values, return_counts=True)
        order = np.argsort(-counts)
        table[code] = [(float(unique[i]), round(float(counts[i] / counts.sum()), 3)) for i in order[:3]]
    return table


def render(relative, path):
    t = np.clip(np.log(relative) / LOG_SCALE, -1, 1)
    rgb = np.zeros(relative.shape + (3,), dtype=np.float64)
    (t0, c0), (t1, c1), (t2, c2) = STOPS
    low = t < 0
    for channel in range(3):
        rgb[..., channel] = np.where(
            low,
            c0[channel] + (t - t0) / (t1 - t0) * (c1[channel] - c0[channel]),
            c1[channel] + (t - t1) / (t2 - t1) * (c2[channel] - c1[channel]))
    rgb[~np.isfinite(relative)] = NO_DATA_RGB
    Image.fromarray(np.round(rgb).astype(np.uint8), "RGB").save(path)


def write_readme(path, width, summary):
    lines = [
        "# Earth rock erodibility from GLiM",
        "",
        "Written by `scripts/refs/glim.py`; reference for the terrain sandbox's procedural hardness map.",
        "",
        "## Sources",
        "",
        "- GLiM v1.0 gridded to 0.5 degree: Hartmann, J. and Moosdorf, N. (2012), PANGAEA, "
        "doi:10.1594/PANGAEA.788537 (CC-BY-3.0), supplement to Hartmann and Moosdorf (2012), "
        "G3 13, Q12004, doi:10.1029/2012GC004370. File `refs/source/glim-2012/hartmann-moosdorf_2012.zip` "
        "(38,670 bytes; receipt beside it). Level-1 classes only (16 classes); the multi-GB polygon "
        "database was not downloaded.",
        "- Erodibility index: Moosdorf, N., Cohen, S. and von Hagke, C. (2018), A global erodibility "
        "index to represent sediment production potential of different rock types, Applied Geography 101, "
        "36-44, doi:10.1016/j.apgeog.2018.10.010. Full text `refs/source/geroid-2018/Moosdorf_et_al_AppGeog2018.pdf`.",
        "- GEroID, the authors' 0.1 degree grid of the index (Mendeley Data, doi:10.17632/y8g64m2v52.1), "
        "`refs/source/geroid-2018/Global_BQART_Lithology_GLiMv1.1_06min_Static.tif`; used only for the two "
        "classes the paper does not list (pyroclastics, evaporites).",
        "",
        "## Classes",
        "",
        "The index is the mean hillslope of acid plutonic rocks divided by that of the class in regions of "
        "equal rock uplift (paper, eq. 1): 1.0 is as hard as granite, higher erodes more easily. "
        f"Relative erodibility = index / {summary['median_index']} (the land-area-weighted median index).",
        "",
        "| Code | GLiM | Class | Index | Relative | Land area share | GEroID over the class (value: share) | Source of the index |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for code, (abbreviation, name, index, source) in CLASSES.items():
        row = summary["classes"][str(code)]
        relative = "" if index is None else f"{index / summary['median_index']:.3f}"
        geroid = ", ".join(f"{value}: {share}" for value, share in row["geroid"])
        lines.append(f"| {code} | {abbreviation} | {name} | {'' if index is None else index} | {relative} | "
                     f"{row['area_share']:.4f} | {geroid} | {source} |")
    percentiles = summary["relative_percentiles"]
    lines += [
        "",
        "Land area share: share of the area with a GLiM class (codes 1-16) on this grid. GEroID column: the "
        "authors' grid (built from GLiM v1.1 at 0.1 degree) sampled at the centres of the 0.5 degree GLiM v1.0 "
        "cells of each class, leading three values; the two maps differ at unit edges, so shares below 1 are "
        "expected. Every class with a paper value carries it in GEroID in over 93 % of its cells except "
        "siliciclastic rocks (ss), where GEroID has 1.4 in about 70 % and 1.5 in about 26 % (probably GLiM's "
        "level-2 grain-size subclasses: Table 1 gives 1.4 for mixed or unknown grain size, 1.5 for fine and "
        "coarse). GEroID gives ice and glaciers 1.0 and water bodies 0; here both are NaN.",
        "",
        "## Processing",
        "",
        "1. Read `glim_wgs84_0point5deg.txt.asc` (ESRI ASCII grid, 720 x 360, lower-left corner 180 W 90 S, "
        "rows from the north, -9999 no data).",
        f"2. Each cell of the {width} x {width // 2} grid takes the class of the 0.5 degree cell containing its "
        "centre (the source is coarser, so this is both the nearest and the majority class; no centre lies on "
        "a source edge). GLiM no data becomes 255; water bodies (11), GLiM no data (15) and ice and "
        "glaciers (16) keep their codes.",
        "3. Index per class from the table above; NaN for 255, 11, 15 and 16.",
        f"4. Median of the index weighted by cell area over the cells with an index: {summary['median_index']}. "
        "Relative erodibility = index / median, float32.",
        f"5. Check image: ln(relative) / {LOG_SCALE} clipped to [-1, 1], stops -1 (196,120,60) brown hard, "
        "0 (232,226,205), +1 (60,120,170) blue soft; NaN (18,22,30).",
        "",
        "Area-weighted percentiles of relative erodibility over land with an index: "
        f"10%: {percentiles['10']}, 50%: {percentiles['50']}, 90%: {percentiles['90']}.",
        "",
        "## Grid",
        "",
        f"Cell-centred equirectangular {width} x {width // 2}: column i spans longitudes "
        f"-180 + 360 i / {width} to -180 + 360 (i + 1) / {width}, row j spans latitudes "
        f"90 - 180 j / {width // 2} to 90 - 180 (j + 1) / {width // 2} (row 0 at the North Pole). "
        "Arrays are numpy `.npy`, shape (rows, columns).",
        "",
        "## Limits",
        "",
        "- One value per lithological class, calibrated on hillslopes in five landscapes; no weathering, "
        "fracturing, climate or soil cover. Level-1 classes only, so siliciclastic rocks of every grain size "
        "share one value.",
        "- 0.5 degree source: small units (dykes, narrow belts) are lost; coastlines are blocky.",
        "- The pyroclastic and evaporite values are the authors' GEroID assignments read from their grid, "
        "not values stated in the paper.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--width", type=int, default=2048)
    arguments = parser.parse_args()
    width = arguments.width

    check_sources()
    glim = read_glim(GLIM_ZIP)
    classes = regrid(glim, width)

    index_of = np.full(256, np.nan)
    for code, (_, _, index, _) in CLASSES.items():
        if index is not None:
            index_of[code] = index
    index = index_of[classes]
    areas = np.broadcast_to(cell_areas(width)[:, None], classes.shape)
    land = np.isfinite(index)
    median = weighted_quantile(index[land], areas[land], 0.5)
    relative = (index / median).astype(np.float32)

    classified = classes != NO_DATA
    total = areas[classified].sum()
    crosstab = geroid_crosstab(glim)
    summary = {
        "grid": {"width": width, "height": width // 2, "rows": "north first", "columns": "from 180 W eastward"},
        "median_index": median,
        "relative_percentiles": {str(q): round(weighted_quantile(relative[land].astype(np.float64), areas[land], q / 100), 4)
                                 for q in (10, 50, 90)},
        "land_with_index_share_of_classified": round(float(areas[land].sum() / total), 4),
        "classes": {str(code): {"abbreviation": abbreviation, "name": name, "index": value, "source": source,
                                "area_share": round(float(areas[classes == code].sum() / total), 5),
                                "geroid": crosstab[code]}
                    for code, (abbreviation, name, value, source) in CLASSES.items()},
        "sources": {str(path.relative_to(ROOT)).replace("\\", "/"): digest for path, digest in RECEIPTS.items()},
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    np.save(OUTPUT / f"glim_classes_{width}.npy", classes)
    np.save(OUTPUT / f"glim_erodibility_{width}.npy", relative)
    render(relative, OUTPUT / f"glim_erodibility_{width}.png")
    (OUTPUT / "glim_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_readme(OUTPUT / "README.md", width, summary)
    print(json.dumps({key: summary[key] for key in ("median_index", "relative_percentiles",
                                                     "land_with_index_share_of_classified")}, indent=2))
    for code, row in summary["classes"].items():
        print(code, row["abbreviation"], row["index"], row["area_share"], row["geroid"])


if __name__ == "__main__":
    main()
