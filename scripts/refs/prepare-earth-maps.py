# /// script
# dependencies = ["rasterio==1.5.1", "fiona==1.10.1", "numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Rebuild Earth maps: uv run scripts/refs/prepare-earth-maps.py [--width 2048] [--verify-only]."""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import fiona
import numpy as np
import rasterio
from PIL import Image
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

ROOT = Path(__file__).resolve().parents[2] / "refs"
SOURCES = ROOT / "source"
WIDTH, HEIGHT = 2048, 1024
SW, SH = 21600, 10800
BOUNDS = (-180, -90, 180, 90)
TRANSFORM = from_bounds(*BOUNDS, WIDTH, HEIGHT)
SOURCE_TRANSFORM = from_bounds(*BOUNDS, SW, SH)
ETOPO = SOURCES / "ETOPO/ETOPO_2022_v1_60s_N90W180_surface.tif"
GSHHG = SOURCES / "GSHHG/gshhg-shp-2.3.7.zip"
HEIGHTMAP = ROOT / "processed/heightmap/2048.tif"
MASKS = {k: ROOT / f"processed/{k}/2048.png" for k in ("land_ocean", "lakes", "rivers")}
METADATA = ROOT / "processed/metadata/2048"
URLS = {
    ETOPO.name: "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/60s/60s_surface_elev_gtif/ETOPO_2022_v1_60s_N90W180_surface.tif",
    GSHHG.name: "https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip",
}


def log(message):
    print(message, flush=True)


def configure_grid(width):
    global WIDTH, HEIGHT, TRANSFORM, HEIGHTMAP, MASKS, METADATA
    WIDTH, HEIGHT = width, width // 2
    TRANSFORM = from_bounds(*BOUNDS, WIDTH, HEIGHT)
    HEIGHTMAP = ROOT / f"processed/heightmap/{WIDTH}.tif"
    MASKS = {k: ROOT / f"processed/{k}/{WIDTH}.png" for k in ("land_ocean", "lakes", "rivers")}
    METADATA = ROOT / f"processed/metadata/{WIDTH}"


def describe(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"file": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size,
            "sha256": digest.hexdigest()}


def shp(member):
    return f"zip://{GSHHG.as_posix()}!{member}"


def burn(member, target, transform, value, all_touched=False):
    with fiona.open(shp(member)) as features:
        assert features.crs.to_epsg() == 4326
        shapes = ((geometry(f, member), value(f) if callable(value) else value) for f in features)
        rasterize(shapes, out=target, transform=transform, all_touched=all_touched)
        log(f"Rasterized {len(features):,} features: {member}")


def geometry(feature, member):
    shape = feature["geometry"]
    if member == "GSHHS_shp/f/GSHHS_f_L5.shp" and feature["properties"]["id"] == "4-E":
        ring = list(shape["coordinates"][0])
        # The distributed eastern Antarctic polygon has duplicated (0,-90)
        # vertices and closes east-to-west across the interior. Keep every
        # coastline vertex and close along the antimeridian and South Pole.
        assert ring[1] == ring[2] == (0.0, -90.0) and ring[0] == ring[3] == ring[-1]
        coast = ring[3:-1]
        assert coast[0][0] == 0 and coast[-1][0] == 180
        return {"type": "Polygon", "coordinates": [coast + [(180, -90), (0, -90), coast[0]]]}
    return shape


def downsample(array):
    result = np.zeros((HEIGHT, WIDTH), dtype=np.float32)
    reproject(array, result, src_transform=SOURCE_TRANSFORM, src_crs="EPSG:4326",
              dst_transform=TRANSFORM, dst_crs="EPSG:4326", resampling=Resampling.average,
              num_threads=2, warp_mem_limit=256)
    return result


def save_png(path, selected):
    with rasterio.open(path, "w", driver="PNG", width=WIDTH, height=HEIGHT,
                       count=1, dtype="uint8", crs="EPSG:4326", transform=TRANSFORM) as output:
        output.write(selected.astype(np.uint8) * 255, 1)
    # Standard world-file order; centre of upper-left pixel, no rotation.
    values = [TRANSFORM.a, 0, 0, TRANSFORM.e,
              TRANSFORM.c + TRANSFORM.a / 2, TRANSFORM.f + TRANSFORM.e / 2]
    path.with_suffix(".pgw").write_text("\n".join(f"{v:.15f}" for v in values) + "\n")
    path.with_suffix(".prj").write_text(rasterio.crs.CRS.from_epsg(4326).to_wkt() + "\n")


def prepare_masks():
    with zipfile.ZipFile(GSHHG) as archive:
        assert archive.testzip() is None, "GSHHG archive CRC failure"
    classes = np.zeros((SH, SW), dtype=np.uint8)
    for level in (1, 5):
        burn(f"GSHHS_shp/f/GSHHS_f_L{level}.shp", classes, SOURCE_TRANSFORM, 1)
    # Negative area identifies river-lake polygons, excluded from the lake mask.
    burn("GSHHS_shp/f/GSHHS_f_L2.shp", classes, SOURCE_TRANSFORM,
         lambda f: 1 if f["properties"]["area"] < 0 else 2)
    burn("GSHHS_shp/f/GSHHS_f_L3.shp", classes, SOURCE_TRANSFORM, 1)
    burn("GSHHS_shp/f/GSHHS_f_L4.shp", classes, SOURCE_TRANSFORM, 2)
    latitude = 90 - (np.arange(SH, dtype=np.float64) + 0.5) * 180 / SH
    weights = np.cos(np.deg2rad(latitude)).astype(np.float32)[:, None]
    total_weight = downsample(np.broadcast_to(weights, (SH, SW)).copy())
    shares = []
    for category in range(3):
        shares.append(downsample((classes == category).astype(np.float32) * weights))
    land = (shares[1] + shares[2]) >= total_weight * 0.5
    lakes = (shares[2] >= total_weight * 0.5) & land
    save_png(MASKS["land_ocean"], land)
    save_png(MASKS["lakes"], lakes)
    rivers = np.zeros((HEIGHT, WIDTH), dtype=np.uint8)
    for level in range(1, 12):
        burn(f"WDBII_shp/f/WDBII_river_f_L{level:02d}.shp", rivers, TRANSFORM, 1, True)
    river_coastal_pixels = int(np.count_nonzero((rivers > 0) & ~land))
    save_png(MASKS["rivers"], (rivers > 0) & land)
    log(f"Saved masks; clipped {river_coastal_pixels} river pixels outside coarse non-ocean mask")
    return classes, weights, shares, land, lakes, river_coastal_pixels


def prepare_height(classes, weights, shares, land, lakes):
    with rasterio.open(ETOPO) as source:
        assert (source.width, source.height) == (SW, SH)
        assert tuple(source.bounds) == BOUNDS
        assert source.transform.almost_equals(SOURCE_TRANSFORM)
        values = source.read(1)
        assert np.isfinite(values).all()
        if source.nodata is not None:
            assert not np.any(values == source.nodata)
        info = {"width": source.width, "height": source.height, "dtype": source.dtypes[0],
                "crs": str(source.crs), "bounds": list(source.bounds),
                "minimum_m": float(values.min()), "maximum_m": float(values.max())}
    output_class = np.where(~land, 0, np.where(lakes, 2, 1))
    height = np.full((HEIGHT, WIDTH), np.nan, dtype=np.float32)
    for category in range(3):
        log(f"Averaging ETOPO category {category} (0=ocean, 1=land, 2=lake)")
        weighted = values.astype(np.float32) * weights
        weighted[classes != category] = 0
        numerator = downsample(weighted)
        selected = output_class == category
        assert np.all(shares[category][selected] > 0)
        height[selected] = numerator[selected] / shares[category][selected]
        del weighted
    assert np.isfinite(height).all()
    with rasterio.open(HEIGHTMAP, "w", driver="GTiff", width=WIDTH, height=HEIGHT,
                       count=1, dtype="float32", crs="EPSG:4326", transform=TRANSFORM,
                       tiled=True, compress="deflate", predictor=3) as output:
        output.write(height, 1)
        output.set_band_description(1, "ETOPO 2022 surface elevation / bathymetry (metres, EGM2008)")
        output.set_band_unit(1, "m")
        output.update_tags(AREA_OR_POINT="Area", VERTICAL_DATUM="EGM2008 (EPSG:3855)",
                           SOURCE=URLS[ETOPO.name], MASK_SOURCE=URLS[GSHHG.name],
                           RESAMPLING="cos(latitude)-weighted mean within output surface class",
                           CLASSIFICATION="ocean, non-lake land, lake; 50 percent coverage masks")
    log("Saved float32 heightmap")
    return info


def verify():
    with rasterio.open(HEIGHTMAP) as source:
        assert (source.width, source.height, source.count) == (WIDTH, HEIGHT, 1)
        assert source.dtypes == ("float32",)
        assert source.crs.to_epsg() == 4326
        assert tuple(source.bounds) == BOUNDS and source.transform == TRANSFORM
        assert source.tags()["VERTICAL_DATUM"] == "EGM2008 (EPSG:3855)"
        assert source.units == ("m",)
        height = source.read(1)
    assert np.isfinite(height).all()
    assert -12000 < height.min() <= height.max() < 9000
    if WIDTH == 2048:
        assert height.min() < -7000 and 4500 < height.max()
    masks = {}
    for name, path in MASKS.items():
        with Image.open(path) as image:
            assert image.mode == "L" and image.size == (WIDTH, HEIGHT)
            values = np.asarray(image)
        assert set(np.unique(values)).issubset({0, 255})
        if name == "land_ocean":
            assert set(np.unique(values)) == {0, 255}
        masks[name] = values > 0
        with rasterio.open(path) as source:
            assert source.transform == TRANSFORM and source.crs.to_epsg() == 4326
    land, lakes, rivers = (masks[k] for k in ("land_ocean", "lakes", "rivers"))
    assert not np.any(lakes & ~land) and not np.any(rivers & ~land)
    assert land[-1].all(), "Antarctic polar cap must cover every longitude"
    assert not land[0].any(), "Northernmost row must be Arctic Ocean"
    row_weight = np.cos(np.deg2rad(90 - (np.arange(HEIGHT) + 0.5) * 180 / HEIGHT))[:, None]
    non_ocean_area = float((land * row_weight).sum() / (WIDTH * row_weight.sum()))
    assert 0.27 < non_ocean_area < 0.33
    samples = []
    points = [
        ("Sahara", 15, 25, True, False, 0, 2000),
        ("Greenland interior", -42, 75, True, False, 1000, 4000),
        ("Antarctic interior", 0, -89, True, False, 1000, 4500),
        ("Arctic Ocean", 0, 89, False, False, -6000, -100),
        ("Pacific Ocean", -140, 0, False, False, -7000, -1000),
        ("Himalayas", 86, 28, True, False, 2000, 9000),
        ("Caspian Sea", 51, 42, True, True, -2000, 100),
        ("Lake Superior", -87, 47.7, True, True, -500, 700),
        ("Lake Victoria", 33, -1, True, True, 500, 2000),
    ]
    # Point-sized features may disappear at other resolutions; keep the original
    # geographic fixtures for their validated 2048-column product.
    for name, lon, lat, expected_land, expected_lake, low, high in (points if WIDTH == 2048 else []):
        x, y = int((lon + 180) / 360 * WIDTH), int((90 - lat) / 180 * HEIGHT)
        assert bool(land[y, x]) == expected_land, f"Land classification: {name}"
        assert bool(lakes[y, x]) == expected_lake, f"Lake classification: {name}"
        elevation = float(height[y, x])
        assert low < elevation < high, f"Elevation: {name} = {elevation}"
        samples.append({"name": name, "longitude": lon, "latitude": lat,
                        "elevation_m": elevation, "non_ocean": bool(land[y, x]),
                        "lake": bool(lakes[y, x])})
    river_points = (("Amazon", -60, -3), ("Nile", 31, 26), ("Ganges", 88, 24)) if WIDTH == 2048 else ()
    for name, lon, lat in river_points:
        x, y = int((lon + 180) / 360 * WIDTH), int((90 - lat) / 180 * HEIGHT)
        assert rivers[y-3:y+4, x-3:x+4].any(), f"Missing river near {name}"
    result = {"dimensions": [WIDTH, HEIGHT], "elevation_dtype": "float32",
              "elevation_minimum_m": float(height.min()), "elevation_maximum_m": float(height.max()),
              "finite_elevation_pixels": int(np.isfinite(height).sum()),
              "white_pixels": {name: int(mask.sum()) for name, mask in masks.items()},
              "non_ocean_spherical_area_fraction": non_ocean_area, "geographic_checks": samples,
              "geographic_fixture_resolution": 2048,
              "geographic_fixtures_checked": WIDTH == 2048}
    METADATA.mkdir(parents=True, exist_ok=True)
    (METADATA / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    log(json.dumps({k: v for k, v in result.items() if k != "geographic_checks"}, indent=2))
    return result


def preview():
    from PIL import ImageDraw, ImageFont
    with rasterio.open(HEIGHTMAP) as source:
        height = source.read(1)
    with Image.open(MASKS["land_ocean"]) as image:
        land = np.asarray(image) > 0
    with Image.open(MASKS["lakes"]) as image:
        lakes = np.asarray(image) > 0
    rgb = np.empty((HEIGHT, WIDTH, 3), dtype=np.uint8)
    for channel in range(3):
        water = np.interp(height, [-11000, 0], [(5, 28, 60)[channel], (86, 155, 177)[channel]])
        ground = np.interp(height, [-500, 0, 1500, 4000, 6500],
                           [(118, 142, 91)[channel], (151, 172, 108)[channel],
                            (158, 143, 105)[channel], (186, 169, 141)[channel], (245, 244, 239)[channel]])
        rgb[:, :, channel] = np.where(land, ground, water).astype(np.uint8)
    rgb[lakes] = (70, 144, 182)
    Image.fromarray(rgb).save(METADATA / "earth_height_preview.png")
    sheet = Image.new("RGB", (2080, 1136), "#151b24")
    draw = ImageDraw.Draw(sheet)
    font_path = Path("C:/Windows/Fonts/segoeui.ttf")
    font = ImageFont.truetype(str(font_path), 22) if font_path.exists() else ImageFont.load_default()
    panels = [("Elevation / bathymetry (colour preview)", METADATA / "earth_height_preview.png"),
              ("White = non-ocean, including inland lakes", MASKS["land_ocean"]),
              ("White = lakes (majority cell coverage)", MASKS["lakes"]),
              ("White = river / waterway presence", MASKS["rivers"])]
    for i, (label, path) in enumerate(panels):
        x, y = 8 + (i % 2) * 1040, 8 + (i // 2) * 564
        draw.text((x, y), label, fill="white", font=font)
        with Image.open(path) as image:
            sheet.paste(image.convert("RGB").resize((1024, 512), Image.Resampling.NEAREST), (x, y + 36))
    sheet.save(METADATA / "preview.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--width", type=int, default=2048, help="Even output width, at least 64; height is width / 2")
    args = parser.parse_args()
    if args.width < 64 or args.width % 2:
        parser.error("--width must be even and at least 64")
    configure_grid(args.width)
    if args.verify_only:
        verify()
        manifest_path = METADATA / "manifest.json"
        if not manifest_path.exists():
            log("Structural/geographic checks passed. Hash verification unavailable: source/output manifest is missing.")
            return
        manifest = json.loads(manifest_path.read_text())
        for record in manifest["sources"] + manifest["outputs"]:
            assert describe(ROOT / record["file"])["sha256"] == record["sha256"], record["file"]
        log("Source and output SHA256 checks passed")
        return
    assert ETOPO.exists() and GSHHG.exists(), "Download the source files listed in README.md first"
    for path in (HEIGHTMAP, *MASKS.values(), METADATA / "manifest.json"):
        path.parent.mkdir(parents=True, exist_ok=True)
    log("Generating masks from full-resolution GSHHG on the ETOPO source grid")
    classes, weights, shares, land, lakes, clipped = prepare_masks()
    source_info = prepare_height(classes, weights, shares, land, lakes)
    del classes, weights, shares
    verification = verify()
    preview()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [dict(describe(path), url=URLS[path.name]) for path in (ETOPO, GSHHG)],
        "source_etopo_metadata": source_info,
        "grid": {"dimensions": [WIDTH, HEIGHT], "bounds": BOUNDS, "crs": "EPSG:4326",
                 "registration": "pixel/cell centre; north up", "vertical_datum": "EGM2008 (EPSG:3855)"},
        "processing": {"mask_source_grid": [SW, SH], "antarctica": "GSHHG L5 ice front",
                       "geometry_repair": "4-E polar closure restored via (180,-90) and (0,-90); coastline vertices unchanged",
                       "mask_threshold": 0.5, "weighting": "cos(source cell-centre latitude)",
                       "height": "weighted arithmetic mean within selected ocean/land/lake class",
                       "lakes": "L2 and L4 minus L3 islands; negative-area river-lakes excluded",
                       "rivers": "all WDBII source levels L01-L11; all_touched line rasterization",
                       "river_pixels_clipped_outside_non_ocean": clipped},
        "software": {"rasterio": rasterio.__version__, "gdal": rasterio.__gdal_version__,
                     "fiona": fiona.__version__, "numpy": np.__version__, "pillow": Image.__version__},
        "outputs": [describe(path) for path in (HEIGHTMAP, *MASKS.values())],
        "verification": verification,
    }
    (METADATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Finished; source URLs and hashes recorded in {METADATA / 'manifest.json'}")


if __name__ == "__main__":
    main()
