# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Validate reference bundle grids, provenance hashes and PNG/float32 TIFF exports."""

import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from PIL import Image


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            result.update(block)
    return result.hexdigest()


def verify(root, source_hashes=False):
    climate = root / "climate"
    bundles = list(climate.glob("*.uwclim"))
    for path in bundles:
        with path.open("rb") as stream:
            header = stream.read(40)
        if header[:8] != b"UWCLIM1\0":
            raise ValueError(f"Invalid bundle magic: {path}")
        version, width, height, layers = struct.unpack("<IIII", header[8:24])
        if version != 1 or width != 2*height or layers < 1 or path.stat().st_size != 40+width*height*layers*4:
            raise ValueError(f"Invalid bundle grid/payload: {path}")
    checked_sources = set()
    for metadata in climate.glob("*-additional-preparation.json"):
        receipt = json.loads(metadata.read_text())
        for product in receipt["products"]:
            path = climate / product["bundle"]
            if digest(path) != product["sha256"]:
                raise ValueError(f"Bundle checksum mismatch: {path}")
        if source_hashes:
            for source in receipt["sources"]:
                path = Path(source["file"])
                if path not in checked_sources and digest(path) != source["sha256"]:
                    raise ValueError(f"Source checksum mismatch: {path}")
                checked_sources.add(path)
    manifest = json.loads((root / "metadata/reference-maps/manifest.json").read_text())
    count = 0
    for resolution in manifest["outputs"]:
        width, height = resolution["width"], resolution["height"]
        if width != 2*height:
            raise ValueError("Manifest has a noncanonical grid")
        for record in (*resolution.get("land_maps", []), *resolution.get("sea_maps", [])):
            path = root / record["file"]
            if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
                raise ValueError(f"Import/preview checksum mismatch: {path}")
        names = set()
        for product in resolution["products"]:
            if product["name"] in names:
                raise ValueError("Duplicate reference product")
            names.add(product["name"])
            for record in (*product["maps"], product["csv"]):
                path = root / record["file"]
                if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
                    raise ValueError(f"Artifact checksum mismatch: {path}")
                if path.suffix not in (".png", ".tif"):
                    continue
                with Image.open(path) as image:
                    if image.size != (width, height):
                        raise ValueError(f"Wrong image dimensions: {path}")
                    if path.suffix == ".tif":
                        if image.mode != "F" or image.tag_v2[258] != (32,) or image.tag_v2[339] != (3,):
                            raise ValueError(f"Expected single-band float32 TIFF: {path}")
                        if image.tag_v2[33550] != (360/width, 180/height, 0):
                            raise ValueError(f"Wrong TIFF registration: {path}")
                        values = np.asarray(image)
                        valid = np.isfinite(values) & (values != np.float32(-9999.9))
                        if int(valid.sum()) != product["valid_cells"]:
                            raise ValueError(f"Incorrect numeric coverage: {path}")
                    elif product["name"] == "land_ocean":
                        if image.mode != "L" or not set(np.unique(image)).issubset({0, 255}):
                            raise ValueError(f"Land/ocean preview is not binary: {path}")
            count += 1
        print(f"Verified {len(names)} product sets at {width}x{height}", flush=True)
    print(f"Verified {len(bundles)} bundles, {count*3} map artifacts, {len(checked_sources)} source hashes", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("refs/processed"))
    parser.add_argument("--source-hashes", action="store_true")
    args = parser.parse_args()
    verify(args.root, args.source_hashes)
