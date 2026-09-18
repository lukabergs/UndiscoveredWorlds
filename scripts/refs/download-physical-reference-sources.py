# /// script
# dependencies = ["requests", "earthaccess"]
# ///
"""Download the bounded physical-reference catalog; retain hashes and provenance."""

import argparse
from functools import lru_cache
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
from urllib.parse import urlsplit, urlunsplit
import zipfile

import requests

REPO = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=2)
def session(authenticated):
    if authenticated:
        import earthaccess
        if not earthaccess.login(strategy="netrc").authenticated:
            raise ValueError("Earthdata local login is not configured")
        return earthaccess.get_requests_https_session()
    return requests.Session()


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def validate(path, kind):
    with path.open("rb") as stream:
        head = stream.read(512)
    if not head:
        raise ValueError("Empty download")
    if kind != "html" and (b"<html" in head.lower() or b"<!doctype html" in head.lower()):
        raise ValueError("Provider returned an HTML page instead of data")
    if kind in ("zip", "xlsx"):
        with zipfile.ZipFile(path) as archive:
            if not archive.namelist():
                raise ValueError("Empty ZIP archive")
            if kind == "xlsx" and "[Content_Types].xml" not in archive.namelist():
                raise ValueError("Not an XLSX workbook")
    elif kind == "tar":
        with tarfile.open(path) as archive:
            if archive.next() is None:
                raise ValueError("Empty TAR archive")
    elif kind == "pdf" and not head.startswith(b"%PDF-"):
        raise ValueError("Invalid PDF signature")
    elif kind == "netcdf" and not head.startswith((b"CDF", b"\x89HDF")):
        raise ValueError("Invalid NetCDF/HDF signature")
    elif kind == "hdf" and not head.startswith((b"\x0e\x03\x13\x01", b"\x89HDF")):
        raise ValueError("Invalid HDF signature")
    elif kind == "tiff" and not head.startswith((b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")):
        raise ValueError("Invalid TIFF signature")
    elif kind in ("json", "geojson"):
        data = json.loads(path.read_text(encoding="utf-8-sig"))
        if kind == "geojson" and not data.get("features"):
            raise ValueError("Empty GeoJSON feature collection")


def fetch(entry, root, budget, minimum_free, verify_only):
    path = root / entry["dataset"] / entry["file"]
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Catalog path escapes output directory")
    receipt = path.with_name(path.name + ".receipt.json")
    if path.exists():
        if not receipt.exists():
            raise ValueError("Existing file has no receipt; refusing to overwrite")
        saved = json.loads(receipt.read_text())
        if saved["request"] != entry or saved["sha256"] != digest(path):
            raise ValueError("Cached source or catalog changed")
        validate(path, entry["kind"])
        print(f"Verified {entry['dataset']}/{path.name}", flush=True)
        return path.stat().st_size
    if verify_only:
        raise FileNotFoundError("Source has not been downloaded")
    path.parent.mkdir(parents=True, exist_ok=True)
    used = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    limit = min(entry.get("max_mib", 512) * 1024**2, budget - used)
    if limit <= 0 or shutil.disk_usage(root).free < minimum_free:
        raise ValueError("Collection budget or free-space floor reached")
    partial = path.with_name(path.name + ".part")
    size = 0
    with session(entry.get("earthdata", False)).get(entry["url"], stream=True, timeout=(20, 60)) as response:
        if response.status_code != 200:
            raise ValueError(f"HTTP {response.status_code}")
        expected = int(response.headers.get("Content-Length", 0))
        if expected > limit:
            raise ValueError(f"Provider file exceeds limit: {expected} > {limit}")
        try:
            with partial.open("wb") as stream:
                for block in response.iter_content(1024**2):
                    size += len(block)
                    if size > limit or shutil.disk_usage(root).free - len(block) < minimum_free:
                        raise ValueError("Transfer exceeds storage limit")
                    stream.write(block)
            if expected and not response.headers.get("Content-Encoding") and size != expected:
                raise ValueError("Incomplete response body")
            validate(partial, entry["kind"])
            checksum = digest(partial)
            if entry.get("sha256") and entry["sha256"] != checksum:
                raise ValueError("Provider checksum mismatch")
            partial.replace(path)
        finally:
            partial.unlink(missing_ok=True)
        url = urlsplit(response.url)
        result = dict(request=entry, bytes=size, sha256=checksum,
                      retrieved_utc=datetime.now(timezone.utc).isoformat(),
                      final_url=urlunsplit((url.scheme, url.netloc, url.path, "", "")),
                      content_type=response.headers.get("Content-Type"),
                      etag=response.headers.get("ETag"))
    receipt.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Downloaded {entry['dataset']}/{path.name}: {size / 1024**2:.2f} MiB", flush=True)
    return size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=REPO / "configs/physical-reference-sources.json")
    parser.add_argument("--output", type=Path, default=REPO / "refs/source/physical")
    parser.add_argument("--dataset", action="append")
    parser.add_argument("--budget-gib", type=float, default=10)
    parser.add_argument("--minimum-free-gib", type=float, default=20)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    catalog = json.loads(args.catalog.read_text(encoding="utf-8"))
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for entry in catalog["files"]:
        if args.dataset and entry["dataset"] not in args.dataset:
            continue
        try:
            size = fetch(entry, args.output, int(args.budget_gib * 1024**3),
                         int(args.minimum_free_gib * 1024**3), args.verify_only)
            results.append(dict(dataset=entry["dataset"], file=entry["file"], status="verified", bytes=size))
        except Exception as exc:
            # Exception URLs may contain temporary provider credentials; omit their text.
            message = str(exc) if isinstance(exc, (ValueError, FileNotFoundError)) else type(exc).__name__
            results.append(dict(dataset=entry["dataset"], file=entry["file"], status="unavailable", error=message))
            print(f"Unavailable {entry['dataset']}/{entry['file']}: {message}", flush=True)
    report = args.output / "download-status.json"
    report.write_text(json.dumps(dict(checked_utc=datetime.now(timezone.utc).isoformat(),
                                    results=results), indent=2) + "\n", encoding="utf-8")
    return int(any(r["status"] != "verified" for r in results))


if __name__ == "__main__":
    raise SystemExit(main())
