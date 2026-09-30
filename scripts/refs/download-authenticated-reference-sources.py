# /// script
# dependencies = ["copernicusmarine>=2.4", "earthaccess", "xarray", "netcdf4"]
# ///
"""Download GLORYS surface ocean fields or CERES EBAF using local provider logins."""

import argparse
import hashlib
import json
from pathlib import Path

import xarray as xr


def record(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            digest.update(block)
    return dict(file=path.name, bytes=path.stat().st_size, sha256=digest.hexdigest())


def marine(root, start, end):
    import copernicusmarine
    directory = root / f"glorys12-{start}-{end}"
    directory.mkdir(parents=True, exist_ok=True)
    for year in range(start, end+1):
        path = directory / f"glorys_surface_monthly_{year}.nc"
        receipt = path.with_suffix(".json")
        request = dict(dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
                       variables=["uo", "vo", "thetao", "mlotst"],
                       minimum_depth=0, maximum_depth=1,
                       start_datetime=f"{year}-01-01", end_datetime=f"{year}-12-31T23:59:59")
        if receipt.exists() and path.exists():
            cached = json.loads(receipt.read_text())
            if cached["request"] != request or cached["output"] != record(path):
                raise ValueError(f"Cached source differs: {path}")
            print(f"Verified {path.name}", flush=True)
            continue
        partial = path.with_suffix(".part.nc")
        copernicusmarine.subset(**request, output_directory=directory,
                               output_filename=partial.name, overwrite=True,
                               disable_progress_bar=True, netcdf_compression_level=1)
        with xr.open_dataset(partial) as data:
            if data.sizes.get("time") != 12 or data.sizes.get("depth") != 1:
                raise ValueError("Expected twelve months at one surface depth")
            if not set(request["variables"]).issubset(data.data_vars):
                raise ValueError("Missing requested ocean variables")
        partial.replace(path)
        receipt.write_text(json.dumps(dict(request=request, output=record(path)), indent=2) + "\n")
        print(f"Downloaded {path}", flush=True)


def ceres(root, start, end):
    import earthaccess
    auth = earthaccess.login(strategy="netrc")
    if not auth.authenticated:
        raise RuntimeError("Earthdata local login is not configured")
    directory = root / "ceres-ebaf-edition4.2.1"
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / "download.json"
    if receipt.exists():
        cached = json.loads(receipt.read_text())
        if cached["output"] != record(directory / cached["output"]["file"]):
            raise ValueError("CERES source checksum mismatch")
        print("Verified cached CERES EBAF", flush=True)
        return
    results = earthaccess.search_data(short_name="CERES_EBAF", version="Edition4.2.1",
                                     temporal=(f"{start}-01-01", f"{end}-12-31"))
    if not results:
        raise ValueError("No CERES EBAF Edition4.2.1 granule found")
    # Each EBAF granule is the complete monthly series. Keep the latest revision.
    selected = max(results, key=lambda g: g["umm"]["TemporalExtent"]["RangeDateTime"]["EndingDateTime"])
    paths = earthaccess.download([selected], local_path=directory, threads=1, show_progress=False)
    if len(paths) != 1:
        raise ValueError("CERES download did not return one source file")
    path = Path(paths[0])
    with xr.open_dataset(path) as data:
        if not data.data_vars:
            raise ValueError("Empty CERES source")
    receipt.write_text(json.dumps(dict(collection="CERES_EBAF", version="Edition4.2.1",
                                      granule=selected["umm"]["GranuleUR"],
                                      output=record(path)), indent=2) + "\n")
    print(f"Downloaded {path}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("marine", "ceres"), required=True)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--output", type=Path, default=Path("refs/source"))
    args = parser.parse_args()
    if args.end_year < args.start_year:
        parser.error("End year precedes start year")
    (marine if args.source == "marine" else ceres)(args.output, args.start_year, args.end_year)
