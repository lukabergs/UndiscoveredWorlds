# /// script
# dependencies = ["cdsapi>=0.7.7", "requests", "xarray", "netcdf4", "filelock"]
# ///
"""Download original monthly climate sources, with resumable CDS request receipts."""

import argparse
import hashlib
import json
from pathlib import Path

import requests

MONTHLY = "reanalysis-era5-single-levels-monthly-means"
GROUPS = {
    "moisture": [
        "vertical_integral_of_eastward_water_vapour_flux",
        "vertical_integral_of_northward_water_vapour_flux",
        "vertical_integral_of_divergence_of_moisture_flux",
    ],
    "surface": [
        "sea_surface_temperature", "2m_temperature", "2m_dewpoint_temperature",
        "total_cloud_cover", "low_cloud_cover", "high_cloud_cover",
        "boundary_layer_height", "surface_pressure",
    ],
    "energy": [
        "mean_surface_net_short_wave_radiation_flux",
        "mean_surface_net_long_wave_radiation_flux",
        "mean_top_net_short_wave_radiation_flux",
        "mean_top_net_long_wave_radiation_flux",
        "mean_surface_sensible_heat_flux", "mean_surface_latent_heat_flux",
        "mean_evaporation_rate",
    ],
    "land": [
        "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2",
        "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4",
        "runoff", "surface_runoff", "sub_surface_runoff", "snow_depth",
        "snow_cover", "total_evaporation",
    ],
}


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def file_record(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return dict(file=path.name, bytes=path.stat().st_size, sha256=digest.hexdigest())


def http_download(url, params, partial):
    state_path = partial.with_suffix(partial.suffix + ".http.json")
    identity = {"url": url, "params": params}
    for attempt in range(4):
        state = json.loads(state_path.read_text()) if state_path.exists() else {}
        etag = state.get("etag") if state.get("request") == identity else None
        offset = partial.stat().st_size if etag and partial.exists() else 0
        headers = {"Range": f"bytes={offset}-", "If-Range": etag} if offset else {}
        try:
            with requests.get(url, params=params, headers=headers, stream=True, timeout=(30, 300)) as response:
                response.raise_for_status()
                resumed = response.status_code == 206
                if resumed and not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                    raise ValueError("Invalid HTTP resume offset")
                if resumed and (not etag or response.headers.get("ETag") != etag):
                    raise ValueError("HTTP source changed during resume")
                new_etag = response.headers.get("ETag", "")
                write_json(state_path, dict(request=identity,
                                           etag=new_etag if not new_etag.startswith("W/") else ""))
                with partial.open("ab" if resumed else "wb") as stream:
                    for block in response.iter_content(1024 * 1024):
                        stream.write(block)
            return
        except requests.RequestException:
            if attempt == 3:
                raise
            print(f"Retrying interrupted transfer ({attempt+1}/3)", flush=True)


def cds_download(group, root, start, end):
    import cdsapi
    dataset = "reanalysis-era5-land-monthly-means" if group == "land" else MONTHLY
    output = root / f"{'era5-land' if group == 'land' else 'era5-cds'}-{start}-{end}"
    output.mkdir(parents=True, exist_ok=True)
    request = dict(product_type="monthly_averaged_reanalysis", variable=GROUPS[group],
                   year=[str(y) for y in range(start, end+1)],
                   month=[f"{m:02d}" for m in range(1, 13)], time="00:00",
                   data_format="netcdf", download_format="unarchived")
    path = output / f"{group}_monthly_{start}-{end}.nc"
    receipt_path = output / f"{group}-download.json"
    receipt = dict(dataset=dataset, request=request, status="pending")
    if receipt_path.exists():
        previous = json.loads(receipt_path.read_text())
        if previous.get("request") != request or previous.get("dataset") != dataset:
            raise ValueError(f"Request differs from retained receipt: {receipt_path}")
        receipt = previous
        if receipt.get("status") == "downloaded":
            path = output / receipt["output"]["file"]
        if receipt.get("status") == "downloaded" and path.exists():
            if file_record(path) != receipt["output"]:
                raise ValueError(f"Source checksum mismatch: {path}")
            print(f"Verified cached {path}", flush=True)
            return
    client = cdsapi.Client(quiet=True, timeout=60, retry_max=3)
    if receipt.get("request_id"):
        remote = client.client.get_remote(receipt["request_id"])
    else:
        # Persist the request ID before waiting; reruns reuse the server job.
        remote = client.client.submit(dataset, request)
        receipt["request_id"] = remote.request_id
        write_json(receipt_path, receipt)
    print(f"CDS {group}: request {receipt['request_id']}", flush=True)
    partial = path.with_suffix(".nc.part")
    remote.download(str(partial))
    # CDS may package mixed step types in a ZIP even when unarchived was requested.
    import zipfile
    if zipfile.is_zipfile(partial):
        with zipfile.ZipFile(partial) as archive:
            if archive.testzip() is not None:
                raise ValueError("CDS archive CRC failure")
        path = path.with_suffix(".zip")
    else:
        import xarray as xr
        with xr.open_dataset(partial, engine="netcdf4") as data:
            if not data.data_vars:
                raise ValueError("Empty CDS dataset")
    partial.replace(path)
    receipt.update(status="downloaded", output=file_record(path))
    write_json(receipt_path, receipt)
    print(f"Downloaded {path} ({path.stat().st_size / 1e6:.1f} MB)", flush=True)


def download_oisst(root, start, end):
    output = root / f"noaa-oisst-v2.1-{start}-{end}"
    output.mkdir(parents=True, exist_ok=True)
    url = "https://psl.noaa.gov/thredds/ncss/grid/Datasets/noaa.oisst.v2.highres/sst.mon.mean.nc"
    # Monthly subsets avoid the server dropping larger annual responses.
    def transfer(date):
        year, month = date
        path = output / f"oisst_sst_monthly_{year}-{month:02d}.nc"
        receipt_path = path.with_suffix(".json")
        if path.exists() and receipt_path.exists():
            previous = json.loads(receipt_path.read_text())
            if file_record(path) == previous["output"]:
                return path, receipt_path, None
            raise ValueError(f"Source checksum mismatch: {path}")
        params = dict(var="sst", north=90, south=-90, west=0, east=360,
                      time=f"{year}-{month:02d}-01T00:00:00Z", accept="netcdf4")
        partial = path.with_suffix(".nc.part")
        http_download(url, params, partial)
        return path, receipt_path, params

    from concurrent.futures import ThreadPoolExecutor
    import xarray as xr
    dates = [(year, month) for year in range(start, end+1) for month in range(1, 13)]
    with ThreadPoolExecutor(max_workers=3) as workers:
        for path, receipt_path, params in workers.map(transfer, dates):
            if params is not None:
                partial = path.with_suffix(".nc.part")
                # netCDF/HDF5 validation stays on the main thread.
                with xr.open_dataset(partial, engine="netcdf4") as data:
                    if data.sizes.get("time") != 1 or "sst" not in data:
                        raise ValueError(f"Incomplete OISST month: {path.name}")
                    if str(data.time.values[0])[:7] != params["time"][:7]:
                        raise ValueError(f"OISST server returned the wrong month: {path.name}")
                partial.replace(path)
                write_json(receipt_path, dict(source=url, request=params, output=file_record(path)))
            print(f"Verified OISST {path.stem}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=(*GROUPS, "oisst"), required=True)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument("--output", type=Path, default=Path("refs/source"))
    parser.add_argument("--prepare", action="store_true",
                        help="After download, ingest at 2048x1024 and refresh this family's existing reference maps.")
    args = parser.parse_args()
    if args.end_year < args.start_year:
        parser.error("End year precedes start year")
    if args.source == "oisst":
        download_oisst(args.output, args.start_year, args.end_year)
    else:
        cds_download(args.source, args.output, args.start_year, args.end_year)
    if args.prepare:
        import os
        import shutil
        import subprocess
        from filelock import FileLock
        repository = Path(__file__).resolve().parents[2]
        family = "oisst" if args.source == "oisst" else "land" if args.source == "land" else "era5"
        prefix = "noaa-oisst-v2.1" if family == "oisst" else "era5-land" if family == "land" else "era5-cds"
        source = args.output.resolve() / f"{prefix}-{args.start_year}-{args.end_year}"
        # uv run exposes its real executable; bare "uv" may only be a Windows
        # .cmd shim, which subprocess cannot resolve like PowerShell does.
        uv = os.environ.get("UV") or shutil.which("uv")
        if not uv:
            raise RuntimeError("Cannot locate uv for reference preparation")
        # Different provider jobs can finish together. Serialize ingestion and
        # manifest replacement; the operating system releases the lock on exit.
        lock = repository / "refs/processed/metadata/reference-maps/refresh.lock"
        lock.parent.mkdir(parents=True, exist_ok=True)
        with FileLock(lock):
            subprocess.run([uv, "run", "scripts/refs/prepare-additional-references.py",
                            "--input", str(source), "--family", family,
                            "--start-year", str(args.start_year), "--end-year", str(args.end_year)],
                           cwd=repository, check=True)
            subprocess.run([uv, "run", "--with", "numpy", "--with", "pillow", "python",
                            "scripts/refs/prepare-reduced-earth-benchmark.py", "--refresh-additional",
                            "--additional-family", family], cwd=repository, check=True)


if __name__ == "__main__":
    main()
