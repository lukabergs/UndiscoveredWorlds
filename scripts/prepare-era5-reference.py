#!/usr/bin/env python3
"""Download ERA5 monthly fields and write 12-month climatology NetCDF files."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import icechunk as ic
import xarray as xr


DEFAULT_VARIABLES = (
    "slp",
    "tcwv",
    "u10m",
    "v10m",
    "u500",
    "v500",
    "u850",
    "v850",
    "w500",
    "pr",
)
BUCKET = "planette-era5"
REGION = "us-east-2"


def store_prefix(variable: str) -> str:
    return (
        f"era5/{variable}/month/0p25latx0p25lon/"
        f"era5_{variable}_month_0p25latx0p25lon.zarr"
    )


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_variable(variable: str, start_year: int, end_year: int, output: Path) -> dict[str, object]:
    prefix = store_prefix(variable)
    output_path = output / f"era5_{variable}_monthly_climatology_{start_year}-{end_year}.nc"
    if output_path.exists():
        return {
            "variable": variable,
            "file": output_path.name,
            "bytes": output_path.stat().st_size,
            "sha256": file_sha256(output_path),
            "source_store": f"s3://{BUCKET}/{prefix}",
        }

    storage = ic.s3_storage(
        bucket=BUCKET,
        prefix=prefix,
        region=REGION,
        anonymous=True,
    )
    repository = ic.Repository.open(storage=storage)
    session = repository.readonly_session("main")
    dataset = xr.open_dataset(
        session.store,
        engine="zarr",
        consolidated=False,
        decode_timedelta=True,
        chunks={},
    )

    data_variables = [
        name
        for name, values in dataset.data_vars.items()
        if {"time", "lat", "lon"}.issubset(values.dims)
    ]
    if len(data_variables) != 1:
        raise RuntimeError(f"Expected one data variable in {prefix}, found {data_variables}")

    source_name = data_variables[0]
    selected = dataset[source_name].sel(
        time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
    )
    if selected.sizes.get("time") != (end_year - start_year + 1) * 12:
        raise RuntimeError(
            f"Unexpected month count for {variable}: {selected.sizes.get('time')}"
        )

    climatology = selected.groupby("time.month").mean("time", keep_attrs=True)
    climatology.name = variable
    climatology.attrs.update(
        {
            "climatology_period": f"{start_year}-01 through {end_year}-12",
            "source_archive": "Planette ERA5 Archive derived from ECMWF/Copernicus ERA5",
            "source_store": f"s3://{BUCKET}/{prefix}",
        }
    )
    result = climatology.to_dataset()
    result.attrs.update(dataset.attrs)
    result.attrs.update(
        {
            "title": f"ERA5 {variable} 12-month climatology, {start_year}-{end_year}",
            "source_store": f"s3://{BUCKET}/{prefix}",
            "processing": "Arithmetic mean of the source monthly means grouped by calendar month",
        }
    )

    result.to_netcdf(
        output_path,
        engine="h5netcdf",
        encoding={
            variable: {
                "dtype": "float32",
                "compression": "gzip",
                "compression_opts": 4,
                "shuffle": True,
            }
        },
    )
    dataset.close()

    return {
        "variable": variable,
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "sha256": file_sha256(output_path),
        "source_store": f"s3://{BUCKET}/{prefix}",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2020)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("extra/reference/climate/sources/era5-planette-2001-2020"),
    )
    parser.add_argument("--variables", nargs="+", default=DEFAULT_VARIABLES)
    args = parser.parse_args()

    if args.end_year < args.start_year:
        parser.error("--end-year must not precede --start-year")

    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    for variable in args.variables:
        print(f"Preparing ERA5 {variable} climatology", flush=True)
        records.append(
            prepare_variable(variable, args.start_year, args.end_year, args.output)
        )

    receipt = {
        "dataset": "planette-era5-monthly-climatology",
        "source": "https://registry.opendata.aws/planette_era5_reanalysis/",
        "underlying_dataset": "ECMWF/Copernicus ERA5",
        "period": f"{args.start_year}-01 through {args.end_year}-12",
        "downloaded_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "files": records,
    }
    receipt_path = args.output / "download-receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
