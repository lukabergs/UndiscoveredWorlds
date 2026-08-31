#!/usr/bin/env python3
"""Build 12-month climatologies from downloaded NSIDC G02202 monthly files."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import xarray as xr


SOURCE = "https://nsidc.org/data/g02202/versions/6"
ERDDAP = "https://coastwatch.pfeg.noaa.gov/erddap/"
VARIABLE = "cdr_seaice_conc_monthly"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def describe(path: Path) -> dict[str, object]:
    return {
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": file_sha256(path),
    }


def prepare(input_path: Path, output_path: Path, hemisphere: str) -> None:
    dataset = xr.open_dataset(input_path, engine="netcdf4", chunks={"time": 12})
    expected_months = 20 * 12
    if dataset.sizes.get("time") != expected_months:
        raise RuntimeError(
            f"Expected {expected_months} months in {input_path}, found {dataset.sizes.get('time')}"
        )

    climatology = dataset[VARIABLE].groupby("time.month").mean("time", keep_attrs=True)
    climatology.name = "sea_ice_concentration"
    climatology.attrs.update(
        {
            "climatology_period": "2001-01 through 2020-12",
            "hemisphere": hemisphere,
            "source_dataset": "NOAA/NSIDC Sea Ice Concentration CDR Version 6",
        }
    )
    result = climatology.to_dataset()
    result.attrs.update(dataset.attrs)
    result.attrs.update(
        {
            "title": f"NSIDC G02202 V6 {hemisphere} 12-month sea-ice concentration climatology",
            "processing": "Arithmetic mean of monthly concentration grouped by calendar month",
            "source": SOURCE,
        }
    )
    result.to_netcdf(
        output_path,
        engine="h5netcdf",
        encoding={
            "sea_ice_concentration": {
                "dtype": "float32",
                "compression": "gzip",
                "compression_opts": 4,
                "shuffle": True,
            }
        },
    )
    dataset.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory",
        type=Path,
        default=Path("extra/reference/nsidc-g02202-v6-2001-2020"),
    )
    args = parser.parse_args()

    files = []
    for hemisphere in ("nh", "sh"):
        input_path = args.directory / f"seaice_concentration_monthly_{hemisphere}_2001-2020.nc"
        output_path = args.directory / f"seaice_concentration_climatology_{hemisphere}_2001-2020.nc"
        prepare(input_path, output_path, hemisphere)
        files.extend((describe(input_path), describe(output_path)))

    for ancillary in sorted(args.directory.glob("G02202-ancillary-*.nc")):
        files.append(describe(ancillary))

    receipt = {
        "dataset": "NOAA/NSIDC G02202 Version 6",
        "doi": "10.7265/b18j-z797",
        "source": SOURCE,
        "download_service": ERDDAP,
        "period": "2001-01 through 2020-12",
        "downloaded_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "files": files,
    }
    (args.directory / "download-receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
