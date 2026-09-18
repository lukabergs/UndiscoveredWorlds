# /// script
# dependencies = ["numpy", "pillow", "xarray", "netcdf4"]
# ///
"""Ingest downloaded monthly NetCDF sources onto the canonical W x W/2 grid.

Reads original coordinates, including native pole samples, rather than resizing
old UWCLIM bundles. Metadata records provenance, units, time coverage and hashes.
"""

import argparse
import hashlib
import json
import struct
import zipfile
from pathlib import Path

import numpy as np
import xarray as xr

from cell_grid import conservative_remap, grid_dimensions

# source variable: (bundle prefix, stored variable, units, scale, offset,
#                   expected source units, display maximum, palette, ocean only)
ERA5 = {
    "viwve": ("era5", "viwve", "kg m^-1 s^-1", 1, 0, ("kg m**-1 s**-1",), 500, "divergence", False),
    "viwvn": ("era5", "viwvn", "kg m^-1 s^-1", 1, 0, ("kg m**-1 s**-1",), 500, "divergence", False),
    "vimd": ("era5", "vimd", "mm/day convergence", -86400, 0, ("kg m**-2 s**-1",), 20, "convergence", False),
    "sst": ("era5", "sst", "degrees C", 1, -273.15, ("K",), 35, "divergence", True),
    "t2m": ("era5", "t2m", "degrees C", 1, -273.15, ("K",), 35, "divergence", False),
    "d2m": ("era5", "d2m", "degrees C", 1, -273.15, ("K",), 35, "divergence", False),
    "tcc": ("era5", "tcc", "fraction", 1, 0, ("(0 - 1)", "1"), 1, "speed", False),
    "lcc": ("era5", "lcc", "fraction", 1, 0, ("(0 - 1)", "1"), 1, "speed", False),
    "hcc": ("era5", "hcc", "fraction", 1, 0, ("(0 - 1)", "1"), 1, "speed", False),
    "blh": ("era5", "blh", "m", 1, 0, ("m",), 3000, "speed", False),
    "sp": ("era5", "sp", "hPa", .01, 0, ("Pa",), 1050, "speed", False),
    **{key: ("era5", key, "W/m^2 downward positive", 1, 0, ("W m**-2",), 300, "convergence", False)
       for key in ("msnswrf", "msnlwrf", "mtnswrf", "mtnlwrf", "msshf", "mslhf")},
    "mer": ("era5", "mer", "mm/day upward positive", -86400, 0, ("kg m**-2 s**-1",), 10, "speed", False),
}
LAND = {
    **{f"swvl{i}": ("era5land", f"swvl{i}", "m^3/m^3", 1, 0, ("m**3 m**-3",), .5, "speed", False)
       for i in range(1, 5)},
    **{key: ("era5land", key, "mm/day", 1000, 0, ("m", "m of water equivalent"), 10, "speed", False)
       for key in ("ro", "sro", "ssro")},
    "e": ("era5land", "e", "mm/day upward positive", -1000, 0, ("m of water equivalent", "m"), 10, "speed", False),
    "sd": ("era5land", "sd", "mm water equivalent", 1000, 0, ("m of water equivalent",), 1000, "speed", False),
    "sde": ("era5land", "snow_depth", "m snow depth", 1, 0, ("m",), 5, "speed", False),
    "snowc": ("era5land", "snowc", "fraction", .01, 0, ("%",), 1, "speed", False),
}
OISST = {"sst": ("oisst", "sst", "degrees C", 1, 0, ("degC", "degree_Celsius", "degrees C"), 35, "divergence", True)}
MARINE = {
    **{key: ("glorys", key, "m/s", 1, 0, ("m s-1", "m s^-1", "m/s"), 2, "divergence", True)
       for key in ("uo", "vo")},
    "thetao": ("glorys", "thetao", "degrees C", 1, 0, ("degrees_C", "degree_Celsius", "degrees_Celsius"), 35, "divergence", True),
    "mlotst": ("glorys", "mlotst", "m", 1, 0, ("m",), 500, "speed", True),
}
ERA5_ALIASES = dict(zip(
    ("avg_snswrf", "avg_snlwrf", "avg_tnswrf", "avg_tnlwrf", "avg_ishf", "avg_slhtf", "avg_ie"),
    ("msnswrf", "msnlwrf", "mtnswrf", "mtnlwrf", "msshf", "mslhf", "mer")))
ERA5_ALIASES["vimdf"] = "vimd"
CDS_GROUP_VARIABLES = {
    "moisture": {"viwve", "viwvn", "vimd"},
    "surface": {"sst", "t2m", "d2m", "tcc", "lcc", "hcc", "blh", "sp"},
    "energy": {"msnswrf", "msnlwrf", "mtnswrf", "mtnlwrf", "msshf", "mslhf", "mer"},
    "land": set(LAND) - {"sd"},  # The CDS snow_depth request supplies geometric depth (sde).
}
CERES = {
    **{key + "_mon": ("ceres", key, "W/m^2 " + direction, 1, 0, ("W m-2",), 400, "convergence", False)
       for key, direction in (
           ("toa_sw_all", "outgoing"), ("toa_lw_all", "outgoing"), ("toa_net_all", "downward positive"),
           ("solar", "incoming"), ("sfc_sw_down_all", "downward"), ("sfc_sw_up_all", "upward"),
           ("sfc_lw_down_all", "downward"), ("sfc_lw_up_all", "upward"),
           ("sfc_net_sw_all", "downward positive"), ("sfc_net_lw_all", "downward positive"),
           ("sfc_net_tot_all", "downward positive"),
           ("toa_cre_sw", "cloud radiative effect"), ("toa_cre_lw", "cloud radiative effect"),
           ("toa_cre_net", "cloud radiative effect"), ("sfc_cre_net_sw", "cloud radiative effect"),
           ("sfc_cre_net_lw", "cloud radiative effect"), ("sfc_cre_net_tot", "cloud radiative effect"))},
    "cldarea_total_daynight_mon": ("ceres", "cloud_fraction", "fraction", .01, 0, ("percent",), 1, "speed", False),
    "cldpress_total_daynight_mon": ("ceres", "cloud_pressure", "hPa", 1, 0, ("hPa",), 1000, "speed", False),
    "cldtemp_total_daynight_mon": ("ceres", "cloud_temp", "degrees C", 1, -273.15, ("K",), 60, "divergence", False),
    "cldtau_total_day_mon": ("ceres", "cloud_tau", "dimensionless", 1, 0, ("dimensionless",), 30, "speed", False),
}


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def normalize(data):
    rename = {old: new for old, new in (("valid_time", "time"), ("latitude", "lat"), ("longitude", "lon"))
              if old in data.dims}
    data = data.rename(rename)
    for dim in set(data.dims) - {"time", "lat", "lon"}:
        if data.sizes[dim] != 1:
            raise ValueError(f"Ambiguous dimension {dim}: select a level/member explicitly")
        data = data.isel({dim: 0}, drop=True)
    if set(data.dims) != {"time", "lat", "lon"}:
        raise ValueError(f"Expected time/latitude/longitude field, got {data.dims}")
    return data.transpose("time", "lat", "lon")


def monthly_climatology(arrays, start, end):
    arrays = [normalize(a).sel(time=slice(f"{start}-01-01", f"{end}-12-31")) for a in arrays]
    arrays = [a for a in arrays if a.sizes["time"]]
    if not arrays:
        raise ValueError("No data within the requested reference period")
    expected = {(y, m) for y in range(start, end+1) for m in range(1, 13)}
    seen = set()
    for a in arrays:
        if not np.array_equal(a.lat, arrays[0].lat) or not np.array_equal(a.lon, arrays[0].lon):
            raise ValueError("Source grids differ between files")
        for y, m in zip(a.time.dt.year.values, a.time.dt.month.values):
            key = (int(y), int(m))
            if key in seen:
                raise ValueError(f"Duplicate monthly sample {key}; hourly/daily data require temporal aggregation first")
            seen.add(key)
    if seen != expected:
        raise ValueError(f"Incomplete monthly coverage; missing {sorted(expected-seen)[:12]}")
    total = np.zeros((12, *arrays[0].shape[1:]), dtype=np.float64)
    count = np.zeros(total.shape, dtype=np.uint16)
    for a in arrays:
        # CDS chunks contain many times. Reading one slice repeatedly decompresses
        # the same chunk; sequential bounded blocks avoid that amplification.
        batch = max(1, min(60, (1024 * 1024**2) // (a.shape[1] * a.shape[2] * 4)))
        chunks = a.encoding.get("preferred_chunks", {})
        time_chunk = chunks.get("time", chunks.get("valid_time", 1))
        if time_chunk <= batch:
            batch = max(time_chunk, batch // time_chunk * time_chunk)
        months = a.time.dt.month.values
        for first in range(0, a.sizes["time"], batch):
            block = np.asarray(a.isel(time=slice(first, first+batch)), dtype=np.float32)
            for index, values in enumerate(block):
                month = int(months[first+index]) - 1
                valid = np.isfinite(values)
                total[month] += np.where(valid, values, 0)
                count[month] += valid
    for month in range(1, 13):
        mean = np.divide(total[month-1], count[month-1],
                         out=np.full_like(total[month-1], np.nan), where=count[month-1] > 0)
        yield month, mean.astype(np.float32), np.asarray(arrays[0].lat), np.asarray(arrays[0].lon)


def netcdf_paths(directory):
    paths = sorted(p for p in directory.glob("*.nc") if not p.name.endswith(".part.nc"))
    for archive_path in sorted(directory.glob("*.zip")):
        destination = directory / archive_path.stem
        destination.mkdir(exist_ok=True)
        with zipfile.ZipFile(archive_path) as archive:
            for info in archive.infolist():
                if not info.filename.endswith(".nc"):
                    continue
                # Flatten only known NetCDF members; never extract archive paths.
                target = destination / Path(info.filename).name
                if not target.exists() or target.stat().st_size != info.file_size:
                    import shutil
                    partial = target.with_suffix(".nc.part")
                    with archive.open(info) as src, partial.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                    partial.replace(target)
                paths.append(target)
    return paths


def prepare(directory, output, family, width, start, end):
    width, height = grid_dimensions(width)
    catalog = {"era5": ERA5, "land": LAND, "oisst": OISST, "marine": MARINE, "ceres": CERES}[family]
    paths = netcdf_paths(directory)
    if not paths:
        raise ValueError(f"No downloaded NetCDF sources in {directory}")
    output.mkdir(parents=True, exist_ok=True)
    from contextlib import ExitStack
    with ExitStack() as stack:
        datasets = [(p, stack.enter_context(xr.open_dataset(p))) for p in paths]
        if family == "era5":
            datasets = [(p, ds.rename({old: new for old, new in ERA5_ALIASES.items() if old in ds}))
                        for p, ds in datasets]
        for group, expected in CDS_GROUP_VARIABLES.items():
            receipt_path = directory / f"{group}-download.json"
            if not receipt_path.exists():
                continue
            receipt = json.loads(receipt_path.read_text())
            if receipt.get("status") != "downloaded":
                continue
            source_name = receipt["output"]["file"]
            source_stem = Path(source_name).stem
            supplied = {key for p, ds in datasets
                        if p.name == source_name or p.parent.name == source_stem for key in ds.data_vars}
            if expected - supplied:
                raise ValueError(f"Downloaded {group} source lacks expected variables: {sorted(expected-supplied)}")
        records = []
        for key, spec in catalog.items():
            selected = [(p, ds[key]) for p, ds in datasets if key in ds]
            if not selected:
                continue
            prefix, variable, unit, scale, offset, expected_units, limit, palette, ocean = spec
            if len(variable.encode("ascii")) > 16:
                raise ValueError(f"UWCLIM variable exceeds sixteen bytes: {variable}")
            for p, a in selected:
                if a.attrs.get("units") not in expected_units:
                    raise ValueError(f"Unexpected units for {key} in {p}: {a.attrs.get('units')}; expected {expected_units}")
            bundle = output / f"{prefix}_{variable}_monthly.uwclim"
            temporary = bundle.with_suffix(".uwclim.part")
            with temporary.open("wb") as stream:
                stream.write(b"UWCLIM1\0" + struct.pack("<IIII", 1, width, height, 12))
                stream.write(variable.encode("ascii").ljust(16, b"\0"))
                for month, values, lat, lon in monthly_climatology([a for _, a in selected], start, end):
                    values = conservative_remap(values * scale + offset, width,
                                                source_latitudes=lat, source_longitudes=lon)
                    stream.write(values.astype("<f4").tobytes())
            if temporary.stat().st_size != 40 + width * height * 12 * 4:
                raise ValueError("Incomplete UWCLIM write")
            temporary.replace(bundle)
            records.append(dict(bundle=bundle.name, variable=variable, units=unit, source_variable=key,
                                description=selected[0][1].attrs.get("long_name", key),
                                source_depth_m=(selected[0][1].depth.values.tolist()
                                                if "depth" in selected[0][1].coords else None),
                                sources=[p.as_posix() for p, _ in selected], sha256=digest(bundle),
                                limit=limit, palette=palette, ocean_only=ocean, land_only=family == "land",
                                period=f"{start}-{end} monthly climatology", width=width, height=height))
            print(f"Prepared {bundle.name}: {width}x{height}, {unit}", flush=True)
    if not records:
        raise ValueError("No supported variables in the downloaded sources")
    metadata = output / f"{family}-additional-preparation.json"
    temporary = metadata.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(dict(
        family=family, width=width, height=height,
        grid="cell-centred W x W/2; north to south; no stored pole rows",
        sources=[dict(file=p.as_posix(), sha256=digest(p)) for p in paths],
        products=records), indent=2) + "\n", encoding="utf-8")
    temporary.replace(metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--family", choices=("era5", "land", "oisst", "marine", "ceres"), required=True)
    parser.add_argument("--output", type=Path, default=Path("refs/processed/climate"))
    parser.add_argument("--width", type=int, default=2048)
    parser.add_argument("--height", type=int)
    parser.add_argument("--start-year", type=int, default=2001)
    parser.add_argument("--end-year", type=int, default=2020)
    args = parser.parse_args()
    grid_dimensions(args.width, args.height)
    if args.end_year < args.start_year:
        parser.error("End year precedes start year")
    prepare(args.input, args.output, args.family, args.width, args.start_year, args.end_year)


if __name__ == "__main__":
    main()
