"""Prepare aligned, masked numerical fields for explicit benchmark interventions."""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import grid, references
from uwtools.hashing import sha256

ROOT = Path(__file__).resolve().parents[2]
A = ROOT / 'runs/reports/climate-reference-experiments-20260915'
REF = ROOT / 'refs/processed/climate'
OUT = A / 'prepared-inputs'
MISSING = np.float32(-9999.9)


def read(name):
    path = REF / f'{name}_monthly.uwclim'
    header, values = references.read(path, memmap=True)
    assert header.shape == (12, 1024, 2048)
    return values


def qsat(temperature, pressure):
    warm = 6.1121 * np.exp((18.678 - temperature / 234.5) * temperature / (257.14 + temperature))
    cold = 6.1115 * np.exp((23.036 - temperature / 333.7) * temperature / (279.82 + temperature))
    vapour = np.minimum(np.where(temperature >= 0, warm, cold), pressure * .95)
    return .622 * vapour / np.maximum(1, pressure - .378 * vapour)


def remap(field, width, extra_mask=None):
    height = width // 2
    ny, nx = field.shape
    factor = nx // width
    assert (ny, nx) == (height * factor, width * factor)
    valid = np.isfinite(field) & (field != MISSING)
    if extra_mask is not None:
        valid &= extra_mask
    weights = grid.band_area(ny)[:, None]
    numerator = (np.where(valid, field, 0) * weights).reshape(height, factor, width, factor).sum(axis=(1, 3))
    denominator = (valid * weights).reshape(height, factor, width, factor).sum(axis=(1, 3))
    total = np.broadcast_to(weights, field.shape).reshape(height, factor, width, factor).sum(axis=(1, 3))
    # At least 5% valid source coverage; native sea/land mask still gates application.
    return np.where(denominator >= .05 * total, numerator / np.maximum(denominator, 1e-30), MISSING).astype('<f4')


def main():
    OUT.mkdir(exist_ok=False)
    source_names = {'qa': ['era5_d2m', 'era5_sp', 'era5_sst'],
        'evap': ['era5_mer', 'era5_sst'], 'sst': ['era5_sst'], 'slp': ['era5_slp_anom'],
        'u850': ['era5_u850'], 'v850': ['era5_v850'], 'u10': ['era5_u10m'], 'v10': ['era5_v10m'],
        'u500': ['era5_u500'], 'v500': ['era5_v500'], 't2m': ['era5_t2m'], 'sp': ['era5_sp'],
        'rad': ['era5_msnswrf', 'era5_msnlwrf']}
    units = {'qa':'kg/kg', 'evap':'mm/day upward', 'sst':'C', 'slp':'hPa anomaly',
        't2m':'C', 'sp':'hPa', 'rad':'W/m2 downward net'}
    receipts = []
    for role, sources in source_names.items():
        arrays = [read(name) for name in sources]
        widths = [1024, 256, 128] if role == 'sst' else [256, 128]
        outputs = {width: np.empty((12, width // 2, width), dtype='<f4') for width in widths}
        for month in range(12):
            field = np.asarray(arrays[0][month], dtype=np.float64)
            mask = None
            if role == 'qa':
                pressure = np.asarray(arrays[1][month], dtype=np.float64)
                field = qsat(field, pressure)
                mask = np.isfinite(arrays[2][month]) & (arrays[2][month] != MISSING)
                assert np.all((field[mask] >= 0) & (field[mask] < .1))
            elif role == 'evap':
                mask = np.isfinite(arrays[1][month]) & (arrays[1][month] != MISSING)
            elif role == 'rad':
                field = field + arrays[1][month]
            for width in widths:
                outputs[width][month] = remap(field, width, mask)
        for width, values in outputs.items():
            target = OUT / str(width) / f'{role}.uwclim'
            target.parent.mkdir(exist_ok=True)
            references.write(target, role, values, mode='xb')
            receipts.append(dict(path=target.relative_to(ROOT).as_posix(), role=role,
                units=units.get(role, 'm/s east' if role.startswith('u') else 'm/s north'),
                sha256=sha256(target),
                sources=[(REF / f'{name}_monthly.uwclim').relative_to(ROOT).as_posix() for name in sources],
                width=width, height=width//2,
                valid_fractions=[float(np.mean(m != MISSING)) for m in values]))
        print(f'Prepared {role}: widths {widths}', flush=True)
    source_hashes = {}
    for name in sorted({name for sources in source_names.values() for name in sources}):
        path = REF / f'{name}_monthly.uwclim'
        source_hashes[path.relative_to(ROOT).as_posix()] = sha256(path)
    (OUT / 'manifest.json').write_text(json.dumps(dict(
        period='ERA5 2001-2020 monthly climatology',
        remap='Spherical-area conservative restriction of valid cells; minimum 5% source coverage; qa and evaporation use ERA5 SST coverage.',
        timing='Inputs retain 12 months; native interventions sample Jan/Apr/Jul/Oct with existing seasonal interpolation.',
        humidity='Buck saturation formula, matching climate_physics; q derived from mean dewpoint and pressure before remapping. Monthly-mean nonlinear approximation remains.',
        source_hashes=source_hashes, products=receipts), indent=2) + '\n')


if __name__ == '__main__':
    main()
