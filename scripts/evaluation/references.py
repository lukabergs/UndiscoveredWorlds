"""Observed reference fields for validation, reduced to model seasons and grids.

References are processed .uwclim files under refs/processed/climate/ (cell-centred,
north row first). They are validation inputs only: nothing here feeds a simulation.
Season rules and unit conversions follow the climate report (climate_report_data.py),
so per-field losses use the same definitions as global-balanced-v1.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIRECTORY = ROOT / 'refs/processed/climate'
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'scripts/refs'))
from cell_grid import conservative_remap  # noqa: E402
from uwtools import references as reference_files  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402

# Mean month lengths used by the climate report for day-weighted quarters and mm/day.
MONTH_DAYS = np.array([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
# Monthly references stored per month that the report converts to per day.
PER_DAY = {'era5_pr'}
SEASON_RULES = ('month', 'quarter', 'annual')


def read_uwclim(path):
    """(variable, values[layers, rows, columns]) as float64 with missing values as NaN."""
    header, values = reference_files.read(path, missing=np.nan, dtype=float)
    return header.variable, values


def remap(values, width):
    """Conservative remap of a 2-D field to width x width/2, as the climate report does."""
    if values.shape[-1] == width:
        return np.asarray(values, dtype=float)
    return conservative_remap(np.asarray(values, dtype=np.float32), width).astype(float)


class References:
    """Loads reference files once and returns season- and grid-reduced fields."""

    def __init__(self, directory=REFERENCE_DIRECTORY):
        self.directory = Path(directory)
        self.files = {}
        self.cache = {}

    def path(self, name, annual):
        return self.directory / (f'{name}.uwclim' if annual else f'{name}_monthly.uwclim')

    def raw(self, name, annual):
        key = (name, annual)
        if key not in self.cache:
            path = self.path(name, annual)
            if not path.is_file():
                raise FileNotFoundError(f'Missing reference {path.relative_to(ROOT)}; see refs/README.md')
            _, values = read_uwclim(path)
            if annual and values.shape[0] != 1:
                raise ValueError(f'{path}: annual reference has {values.shape[0]} layers')
            if not annual:
                if values.shape[0] != 12:
                    raise ValueError(f'{path}: monthly reference has {values.shape[0]} layers')
                if name in PER_DAY:
                    values = values / MONTH_DAYS[:, None, None]
            self.cache[key] = values
            self.files[path.relative_to(ROOT).as_posix()] = None
        return self.cache[key]

    def season(self, name, rule, season, width):
        """Reference field for one model season on a width x width/2 grid.

        month: the season's first month (January, April, July, October); quarter: the
        day-weighted mean of that month and the next two; annual: the annual file.
        """
        if rule not in SEASON_RULES:
            raise ValueError(f'Unknown season rule {rule!r}')
        key = (name, rule, season, width)
        if key not in self.cache:
            values = self.raw(name, rule == 'annual')
            if rule == 'annual':
                field = values[0]
            elif rule == 'month':
                field = values[season * 3]
            else:
                months = np.arange(season * 3, season * 3 + 3) % 12
                field = np.average(values[months], axis=0, weights=MONTH_DAYS[months])
            self.cache[key] = remap(field, width)
        return self.cache[key]

    def provenance(self):
        """sha256 of every reference file read so far."""
        for relative in self.files:
            if self.files[relative] is None:
                self.files[relative] = sha256(ROOT / relative)
        return dict(sorted(self.files.items()))
