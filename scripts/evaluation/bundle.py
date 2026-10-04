"""Read run bundles written by world_pipeline (uw-run-bundle/v1).

A bundle is a directory with manifest.json, one NumPy file per field under fields/
and text tables under tables/. Grids are cell-centred, W x W/2, row-major with the
north row first (docs/architecture/world-model.md).
"""
import csv
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

SCHEMA = 'uw-run-bundle/v1'
SEASONS = ('January', 'April', 'July', 'October')


@dataclass
class FieldInfo:
    id: str
    component: str
    layout: str
    dtype: str
    shape: tuple
    units: str
    file: str
    hash: str
    written_by: str
    description: str


@dataclass
class Bundle:
    path: Path
    manifest: dict
    fields: dict = field(default_factory=dict)
    tables: dict = field(default_factory=dict)

    @classmethod
    def open(cls, path):
        path = Path(path)
        manifest = json.loads((path / 'manifest.json').read_text(encoding='utf-8'))
        if manifest.get('schema') != SCHEMA:
            raise ValueError(f'{path}: unsupported bundle schema {manifest.get("schema")!r}')
        bundle = cls(path, manifest)
        for entry in manifest['fields']:
            bundle.fields[entry['id']] = FieldInfo(entry['id'], entry['component'], entry['layout'], entry['dtype'],
                                                   tuple(entry['shape']), entry['units'], entry['file'], entry['hash'],
                                                   entry.get('written_by', ''), entry.get('description', ''))
        for entry in manifest['tables']:
            bundle.tables[entry['id']] = entry
        return bundle

    @property
    def seed(self):
        return self.manifest['seed']

    @property
    def columns(self):
        return self.manifest['grid']['columns']

    @property
    def rows(self):
        return self.manifest['grid']['rows']

    @property
    def recipe(self):
        return self.manifest.get('recipe', '')

    def has(self, identifier):
        return identifier in self.fields

    def load(self, identifier):
        """Field values as a NumPy array (scalars as 0-d arrays)."""
        info = self.fields[identifier]
        return np.load(self.path / info.file)

    def scalar(self, identifier):
        return self.load(identifier).item()

    def table_text(self, identifier):
        return (self.path / self.tables[identifier]['file']).read_text(encoding='utf-8')

    def table_rows(self, identifier):
        return list(csv.DictReader(io.StringIO(self.table_text(identifier))))

    def latitudes(self):
        """Cell-centre latitudes in degrees, north row first."""
        return 90.0 - (np.arange(self.rows) + 0.5) * 180.0 / self.rows

    def longitudes(self):
        return -180.0 + (np.arange(self.columns) + 0.5) * 360.0 / self.columns

    def area_weights(self):
        """Cell area weights (rows x columns), proportional to the cosine of latitude, summing to one."""
        weights = np.cos(np.deg2rad(self.latitudes()))[:, None] * np.ones((1, self.columns))
        return weights / weights.sum()

    def sea_mask(self):
        """True for sea cells of the base elevation (lakes count as land), as planet::sea without volcano islands."""
        elevation = self.load('geosphere.base_elevation')
        sea_level = self.scalar('geosphere.sea_level')
        lakes = self.load('hydrosphere.lake_surface') if self.has('hydrosphere.lake_surface') else 0
        return (elevation <= sea_level) & (np.asarray(lakes) == 0)

    def sea_fraction(self, rows, columns):
        """Area fraction of surface sea cells in each cell of a coarser rows x columns grid
        that divides the surface grid (solver grids), as the solvers' land fraction."""
        sea = self.sea_mask().astype(float)
        if (rows, columns) == sea.shape:
            return sea
        if self.rows % rows or self.columns % columns:
            raise ValueError(f'{columns} x {rows} does not divide the {self.columns} x {self.rows} surface grid')
        fy, fx = self.rows // rows, self.columns // columns
        weights = np.cos(np.deg2rad(self.latitudes()))[:, None] * np.ones((1, self.columns))
        total = (weights * sea).reshape(rows, fy, columns, fx).sum(axis=(1, 3))
        return total / weights.reshape(rows, fy, columns, fx).sum(axis=(1, 3))
