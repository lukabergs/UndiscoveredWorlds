"""Read run bundles written by world_pipeline (uw-run-bundle/v2; v1 reads as grid kind latlon).

A bundle is a directory with manifest.json, one NumPy file per field under fields/
and text tables under tables/ (docs/architecture/world-model.md). Bundle v2
(docs/architecture/sphere-mesh.md section 7) has a grid kind:
- latlon: cell-centred W x W/2, row-major with the north row first; fields [layer,] rows, columns.
- cube: N cells per face edge (orientation uw-cube-1); native fields [layer,] 6, N, N in id
  order, a cell table under grid/ (`cells()`) and the lat-lon export of every field on the
  world grid under latlon/ at export_columns W = 4 N.
`load` returns the lat-lon array in every kind (the export on a cube), `load_native` the
stored one; `columns` and `rows` are the lat-lon (export) dimensions.
"""
import csv
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import grid

SCHEMA = 'uw-run-bundle/v2'
SCHEMA_V1 = 'uw-run-bundle/v1'
SCHEMAS = (SCHEMA_V1, SCHEMA)
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


def read_manifest(path):
    """manifest.json of a bundle directory, whatever its schema."""
    return json.loads((Path(path) / 'manifest.json').read_text(encoding='utf-8'))


def grid_kind(manifest):
    """The grid kind of a manifest: v1 bundles are lat-lon."""
    return manifest.get('grid', {}).get('kind', 'latlon')


@dataclass
class Bundle:
    path: Path
    manifest: dict
    fields: dict = field(default_factory=dict)
    tables: dict = field(default_factory=dict)
    exports: dict = field(default_factory=dict)
    cell_entries: dict = field(default_factory=dict)

    @classmethod
    def open(cls, path):
        path = Path(path)
        manifest = read_manifest(path)
        if manifest.get('schema') not in SCHEMAS:
            raise ValueError(f'{path}: unsupported bundle schema {manifest.get("schema")!r}')
        kind = grid_kind(manifest)
        if kind not in ('latlon', 'cube'):
            raise ValueError(f'{path}: unsupported grid kind {kind!r}')
        bundle = cls(path, manifest)
        for entry in manifest['fields']:
            bundle.fields[entry['id']] = FieldInfo(entry['id'], entry['component'], entry['layout'], entry['dtype'],
                                                   tuple(entry['shape']), entry['units'], entry['file'], entry['hash'],
                                                   entry.get('written_by', ''), entry.get('description', ''))
        for entry in manifest['tables']:
            bundle.tables[entry['id']] = entry
        for entry in manifest.get('exports', []):
            bundle.exports[entry['id']] = entry
        for entry in manifest.get('cells', []):
            bundle.cell_entries[entry['id']] = entry
        return bundle

    @property
    def seed(self):
        return self.manifest['seed']

    @property
    def grid_kind(self):
        return grid_kind(self.manifest)

    @property
    def face_cells(self):
        """N cells per face edge on a cube; None on lat-lon."""
        return self.manifest['grid']['face_cells'] if self.grid_kind == 'cube' else None

    @property
    def columns(self):
        """Columns of the lat-lon grid, or of the cube's lat-lon export."""
        grid_block = self.manifest['grid']
        return grid_block['export_columns'] if self.grid_kind == 'cube' else grid_block['columns']

    @property
    def rows(self):
        return self.columns // 2 if self.grid_kind == 'cube' else self.manifest['grid']['rows']

    @property
    def radius_km(self):
        """The radius of the grid block; None when not recorded (v1)."""
        return self.manifest['grid'].get('radius_km')

    @property
    def recipe(self):
        return self.manifest.get('recipe', '')

    def has(self, identifier):
        return identifier in self.fields

    def load(self, identifier):
        """Field values as a lat-lon NumPy array (scalars as 0-d arrays): the stored field on
        lat-lon, its lat-lon export on a cube (stored fields without an export, such as scalars,
        zonal profiles and solver grids, as stored)."""
        if identifier in self.exports:
            return np.load(self.path / self.exports[identifier]['file'])
        return self.load_native(identifier)

    def load_native(self, identifier):
        """The field as stored: [layer,] rows, columns on lat-lon, [layer,] 6, N, N on a cube."""
        info = self.fields[identifier]
        return np.load(self.path / info.file)

    def export_rule(self, identifier):
        """The remap rule of a field's lat-lon export ('mean', 'total', 'nearest', 'direction'), None without one."""
        entry = self.exports.get(identifier)
        return entry['rule'] if entry else None

    def cells(self):
        """The cube's cell table as a dict of arrays in id order: centre (n x 3), latitude_deg,
        longitude_deg, area_m2, neighbours (n x 8, slots 1..8, -1 none), reverse (n x 8, 0 none).
        Lat-lon bundles have none (the table is implicit: uwtools.grid)."""
        if self.grid_kind == 'latlon':
            raise ValueError(f'{self.path}: a lat-lon bundle has no cell table (uwtools.grid gives it)')
        return {identifier: np.load(self.path / entry['file']) for identifier, entry in self.cell_entries.items()}

    def scalar(self, identifier):
        return self.load(identifier).item()

    def table_text(self, identifier):
        return (self.path / self.tables[identifier]['file']).read_text(encoding='utf-8')

    def table_rows(self, identifier):
        return list(csv.DictReader(io.StringIO(self.table_text(identifier))))

    def latitudes(self):
        """Cell-centre latitudes in degrees, north row first."""
        return grid.latitudes(self.rows)

    def longitudes(self):
        return grid.longitudes(self.columns)

    def area_weights(self):
        """Cell area weights (rows x columns), proportional to the cosine of latitude, summing to one."""
        weights = grid.cosine_latitude(self.rows)[:, None] * np.ones((1, self.columns))
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
        weights = grid.cosine_latitude(self.rows)[:, None] * np.ones((1, self.columns))
        total = (weights * sea).reshape(rows, fy, columns, fx).sum(axis=(1, 3))
        return total / weights.reshape(rows, fy, columns, fx).sum(axis=(1, 3))
