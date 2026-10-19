"""Read run bundles written by world_pipeline (uw-run-bundle/v2; v1 reads as grid kind latlon).

A bundle is a directory with manifest.json, one NumPy file per field under fields/
and text tables under tables/ (docs/architecture/world-model.md). Bundle v2
(docs/architecture/sphere-mesh.md section 7) has a grid kind:
- latlon: cell-centred W x W/2, row-major with the north row first; fields [layer,] rows, columns.
- cube: N cells per face edge (orientation uw-cube-1); native fields [layer,] 6, N, N in id
  order, a cell table under grid/ (`cells()`) and the lat-lon export of every field on the
  world grid under latlon/ at export_columns W = 4 N.
- cube-window: cell_count level cells of the cube level with face_cells N_L (level L of a world
  of base_face_cells N); native fields [layer,] n in ascending key order, the cell table with the
  u8 level keys grid/cell_ids.npy (`cell_ids()`), and the lat-lon export of the window's box
  under latlon/ (export_columns x export_rows; the window block gives its place).
`load` returns the lat-lon array in every kind (the export on a cube or a cube window),
`load_native` the stored one; `columns` and `rows` are the lat-lon (export) dimensions.
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
GRID_KINDS = ('latlon', 'cube', 'cube-window')


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
        if kind not in GRID_KINDS:
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
        """N cells per face edge on a cube (N_L, the level's, on a cube window); None on lat-lon."""
        return self.manifest['grid']['face_cells'] if self.grid_kind != 'latlon' else None

    @property
    def level(self):
        """A cube window's refinement level; None for the other kinds."""
        return self.manifest['grid']['level'] if self.grid_kind == 'cube-window' else None

    @property
    def base_face_cells(self):
        """A cube window's world N (level 0); None for the other kinds."""
        return self.manifest['grid']['base_face_cells'] if self.grid_kind == 'cube-window' else None

    @property
    def cell_count(self):
        """The cells of a cube window; None for the other kinds."""
        return self.manifest['grid']['cell_count'] if self.grid_kind == 'cube-window' else None

    @property
    def columns(self):
        """Columns of the lat-lon grid, or of the cube's (or cube window's) lat-lon export."""
        grid_block = self.manifest['grid']
        return grid_block['columns'] if self.grid_kind == 'latlon' else grid_block['export_columns']

    @property
    def rows(self):
        kind = self.grid_kind
        if kind == 'cube':
            return self.columns // 2
        return self.manifest['grid']['export_rows' if kind == 'cube-window' else 'rows']

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
        """The field as stored: [layer,] rows, columns on lat-lon, [layer,] 6, N, N on a cube,
        [layer,] n on a cube window (cells in the order of `cell_ids()`)."""
        info = self.fields[identifier]
        return np.load(self.path / info.file)

    def export_rule(self, identifier):
        """The remap rule of a field's lat-lon export ('mean', 'total', 'nearest', 'direction'), None without one."""
        entry = self.exports.get(identifier)
        return entry['rule'] if entry else None

    def cells(self):
        """The cube's cell table as a dict of arrays in id order: centre (n x 3), latitude_deg,
        longitude_deg, area_m2, neighbours (n x 8, slots 1..8, -1 none), reverse (n x 8, 0 none);
        a cube window's also cell_ids (n, u8 level keys) and neighbours as its own local ids.
        Lat-lon bundles have none (the table is implicit: uwtools.grid)."""
        if self.grid_kind == 'latlon':
            raise ValueError(f'{self.path}: a lat-lon bundle has no cell table (uwtools.grid gives it)')
        return {identifier: np.load(self.path / entry['file']) for identifier, entry in self.cell_entries.items()}

    def cell_ids(self):
        """A cube window's level keys f N_L^2 + j N_L + i (u8, ascending), the order of its native fields."""
        if self.grid_kind != 'cube-window':
            raise ValueError(f'{self.path}: only a cube-window bundle lists its cells')
        return np.load(self.path / self.cell_entries['cell_ids']['file'])

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
