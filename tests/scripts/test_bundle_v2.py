"""Bundle v2 in scripts/uwtools/bundle.py, compare_bundles.py and check_grid (docs/architecture/sphere-mesh.md
section 7) on tiny synthetic bundles: v1 and v2 lat-lon (8 x 4) and a cube of N 2 with its export and cell table."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

from uwtools.bundle import SCHEMA, SCHEMA_V1, Bundle  # noqa: E402
from checks import check_grid  # noqa: E402

N = 2
COLUMNS, ROWS = 4 * N, 2 * N
LATLON_GRID = {"columns": COLUMNS, "rows": ROWS, "registration": "cell-centred", "layout": "row-major, north row first"}


def save(directory, relative, values):
    path = Path(directory) / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, values)
    return {"dtype": values.dtype.str, "shape": list(values.shape), "file": relative, "hash": "0"}


def field_entry(directory, identifier, values, layout):
    component, name = identifier.split(".")
    entry = save(directory, f"fields/{component}/{name}.npy", values)
    entry.update({"id": identifier, "component": component, "layout": layout, "units": "", "written_by": "test",
                  "description": ""})
    return entry


def latlon_fields(seed=1):
    random = np.random.default_rng(seed)
    return {"geosphere.base_elevation": (random.normal(0.0, 1000.0, (ROWS, COLUMNS)).astype(np.float32), "grid"),
            "geosphere.sea_level": (np.array(0, dtype=np.int32), "scalar"),
            "hydrosphere.river_direction": (random.integers(0, 9, (ROWS, COLUMNS)).astype(np.int16), "grid")}


def write_latlon(directory, fields, schema):
    """A lat-lon bundle of schema v1 (grid without kind) or v2 (kind latlon, radius)."""
    directory = Path(directory)
    entries = [field_entry(directory, identifier, values, layout) for identifier, (values, layout) in fields.items()]
    grid = dict(LATLON_GRID) if schema == SCHEMA_V1 else {"kind": "latlon", **LATLON_GRID, "radius_km": 6371.0}
    manifest = {"schema": schema, "seed": 1, "recipe": "procedural-world", "grid": grid, "fields": entries, "tables": []}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def write_cube(directory):
    """A cube bundle of N 2: native (6, N, N) and (4, 6, N, N) fields, their exports and a cell table."""
    directory = Path(directory)
    random = np.random.default_rng(1)
    count = 6 * N * N
    base = random.normal(0.0, 1000.0, (6, N, N)).astype(np.float32)
    rain = random.uniform(0.0, 300.0, (4, 6, N, N)).astype(np.float32)
    fields = [field_entry(directory, "geosphere.base_elevation", base, "grid"),
              field_entry(directory, "atmosphere.precipitation", rain, "seasonal"),
              field_entry(directory, "geosphere.sea_level", np.array(0, dtype=np.int32), "scalar")]
    exports = []
    for identifier, values, rule in (("geosphere.base_elevation", np.full((ROWS, COLUMNS), 5.0, np.float32), "mean"),
                                     ("atmosphere.precipitation", np.full((4, ROWS, COLUMNS), 7.0, np.float32), "mean")):
        component, name = identifier.split(".")
        entry = save(directory, f"latlon/{component}/{name}.npy", values)
        entry.update({"id": identifier, "rule": rule})
        exports.append(entry)
    cells = []
    for identifier, values in (("centre", np.zeros((count, 3))), ("latitude_deg", np.zeros(count)),
                               ("longitude_deg", np.zeros(count)), ("area_m2", np.ones(count)),
                               ("neighbours", np.full((count, 8), -1, np.int32)), ("reverse", np.zeros((count, 8), np.int8))):
        entry = save(directory, f"grid/{identifier}.npy", values)
        entry["id"] = identifier
        cells.append(entry)
    grid = {"kind": "cube", "face_cells": N, "orientation": "uw-cube-1", "layout": "face, row (top first), column",
            "radius_km": 6371.0, "export_columns": COLUMNS}
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "procedural-world", "grid": grid, "fields": fields,
                "exports": exports, "cells": cells, "tables": []}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def write_cube_window(directory, keys, level=1, base=N, factor=4, export=(8, 6)):
    """A cube-window bundle of the level keys `keys` (ascending) at N_L = base factor^level: native [n] and
    [4, n] fields whose values are their keys, their exports of the box and the cell table with cell_ids."""
    directory = Path(directory)
    keys = np.asarray(keys, dtype=np.uint64)
    count = keys.size
    rows, columns = export[1], export[0]
    fields = [field_entry(directory, "geosphere.base_elevation", keys.astype(np.float32), "grid"),
              field_entry(directory, "atmosphere.precipitation", np.stack([keys.astype(np.float32) + s for s in range(4)]),
                          "seasonal"),
              field_entry(directory, "geosphere.sea_level", np.array(0, dtype=np.int32), "scalar")]
    exports = []
    for identifier, values in (("geosphere.base_elevation", np.full((rows, columns), 5.0, np.float32)),
                               ("atmosphere.precipitation", np.full((4, rows, columns), 7.0, np.float32))):
        component, name = identifier.split(".")
        entry = save(directory, f"latlon/{component}/{name}.npy", values)
        entry.update({"id": identifier, "rule": "mean"})
        exports.append(entry)
    cells = []
    for identifier, values in (("cell_ids", keys), ("centre", np.zeros((count, 3))), ("area_m2", np.ones(count)),
                               ("neighbours", np.full((count, 8), -1, np.int32)), ("reverse", np.zeros((count, 8), np.int8))):
        entry = save(directory, f"grid/{identifier}.npy", values)
        entry["id"] = identifier
        cells.append(entry)
    grid = {"kind": "cube-window", "face_cells": base * factor ** level, "level": level, "base_face_cells": base,
            "orientation": "uw-cube-1", "layout": "cell list, ascending id", "cell_count": int(count), "radius_km": 6371.0,
            "export_columns": columns, "export_rows": rows}
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "region", "grid": grid, "fields": fields, "exports": exports,
                "cells": cells, "tables": [],
                "window": {"scale": factor ** level, "world_columns": 4 * base, "world_rows": 2 * base, "column0": 0,
                           "columns": 2, "row0": 1, "rows": 2, "fine_column0": 0, "fine_row0": 4, "halo": 2}}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def compare(left, right):
    script = SCRIPTS / "evaluation" / "compare_bundles.py"
    return subprocess.run([sys.executable, str(script), str(left), str(right)], capture_output=True, text=True)


class BundleV2Tests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_v1_and_v2_latlon_read_alike(self):
        fields = latlon_fields()
        for schema in (SCHEMA_V1, SCHEMA):
            bundle = Bundle.open(write_latlon(self.root / schema.replace("/", "-"), fields, schema))
            self.assertEqual(bundle.grid_kind, "latlon")
            self.assertEqual((bundle.columns, bundle.rows), (COLUMNS, ROWS))
            self.assertIsNone(bundle.face_cells)
            np.testing.assert_array_equal(bundle.load("geosphere.base_elevation"), fields["geosphere.base_elevation"][0])
            np.testing.assert_array_equal(bundle.load_native("geosphere.base_elevation"),
                                          fields["geosphere.base_elevation"][0])
            self.assertEqual(bundle.scalar("geosphere.sea_level"), 0)
            with self.assertRaises(ValueError):
                bundle.cells()
        self.assertIsNone(Bundle.open(self.root / SCHEMA_V1.replace("/", "-")).radius_km)
        self.assertEqual(Bundle.open(self.root / SCHEMA.replace("/", "-")).radius_km, 6371.0)

    def test_cube_load_is_the_export_and_load_native_the_faces(self):
        bundle = Bundle.open(write_cube(self.root / "cube"))
        self.assertEqual(bundle.grid_kind, "cube")
        self.assertEqual(bundle.face_cells, N)
        self.assertEqual((bundle.columns, bundle.rows), (COLUMNS, ROWS))
        self.assertEqual(bundle.load("geosphere.base_elevation").shape, (ROWS, COLUMNS))
        self.assertTrue(np.all(bundle.load("geosphere.base_elevation") == 5.0))
        self.assertEqual(bundle.load_native("geosphere.base_elevation").shape, (6, N, N))
        self.assertEqual(bundle.load("atmosphere.precipitation").shape, (4, ROWS, COLUMNS))
        self.assertEqual(bundle.load_native("atmosphere.precipitation").shape, (4, 6, N, N))
        self.assertEqual(bundle.export_rule("geosphere.base_elevation"), "mean")
        self.assertEqual(bundle.scalar("geosphere.sea_level"), 0)   # no export: as stored
        self.assertEqual(bundle.latitudes().shape, (ROWS,))
        cells = bundle.cells()
        self.assertEqual(set(cells), {"centre", "latitude_deg", "longitude_deg", "area_m2", "neighbours", "reverse"})
        self.assertEqual(cells["centre"].shape, (6 * N * N, 3))
        self.assertEqual(cells["neighbours"].shape, (6 * N * N, 8))
        self.assertEqual(cells["neighbours"].dtype, np.int32)
        self.assertEqual(cells["reverse"].dtype, np.int8)

    def test_cube_window_load_is_the_export_and_load_native_the_cells(self):
        keys = np.array([3, 4, 9, 70, 71, 300], dtype=np.uint64)
        bundle = Bundle.open(write_cube_window(self.root / "window", keys))
        self.assertEqual(bundle.grid_kind, "cube-window")
        self.assertEqual((bundle.face_cells, bundle.level, bundle.base_face_cells, bundle.cell_count), (8, 1, N, 6))
        self.assertEqual((bundle.columns, bundle.rows), (8, 6))
        self.assertEqual(bundle.load("geosphere.base_elevation").shape, (6, 8))
        self.assertTrue(np.all(bundle.load("atmosphere.precipitation") == 7.0))
        np.testing.assert_array_equal(bundle.load_native("geosphere.base_elevation"), keys.astype(np.float32))
        self.assertEqual(bundle.load_native("atmosphere.precipitation").shape, (4, 6))
        np.testing.assert_array_equal(bundle.cell_ids(), keys)
        self.assertEqual(bundle.cell_ids().dtype, np.uint64)
        self.assertIn("cell_ids", bundle.cells())
        self.assertEqual(bundle.scalar("geosphere.sea_level"), 0)
        with self.assertRaises(ValueError):
            Bundle.open(write_cube(self.root / "cube")).cell_ids()

    def test_coordinates_describe_the_box_of_a_window(self):
        # A cube window: its export box starts at fine row 4 of the 16-row level grid (scale 4, world 8 x 4).
        window = Bundle.open(write_cube_window(self.root / "window", [1, 2, 3]))
        np.testing.assert_allclose(window.latitudes(), 90.0 - (4 + np.arange(6) + 0.5) * 180.0 / 16)
        np.testing.assert_allclose(window.longitudes(), -180.0 + (np.arange(8) + 0.5) * 360.0 / 32)
        # A lat-lon regional bundle across the date line: fine columns 14 .. 21 of 16 at scale 2.
        regional = write_latlon(self.root / "regional", latlon_fields(), SCHEMA)
        manifest = json.loads((regional / "manifest.json").read_text(encoding="utf-8"))
        manifest["window"] = {"scale": 2, "world_columns": 8, "world_rows": 4, "column0": 7, "columns": 4, "row0": 1,
                              "rows": 2, "fine_column0": 14, "fine_row0": 2, "halo": 0}
        (regional / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        bundle = Bundle.open(regional)
        np.testing.assert_allclose(bundle.latitudes(), 90.0 - (2 + np.arange(ROWS) + 0.5) * 22.5)
        expected = -180.0 + ((14 + np.arange(COLUMNS)) % 16 + 0.5) * 22.5
        np.testing.assert_allclose(bundle.longitudes(), expected)
        self.assertEqual(bundle.longitudes()[2], -180.0 + 0.5 * 22.5)  # past the date line
        weights = bundle.area_weights()
        self.assertEqual(weights.shape, (ROWS, COLUMNS))
        self.assertAlmostEqual(weights.sum(), 1.0)
        np.testing.assert_allclose(weights[:, 0] / weights[0, 0], np.cos(np.deg2rad(bundle.latitudes())) /
                                   np.cos(np.deg2rad(bundle.latitudes()[0])))
        # A world bundle keeps the whole grid.
        world = Bundle.open(write_latlon(self.root / "world", latlon_fields(), SCHEMA))
        np.testing.assert_allclose(world.latitudes(), 90.0 - (np.arange(ROWS) + 0.5) * 45.0)
        self.assertIsNone(world.export_box())

    def test_check_grid_accepts_both_kinds(self):
        latlon = Bundle.open(write_latlon(self.root / "latlon", latlon_fields(), SCHEMA))
        cube = Bundle.open(write_cube(self.root / "cube"))
        self.assertEqual(check_grid(latlon, {}, {})["status"], "pass")
        outcome = check_grid(cube, {}, {})
        self.assertEqual(outcome["status"], "pass", outcome["details"])

    def test_compare_v1_and_v2_latlon_with_equal_fields_is_identical(self):
        fields = latlon_fields()
        left = write_latlon(self.root / "v1", fields, SCHEMA_V1)
        right = write_latlon(self.root / "v2", fields, SCHEMA)
        same = compare(left, right)
        self.assertEqual(same.returncode, 0, same.stdout + same.stderr)
        self.assertIn("3 of 3 fields and tables identical", same.stdout)
        different = compare(left, write_latlon(self.root / "other", latlon_fields(seed=2), SCHEMA))
        self.assertEqual(different.returncode, 1)

    def test_compare_cube_bundles_includes_exports_and_cells(self):
        left, right = write_cube(self.root / "a"), write_cube(self.root / "b")
        same = compare(left, right)
        self.assertEqual(same.returncode, 0, same.stdout + same.stderr)
        self.assertIn("11 of 11 fields, exports and cell arrays and tables identical", same.stdout)
        np.save(right / "latlon/geosphere/base_elevation.npy", np.full((ROWS, COLUMNS), 6.0, np.float32))
        different = compare(left, right)
        self.assertEqual(different.returncode, 1)
        self.assertIn("latlon:geosphere.base_elevation", different.stdout)

    def test_compare_reports_grid_and_one_sided_exports(self):
        cube = write_cube(self.root / "cube")
        latlon = compare(write_latlon(self.root / "latlon", latlon_fields(), SCHEMA), cube)
        self.assertEqual(latlon.returncode, 1)
        self.assertIn("grid.kind  left=latlon, right=cube", latlon.stdout)
        self.assertIn("grid.face_cells", latlon.stdout)
        # A cube bundle without its exports and cells: each entry is missing on that side.
        stripped = write_cube(self.root / "stripped")
        manifest = json.loads((stripped / "manifest.json").read_text(encoding="utf-8"))
        del manifest["exports"], manifest["cells"]
        (stripped / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        missing = compare(cube, stripped)
        self.assertEqual(missing.returncode, 1)
        self.assertIn("missing  latlon:geosphere.base_elevation  side=right", missing.stdout)
        self.assertIn("missing  grid:centre  side=right", missing.stdout)


if __name__ == "__main__":
    unittest.main()
