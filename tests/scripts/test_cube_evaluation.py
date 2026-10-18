"""Evaluation tools on bundle v2 (W2-EV; docs/architecture/sphere-mesh.md sections 2, 3, 7, 12): v1 and v2 schema
acceptance (tectonics, evaluate, regional_overlap), the cell size without world_width_km, plate_metrics' true
areas on a cube, striping_ratio --native and cube_seams on synthetic cube bundles.

The cube fixture's cell table is built here from the frozen geometry of section 2 (centres, exact areas, point
location; neighbours by locating a point just past each cell edge or corner), test data only: production
Python reads the table C++ writes."""

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_bundle_v2 import field_entry, latlon_fields, save, write_latlon  # noqa: E402
from uwtools import grid, plate_metrics, tectonics  # noqa: E402
from uwtools.bundle import SCHEMA, SCHEMA_V1, Bundle  # noqa: E402
import cube_seams  # noqa: E402
import regional_overlap  # noqa: E402

RADIUS_KM = 6371.0
# uw-cube-1: centre c, right axis u, up axis v of faces 0-5.
FACES = np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
                  [[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                  [[0, 0, 1], [0, 1, 0], [-1, 0, 0]], [[0, 0, -1], [0, 1, 0], [1, 0, 0]]], dtype=np.float64)
# Slots 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW, 7 W, 8 NW as steps in (a, b).
STEPS = ((0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1))


def face_coordinates(n):
    """(face, a, b) of every cell centre in id order."""
    delta = math.pi / (2 * n)
    face, j, i = np.meshgrid(np.arange(6), np.arange(n), np.arange(n), indexing="ij")
    return face.ravel(), -math.pi / 4 + (i.ravel() + 0.5) * delta, math.pi / 4 - (j.ravel() + 0.5) * delta


def point(face, a, b):
    c, u, v = FACES[face, 0], FACES[face, 1], FACES[face, 2]
    p = c + np.tan(a)[:, None] * u + np.tan(b)[:, None] * v
    return p / np.linalg.norm(p, axis=1, keepdims=True)


def locate(p, n):
    delta = math.pi / (2 * n)
    face = np.argmax(p @ FACES[:, 0].T, axis=1)   # the largest coordinate with its sign; ties to the lower face
    c, u, v = FACES[face, 0], FACES[face, 1], FACES[face, 2]
    pc = np.einsum("ij,ij->i", p, c)
    a = np.arctan(np.einsum("ij,ij->i", p, u) / pc)
    b = np.arctan(np.einsum("ij,ij->i", p, v) / pc)
    i = np.clip(np.floor((a + math.pi / 4) / delta), 0, n - 1).astype(np.int64)
    j = np.clip(np.floor((math.pi / 4 - b) / delta), 0, n - 1).astype(np.int64)
    return face * n * n + j * n + i


def cube_cells(n, radius_km=RADIUS_KM):
    """The cell table of a cube of N n (section 3 layout: neighbours n x 8 int32, -1 none; reverse n x 8 int8)."""
    delta = math.pi / (2 * n)
    face, a, b = face_coordinates(n)
    centre = point(face, a, b)
    count = centre.shape[0]
    neighbours = np.full((count, 8), -1, dtype=np.int32)
    for slot, (da, db) in enumerate(STEPS):
        reach = 0.5 + 1.0e-3
        target = locate(point(face, a + da * reach * delta, b + db * reach * delta), n)
        if da and db:   # the diagonal pointing at a cube corner does not exist
            corner = np.isclose(np.abs(a + da * 0.5 * delta), math.pi / 4) & np.isclose(np.abs(b + db * 0.5 * delta),
                                                                                        math.pi / 4)
            target = np.where(corner, -1, target)
        neighbours[:, slot] = target
    reverse = np.zeros((count, 8), dtype=np.int8)
    for cell in range(count):
        for slot in range(8):
            other = neighbours[cell, slot]
            if other >= 0:
                back = np.flatnonzero(neighbours[other] == cell)
                reverse[cell, slot] = back[0] + 1 if back.size else 0

    def f(x, y):
        return np.arctan(x * y / np.sqrt(1 + x * x + y * y))
    x0, x1 = np.tan(a - delta / 2), np.tan(a + delta / 2)
    y0, y1 = np.tan(b - delta / 2), np.tan(b + delta / 2)
    area = (radius_km * 1000.0) ** 2 * (f(x1, y1) - f(x0, y1) - f(x1, y0) + f(x0, y0))
    latitude = np.degrees(np.arcsin(np.clip(centre[:, 2], -1, 1)))
    longitude = np.degrees(np.arctan2(centre[:, 1], centre[:, 0]))
    return {"centre": centre, "latitude_deg": latitude, "longitude_deg": longitude, "area_m2": area,
            "neighbours": neighbours, "reverse": reverse}


def latlon_points(columns):
    rows = columns // 2
    lat = np.radians(grid.latitudes(rows))[:, None]
    lon = np.radians(grid.longitudes(columns))[None, :]
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)], axis=-1)


def isotropic(p, waves=300, wavenumber=20.0, seed=1):
    """A smooth, statistically isotropic field of unit vectors p (..., 3), in metres about 0."""
    random = np.random.default_rng(seed)
    k = random.normal(size=(waves, 3))
    k *= wavenumber / np.linalg.norm(k, axis=1, keepdims=True)
    phase = random.uniform(0, 2 * math.pi, waves)
    return 1000.0 * np.cos(p @ k.T + phase).sum(axis=-1) / math.sqrt(waves / 2)


def geodesic(p, pole=(0.3, -0.5, 0.81)):
    """Metres proportional to the angle from a pole off every cube symmetry: |grad| is the same everywhere
    (1000 m per radian, away from the pole and its antipode), so the field has no seam by construction."""
    pole = np.asarray(pole) / np.linalg.norm(pole)
    return 1000.0 * np.arccos(np.clip(p @ pole, -1.0, 1.0)) - 1500.0


def write_cube(directory, n, fields, tables=()):
    """A cube bundle of N n: fields {id: (native, export or None, layout, units)}, the fixture's cell table."""
    directory = Path(directory)
    entries, exports = [], []
    for identifier, (native, export, layout, units) in fields.items():
        entry = field_entry(directory, identifier, native, layout)
        entry["units"] = units
        entries.append(entry)
        if export is not None:
            component, name = identifier.split(".")
            item = save(directory, f"latlon/{component}/{name}.npy", export)
            item.update({"id": identifier, "rule": "mean"})
            exports.append(item)
    cells = []
    for identifier, values in cube_cells(n).items():
        item = save(directory, f"grid/{identifier}.npy", values)
        item["id"] = identifier
        cells.append(item)
    table_entries = []
    for identifier, relative, document in tables:
        path = directory / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(document), encoding="utf-8")
        table_entries.append({"id": identifier, "file": relative, "hash": "0"})
    grid_block = {"kind": "cube", "face_cells": n, "orientation": "uw-cube-1", "layout": "face, row (top first), column",
                  "radius_km": RADIUS_KM, "export_columns": 4 * n}
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "procedural-world", "grid": grid_block, "fields": entries,
                "exports": exports, "cells": cells, "tables": table_entries}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def elevation_bundle(directory, n, native, export):
    return write_cube(directory, n, {"geosphere.base_elevation": (native.astype(np.float32), export.astype(np.float32),
                                                                  "grid", "m"),
                                     "geosphere.sea_level": (np.array(0, dtype=np.int32), None, "scalar", "m")})


class CubeFixtureTests(unittest.TestCase):
    def test_cell_table_is_consistent(self):
        n = 8
        cells = cube_cells(n)
        self.assertAlmostEqual(cells["area_m2"].sum() / (4 * math.pi * (RADIUS_KM * 1000.0) ** 2), 1.0, places=9)
        neighbours, reverse = cells["neighbours"], cells["reverse"]
        self.assertEqual(int((neighbours < 0).sum()), 24)   # one missing diagonal at each of the 24 corner cells
        for cell in range(neighbours.shape[0]):
            for slot in (0, 2, 4, 6):
                other = neighbours[cell, slot]
                self.assertGreaterEqual(other, 0)
                self.assertEqual(neighbours[other, reverse[cell, slot] - 1], cell)
        # Interior frame: N is j - 1, E is i + 1.
        cell = 0 * n * n + 3 * n + 3
        self.assertEqual(neighbours[cell, 0], cell - n)
        self.assertEqual(neighbours[cell, 2], cell + 1)


class CubeMeasureTests(unittest.TestCase):
    N = 64

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        n = cls.N
        cells = cube_cells(n)
        cls.cells = cells
        smooth = isotropic(cells["centre"]).reshape(6, n, n)
        export = isotropic(latlon_points(4 * n))
        cls.smooth = elevation_bundle(cls.root / "smooth", n, smooth, export)
        cls.geodesic = elevation_bundle(cls.root / "geodesic", n, geodesic(cells["centre"]).reshape(6, n, n),
                                        geodesic(latlon_points(4 * n)))
        face = np.repeat(np.arange(6), n * n).reshape(6, n, n)
        cls.seam = elevation_bundle(cls.root / "seam", n, smooth + 400.0 * face, export)
        _, a, _ = face_coordinates(n)
        aligned = (1000.0 * np.cos(a * 4 * n / 3.0) + 200.0).reshape(6, n, n)   # ridges along each face's v axis
        cls.aligned = elevation_bundle(cls.root / "aligned", n, aligned, export)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_seams_of_a_field_without_seam_are_within_five_percent(self):
        record = cube_seams.seams(Bundle.open(self.geodesic), "geosphere.base_elevation")
        self.assertAlmostEqual(record["ratio_median"], 1.0, delta=0.05, msg=record)
        self.assertAlmostEqual(record["interior"]["median_grad"], 1000.0 / RADIUS_KM, delta=0.01 * 1000.0 / RADIUS_KM)
        self.assertEqual(record["edge"]["cells"] + record["interior"]["cells"], 6 * self.N * self.N)
        # An isotropic random field: no face-aligned aspect at the edges or inside (its gradient ratio between
        # the band and the interior varies by a few per cent with the realisation, so it is not gated here).
        record = cube_seams.seams(Bundle.open(self.smooth), "geosphere.base_elevation")
        self.assertLess(abs(record["interior"]["a4"]), 0.05, record)
        self.assertLess(abs(record["edge"]["a4"]), 0.05, record)

    def test_seams_find_a_jump_at_face_edges(self):
        record = cube_seams.seams(Bundle.open(self.seam), "geosphere.base_elevation")
        self.assertGreater(record["ratio_median"], 1.3, record)

    def test_seams_cli_writes_json(self):
        out = self.root / "seams.json"
        script = SCRIPTS / "evaluation" / "cube_seams.py"
        run = subprocess.run([sys.executable, str(script), str(self.smooth), "--field", "geosphere.base_elevation",
                              "--json", str(out)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        document = json.loads(out.read_text(encoding="utf-8"))
        self.assertEqual(document["grid"]["face_cells"], self.N)
        self.assertIn("geosphere.base_elevation", document["fields"])

    def striping(self, bundle, *extra):
        out = self.root / f"striping-{Path(bundle).name}.json"
        script = SCRIPTS / "sandbox" / "striping_ratio.py"
        return subprocess.run([sys.executable, str(script), "--bundle", str(bundle), "--field",
                               "geosphere.base_elevation", *extra, "--json", str(out)],
                              capture_output=True, text=True), out

    def test_striping_native_aspect_in_face_coordinates(self):
        run, out = self.striping(self.smooth, "--native")
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        record = json.loads(out.read_text(encoding="utf-8"))["bundle"]["fields"]["geosphere.base_elevation"]
        self.assertEqual(record["shape"], [2 * self.N, 4 * self.N])   # the lat-lon export
        self.assertEqual(record["native"]["shape"], [6, self.N, self.N])
        self.assertLess(abs(record["native"]["bands"]["all-0-80"]["a4"]), 0.05, record["native"])
        run, out = self.striping(self.aligned, "--native")
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        record = json.loads(out.read_text(encoding="utf-8"))["bundle"]["fields"]["geosphere.base_elevation"]
        # Below 1 because the frame is orthonormal: near face corners the v grid lines are not perpendicular to u.
        self.assertGreater(record["native"]["bands"]["all-0-80"]["a4"], 0.6, record["native"])

    def test_striping_native_refuses_latlon(self):
        directory = self.root / "latlon-v1"
        if not directory.exists():
            write_latlon(directory, latlon_fields(), SCHEMA_V1)
        run, _ = self.striping(directory, "--native")
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("--native needs a cube bundle", run.stdout + run.stderr)


class PlateMetricsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_cube_areas_from_the_cell_table(self):
        n = 8
        plate = np.repeat(np.arange(6, dtype=np.int32), n * n).reshape(6, n, n)
        snapshot = {"width": 4 * n, "height": 2 * n, "delta_time_myr": 1.0, "time_myr": 10.0,
                    "run_scenario": {"delta_time_myr": 1.0},
                    "plates": [{"plate_id": i, "linear_velocity": [0.1, 0.0], "angular_velocity": 0.0}
                               for i in range(6)]}
        bundle = write_cube(self.root / "cube", n, {
            "geosphere.plate_id": (plate, None, "grid", ""),
            "planetbody.radius": (np.array(RADIUS_KM * 1000.0), None, "scalar", "m")},
            tables=[("geosphere.tectonic_snapshot", "tables/geosphere/tectonic_snapshot.json", snapshot)])
        record = plate_metrics.measure(bundle, 0.05)
        self.assertEqual(record["grid"], {"kind": "cube", "face_cells": n, "export": [4 * n, 2 * n]})
        self.assertEqual(record["plate_count"], 6)
        self.assertAlmostEqual(record["surface_km2"] / (4 * math.pi * RADIUS_KM ** 2), 1.0, places=6)
        for p in record["plates"]:
            self.assertAlmostEqual(p["share"], 1 / 6, places=4)   # the six faces have equal areas
        cell_km = 2 * math.pi * RADIUS_KM / (4 * n)   # no world_width_km: the sphere's equatorial cell
        self.assertAlmostEqual(record["speed"]["km_per_cell"], cell_km)
        self.assertAlmostEqual(record["speed"]["rms_cm_yr_engine"], round(0.1 * cell_km / 10.0, 3))

    def test_latlon_v1_and_v2_agree(self):
        fields = latlon_fields()
        fields["geosphere.plate_id"] = ((np.arange(32).reshape(4, 8) % 3).astype(np.int32), "grid")
        fields["planetbody.radius"] = (np.array(RADIUS_KM * 1000.0), "scalar")
        records = []
        for schema in (SCHEMA_V1, SCHEMA):
            record = plate_metrics.measure(write_latlon(self.root / schema.replace("/", "-"), fields, schema), 0.05)
            record.pop("bundle")
            records.append(record)
        self.assertEqual(records[0], records[1])
        self.assertEqual(records[0]["grid"], [8, 4])

    def test_cell_km(self):
        self.assertEqual(grid.cell_km({"world_width_km": 40000}, 256), 156.25)
        self.assertAlmostEqual(grid.cell_km({}, 256), 2 * math.pi * 6371.0 / 256)
        self.assertAlmostEqual(grid.cell_km(None, 256, 3000.0), 2 * math.pi * 3000.0 / 256)


class SchemaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_v1_and_v2_run_bundles_are_run_bundles(self):
        fields = {"geosphere.plate_id": (np.zeros((4, 8), np.int32), "grid")}
        for schema in (SCHEMA_V1, SCHEMA):
            directory = write_latlon(self.root / schema.replace("/", "-"), fields, schema)
            self.assertTrue(tectonics.is_run_bundle(directory), schema)
            self.assertEqual(tectonics.bundle_files(directory), {"plate_id": "fields/geosphere/plate_id.npy"})
        former = self.root / "former"
        former.mkdir()
        (former / "manifest.json").write_text(json.dumps({"schema": "tectonic-snapshot-bundle/v1",
                                                          "files": {"plate_id": "plate_id.bin"}}))
        self.assertFalse(tectonics.is_run_bundle(former))
        self.assertEqual(tectonics.bundle_files(former), {"plate_id": "plate_id.bin"})

    def test_evaluate_finds_checkpoints_of_both_schemas(self):
        import evaluate
        checkpoints = self.root / "checkpoints"
        write_latlon(checkpoints / "a", latlon_fields(), SCHEMA_V1)
        write_latlon(checkpoints / "b", latlon_fields(), SCHEMA)
        self.assertEqual([b.path.name for b in evaluate.checkpoint_bundles(checkpoints)], ["a", "b"])

    def test_regional_overlap_reads_both_schemas_and_refuses_a_cube(self):
        for schema in (SCHEMA_V1, SCHEMA):
            directory = write_latlon(self.root / schema.replace("/", "-"), latlon_fields(), schema)
            manifest = json.loads((directory / "manifest.json").read_text())
            manifest["window"] = {"scale": 1}
            (directory / "manifest.json").write_text(json.dumps(manifest))
            self.assertEqual(regional_overlap.window_of(directory)[1], {"scale": 1})
        cube = write_cube(self.root / "cube", 2, {"geosphere.sea_level": (np.array(0, np.int32), None, "scalar", "")})
        with self.assertRaises(regional_overlap.Refused):
            regional_overlap.window_of(cube)


if __name__ == "__main__":
    unittest.main()
