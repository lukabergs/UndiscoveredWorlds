"""The terrain sandbox on the sphere (W2-SB1; docs/architecture/sphere-mesh.md sections 7, 11, 12):
- overlaps: two lat-lon windows of the same scale that share cells compute every pointwise layer of the procedural
  pass (the base, the warp, the strike, the relief noise, the rock, the provinces, the rain, the fill noise) with the
  same bits there, now that the noise is evaluated at the cells' unit vectors (a synthetic world; base-v9 too when
  the main tree holds it);
- the cube grid of a synthetic cube base (`terrain_sandbox synthetic-base`): Python's point location agrees with the
  C++ table, the bilinear stencil returns the cells' values at their centres, the smoothing keeps constants, the
  labels, distances and gradients behave, the noise shows no face imprint;
- one global run on the cube (cube.enabled) of the synthetic base: finite native and exported elevations, the maps.
The cube tests need terrain_sandbox.exe (out/build/headless) and are skipped without it."""

import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "sandbox"))

import cube_domain  # noqa: E402
import terrain_sandbox as ts  # noqa: E402
from uwtools.bundle import Bundle  # noqa: E402

EXE = ts.PHYSICS
BASE_V9 = Path("C:/dev/undiscovered-worlds/out/scratch/sandbox/base-v9")
RADIUS_KM = 6371.0
POINTWISE = ("base", "relief_noise", "hardness", "precipitation")


def synthetic_latlon_world(columns=32):
    """World fields on a W x W/2 lat-lon grid from smooth functions of position (seed 1): a continent, a belt, oceans."""
    rows = columns // 2
    latitude = np.radians(90.0 - (np.arange(rows) + 0.5) * 180.0 / rows)[:, None]
    longitude = np.radians(-180.0 + (np.arange(columns) + 0.5) * 360.0 / columns)[None, :]
    land = np.cos(latitude - 0.3) * np.cos(longitude - 0.5) - 0.45
    belt = np.exp(-((longitude - 0.6) / 0.15) ** 2) * (land > 0)
    elevation = np.where(land > 0, 300.0 + 2000.0 * land + 2500.0 * belt, -1500.0 + 6000.0 * np.minimum(land, 0.0))
    shape = elevation.shape
    world = {"elevation": elevation, "sea": elevation <= 0.0, "shortening": 0.4 * belt,
             "direction": np.zeros(shape), "uplift": 200.0 * belt, "sediment": np.where(land > 0, 300.0, 500.0),
             "continental_crust": np.where(elevation > -1000.0, 35000.0 + 15000.0 * belt, 0.0),
             "crust_class": np.where(elevation > -1000.0, 2, 1).astype(np.uint8),
             "crust_age": np.minimum(np.abs(np.degrees(longitude + 2.0)) * 2.0, 180.0) * np.ones(shape),
             "boundary_type": np.where(belt > 0.8, 1, 0).astype(np.uint8), "accumulated_strain": 0.3 * belt,
             "deforming_region_type": np.where(belt > 0.5, 2, 0).astype(np.uint8), "grid": None}
    return world


def prepare(world, config, radius_km=RADIUS_KM):
    """The world-scale fields main() adds (rain, abyss, sea floor, provinces)."""
    world["precipitation"] = ts.precipitation(world, config["precipitation"], radius_km * 1000.0)
    if config["seafloor"].get("abyss", {}).get("enabled", False):
        world["elevation"] = ts.damp_abyss(world, config["seafloor"]["abyss"])
    world["ocean"] = ts.ocean_world(world, config["seafloor"], config["seed"], radius_km)
    columns = world["elevation"].shape[-1]
    world["provinces"] = ts.provinces_world(world, config["platform"], 2.0 * math.pi * radius_km / columns)
    return world


def layers(fields, domain, world, config):
    """The pointwise layers of a procedural pass, by name, and the sea floor's and the coast's position functions:
    the segments (shifted positions and trough), the seamounts, the hills bank on the shifted positions, the coast
    term's displacement noise and the micro texture's."""
    stages = fields["stages"]
    out = {name: stages[name] for name in POINTWISE}
    out.update({f"warp{axis}": stages["warp"][axis] for axis in range(3)}, theta=stages["strike"][0],
               alignment=stages["strike"][1], amplitude=fields["amplitude"], stability=fields["stability"],
               province=fields["province"], fill_noise=fields["fill_noise"])
    seed = config["seed"]
    spec = ts.seafloor_spec(config)
    sampler = ts.make_sampler(domain, world)
    points = domain.points()
    oceanic = np.clip(sampler(world["ocean"]["oceanic"]), 0.0, 1.0)
    doubled = sampler(world["ocean"]["doubled"].real) + 1j * sampler(world["ocean"]["doubled"].imag)
    theta = 0.5 * np.angle(doubled)
    shifted, trough = ts.segment_shift(domain, points, theta, np.ones_like(oceanic), spec, seed)
    out.update({f"segments{axis}": shifted[axis] for axis in range(3)}, trough=trough,
               seamounts=ts.seamounts(points, spec["seamounts"], seed, domain.radius_km),
               hills=ts.oriented_fbm(shifted, theta, np.ones_like(oceanic), {**spec["hills"], **ts.bank_options(spec)},
                                     seed + 121, domain))
    coast = config["coast"]
    warp_spec = {**coast, "wavelength_km": coast["warp_wavelength_km"], "octaves": coast["warp_octaves"]}
    out["coast_displacement"] = ts.fbm(fields["coords"], warp_spec, seed + 71, domain)
    micro = config["micro"]
    moved = fields["plain"] + micro["warp_fraction"] * (fields["warped"] - fields["plain"])
    out["micro_displacement"] = ts.fbm(moved, micro, seed + 75, domain)
    return out


class OverlapTest(unittest.TestCase):
    """Two windows of one scale sharing cells: the pointwise layers are bit-identical on the shared cells."""

    def check_overlap(self, config, world, bundle_columns, bundle_rows, first, second, scale):
        config = json.loads(json.dumps(config))
        config["cube"]["enabled"] = False   # lat-lon worlds (the default base is the cube base-v10)
        config["domains"] = {"A": {"pixels": first, "scale": scale, "margin": 1},
                             "B": {"pixels": second, "scale": scale, "margin": 1}}
        config["physics"]["use_levels"] = False

        class Grid:  # what make_domains reads of a bundle
            columns, rows, manifest, grid_kind = bundle_columns, bundle_rows, {}, "latlon"
        domains = ts.make_domains(config, Grid, ["A", "B"])
        results = [layers(ts.procedural(domain, world, config, config["seed"])[0], domain, world, config)
                   for domain in domains]
        a, b = domains
        trim = 2  # cells away from either window's edge
        r0, r1 = max(a.row0, b.row0) + trim, min(a.row0 + a.rows, b.row0 + b.rows) - trim
        c0, c1 = max(a.column0, b.column0) + trim, min(a.column0 + a.columns, b.column0 + b.columns) - trim
        self.assertGreater(r1 - r0, 0)
        self.assertGreater(c1 - c0, 0)
        for name in results[0]:
            first_cells = results[0][name][r0 - a.row0:r1 - a.row0, c0 - a.column0:c1 - a.column0]
            second_cells = results[1][name][r0 - b.row0:r1 - b.row0, c0 - b.column0:c1 - b.column0]
            self.assertTrue(np.array_equal(first_cells, second_cells), f"{name} differs on the shared cells")
            if name in ("relief_noise", "warp0", "hardness", "segments0", "hills", "coast_displacement"):
                self.assertGreater(float(np.std(first_cells)), 0.0, name)   # not a constant overlap

    def test_synthetic_world(self):
        config = ts.load_config(ts.DEFAULT_CONFIG, [])
        world = prepare(synthetic_latlon_world(128), config)
        # Scale 64 on 128 columns: 4.9 km cells, so the sea-floor hills (20 km) are drawn too (none fade out).
        self.check_overlap(config, world, 128, 64, [40, 16, 41, 17], [41, 17, 42, 18], 64)

    @unittest.skipUnless((BASE_V9 / "manifest.json").is_file(), "base-v9 is not in the main tree")
    def test_base_v9(self):
        config = ts.load_config(ts.DEFAULT_CONFIG, [f"base_bundle={BASE_V9.as_posix()}"])
        bundle = Bundle.open(BASE_V9)
        world = prepare(ts.world_fields(bundle), config)
        self.check_overlap(config, world, bundle.columns, bundle.rows, [160, 64, 161, 65], [161, 65, 162, 66], 64)


@unittest.skipUnless(EXE.is_file(), f"{EXE} is not built")
class CubeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.base = cls.root / "synthetic-n8"
        run = subprocess.run([str(EXE), "synthetic-base", "--face-cells", "8", "--out", str(cls.base)],
                             capture_output=True, text=True)
        assert run.returncode == 0, run.stderr
        cls.bundle = Bundle.open(cls.base)
        cls.grid = cube_domain.CubeGrid.from_bundle(cls.bundle)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_bundle_is_v2_cube(self):
        self.assertEqual(self.bundle.grid_kind, "cube")
        self.assertEqual(self.bundle.load_native("geosphere.base_elevation").shape, (6, 8, 8))
        self.assertEqual(self.bundle.load("geosphere.base_elevation").shape, (16, 32))
        self.assertEqual(self.bundle.cells()["neighbours"].shape, (384, 8))

    def test_location_matches_the_table(self):
        ids, ta, tb = self.grid.locate(self.grid.centre.T)
        self.assertTrue(np.array_equal(ids, np.arange(self.grid.n)))
        self.assertLess(np.abs(ta).max(), 1e-9)
        self.assertLess(np.abs(tb).max(), 1e-9)

    def test_bilinear_and_smooth(self):
        values = np.random.default_rng(1).normal(size=self.grid.n)
        self.assertLess(np.abs(self.grid.bilinear(values, self.grid.centre.T) - values).max(), 1e-12)
        constant = np.full(self.grid.n, 3.5)
        self.assertLess(np.abs(self.grid.smooth(constant, 2.0) - 3.5).max(), 1e-12)
        points = np.random.default_rng(2).normal(size=(3, 1000))
        sampled = self.grid.bilinear(constant, points)
        self.assertLess(np.abs(sampled - 3.5).max(), 1e-12)   # corners included: the weights sum to 1

    def test_labels_distances_gradient(self):
        labels, count = self.grid.label(np.ones(self.grid.n, bool))
        self.assertEqual(count, 1)
        face = np.arange(self.grid.n) // 64
        labels, count = self.grid.label((face == 0) | (face == 2))
        self.assertEqual(count, 2)   # faces 0 and 2 share no edge
        source = np.zeros(self.grid.n, bool)
        source[0] = True
        distance = self.grid.distance_km(source)
        neighbour = self.grid.neighbours[0, 2]
        self.assertAlmostEqual(distance[neighbour], self.grid.slot_km[0, 2], places=4)
        self.assertTrue(np.isfinite(distance).all())
        # The great-circle distance to the nearest source, exact for one source.
        exact = self.grid.great_circle_km(np.arange(self.grid.n), np.zeros(self.grid.n, dtype=np.int64))
        self.assertLess(np.abs(distance - exact).max(), 1e-6)
        bounded = self.grid.distance_km(source, cutoff_km=1500.0)
        self.assertTrue(np.array_equal(np.isfinite(bounded), exact <= 1500.0))
        self.assertLess(np.abs(bounded - exact)[np.isfinite(bounded)].max(), 1e-6)
        two = source.copy()
        two[200] = True
        nearest = np.minimum(exact, self.grid.great_circle_km(np.arange(self.grid.n), np.full(self.grid.n, 200)))
        self.assertLess(np.abs(self.grid.distance_km(two) - nearest).max(), 1e-6)
        east, north = self.grid.gradient(RADIUS_KM * self.grid.centre[:, 2])   # z: d/dnorth = cos(latitude)
        interior = np.abs(np.degrees(self.grid.latitude)) < 60.0
        self.assertLess(np.abs(north - np.cos(self.grid.latitude))[interior].max(), 0.05)
        self.assertLess(np.abs(east)[interior].max(), 0.05)

    def test_wide_smoothing(self):
        """smooth_wide (dilated 3 x 3 means along straight walks) keeps a constant (a delta's total within 3 %) and
        spreads a delta as the 3 x 3 passes of the same variance do; the walks cross face edges by the S rule."""
        table = cube_domain.cube_table(EXE, 32, RADIUS_KM, self.root)
        constant = np.full(table.n, 2.5)
        self.assertLess(np.abs(table.smooth_wide(constant, 20.0) - 2.5).max(), 1e-12)
        spacing = 2.0 * math.pi * RADIUS_KM / 128.0
        for cell in (0, 32 * 16 + 16):   # a face corner cell and a face centre
            delta = np.zeros(table.n)
            delta[cell] = 1.0
            d2 = (np.linalg.norm(table.centre - table.centre[cell], axis=1) * RADIUS_KM / spacing) ** 2
            wide = table.smooth_wide(delta, 20.0)
            passes = table.smooth(delta, math.sqrt(20.0))
            self.assertLess(abs(wide.sum() - 1.0), 0.03, cell)   # a gather mean keeps constants, sums only about
            self.assertLess(abs((wide * d2).sum() / (passes * d2).sum() - 1.0), 0.15, cell)

    def test_noise_has_no_face_imprint(self):
        table = cube_domain.cube_table(EXE, 64, RADIUS_KM, self.root)
        domain = cube_domain.CubeDomain("noise", table, 8, RADIUS_KM)
        spec = {"wavelength_km": 600, "octaves": 4, "lacunarity": 2.0, "gain": 0.5, "min_cells": 2.5}
        noise = ts.fbm(domain.points(), spec, 1, domain)
        column = np.arange(table.n) % 64
        row = (np.arange(table.n) // 64) % 64
        edge = np.minimum.reduce([column, row, 63 - column, 63 - row]) < 2
        self.assertLess(abs(noise[edge].std() / noise[~edge].std() - 1.0), 0.1)
        self.assertLess(abs(noise[edge].mean() - noise[~edge].mean()), 0.1)

    def test_physics_checks_the_table(self):
        n = self.grid.n
        for name, values in (("elevation", self.bundle.load_native("geosphere.base_elevation").reshape(-1)),
                             ("sea", (self.bundle.load_native("geosphere.base_elevation") <= 0).reshape(-1))):
            np.asarray(values, dtype="<f4").tofile(self.root / f"{name}.f32")
        wrong = self.grid.neighbours.astype("<i4").copy()
        wrong[5, 2] = 0
        wrong.tofile(self.root / "wrong.i32")
        command = [str(EXE), "--face-cells", "8", "--elevation", str(self.root / "elevation.f32"),
                   "--sea", str(self.root / "sea.f32"), "--out-elevation", str(self.root / "out.f32"),
                   "--out-discharge", str(self.root / "discharge.f32"), "--duration", "1", "--step", "1"]
        run = subprocess.run(command + ["--cells", str(self.root / "wrong.i32")], capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("neighbour table differs", run.stderr)
        self.grid.neighbours.astype("<i4").tofile(self.root / "right.i32")
        run = subprocess.run(command + ["--cells", str(self.root / "right.i32")], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        out = np.fromfile(self.root / "out.f32", dtype="<f4")
        self.assertEqual(out.size, n)
        self.assertTrue(np.isfinite(out).all())

    def test_global_run_on_the_cube(self):
        label = self.root / "run"
        run = subprocess.run([sys.executable, str(SCRIPTS / "sandbox" / "terrain_sandbox.py"), "--label", str(label),
                              "--domains", "global", "--set", f"base_bundle={self.base.as_posix()}",
                              "--set", "cube.enabled=true", "--set", "domains.global.columns=64"],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout[-2000:] + run.stderr[-4000:])
        native = np.load(label / "global-elevation-native.npy")
        exported = np.load(label / "global-elevation.npy")
        self.assertEqual(native.shape, (6, 16, 16))
        self.assertEqual(exported.shape, (32, 64))
        self.assertTrue(np.isfinite(native).all() and np.isfinite(exported).all())
        self.assertTrue((native > 0).any() and (native < 0).any())
        for name in ("index.html", "global-relief.png", "global-stage1-base.png", "global-diag-strike.png"):
            self.assertTrue((label / name).is_file(), name)
        with self.assertRaises(ValueError):   # nesting (rule B) has no cube form
            config = ts.load_config(ts.DEFAULT_CONFIG, ["cube.enabled=true", "nesting.enabled=true"])
            ts.make_domains(config, self.bundle, ["B"], self.root)

    def cube_world(self, config):
        world = ts.world_fields(self.bundle, native=True)
        world["precipitation"] = ts.precipitation(world, config["precipitation"], RADIUS_KM * 1000.0)
        world["elevation"] = ts.damp_abyss(world, config["seafloor"]["abyss"])
        world["ocean"] = ts.ocean_world(world, ts.seafloor_spec(config), config["seed"], RADIUS_KM, self.root)
        world["provinces"] = ts.provinces_world(world, config["platform"], 2.0 * math.pi * RADIUS_KM / self.bundle.columns)
        return world

    def test_cube_window_overlap(self):
        """Two cube windows sharing cells, one across the face 0 / 1 edge (longitude 45, pixel column 20 of 32):
        every pointwise layer bit-identical on the shared level keys."""
        config = ts.load_config(ts.DEFAULT_CONFIG, ["cube.enabled=true", "physics.use_levels=false"])
        config["domains"] = {"T1": {"pixels": [18, 6, 20, 8], "scale": 64, "margin": 1},
                             "T2": {"pixels": [19, 7, 21, 9], "scale": 64, "margin": 1}}
        world = self.cube_world(config)
        domains = ts.make_domains(config, self.bundle, ["T1", "T2"], self.root)
        faces = [np.unique(domain.grid.keys // np.uint64(domain.grid.N * domain.grid.N)) for domain in domains]
        self.assertTrue(set(faces[0].tolist()) >= {0, 1}, faces[0])   # the first window crosses the face edge
        results = [layers(ts.procedural(domain, world, config, config["seed"])[0], domain, world, config)
                   for domain in domains]
        # Away from either window's edge by more than the widest filter (the trench crest's maximum over k / 2 rings).
        interior = []
        for domain in domains:
            edge = ~domain.grid.valid[:, list(cube_domain.CARDINAL)].all(1)
            interior.append(~domain.grid.dilate(edge, domain.k // 2 + 2))
        shared, first, second = np.intersect1d(domains[0].grid.keys, domains[1].grid.keys, return_indices=True)
        keep = interior[0][first] & interior[1][second]
        shared, first, second = shared[keep], first[keep], second[keep]
        self.assertGreater(shared.size, 1000)
        shared_faces = np.unique(shared // np.uint64(domains[0].grid.N * domains[0].grid.N))
        self.assertTrue(set(shared_faces.tolist()) >= {0, 1}, shared_faces)   # shared cells on both faces
        for name in results[0]:
            a, b = np.asarray(results[0][name])[..., first], np.asarray(results[1][name])[..., second]
            self.assertTrue(np.array_equal(a, b), f"{name} differs on the shared cells")
            if name in ("relief_noise", "warp0", "hardness", "base"):
                self.assertGreater(float(np.std(a)), 0.0, name)

    def test_window_run_on_the_cube(self):
        label = self.root / "window-run"
        run = subprocess.run([sys.executable, str(SCRIPTS / "sandbox" / "terrain_sandbox.py"), "--label", str(label),
                              "--domains", "W", "--set", f"base_bundle={self.base.as_posix()}",
                              "--set", "cube.enabled=true",
                              "--set", 'domains.W={"pixels": [18, 6, 20, 8], "scale": 16, "margin": 1}'],
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout[-2000:] + run.stderr[-4000:])
        exported = np.load(label / "W-elevation.npy")
        self.assertEqual(exported.shape, (3 * 16, 3 * 16))   # the box at k pixels per pixel
        native = np.load(label / "W-elevation-native.npy")
        self.assertEqual(native.shape, np.load(label / "W-keys.npy").shape)
        self.assertTrue(np.isfinite(native).all() and np.isfinite(exported).all())
        self.assertTrue((label / "W-relief.png").is_file())


if __name__ == "__main__":
    unittest.main()
