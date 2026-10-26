"""The sandbox's sloped fill (terrain_sandbox.sloped_fill, CubeGrid.sloped_fill; master plan 5.3):
- metric "chamfer" (the default) gives exactly the output of the function as it was before the metric existed (frozen
  copies below), on a lat-lon window and on a cube, with noise;
- metric "euclid" raises a basin drained by one spill cell to g x the great-circle distance to the spill (a round cone,
  not the chamfer's octagon with up to 8.2 % long diagonals), on a lat-lon window and across a cube face edge;
- with either metric every raised cell has a lower neighbour (the fill leaves no closed depression);
- the exe's threshold-slope limiter: --slope-limit-neighbours 8 is the default, 16 cuts a rounder cone.
The cube and limiter tests need terrain_sandbox.exe (out/build/headless) and are skipped without it.
Run: python -m unittest tests/scripts/test_sloped_fill.py"""

import heapq
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "sandbox"))

import cube_domain  # noqa: E402
import terrain_sandbox as ts  # noqa: E402

RADIUS_KM = 6371.0
GRADIENT = 1e-3
# euclid: |fill - g D| at most this share of g D (the 0.1 g d floor binds nowhere on these cones).
EUCLID_TOLERANCE = 1e-9


def legacy_latlon(domain, elevation, sea, gradient, noise=None):
    """terrain_sandbox.sloped_fill's lat-lon path at 18d0233 (before the metric)."""
    rows, columns = elevation.shape
    wrap = domain.whole
    pad = columns // 8 if wrap else 0
    padding = ((1, 1), (pad + 1, pad + 1))
    mode = 'wrap' if wrap else 'edge'
    z = np.pad(np.asarray(elevation, dtype=float), padding, mode=mode)
    water = np.pad(sea, padding, mode=mode)
    n = np.pad(np.zeros(elevation.shape) if noise is None else np.asarray(noise, dtype=float), padding, mode=mode)
    done = water.copy()
    done[0, :] = done[-1, :] = done[:, 0] = done[:, -1] = True
    width = z.shape[1]
    dy = domain.fine_km * 1000.0
    dx = dy * np.pad(domain.coslat, (1, 1), mode='edge')
    offsets = (-width - 1, -width, -width + 1, -1, 1, width - 1, width, width + 1)
    rises = [[gradient * math.hypot(dx[row] if dc else 0.0, dy if dr else 0.0)
              for dr, dc in ((1, 1), (1, 0), (1, 1), (0, 1), (0, 1), (1, 1), (1, 0), (1, 1))] for row in range(z.shape[0])]
    shore = water & ndimage.binary_dilation(~water, structure=np.ones((3, 3)))
    shore[0, :] = shore[-1, :] = shore[:, 0] = shore[:, -1] = False
    height = z.ravel().tolist()
    lift = n.ravel().tolist()
    filled = list(height)
    closed = done.ravel().tolist()
    heap = [(height[i], i) for i in np.flatnonzero(shore.ravel()).tolist()]
    heapq.heapify(heap)
    while heap:
        level, cell = heapq.heappop(heap)
        run = rises[cell // width]
        for k, offset in enumerate(offsets):
            j = cell + offset
            if closed[j]:
                continue
            closed[j] = True
            if height[j] <= level:
                step = run[k] + lift[j] - lift[cell]
                filled[j] = level + max(0.1 * run[k], step)
                heapq.heappush(heap, (filled[j], j))
            else:
                heapq.heappush(heap, (height[j], j))
    filled = np.asarray(filled).reshape(z.shape)[1:-1, pad + 1:pad + 1 + columns]
    return np.where(sea, elevation, filled)


def legacy_cube(grid, elevation, sea, gradient, noise=None):
    """CubeGrid.sloped_fill at 18d0233 (before the metric)."""
    elevation = np.asarray(elevation, dtype=np.float64)
    sea = np.asarray(sea, dtype=bool)
    neighbours = grid.neighbours.tolist()
    rises = (gradient * np.where(grid.valid, grid.slot_km, 0.0) * 1000.0).tolist()
    height = elevation.tolist()
    lift = (np.zeros(grid.n) if noise is None else np.asarray(noise, dtype=np.float64)).tolist()
    filled = list(height)
    closed = sea.tolist()
    shore = sea & ~grid.minimum(sea, 1)
    heap = [(height[i], i) for i in np.flatnonzero(shore).tolist()]
    heapq.heapify(heap)
    while heap:
        level, cell = heapq.heappop(heap)
        run = rises[cell]
        for k, j in enumerate(neighbours[cell]):
            if j < 0 or closed[j]:
                continue
            closed[j] = True
            if height[j] <= level:
                filled[j] = level + max(0.1 * run[k], run[k] + lift[j] - lift[cell])
                heapq.heappush(heap, (filled[j], j))
            else:
                heapq.heappush(heap, (height[j], j))
    return np.where(sea, elevation, np.asarray(filled))


def window(size=41):
    """A size x size window at k 8 on a 256 x 128 world, centred on the equator (cells about 19.5 km)."""
    k, half = 8, size // 2
    x0, y0 = (128 * k - half) / k, (64 * k - half) / k
    return ts.Domain('fill-test', size, size, x0, y0, k, 256, 128, RADIUS_KM, crop=(0, size, 0, size))


def great_circle_m(lat1, lon1, lat2, lon2):
    """Haversine distance (m) between points in radians."""
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2.0 * RADIUS_KM * 1000.0 * np.arcsin(np.sqrt(np.minimum(a, 1.0)))


def rough_terrain(shape, seed=1):
    """Random land with pits (seed 1): smoothed white noise, a sea margin on the west."""
    rng = np.random.default_rng(seed)
    elevation = 200.0 + 150.0 * ndimage.gaussian_filter(rng.standard_normal(shape), 2.0) * 4.0
    noise = 5.0 * ndimage.gaussian_filter(rng.standard_normal(shape), 1.5) * 3.0
    return elevation, noise


class LatLonTest(unittest.TestCase):
    def setUp(self):
        self.domain = window()
        self.shape = (self.domain.rows, self.domain.columns)

    def cone(self, metric):
        """A flat basin at -100 m drained by one sea cell (0 m) in the window's centre."""
        elevation = np.full(self.shape, -100.0)
        sea = np.zeros(self.shape, dtype=bool)
        centre = self.shape[0] // 2
        sea[centre, centre] = True
        elevation[centre, centre] = 0.0
        filled = ts.sloped_fill(self.domain, elevation, sea, GRADIENT, metric=metric)
        latitude = np.radians(self.domain.latitude_deg)[:, None] * np.ones(self.shape)
        longitude = np.radians(self.domain.longitude_deg)[None, :] * np.ones(self.shape)
        distance = great_circle_m(latitude, longitude, latitude[centre, centre], longitude[centre, centre])
        # The window's frame is closed (never filled): compare the interior.
        inner = np.zeros(self.shape, dtype=bool)
        inner[1:-1, 1:-1] = True
        inner &= ~sea
        return filled, GRADIENT * distance, inner

    def test_chamfer_is_the_legacy_function(self):
        elevation, noise = rough_terrain(self.shape)
        sea = np.zeros(self.shape, dtype=bool)
        sea[:, :3] = True
        elevation[sea] = -50.0
        legacy = legacy_latlon(self.domain, elevation, sea, 2e-5, noise)
        self.assertTrue((legacy - elevation > 1.0).any())  # the terrain has pits that fill
        for kwargs in ({}, {'metric': 'chamfer'}):
            np.testing.assert_array_equal(ts.sloped_fill(self.domain, elevation, sea, 2e-5, noise, **kwargs), legacy)
        whole = ts.Domain('fill-whole', 64, 32, 0.0, 0.0, 1, 64, 32, RADIUS_KM)
        elevation, noise = rough_terrain((32, 64), seed=2)
        sea = elevation < 150.0
        np.testing.assert_array_equal(ts.sloped_fill(whole, elevation, sea, 2e-5, noise),
                                      legacy_latlon(whole, elevation, sea, 2e-5, noise))

    def test_euclid_is_a_round_cone(self):
        filled, expected, inner = self.cone('euclid')
        error = np.abs(filled - expected)[inner] / expected[inner]
        self.assertLessEqual(float(error.max()), EUCLID_TOLERANCE)
        # The chamfer cone is octagonal: its diagonals' rise is up to 8.2 % above g D (sqrt(4 - 2 sqrt 2) - 1).
        chamfer, _, _ = self.cone('chamfer')
        excess = (chamfer - expected)[inner] / expected[inner]
        self.assertGreater(float(excess.max()), 0.07)
        self.assertLess(float(excess.max()), 0.0825)

    def test_every_raised_cell_drains(self):
        elevation, noise = rough_terrain(self.shape, seed=3)
        sea = np.zeros(self.shape, dtype=bool)
        sea[:, :3] = True
        for metric in ('chamfer', 'euclid'):
            filled = ts.sloped_fill(self.domain, elevation, sea, 2e-5, noise, metric=metric)
            raised = (filled > elevation) & ~sea
            raised[[0, -1], :] = raised[:, [0, -1]] = False
            lowest = ndimage.minimum_filter(np.where(sea, -np.inf, filled), size=3, mode='nearest')
            self.assertTrue((lowest[raised] < filled[raised]).all(), metric)

    def test_unknown_metric(self):
        with self.assertRaises(ValueError):
            ts.sloped_fill(self.domain, np.zeros(self.shape), np.zeros(self.shape, dtype=bool), 2e-5, metric='taxicab')


@unittest.skipUnless(ts.PHYSICS.is_file(), 'needs terrain_sandbox.exe for the cube table')
class CubeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.grid = cube_domain.cube_table(ts.PHYSICS, 16, RADIUS_KM, Path(cls.temporary.name))

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_chamfer_is_the_legacy_function(self):
        grid = self.grid
        rng = np.random.default_rng(1)
        elevation = 200.0 + 300.0 * grid.smooth(rng.standard_normal(grid.n), 1.5)
        noise = 5.0 * grid.smooth(rng.standard_normal(grid.n), 1.0)
        sea = elevation < 150.0
        legacy = legacy_cube(grid, elevation, sea, 2e-5, noise)
        self.assertTrue((legacy - elevation > 1.0).any())
        for kwargs in ({}, {'metric': 'chamfer'}):
            np.testing.assert_array_equal(grid.sloped_fill(elevation, sea, 2e-5, noise, **kwargs), legacy)

    def test_euclid_is_a_round_cone_across_a_face_edge(self):
        grid, N = self.grid, self.grid.N
        spill = 0 * N * N + (N // 2) * N + (N - 1)  # face 0, the cell beside its east edge
        sea = np.zeros(grid.n, dtype=bool)
        sea[spill] = True
        elevation = np.where(sea, 0.0, -100.0)
        filled = grid.sloped_fill(elevation, sea, GRADIENT, metric='euclid')
        distance = grid.great_circle_km(np.arange(grid.n), np.full(grid.n, spill)) * 1000.0
        near = ~sea & (distance < 6 * 2 * math.pi * RADIUS_KM * 1000.0 / (4 * N))  # within about 6 cells
        across = near & (np.arange(grid.n) // (N * N) != 0)
        self.assertGreater(int(across.sum()), 20)  # the cone reaches over the edge onto the next face
        expected = GRADIENT * distance
        error = np.abs(filled - expected)[near] / expected[near]
        self.assertLessEqual(float(error.max()), EUCLID_TOLERANCE)
        chamfer = grid.sloped_fill(elevation, sea, GRADIENT)
        self.assertGreater(float(((chamfer - expected)[near] / expected[near]).max()), 0.05)

    def test_every_raised_cell_drains(self):
        grid = self.grid
        rng = np.random.default_rng(3)
        elevation = 200.0 + 300.0 * grid.smooth(rng.standard_normal(grid.n), 1.5)
        noise = 5.0 * grid.smooth(rng.standard_normal(grid.n), 1.0)
        sea = elevation < 150.0
        for metric in ('chamfer', 'euclid'):
            filled = grid.sloped_fill(elevation, sea, 2e-5, noise, metric=metric)
            raised = (filled > elevation) & ~sea
            around = np.where(grid.valid, np.where(sea[grid.safe], -np.inf, filled[grid.safe]), np.inf).min(1)
            self.assertTrue((around[raised] < filled[raised]).all(), metric)


@unittest.skipUnless(ts.PHYSICS.is_file(), 'needs terrain_sandbox.exe')
class SlopeLimitTest(unittest.TestCase):
    """The exe's threshold-slope limiter (--slope-limit-neighbours) on a cube of N 32: a high plateau around one sea
    cell at a face centre, one step without erosion or diffusion, is cut to a cone of slope S; with 8 neighbours its
    slope runs up to 8.2 % above S between the axes and the diagonals (an octagon), with 16 up to about 2.7 %."""

    N, SLOPE, HEIGHT = 32, 0.001, 1e5

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.grid = cube_domain.cube_table(ts.PHYSICS, cls.N, RADIUS_KM, cls.root)
        N = cls.N
        cls.spill = 0 * N * N + (N // 2) * N + N // 2
        sea = np.zeros(cls.grid.n, dtype='<f4')
        sea[cls.spill] = 1.0
        elevation = np.where(sea > 0, -100.0, cls.HEIGHT).astype('<f4')
        sea.tofile(cls.root / 'sea.f32')
        elevation.tofile(cls.root / 'elevation.f32')

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def cone(self, neighbours):
        out = self.root / f'out-{neighbours}.f32'
        command = [str(ts.PHYSICS), '--face-cells', str(self.N), '--elevation', str(self.root / 'elevation.f32'),
                   '--sea', str(self.root / 'sea.f32'), '--duration', '1', '--step', '1', '--k', '0',
                   '--diffusivity', '0', '--threshold-slope', str(self.SLOPE), '--radius-km', str(RADIUS_KM),
                   '--out-elevation', str(out), '--out-discharge', str(self.root / 'discharge.f32')]
        if neighbours is not None:
            command += ['--slope-limit-neighbours', str(neighbours)]
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        return np.fromfile(out, dtype='<f4').astype(np.float64)

    def slopes(self, surface):
        """The cone's slope over S by direction: a line fitted to height against the great-circle distance D on the
        spill's face, 4 to 14 cells out (the apex, rounded by the limiter's cut smoothing, left out), in 5-degree
        bins of the bearing folded into 0 to 45 degrees."""
        N, c = self.N, self.N // 2
        distance = self.grid.great_circle_km(np.arange(N * N), np.full(N * N, self.spill)) * 1000.0
        spacing = 2.0 * math.pi * RADIUS_KM * 1000.0 / (4 * N)
        row, column = np.divmod(np.arange(N * N), N)
        angle = np.degrees(np.arctan2(np.abs(row - c), np.abs(column - c)))
        angle = np.where(angle > 45.0, 90.0 - angle, angle)
        ring = (distance > 4 * spacing) & (distance < 14 * spacing)
        height = surface[:N * N]
        return np.array([np.polyfit(distance[ring & (angle >= lo) & (angle < lo + 5)],
                                    height[ring & (angle >= lo) & (angle < lo + 5)], 1)[0] / self.SLOPE
                         for lo in range(0, 45, 5)])

    def test_default_is_eight(self):
        np.testing.assert_array_equal(self.cone(None), self.cone(8))

    def test_sixteen_is_rounder(self):
        # Measured (master plan 5.3 review): 8 gives slopes 0.993 to 1.059 (spread 0.066), 16 0.994 to 1.009 (0.014).
        eight, sixteen = self.slopes(self.cone(8)), self.slopes(self.cone(16))
        self.assertGreater(float(eight.max() - eight.min()), 0.05)
        self.assertLess(float(sixteen.max() - sixteen.min()), 0.025)
        self.assertLess(float(np.abs(sixteen - 1.0).max()), 0.02)

    def test_refuses_other_counts(self):
        run = subprocess.run([str(ts.PHYSICS), '--face-cells', str(self.N), '--elevation', str(self.root / 'elevation.f32'),
                              '--sea', str(self.root / 'sea.f32'), '--slope-limit-neighbours', '12',
                              '--out-elevation', str(self.root / 'x.f32'), '--out-discharge', str(self.root / 'y.f32')],
                             capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)


if __name__ == '__main__':
    unittest.main()
