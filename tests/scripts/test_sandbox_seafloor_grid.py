"""The sea-floor grid-imprint fixes of the terrain sandbox (master plan 6.1b), on a synthetic lat-lon domain (no
executable needed):
- shelf_warp (M1): zero where the base is land or at the coast, the full warp below the shelf break, no sea cell read
  above the shelf edge or its own depth, deterministic, rows not wrapped across the poles;
- signed_distance_km: the distance of a cell's centre to a straight boundary, clipped;
- oceanic_weight (M2): in [0, 1], 0 and 1 beyond half the ramp from the indicator's boundary, its 0.5 contour off a
  one-cell staircase of the indicator; the column mode reads the oceanic share of the crust and refuses a bundle
  without it;
- volcano_cones (M3): the cell's height at its centre, nothing beyond radius_km, the tallest where two overlap."""

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "sandbox"))

import terrain_sandbox as ts  # noqa: E402

RADIUS_KM = 6371.0
COLUMNS, ROWS, K = 32, 16, 8
CELL_KM = 2.0 * 3.141592653589793 * RADIUS_KM / COLUMNS  # 1251 km: the test widths are in cells


def domain():
    """The whole lat-lon world of 32 x 16 cells at k 8 (256 x 128 fine cells)."""
    return ts.Domain("test", COLUMNS * K, ROWS * K, 0.0, 0.0, K, COLUMNS, ROWS, RADIUS_KM)


def sampler(d):
    """World fields (ROWS, COLUMNS) to the domain by nearest cell: the indicator's staircase kept exactly."""
    rows = np.minimum((d.y).astype(int), ROWS - 1)
    columns = np.minimum((d.x).astype(int), COLUMNS - 1)
    return lambda field: np.asarray(field, dtype=float)[rows[:, None], columns[None, :]]  # noqa: E731


WARP = {"amplitude_cells": 0.5, "wavelength_km": 4 * CELL_KM, "octaves": 2, "lacunarity": 2.0, "gain": 0.5, "min_cells": 2.5}
WEIGHT = {"mode": "distance", "ramp_km": CELL_KM, "smooth_km": 0.5 * CELL_KM, "column_share": 0.5}


class ShelfWarp(unittest.TestCase):
    def test_land_and_coast_keep_the_base(self):
        d = domain()
        x = np.broadcast_to(d.x[None, :], (d.rows, d.columns))
        base = np.where(x < 16, 500.0, -3000.0) + 0.0 * x  # land west, deep sea east
        base = base + np.where((x > 15.5) & (x < 17), 2500.0, 0.0)  # a shallow strip at the coast
        warped = ts.shelf_warp(d, base, WARP, 200.0, 1)
        land = base >= 0.0
        self.assertTrue(np.array_equal(warped[land], base[land]))
        self.assertTrue(np.all(np.isfinite(warped)))
        self.assertTrue(np.all(warped <= np.maximum(base, -200.0)))  # no sea cell read above the shelf edge or itself
        deep = base < -200.0
        self.assertGreater(np.abs(warped[deep] - base[deep]).max(), 100.0)  # the warp acts below the break
        self.assertTrue(np.array_equal(warped, ts.shelf_warp(d, base, WARP, 200.0, 1)))  # deterministic

    def test_shelf_profile_reads_the_warped_depth(self):
        spec = {"shelf_m": 150, "break_m": 200, "slope_m": 1000, "foot_m": 2200, "bottom_m": 6000}
        base = np.array([300.0, -100.0, -100.0, -3000.0])
        depth_base = np.array([300.0, -100.0, -1500.0, -3000.0])
        oceanic = np.zeros(4)
        plain = ts.shelf_profile(base, oceanic, spec)
        moved = ts.shelf_profile(base, oceanic, spec, depth_base)
        self.assertEqual(moved[0], 300.0)  # land stays the base's
        self.assertEqual(moved[1], plain[1])
        self.assertLess(moved[2], plain[2] - 1000.0)  # the profile of the warped depth
        self.assertTrue(np.array_equal(ts.shelf_profile(base, oceanic, spec, None), plain))

    def test_rows_do_not_wrap_across_the_poles(self):
        """A whole lat-lon world: the bottom rows read no value of the top rows (rows clamped, columns wrapped). With
        seed 1 the bottom row moves up to 0.6 fine rows south: wrapped rows read -4009 m there (checked on the
        former path, displaced)."""
        d = domain()
        base = np.full((d.rows, d.columns), -1000.0)
        base[:2 * K, :] = -6000.0
        warped = ts.shelf_warp(d, base, WARP, 200.0, 1)
        self.assertTrue(np.allclose(warped[-K // 2:], -1000.0))  # bilinear of a constant: up to rounding
        self.assertLess(warped[2 * K:].min(), -1000.0)  # the warp does reach across the step


class SignedDistance(unittest.TestCase):
    def test_distance_to_a_straight_boundary(self):
        """Inside from fine column 128 on: column 131's centre lies 3.5 fine cells inside, column 126's 1.5 outside."""
        d = domain()
        inside = np.zeros((d.rows, d.columns), dtype=bool)
        inside[:, 128:] = True
        distance = ts.signed_distance_km(d, inside, 1e9)
        self.assertTrue(np.allclose(distance[:, 131], 3.5 * d.fine_km))
        self.assertTrue(np.allclose(distance[:, 126], -1.5 * d.fine_km))
        self.assertTrue(np.allclose(distance[:, 128], 0.5 * d.fine_km))
        self.assertTrue(np.allclose(distance[:, 127], -0.5 * d.fine_km))
        clipped = ts.signed_distance_km(d, inside, 2.0 * d.fine_km)
        self.assertTrue(np.allclose(clipped[:, 131], 1.5 * d.fine_km))  # min(4, 2) cells less half a cell
        self.assertTrue(np.all(np.abs(clipped) <= 2.0 * d.fine_km))


class OceanicWeight(unittest.TestCase):
    def world(self):
        """Crust class 1 east of a one-cell staircase (column 16 + row // 2); columns 9 to 11 class 3 with an
        oceanic column."""
        rows, columns = np.indices((ROWS, COLUMNS))
        oceanic = columns >= 16 + rows // 2
        strip = (columns >= 9) & (columns <= 11)
        return {"crust_class": np.where(oceanic, 1, 3).astype(np.uint8),
                "oceanic_crust": np.where(oceanic, 7000.0, 0.0) + np.where(strip, 9000.0, 0.0),
                "continental_crust": np.where(oceanic, 1e-18, np.where(strip, 0.0, 35000.0))}

    def test_range_and_reach(self):
        d = domain()
        weight = ts.oceanic_weight(d, self.world(), sampler(d), WEIGHT)
        self.assertGreaterEqual(weight.min(), 0.0)
        self.assertLessEqual(weight.max(), 1.0)
        x = np.broadcast_to(d.x[None, :], weight.shape)
        boundary = 16 + np.floor(np.broadcast_to(d.y[:, None], weight.shape)).astype(int) // 2
        far = (np.abs(x - boundary) * d.cell_km > 0.5 * WEIGHT["ramp_km"] + 3 * WEIGHT["smooth_km"])
        self.assertTrue(np.all(weight[far & (x > boundary)] == 1.0))
        self.assertTrue(np.all(weight[far & (x < boundary)] == 0.0))

    def test_contour_leaves_the_staircase(self):
        """Along each fine row the 0.5 crossing moves smoothly with the row: its steps between adjacent rows are far
        below the staircase's one cell (K fine cells) every two cells of rows."""
        d = domain()
        weight = ts.oceanic_weight(d, self.world(), sampler(d), WEIGHT)
        band = slice(2 * K, (ROWS - 2) * K)  # rows away from the poles
        crossing = np.array([np.argmax(row >= 0.5) for row in weight[band, : 30 * K]])
        self.assertLess(np.abs(np.diff(crossing)).max(), K / 2)

    def test_column_mode_reads_the_crust(self):
        d = domain()
        world = self.world()
        cell = ts.oceanic_weight(d, world, sampler(d), WEIGHT)
        column = ts.oceanic_weight(d, world, sampler(d), {**WEIGHT, "mode": "column"})
        x = np.broadcast_to(d.x[None, :], cell.shape)
        strip = (x > 9.9) & (x < 10.1)  # class 3 with an oceanic column: oceanic only by the column
        self.assertTrue(np.all(cell[strip] == 0.0))
        self.assertGreater(column[strip].min(), 0.4)
        self.assertTrue(np.array_equal(cell[x > 28], column[x > 28]))

    def test_column_mode_needs_the_oceanic_crust(self):
        d = domain()
        world = {**self.world(), "oceanic_crust": None}
        with self.assertRaises(ValueError):
            ts.oceanic_weight(d, world, sampler(d), {**WEIGHT, "mode": "column"})

    def test_unknown_mode(self):
        d = domain()
        with self.assertRaises(ValueError):
            ts.oceanic_weight(d, self.world(), sampler(d), {**WEIGHT, "mode": "nearest"})


class VolcanoCones(unittest.TestCase):
    def test_cone_shape(self):
        d = domain()
        heights = np.zeros((ROWS, COLUMNS))
        heights[8, 5] = 3000.0
        heights[8, 6] = 1000.0  # one cell east: overlapping cones keep the tallest
        world = {"volcano": heights, "elevation": np.zeros((ROWS, COLUMNS)), "grid": None}
        spec = {"radius_km": 900.0, "profile": 1.3}
        cones = ts.volcano_cones(d, world, spec)
        centres = ts.world_centres(world) * RADIUS_KM
        points = d.points().reshape(3, -1)
        chord = np.linalg.norm(points - centres[:, 8 * COLUMNS + 5, None], axis=0)
        other = np.linalg.norm(points - centres[:, 8 * COLUMNS + 6, None], axis=0)
        flat = cones.reshape(-1)
        expected = np.maximum(3000.0 * np.clip(1 - chord / 900.0, 0, 1) ** 1.3,
                              1000.0 * np.clip(1 - other / 900.0, 0, 1) ** 1.3)
        self.assertTrue(np.allclose(flat, expected))
        self.assertTrue(np.all(flat[(chord >= 900.0) & (other >= 900.0)] == 0.0))
        self.assertLessEqual(flat.max(), 3000.0)
        self.assertGreater(flat.max(), 3000.0 * (1.0 - d.fine_km / 900.0) ** 1.3)  # a fine cell near the centre


if __name__ == "__main__":
    unittest.main()
