"""nested_window.parent_grid (master plan 4.2): a cube parent's pixels are its window's lat-lon export, the same
pixels (centres, latitudes, longitudes) a lat-lon run of the window keeps, also across the date line; the global
domain is refused as a cube parent; a lat-lon parent is terrain_sandbox.make_domains' window.

  python -m unittest tests/scripts/test_nested_window_cube.py"""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts' / 'sandbox'))
import nested_window  # noqa: E402
import terrain_sandbox as ts  # noqa: E402


class FakeBundle:
    """The attributes parent_grid and make_domains read of a cube base bundle of N 64 (256 export columns)."""
    face_cells, columns, rows, grid_kind, manifest = 64, 256, 128, 'cube', {}


def config(cube, pixels, scale=64, margin=1):
    return {'cube': {'enabled': cube}, 'physics': {'window': {'k': 1.0}},
            'domains': {'global': {'columns': 2048}, 'W': {'pixels': pixels, 'scale': scale, 'margin': margin}}}


class ParentGrid(unittest.TestCase):
    def check_same_pixels(self, pixels, scale):
        cube = nested_window.parent_grid(config(True, pixels, scale), FakeBundle(), 'W')
        latlon = nested_window.parent_grid(config(False, pixels, scale), FakeBundle(), 'W')
        r0, r1, c0, c1 = latlon.crop
        self.assertEqual(cube.crop, (0, cube.rows, 0, cube.columns))
        self.assertEqual((cube.rows, cube.columns), (r1 - r0, c1 - c0))
        np.testing.assert_array_equal(cube.x, latlon.x[c0:c1])
        np.testing.assert_array_equal(cube.y, latlon.y[r0:r1])
        np.testing.assert_allclose(cube.latitude_deg, latlon.latitude_deg[r0:r1], rtol=0, atol=1e-12)
        np.testing.assert_allclose(cube.longitude_deg, latlon.longitude_deg[c0:c1], rtol=0, atol=1e-12)
        self.assertEqual((cube.k, cube.world_columns, cube.world_rows), (scale, 256, 128))
        self.assertAlmostEqual(cube.fine_km, latlon.fine_km, places=12)
        return cube

    def test_cube_parent_is_the_lat_lon_window_kept(self):
        cube = self.check_same_pixels([185, 15, 188, 18], 256)
        self.assertEqual((cube.rows, cube.columns), (1024, 1024))
        self.assertAlmostEqual(cube.x[0], 185 + 0.5 / 256, places=12)

    def test_across_the_date_line(self):
        cube = self.check_same_pixels([254, 60, 1, 63], 64)
        self.assertEqual(cube.columns, 4 * 64)
        self.assertGreater(cube.x[-1], 256.0)  # past column W - 1, sampled with the columns wrapped

    def test_global_is_no_cube_parent(self):
        with self.assertRaises(SystemExit):
            nested_window.parent_grid(config(True, [0, 0, 1, 1]), FakeBundle(), 'global')

    def test_lat_lon_parent_is_make_domains(self):
        spec = config(False, [45, 25, 48, 28], 256)
        got = nested_window.parent_grid(spec, FakeBundle(), 'W')
        expected = ts.make_domains(spec, FakeBundle(), ['W'])[0]
        self.assertEqual((got.columns, got.rows, got.crop, got.k), (expected.columns, expected.rows, expected.crop,
                                                                    expected.k))
        np.testing.assert_array_equal(got.x, expected.x)
        np.testing.assert_array_equal(got.y, expected.y)


if __name__ == '__main__':
    unittest.main()
