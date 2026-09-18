"""Deterministic checks for processed-reference spherical grids."""

import sys
from pathlib import Path
import unittest

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "refs"))
import cell_grid


class CellGridTests(unittest.TestCase):
    def test_dimensions_and_latitude_centres(self):
        self.assertEqual(cell_grid.grid_dimensions(512), (512, 256))
        with self.assertRaises(ValueError):
            cell_grid.grid_dimensions(512, 255)
        latitudes = cell_grid.latitude_centres(256)
        self.assertLess(latitudes[0], 90.0)
        self.assertGreater(latitudes[-1], -90.0)
        np.testing.assert_allclose(latitudes, -latitudes[::-1], atol=1.0e-12)

    def test_area_and_conservative_remapping(self):
        source = np.arange(32 * 16, dtype=np.float32).reshape(16, 32) ** 2
        destination = cell_grid.conservative_remap(source, 8)
        source_integral = float(np.sum(source * cell_grid.spherical_area_weights(32)))
        destination_integral = float(
            np.sum(destination * cell_grid.spherical_area_weights(8))
        )
        self.assertAlmostEqual(destination_integral / source_integral, 1.0, places=7)

    def test_native_coordinate_ingestion_isolated_at_boundary(self):
        source = np.full((3, 4), 7.25, dtype=np.float32)
        result = cell_grid.conservative_remap(
            source,
            8,
            source_latitudes=np.array([90.0, 0.0, -90.0]),
            source_longitudes=np.array([0.0, 90.0, 180.0, 270.0]),
        )
        self.assertEqual(result.shape, (4, 8))
        np.testing.assert_allclose(result, 7.25, atol=1.0e-6)

    def test_partial_latitude_coverage_does_not_invent_polar_data(self):
        result = cell_grid.conservative_remap(
            np.ones((2, 4), dtype=np.float32), 8,
            source_latitudes=np.array([20., -20.]),
            source_longitudes=np.array([0., 90., 180., 270.]),
        )
        self.assertTrue(np.isnan(result[[0, -1]]).all())
        np.testing.assert_array_equal(result[1:3], 1)

    def test_float32_native_longitudes_are_accepted(self):
        width = 4320
        longitudes = (np.arange(width) * 360 / width - 180).astype(np.float32)
        result = cell_grid.conservative_remap(np.ones((2, width)), 8,
                                               source_longitudes=longitudes)
        np.testing.assert_allclose(result, 1)


if __name__ == "__main__":
    unittest.main()
