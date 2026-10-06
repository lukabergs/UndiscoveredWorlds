"""The drainage_pattern and ocean_depths validation kinds of scripts/evaluation/landscape_statistics.py
on a tiny synthetic bundle (8 x 4; erosion plan, Phase C, step C0)."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))
sys.path.insert(0, str(TESTS))

from landscape_statistics import depth_shares, measure_drainage_pattern  # noqa: E402
from test_river_mouths import COLUMNS, ROWS, write_bundle  # noqa: E402


def drainage_world():
    """Column 0 is sea, columns 1..7 land (28 land cells). Every land cell drains north (code 1)
    except row 1, which drains west (7) to the sea, and the cell x 3 of row 2, which drains
    north-west (8) into row 1."""
    elevation = np.full((ROWS, COLUMNS), 200.0, dtype=np.float32)
    elevation[:, 0] = -100.0
    direction = np.zeros((ROWS, COLUMNS), dtype=np.uint8)
    direction[:, 1:] = 1
    direction[1, 1:] = 7
    direction[2, 3] = 8
    lake = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    lake[3, 5] = 250.0
    special = np.zeros((ROWS, COLUMNS), dtype=np.int16)
    special[3, 5], special[3, 6], special[2, 6] = 100, 110, 131
    basin = np.zeros((ROWS, COLUMNS), dtype=np.uint8)
    basin[3, 4:7] = 2
    delta = np.zeros((ROWS, COLUMNS), dtype=np.int16)
    delta[1, 0] = 7  # deltas lie on sea cells as well
    rate = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    rate[:, 1:] = 2.0
    rate[1, 1:] = 9.0
    return {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(0, dtype=np.int32),
            "hydrosphere.flow_direction": direction, "hydrosphere.lake_surface": lake,
            "hydrosphere.special_features": special, "hydrosphere.basin_class": basin,
            "hydrosphere.delta_direction": delta, "geosphere.erosion_rate": rate}


class DrainagePatternTests(unittest.TestCase):
    def test_counts_and_shares(self):
        with tempfile.TemporaryDirectory() as directory:
            metrics = measure_drainage_pattern(write_bundle(directory, drainage_world()))
        self.assertEqual(metrics["land_cells"], 28)
        self.assertAlmostEqual(metrics["share_east_west"], 7 / 28)
        # Land receivers: row 0 has none (north of the pole), the outlet x 1 of row 1 drains to the
        # sea: 28 - 7 - 1 = 20. Straight: row 1 x 2..7 (6); row 2 drains into row 1's westward cells
        # (0); row 3 drains north into row 2's northward cells except x 3, whose receiver drains
        # north-west (6).
        self.assertEqual(metrics["cells_with_land_receiver"], 20)
        self.assertAlmostEqual(metrics["share_straight"], 12 / 20)
        self.assertEqual(metrics["lake_cells"], 1)
        self.assertEqual(metrics["salt_lake_cells"], 1)
        self.assertEqual(metrics["salt_pan_cells"], 1)
        self.assertEqual(metrics["wetland_cells"], 1)
        self.assertEqual(metrics["endorheic_cells"], 3)
        self.assertEqual(metrics["delta_cells"], 1)
        self.assertAlmostEqual(metrics["land_mean_erosion_rate_m_myr"], (7 * 9.0 + 21 * 2.0) / 28)


class OceanDepthTests(unittest.TestCase):
    def test_bands_mean_and_median(self):
        depth = np.array([100.0, 500.0, 4500.0, 7000.0])
        weights = np.array([1.0, 1.0, 2.0, 4.0])
        shares = depth_shares(depth, weights)
        self.assertAlmostEqual(shares["ocean 0-0.2 km percent"], 12.5)
        self.assertAlmostEqual(shares["ocean 0.2-1 km percent"], 12.5)
        self.assertAlmostEqual(shares["ocean 4-5 km percent"], 25.0)
        self.assertAlmostEqual(shares["ocean over 6 km percent"], 50.0)
        self.assertAlmostEqual(shares["ocean 1-2 km percent"], 0.0)
        self.assertAlmostEqual(shares["mean_depth_m"], (100 + 500 + 9000 + 28000) / 8)
        self.assertEqual(shares["median_depth_m"], 4500.0)


if __name__ == "__main__":
    unittest.main()
