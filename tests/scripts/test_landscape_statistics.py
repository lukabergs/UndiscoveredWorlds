"""The drainage_pattern, ocean_depths and glaciation validation kinds of
scripts/evaluation/landscape_statistics.py on a tiny synthetic bundle (8 x 4; erosion plan, Phase C,
step C0, and Phase G, step G3)."""

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

from landscape_statistics import depth_shares, measure_drainage_pattern, measure_glaciation  # noqa: E402
from test_river_mouths import COLUMNS, RADIUS_M, ROWS, write_bundle  # noqa: E402


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
        # Rows are centred at 67.5, 22.5, -22.5 and -67.5 degrees: rows 1 and 2 are the 0-30 band
        # (row 1's 7 westward cells of 14), rows 0 and 3 the 60-90 band (northward), none 30-60.
        self.assertAlmostEqual(metrics["share_east_west_0_30"], 7 / 14)
        self.assertIsNone(metrics["share_east_west_30_60"])
        self.assertAlmostEqual(metrics["share_east_west_60_90"], 0.0)
        self.assertNotIn("direction_deviation_deg", metrics)  # no filled elevation

    def test_direction_deviation(self):
        """The filled surface rises 100 m per column east of the sea column, so the downhill
        direction is west, except at column 7, whose east neighbour across the date line is the
        sea (east). Rows 0, 2 and 3 drain west: 0 degrees, 180 at column 7. Row 1 drains north: 90
        everywhere. The cell x 3 of row 2 drains north-west: atan2(dy, mean dx of rows 1 and 2) =
        atan2(pi / 4, cos(22.5) pi / 4) from west."""
        world = drainage_world()
        filled = np.zeros((ROWS, COLUMNS), dtype=np.float32)
        filled[:, 0] = -100.0
        filled[:, 1:] = 100.0 * np.arange(1, COLUMNS)
        direction = np.full((ROWS, COLUMNS), 7, dtype=np.uint8)
        direction[:, 0] = 0
        direction[1, 1:] = 1
        direction[2, 3] = 8
        world["hydrosphere.filled_elevation"] = filled
        world["hydrosphere.flow_direction"] = direction
        with tempfile.TemporaryDirectory() as directory:
            metrics = measure_drainage_pattern(write_bundle(directory, world))
        diagonal = np.degrees(np.arctan2(1.0, np.cos(np.radians(22.5))))
        self.assertAlmostEqual(metrics["direction_deviation_deg"], (3 * 180.0 + 7 * 90.0 + diagonal) / 28)
        self.assertAlmostEqual(metrics["share_east_west_0_30"], 6 / 14)
        self.assertAlmostEqual(metrics["share_east_west_60_90"], 1.0)


def glacial_world():
    """Sea level 10 m; column 0 is sea, columns 1..7 land (28 cells). Row 0 (67.5 N) is glaciated
    (0.5) with 2 m of glacial denudation and holds a lake at x 2 and present-day ice at x 1..3;
    every land cell has 4 m of denudation. Sea cells: (0, 0) has 5 m of glacial denudation at 110 m
    below the sea level, (1, 0) 1 m at 40 m below it and glaciation 0.25, (3, 0) 5 m but a lake
    surface (not sea). The ELA is 1000 m more per row from 1010."""
    elevation = np.full((ROWS, COLUMNS), 200.0, dtype=np.float32)
    elevation[:, 0] = -100.0
    elevation[1, 0] = -30.0
    glaciation = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    glaciation[0, 1:] = 0.5
    glaciation[1, 0] = 0.25
    ice = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    ice[0, 1:4] = 100.0
    glacial = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    glacial[0, 1:] = 2.0
    glacial[0, 0] = glacial[3, 0] = 5.0
    glacial[1, 0] = 1.0
    denudation = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    denudation[:, 1:] = 4.0
    lake = np.zeros((ROWS, COLUMNS), dtype=np.float32)
    lake[0, 2] = 250.0
    lake[3, 0] = -50.0
    ela = (1010.0 + 1000.0 * np.arange(ROWS, dtype=np.float32))[:, None] * np.ones((1, COLUMNS), dtype=np.float32)
    return {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(10, dtype=np.int32),
            "planetbody.radius": np.array(RADIUS_M, dtype=np.float64),
            "cryosphere.glaciation": glaciation, "cryosphere.ice_thickness": ice,
            "cryosphere.equilibrium_line_altitude": ela.astype(np.float32),
            "geosphere.glacial_denudation": glacial, "geosphere.denudation": denudation,
            "hydrosphere.lake_surface": lake}


class GlaciationTests(unittest.TestCase):
    def test_shares_sea_cells_bands_and_denudation(self):
        with tempfile.TemporaryDirectory() as directory:
            metrics = measure_glaciation(write_bundle(directory, glacial_world()))
        polar, tropical = np.cos(np.radians(67.5)), np.cos(np.radians(22.5))  # row weights
        land = 14 * (polar + tropical)
        self.assertEqual(metrics["land_cells"], 28)
        self.assertEqual(metrics["glaciated_land_cells"], 7)
        self.assertEqual(metrics["ice_covered_land_cells"], 3)
        self.assertAlmostEqual(metrics["glaciated_land_share"], 7 * polar / land)
        self.assertAlmostEqual(metrics["land_mean_glaciation"], 0.5 * 7 * polar / land)
        self.assertAlmostEqual(metrics["ice_covered_land_share"], 3 * polar / land)
        self.assertAlmostEqual(metrics["maximum_ice_thickness_m"], 100.0)
        self.assertEqual(metrics["glaciated_sea_cells"], 2)
        sphere_km2 = 4 * np.pi * (RADIUS_M / 1000) ** 2
        self.assertAlmostEqual(metrics["glaciated_sea_area_km2"], sphere_km2 / 16, places=3)
        self.assertAlmostEqual(metrics["glaciated_sea_min_depth_m"], 40.0)
        self.assertAlmostEqual(metrics["glaciated_sea_median_depth_m"], 75.0)
        self.assertEqual(metrics["glaciated_sea_cells_deeper_than_100_m"], 1)
        self.assertEqual(metrics["glaciated_cells_below_sea_level"], 1)
        self.assertAlmostEqual(metrics["maximum_glacial_denudation_m"], 5.0)
        self.assertEqual(metrics["glaciated_lake_cells"], 1)
        self.assertEqual(metrics["ice_on_lake_cells"], 1)
        # Rows are centred at 67.5, 22.5, -22.5 and -67.5 degrees; the ELA less the sea level.
        self.assertAlmostEqual(metrics["ela_above_sea_m_60N_75N"], 1000.0)
        self.assertAlmostEqual(metrics["ela_above_sea_m_15N_30N"], 2000.0)
        self.assertAlmostEqual(metrics["ela_above_sea_m_30S_15S"], 3000.0)
        self.assertAlmostEqual(metrics["ela_above_sea_m_75S_60S"], 4000.0)
        self.assertIsNone(metrics["ela_above_sea_m_0N_15N"])
        self.assertIsNone(metrics["ela_above_sea_m_75N_90N"])
        self.assertAlmostEqual(metrics["land_mean_glacial_denudation_m"], 2.0 * 7 * polar / land)
        self.assertAlmostEqual(metrics["land_mean_denudation_m"], 4.0)
        self.assertAlmostEqual(metrics["glacial_denudation_share"], 0.5 * 7 * polar / land)
        # An ELA of zero everywhere is the term switched off: no band values.
        world = glacial_world()
        world["cryosphere.equilibrium_line_altitude"] = np.zeros((ROWS, COLUMNS), dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            metrics = measure_glaciation(write_bundle(directory, world))
        self.assertIsNone(metrics["ela_above_sea_m_60N_75N"])


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
