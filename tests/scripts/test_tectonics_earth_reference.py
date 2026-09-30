import sys
from pathlib import Path
import unittest
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from PIL import Image, TiffImagePlugin

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/benchmarks"))
from tectonics_earth_reference import measure_earth, weighted_distribution, compare_to_earth, load_earth_reference


class EarthHypsographyTests(unittest.TestCase):
    def test_cache_reuses_measures_and_rejects_changed_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            heightmap = root / "refs/processed/heightmap/8.tif"
            mask = root / "refs/processed/land_ocean/maps/earth_land_ocean_8x4.png"
            heightmap.parent.mkdir(parents=True)
            mask.parent.mkdir(parents=True)
            values = np.full((4, 8), -4000, dtype=np.float32)
            values[1] = 1000
            tags = TiffImagePlugin.ImageFileDirectory_v2()
            tags[33550] = (45., 45., 0.)
            tags[33922] = (0., 0., 0., -180., 90., 0.)
            Image.fromarray(values).save(heightmap, tiffinfo=tags)
            Image.fromarray(np.uint8(values > 0) * 255).save(mask)
            with patch("tectonics_earth_reference.measure_earth", wraps=measure_earth) as measure:
                first = load_earth_reference(root, 8)
                self.assertEqual(load_earth_reference(root, 8), first)
                self.assertEqual(measure.call_count, 1)
                (root / "refs/README.md").write_text("changed provenance")
                with self.assertRaisesRegex(ValueError, "Stale Earth"):
                    load_earth_reference(root, 8)
                refreshed = load_earth_reference(root, 8, refresh=True)
                self.assertEqual(measure.call_count, 2)
                self.assertNotEqual(first["cache"]["sha256"], refreshed["cache"]["sha256"])
                self.assertEqual(first["metrics"], refreshed["metrics"])

    def test_weights_and_mask_preserve_below_sea_level_land(self):
        z = np.full((4, 8), -4000.)
        land = np.zeros((4, 8), dtype=bool)
        land[0] = True
        z[0] = -100.  # Depression still belongs to the non-ocean mask.
        m = measure_earth(z, land)
        self.assertAlmostEqual(m["land_fraction"], (1 - np.sqrt(.5)) / 2)
        self.assertEqual(m["land_elevation_m"]["mean"], -100)
        self.assertEqual(m["ocean_depth_m"]["mean"], 4000)
        self.assertEqual(m["hypsography"]["elevation_m"][100], -4000)
        self.assertAlmostEqual(compare_to_earth(m, {"metrics": m})["hypsography_mean_absolute_quantile_gap_m"], 0)

    def test_invalid_grids_and_values_fail(self):
        with self.assertRaisesRegex(ValueError, "2:1"):
            measure_earth(np.ones((3, 4)), np.ones((3, 4)))
        for values, weights in (([np.nan], [1]), ([1], [0]), ([1, 2], [1])):
            with self.assertRaises(ValueError):
                weighted_distribution(values, weights)

    def test_weighted_quantiles_do_not_count_pixels_equally(self):
        statistics, curve = weighted_distribution([0, 1000], [9, 1])
        self.assertEqual(statistics["mean"], 100)
        self.assertEqual(statistics["median"], 0)
        self.assertEqual(curve["elevation_m"][-1], 1000)


if __name__ == "__main__":
    unittest.main()
