"""scripts/evaluation/maps.py: the checkpoint time (the evolver's property before the tectonic snapshot) and
the outlines of the regional terrain map."""

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

import numpy as np  # noqa: E402

from uwtools.bundle import Bundle, FieldInfo  # noqa: E402
from maps import (FADE_FLOW_OPAQUE, FADE_RIVER, LAKE_OUTLINE, SNAPSHOT_TABLE, checkpoint_time,  # noqa: E402
                  pixel_geometry, regional_terrain, regional_terrain_fade)


def bundle(directory, properties, snapshot_time=None):
    manifest = {"properties": properties}
    result = Bundle(Path(directory), manifest)
    if snapshot_time is not None:
        (Path(directory) / "snapshot.json").write_text(json.dumps({"time_myr": snapshot_time}), encoding="utf-8")
        result.tables[SNAPSHOT_TABLE] = {"id": SNAPSHOT_TABLE, "file": "snapshot.json"}
    return result


class CheckpointTimeTests(unittest.TestCase):
    def test_property_comes_before_the_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            item = bundle(directory, {"checkpoint_process": "terrain.fluvial_erosion", "checkpoint_time_myr": "13"}, 50.0)
            self.assertEqual(checkpoint_time(item), 13.0)

    def test_snapshot_without_the_property(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(checkpoint_time(bundle(directory, {"command": "run"}, 50.0)), 50.0)

    def test_neither(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertIsNone(checkpoint_time(bundle(directory, {})))


def regional_bundle(directory, values):
    """A bundle of the given fields (id -> array) as `world_pipeline region` names them."""
    result = Bundle(Path(directory), {})
    for identifier, array in values.items():
        relative = f"{identifier}.npy"
        np.save(Path(directory) / relative, array)
        result.fields[identifier] = FieldInfo(identifier, identifier.split(".")[0], "grid", str(array.dtype),
                                              array.shape, "", relative, "", "", "")
    return result


class PixelGeometryTests(unittest.TestCase):
    def test_a_resampled_bundle_takes_metres_per_pixel(self):
        """44 rows of 20 km pixels in an 8-degree box: the pixel is 20 km, not 8/44 degrees."""
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, {"geosphere.base_elevation": np.zeros((44, 44), dtype=np.float32)})
            item.manifest["window"] = {"north": 18.0, "south": 10.0, "west": 0.0, "east": 8.0, "pixel_columns": 44,
                                       "pixel_rows": 44, "pixel_cells": 1.11, "resampled": True, "metres_per_pixel": 20000.0}
            north, pixel_deg, radius_m = pixel_geometry(item)
            self.assertEqual((north, radius_m), (18.0, 6371000.0))
            self.assertAlmostEqual(pixel_deg, math.degrees(20000.0 / 6371000.0), places=12)
            del item.manifest["window"]["pixel_columns"]
            self.assertAlmostEqual(pixel_geometry(item)[1], 8.0 / 44, places=12)


class RegionalTerrainOutlineTests(unittest.TestCase):
    """Erosion plan R1: the landmass outline (black, on the land side) and the lake shores (on the
    lake side), 4-neighbours, the window's edge not a boundary, nothing wrapping, over the rivers."""

    def test_outline_cells(self):
        rows, columns = 6, 8
        sea = np.zeros((rows, columns), dtype=np.uint8)
        sea[:, :2] = 1                                   # columns 0 and 1: sea
        surface = np.zeros((rows, columns), dtype=np.float32)
        surface[1:4, 5:8] = 50.0                         # a 3 x 3 lake against the east edge
        direction = np.zeros((rows, columns), dtype=np.uint8)
        direction[4, 2] = 7                              # a river on the coast cell (4, 2)
        direction[4, 4] = 7                              # and on the inland cell (4, 4)
        direction[2, 6] = 7                              # and through the lake's interior cell (2, 6)
        flow = (direction > 0).astype(np.int32)
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, {
                "geosphere.base_elevation": np.where(sea != 0, -100.0, 100.0).astype(np.float32),
                "geosphere.sea_level": np.array(0.0, dtype=np.float32),
                "geosphere.sea": sea,
                "cryosphere.ice_thickness": np.zeros((rows, columns), dtype=np.float32),
                "hydrosphere.lake_surface": surface,
                "hydrosphere.channel_direction": direction,
                "hydrosphere.discharge": np.where(direction > 0, 100.0, 0.0).astype(np.float32),
                "hydrosphere.river_flow_january": flow,
                "hydrosphere.river_flow_july": flow,
            })
            rgb, note = regional_terrain(item)
        black = np.all(rgb == 0, axis=2)
        expected_black = np.zeros((rows, columns), dtype=bool)
        expected_black[:, 2] = True                      # every land cell beside the sea, (4, 2)'s river included
        np.testing.assert_array_equal(black, expected_black)
        shore = np.all(rgb == np.array(LAKE_OUTLINE, dtype=np.uint8), axis=2)
        expected_shore = surface > 0
        expected_shore[2, 6] = False                     # the lake's interior cell
        expected_shore[2, 7] = False                     # beside the window's east edge only: not a shore
        np.testing.assert_array_equal(shore, expected_shore)
        self.assertFalse(black[4, 4] or shore[4, 4])     # the inland river keeps its colour
        self.assertFalse(np.array_equal(rgb[2, 6], rgb[4, 4]))  # lakes are drawn over the rivers

    def test_minimum_flow_option_culls_small_rivers(self):
        rows, columns = 3, 6
        direction = np.zeros((rows, columns), dtype=np.uint8)
        direction[1, 1:5] = 3                            # four river cells, annual flows 5, 100, 150 and 1000
        january = np.zeros((rows, columns), dtype=np.int32)
        july = np.zeros((rows, columns), dtype=np.int32)
        january[1, 1:5] = [5, 100, 100, 1000]
        july[1, 1:5] = [5, 100, 200, 1000]
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, {
                "geosphere.base_elevation": np.full((rows, columns), 100.0, dtype=np.float32),
                "geosphere.sea_level": np.array(0.0, dtype=np.float32),
                "geosphere.sea": np.zeros((rows, columns), dtype=np.uint8),
                "cryosphere.ice_thickness": np.zeros((rows, columns), dtype=np.float32),
                "hydrosphere.lake_surface": np.zeros((rows, columns), dtype=np.float32),
                "hydrosphere.channel_direction": direction,
                "hydrosphere.discharge": np.where(direction > 0, 100.0, 0.0).astype(np.float32),
                "hydrosphere.river_flow_january": january,
                "hydrosphere.river_flow_july": july,
            })
            every, _ = regional_terrain(item)
            culled, note = regional_terrain(item, minimum_flow=150)
        land = every[0, 1]
        for column in range(1, 5):
            self.assertFalse(np.array_equal(every[1, column], land))   # without the option every flowing channel is drawn
        np.testing.assert_array_equal(culled[1, 1], land)             # 5: culled
        np.testing.assert_array_equal(culled[1, 2], land)             # 100: culled
        self.assertFalse(np.array_equal(culled[1, 3], land))          # (100 + 200) / 2 = 150: drawn
        self.assertFalse(np.array_equal(culled[1, 4], land))          # 1000: drawn
        self.assertIn('at least 150', note)

    def test_fade_rivers_grow_opaque_with_flow(self):
        rows, columns = 3, 6
        direction = np.zeros((rows, columns), dtype=np.uint8)
        direction[1, 1:5] = 3                            # four river cells, flows 5, 10, 300 and 1000
        flow = np.zeros((rows, columns), dtype=np.int32)
        flow[1, 1:5] = [5, 10, 300, 1000]
        surface = np.zeros((rows, columns), dtype=np.float32)
        surface[1, 5] = 150.0                            # a lake cell on the river's last column
        direction[1, 5] = 3
        flow[1, 5] = 1000
        with tempfile.TemporaryDirectory() as directory:
            item = regional_bundle(directory, {
                "geosphere.base_elevation": np.full((rows, columns), 100.0, dtype=np.float32),
                "geosphere.sea_level": np.array(0.0, dtype=np.float32),
                "geosphere.sea": np.zeros((rows, columns), dtype=np.uint8),
                "cryosphere.ice_thickness": np.zeros((rows, columns), dtype=np.float32),
                "hydrosphere.lake_surface": surface,
                "hydrosphere.channel_direction": direction,
                "hydrosphere.discharge": np.where(direction > 0, 100.0, 0.0).astype(np.float32),
                "hydrosphere.river_flow_january": flow,
                "hydrosphere.river_flow_july": flow,
            })
            rgb, note = regional_terrain_fade(item)
        land = rgb[0, 1]
        np.testing.assert_array_equal(rgb[1, 1], land)   # below the transparent flow: the terrain
        np.testing.assert_array_equal(rgb[1, 2], land)   # at it: still the terrain
        river = np.array(FADE_RIVER, dtype=np.uint8)
        np.testing.assert_array_equal(rgb[1, 3], river)  # at the opaque flow: the river colour
        np.testing.assert_array_equal(rgb[1, 4], river)
        self.assertFalse(np.array_equal(rgb[1, 5], river))  # the lake is drawn over the river
        self.assertIn(f"{FADE_FLOW_OPAQUE:g}", note)
        self.assertIn("outline", note)


if __name__ == "__main__":
    unittest.main()
