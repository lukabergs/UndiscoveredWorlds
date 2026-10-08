"""scripts/evaluation/maps.py: the checkpoint time (the evolver's property before the tectonic snapshot) and
the outlines of the regional terrain map."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

import numpy as np  # noqa: E402

from uwtools.bundle import Bundle, FieldInfo  # noqa: E402
from maps import LAKE_OUTLINE, SNAPSHOT_TABLE, checkpoint_time, regional_terrain  # noqa: E402


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
        self.assertIn("outline", note)


if __name__ == "__main__":
    unittest.main()
