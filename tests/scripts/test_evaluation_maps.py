"""The checkpoint time of scripts/evaluation/maps.py: the evolver's property before the tectonic snapshot."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

from uwtools.bundle import Bundle  # noqa: E402
from maps import SNAPSHOT_TABLE, checkpoint_time  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
