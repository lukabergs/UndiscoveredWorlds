"""scripts/uwtools/naming.py: tags, ids, the scratch registry, find and previous (master plan 2.2).

The main tree and the storage root are temporary directories: nothing is written to the real
registry, to out/ or to the storage drive, and nothing is routed or linked.
"""

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from uwtools import naming, runs, storage  # noqa: E402

CODE = {"revision": "0" * 40, "dirty": False}


class NamingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.store = Path(self.temporary.name)
        self.root = self.store / "main"
        (self.root / "runs/registry").mkdir(parents=True)
        (self.root / "runs/registry/scratch.json").write_text(
            json.dumps({"schema": "scratch-run-registry/v1", "runs": []}, indent=2) + "\n", encoding="utf-8")
        for patcher in (patch.object(naming, "_main_root", self.root),
                        patch.object(naming, "_code_state", return_value=CODE),
                        patch.object(storage, "storage_root", return_value=self.store),
                        patch.object(runs, "ROOT", self.root)):
            patcher.start()
            self.addCleanup(patcher.stop)

    def registry(self):
        return json.loads((self.root / "runs/registry/scratch.json").read_text(encoding="utf-8"))["runs"]

    def test_tag_rule(self):
        for tag in ("base", "seed1-w128", "t8-acceptance", "parity", "a" * 40):
            self.assertEqual(naming.check_tag(tag), tag)
        for tag in ("Base", "base-v10", "try-2", "1abc", "a_b", "a--b", "-a", "a-", "", "a" * 41, "base v"):
            with self.assertRaises(ValueError, msg=tag):
                naming.check_tag(tag)

    def test_one_sequence_across_kinds_and_records(self):
        first = naming.new_directory("sandbox", "base", label="the base")
        second = naming.new_directory("sandbox", "base")
        check = naming.new_directory("step-check", "parity")
        page = naming.new_directory("page", "t8-acceptance")
        experiment = naming.new_directory("experiment", "seed1-w128", module="climate")
        self.assertEqual([first.id, second.id, check.id, page.id, experiment.id], [1, 2, 3, 4, 5])
        self.assertEqual(first.path, self.root / "out/scratch/sandbox/1-base")
        self.assertEqual(check.path, self.root / "out/scratch/step-check/3-parity")
        self.assertEqual(page.path, self.root / "out/evaluation/4-t8-acceptance")
        self.assertEqual(experiment.path, self.root / "out/evaluation/climate/5-seed1-w128")
        for named in (first, second, check, page, experiment):
            self.assertTrue(named.path.is_dir())
            self.assertEqual(list(named.path.iterdir()), [])
        records = self.registry()
        self.assertEqual([entry["id"] for entry in records], [1, 2, 3, 4, 5])
        self.assertEqual(records[0], {"id": 1, "kind": "sandbox", "tag": "base", "run": "out/scratch/sandbox/1-base",
                                      "date": records[0]["date"], "label": "the base", "tool": records[0]["tool"],
                                      "code": CODE})
        self.assertEqual(records[4]["module"], "climate")
        self.assertNotIn("module", records[0])

    def test_next_id_counts_numbered_folders_only(self):
        for folder in ("out/scratch/sandbox/7-old", "out/scratch/sandbox/v9", "out/evaluation/climate/9-x",
                       "out/evaluation/climate/20261005T141419Z-seed1-128", "out/scratch/geology/3-seed1",
                       "out/evaluation/2-page"):
            (self.root / folder).mkdir(parents=True)
        self.assertEqual(naming.new_directory("geology", "seed1").id, 10)

    def test_invalid_requests_create_nothing(self):
        with self.assertRaises(ValueError):
            naming.new_directory("sandbox", "base-v11")
        with self.assertRaises(ValueError):
            naming.new_directory("experiment", "seed1-w128")
        with self.assertRaises(ValueError):
            naming.new_directory("bench", "base")
        self.assertFalse((self.root / "out").exists())
        self.assertEqual(self.registry(), [])

    def test_refuses_a_parent_outside_the_storage_root(self):
        with tempfile.TemporaryDirectory() as other, \
                patch.object(storage, "storage_root", return_value=Path(other)):
            with self.assertRaises(SystemExit):
                naming.new_directory("sandbox", "base")
        self.assertFalse((self.root / "out/scratch/sandbox").exists())
        self.assertEqual(self.registry(), [])

    def test_find(self):
        named = naming.new_directory("sandbox", "base")
        legacy = self.root / "out/scratch/sandbox/v9"
        legacy.mkdir()
        self.assertEqual(naming.find("sandbox", "1"), named.path)
        self.assertEqual(naming.find("sandbox", "1-base"), named.path)
        self.assertEqual(naming.find("sandbox", "v9"), legacy)
        self.assertEqual(naming.find("sandbox", str(legacy)), legacy)
        with self.assertRaises(FileNotFoundError):
            naming.find("sandbox", "v8")
        with self.assertRaises(FileNotFoundError):
            naming.find("sandbox", "2")

    def test_find_step_check_of_the_earlier_layout(self):
        (self.root / "out/scratch/phase5/step4a2-tectonics").mkdir(parents=True)
        (self.root / "out/scratch/sphere/step-check/a15-cube-tectonics").mkdir(parents=True)
        self.assertEqual(naming.find("step-check", "step4a2"), self.root / "out/scratch/phase5/step4a2")
        self.assertEqual(naming.find("step-check", "a15-cube"), self.root / "out/scratch/sphere/step-check/a15-cube")
        with self.assertRaises(FileNotFoundError):
            naming.find("step-check", "a15")
        check = naming.new_directory("step-check", "parity")
        self.assertEqual(naming.find("step-check", check.name), check.path)

    def test_previous_is_the_topic_s_newest_earlier_record(self):
        first = naming.new_directory("sandbox", "base")
        naming.new_directory("sandbox", "other")
        naming.new_directory("experiment", "base", module="climate")
        third = naming.new_directory("sandbox", "base")
        self.assertIsNone(naming.previous(first))
        self.assertEqual(naming.previous(third)["id"], first.id)

    def test_registered_run_keeps_its_module_sequence(self):
        (self.root / "runs/registry/world.json").write_text(json.dumps(
            {"schema": "world-run-registry/v1", "runs": [{"id": 1, "run": "runs/world/1", "tag": "world"},
                                                          {"id": 2, "run": "runs/world/2", "tag": "other"}]}),
            encoding="utf-8")

        def run_directory(module, identifier):
            path = self.root / "runs" / module / str(identifier)
            path.mkdir(parents=True)
            return path

        with patch.object(storage, "run_directory", side_effect=run_directory):
            named = naming.new_directory("run", "world", module="world")
        self.assertEqual((named.id, named.name, named.path), (3, "3", self.root / "runs/world/3"))
        self.assertEqual(self.registry(), [])
        self.assertEqual(naming.previous(named)["id"], 1)

    def test_lock_file_stays_and_a_held_lock_times_out(self):
        lock = self.root / naming.LOCK
        lock.parent.mkdir(parents=True)
        lock.write_text("left by an earlier process", encoding="utf-8")
        self.assertEqual(naming.new_directory("sandbox", "base").id, 1)  # a leftover file does not block
        self.assertTrue(lock.exists())
        with naming._Lock(self.root), patch.object(naming, "LOCK_TIMEOUT_SECONDS", 0.3):
            started = time.monotonic()
            with self.assertRaises(TimeoutError):
                naming.new_directory("sandbox", "base")
            self.assertGreaterEqual(time.monotonic() - started, 0.25)
        self.assertEqual([entry["id"] for entry in self.registry()], [1])
        self.assertEqual(naming.new_directory("sandbox", "base").id, 2)  # released on exit

    def test_command_line_prints_the_path(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(naming.main(["new", "sandbox", "base", "--label", "base bundle"]), 0)
        self.assertEqual(output.getvalue().strip(), str(self.root / "out/scratch/sandbox/1-base"))
        self.assertEqual(self.registry()[0]["label"], "base bundle")


if __name__ == "__main__":
    unittest.main()
