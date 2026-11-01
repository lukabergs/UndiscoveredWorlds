"""scripts/uwtools/naming.py: tags, ids, the scratch registry, find and previous (master plan 2.2).

The main tree and the storage root are temporary directories: nothing is written to the real
registry, to out/ or to the storage drive, and nothing is linked. Named directories are created on the
storage root at their logical path (storage.data_path); the registry keeps the logical path.
"""

import contextlib
import io
import json
import os
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
        self.store = Path(os.path.realpath(self.temporary.name))
        self.root = self.store / "main"
        (self.root / "runs/registry").mkdir(parents=True)
        (self.root / "runs/registry/scratch.json").write_text(
            json.dumps({"schema": "scratch-run-registry/v1", "runs": []}, indent=2) + "\n", encoding="utf-8")
        self.name_storage(self.store)
        (self.root / "configs/evaluation").mkdir(parents=True)
        (self.root / "configs/evaluation/climate.json").write_text("{}", encoding="utf-8")
        # A worktree: this checkout's runs/ and storage.json are not the main tree's and are never used.
        self.worktree = self.store / "worktree"
        self.worktree.mkdir()
        for patcher in (patch.object(naming, "_main_root", self.root),
                        patch.object(storage, "_main", return_value=self.root),
                        patch.object(naming, "_code_state", return_value=CODE),
                        patch.object(storage, "CONFIG", self.worktree / "runs/storage.json"),
                        patch.object(runs, "ROOT", self.worktree)):
            patcher.start()
            self.addCleanup(patcher.stop)

    def name_storage(self, store):
        """The main tree's runs/storage.json names store."""
        (self.root / "runs/storage.json").write_text(json.dumps({"storage_root": str(store)}), encoding="utf-8")

    def registry(self):
        return json.loads((self.root / "runs/registry/scratch.json").read_text(encoding="utf-8"))["runs"]

    def test_tag_rule(self):
        for tag in ("base", "seed1-w128", "t8-acceptance", "parity", "step-1-try", "a" * 40):
            self.assertEqual(naming.check_tag(tag), tag)
        for tag in ("Base", "base-v10", "try-2", "step-1-try-2", "1abc", "a_b", "a--b", "-a", "a-", "", "a" * 41, "base v"):
            with self.assertRaises(ValueError, msg=tag):
                naming.check_tag(tag)

    def test_one_sequence_across_kinds_and_records(self):
        first = naming.new_directory("sandbox", "base", label="the base")
        second = naming.new_directory("sandbox", "base")
        check = naming.new_directory("step-check", "parity")
        page = naming.new_directory("page", "t8-acceptance")
        experiment = naming.new_directory("experiment", "seed1-w128", module="climate")
        self.assertEqual([first.id, second.id, check.id, page.id, experiment.id], [1, 2, 3, 4, 5])
        self.assertEqual(first.path, self.store / "out/scratch/sandbox/1-base")
        self.assertEqual(check.path, self.store / "out/scratch/step-check/3-parity")
        self.assertEqual(page.path, self.store / "out/evaluation/4-t8-acceptance")
        self.assertEqual(experiment.path, self.store / "out/evaluation/climate/5-seed1-w128")
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
                       "out/evaluation/2-page", "out/evaluation/notamodule/50-x"):
            (self.store / folder).mkdir(parents=True)
        self.assertEqual(naming.new_directory("geology", "seed1").id, 10)  # only configs/evaluation modules count

    def test_invalid_requests_create_nothing(self):
        with self.assertRaises(ValueError):
            naming.new_directory("sandbox", "base-v11")
        with self.assertRaises(ValueError):
            naming.new_directory("experiment", "seed1-w128")
        with self.assertRaises(ValueError):
            naming.new_directory("bench", "base")
        self.assertFalse((self.store / "out").exists())
        self.assertEqual(self.registry(), [])

    def test_directories_go_to_the_storage_root_never_the_main_tree(self):
        with tempfile.TemporaryDirectory() as other:
            other = Path(os.path.realpath(other))
            self.name_storage(other)
            named = naming.new_directory("sandbox", "base")
            run = naming.new_directory("run", "world", module="world")
            self.assertEqual(named.path, other / "out/scratch/sandbox/1-base")
            self.assertEqual(run.path, other / "runs/world/1")
            self.assertTrue(named.path.is_dir() and run.path.is_dir())
        self.assertFalse((self.root / "out").exists())
        self.assertFalse((self.root / "runs/world").exists())
        self.assertEqual(self.registry()[0]["run"], "out/scratch/sandbox/1-base")  # the logical path

    def test_refuses_a_parent_outside_the_storage_root(self):
        # A parent whose logical path is not a storage place resolves in the main tree: refused, nothing created
        # (a storage root beside the main tree: in setUp the main tree lies inside the temporary storage root).
        with tempfile.TemporaryDirectory() as other, patch.dict(naming.PARENTS, {"sandbox": "docs/sandbox"}):
            self.name_storage(Path(os.path.realpath(other)))
            with self.assertRaises(SystemExit):
                naming.new_directory("sandbox", "base")
        self.assertFalse((self.root / "docs").exists())
        self.assertEqual(self.registry(), [])

    def test_storage_root_is_read_from_the_main_tree(self):
        # A main tree without runs/storage.json; the worktree's copy names a valid store but is not read.
        bare = self.store / "bare"
        (bare / "runs/registry").mkdir(parents=True)
        (self.worktree / "runs").mkdir()
        (self.worktree / "runs/storage.json").write_text(json.dumps({"storage_root": str(self.store)}), encoding="utf-8")
        with patch.object(naming, "_main_root", bare), patch.object(storage, "_main", return_value=bare):
            with self.assertRaises(SystemExit):
                naming.new_directory("sandbox", "base")
        self.assertFalse((bare / "out").exists())
        self.assertFalse((self.store / "out").exists())

    def test_find(self):
        named = naming.new_directory("sandbox", "base")
        legacy = self.store / "out/scratch/sandbox/v9"
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
        (self.store / "out/scratch/phase5/step4a2-tectonics").mkdir(parents=True)
        (self.store / "out/scratch/sphere/step-check/a15-cube-tectonics").mkdir(parents=True)
        self.assertEqual(naming.find("step-check", "step4a2"), self.store / "out/scratch/phase5/step4a2")
        self.assertEqual(naming.find("step-check", "a15-cube"), self.store / "out/scratch/sphere/step-check/a15-cube")
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

    def test_registered_run_uses_the_main_tree_from_a_worktree(self):
        (self.root / "runs/registry/world.json").write_text(json.dumps(
            {"schema": "world-run-registry/v1", "runs": [{"id": 1, "run": "runs/world/1", "tag": "world"},
                                                          {"id": 2, "run": "runs/world/2", "tag": "other"}]}),
            encoding="utf-8")
        (self.store / "runs/world/4").mkdir(parents=True)  # a directory without a record still takes its id
        named = naming.new_directory("run", "world", module="world")
        self.assertEqual((named.id, named.name, named.path), (5, "5", self.store / "runs/world/5"))
        self.assertTrue(named.path.is_dir())
        self.assertEqual(list(self.worktree.iterdir()), [])  # nothing in the worktree
        self.assertEqual(self.registry(), [])
        self.assertEqual(naming.previous(named)["id"], 1)
        self.assertEqual(naming.find("run", "2", module="world"), self.store / "runs/world/2")

    def test_lock_file_stays_and_a_held_lock_times_out(self):
        lock = self.store / naming.LOCK
        lock.parent.mkdir(parents=True)
        lock.write_text("left by an earlier process", encoding="utf-8")
        self.assertEqual(naming.new_directory("sandbox", "base").id, 1)  # a leftover file does not block
        self.assertTrue(lock.exists())
        with naming._Lock(), patch.object(naming, "LOCK_TIMEOUT_SECONDS", 0.3):
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
        self.assertEqual(output.getvalue().strip(), str(self.store / "out/scratch/sandbox/1-base"))
        self.assertEqual(self.registry()[0]["label"], "base bundle")


if __name__ == "__main__":
    unittest.main()
