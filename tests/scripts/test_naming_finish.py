"""scripts/uwtools/naming.py: finish, pinned, tool_environment and page_directory (master plan 2.4).

The main tree and the storage root are temporary directories, as in test_naming.py: nothing is written to the real
registry, to out/ or to the storage drive. Most tests stub archive.archive; the last ones run the real archive on
the temporary tree.
"""

import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from uwtools import archive, naming, paths, runs, storage  # noqa: E402

CODE = {"revision": "0" * 40, "dirty": False}


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def receipt(dropped=2, size=2_000_000):
    return {"dropped": [{"path": f"evaluation/{index}"} for index in range(dropped)], "dropped_bytes": size}


class Tree(unittest.TestCase):
    """A temporary main tree inside a temporary storage root; finish's output captured."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.store = Path(os.path.realpath(self.temporary.name))
        self.root = self.store / "main"
        write(self.root / "runs/registry/scratch.json",
              json.dumps({"schema": "scratch-run-registry/v1", "runs": []}, indent=2) + "\n")
        write(self.root / "runs/storage.json", json.dumps({"storage_root": str(self.store)}))
        write(self.root / "configs/evaluation/climate.json", "{}")
        self.worktree = self.store / "worktree"
        self.worktree.mkdir()
        for patcher in (patch.object(naming, "_main_root", self.root),
                        patch.object(naming, "_code_state", return_value=CODE),
                        patch.object(storage, "CONFIG", self.worktree / "runs/storage.json"),
                        patch.object(runs, "ROOT", self.worktree),
                        patch.object(paths, "ROOT", self.worktree)):
            patcher.start()
            self.addCleanup(patcher.stop)
        self.output = io.StringIO()
        redirect = contextlib.redirect_stdout(self.output)
        redirect.__enter__()
        self.addCleanup(redirect.__exit__, None, None, None)

    def registry(self):
        return json.loads((self.root / "runs/registry/scratch.json").read_text(encoding="utf-8"))["runs"]

    def stub(self, **behaviour):
        patcher = patch.object(archive, "archive", **behaviour)
        mock = patcher.start()
        self.addCleanup(patcher.stop)
        return mock

    def new(self, kind, tag, module=None, grid="cube", social=False):
        """A new named directory holding its kind's completion marker (a step check's says grid and social)."""
        named = naming.new_directory(kind, tag, module=module)
        marker = {"grid": grid, "social": social} if kind == "step-check" else {}
        write(named.path / naming.COMPLETE[kind], json.dumps(marker))
        return named


def only(outcomes):
    """The one outcome of a finish that saw one earlier run."""
    assert outcomes is not None and len(outcomes) == 1, outcomes
    return outcomes[0]


class FinishTests(Tree):
    def test_a_failed_run_or_a_first_run_archives_nothing(self):
        mock = self.stub(return_value=receipt())
        first = self.new("sandbox", "base")
        self.assertEqual(naming.finish(first, True), [])
        second = self.new("sandbox", "base")
        self.assertIsNone(naming.finish(second, False))
        mock.assert_not_called()
        self.assertNotIn("superseded_by", self.registry()[0])

    def test_a_passed_run_archives_the_previous_run_of_its_topic(self):
        mock = self.stub(return_value=receipt())
        first = self.new("experiment", "seed1-w128", module="climate")
        self.new("experiment", "other", module="climate")
        self.new("sandbox", "seed1-w128")
        third = self.new("experiment", "seed1-w128", module="climate")
        outcome = only(naming.finish(third, True))
        mock.assert_called_once_with(first.path, f"superseded by {third.name}", apply=True)
        self.assertEqual((outcome["action"], outcome["previous"], outcome["id"]),
                         ("archived", "out/evaluation/climate/1-seed1-w128", 1))
        self.assertIn("2 files, 2.0 MB dropped", self.output.getvalue())

    def test_every_earlier_unarchived_run_is_archived_after_a_failed_one(self):
        mock = self.stub(return_value=receipt())
        first = self.new("experiment", "seed1-w128", module="climate")
        failed = self.new("experiment", "seed1-w128", module="climate")
        self.assertIsNone(naming.finish(failed, False))
        pinned = self.new("experiment", "seed1-w128", module="climate")
        write(self.root / "configs/evaluation/world.json", json.dumps({"baseline": f"out/evaluation/climate/{pinned.name}"}))
        done = self.new("experiment", "seed1-w128", module="climate")
        write(done.path / archive.RECEIPT, "{}")
        last = self.new("experiment", "seed1-w128", module="climate")
        outcomes = naming.finish(last, True)
        self.assertEqual([(outcome["id"], outcome["action"]) for outcome in outcomes],
                         [(4, "already-archived"), (3, "pinned"), (2, "archived"), (1, "archived")])
        self.assertEqual([call.args[0] for call in mock.call_args_list], [failed.path, first.path])
        self.assertIn("1 already-archived", self.output.getvalue())

    def test_a_folder_that_is_not_a_run_is_marked_superseded(self):
        self.stub(side_effect=archive.NotARun("no evaluation.json"))
        first = self.new("step-check", "parity")
        other = self.new("sandbox", "base")
        second = self.new("step-check", "parity")
        before = self.registry()
        self.assertEqual(only(naming.finish(second, True))["action"], "superseded")
        after = self.registry()
        self.assertEqual(after[0]["superseded_by"], second.name)
        self.assertIn("superseded_utc", after[0])
        marked = dict(after[0])
        del marked["superseded_by"], marked["superseded_utc"]
        self.assertEqual([marked, *after[1:]], before)  # nothing else changed
        self.assertTrue(first.path.is_dir() and other.path.is_dir())
        # Once marked, a later finish of the same run leaves it as it is.
        self.assertEqual(only(naming.finish(second, True))["action"], "already-superseded")
        self.assertEqual(self.registry(), after)
        # The next run of the topic marks the run before it, not the first one again.
        third = self.new("step-check", "parity")
        self.assertEqual([outcome["action"] for outcome in naming.finish(third, True)],
                         ["superseded", "already-superseded"])
        self.assertEqual(self.registry()[2]["superseded_by"], third.name)
        self.assertEqual(self.registry()[0]["superseded_by"], second.name)

    def test_a_refusal_or_an_error_is_reported_and_never_raised(self):
        self.new("sandbox", "base")
        second = self.new("sandbox", "base")
        before = self.registry()
        mock = self.stub(side_effect=archive.Refused("holds tracked files"))
        outcome = only(naming.finish(second, True))
        self.assertEqual((outcome["action"], outcome["reason"]), ("refused", "holds tracked files"))
        self.assertIn("refused (holds tracked files)", self.output.getvalue())
        mock.side_effect = RuntimeError("disk gone")
        self.assertEqual(only(naming.finish(second, True))["action"], "error")
        self.assertEqual(self.registry(), before)

    def test_pinned_runs_are_left_alone(self):
        mock = self.stub(return_value=receipt())
        # Named in a tracked sandbox configuration (a path inside the run).
        base = self.new("sandbox", "base")
        write(self.root / "configs/sandbox/terrain.json", json.dumps({"base_bundle": f"out/scratch/sandbox/{base.name}"}))
        outcome = only(naming.finish(self.new("sandbox", "base"), True))
        self.assertEqual(outcome["action"], "pinned")
        self.assertIn("configs/sandbox/terrain.json", outcome["reason"])
        # An input of the new run: the step check's reference or --parent-run (an absolute path inside the run).
        reference = self.new("step-check", "parity")
        outcome = only(naming.finish(self.new("step-check", "parity"), True, [reference.path / "tectonics"]))
        self.assertEqual(outcome["action"], "pinned")
        self.assertIn("an input of", outcome["reason"])
        # The baseline of a module registry.
        write(self.root / "runs/registry/world.json", json.dumps(
            {"schema": "world-run-registry/v1", "baseline": 1, "runs": [{"id": 1, "run": "runs/world/1", "tag": "world"}]}))
        write(self.root / "runs/world/1/run.json", "{}")
        outcome = only(naming.finish(self.new("run", "world", module="world"), True))
        self.assertEqual(outcome["action"], "pinned")
        self.assertIn("baseline", outcome["reason"])
        mock.assert_not_called()
        self.assertTrue(all("superseded_by" not in entry for entry in self.registry()))

    def test_a_key_value_setting_in_a_configuration_pins_its_path(self):
        mock = self.stub(return_value=receipt())
        base = self.new("sandbox", "base")
        write(self.root / "configs/sandbox/presets.json", json.dumps(
            {"old-base": {"set": [f'base_bundle="out/scratch/sandbox/{base.name}"', "cube.enabled=false"]}}))
        outcome = only(naming.finish(self.new("sandbox", "base"), True))
        self.assertEqual(outcome["action"], "pinned")
        self.assertIn("configs/sandbox/presets.json", outcome["reason"])
        mock.assert_not_called()

    def test_a_run_awaiting_its_review_is_pinned(self):
        mock = self.stub(return_value=receipt())
        pending = self.new("experiment", "seed1-w128", module="climate")
        write(pending.path / "evaluation/evaluation.json", json.dumps({"verdict": {"checks": "pass", "review": "pending"}}))
        requested = self.new("experiment", "seed1-w128", module="climate")
        write(requested.path / "evaluation/evaluation.json", json.dumps({"verdict": {"review": "optional"}}))
        write(requested.path / "evaluation/review-request.json", "{}")
        write(requested.path / "evaluation/review.json", json.dumps({"entries": []}))
        last = self.new("experiment", "seed1-w128", module="climate")
        outcomes = naming.finish(last, True)
        self.assertEqual([outcome["action"] for outcome in outcomes], ["pinned", "pinned"])
        self.assertTrue(all("awaits a review" in outcome["reason"] for outcome in outcomes))
        mock.assert_not_called()
        # A verdict releases both: they are archived by the next passed run.
        for run in (pending, requested):
            write(run.path / "evaluation/review.json", json.dumps({"entries": [{"verdict": "approved", "reviewer": "user"}]}))
        outcomes = naming.finish(self.new("experiment", "seed1-w128", module="climate"), True)
        self.assertEqual([(outcome["id"], outcome["action"]) for outcome in outcomes],
                         [(3, "archived"), (2, "archived"), (1, "archived")])

    def test_a_step_check_of_another_grid_or_social_setting_is_left_alone(self):
        self.stub(side_effect=archive.NotARun("no evaluation.json"))
        latlon = self.new("step-check", "parity", grid="latlon")
        social = self.new("step-check", "parity", social=True)
        cube = self.new("step-check", "parity")
        outcomes = naming.finish(self.new("step-check", "parity"), True)
        self.assertEqual([(outcome["id"], outcome["action"]) for outcome in outcomes],
                         [(cube.id, "superseded"), (social.id, "other-setting"), (latlon.id, "other-setting")])
        self.assertEqual([entry.get("superseded_by") is not None for entry in self.registry()], [False, False, True, False])

    def test_a_run_without_its_completion_marker_is_never_touched(self):
        mock = self.stub(return_value=receipt())
        running = naming.new_directory("sandbox", "base")  # no config.json yet: still running
        write(running.path / "global-relief.png", "png")
        outcome = only(naming.finish(self.new("sandbox", "base"), True))
        self.assertEqual(outcome["action"], "incomplete")
        experiment = naming.new_directory("experiment", "seed1-w128", module="climate")  # no run.json
        self.assertEqual(only(naming.finish(self.new("experiment", "seed1-w128", module="climate"), True))["action"],
                         "incomplete")
        mock.assert_not_called()
        self.assertTrue(all("superseded_by" not in entry for entry in self.registry()))
        self.assertTrue(experiment.path.is_dir())

    def test_already_archived_and_missing_predecessors_are_skipped(self):
        mock = self.stub(return_value=receipt())
        first = self.new("experiment", "seed1-w128", module="climate")
        write(first.path / archive.RECEIPT, "{}")
        second = self.new("experiment", "seed1-w128", module="climate")
        self.assertEqual(only(naming.finish(second, True))["action"], "already-archived")
        # A record whose folder never existed (as after a folder went to the trash).
        runs.append_record("scratch", {"id": 3, "kind": "geology", "tag": "seed1", "run": "out/scratch/geology/3-seed1"},
                           root=self.root)
        self.assertEqual(only(naming.finish(self.new("geology", "seed1"), True))["action"], "missing")
        mock.assert_not_called()

    def test_a_registered_run_without_an_evaluation_is_not_marked(self):
        self.stub(side_effect=archive.NotARun("no evaluation.json"))
        write(self.root / "runs/registry/world.json", json.dumps(
            {"schema": "world-run-registry/v1", "runs": [{"id": 1, "run": "runs/world/1", "tag": "world"}]}))
        write(self.root / "runs/world/1/run.json", "{}")
        registry =(self.root / "runs/registry/world.json").read_bytes()
        outcome = only(naming.finish(self.new("run", "world", module="world"), True))
        self.assertEqual(outcome["action"], "not-a-run")
        self.assertEqual((self.root / "runs/registry/world.json").read_bytes(), registry)
        self.assertEqual(self.registry(), [])

    def test_tool_environment_points_below_the_directory_on_the_storage_drive(self):
        directory = self.root / "out/scratch/step-check/1-parity"
        environment = naming.tool_environment(directory, {"PATH": "x"})
        self.assertEqual(environment["PATH"], "x")
        self.assertEqual(environment["UW_CLIMATE_DIAGNOSTICS_DIR"], str(directory / "diagnostics"))
        self.assertEqual(environment["UW_DEFAULT_IMAGE_DIR"], str(directory / "images"))
        self.assertEqual(environment["UW_DEFAULT_WORLD_DIR"], str(directory / "worlds"))
        self.assertEqual(environment["UW_DEFAULT_UWS_DIR"], str(directory / "worlds"))
        with tempfile.TemporaryDirectory() as other, self.assertRaises(SystemExit):
            naming.tool_environment(Path(other) / "outputs")

    def test_page_directory_is_named_like_the_run(self):
        named = self.new("sandbox", "base")
        page = naming.page_directory(named)
        self.assertEqual(page, self.root / "out/evaluation/1-base")
        self.assertTrue(page.is_dir())
        self.assertEqual(self.new("page", "t8-acceptance").id, 2)  # the id is not drawn twice


class RealArchiveTests(Tree):
    """finish with the real archive.archive on the temporary tree (a git workspace, the storage root around it)."""

    def setUp(self):
        super().setUp()
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        for target, name, value in ((archive, "ROOT", self.root), (runs, "ARCHIVE", self.store / "runs" / "archive"),
                                    (storage, "CONFIG", self.root / "runs/storage.json")):
            patcher = patch.object(target, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        patcher = patch.object(archive, "main_root", return_value=self.root)
        patcher.start()
        self.addCleanup(patcher.stop)

    def evaluated(self, named):
        folder = named.path / "evaluation"
        write(folder / "evaluation.json", json.dumps({"schema": "uw-evaluation/v1", "profile": {"module": "climate"},
                                                      "code": CODE, "maps": []}))
        write(folder / "index.html", "<html></html>")
        write(folder / "report.md", "# report")
        write(named.path / "run.json", "{}")

    def test_the_real_archive_drops_the_pages_and_marks_a_sandbox(self):
        first = self.new("experiment", "seed1-w128", module="climate")
        self.evaluated(first)
        second = self.new("experiment", "seed1-w128", module="climate")
        self.evaluated(second)
        outcome = only(naming.finish(second, True))
        self.assertEqual(outcome["action"], "archived", self.output.getvalue())
        self.assertFalse((first.path / "evaluation/index.html").exists())
        self.assertFalse((first.path / "evaluation/report.md").exists())
        self.assertTrue((first.path / "evaluation/evaluation.json").is_file())
        self.assertTrue((first.path / "run.json").is_file())
        document = json.loads((first.path / archive.RECEIPT).read_text(encoding="utf-8"))
        self.assertEqual((document["status"], document["reason"]), ("complete", f"superseded by {second.name}"))
        self.assertTrue((second.path / "evaluation/index.html").is_file())  # the new run is untouched
        sandbox = self.new("sandbox", "base")
        write(sandbox.path / "global-relief.png", "png")
        later = self.new("sandbox", "base")
        self.assertEqual(only(naming.finish(later, True))["action"], "superseded")
        self.assertTrue((sandbox.path / "global-relief.png").is_file())
        self.assertEqual(self.registry()[2]["superseded_by"], later.name)


if __name__ == "__main__":
    unittest.main()
