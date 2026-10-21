"""archive.py drops only regenerable map PNGs (R1) and pages (R2), with a receipt (master plan 2.3).

Everything happens in a temporary folder: a git workspace holding a fake storage root, with the repository root,
the main tree, the storage root and runs/archive patched to it. Nothing outside the temporary folder is read or
written.
"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from uwtools import archive, runs, storage, trash  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402

DROPPED = {
    "evaluation/maps/terrain.png": "R1",
    "evaluation/maps/relief.png": "R1",
    "evaluation/index.html": "R2",
    "evaluation/report.md": "R2",
}
KEPT = [
    "run.json", "log.txt", "bundle/manifest.json", "bundle/fields/elevation.npy", "evaluation/evaluation.json",
    "evaluation/review-request.json", "evaluation/maps/palette.json", "evaluation/maps/unlisted.png",
]
REVISION = "0123456789abcdef0123456789abcdef01234567"


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content, encoding="utf-8")


def make_run(run, manifest=True):
    write(run / "run.json", "{}")
    write(run / "log.txt", "log")
    write(run / "bundle/fields/elevation.npy", b"\x93NUMPY" + bytes(64))
    if manifest:
        write(run / "bundle/manifest.json", '{"schema": "uw-bundle/v2"}')
    evaluation = {
        "schema": "uw-evaluation/v1",
        "profile": {"module": "regional"},
        "run": {"bundle": (run / "bundle").as_posix(),
                "manifest_sha256": sha256(run / "bundle/manifest.json") if manifest else "0" * 64},
        "baseline": None,
        "code": {"revision": REVISION, "dirty": False},
        "maps": [{"map": "terrain", "kind": "model", "path": "maps/terrain.png"},
                 {"map": "relief", "kind": "model", "path": "maps/relief.png"}],
    }
    write(run / "evaluation/evaluation.json", json.dumps(evaluation))
    write(run / "evaluation/review-request.json", '{"schema": "uw-review-request/v1"}')
    write(run / "evaluation/index.html", "<html></html>")
    write(run / "evaluation/report.md", "# report")
    write(run / "evaluation/maps/terrain.png", b"\x89PNG terrain")
    write(run / "evaluation/maps/relief.png", b"\x89PNG relief")
    write(run / "evaluation/maps/unlisted.png", b"\x89PNG unlisted")
    write(run / "evaluation/maps/palette.json", "{}")
    return {path: sha256(run / path) for path in [*DROPPED, *KEPT] if (run / path).is_file()}


def edit_evaluation(run, change):
    path = run / "evaluation/evaluation.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    change(document)
    path.write_text(json.dumps(document), encoding="utf-8")


def classes(receipt):
    return {item["path"]: item["class"] for item in receipt["dropped"]}


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.workspace = Path(os.path.realpath(self.temporary.name)) / "workspace"
        self.store = self.workspace / "store"
        self.store.mkdir(parents=True)
        subprocess.run(["git", "init", "-q", str(self.workspace)], check=True)
        for target, name, value in (
                (archive, "ROOT", self.workspace), (trash, "ROOT", self.workspace),
                (runs, "ARCHIVE", self.store / "runs" / "archive")):
            patcher = patch.object(target, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        for target, name, value in ((storage, "storage_root", self.store), (archive, "main_root", self.workspace)):
            patcher = patch.object(target, name, return_value=value)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.run_dir = self.store / "runs" / "regional" / "44"
        self.digests = make_run(self.run_dir)

    def receipt(self):
        return json.loads((self.run_dir / archive.RECEIPT).read_text(encoding="utf-8"))

    def copy(self):
        return self.store / "receipts" / "archive" / "runs__regional__44.json"

    def assert_untouched(self, run=None):
        run = run or self.run_dir
        for path in [*DROPPED, *KEPT]:
            self.assertTrue((run / path).is_file(), path)
        self.assertFalse((run / archive.RECEIPT).exists())

    def assert_refused(self, path):
        with self.assertRaises(archive.Refused):
            archive.archive(path, "x", apply=True)

    def test_plan_lists_both_classes_and_drops_nothing(self):
        receipt = archive.archive(self.run_dir, "", apply=False)
        self.assertEqual(classes(receipt), DROPPED)
        self.assertEqual(receipt["status"], "plan")
        self.assertEqual(receipt["kept_files"], len(KEPT))
        self.assert_untouched()
        self.assertFalse(self.copy().exists())

    def test_unlisted_png_is_kept(self):
        receipt = archive.archive(self.run_dir, "", apply=False)
        self.assertNotIn("evaluation/maps/unlisted.png", classes(receipt))

    def test_dropped_items_record_the_code_and_a_bare_command(self):
        receipt = archive.archive(self.run_dir, "", apply=False)
        for item in receipt["dropped"]:
            self.assertEqual(item["code"], {"revision": REVISION, "dirty": False}, item["path"])
            self.assertNotIn(REVISION, item["regenerate"], item["path"])
            self.assertTrue(item["regenerate"].startswith("uv run --offline python scripts/"), item["path"])

    def test_apply_drops_exactly_the_regenerable_files(self):
        receipt = archive.archive(self.run_dir, "superseded", apply=True)
        for path in DROPPED:
            self.assertFalse((self.run_dir / path).exists(), path)
        for path in KEPT:
            self.assertEqual(sha256(self.run_dir / path), self.digests[path], path)
        stored = self.receipt()
        self.assertEqual(stored, receipt)
        self.assertEqual(stored["schema"], "uw-archive-receipt/v1")
        self.assertEqual(stored["status"], "complete")
        self.assertEqual(stored["reason"], "superseded")
        self.assertEqual({item["path"]: item["sha256"] for item in stored["dropped"]},
                         {path: self.digests[path] for path in DROPPED})
        self.assertEqual(stored["dropped_bytes"], sum(item["bytes"] for item in stored["dropped"]))
        self.assertEqual(json.loads(self.copy().read_text(encoding="utf-8")), stored)
        # Applying again finds nothing new and returns the same receipt.
        self.assertEqual(archive.archive(self.run_dir, "again", apply=True), stored)

    def test_run_without_bundle_manifest_keeps_its_pngs(self):
        run = self.store / "runs" / "regional" / "45"
        make_run(run, manifest=False)
        receipt = archive.archive(run, "", apply=False)
        self.assertEqual(classes(receipt), {path: kind for path, kind in DROPPED.items() if kind == "R2"})

    def test_an_animation_keeps_every_map(self):
        write(self.run_dir / "evaluation/maps/animation.png", b"\x89PNG frames")
        edit_evaluation(self.run_dir, lambda document: document["maps"].append(
            {"map": "terrain", "kind": "animation", "path": "maps/animation.png"}))
        receipt = archive.archive(self.run_dir, "", apply=False)
        self.assertEqual(classes(receipt), {path: kind for path, kind in DROPPED.items() if kind == "R2"})

    def test_a_parent_map_needs_the_matching_baseline_bundle(self):
        write(self.run_dir / "evaluation/maps/terrain-parent.png", b"\x89PNG parent")
        edit_evaluation(self.run_dir, lambda document: document["maps"].append(
            {"map": "terrain", "kind": "parent", "path": "maps/terrain-parent.png"}))
        self.assertEqual(classes(archive.archive(self.run_dir, "", apply=False)), DROPPED)
        parent = self.store / "runs" / "regional" / "43"
        make_run(parent)
        baseline = {"bundle": (parent / "bundle").as_posix(), "manifest_sha256": "f" * 64}
        edit_evaluation(self.run_dir, lambda document: document.update(baseline=baseline))
        self.assertEqual(classes(archive.archive(self.run_dir, "", apply=False)), DROPPED)
        baseline["manifest_sha256"] = sha256(parent / "bundle/manifest.json")
        edit_evaluation(self.run_dir, lambda document: document.update(baseline=baseline))
        self.assertEqual(classes(archive.archive(self.run_dir, "", apply=False)),
                         dict(DROPPED, **{"evaluation/maps/terrain-parent.png": "R1"}))

    def test_a_judged_run_keeps_its_maps(self):
        write(self.run_dir / "evaluation/review.json",
              json.dumps({"schema": "uw-review/v1", "entries": [{"verdict": "approved", "reviewer": "user"}]}))
        receipt = archive.archive(self.run_dir, "", apply=False)
        self.assertEqual(classes(receipt), {path: kind for path, kind in DROPPED.items() if kind == "R2"})

    def test_archived_root_with_an_evaluation_is_an_evaluation_folder(self):
        source = self.run_dir / "evaluation"
        for schema, expected in (("uw-evaluation/v1", {"maps/terrain.png": "R1", "maps/relief.png": "R1",
                                                       "index.html": "R2", "report.md": "R2"}),
                                 ("other/v1", {})):
            folder = self.store / "out" / "evaluation" / "regional" / f"3-{schema.split('/')[0]}"
            for name in ("index.html", "report.md", "maps/terrain.png", "maps/relief.png", "maps/unlisted.png"):
                write(folder / name, (source / name).read_bytes())
            document = json.loads((source / "evaluation.json").read_text(encoding="utf-8"))
            write(folder / "evaluation.json", json.dumps(dict(document, schema=schema)))
            receipt = archive.archive(folder, "", apply=False)
            self.assertEqual(classes(receipt), expected, schema)

    def test_a_named_evaluation_folder_needs_the_schema(self):
        edit_evaluation(self.run_dir, lambda document: document.update(schema="other/v1"))
        self.assertEqual(classes(archive.archive(self.run_dir, "", apply=False)), {})

    def test_refuses_anything_but_one_run_level_folder(self):
        sandbox = self.store / "out" / "scratch" / "sandbox" / "7-base"
        write(sandbox / "terrain.png", b"\x89PNG sandbox")
        for path in (sandbox, self.store / "runs" / "regional"):  # no evaluation.json of its own
            with self.assertRaises(archive.NotARun):
                archive.archive(path, "x", apply=True)
        self.assertTrue((sandbox / "terrain.png").is_file())
        make_run(self.run_dir / "inner")  # a second run inside the first
        with self.assertRaises(archive.Refused) as raised:
            archive.archive(self.run_dir, "x", apply=True)
        self.assertNotIsInstance(raised.exception, archive.NotARun)
        self.assert_untouched()

    def test_a_comparison_nested_in_an_evaluation_folder_is_kept_whole(self):
        comparison = self.run_dir / "evaluation" / "comparison-8"
        make_run(comparison)  # its evaluation.json, pages and maps lie below evaluation/comparison-8/
        nested = {path: sha256(comparison / path) for path in [*DROPPED, *KEPT]}
        receipt = archive.archive(self.run_dir, "superseded", apply=True)
        self.assertEqual(classes(receipt), DROPPED)
        for path, digest in nested.items():
            self.assertEqual(sha256(comparison / path), digest, path)

    def test_refuses_a_link_inside(self):
        import _winapi
        target = self.workspace.parent / "elsewhere"
        target.mkdir()
        _winapi.CreateJunction(str(target), str(self.run_dir / "linked"))
        self.assert_refused(self.run_dir)
        self.assert_untouched()

    def test_refuses_a_junction_target(self):
        import _winapi
        target = self.store / "out" / "scratch-target"
        make_run(target)
        (self.workspace / "out").mkdir()
        _winapi.CreateJunction(str(target), str(self.workspace / "out" / "scratch"))
        self.assert_refused(target)
        self.assert_refused(self.workspace / "out" / "scratch")
        self.assert_untouched(target)

    def test_refuses_a_path_outside_the_storage_root(self):
        run = self.workspace / "out" / "scratch" / "run"
        make_run(run)
        self.assert_refused(run)
        self.assert_untouched(run)

    def test_refuses_tracked_files(self):
        subprocess.run(["git", "-C", str(self.workspace), "add", "--", (self.run_dir / "run.json").as_posix()],
                       check=True)
        self.assert_refused(self.run_dir)
        self.assert_untouched()

    def test_refuses_files_tracked_at_the_main_tree_path(self):
        # The storage path maps back to runs/regional/44 in the main tree, which holds a tracked file.
        write(self.workspace / "runs/regional/44/run.json", "{}")
        subprocess.run(["git", "-C", str(self.workspace), "add", "--", "runs/regional/44/run.json"], check=True)
        self.assert_refused(self.run_dir)
        self.assert_untouched()

    def test_refuses_runs_archive(self):
        run = self.store / "runs" / "archive" / "old"
        make_run(run)
        self.assert_refused(run)
        self.assert_untouched(run)

    def test_refuses_a_protected_place(self):
        run = self.store / "receipts" / "run"
        make_run(run)
        self.assert_refused(run)
        self.assert_untouched(run)
        self.assert_refused(self.store / "runs")
        self.assert_untouched()

    def test_refuses_a_corrupt_receipt_and_an_unavailable_storage_root(self):
        write(self.run_dir / archive.RECEIPT, "{not json")
        self.assert_refused(self.run_dir)
        for path in [*DROPPED, *KEPT]:
            self.assertTrue((self.run_dir / path).is_file(), path)
        with patch.object(storage, "storage_root", side_effect=OSError("drive gone")):
            with self.assertRaises(archive.Refused):
                archive.archive(self.run_dir, "", apply=False)

    def test_several_paths_are_all_checked_before_any_drop(self):
        other = self.store / "runs" / "regional" / "45"
        make_run(other)
        write(other / archive.RECEIPT, json.dumps({"status": "unknown", "dropped": []}))
        with self.assertRaises(archive.Refused):
            archive.archive_all([self.run_dir, other], "superseded")
        self.assert_untouched()

    def test_rerun_finishes_a_stopped_job(self):
        removed = []

        def interrupted(path):
            if removed:
                raise KeyboardInterrupt
            os.unlink(path)
            removed.append(path)

        with patch.object(archive, "unlink", side_effect=interrupted):
            with self.assertRaises(KeyboardInterrupt):
                archive.archive(self.run_dir, "superseded", apply=True)
        self.assertEqual(self.receipt()["status"], "dropping")
        self.assertEqual(sum((self.run_dir / path).exists() for path in DROPPED), len(DROPPED) - 1)
        receipt = archive.archive(self.run_dir, "ignored", apply=True)
        self.assertEqual(receipt["status"], "complete")
        self.assertEqual(receipt["reason"], "superseded")
        self.assertIn("resumed_utc", receipt)
        for path in DROPPED:
            self.assertFalse((self.run_dir / path).exists(), path)
        for path in KEPT:
            self.assertTrue((self.run_dir / path).is_file(), path)

    def test_a_file_changed_since_planned_is_refused(self):
        original = archive.write_receipt

        def changing(receipt, real, root):
            original(receipt, real, root)
            if receipt["status"] == "dropping":
                (self.run_dir / "evaluation/maps/relief.png").write_bytes(b"\x89PNG changed")

        with patch.object(archive, "write_receipt", side_effect=changing):
            with self.assertRaises(archive.Refused):
                archive.archive(self.run_dir, "superseded", apply=True)
        for path in DROPPED:
            self.assertTrue((self.run_dir / path).is_file(), path)
        self.assertEqual(self.receipt()["status"], "dropping")
        self.assert_refused(self.run_dir)
        self.assertTrue((self.run_dir / "evaluation/maps/terrain.png").is_file())


if __name__ == "__main__":
    unittest.main()
