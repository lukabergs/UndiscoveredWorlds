"""storage.route refuses a non-empty storage destination that has no receipt (master plan 2.6).

Only the refusal path is exercised: the workspace and the storage root are temporary
directories, nothing is moved, linked or deleted.
"""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from uwtools import storage  # noqa: E402


class RouteRefusalTests(unittest.TestCase):
    def test_non_empty_destination_without_receipt_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "workspace"
            store = Path(directory) / "store"
            source = workspace / "runs" / "demo"
            source.mkdir(parents=True)
            (source / "result.txt").write_text("source", encoding="utf-8")
            destination = store / "runs" / "demo"
            destination.mkdir(parents=True)
            (destination / "other.txt").write_text("stored", encoding="utf-8")

            with patch.object(storage, "ROOT", workspace), \
                    patch.object(storage, "CONFIG", workspace / "runs" / "storage.json"), \
                    patch.object(storage, "storage_root", return_value=store):
                with self.assertRaises(FileExistsError):
                    storage.route(source)

            self.assertFalse(storage.junction_target(source))
            self.assertTrue(source.is_dir())
            self.assertEqual([item.name for item in source.iterdir()], ["result.txt"])
            self.assertEqual((source / "result.txt").read_text(encoding="utf-8"), "source")
            self.assertEqual([item.name for item in destination.iterdir()], ["other.txt"])
            self.assertFalse((store / "receipts").exists())
            self.assertFalse(source.with_name("demo.storage-migration-backup").exists())


if __name__ == "__main__":
    unittest.main()
