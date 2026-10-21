"""scripts/evaluation/process_access.py fails on an unaudited report and on a stale file (master plan 2.6).

OUTPUT is patched to a temporary file, so src/pipeline/world_process_access.inc is never written.
"""

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "evaluation"))
import process_access  # noqa: E402


def step(measured):
    return {"label": "demo", "process": "test.unknown_process", "measured": measured, "reads": [], "writes": [],
            "table_reads": [], "table_writes": [], "drew_random": False}


def report(measured):
    return {"schema": "uw-access-report/v1", "tracking_reads": True, "recipe": "procedural-world",
            "steps": [step(measured)]}


class ProcessAccessTests(unittest.TestCase):
    def run_main(self, directory, content, *flags):
        path = Path(directory) / "report.json"
        path.write_text(json.dumps(content), encoding="utf-8")
        argv = ["process_access.py", str(path), *flags]
        with patch.object(sys, "argv", argv), self.assertRaises(SystemExit) as raised:
            process_access.main()
        return raised.exception.code

    def test_unmeasured_step_exits_non_zero_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "world_process_access.inc"
            with patch.object(process_access, "OUTPUT", output):
                code = self.run_main(directory, report(measured=False))
            self.assertNotIn(code, (0, None))
            self.assertIn("was not audited", str(code))
            self.assertFalse(output.exists())

    def test_check_exits_one_when_the_generated_file_differs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "world_process_access.inc"
            output.write_text("stale\n", encoding="utf-8")
            processes = Path(directory) / "world_processes.cpp"
            processes.write_text("// no process table\n", encoding="utf-8")
            # ROOT is the temporary directory so the script's relative paths in its messages resolve.
            with patch.object(process_access, "OUTPUT", output), \
                    patch.object(process_access, "ROOT", Path(directory)), \
                    patch.object(process_access, "PROCESSES", processes):
                code = self.run_main(directory, report(measured=True), "--check")
            self.assertEqual(code, 1)
            self.assertEqual(output.read_text(encoding="utf-8"), "stale\n")


if __name__ == "__main__":
    unittest.main()
