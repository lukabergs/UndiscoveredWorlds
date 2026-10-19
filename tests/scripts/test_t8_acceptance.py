"""The T8 acceptance collector (scripts/evaluation/acceptance.py): missing or unreadable inputs give placeholders,
the amended engine gates, the window-log parser and the page written from the document."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

import acceptance  # noqa: E402


def entry(value, limit="limit", result="recorded"):
    return {"value": value, "limit": limit, "result": result}


def engine_fixture():
    return {"scenario": "fixture", "n64_after": {
        "bundle": "no-such-bundle",
        "table": {
            "plates": {"count": entry(10), "largest_share": entry(0.15), "plates_above_5_percent": entry(8),
                       "rms_cm_yr_sphere": entry(1.348, result="fail"), "spin_share_mean": entry(0.58),
                       "rms_speed_mean_over_run_cm_yr": entry(1.358)},
            "max_displacement_over_edge": entry(0.37, "below 1", "pass"),
            "t412": {"land": entry(0.279, "W3 0.29 +- 0.02", "pass"), "band_3000": entry(0.001, "Earth", "fail")},
            "striping_engine_surface": {"ratio_70-80": entry(0.46, "0.8-1.4", "fail")},
            "cube_imprint": {"a4_face_coordinates_all_0_80": entry(0.0279, "at most 0.028", "pass"),
                             "seam_rank_median_all": entry(0.11, "at most 0.95", "pass")},
            "time": entry(761.425, "at most 715.5 s", "fail")}}}


class CollectorTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.scratch = Path(self.directory.name) / "sphere"
        self.output = Path(self.directory.name) / "page"
        self.scratch.mkdir()

    def tearDown(self):
        self.directory.cleanup()

    def section(self, document, identifier):
        return next(section for section in document["sections"] if section["id"] == identifier)

    def test_missing_inputs_give_placeholders_and_a_page(self):
        document = acceptance.collect(self.scratch, self.output)
        for section in document["sections"]:
            self.assertEqual(section["status"], "placeholder", section["id"])
            self.assertTrue(section["missing"], section["id"])
        self.assertEqual(document["summary"]["fail"], 0)
        self.assertEqual(len(document["summary"]["placeholders"]), len(acceptance.SECTIONS))
        page = acceptance.write(document, self.scratch / "acceptance/t8-acceptance.json", self.output)
        self.assertTrue(page.is_file())
        self.assertIn(acceptance.SCHEMA, page.read_text(encoding="utf-8"))
        written = json.loads((self.scratch / "acceptance/t8-acceptance.json").read_text(encoding="utf-8"))
        self.assertEqual(written["schema"], acceptance.SCHEMA)

    def test_unreadable_input_is_a_placeholder(self):
        (self.scratch / "E6").mkdir()
        (self.scratch / "E6/acceptance.json").write_text(json.dumps({"n64_after": {}}), encoding="utf-8")
        engine = self.section(acceptance.collect(self.scratch, self.output, render=False), "engine")
        self.assertEqual(engine["status"], "placeholder")
        self.assertTrue(any("unreadable input" in item for item in engine["missing"]))

    def test_engine_gates_as_amended(self):
        (self.scratch / "E6").mkdir()
        (self.scratch / "E6/acceptance.json").write_text(json.dumps(engine_fixture()), encoding="utf-8")
        document = acceptance.collect(self.scratch, self.output, render=False)
        engine = self.section(document, "engine")
        rows = {row["id"]: row for row in engine["gates"]}
        self.assertEqual(rows["plates.rms_speed_mean_over_run"]["result"], "pass")
        self.assertEqual(rows["plates.rms_end_of_run"]["result"], "recorded")
        self.assertEqual(rows["time"]["result"], "accepted")
        self.assertEqual(rows["t412.band_3000"]["result"], "fail")
        self.assertEqual(rows["cells.area_ratio"]["result"], "missing")
        self.assertEqual(engine["status"], "partial")
        striping = {row["id"]: row for row in self.section(document, "striping")["gates"]}
        self.assertFalse(striping["engine.latlon.ratio_70-80"]["gated"])
        self.assertTrue(striping["engine.native.a4_face_coordinates_all_0_80"]["gated"])
        self.assertIn("engine: t412.band_3000", document["summary"]["failing"])
        self.assertNotIn("striping: engine.latlon.ratio_70-80", document["summary"]["failing"])

    def test_window_log_pairs(self):
        text = "\n".join([
            "interior shared cells 256 fields 14 identical 14",
            "edge01 scale 16: overlap exit 0 (6144 shared cells, 13/14 fields identical); a 6.44 s / 144 MB; "
            "ratio a 1.66 x time 1.19 x memory, b 2.62 x time 1.19 x memory; levels a",
            "== Baseline 1c002f0 for comparison",
            "corner scale 64: overlap exit 0 (102400 shared cells, 14/14 fields identical); a 105 s"])
        rows = acceptance.parse_windows("W2-W3x", text)
        self.assertEqual(len(rows), 2)
        self.assertEqual((rows[0]["scale"], rows[0]["identical"]), (4, 14))
        self.assertEqual((rows[1]["identical"], rows[1]["fields"], rows[1]["time_ratio"]), (13, 14, 2.62))
        (self.scratch / "logs").mkdir()
        (self.scratch / "logs/W2-W3b.txt").write_text(text, encoding="utf-8")
        windows = self.section(acceptance.collect(self.scratch, self.output, render=False), "windows")
        results = {row["id"]: row["result"] for row in windows["gates"]}
        self.assertEqual(results["W2-W3b.interior.4"], "pass")
        self.assertEqual(results["W2-W3b.edge01.16"], "fail")
        self.assertEqual(results["W2-W3b.edge01.16.time"], "fail")
        self.assertEqual(windows["status"], "partial")


if __name__ == "__main__":
    unittest.main()
