"""The erosion check kinds of scripts/evaluation/checks.py on a tiny synthetic bundle (4 x 2)."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

from uwtools.bundle import SCHEMA, Bundle  # noqa: E402
from checks import CHECKS, land_mask  # noqa: E402

COLUMNS, ROWS = 4, 2

BUDGET = {"check": "volume-budget", "total": "geosphere.eroded_volume",
          "parts": ["geosphere.land_deposited_volume", "geosphere.marine_deposited_volume"],
          "max_relative_residual": 1e-6}
DRAINED = {"check": "drained", "direction": "hydrosphere.flow_direction", "depression": "hydrosphere.depression_id",
           "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
MONOTONE = {"check": "monotone", "direction": "hydrosphere.flow_direction", "discharge": "hydrosphere.discharge",
            "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
LAKES = {"check": "lake-above-terrain", "lake": "hydrosphere.lake_surface", "elevation": "geosphere.base_elevation"}


def write_bundle(directory, fields):
    """A minimal uw-run-bundle/v1 directory holding the given fields ({id: array})."""
    directory = Path(directory)
    entries = []
    for identifier, values in fields.items():
        values = np.asarray(values)
        component, name = identifier.split(".")
        relative = f"fields/{component}/{name}.npy"
        (directory / "fields" / component).mkdir(parents=True, exist_ok=True)
        np.save(directory / relative, values)
        entries.append({"id": identifier, "component": component, "layout": "scalar" if values.ndim == 0 else "grid",
                        "dtype": values.dtype.str, "shape": list(values.shape), "units": "", "file": relative,
                        "hash": "0", "written_by": "test", "description": ""})
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "eroded-world",
                "grid": {"columns": COLUMNS, "rows": ROWS, "registration": "cell-centred",
                         "layout": "row-major, north row first"},
                "fields": entries, "tables": []}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return Bundle.open(directory)


def routing(direction, discharge, depression=None):
    """Routing fields on a world whose column 0 is sea and columns 1..3 are land."""
    elevation = np.array([[-100.0, 100.0, 200.0, 300.0], [-100.0, 100.0, 200.0, 300.0]], dtype=np.float32)
    return {
        "geosphere.base_elevation": elevation,
        "geosphere.sea_level": np.array(0, dtype=np.int32),
        "hydrosphere.flow_direction": np.array(direction, dtype=np.uint8),
        "hydrosphere.discharge": np.array(discharge, dtype=np.float32),
        "hydrosphere.depression_id": np.zeros((ROWS, COLUMNS), dtype=np.int32) if depression is None
        else np.array(depression, dtype=np.int32),
    }


class EvaluationCheckTests(unittest.TestCase):
    def run_check(self, fields, spec):
        with tempfile.TemporaryDirectory() as directory:
            return CHECKS[spec["check"]](write_bundle(directory, fields), spec, {})["status"]

    def budget(self, total, land, marine):
        return {"geosphere.eroded_volume": np.array(total, dtype=np.float32),
                "geosphere.land_deposited_volume": np.array(land, dtype=np.float32),
                "geosphere.marine_deposited_volume": np.array(marine, dtype=np.float32)}

    def test_volume_budget(self):
        self.assertEqual(self.run_check(self.budget(10, 4, 6), BUDGET), "pass")
        self.assertEqual(self.run_check(self.budget(10, 4, 5), BUDGET), "fail")

    def test_drained(self):
        west = [[0, 7, 7, 7], [0, 7, 7, 7]]
        flow = [[0, 3, 2, 1], [0, 3, 2, 1]]
        self.assertEqual(self.run_check(routing(west, flow), DRAINED), "pass")
        undrained = [[0, 7, 0, 7], [0, 7, 7, 7]]
        self.assertEqual(self.run_check(routing(undrained, flow), DRAINED), "fail")
        in_depression = [[0, 0, 1, 0], [0, 0, 0, 0]]
        self.assertEqual(self.run_check(routing(undrained, flow, in_depression), DRAINED), "pass")

    def test_monotone(self):
        west = [[0, 7, 7, 7], [0, 7, 7, 7]]
        self.assertEqual(self.run_check(routing(west, [[0, 3, 2, 1], [0, 3, 2, 1]]), MONOTONE), "pass")
        self.assertEqual(self.run_check(routing(west, [[0, 3, 2, 1], [0, 1, 2, 1]]), MONOTONE), "fail")

    def test_lake_above_terrain(self):
        fields = routing([[0, 7, 7, 7], [0, 7, 7, 7]], [[0, 3, 2, 1], [0, 3, 2, 1]])
        lake = np.zeros((ROWS, COLUMNS), dtype=np.float32)
        lake[0, 1], lake[1, 2] = 100.0, 250.0  # at its terrain (100) and above it (200)
        fields["hydrosphere.lake_surface"] = lake
        self.assertEqual(self.run_check(fields, LAKES), "pass")
        lake[0, 3] = 299.5  # below its terrain (300)
        self.assertEqual(self.run_check(fields, LAKES), "fail")
        del fields["hydrosphere.lake_surface"]
        self.assertEqual(self.run_check(fields, LAKES), "skip")

    def test_land_rule_rounds_to_whole_metres(self):
        fields = routing([[0, 7, 7, 7], [0, 7, 7, 7]], [[0, 3, 2, 1], [0, 3, 2, 1]])
        fields["geosphere.base_elevation"][0, 1] = 0.4
        with tempfile.TemporaryDirectory() as directory:
            land = land_mask(write_bundle(directory, fields), DRAINED)
        self.assertFalse(land[0, 1])
        self.assertTrue(land[1, 1])


if __name__ == "__main__":
    unittest.main()
