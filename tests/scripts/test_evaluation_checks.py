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
from checks import CHECKS, land_mask, ocean_mask  # noqa: E402

COLUMNS, ROWS = 4, 2

BUDGET = {"check": "volume-budget", "total": "geosphere.eroded_volume",
          "parts": ["geosphere.land_deposited_volume", "geosphere.marine_deposited_volume"],
          "max_relative_residual": 1e-6}
DRAINED = {"check": "drained", "direction": "hydrosphere.flow_direction", "depression": "hydrosphere.depression_id",
           "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
MONOTONE = {"check": "monotone", "direction": "hydrosphere.flow_direction", "discharge": "hydrosphere.discharge",
            "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
LAKES = {"check": "lake-above-terrain", "lake": "hydrosphere.lake_surface", "elevation": "geosphere.base_elevation"}


def write_bundle(directory, fields, columns=COLUMNS, rows=ROWS, parameters=None):
    """A minimal uw-run-bundle/v1 directory holding the given fields ({id: array}) and, when given,
    the manifest's parameters ({"module.group.key": value})."""
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
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "procedural-world",
                "grid": {"columns": columns, "rows": rows, "registration": "cell-centred",
                         "layout": "row-major, north row first"},
                "fields": entries, "tables": []}
    if parameters is not None:
        manifest["parameters"] = parameters
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

    def test_enclosed_below_sea_cell_is_land(self):
        """Erosion plan H3: on 8 x 4 cells an ocean in column 0 (area 2 of the planet's 16) and one
        enclosed below-sea cell (row 2, column 4; area 0.71, 4.4 percent). Under a 10 percent
        threshold the cell is land, so without a direction or a depression it is undrained; in a
        depression it passes; under the default 0.2 percent it is sea."""
        columns, rows = 8, 4
        elevation = np.full((rows, columns), 100.0, dtype=np.float32)
        elevation[:, 0] = -100.0
        elevation[2, 4] = -50.0
        direction = np.full((rows, columns), 7, dtype=np.uint8)
        direction[:, 0] = 0
        direction[2, 4] = 0
        fields = {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(0, dtype=np.int32),
                  "hydrosphere.flow_direction": direction,
                  "hydrosphere.depression_id": np.zeros((rows, columns), dtype=np.int32)}

        def status(parameters):
            with tempfile.TemporaryDirectory() as directory:
                bundle = write_bundle(directory, fields, columns, rows, parameters)
                return CHECKS["drained"](bundle, DRAINED, {})["status"], ocean_mask(bundle, DRAINED)

        threshold = {"landscape.routing.ocean_minimum_area_fraction": 0.1}
        outcome, ocean = status(threshold)
        self.assertEqual(outcome, "fail")
        self.assertTrue(ocean[:, 0].all())
        self.assertFalse(ocean[2, 4])
        self.assertEqual(status(None)[0], "pass")
        fields["hydrosphere.depression_id"][2, 4] = 1
        self.assertEqual(status(threshold)[0], "pass")

    def test_ocean_is_largest_by_area(self):
        """A polar row of 8 cells (area 2.34) loses to 4 cells near the equator (2.83) when only the
        largest region is sea."""
        columns, rows = 8, 4
        elevation = np.full((rows, columns), 100.0, dtype=np.float32)
        elevation[0, :] = -100.0
        elevation[2, 2:6] = -100.0
        fields = {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(0, dtype=np.int32)}
        with tempfile.TemporaryDirectory() as directory:
            ocean = ocean_mask(write_bundle(directory, fields, columns, rows,
                                            {"landscape.routing.ocean_minimum_area_fraction": 1.0}), DRAINED)
        self.assertFalse(ocean[0, :].any())
        self.assertTrue(ocean[2, 2:6].all())


if __name__ == "__main__":
    unittest.main()
