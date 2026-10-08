"""scripts/evaluation/regional_overlap.py on tiny synthetic regional bundles written by the test."""

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "evaluation"))

from uwtools.bundle import SCHEMA  # noqa: E402
import regional_overlap  # noqa: E402

WORLD_COLUMNS, WORLD_ROWS = 8, 4


def planet_field(scale, dtype):
    """A field on the planet-wide fine grid (kH x kW) whose every value differs."""
    rows, columns = WORLD_ROWS * scale, WORLD_COLUMNS * scale
    return (np.arange(rows * columns).reshape(rows, columns) * 0.5 + 1).astype(dtype)


def write_regional(directory, column0, columns, row0, rows, scale=2, source="a" * 64, change=None, halo=0,
                   parameters=None):
    """A regional bundle of global cells column0 .. column0 + columns - 1 (modulo W) and
    row0 .. row0 + rows - 1, cut from the planet-wide fields; change = (field, row, column)
    adds 1 to one fine cell of the bundle."""
    directory = Path(directory)
    absolute_columns = (column0 * scale + np.arange(columns * scale)) % (WORLD_COLUMNS * scale)
    absolute_rows = row0 * scale + np.arange(rows * scale)
    fields = {"geosphere.base_elevation": planet_field(scale, np.float32),
              "hydrosphere.channel_direction": planet_field(scale, np.int16) % 9,
              "geosphere.sea": (planet_field(scale, np.int64) % 2).astype(np.uint8)}
    entries = []
    for identifier, values in fields.items():
        values = values[np.ix_(absolute_rows, absolute_columns)].copy()
        if change and change[0] == identifier:
            values[change[1], change[2]] += 1
        component, name = identifier.split(".")
        relative = f"fields/{component}/{name}.npy"
        (directory / "fields" / component).mkdir(parents=True, exist_ok=True)
        np.save(directory / relative, values)
        entries.append({"id": identifier, "component": component, "layout": "grid", "dtype": values.dtype.str,
                        "shape": list(values.shape), "units": "", "file": relative, "hash": "0",
                        "written_by": "test", "description": ""})
    level = np.array(0, dtype=np.int32)
    (directory / "fields/geosphere").mkdir(parents=True, exist_ok=True)
    np.save(directory / "fields/geosphere/sea_level.npy", level)
    entries.append({"id": "geosphere.sea_level", "component": "geosphere", "layout": "scalar", "dtype": "<i4",
                    "shape": [], "units": "", "file": "fields/geosphere/sea_level.npy", "hash": "0",
                    "written_by": "test", "description": ""})
    manifest = {"schema": SCHEMA, "seed": 1, "recipe": "region",
                "grid": {"columns": columns * scale, "rows": rows * scale, "registration": "cell-centred",
                         "layout": "row-major, north row first"},
                "window": {"scale": scale, "world_columns": WORLD_COLUMNS, "world_rows": WORLD_ROWS,
                           "column0": column0, "columns": columns, "row0": row0, "rows": rows,
                           "fine_column0": column0 * scale, "fine_row0": row0 * scale, "halo": halo},
                "source": {"bundle": "world", "manifest_sha256": source, "seed": 1},
                "parameters": parameters or {"landscape.regional.elapsed_myr": 1.0},
                "fields": entries, "tables": []}
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def run(*arguments):
    output = io.StringIO()
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(io.StringIO()):
        code = regional_overlap.main([str(argument) for argument in arguments])
    return code, output.getvalue()


class RegionalOverlapTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_identical_overlap(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3)
        b = write_regional(self.root / "b", 3, 4, 1, 3)
        code, text = run(a, b)
        self.assertEqual(code, 0, text)
        result = regional_overlap.overlap(a, b)
        # Global columns 3..4 and rows 1..2 at scale 2.
        self.assertEqual((result["columns"], result["rows"]), (4, 4))
        self.assertEqual(result["fields"]["geosphere.base_elevation"],
                         {"status": "identical", "cells": 16, "differing": 0, "max_abs_difference": 0.0})
        self.assertEqual(result["fields"]["geosphere.sea_level"]["cells"], 1)

    def test_changed_cell(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3)
        # Fine cell (row 0, column 1) of B is global cell (3, 1): inside the overlap.
        b = write_regional(self.root / "b", 3, 4, 1, 3, change=("geosphere.base_elevation", 0, 1))
        code, text = run(a, b, "--json")
        self.assertEqual(code, 1)
        outcome = json.loads(text)["fields"]["geosphere.base_elevation"]
        self.assertEqual((outcome["status"], outcome["differing"], outcome["max_abs_difference"]), ("different", 1, 1.0))
        self.assertEqual(json.loads(text)["fields"]["geosphere.sea"]["status"], "identical")

    def test_change_outside_overlap_is_not_counted(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3, change=("geosphere.sea", 0, 0))
        b = write_regional(self.root / "b", 3, 4, 1, 3)
        self.assertEqual(run(a, b)[0], 0)

    def test_date_line_crossing(self):
        # A: global columns 6, 7, 0, 1 (across the date line); B: columns 0..2. Overlap: columns 0 and 1.
        a = write_regional(self.root / "a", 6, 4, 1, 2)
        b = write_regional(self.root / "b", 0, 3, 0, 4)
        code, text = run(a, b)
        self.assertEqual(code, 0, text)
        result = regional_overlap.overlap(a, b)
        self.assertEqual((result["columns"], result["rows"]), (4, 4))
        changed = write_regional(self.root / "c", 0, 3, 0, 4, change=("hydrosphere.channel_direction", 2, 3))
        code, text = run(a, changed, "--json")
        self.assertEqual(code, 1)
        self.assertEqual(json.loads(text)["fields"]["hydrosphere.channel_direction"]["differing"], 1)

    def test_refusals(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3)
        self.assertEqual(run(a, write_regional(self.root / "scale", 1, 4, 0, 3, scale=3))[0], 2)
        self.assertEqual(run(a, write_regional(self.root / "source", 1, 4, 0, 3, source="b" * 64))[0], 2)
        self.assertEqual(run(a, write_regional(self.root / "apart", 6, 1, 0, 3))[0], 2)

    def test_refuses_another_halo(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3)
        self.assertEqual(run(a, write_regional(self.root / "halo", 1, 4, 0, 3, halo=1))[0], 2)

    def test_refuses_other_parameters(self):
        a = write_regional(self.root / "a", 1, 4, 0, 3)
        other = write_regional(self.root / "parameters", 1, 4, 0, 3, parameters={"landscape.regional.elapsed_myr": 2.0})
        self.assertEqual(run(a, other)[0], 2)
        with self.assertRaisesRegex(regional_overlap.Refused, "elapsed_myr"):
            regional_overlap.overlap(a, other)


if __name__ == "__main__":
    unittest.main()
