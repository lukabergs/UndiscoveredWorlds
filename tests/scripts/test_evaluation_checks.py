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
from checks import CHECKS, box_average, land_mask, ocean_mask, run_checks  # noqa: E402

COLUMNS, ROWS = 4, 2

BUDGET = {"check": "volume-budget", "total": "geosphere.eroded_volume",
          "parts": ["geosphere.land_deposited_volume", "geosphere.marine_deposited_volume"],
          "max_relative_residual": 1e-6}
DRAINED = {"check": "drained", "direction": "hydrosphere.flow_direction", "depression": "hydrosphere.depression_id",
           "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
MONOTONE = {"check": "monotone", "direction": "hydrosphere.flow_direction", "discharge": "hydrosphere.discharge",
            "elevation": "geosphere.base_elevation", "sea_level": "geosphere.sea_level"}
LAKES = {"check": "lake-above-terrain", "lake": "hydrosphere.lake_surface", "elevation": "geosphere.base_elevation"}


def write_bundle(directory, fields, columns=COLUMNS, rows=ROWS, parameters=None, blocks=None):
    """A minimal uw-run-bundle/v1 directory holding the given fields ({id: array}) and, when given,
    the manifest's parameters ({"module.group.key": value}) and further manifest blocks."""
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
    manifest.update(blocks or {})
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


BOUND = 0.001953125  # half a float step at 32 to 64 km, run 9's meanbound_m
WORLD_COLUMNS, WORLD_ROWS = 8, 4


def level_entries(count, errors=None, overall=None):
    """diagnostics.levels of a hierarchy of `count` levels, every bound BOUND; parent-mean errors (default
    the bound itself, a tie) and global-mean errors (default the parent-mean error) per level."""
    errors = errors or [BOUND] * count
    overall = overall or errors
    return [{"level": level + 1, "cellsize_km": 39.1 / 4 ** level, "maxparentmeanerror_m": errors[level],
             "meanbound_m": BOUND, "globalmeanerror_m": overall[level], "seconds": 1.0} for level in range(count)]


def hierarchy_blocks(level, levels=None, column0=2, columns=2, row0=1, rows=1, factor=4, digest="a" * 64):
    """Manifest blocks of a regional hierarchy bundle at `level` (scale factor^level) of a window of global
    cells on an 8 x 4 world; levels None leaves out diagnostics.levels (a single pass)."""
    scale = factor ** level
    window = {"scale": scale, "world_columns": WORLD_COLUMNS, "world_rows": WORLD_ROWS, "column0": column0, "columns": columns,
              "row0": row0, "rows": rows, "fine_column0": column0 * scale, "fine_row0": row0 * scale, "halo": 2,
              "level": level, "factor": factor}
    blocks = {"window": window, "source": {"bundle": "world", "manifest_sha256": digest, "seed": 1}, "diagnostics": {}}
    if levels is not None:
        blocks["diagnostics"]["levels"] = levels
    return blocks


def regional_fields(elevation, sea=None, channel=None):
    elevation = np.asarray(elevation, dtype=np.float32)
    return {"geosphere.base_elevation": elevation, "geosphere.sea_level": np.array(1000, dtype=np.int32),
            "geosphere.sea": np.zeros(elevation.shape, dtype=np.uint8) if sea is None else np.asarray(sea, dtype=np.uint8),
            "hydrosphere.channel_direction": np.zeros(elevation.shape, dtype=np.int16) if channel is None
            else np.asarray(channel, dtype=np.int16)}


def fine_elevation(global_columns, global_rows, scale, low=30000.0, high=40000.0, seed=1):
    """Elevation of `scale` fine cells per global cell, by default 30 to 40 km (where BOUND is at most half
    a float step, so a float32 box average lies within BOUND of the float64 one)."""
    generator = np.random.default_rng(seed)
    return generator.uniform(low, high, (global_rows * scale, global_columns * scale)).astype(np.float32)


class HierarchyCheckTests(unittest.TestCase):
    """The refinement hierarchy checks (multiscale plan S2.5) on synthetic manifests and arrays."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.count = 0

    def tearDown(self):
        self.directory.cleanup()

    def bundle(self, fields, blocks):
        self.count += 1
        rows, columns = fields["geosphere.base_elevation"].shape
        return write_bundle(Path(self.directory.name) / str(self.count), fields, columns, rows, blocks=blocks)

    def check(self, name, bundle, coarser=None):
        return CHECKS[name](bundle, {"check": name}, {"coarser": coarser})

    def test_parent_mean_exact(self):
        fields = regional_fields(fine_elevation(2, 1, 16))
        tie = self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2, level_entries(2))))
        self.assertEqual(tie["status"], "pass", tie)
        above = self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2, level_entries(2, [BOUND, BOUND * 1.001]))))
        self.assertEqual(above["status"], "fail")
        self.assertIn("level 2: parent-mean error", above["details"][0])
        # Level 2's global-mean error may reach the sum of both bounds, not more.
        summed = level_entries(2, [BOUND, BOUND], [BOUND, 2 * BOUND])
        self.assertEqual(self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2, summed)))["status"], "pass")
        beyond = level_entries(2, [BOUND, BOUND], [BOUND, 2.01 * BOUND])
        self.assertEqual(self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2, beyond)))["status"], "fail")
        single = self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2)))
        self.assertEqual(single["status"], "skip")
        self.assertTrue(single["summary"].startswith("not applicable"))

    def test_border_step(self):
        """Each level's borderstep99_m against p99_max_m (default 50 m); a missing value fails; a single pass
        is not applicable; the profile's 'report' turns a miss into 'reported'."""
        fields = regional_fields(fine_elevation(2, 1, 16))
        levels = level_entries(2)
        for entry, p99 in zip(levels, (12.0, 49.5)):
            entry.update({"borderstep99_m": p99, "maxborderstep_m": 4 * p99, "maxlambdabefore_m": 3000.0,
                          "maxlambda_m": 40.0, "smoothpasses": 6, "forcedparents": 0})
        passing = self.check("border-step", self.bundle(fields, hierarchy_blocks(2, levels)))
        self.assertEqual(passing["status"], "pass", passing)
        self.assertIn("level 2 49.5 m (max 198 m)", passing["summary"])
        self.assertIn("3000 m before 6 smoothing pass(es), 40 m after", passing["details"][1])
        levels[1]["borderstep99_m"] = 51.0
        run = self.bundle(fields, hierarchy_blocks(2, levels))
        above = self.check("border-step", run)
        self.assertEqual(above["status"], "fail")
        self.assertIn("level 2: border step p99 51 m above 50 m", above["details"][0])
        self.assertEqual(CHECKS["border-step"](run, {"check": "border-step", "p99_max_m": 60}, {})["status"], "pass")
        reported = run_checks(run, [{"check": "border-step", "p99_max_m": 50, "report": True}], {})[0]
        self.assertEqual(reported["status"], "reported")
        levels[0]["borderstep99_m"] = None
        self.assertEqual(self.check("border-step", self.bundle(fields, hierarchy_blocks(2, levels)))["status"], "fail")
        self.assertEqual(self.check("border-step", self.bundle(fields, hierarchy_blocks(2)))["status"], "skip")

    def test_box_average_factors_4_and_16(self):
        """Blocks of 4 and of 16 cells per side holding a block value plus a zero-mean pattern average to the
        block value; the 16-cell average is the 4-cell average averaged again."""
        blocks = np.arange(2.0, dtype=np.float64).reshape(1, 2) * 1000.0 + 31000.0
        pattern = np.tile(np.array([[1.0, -1.0], [-3.0, 3.0]]), (8, 8))
        values = (np.kron(blocks, np.ones((16, 16))) + np.tile(pattern, (1, 2))).astype(np.float32)
        np.testing.assert_array_equal(box_average(values, 16), blocks)
        np.testing.assert_array_equal(box_average(values, 4), np.kron(blocks, np.ones((4, 4))))
        np.testing.assert_array_equal(box_average(box_average(values, 4), 4), box_average(values, 16))
        self.assertEqual(box_average(values, 4).dtype, np.float64)

    def coarsened(self, fine, ratio):
        return box_average(fine, ratio).astype(np.float32)

    def test_level_consistency_factor_4(self):
        fine = fine_elevation(2, 1, 16)
        run = self.bundle(regional_fields(fine), hierarchy_blocks(2, level_entries(2)))
        coarse = self.coarsened(fine, 4)
        passing = self.check("level-consistency", run, self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1))))
        self.assertEqual(passing["status"], "pass", passing)
        self.assertIn("allowed 0.00195312 m (levels 2 to 2)", passing["summary"])
        coarse[2, 5] += 0.01
        failing = self.check("level-consistency", run, self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1))))
        self.assertEqual(failing["status"], "fail")
        self.assertIn("row 2, column 5 (global cell 3, 1)", failing["summary"])

    def test_level_consistency_factor_16(self):
        """Level 3 against level 1: 16 x 16 fine cells per coarse cell, allowed the bounds of levels 2 and 3.
        Heights of 1 to 2 km (float step 1.2e-4 m) so the offsets below are stored nearly exactly."""
        fine = fine_elevation(2, 1, 64, 1000.0, 2000.0)
        run = self.bundle(regional_fields(fine), hierarchy_blocks(3, level_entries(3)))
        coarse = self.coarsened(fine, 16) + np.float32(1.5 * BOUND)  # inside 2 x BOUND, outside one bound
        outcome = self.check("level-consistency", run, self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1))))
        self.assertEqual(outcome["status"], "pass", outcome)
        self.assertIn("factor 16", outcome["summary"])
        coarse = self.coarsened(fine, 16) + np.float32(2.5 * BOUND)
        outcome = self.check("level-consistency", run, self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1))))
        self.assertEqual(outcome["status"], "fail")

    def test_level_consistency_shared_window_across_the_date_line(self):
        """The run keeps global columns 7 and 0 (across the date line), the coarser bundle columns 0 and 1:
        compared on column 0 only; a difference on column 1 is outside the shared window."""
        fine = fine_elevation(2, 1, 16)
        run = self.bundle(regional_fields(fine), hierarchy_blocks(2, level_entries(2), column0=7))
        shared = self.coarsened(fine, 4)[:, 4:]
        coarse = np.concatenate([shared, np.full_like(shared, 50000.0)], axis=1)
        outcome = self.check("level-consistency", run,
                             self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1), column0=0)))
        self.assertEqual(outcome["status"], "pass", outcome)
        self.assertIn("1 x 1 global cells shared", outcome["summary"])

    def test_level_consistency_not_applicable_or_mismatched(self):
        fine = fine_elevation(2, 1, 16)
        run = self.bundle(regional_fields(fine), hierarchy_blocks(2, level_entries(2)))
        coarse = regional_fields(self.coarsened(fine, 4))
        absent = self.check("level-consistency", run)
        self.assertEqual(absent["status"], "skip")
        self.assertTrue(absent["summary"].startswith("not applicable"))
        single = self.check("level-consistency", run, self.bundle(coarse, hierarchy_blocks(1)))
        self.assertEqual(single["status"], "skip")
        other = self.check("level-consistency", run, self.bundle(coarse, hierarchy_blocks(1, level_entries(1), digest="b" * 64)))
        self.assertEqual(other["status"], "fail")
        apart = self.check("level-consistency", run, self.bundle(coarse, hierarchy_blocks(1, level_entries(1), column0=5)))
        self.assertEqual(apart["status"], "fail")

    def test_parent_mean_exact_null_error_fails(self):
        fields = regional_fields(fine_elevation(2, 1, 16))
        levels = level_entries(2)
        levels[1]["globalmeanerror_m"] = None
        outcome = self.check("parent-mean-exact", self.bundle(fields, hierarchy_blocks(2, levels)))
        self.assertEqual(outcome["status"], "fail")
        self.assertIn("level 2: globalmeanerror_m missing or null", outcome["details"][0])

    def test_level_consistency_mismatches_fail(self):
        """Another factor, another world grid, a scale that is not the levels' and a grid that is not the
        window's fail; a run without levels is not applicable."""
        fine = fine_elevation(2, 1, 16)
        coarse = self.bundle(regional_fields(self.coarsened(fine, 4)), hierarchy_blocks(1, level_entries(1)))
        factor2 = hierarchy_blocks(2, level_entries(2))
        factor2["window"]["factor"] = 2
        world = hierarchy_blocks(2, level_entries(2))
        world["window"]["world_columns"], world["window"]["world_rows"] = 16, 8
        scale = hierarchy_blocks(2, level_entries(2))
        scale["window"]["scale"] = 8
        for blocks in (factor2, world, scale):
            outcome = self.check("level-consistency", self.bundle(regional_fields(fine), blocks), coarse)
            self.assertEqual(outcome["status"], "fail", outcome)
        cut = self.check("level-consistency", self.bundle(regional_fields(fine[:, :16]), hierarchy_blocks(2, level_entries(2))), coarse)
        self.assertEqual(cut["status"], "fail")
        self.assertEqual(cut["summary"], "grid does not match the window")
        single = self.check("level-consistency", self.bundle(regional_fields(fine), hierarchy_blocks(2)), coarse)
        self.assertEqual(single["status"], "skip")
        self.assertIn("not applicable: the run is a single pass", single["summary"])

    def test_level_consistency_coarser_window_inside(self):
        """The run keeps global columns 1 to 4, the coarser bundle columns 2 and 3: compared on those two."""
        fine = fine_elevation(4, 1, 16)
        run = self.bundle(regional_fields(fine), hierarchy_blocks(2, level_entries(2), column0=1, columns=4))
        coarse = self.coarsened(fine, 4)[:, 4:12]
        outcome = self.check("level-consistency", run,
                             self.bundle(regional_fields(coarse), hierarchy_blocks(1, level_entries(1), column0=2, columns=2)))
        self.assertEqual(outcome["status"], "pass", outcome)
        self.assertIn("1 x 2 global cells shared", outcome["summary"])

    def test_level_consistency_pixel_grid_not_applicable(self):
        fine = fine_elevation(2, 1, 16)
        pixels = hierarchy_blocks(2, level_entries(2))
        pixels["window"].update({"pixel_columns": 20, "pixel_rows": 10, "pixel_cells": 1.6, "resampled": True})
        coarse = self.bundle(regional_fields(self.coarsened(fine, 4)), hierarchy_blocks(1, level_entries(1)))
        outcome = self.check("level-consistency", self.bundle(regional_fields(np.zeros((10, 20))), pixels), coarse)
        self.assertEqual(outcome["status"], "skip")
        self.assertIn("not applicable: the run is on a pixel grid", outcome["summary"])

    def test_regional_window_pixel_grid(self):
        """A window of 2 x 1 global cells at scale 16 (32 x 16 cells) in pixels of 1.5 cells: 21 x 10 pixels
        (partial edge pixels dropped). Another grid or other pixel counts fail; a cell-grid bundle still passes."""
        cells = self.check("regional-window", self.bundle(regional_fields(fine_elevation(2, 1, 16)), hierarchy_blocks(2, level_entries(2))))
        self.assertEqual(cells["status"], "pass", cells)

        def pixels(columns, rows, grid):
            blocks = hierarchy_blocks(2, level_entries(2))
            blocks["window"].update({"pixel_columns": columns, "pixel_rows": rows, "pixel_cells": 1.5, "resampled": True,
                                     "metres_per_pixel": 15000.0})
            return self.check("regional-window", self.bundle(regional_fields(np.zeros(grid)), blocks))

        self.assertEqual(pixels(21, 10, (10, 21))["status"], "pass")
        self.assertEqual(pixels(21, 10, (10, 22))["status"], "fail")
        self.assertEqual(pixels(22, 10, (10, 22))["status"], "fail")

    def test_resolution_invariance_reports(self):
        fine = fine_elevation(2, 1, 4)
        sea = np.zeros(fine.shape, dtype=np.uint8)
        sea[:, :2] = 1
        channel = np.zeros(fine.shape, dtype=np.int16)
        channel[0, 2:5] = 3  # 3 channel cells of 24 land cells
        channel[1, 0] = 3  # a channel code on a sea cell is not counted
        run = self.bundle(regional_fields(fine, sea, channel), hierarchy_blocks(1, level_entries(1)))
        alone = self.check("resolution-invariance", run)
        self.assertEqual(alone["status"], "reported")
        self.assertEqual(list(alone["values"]), ["run"])
        self.assertAlmostEqual(alone["values"]["run"]["channel_share"], 3 / 24)
        expected = np.percentile(fine[:, 2:].astype(np.float64) - 1000.0, [5, 25, 50, 75, 95])
        np.testing.assert_allclose(list(alone["values"]["run"]["hypsometry_m"].values()), expected)
        coarse = self.bundle(regional_fields(self.coarsened(fine, 4)), hierarchy_blocks(0, level_entries(1)))
        beside = self.check("resolution-invariance", run, coarse)
        self.assertEqual(beside["status"], "reported")
        self.assertEqual(list(beside["values"]), ["run", "coarser"])
        self.assertEqual(beside["values"]["coarser"]["channel_share"], 0.0)
        fields = regional_fields(fine)
        del fields["hydrosphere.channel_direction"]
        self.assertEqual(self.check("resolution-invariance", self.bundle(fields, hierarchy_blocks(1)))["status"], "skip")


if __name__ == "__main__":
    unittest.main()
