"""Fast semantic tests for the climate reference regression gate."""
import copy
import runpy
import unittest
from pathlib import Path

compare = runpy.run_path(str(Path(__file__).resolve().parents[2] /
    "scripts/benchmarks/check-climate-baseline.py"))["compare_reference_skill"]


class ReferenceGateTests(unittest.TestCase):
    def setUp(self):
        metrics = dict.fromkeys(("area_weighted_land_precipitation_correlation",
            "era5_column_water_correlation", "era5_surface_eastward_wind_correlation",
            "era5_surface_northward_wind_correlation", "era5_global_pressure_correlation"), 0.5)
        self.baseline = dict(run_id=1, historical_metrics={"1": metrics},
            tropical_structure={s: dict(vector_rmse_mps=6.0) for s in ("1", "3")},
            ocean_reference_scores={"mean_sst_c_quarter_1": dict(rmse=5.0)})
        self.current = copy.deepcopy(self.baseline)
        self.current["run_id"] = 2
        self.current["historical_metrics"] = {"2": copy.deepcopy(metrics)}

    def test_equal_or_improved_targets_pass(self):
        self.assertEqual(compare(self.current, self.baseline)["status"], "no_detected_regression")
        self.current["tropical_structure"]["1"]["vector_rmse_mps"] = 3.0
        self.assertEqual(compare(self.current, self.baseline)["status"], "no_detected_regression")

    def test_better_wind_cannot_hide_worse_rain(self):
        self.current["tropical_structure"]["1"]["vector_rmse_mps"] = 3.0
        self.current["historical_metrics"]["2"]["area_weighted_land_precipitation_correlation"] = 0.05
        result = compare(self.current, self.baseline)
        self.assertEqual(result["status"], "regressed")
        self.assertTrue(result["metrics"]["area_weighted_land_precipitation_correlation"]["regressed"])

    def test_missing_metric_is_not_zero_error(self):
        self.current["ocean_reference_scores"]["mean_sst_c_quarter_1"]["rmse"] = None
        with self.assertRaisesRegex(ValueError, "Missing reference metric"):
            compare(self.current, self.baseline)

    def test_historical_scores_do_not_change_previous_run_gate(self):
        self.current["historical_metrics"]["145"] = {"area_weighted_land_precipitation_correlation": 1.0}
        self.current["historical_comparison"] = {"role": "Informational only"}
        self.assertEqual(compare(self.current, self.baseline)["status"], "no_detected_regression")


if __name__ == "__main__":
    unittest.main()
