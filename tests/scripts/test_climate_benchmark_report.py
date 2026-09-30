"""Historical report comparisons preserve support and original map grids."""
import importlib.util
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("benchmark_report", Path(__file__).resolve().parents[2] /
    "scripts/benchmarks/write-benchmark-report.py")
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


class HistoricalReportTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        previous = REPORT.ROOT
        REPORT.ROOT = Path(self.directory.name)
        self.addCleanup(setattr, REPORT, "ROOT", previous)

    def test_unsupported_scores_do_not_become_leaders(self):
        for metrics, key in (({"weighted_relative_error": 0}, "weighted_relative_error"),
            ({"spatial_compared_cells": 0, "spatial_exact_accuracy": 1}, "spatial_exact_accuracy"),
            ({"correlation": math.nan}, "correlation"), ({"correlation": 2}, "correlation")):
            self.assertIsNone(REPORT.recorded_metric({"metrics": metrics}, key))
        records = [{"id": 1, "metrics": {"weighted_relative_error": 0}},
                   {"id": 2, "metrics": {"weighted_relative_error": 0.4}}]
        self.assertEqual(REPORT.record_leader(records, "weighted_relative_error", False),
                         "0.400 · runs 2 (output width unavailable)")

    def test_historical_set_and_previous_are_explicit_without_run_161(self):
        records = [{"id": i, "metrics": {"era5_column_water_correlation": 0.7}} for i in (142, 143, 144, 145, 178)]
        table, markdown = REPORT.baseline_comparison(records, 178, 145)
        for i in (142, 143, 144, 145, 178):
            self.assertIn(f"Run {i}", table)
        self.assertNotIn("Run 161", table)
        self.assertIn("1 mm/month mean (12 mm/year)", markdown)
        self.assertIn("No monotonic target", markdown)

    def test_temperature_requires_compared_cells_and_nonplaceholder_error(self):
        path = REPORT.ROOT / "runs/diagnostics/climate/145/monthly_climate_reference_comparison.csv"
        path.parent.mkdir(parents=True)
        for cells, rmse, expected in ((0, 2, None), (10, 0, None), (10, 2, 2)):
            path.write_text(f"variable,period,compared_cells,area_weighted_rmse\ntemperature_c,annual_mean,{cells},{rmse}\n")
            self.assertEqual(REPORT.annual_temperature_rmse(145), expected)

    def png_header(self, path, width, height):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\0" * 8 + struct.pack(">II", width, height))

    def test_original_historical_dimensions_and_missing_maps_are_visible(self):
        path = REPORT.ROOT / "145.png"
        self.png_header(path, 512, 257)
        panel = REPORT.map_panel(path, REPORT.ROOT, "Historical run 145", "rain")
        self.assertIn("512×257 original pixels", panel)
        self.assertIn('height="257"', panel)
        missing = REPORT.map_panel(REPORT.ROOT / "absent.png", REPORT.ROOT, "Historical run 145", "rain")
        self.assertIn("Historical run 145: map unavailable", missing)

    def test_gallery_has_reference_previous_and_historical_panels(self):
        registry = REPORT.ROOT / "runs/registry/climate.json"
        registry.parent.mkdir(parents=True)
        registry.write_text(json.dumps({"runs": [{"id": 178, "horizontal_resolution": 4}]}))
        maps = REPORT.ROOT / "runs/maps/climate/rain/annual/s"
        for i, width, height in ((178, 4, 2), (174, 4, 2), (145, 512, 257)):
            self.png_header(maps / f"{i}.png", width, height)
        self.png_header(REPORT.ROOT / "refs/processed/climate/maps/climate/rain/annual/s/4.png", 4, 2)
        REPORT.write_report(178, compare_run_id=174)
        document = (REPORT.ROOT / "runs/reports/178.html").read_text(encoding="utf-8")
        self.assertEqual(document.count("<figure>"), 4)
        self.assertIn("Previous comparison run 174", document)
        self.assertIn("Historical run 145 · 512×257 original pixels", document)
        REPORT.write_report(178, compare_run_id=145)
        document = (REPORT.ROOT / "runs/reports/178.html").read_text(encoding="utf-8")
        self.assertEqual(document.count("<figure>"), 3)

    def test_mixed_wind_resolution_uses_product_metadata_and_matching_reference(self):
        registry = REPORT.ROOT / "runs/registry/climate.json"
        registry.parent.mkdir(parents=True)
        registry.write_text(json.dumps({"runs": [{"id": 273, "horizontal_resolution": 128}]}))
        metadata = REPORT.ROOT / "runs/diagnostics/climate/273/maps/wind_render_metadata.json"
        metadata.parent.mkdir(parents=True)
        data = dict(version=1, columns=2048, rows=1024, source_columns=128, source_rows=64,
                    products=['lic', 'particles', 'speed'])
        metadata.write_text(json.dumps(data))
        wind = REPORT.ROOT / "runs/maps/climate/wind/speed/3/u/273.png"
        self.png_header(wind, 2048, 1024)
        self.png_header(REPORT.ROOT / "runs/maps/climate/rain/annual/s/273.png", 128, 64)
        self.png_header(REPORT.ROOT / "refs/processed/climate/maps/climate/wind/speed/3/u/2048.png", 2048, 1024)
        REPORT.write_report(273)
        document = (REPORT.ROOT / "runs/reports/273.html").read_text(encoding='utf-8')
        self.assertIn('Run 273 · 2048×1024', document)
        self.assertIn('Reference · 2048×1024', document)
        self.assertIn('Run 273 · 128×64', document)
        # New metadata must never silently excuse malformed or old-size wind output.
        self.png_header(wind, 128, 64)
        with self.assertRaisesRegex(ValueError, 'Unexpected map dimensions'):
            REPORT.write_report(273)
        data['rows'] = 1025
        metadata.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'Invalid wind render geometry'):
            REPORT.write_report(273)


if __name__ == "__main__":
    unittest.main()
