"""Path contracts and collision-safe migration, using isolated fixtures."""

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts/refs"
sys.path.insert(0, str(SCRIPTS))
from reference_map_layout import product_directory, product_path

spec = importlib.util.spec_from_file_location("organize", SCRIPTS / "organize-climate-maps.py")
organize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(organize)


class MapLayoutTests(unittest.TestCase):
    def test_july_is_three_and_width_is_the_filename(self):
        self.assertEqual(str(product_path("jul_surface_wind_lic", 2048)),
                         "climate/maps/climate/wind/lic/3/s/2048.png")
        self.assertEqual(str(product_path("oct_upper_divergence", 512, "fields")),
                         "climate/fields/climate/wind/divergence/4/u/512.tif")

    def test_providers_and_layer_components_do_not_collide(self):
        names = ["jul_oisst_sst", "jul_era5_sst", "jul_glorys_thetao",
                 "jul_column_moisture_flux", "jul_column_moisture_flux_east",
                 "jul_column_moisture_flux_north", "jul_moisture_flux_convergence", "jul_era5_vimd"]
        self.assertEqual(len({product_directory(n) for n in names}), len(names))
        with self.assertRaises(KeyError):
            product_directory("jul_unknown_quantity")

    def test_legacy_run_categories_and_months(self):
        for old, new in {
            "wind/lic/2/s": "wind/lic/3/s",
            "wind/cons/0_u": "wind/consistency/1/u",
            "wind/part/3/u": "wind/particles/4/u",
            "wind/v/0/s": "wind/north/1/s",
            "moisture/flux/0_b": "moisture/flux/1/b",
            "s_div/0": "wind/divergence/1/s",
            "air_temperature": "air_temp/annual/s",
        }.items():
            self.assertEqual(organize.legacy_run_directory(old.split("/")), new)

    def test_staging_preserves_files_when_destination_is_another_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            manifest = root / "refs/processed/metadata/reference-maps/manifest.json"
            manifest.parent.mkdir(parents=True)
            manifest.write_text('{"outputs": []}', encoding="utf-8")
            journal_path = root / "runs/manifests/climate-map-layout.json"
            journal_path.parent.mkdir(parents=True)
            moves = {"runs/maps/climate/wind/lic/0/s/7.png": "runs/maps/climate/wind/lic/1/s/7.png",
                     "runs/maps/climate/wind/lic/1/s/7.png": "runs/maps/climate/wind/lic/2/s/7.png"}
            for i, source in enumerate(moves):
                path = root / source
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(bytes([i]))
            journal = {"version": 1, "phase": "staging", "moves": moves}
            with patch.object(organize, "ROOT", root), patch.object(organize, "JOURNAL", journal_path):
                organize.apply(journal)
                for i, destination in enumerate(moves.values()):
                    self.assertEqual((root / destination).read_bytes(), bytes([i]))
                self.assertEqual(json.loads(journal_path.read_text())["phase"], "complete")
                with self.assertRaises(ValueError):
                    organize.checked(root / "../outside")


if __name__ == "__main__":
    unittest.main()
