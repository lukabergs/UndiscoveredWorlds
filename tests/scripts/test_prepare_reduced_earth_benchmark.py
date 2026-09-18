"""Deterministic checks for arbitrary-resolution Earth benchmark maps."""

import importlib.util
from pathlib import Path
import csv
import hashlib
import json
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image


path = (
    Path(__file__).resolve().parents[2]
    / "scripts/refs/prepare-reduced-earth-benchmark.py"
)
sys.path.insert(0, str(path.parent))
spec = importlib.util.spec_from_file_location("earth_benchmark_maps", path)
maps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(maps)
import reference_map_rendering as rendering
import reference_climate_maps as climate_maps


class EarthBenchmarkMapTests(unittest.TestCase):
    def test_topography_has_exact_dimensions_and_complete_coverage(self):
        elevation = np.array(
            [
                [100, 200, -300, -400, 100, 200, -300, -400],
                [100, 200, -300, -400, 100, 200, -300, -400],
                [100, 200, -300, -400, 100, 200, -300, -400],
                [100, 200, -300, -400, 100, 200, -300, -400],
            ],
            dtype=np.float32,
        )
        land_mask = elevation > 0
        land, sea = maps.prepare_topography(elevation, land_mask, 4)
        self.assertEqual(land.shape, (2, 4))
        self.assertEqual(sea.shape, (2, 4))
        self.assertFalse(np.any((land > 0) & (sea > 0)))
        self.assertTrue(np.all((land > 0) | (sea > 0)))

    def test_writer_emits_importer_compatible_uint16_tiff(self):
        expected = np.arange(8, dtype=np.uint16).reshape(2, 4)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "map.tif"
            maps.write_tiff16(output, expected)
            encoded = output.read_bytes()
        self.assertEqual(encoded[:4], b"II*\x00")
        self.assertEqual(encoded[146:], expected.astype("<u2").tobytes())

    def test_csv_and_visual_map_match_the_grid(self):
        expected = np.arange(8, dtype=np.float32).reshape(2, 4)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            csv_path = directory / "map.csv"
            image_path = directory / "map.png"
            maps.write_grid_csv(csv_path, expected)
            maps.write_visual_map(image_path, expected, "precipitation")
            with csv_path.open(newline="", encoding="utf-8") as source:
                rows = list(csv.reader(source))
            with Image.open(image_path) as image:
                self.assertEqual(image.size, (4, 2))
        self.assertEqual(len(rows), 3)
        self.assertEqual([float(rows[1][1]), float(rows[2][1])], [45.0, -45.0])
        self.assertEqual(len(rows[1]) - 2, 4)

    def test_float_tiff_keeps_unclamped_signed_values_and_georeferencing(self):
        expected = np.array([[-54.125, 9000.5, np.nan, 0], [1, -1, 1e6, 1e-5]], dtype=np.float32)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "field.tif"
            rendering.write_float32_geotiff(output, expected)
            with Image.open(output) as image:
                self.assertEqual(image.mode, "F")
                self.assertEqual(image.tag_v2[258], (32,))
                self.assertEqual(image.tag_v2[339], (3,))
                self.assertEqual(image.tag_v2[277], 1)
                self.assertEqual(image.tag_v2[33550], (90, 90, 0))
                self.assertEqual(image.tag_v2[33922], (0, 0, 0, -180, 90, 0))
                self.assertEqual(image.tag_v2[34735][-4:], (2048, 0, 1, 4326))
                self.assertAlmostEqual(float(image.tag_v2[42113]), -9999.9)
                stored = np.asarray(image)
                np.testing.assert_array_equal(stored[np.isfinite(expected)], expected[np.isfinite(expected)])
                self.assertEqual(stored[0, 2], np.float32(-9999.9))

    def test_divergence_is_periodic_conservative_and_has_northward_sign(self):
        from cell_grid import spherical_area_weights, latitude_centres
        east = np.full((8, 16), 12, dtype=np.float32)
        north = np.zeros_like(east)
        np.testing.assert_array_equal(rendering.divergence(east, north), 0)
        north[:] = np.cos(np.deg2rad(latitude_centres(8)))[:, None]
        result = rendering.divergence(east, north)
        self.assertTrue(np.all(result[:4] < 0))
        self.assertTrue(np.all(result[4:] > 0))
        rng = np.random.default_rng(17)
        result = rendering.divergence(rng.normal(size=(8, 16)), rng.normal(size=(8, 16)))
        self.assertAlmostEqual(float(np.sum(result * spherical_area_weights(16))), 0, places=12)

    def test_flow_is_deterministic_and_northward_points_up(self):
        east = np.full((4, 8), 3.0)
        north = np.full((4, 8), 4.0)
        dx, dy, valid = rendering.direction(east, north, np.array(2.), np.array(2.))
        self.assertTrue(valid)
        self.assertGreater(dx, 0)
        self.assertLess(dy, 0)
        for render in (rendering.lic_luminance, rendering.particle_intensity):
            result = render(east, north, 17)
            np.testing.assert_array_equal(result, render(east, north, 17))
            self.assertEqual(result.shape, east.shape)
            self.assertTrue(np.all((result >= 0) & (result <= 1)))
        east[0, 0] = np.nan
        self.assertTrue(np.isnan(rendering.lic_luminance(east, north, 17)[0, 0]))
        row = np.arange(8, dtype=float)[None, :]
        self.assertEqual(rendering.sample(row, np.array(7.5), np.array(0.)), 3.5)

    def test_partial_quarter_and_mismatched_vector_coverage_stay_missing(self):
        months = np.ones((12, 2, 4))
        months[1, 0, 0] = np.nan
        self.assertTrue(np.isnan(climate_maps.quarter_mean(months, 0)[0, 0]))
        east, north = np.ones((2, 4)), np.ones((2, 4))
        east[:, ::2] = np.nan
        north[:, 1::2] = np.nan
        east, north = climate_maps.paired_remap(east, north, 2)
        self.assertTrue(np.isnan(east).all() and np.isnan(north).all())

    def test_full_cli_exports_mask_and_available_reference_families(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources, output = root / "sources", root / "output"
            sources.mkdir()
            elevation = np.array([[100, 200, -300, -400]] * 2, dtype=np.float32)
            Image.fromarray(elevation).save(sources / "height.tif")
            Image.fromarray((elevation > 0).astype(np.uint8)*255).save(sources / "mask.png")

            def bundle(filename, variable, value, layers=12):
                values = np.full((layers, 2, 4), value, dtype="<f4")
                header = maps.UWCLIM_MAGIC + struct.pack("<IIII", 1, 4, 2, layers)
                (sources / filename).write_bytes(header + variable.encode().ljust(16, b"\0") + values.tobytes())

            bundle("rain.uwclim", "prec_annual", 1200, 1)
            bundle("temperature.uwclim", "tavg", 25)
            bundle("monthly.uwclim", "prec", 100)
            for variable, value in (("u10m", 3), ("v10m", 4), ("u500", 0), ("v500", 5),
                                    ("u850", 5), ("v850", 0), ("tcwv", 30), ("slp_anom", -10),
                                    ("w500_ascent", -150), ("pr", 100), ("viwve", 300), ("viwvn", 400)):
                bundle(f"era5_{variable}_monthly.uwclim", variable, value)
            argv = [str(path), "--width", "4", "6", "--elevation", str(sources / "height.tif"),
                    "--land-mask", str(sources / "mask.png"), "--precipitation", str(sources / "rain.uwclim"),
                    "--temperature", str(sources / "temperature.uwclim"), "--diagnostics", str(sources),
                    "--monthly-precipitation", str(sources / "monthly.uwclim"), "--output-root", str(output)]
            with patch.object(sys, "argv", argv):
                maps.main()
            receipt = json.loads((output / "metadata/reference-maps/manifest.json").read_text())
            for resolution in receipt["outputs"]:
                w, h = resolution["width"], resolution["height"]
                products = {p["name"]: p for p in resolution["products"]}
                self.assertEqual(products["jan_surface_wind_speed"]["maps"][0]["file"],
                                 f"climate/maps/climate/wind/speed/1/s/{w}.png")
                self.assertEqual(products["jul_surface_divergence"]["maps"][1]["file"],
                                 f"climate/fields/climate/wind/divergence/3/s/{w}.tif")
                self.assertTrue({"koppen", "land_ocean", "jan_surface_wind_speed", "oct_upper_wind_lic",
                                 "apr_surface_wind_particles", "jul_surface_divergence", "jan_ascent",
                                 "jan_column_water", "jan_column_moisture_flux", "jan_moisture_flux_convergence"} <= products.keys())
                for product in products.values():
                    for record in (*product["maps"], product["csv"]):
                        file = output / record["file"]
                        self.assertEqual(hashlib.sha256(file.read_bytes()).hexdigest(), record["sha256"])
                    with Image.open(output / product["maps"][1]["file"]) as image:
                        self.assertEqual(image.size, (w, h))
                        self.assertEqual(image.mode, "F")
                with Image.open(output / products["land_ocean"]["maps"][0]["file"]) as image:
                    self.assertEqual(image.mode, "L")
                    self.assertEqual(set(np.unique(image)), {0, 255})
                with Image.open(output / products["jan_surface_wind_speed"]["maps"][1]["file"]) as image:
                    np.testing.assert_array_equal(np.asarray(image), 5)
                with Image.open(output / products["jan_column_moisture_flux"]["maps"][1]["file"]) as image:
                    np.testing.assert_array_equal(np.asarray(image), 500)
                with Image.open(output / products["jan_ascent"]["maps"][1]["file"]) as image:
                    np.testing.assert_array_equal(np.asarray(image), -150)
                self.assertTrue(any("consistency" in p["product"] for p in resolution["unavailable"]))

            # Regenerating one width must retain the other width and its hashes.
            with patch.object(sys, "argv", argv[:3] + argv[4:]):
                maps.main()
            partial = json.loads((output / "metadata/reference-maps/manifest.json").read_text())
            self.assertEqual(partial["outputs"], receipt["outputs"])
            bundle("oisst_sst_monthly.uwclim", "sst", 20)
            (sources / "oisst-additional-preparation.json").write_text(json.dumps(dict(
                family="oisst", width=4, height=2, products=[dict(bundle="oisst_sst_monthly.uwclim", variable="sst",
                units="degrees C", ocean_only=True, limit=35, palette="divergence",
                sources=["fixture"], period="2001-2020 monthly climatology")])) )
            with patch.object(sys, "argv", argv + ["--refresh-additional", "--additional-family", "oisst"]):
                maps.main()
            updated = json.loads((output / "metadata/reference-maps/manifest.json").read_text())
            for before, after in zip(receipt["outputs"], updated["outputs"]):
                products = {p["name"]: p for p in after["products"]}
                self.assertEqual(len(products), len(before["products"])+4)
                self.assertEqual(next(p for p in before["products"] if p["name"] == "jan_surface_wind_lic"),
                                 products["jan_surface_wind_lic"])
                self.assertIn("jan_oisst_sst", products)
                self.assertFalse(any(p["product"] == "sea_surface_temperature" for p in after["unavailable"]))
                self.assertTrue(any(p["product"] == "ocean_currents" for p in after["unavailable"]))

    def test_missing_references_are_reported_without_fabricated_maps(self):
        emitted = []
        missing = climate_maps.export_diagnostics({}, {}, 4, np.zeros((2, 4), dtype=bool),
                                                   lambda *args, **kwargs: emitted.append(args))
        self.assertEqual(emitted, [])
        self.assertTrue(any("moisture_flux_convergence" in item["product"] for item in missing))


if __name__ == "__main__":
    unittest.main()
