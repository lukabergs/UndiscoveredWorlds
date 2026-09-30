import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import rasterio
from rasterio.transform import from_bounds

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("geology_benchmark", ROOT / "scripts/benchmarks/benchmark-geology.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class GeologyBenchmarkTests(unittest.TestCase):
    def test_area_weights_and_known_distribution_distance(self):
        coarse = benchmark.area_weights(4, 8)
        fine = benchmark.area_weights(8, 16)
        self.assertAlmostEqual(float(coarse[:2].sum() / coarse.sum()), 0.5)
        self.assertAlmostEqual(float(coarse[0].sum() / coarse.sum()),
                               float(fine[:2].sum() / fine.sum()))
        self.assertAlmostEqual(benchmark.age_distribution_distance(
            np.array([0., 10.]), np.array([1., 1.]),
            np.array([3., 13.]), np.array([1., 1.])), 3.)
        self.assertAlmostEqual(benchmark.age_distribution_distance(
            np.array([0., 10.]), np.array([3., 1.]),
            np.array([0., 10.]), np.array([1., 3.])), 5.)

    def fixture(self, root):
        refs = root / "refs"
        products = []
        for key, (units, _, _) in benchmark.REFERENCES.items():
            products.append({"key": key, "units": units, "dataset": "test"})
            path = refs / "fields" / key / "8.tif"
            path.parent.mkdir(parents=True, exist_ok=True)
            data = np.full((4, 8), 20, dtype="float32")
            data[:, 0] = -9999
            with rasterio.open(path, "w", driver="GTiff", width=8, height=4, count=1,
                               dtype="float32", crs="EPSG:4326", nodata=-9999,
                               transform=from_bounds(-180, -90, 180, 90, 8, 4)) as ds:
                ds.write(data, 1)
        (refs / "manifest.json").write_text(json.dumps({"products": products}))
        bundle = root / "bundle"
        bundle.mkdir()
        manifest = {"schema": "tectonic-snapshot-bundle/v1", "contract_schema_version": 6,
                    "width": 8, "height": 4, "seed": 1, "time_myr": 5,
                    "run_scenario": {}, "files": {}}
        for name, dtype, value in [("heightmap", "<f4", .5), ("crust_age_myr", "<f4", 10),
                                    ("crust_class", "u1", 1), ("geologic_regime", "u1", 0),
                                    ("convergence_score", "u1", 0), ("divergence_score", "u1", 0),
                                    ("shear_score", "u1", 0)]:
            manifest["files"][name] = name + ".bin"
            np.full(32, value, dtype=dtype).tofile(bundle / (name + ".bin"))
        (bundle / "manifest.json").write_text(json.dumps(manifest))
        return refs, bundle

    def test_nodata_is_excluded_and_unit_mismatches_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            refs, bundle = self.fixture(Path(temp))
            result = benchmark.benchmark([bundle], refs, 8)
            age = result["references"]["tectonics/seafloor_age"]
            self.assertAlmostEqual(age["statistics_over_valid_cells"]["mean"], 20)
            self.assertAlmostEqual(age["valid_area_fraction"], 7 / 8)
            self.assertAlmostEqual(result["runs"][0]["age_wasserstein_distance_ma"], 10)
            self.assertTrue(result["runs"][0]["invariants_pass"])
            catalog = json.loads((refs / "manifest.json").read_text())
            catalog["products"][0]["units"] = "years"
            (refs / "manifest.json").write_text(json.dumps(catalog))
            with self.assertRaisesRegex(ValueError, "Unexpected units"):
                benchmark.benchmark([bundle], refs, 8)

    def test_collision_invariant_and_invalid_binary_shape(self):
        with tempfile.TemporaryDirectory() as temp:
            refs, bundle = self.fixture(Path(temp))
            np.full(32, 2, dtype="u1").tofile(bundle / "geologic_regime.bin")
            run = benchmark.benchmark([bundle], refs, 8)["runs"][0]
            self.assertFalse(run["invariants_pass"])
            self.assertEqual(run["oceanic_cells_with_continental_collision"], 32)
            (bundle / "crust_class.bin").write_bytes(b"\x01")
            with self.assertRaisesRegex(ValueError, "Invalid shape"):
                benchmark.read_bundle(bundle)


if __name__ == "__main__":
    unittest.main()
