"""Physical input contracts used by the reference-substitution experiments."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('reference_experiment_prepare',
    ROOT / 'runs/reports/climate-reference-experiment-prepare.py')
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


class ReferenceExperimentInputs(unittest.TestCase):
    def test_masked_restriction_retains_valid_values_and_missing_cells(self):
        field = np.full((4, 8), p.MISSING, dtype=np.float32)
        field[:2, :2] = 7.0
        field[0, 0] = p.MISSING
        result = p.remap(field, 4)
        self.assertEqual(result[0, 0], 7.0)
        self.assertTrue(np.all(result[:, 1:] == p.MISSING))
        self.assertEqual(result[1, 0], p.MISSING)

    def test_spherical_area_integral_is_conserved(self):
        field = np.arange(32, dtype=float).reshape(4, 8)
        coarse = p.remap(field, 4)
        def mean(a):
            edges = np.sin(np.linspace(np.pi / 2, -np.pi / 2, a.shape[0] + 1))
            weights = np.broadcast_to((edges[:-1] - edges[1:])[:, None], a.shape)
            return np.average(a, weights=weights)
        self.assertAlmostEqual(mean(field), mean(coarse), places=5)

    def test_dewpoint_humidity_is_specific_not_relative(self):
        humidity = p.qsat(np.array([20.0]), np.array([1000.0]))[0]
        self.assertGreater(humidity, .014)
        self.assertLess(humidity, .015)
        self.assertLess(humidity, p.qsat(np.array([28.0]), np.array([1000.0]))[0])


if __name__ == '__main__':
    unittest.main()
