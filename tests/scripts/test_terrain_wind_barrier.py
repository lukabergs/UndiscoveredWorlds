"""Opt-in global barrier checks: uv run --with numpy --with scipy --with pyamg ..."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
try:
    from terrain_wind_barrier import GlobalBarrier
except ModuleNotFoundError:
    GlobalBarrier = None


@unittest.skipIf(GlobalBarrier is None, 'Optional scipy/pyamg experiment dependencies')
class GlobalBarrierTests(unittest.TestCase):
    def test_flat_zonal_and_polar_faces(self):
        z = np.zeros((32, 64))
        u = np.broadcast_to((5+np.sin(np.arange(32)))[:, None], z.shape)
        r = GlobalBarrier(z).solve(u, z)
        np.testing.assert_allclose(r['east'], u, atol=1e-12)
        np.testing.assert_array_equal(r['north'], z)
        self.assertTrue(np.all(r['qy'][[0, -1]] == 0))

    def test_gap_fluxes_and_longitude_registration(self):
        y, x = np.indices((48, 96))
        z = np.zeros_like(x, dtype=float)
        z[((x-48)**2+(y-20)**2 < 9) | ((x-48)**2+(y-28)**2 < 9)] = 1000
        u = z*0+8
        v = z*0
        r = GlobalBarrier(z).solve(u, v)
        self.assertGreater(r['east'][24, 48], 8)
        self.assertLess(r['stats']['relative_flux_residual'], 1.2e-8)
        self.assertLess(abs(r['stats']['budget_error_m3s'])/r['stats']['absolute_flux_scale_m3s'], 1e-9)
        self.assertTrue(np.isfinite(r['east']).all())
        self.assertTrue(np.all(r['east'][z >= 500] == 0))
        shifted = GlobalBarrier(np.roll(z, 17, axis=1)).solve(u, v)
        np.testing.assert_allclose(np.roll(r['east'], 17, axis=1), shifted['east'], atol=2e-6)
        np.testing.assert_allclose(np.roll(r['north'], 17, axis=1), shifted['north'], atol=2e-6)
        mirrored = GlobalBarrier(z[:, ::-1]).solve(-u, v)
        np.testing.assert_allclose(r['east'], -mirrored['east'][:, ::-1], atol=2e-6)
        np.testing.assert_allclose(r['north'], mirrored['north'][:, ::-1], atol=2e-6)

    def test_invalid_input(self):
        for z in [np.zeros((3, 6)), np.zeros((4, 7)), np.full((4, 8), np.nan)]:
            with self.assertRaises(ValueError): GlobalBarrier(z)
        with self.assertRaises(ValueError): GlobalBarrier(np.zeros((4, 8)), vertical_mobility=0)
        grid = GlobalBarrier(np.zeros((4, 8)))
        with self.assertRaises(ValueError): grid.solve(np.zeros((4, 8)), np.full((4, 8), np.inf))


if __name__ == '__main__': unittest.main()
