import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
from terrain_wind_mass import TerrainMassGrid, MassConfig


class TerrainMassTests(unittest.TestCase):
    def grid(self, h):
        return TerrainMassGrid(h, 46-(np.arange(h.shape[0])+.5)*.17578125, .17578125)

    def test_flat_identity_and_calm(self):
        h = np.zeros((21,31)); g = self.grid(h)
        for speed in (0.,8.):
            u = h+speed; r = g.solve(u,h,u,h)
            np.testing.assert_array_equal(r['east'],u)
            np.testing.assert_array_equal(r['north'],h)
            np.testing.assert_array_equal(r['vertical_layers'],0)
            self.assertEqual(r['stats']['iterations'],0)

    def test_ridge_vertical_flow_and_local_global_closure(self):
        yy,xx = np.indices((27,41)); h = 2000*np.exp(-((xx-20)/3)**2-((yy-13)/9)**2)
        g = self.grid(h); u = h*0+8; v = h*0
        before = g.initialize(u,v,u,v); r = g.solve(u,v,u,v)
        self.assertLess(r['stats']['relative_flux_residual'],1.1e-8)
        self.assertLess(r['stats']['maximum_divergence_per_second'],1e-10)
        self.assertGreater(np.abs(r['vertical_layers']).max(),.01)
        np.testing.assert_array_equal(r['qz'][0],0)
        np.testing.assert_array_equal(r['qx'][:,:,[0,-1]],before[0][:,:,[0,-1]])
        np.testing.assert_array_equal(r['qy'][:,[0,-1]],before[1][:,[0,-1]])
        self.assertLess(abs(r['stats']['budget_error_m3s']),1e-7*np.sum(np.abs(before[0])))
        for key in ('east','north','qx','qy','qz'):
            np.testing.assert_array_equal(r[key],g.solve(u,v,u,v)[key])
        mirror = self.grid(h[:,::-1]).solve(-u[:,::-1],v[:,::-1],-u[:,::-1],v[:,::-1])
        np.testing.assert_allclose(r['east'],-mirror['east'][:,::-1],atol=1e-8)
        np.testing.assert_allclose(r['north'],mirror['north'][:,::-1],atol=1e-8)

    def test_divergent_input_top_exchange_and_spd(self):
        yy,xx = np.indices((23,35)); h = xx*0.; g = self.grid(h)
        u = (xx-17)/3; v = -(yy-11)/3
        r = g.solve(u,v,u,v)
        self.assertLess(r['stats']['relative_flux_residual'],1.1e-8)
        self.assertLess(r['stats']['top_outflow_m3s'],0)
        self.assertGreater(r['stats']['lateral_outflow_m3s'],0)
        rng = np.random.default_rng(1729); a,b = rng.normal(size=(2,*g.volume.shape))
        self.assertAlmostEqual(float(np.sum(a*g.operator(b))/np.sum(b*g.operator(a))),1,places=12)
        self.assertGreater(np.sum(a*g.operator(a)),0)

    def test_invalid_input_and_failed_solve(self):
        h = np.zeros((7,9)); g = self.grid(h)
        with self.assertRaises(ValueError):g.solve(h+np.nan,h,h,h)
        with self.assertRaises(ValueError):TerrainMassGrid(h,np.arange(7),1)
        with self.assertRaises(ValueError):TerrainMassGrid(h,90-np.arange(7),1)
        h[3,4] = 2000
        g = TerrainMassGrid(h,45-np.arange(7)*.2,.2,MassConfig(maximum_iterations=1))
        with self.assertRaises(RuntimeError):g.solve(h*0+8,h*0,h*0+8,h*0)


if __name__ == '__main__':unittest.main()
