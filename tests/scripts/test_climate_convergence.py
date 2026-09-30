import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
from climate_convergence import surface_convergence,tropical_band,smooth_vectors

class ConvergenceTests(unittest.TestCase):
    def test_closed_budget_and_cyclic_longitude(self):
        rng=np.random.default_rng(20260913);u,v=rng.normal(size=(2,32,64));c=surface_convergence(u,v)
        edges=np.linspace(np.pi/2,-np.pi/2,33);area=(np.sin(edges[:-1])-np.sin(edges[1:]))[:,None]
        self.assertLess(abs(np.sum(c*area)),1e-18)
        np.testing.assert_allclose(surface_convergence(np.roll(u,7,1),np.roll(v,7,1)),np.roll(c,7,1),atol=1e-20)
    def test_uniform_zonal_and_no_flow_have_no_band(self):
        u=np.ones((32,64))*8;v=np.zeros_like(u)
        np.testing.assert_array_equal(surface_convergence(u,v),v)
        self.assertFalse(tropical_band(v)[0].any())
    def test_equatorial_convergence_and_multiple_branches(self):
        lat=np.deg2rad(90-(np.arange(64)+.5)*180/64);u=np.zeros((64,128));v=np.broadcast_to(-20*np.sin(lat[:,None]),u.shape)
        c=surface_convergence(u,v);self.assertTrue((c[30:34]>0).all());self.assertTrue((c[0]<0).all())
        a=np.zeros((64,128));a[27:29]=1/86400;a[35:37]=.8/86400
        mask,_=tropical_band(a);self.assertTrue(mask[27:29].all() and mask[35:37].all());self.assertFalse(mask[30:34].any())
    def test_smoothing_is_shift_equivalent_and_invalid_inputs_rejected(self):
        rng=np.random.default_rng(14);u,v=rng.normal(size=(2,16,32))
        for a,b in zip(smooth_vectors(np.roll(u,3,1),np.roll(v,3,1)),smooth_vectors(u,v)):
            np.testing.assert_array_equal(a,np.roll(b,3,1))
        with self.assertRaises(ValueError):surface_convergence(u,v,0)
        with self.assertRaises(ValueError):surface_convergence(u,np.ones((2,2)))
        with self.assertRaises(ValueError):tropical_band(u,minimum_per_day=0)

if __name__=='__main__':unittest.main()
