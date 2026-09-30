import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('ocean_flow',ROOT/'scripts/benchmarks/ocean_flow_rendering.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class OceanFlowRenderingTests(unittest.TestCase):
    def test_physical_velocity_and_periodic_seam(self):
        e=np.full((8,16),.5);n=np.zeros_like(e)
        x=np.array([15.999]);y=np.array([3.5]);remaining=np.array([21600.])
        nx,ny,nr,valid=m.advance(e,n,x,y,remaining)
        self.assertTrue(valid[0]);self.assertLess(nx[0],1.)
        expected=.5*21600/(2*np.pi*6371000/16)
        self.assertAlmostEqual((nx[0]-x[0])%16,expected,places=10)
        self.assertEqual(ny[0],y[0]);self.assertEqual(nr[0],0.)

    def test_northward_sign_and_coast_stop(self):
        e=np.zeros((8,16));n=np.ones_like(e)
        nx,ny,_,valid=m.advance(e,n,np.array([4.]),np.array([4.]),np.array([21600.]))
        self.assertTrue(valid[0]);self.assertLess(ny[0],4.)
        e[:]=.5;n[:]=0.;e[:,8]=np.nan;n[:,8]=np.nan
        _,_,_,valid=m.advance(e,n,np.array([7.999]),np.array([4.]),np.array([21600.]))
        self.assertFalse(valid[0])

    def test_masked_interpolation_does_not_cross_land(self):
        e=np.ones((8,16));wet=np.ones_like(e,dtype=bool);wet[:,8]=False
        east,north,valid=m.resample(e,e*0,wet,64)
        self.assertFalse(valid[:,32:36].any())
        self.assertTrue(np.isnan(east[:,32:36]).all())
        self.assertTrue(np.all(east[valid]==1))

    def test_determinism_calm_and_small_artifact(self):
        e=np.full((8,16),.5);n=np.zeros_like(e)
        a=m.particles(e,n,7,days=1,count=100)
        b=m.particles(e,n,7,days=1,count=100)
        np.testing.assert_array_equal(a[0],b[0]);self.assertEqual(a[1:],b[1:])
        calm=m.particles(n,n,7,days=1,count=100)
        self.assertEqual(calm[0].max(),0.)
        with tempfile.TemporaryDirectory() as d:
            result=m.render_pair(e,n,e,n,np.ones_like(e,dtype=bool),Path(d)/'flow',width=64,days=1)
            self.assertEqual(result['native_grid'],[16,8])
            self.assertEqual(result['panels'][0]['valid_pixels'],result['panels'][1]['valid_pixels'])
            self.assertEqual(len(result['output_sha256']),2)


if __name__=='__main__':unittest.main()
