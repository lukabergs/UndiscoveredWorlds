"""Deterministic invariants for the workbook's spatial aggregation inputs."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

path = Path(__file__).resolve().parents[1]/'scripts/summarize-climate-validation.py'
spec = importlib.util.spec_from_file_location('validation_summary',path)
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)

class SpatialAggregationTests(unittest.TestCase):
    def test_spherical_area_is_independent_of_resolution_and_layout(self):
        expected=4*np.pi*summary.RADIUS**2/1e6
        for ny,nx,pole in [(33,64,True),(65,128,True),(32,64,False)]:
            self.assertAlmostEqual(summary.grid(ny,nx,pole)[2].sum()/expected,1,places=12)

    def test_conservative_regrid_preserves_integral(self):
        source=np.arange(7*12,dtype=float).reshape(7,12)**2
        old_area=summary.grid(7,12,True)[2]
        for ny,nx,pole in [(4,4,True),(6,8,False),(9,18,True)]:
            result=summary.remap(source,ny,nx,pole)
            new_area=summary.grid(ny,nx,pole)[2]
            self.assertAlmostEqual(np.sum(result*new_area)/np.sum(source*old_area),1,places=12)

    def test_paired_mask_excludes_missing_reference_and_zero_weight(self):
        result=summary.moments(np.array([1.,3.,100.,50.]),np.array([2.,4.,np.nan,0.]),np.array([1.,3.,9.,0.]))
        self.assertEqual(result['cells'],2)
        self.assertEqual(result['sim_mean'],2.5)
        self.assertEqual(result['ref_mean'],3.5)
        self.assertEqual(result['mse'],1)
        self.assertAlmostEqual(result['covariance']/np.sqrt(result['sim_variance']*result['ref_variance']),1)

    def test_constant_field_has_no_invented_variability(self):
        result=summary.moments(np.ones(4),np.zeros(4),np.ones(4))
        self.assertEqual(result['ref_variance'],0)
        self.assertEqual(result['sim_variance'],0)
        self.assertEqual(result['mse'],1)

    def test_ocean_mask_does_not_infer_water_from_temperature(self):
        fields={'east_current_mps':np.array([[0.,.1,0.,0.],[0.,0.,0.,0.]]),
                'south_current_mps':np.array([[0.,0.,0.,0.],[0.,0.,-.1,0.]]),
                'ice_thickness_m':np.array([[0.,0.,0.,1.],[0.,0.,0.,0.]]),
                'sst_c':np.array([[-20.,10.,15.,-1.8],[-10.,11.,16.,-1.8]])}
        np.testing.assert_array_equal(summary.identified_ocean_mask(fields),[False,True,True,True])

if __name__=='__main__':unittest.main()
