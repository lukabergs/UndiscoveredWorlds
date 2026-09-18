import sys
import unittest
from dataclasses import replace
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
from terrain_wind import coordinates
from terrain_wind_obstacles import ObstacleWind,ObstacleConfig
from terrain_wind_wake import SpreadTerrainWind,SpreadObstacleWind,WakeConfig,Stencil


class WakeTests(unittest.TestCase):
    def fixture(self):
        yy,xx=np.indices((128,256));h=1500*np.exp(-((xx-128)/1.2)**2-((yy-64)/1.4)**2)
        return h.astype(np.float32),np.zeros_like(h,dtype=np.float32)

    def test_flat_same_terrain_calm_and_invalid(self):
        h,z=self.fixture();u=z+8
        for fine,coarse in ((z,z),(h,h)):
            a,b=SpreadObstacleWind(fine,coarse).apply(u,z)
            np.testing.assert_allclose(a,u,atol=1e-6);np.testing.assert_allclose(b,z,atol=1e-6)
        a,b=SpreadObstacleWind(h,z).apply(z,z)
        np.testing.assert_array_equal(a,z);np.testing.assert_array_equal(b,z)
        with self.assertRaises(ValueError):SpreadObstacleWind(h,z,wake=WakeConfig(crosswind_mixing_m=-1))
        with self.assertRaises(ValueError):SpreadObstacleWind(h,z,wake=WakeConfig(step_m=np.nan))
        with self.assertRaises(ValueError):SpreadObstacleWind(h,z).apply(z+np.nan,z)

    def test_wake_spreads_recovers_and_is_bounded(self):
        h,z=self.fixture();config=ObstacleConfig().terrain
        model=SpreadTerrainWind(h,z,config,WakeConfig(crosswind_mixing_m=6000),radius=1e6)
        wake=model.shelter(z+1,z)[0]
        narrow=SpreadTerrainWind(h,z,config,WakeConfig(crosswind_mixing_m=0),radius=1e6).shelter(z+1,z)[0]
        self.assertGreaterEqual(wake.min(),0);self.assertLessEqual(wake.max(),h.max())
        def width(column):
            y=np.arange(128);mean=np.average(y,weights=column)
            return np.sqrt(np.average((y-mean)**2,weights=column))
        self.assertGreater(width(wake[:,134]),width(wake[:,130]))
        self.assertGreater(width(wake[:,134]),width(narrow[:,134]))
        self.assertLess(wake[:,134].max(),wake[:,130].max())
        self.assertLess(wake[:,:123].max(),1e-4)
        self.assertGreater(wake[67,134],narrow[67,134])

    def test_unchanged_nonshelter_algebra(self):
        h,z=self.fixture();config=replace(ObstacleConfig(),terrain=replace(ObstacleConfig().terrain,shelter=0))
        u=z+8;v=z+2
        old=ObstacleWind(h,z,config,radius=1e6).apply(u,v)
        new=SpreadObstacleWind(h,z,config,radius=1e6).apply(u,v)
        for a,b in zip(old,new):np.testing.assert_array_equal(a,b)

    def test_symmetry_seam_poles_repeatability_and_bounds(self):
        h,z=self.fixture();u=z+8;v=z+1
        model=SpreadObstacleWind(h,z,radius=1e6);a,b=model.apply(u,v)
        for x,y in zip((a,b),model.apply(u,v)):np.testing.assert_array_equal(x,y)
        speed=np.hypot(a,b);parent=np.hypot(u,v)
        self.assertTrue(np.isfinite(speed).all());self.assertTrue((speed>=.1*parent-1e-5).all())
        self.assertTrue((speed<=1.6*parent+1e-5).all())
        aa,bb=SpreadObstacleWind(np.roll(h,128,axis=1),z,radius=1e6).apply(u,v)
        np.testing.assert_allclose(np.roll(a,128,axis=1),aa,atol=2e-4)
        np.testing.assert_allclose(np.roll(b,128,axis=1),bb,atol=2e-4)
        aa,bb=SpreadObstacleWind(h[:,::-1],z,radius=1e6).apply(-u,v)
        np.testing.assert_allclose(a,-aa[:,::-1],atol=2e-4)
        np.testing.assert_allclose(b,bb[:,::-1],atol=2e-4)
        lat,lon=coordinates(h.shape)
        np.testing.assert_allclose(Stencil(h.shape,lat,lon).apply(z+1),1,atol=2e-7)


if __name__=='__main__':unittest.main()
