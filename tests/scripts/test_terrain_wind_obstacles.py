import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
from terrain_wind import TerrainWindConfig,coordinates
from terrain_wind_obstacles import ObstacleWind,ObstacleConfig


class ObstacleWindTests(unittest.TestCase):
    def test_identity_calm_and_repeatability(self):
        rng=np.random.default_rng(252)
        h=rng.uniform(0,3000,(32,64)).astype(np.float32)
        u,v=rng.normal(size=(2,32,64)).astype(np.float32)
        a,b=ObstacleWind(h,h).apply(u,v)
        np.testing.assert_allclose(a,u,atol=2e-6);np.testing.assert_allclose(b,v,atol=2e-6)
        z=np.zeros_like(h)
        a,b=ObstacleWind(h,z).apply(z,z)
        self.assertFalse(np.any(a) or np.any(b))
        model=ObstacleWind(h,z)
        first=model.apply(u,v);second=model.apply(u,v)
        for a,b in zip(first,second):np.testing.assert_array_equal(a,b)

    def test_diversion_away_from_upstream_barrier(self):
        lat,lon=coordinates((128,256));radius=1e6
        h=(2500*np.exp(-((lon*radius/35000)**2+(lat*radius/60000)**2))).astype(np.float32)
        u=np.full_like(h,8);v=np.zeros_like(h)
        # Isolate look-ahead steering from the old local slope correction.
        c=ObstacleConfig(terrain=TerrainWindConfig(shelter=0,exposure=0,turning=0,channel=0))
        a,b=ObstacleWind(h,np.zeros((16,32)),c,radius).apply(u,v)
        north=(np.broadcast_to(lat,h.shape)>.015)&(np.broadcast_to(lat,h.shape)<.06)&(np.broadcast_to(lon,h.shape)<-.04)&(np.broadcast_to(lon,h.shape)>-.1)
        south=north[::-1]
        self.assertGreater(float(b[north].mean()),.1)
        self.assertLess(float(b[south].mean()),-.1)
        np.testing.assert_allclose(np.hypot(a,b),8,atol=2e-6)

    def test_bounds_longitude_symmetry_and_poles(self):
        lat,lon=coordinates((64,128));h=(3000*np.exp(-((lon/.2)**2+(lat/.15)**2))).astype(np.float32)
        u=np.full_like(h,8);v=np.full_like(h,1)
        config=ObstacleConfig(diversion_degrees=80,lookahead_m=150000,flank_distance_m=75000)
        model=ObstacleWind(h,np.zeros((16,32)),config,radius=1e6)
        a,b=model.apply(u,v);speed=np.hypot(a,b);original=np.hypot(u,v)
        self.assertTrue(np.isfinite(a).all() and np.isfinite(b).all())
        self.assertTrue(np.all(speed>=original*.09999) and np.all(speed<=original*1.60001))
        turn=np.abs(np.rad2deg(np.arctan2(u*b-v*a,u*a+v*b)))
        self.assertLessEqual(float(turn.max()),75.001)
        c,d=ObstacleWind(h[:,::-1],np.zeros((16,32)),config,radius=1e6).apply(-u,v)
        np.testing.assert_allclose(a[:,::-1],-c,atol=5e-4);np.testing.assert_allclose(b[:,::-1],d,atol=5e-4)
        c,d=ObstacleWind(np.roll(h,64,axis=1),np.zeros((16,32)),config,radius=1e6).apply(u,v)
        np.testing.assert_allclose(np.roll(a,64,axis=1),c,atol=5e-4)
        np.testing.assert_allclose(np.roll(b,64,axis=1),d,atol=5e-4)


if __name__=='__main__':unittest.main()
