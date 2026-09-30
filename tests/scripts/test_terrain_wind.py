import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts/benchmarks"))
from terrain_wind import TerrainWind,TerrainWindConfig,coordinates,destination,sample


class TerrainWindTests(unittest.TestCase):
    def test_flat_and_identical_terrain(self):
        rng=np.random.default_rng(250)
        u,v=rng.normal(size=(2,32,64)).astype(np.float32)
        flat=np.full((32,64),500,dtype=np.float32)
        a,b=TerrainWind(flat,np.full((8,16),500)).apply(u,v)
        np.testing.assert_array_equal(a,u)
        np.testing.assert_array_equal(b,v)
        h=rng.uniform(0,3000,(32,64)).astype(np.float32)
        a,b=TerrainWind(h,h).apply(u,v)
        np.testing.assert_allclose(a,u,atol=2e-5)
        np.testing.assert_allclose(b,v,atol=2e-5)

    def test_island_wake_calm_bound_and_wrap(self):
        rows,cols=64,128
        y,x=np.indices((rows,cols))
        # Small planet fixes physical spacing to ~10 km for this fixture.
        h=(2500*np.exp(-((x-64)**2+(y-32)**2)/5)).astype(np.float32)
        model=TerrainWind(h,np.zeros((16,32)),radius=200000.)
        u=np.full(h.shape,8,dtype=np.float32); v=np.zeros_like(u)
        a,b=model.apply(u,v)
        self.assertLess(a[32,70],u[32,70])
        self.assertGreater(a[32,64],u[32,64])
        self.assertTrue(np.isfinite(a).all() and np.isfinite(b).all())
        self.assertLessEqual(float(np.hypot(a-u,b-v).max()),4.0001)
        z1,z2=model.apply(v,v)
        self.assertFalse(np.any(z1) or np.any(z2))
        shifted=TerrainWind(np.roll(h,64,axis=1),np.zeros((16,32)),radius=200000.)
        c,d=shifted.apply(u,v)
        np.testing.assert_allclose(np.roll(a,64,axis=1),c,atol=2e-4)
        np.testing.assert_allclose(np.roll(b,64,axis=1),d,atol=2e-4)

    def test_mirror_and_polar_coordinates(self):
        y,x=np.indices((64,128))
        h=(2000*np.exp(-((x-64)**2+(y-32)**2)/7)).astype(np.float32)
        u=np.full(h.shape,8,dtype=np.float32); v=np.zeros_like(u)
        a,b=TerrainWind(h,np.zeros((16,32)),radius=200000.).apply(u,v)
        c,d=TerrainWind(h[:,::-1],np.zeros((16,32)),radius=200000.).apply(-u,v)
        np.testing.assert_allclose(a[:,::-1],-c,atol=2e-4)
        np.testing.assert_allclose(b[:,::-1],d,atol=2e-4)
        lat,lon=destination(np.deg2rad(89.9),np.deg2rad(30.),0,1,50000,6371000)
        self.assertTrue(np.isfinite(lat) and np.isfinite(lon))
        self.assertLess(np.rad2deg(lat),89.9)

    def test_gap_and_resolution(self):
        y,x=np.indices((64,128))
        h=(2500*np.exp(-(x-64)**2/5)*(1-np.exp(-(y-32)**2/2))).astype(np.float32)
        u=np.full(h.shape,8,dtype=np.float32); v=np.zeros_like(u)
        a,b=TerrainWind(h,np.zeros((16,32)),radius=200000.).apply(u,v)
        self.assertGreater(a[32,64],8)
        self.assertTrue(np.any(np.abs(b)>0.05))
        fields=[]
        for width in (128,256,512):
            lat,lon=coordinates((width//2,width))
            h=2000*np.exp(-((lon/.08)**2+(lat/.08)**2))
            u=np.full(h.shape,8,dtype=np.float32); v=np.zeros_like(u)
            a,b=TerrainWind(h,np.zeros((16,32)),radius=200000.).apply(u,v)
            fields.append(sample(a,np.array([0.,.03,-.03]),np.array([.15,.1,.1])))
        coarse_error=np.max(np.abs(fields[0]-fields[1]))
        fine_error=np.max(np.abs(fields[1]-fields[2]))
        self.assertLess(fine_error,coarse_error)
        self.assertLess(fine_error,.2)


if __name__=="__main__": unittest.main()
