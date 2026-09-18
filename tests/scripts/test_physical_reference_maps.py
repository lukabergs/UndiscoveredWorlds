"""Deterministic public data-product invariants for the physical reference atlas."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/refs'))
import numpy as np
import pandas as pd
from PIL import Image
import rasterio
from physical_map_core import Atlas, Layer, WORLD_COVER, aggregate_points, regular_grid, rgba, provider_cpt, band_rgb, atlas_lock
from physical_map_fields import magnetic_local_components
from physical_map_points import whole_rock_rows
from physical_map_rasters import sinusoidal_pixel_cells, accumulate_equal_area


class PointAggregationTests(unittest.TestCase):
    def test_whole_rock_properties_keep_semantic_column_order(self):
        data=pd.DataFrame({'CAO(WT%)':[10,20,30],'MGO(WT%)':[7,8,9],
                           'SIO2(WT%)':[50,60,70],'MATERIAL':['WR','MIN','WR'],
                           'LONGITUDE MAX':[12,12,14],'LONGITUDE MIN':[12,12,12],
                           'LATITUDE MAX':[40,40,40],'LATITUDE MIN':[40,40,40]})
        np.testing.assert_array_equal(whole_rock_rows(data),[[40,40,12,12,50,7,10]])

    def test_dateline_and_poles_conserve_records(self):
        data=aggregate_points([-180,180,0,0],[0,0,90,-90],128)
        self.assertEqual(data.sum(),4)
        self.assertEqual(data[32,0],2)
        self.assertEqual(data[0,64],1)
        self.assertEqual(data[-1,64],1)

    def test_invalid_locations_do_not_wrap_into_valid_data(self):
        data=aggregate_points([0,190,np.nan,0],[0,0,0,91],128)
        self.assertEqual(data.sum(),1)

    def test_mean_and_missing_are_distinct_from_zero(self):
        data=aggregate_points([0,0,90],[0,0,0],128,[2,4,0])
        self.assertEqual(data[32,64],3)
        self.assertEqual(data[32,96],0)
        self.assertTrue(np.isnan(data[0,0]))

    def test_counts_are_rebinned_at_each_resolution(self):
        lon=[.1,.2,1.5];lat=[0,0,0]
        for w in (128,256,512,1024,2048):
            self.assertEqual(aggregate_points(lon,lat,w).sum(),3)


class GridTests(unittest.TestCase):
    def test_sinusoidal_aggregation_preserves_valid_pixel_sum_and_weight(self):
        cells=sinusoidal_pixel_cells(18,9,(20,20),128)
        values=np.full((20,20),3.);values[0,0]=np.nan
        total=np.zeros((64,128));weights=np.zeros_like(total)
        accumulate_equal_area(total,weights,values,cells)
        self.assertEqual(weights.sum(),399)
        self.assertEqual(total.sum(),1197)
        np.testing.assert_allclose(total[weights>0]/weights[weights>0],3)
        self.assertTrue((cells//128>=32).all())
        self.assertTrue((cells%128>=64).all())
        self.assertTrue((sinusoidal_pixel_cells(0,0,(20,20),128)==-1).all())

    def test_partial_latitude_coverage_is_not_filled_to_poles(self):
        data=regular_grid(np.ones((4,8)),[45,15,-15,-45],np.arange(-157.5,180,45),128)
        self.assertTrue(np.isnan(data[0]).all())
        self.assertTrue(np.isnan(data[-1]).all())
        self.assertTrue(np.nanmax(abs(data-1))<1e-6)

    def test_source_orientation_and_periodic_seam(self):
        data=np.array([[10]*9,[20]*9,[30]*9],dtype=float)
        out=regular_grid(data,[-60,0,60],np.arange(-180,181,45),128)
        self.assertGreater(np.nanmean(out[:10]),np.nanmean(out[-10:]))
        np.testing.assert_allclose(out[:,0],out[:,-1])

    def test_conflicting_duplicate_seam_is_rejected(self):
        data=np.ones((3,9));data[1,-1]=2
        with self.assertRaisesRegex(ValueError,'Conflicting'):
            regular_grid(data,[-60,0,60],np.arange(-180,181,45),128)

    def test_classes_remain_integer_codes(self):
        data=np.tile([10,80,10,80,10,80,10,80],(4,1))
        out=regular_grid(data,[67.5,22.5,-22.5,-67.5],np.arange(-157.5,180,45),128,categorical=True)
        self.assertEqual(set(out[np.isfinite(out)]),{10,80})


class ColourAndExportTests(unittest.TestCase):
    def test_atlas_rejects_concurrent_writer_and_releases_lock(self):
        with tempfile.TemporaryDirectory() as root:
            with atlas_lock(root):
                with self.assertRaisesRegex(RuntimeError,'Another process'):
                    with atlas_lock(root):pass
            with atlas_lock(root):pass

    def test_dipole_north_down_signs_and_rotation_preserves_intensity(self):
        latitude=np.array([-30.,0.,30.])
        radial=(-60000*np.sin(np.deg2rad(latitude)))[:,None]
        theta=(-30000*np.cos(np.deg2rad(latitude)))[:,None]
        east=np.zeros_like(theta)
        _,n,e,d=magnetic_local_components(radial,theta,east,latitude,0)
        np.testing.assert_allclose(n,-theta)
        np.testing.assert_allclose(d,-radial)
        self.assertLess(d[0,0],0)
        self.assertGreater(d[2,0],0)
        _,n,e,d=magnetic_local_components(radial,theta,east,latitude,1/298.257223563)
        np.testing.assert_allclose(n*n+e*e+d*d,radial*radial+theta*theta,rtol=1e-12)

    def test_worldcover_official_colours_and_missing_alpha(self):
        layer=Layer('land','Land','class','test',(10,100),classes=WORLD_COVER)
        out=rgba(np.array([[10,80,np.nan]]),layer)
        self.assertEqual(out[0,0].tolist(),[0,100,0,255])
        self.assertEqual(out[0,1].tolist(),[0,100,200,255])
        self.assertEqual(out[0,2].tolist(),[0,0,0,0])

    def test_provider_cpt_thresholds_are_not_smoothed(self):
        colours=band_rgb(np.array([1,49,51,199]),provider_cpt('fullrategrid.cpt'))
        self.assertEqual(colours[0].tolist(),[255,0,255])
        self.assertEqual(colours[1].tolist(),[0,85,255])
        self.assertEqual(colours[2].tolist(),[0,216,255])
        self.assertEqual(colours[3].tolist(),[255,73,0])

    def test_export_keeps_unclamped_values_and_body_crs(self):
        with tempfile.TemporaryDirectory() as root:
            atlas=Atlas(root,Path(root)/'output',[128])
            values=np.full((64,128),500.,dtype='float32');values[0,0]=np.nan
            layer=Layer('moon/test','Moon test','m','test',(0,100),body='Moon')
            atlas.emit(layer,lambda w:values)
            path=atlas.output/'fields/moon/test/128.tif'
            with rasterio.open(path) as data:
                self.assertEqual(data.shape,(64,128))
                self.assertEqual(data.bounds.left,-180)
                self.assertEqual(data.bounds.top,90)
                self.assertEqual(data.tags()['units'],'m')
                self.assertEqual(data.read(1)[1,1],500)
                self.assertNotEqual(data.crs.to_epsg(),4326)
                self.assertTrue(data.read(1,masked=True).mask[0,0])
            image=np.asarray(Image.open(atlas.output/'maps/moon/test/128.png'))
            self.assertEqual(image[0,0,3],0)
            manifest=json.loads((atlas.output/'manifest.json').read_text())
            self.assertEqual(manifest['products'][0]['outputs'][0]['max'],500)


if __name__=='__main__':unittest.main()
