import sys
import unittest
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts/benchmarks'))
from climate_palettes import rgb,DEFINITION
from climate_report_metrics import evaluate,rank
from climate_report_data import FIELDS,weights,remap
from climate_report_render import resize_scalar,raster
import climate_report
import climate_run_record
import climate_flow_texture
import reference_map_rendering

class StandardReportTests(unittest.TestCase):
    @unittest.skipIf(climate_flow_texture._library is None,'Build climate_flow_texture for accelerator equivalence')
    def test_compiled_flow_textures_match_numpy_reference(self):
        yy,xx=np.indices((16,32));east=10+8*np.sin(xx/5);north=3*np.cos(yy/3)+np.sin(xx/8)
        east[5:8,6:9]=np.nan;north[5:8,6:9]=np.nan
        for name in ('lic_luminance','particle_intensity'):
            expected=getattr(reference_map_rendering,name)(east,north,43)
            actual=getattr(climate_flow_texture,name)(east,north,43)
            np.testing.assert_allclose(actual,expected,atol=2e-6,rtol=0,equal_nan=True)
    def test_visual_baseline_choice_survives_a_report_rebuild(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);reports=root/'runs/reports';reports.mkdir(parents=True)
            selection=reports/'climate-global-current-selection.json'
            selection.write_text(json.dumps({'selected_run':10}))
            archive=SimpleNamespace(path=reports/'source',ids=[10,11],entries={10:('10',{}),11:('11',{})})
            for run in archive.ids:
                path=archive.path/str(run);path.mkdir(parents=True)
                (path/'manifest.json').write_text(json.dumps({'source_zip_sha256':f'source-{run}','executable_sha256':f'exe-{run}'}))
            ranking={'recommended_baseline':11,'best_in_batch':11,'baseline':10,'version':'fixture'}
            out=reports/'batch';out.mkdir()
            with patch.object(climate_report,'ROOT',root):
                climate_report.apply_selection(out,archive,ranking)
                self.assertEqual(json.loads(selection.read_text())['selected_run'],11)
                climate_report.apply_selection(out,archive,ranking,10,'User prefers rainfall pattern')
                climate_report.apply_selection(out,archive,ranking)
            self.assertEqual(json.loads(selection.read_text())['selected_run'],10)
            self.assertEqual(json.loads(selection.read_text())['selection_mode'],'visual_override')
            self.assertEqual(json.loads((out/'selection-before.json').read_text()),{'selected_run':10})
    def test_default_run_record_uses_numeric_id_and_raw_log_link(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);registry=root/'runs/registry';registry.mkdir(parents=True)
            (registry/'climate.json').write_text(json.dumps({'runs':[{'id':42,'horizontal_resolution':1024}]}))
            with patch.object(climate_run_record,'ROOT',root):climate_run_record.write_run_record(42,300)
            page=(root/'runs/reports/42.html').read_text()
            self.assertIn('<title>42</title>',page)
            self.assertIn('../logs/climate/42.log',page)
            self.assertIn(climate_run_record.MARKER,page)
    def test_shared_palette_endpoints_zero_missing_and_saturation(self):
        a=rgb(np.array([[-100.,-30.,0.,30.,100.,np.nan]]),-30,30,True)[0]
        np.testing.assert_array_equal(a[0],a[1]);np.testing.assert_array_equal(a[3],a[4])
        np.testing.assert_array_equal(a[2],[0,220,230]);np.testing.assert_array_equal(a[-1],[0,0,0])
        self.assertTrue(np.any(rgb(np.array([[0.]]))[0,0]!=0))
    def test_correct_histogram_does_not_hide_displaced_rainfall(self):
        ref=np.tile([0.,1.,4.,12.,0.,1.,4.,12.],(4,1));mask=np.ones_like(ref,bool)
        perfect=evaluate(ref,ref,mask,'rain');shifted=evaluate(np.roll(ref,2,axis=1),ref,mask,'rain')
        self.assertEqual(perfect['loss'],0)
        np.testing.assert_allclose(shifted['model_hist'],shifted['reference_hist'])
        self.assertGreater(shifted['loss'],.1)
    def test_reversed_and_calm_vectors_are_not_good_directions(self):
        ref=np.zeros((2,4,8));ref[0]=5;mask=np.ones((4,8),bool)
        flipped=evaluate(-ref,ref,mask,'surface_wind');calm=evaluate(np.zeros_like(ref),ref,mask,'surface_wind')
        rows=lambda d:{k:v for k,v,_ in d['rows']}
        self.assertAlmostEqual(rows(flipped)['Direction error (°)'],180)
        self.assertAlmostEqual(rows(calm)['Direction error (°)'],90)
        self.assertGreater(flipped['loss'],calm['loss'])
    def test_area_weights_and_conservative_registration(self):
        w=weights((4,8));self.assertGreater(w[1,0],w[0,0])
        np.testing.assert_allclose(remap(np.full((4,8),3.),4),3.)
    def test_scalar_export_interpolates_across_longitude_seam(self):
        a=np.zeros((2,4));a[:,-1]=8
        enlarged=resize_scalar(a,8)
        self.assertAlmostEqual(enlarged[1,0],2.)
        np.testing.assert_allclose(resize_scalar(np.full((2,4),7.),8),7.)
    def test_unsmoothed_export_preserves_native_cells_and_mask(self):
        values=np.array([[1,2,3,4],[5,6,7,8]],dtype=float);mask=np.ones_like(values,bool);mask[0,0]=False
        out,valid=raster('rain',values,mask,8,False)
        np.testing.assert_allclose(out[2:,2:4],6)
        self.assertFalse(valid[:2,:2].any())
        self.assertEqual(set(out[valid]),{2,3,4,5,6,7,8})
    def test_rank_has_fixed_weights_and_retains_baseline_on_tie(self):
        def values(loss):return {f:{str(s):{'loss':loss} for s in ([0] if f=='annual_rain' else range(4))} for f,v in FIELDS.items() if v[3]}
        data={'10':values(1),'11':values(1),'12':values(2)}
        result=rank(data,10,[10,11,12]);self.assertEqual(result['recommended_baseline'],10)
        data['11']=values(.9);result=rank(data,10,[10,11,12]);self.assertEqual(result['recommended_baseline'],11)
        self.assertAlmostEqual(result['improvement_percent'],10)
        del data['11']['pressure'];self.assertRaises(KeyError,rank,data,10,[10,11,12])
    def test_percentiles_include_out_of_palette_values(self):
        ref=np.ones((4,8));model=ref.copy();model[1:3,:]=100
        result=evaluate(model,ref,np.ones_like(ref,bool),'rain')
        self.assertGreater(result['model_percentiles']['0.99'],20)
        self.assertGreater(result['model_hist'][-1],0)
        self.assertAlmostEqual(sum(result['model_hist']),100)

if __name__=='__main__':unittest.main()
