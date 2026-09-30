"""Offline atlas index and independent output verification."""
import dataclasses
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
import rasterio

from physical_map_core import REPO, Layer, colour_table, legend, rgba, sha256, coverage_bounds

NONSPATIAL={
    'phreeqc-3.8.6':'Thermodynamic reaction databases; water-chemistry maps require spatial water compositions and a reaction/transport model.',
    'usgs-spectral-v7-ascii':'Laboratory/field spectra and material descriptions; no continuous geographic distribution to rasterise.',
    'usgs-deposit-models':'Descriptive and grade-tonnage models, not spatial observations.',
    'crism-mica-v1':'Mars mineral type spectra; downloaded archive is not a planetary mineral-distribution map.',
}


def layer_from_record(record):
    names={f.name for f in dataclasses.fields(Layer)}
    return Layer(**{k:v for k,v in record.items() if k in names})


def render_existing(root):
    root=Path(root);manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    for record in manifest['products']:
        layer=layer_from_record(record)
        record['palette_source']=layer.palette_source
        if layer.method.startswith('Record count'):
            layer.zero_colour='#111820';record['zero_colour']=layer.zero_colour
        legend(root/record['legend'],layer)
        for output in record['outputs']:
            with rasterio.open(root/output['tiff']) as src:
                values=src.read(1,masked=True).filled(np.nan)
            path=root/output['png']
            Image.fromarray(rgba(values,layer)).save(path)
            output['png_sha256']=sha256(path)
            output['coverage_pixel_bounds']=coverage_bounds(values)
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')


def verify(root,expected_widths):
    root=Path(root);manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    results=[];errors=[]
    keys=[p['key'] for p in manifest['products']]
    if len(set(keys))!=len(keys):errors.append('Duplicate map keys')
    for record in manifest['products']:
        layer=layer_from_record(record)
        if sorted(o['width'] for o in record['outputs'])!=sorted(expected_widths):
            errors.append(record['key']+': incomplete resolution set')
        if not (root/record['legend']).exists():errors.append(record['key']+': missing legend')
        for out in record['outputs']:
            try:
                width=out['width'];png=root/out['png'];tif=root/out['tiff']
                if sha256(png)!=out['png_sha256'] or sha256(tif)!=out['tiff_sha256']:
                    raise ValueError('Output checksum mismatch')
                with rasterio.open(tif) as data:
                    values=data.read(1,masked=True).filled(np.nan)
                    assert data.shape==(width//2,width),'Incorrect field shape'
                    assert np.allclose(tuple(data.bounds),[-180,-90,180,90]),'Incorrect geographic bounds'
                    assert data.tags()['units']==record['units'],'Units metadata mismatch'
                    assert data.tags()['AREA_OR_POINT']=='Area','Not a cell grid'
                    assert data.transform.e<0,'Not north-up'
                    assert data.dtypes==('float32',),'Not a float32 field'
                    if record['body']=='Earth':assert data.crs.to_epsg()==4326,'Earth CRS mismatch'
                    else:assert data.crs.to_epsg()!=4326,'Lunar data labelled as Earth'
                image=np.asarray(Image.open(png))
                assert image.shape==(width//2,width,4),'Incorrect image shape'
                np.testing.assert_array_equal(image,rgba(values,layer),err_msg='PNG does not correspond to numeric field / declared palette')
                count=int(np.isfinite(values).sum())
                assert count==out['valid_cells'],'Valid-cell count mismatch'
                assert coverage_bounds(values)==out['coverage_pixel_bounds'],'Coverage bounds mismatch'
                if count==0:raise ValueError('Map has no valid cells')
                if record['key'].startswith(('crust/','gravity/','moon/')) or (record['key'].startswith('magnetism/') and record['key']!='magnetism/declination'):
                    assert count==width*width//2,'Unexpected gaps in a global model field'
                if layer.classes:
                    assert set(np.unique(values[np.isfinite(values)]))<=set(map(int,layer.classes)),'Invented category codes'
                assert png.with_suffix('.pgw').exists() and png.with_suffix('.prj').exists(),'Missing PNG georeferencing'
                results.append(dict(key=record['key'],width=width,status='passed'))
            except Exception as exc:
                errors.append(f"{record['key']} / {out['width']}: {exc}")
    errors.extend(x['family']+': '+x.get('error','failed') for x in manifest['datasets'] if x['status']!='complete')
    report=dict(products=len(manifest['products']),outputs_checked=len(results),errors=errors,
                checks='Hashes; field/image agreement; dimensions; CRS/body; units; masks; valid category codes; complete requested widths. Not scientific accuracy or simulator validation.',results=results)
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='results'}),flush=True)
    return int(bool(errors))


def index(root):
    root=Path(root);manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    manifest['nonspatial_inputs']=NONSPATIAL
    manifest['model_dependent_outputs']=['Weathering/erosion rates','Soil development','Spatial water chemistry','Biome/habitat suitability','Undiscovered mineral potential']
    manifest['model_dependent_explanation']='These require additional process models and boundary conditions; they are not generated as observed reference maps.'
    manifest['palette_assets']={str(p.relative_to(REPO)).replace('\\','/'):sha256(p) for p in (REPO/'assets/palettes/physical').glob('*') if p.is_file()}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    backdrop=root/'earth-context.png'
    source=REPO/'refs/processed/land_ocean/maps/earth_land_ocean_2048x1024.png'
    if source.exists():
        mask=np.asarray(Image.open(source).convert('L').resize((1024,512),Image.Resampling.NEAREST))>0
        image=np.empty((512,1024,3),'uint8');image[:]=[9,18,27];image[mask]=[40,52,60]
        Image.fromarray(image).save(backdrop)
    data=json.dumps(manifest,ensure_ascii=False).replace('</',r'<\/')
    template=Path(__file__).with_name('physical_map_index.html').read_text(encoding='utf-8')
    (root/'index.html').write_text(template.replace('__ATLAS_DATA__',data),encoding='utf-8')
    # A bounded overview for visual inspection and sharing; numeric maps remain separate.
    desired=['geology/lithology','crust/solid_crust_thickness','tectonics/seafloor_age','mineralogy/hematite',
             'soils/phh2o','ocean/o/1000m','hydrology/rivers_discharge','cryosphere/glacier_presence',
             'biosphere/npp','gravity/earth_geoid','magnetism/intensity','moon/elevation']
    products={p['key']:p for p in manifest['products']}
    selected=[products[k] for k in desired if k in products]
    canvas=Image.new('RGB',(1536,300*((len(selected)+2)//3)),'#111820');draw=ImageDraw.Draw(canvas)
    from physical_map_core import font
    for i,product in enumerate(selected):
        x=(i%3)*512;y=(i//3)*300
        output=min(product['outputs'],key=lambda o:abs(o['width']-512))
        image=Image.open(root/output['png']).convert('RGBA').resize((496,248))
        if product['body']=='Earth' and backdrop.exists():
            base=Image.open(backdrop).convert('RGBA').resize((496,248))
        else:base=Image.new('RGBA',(496,248),'#09121b')
        base.alpha_composite(image);canvas.paste(base.convert('RGB'),(x+8,y+8))
        draw.text((x+10,y+263),product['title'][:58],font=font(13),fill='white')
        draw.text((x+10,y+282),product['units'],font=font(11),fill='#a8bbc9')
    canvas.save(root/'overview.png')
