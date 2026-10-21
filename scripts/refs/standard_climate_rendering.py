"""Standard climate reference images: field table, reference readers and the shared map, LIC and particle renderers.

Used by render_standard_climate_maps.py and render_standard_flow_references.py (refs/processed/climate/maps/standard-v2).
Moved unchanged from the retired climate gallery (2026-10-05); only its native-cell export was dropped.
No measurement is extracted from colours.
"""
import html
import json
import os
from pathlib import Path
import sys
import numpy as np
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'scripts/refs'))
sys.path.insert(0,str(ROOT/'scripts/benchmarks'))
from uwtools.cell_grid import conservative_remap
from uwtools import references
from uwtools.protected_io import install
from climate_palettes import field_rgb,anchors,DEFINITION
from ocean_flow_rendering import resample,particles,arrows
from reference_map_rendering import write_float32_geotiff
from climate_flow_texture import lic_luminance,particle_intensity

FIELDS={
 'surface_wind':('Surface wind', 'm/s', 'vector', .10, 'wind/speed/{s}/s'),
 '850_wind':('850 hPa wind', 'm/s', 'vector', .06, 'wind/speed/{s}/850'),
 '500_wind':('500 hPa wind', 'm/s', 'vector', .06, 'wind/speed/{s}/u'),
 'pressure':('Sea-level pressure anomaly', 'hPa', 'scalar', .12, 'pressure/anomaly/{s}/slp'),
 'rain':('Seasonal rainfall · ERA5', 'mm/day', 'scalar', .10, 'rain/{s}/s/era5'),
 'annual_rain':('Annual rainfall · IMERG', 'mm/year', 'scalar', .15, 'rain/annual/s'),
 'water':('Column water', 'kg/m²', 'scalar', .06, 'moisture/water/{s}/column'),
 'moisture_flux':('Column moisture transport', 'kg/(m·s)', 'vector', .09, 'moisture/flux/{s}/column'),
 'convergence':('Moisture convergence', 'mm/day', 'scalar', .05, 'moisture/convergence/{s}/column/era5_native'),
 'sst_oisst':('Bulk SST · OISST', '°C', 'scalar', .11, 'sea_temp/{s}/s'),
 'sst':('Bulk SST · ERA5', '°C', 'scalar', 0, 'sea_temp/{s}/s/era5'),
 'current':('Ocean currents · GLORYS', 'm/s', 'vector', .10, 'ocean/current/{s}/s'),
 'storage_depth':('Thermal storage / density mixed-layer depth', 'm', 'scalar', 0, 'ocean/mixed_layer_depth/{s}/s'),
 'skin_sst':('Ocean / ice skin temperature', '°C', 'scalar', 0, 'ocean/skin_temperature/{s}/s'),
 'evaporation':('Ocean latent heat loss', 'W/m²', 'scalar', 0, 'energy/latent_upward/{s}/s'),
 'surface_radiation':('Net surface radiation', 'W/m²', 'scalar', 0, 'energy/net_radiation/{s}/s')}
MONTH_DAYS=np.array([31,28.25,31,30,31,30,31,31,30,31,30,31])
SOURCES={
 'surface_wind':('era5_u10m','era5_v10m'), '850_wind':('era5_u850','era5_v850'),
 '500_wind':('era5_u500','era5_v500'), 'pressure':('era5_slp_anom',),
 'rain':('era5_pr',), 'water':('era5_tcwv',), 'moisture_flux':('era5_viwve','era5_viwvn'),
 'convergence':('era5_vimd',), 'sst':('era5_sst',), 'sst_oisst':('oisst_sst',),
 'current':('glorys_uo','glorys_vo'), 'storage_depth':('glorys_mlotst',),
 'evaporation':('era5_mer',), 'surface_radiation':('era5_msnswrf','era5_msnlwrf')}

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def relative(path,parent):return Path(os.path.relpath(path,parent)).as_posix()
def quarter(values, season):
    indices=np.arange(season*3,season*3+3)%12
    return np.average(values[indices],axis=0,weights=MONTH_DAYS[indices])
def magnitude(a):return np.hypot(*a) if a.ndim==3 else a
def remap(a,width):
    if a.shape[-1]==width:return np.asarray(a,dtype=float)
    if a.ndim==3:return np.stack([remap(x,width) for x in a])
    return conservative_remap(np.asarray(a,dtype=np.float32),width).astype(float)
def bundle(name):
    path=references.path(name)
    header,a=references.read(path,missing=np.nan,dtype=float)
    assert header.rows*2==header.columns and header.layers==12,path
    if name=='era5_pr':a/=MONTH_DAYS[:,None,None]
    return a
def annual_rain(width):
    p=references.path('imerg_prec_annual',annual=True)
    return remap(references.read(p,missing=np.nan,dtype=float)[1][0],width)

def resize_scalar(a,width):
    if a.shape[-1]==width:return np.array(a,dtype=float)
    valid=np.isfinite(a)
    def transform(x):
        # Longitude is periodic. Padding avoids clamping the two seam half-cells.
        if width % a.shape[-1] == 0:
            scale=width//a.shape[-1]
            padded=np.pad(np.asarray(x,dtype=np.float32),((0,0),(1,1)),mode='wrap')
            result=np.asarray(Image.fromarray(padded).resize((width+2*scale,width//2),Image.Resampling.BILINEAR),dtype=float)
            return result[:,scale:-scale]
        from reference_map_rendering import sample
        yy,xx=np.indices((width//2,width),dtype=float)
        return sample(np.asarray(x,dtype=float),(xx+.5)*a.shape[1]/width-.5,(yy+.5)*a.shape[0]/(width//2)-.5)
    fraction=transform(valid)
    value=transform(np.where(valid,a,0))/np.maximum(fraction,1e-20)
    value[fraction<.999]=np.nan
    return value

def flow_path(prefix,style,suffix='.png'):
    prefix=Path(prefix)
    # .../ocean/current/season/layer/id -> .../ocean/current_STYLE/season/layer/id
    quantity=prefix.parents[2].name
    target=style if quantity=='speed' else 'current_'+style
    path=prefix.parents[2].with_name(target)/prefix.parents[1].name/prefix.parent.name/(prefix.name+suffix)
    if suffix=='.tif':
        reference_maps=ROOT/'refs/processed/climate/maps'
        if path.is_relative_to(reference_maps):path=ROOT/'refs/processed/climate/fields'/path.relative_to(reference_maps)
        else:
            parts=list(path.parts);index=parts.index('maps');parts[index]='fields';path=Path(*parts)
    return path

def flow_complete(prefix,adopt_legacy=True):
    prefix=Path(prefix)
    try:
        for style in ('lic','particles'):
            for suffix in ('.png','.tif'):
                destination=flow_path(prefix,style,suffix)
                old=Path(str(prefix)+f'-{style}'+suffix)
                if not destination.exists() and adopt_legacy and old.exists():
                    # These are unsealed products generated by this renderer. Copied, not moved:
                    # refs/ forbids moving a file away, so the legacy file stays.
                    roots=[ROOT/'refs/processed/climate']
                    reports=ROOT/'runs/reports'
                    if prefix.is_relative_to(reports):
                        roots.append(reports/prefix.relative_to(reports).parts[0])
                    assert any(old.resolve().is_relative_to(r.resolve()) and destination.resolve().is_relative_to(r.resolve()) for r in roots)
                    install(old,destination)
                with Image.open(destination) as im:
                    if im.size!=(2048,1024):return False
                    im.verify()
        return True
    except (OSError,ValueError):return False

def raster(field,value,mask,width=2048):
    if value.ndim==3:
        east,north,valid=resample(*value,mask,width)
        return np.array([east,north]),valid
    a=resize_scalar(np.where(mask,value,np.nan),width)
    return a,np.isfinite(a)

def save_map(path,field,value,mask,width=2048):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    data,valid=raster(field,value,mask,width)
    im=Image.fromarray(field_rgb(magnitude(data),field))
    if data.ndim==3:arrows(im,*data,valid)
    im.save(path)
    return data

def compose_flow(data,texture,field,style):
    """Keep outlines separate so a page toggle removes them completely."""
    east,north=data;valid=np.isfinite(east)&np.isfinite(north)
    base=field_rgb(np.hypot(east,north),field).astype(float)
    light=np.nan_to_num(texture)[...,None]
    rgb=base*(.65+.7*light) if style=='lic' else base*.7+(255-base*.7)*light
    rgb[~valid]=0
    im=Image.fromarray(np.clip(rgb,0,255).astype(np.uint8));arrows(im,east,north,valid)
    return im

def save_wind_flow(prefix,field,value,mask,land,width=2048):
    data,valid=raster(field,value,mask,width);east,north=data
    seed=20260906
    textures={'lic':lic_luminance(east,north,seed),'particles':particle_intensity(east,north,seed)}
    for style,texture in textures.items():
        image=compose_flow(data,texture,field,style)
        path=flow_path(prefix,style);path.parent.mkdir(parents=True,exist_ok=True);image.save(path)
        write_float32_geotiff(flow_path(prefix,style,'.tif'),texture)
    return dict(field=field,seed=seed,particle_days=4.5,width=width)

def save_ocean_flow(prefix,value,mask,width=2048):
    prefix=Path(prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    data,valid=raster('current',value,mask,width);east,north=data
    light=lic_luminance(east,north,20260906)
    intensity,heads,steps=particles(east,north,20260906,45)
    for name,texture in [('lic',light),('particles',intensity)]:
        im=compose_flow(data,texture,'current',name)
        if name=='particles':
            draw=ImageDraw.Draw(im)
            for x,y in heads:draw.ellipse((x-1,y-1,x+1,y+1),fill=(255,246,184))
        path=flow_path(prefix,name);path.parent.mkdir(parents=True,exist_ok=True);im.save(path)
        write_float32_geotiff(flow_path(prefix,name,'.tif'),texture)
    return dict(days=45,seed=20260906,steps=steps,terminal_particles=len(heads),width=width)

def full_reference(field,s,cache):
    if field=='skin_sst':return None
    if field=='annual_rain':return annual_rain(2048)
    data=[]
    for name in SOURCES[field]:
        if name not in cache:cache[name]=bundle(name)
        source=cache[name]
        a=source[s*3] if field in ('surface_wind','850_wind','500_wind','pressure') else quarter(source,s)
        data.append(resize_scalar(a,2048))
    if field=='evaporation':return data[0]*2.5e6/86400
    if field=='surface_radiation':return sum(data)
    return np.array(data) if len(data)>1 else data[0]

def write_reference_index(root,records,audit,parent):
    extras=root/'wind-flow-manifest.json'
    records=dict(records)
    if extras.exists():
        records.update({key:relative(ROOT/item['path'],parent) for key,item in read_json(extras)['products'].items()})
    rows=[]
    for key,path in records.items():
        rows.append(f'<tr><td>{html.escape(key)}</td><td><a href="{relative(parent/path,root)}">PNG</a></td></tr>')
    legends=[]
    for field,(label,unit,*_) in FIELDS.items():
        lo,hi=DEFINITION['scales'][field];stops=anchors(lo,hi,lo<0)
        gradient=','.join(f'rgb({",".join(map(str,c))}) {(v-lo)/(hi-lo)*100:.4f}%' for v,c in stops)
        legends.append(f'<tr><td>{html.escape(label)} ({unit})</td><td><div style="width:400px;height:18px;background:linear-gradient(to right,{gradient})"></div>{lo:g} → {hi:g}</td></tr>')
    text='<html><meta charset="utf-8"><title>Standard climate references</title><style>body{font:16px system-ui;margin:2rem;background:#101923;color:#eee}a{color:#8bd1ff}td{padding:8px}</style><h1>Standard climate references</h1><p>Shared climate-colour-v2 palettes. Maps preserve geographic registration. LIC and particles are visual diagnostics, not scalar measurements. Zero is cyan on signed scales; black is missing data.</p><table>'+''.join(legends)+'</table><h2>Maps</h2><table>'+''.join(rows)+'</table><h2>Reference coverage</h2><ul>'+''.join(f'<li>{html.escape(k)}: {html.escape(v)}</li>' for k,v in audit.items())+'</ul><a href="manifest.json">Provenance and image hashes</a> · <a href="wind-flow-manifest.json">Wind texture provenance</a></html>'
    (root/'index.html').write_text(text,encoding='utf-8')
