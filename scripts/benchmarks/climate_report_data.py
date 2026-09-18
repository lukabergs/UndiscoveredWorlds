"""Numerical input adapter. Historical archive aliases never become UI identifiers."""
import json
from pathlib import Path
import struct
import sys
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/refs'))
from cell_grid import conservative_remap

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
REGIONS={
 'Global':[-90,90,-180,180], 'South Asian monsoon':[-10,35,35,110],
 'Central / eastern Pacific':[-25,25,-160,-80], 'Western Pacific':[-25,25,120,180],
 'North Atlantic / Greenland':[30,90,-85,45], 'North Pacific':[15,70,120,-110],
 'Central Asia':[15,65,35,120], 'Southern Ocean / Antarctica':[-90,-35,-180,180],
 'Tropical Atlantic / Africa':[-25,30,-65,55], 'Kuroshio':[15,45,115,160],
 'Agulhas':[-50,-15,10,50], 'South America':[-60,15,-95,-25]}
MONTH_DAYS=np.array([31,28.25,31,30,31,30,31,31,30,31,30,31])
SOURCES={
 'surface_wind':('era5_u10m','era5_v10m'), '850_wind':('era5_u850','era5_v850'),
 '500_wind':('era5_u500','era5_v500'), 'pressure':('era5_slp_anom',),
 'rain':('era5_pr',), 'water':('era5_tcwv',), 'moisture_flux':('era5_viwve','era5_viwvn'),
 'convergence':('era5_vimd',), 'sst':('era5_sst',), 'sst_oisst':('oisst_sst',),
 'current':('glorys_uo','glorys_vo'), 'storage_depth':('glorys_mlotst',),
 'evaporation':('era5_mer',), 'surface_radiation':('era5_msnswrf','era5_msnlwrf')}
CACHE_KEYS={'surface_wind':('u','v'),'850_wind':('u850','v850'), '500_wind':('u500','v500'),
 'pressure':('pressure',),'rain':('rain',),'water':('water',),'moisture_flux':('flux_u','flux_v'),
 'convergence':('convergence',),'sst':('sst',),'sst_oisst':('oisst',),'current':('current_u','current_v'),
 'evaporation':('evap',),'surface_radiation':('rad_sw','rad_lw')}

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def quarter(values, season):
    indices=np.arange(season*3,season*3+3)%12
    return np.average(values[indices],axis=0,weights=MONTH_DAYS[indices])
def magnitude(a):return np.hypot(*a) if a.ndim==3 else a
def weights(shape):
    lat=np.linspace(np.pi/2,-np.pi/2,shape[0]+1)
    return np.broadcast_to((np.sin(lat[:-1])-np.sin(lat[1:]))[:,None],shape)
def remap(a,width):
    if a.shape[-1]==width:return np.asarray(a,dtype=float)
    if a.ndim==3:return np.stack([remap(x,width) for x in a])
    return conservative_remap(np.asarray(a,dtype=np.float32),width).astype(float)
def table(path):
    x=np.genfromtxt(path,delimiter=',',names=True)
    w=len(np.unique(x['longitude']));h=len(np.unique(x['latitude']))
    return {n:x[n].reshape(4,h,w) for n in x.dtype.names}
def bundle(name):
    path=ROOT/f'refs/processed/climate/{name}_monthly.uwclim'
    with path.open('rb') as stream:
        magic=stream.read(8);version,w,h,n=struct.unpack('<IIII',stream.read(16))
    assert version==1 and h*2==w and n==12,(path,magic)
    a=np.array(np.memmap(path,dtype='<f4',mode='r',offset=40,shape=(n,h,w)),dtype=float)
    a[a < -9000]=np.nan
    if name=='era5_pr':a/=MONTH_DAYS[:,None,None]
    return a

class Archive:
    def __init__(self,path,ids):
        self.path=Path(path).absolute();self.entries={};self.ids=ids
        for p in self.path.glob('*-analysis.json'):
            d=read_json(p)
            if d['run'] in ids:self.entries[d['run']]=(p.name[:-len('-analysis.json')],d)
        missing=set(ids)-set(self.entries)
        if missing:raise ValueError(f'Missing analyzed numerical runs: {sorted(missing)}')
        self.refs=dict(np.load(self.path/'references.npz'))
        self.monthly_full={}
    def load(self,run):
        key,analysis=self.entries[run]
        fields=dict(np.load(self.path/f'{key}-fields.npz'))
        cmp=dict(np.load(self.path/f'{key}-comparison.npz'))
        masks={}
        for s,capture in enumerate(analysis['capture_paths']):
            raw=np.genfromtxt(ROOT/capture,delimiter=',',names=True)
            h,w=fields[f'pressure_{s}'].shape
            masks[f'850_{s}']=raw['low_level_850_available'].reshape(h,w)>.5
            masks[f'wet_{s}']=raw['land_fraction'].reshape(h,w)<.1
        storage=table(ROOT/f'runs/diagnostics/climate/{run}/maps/ocean_storage_fields.csv')
        for s in range(4):fields[f'storage_depth_{s}']=storage['mean_storage_depth_m'][s]
        with Image.open(ROOT/f'runs/fields/climate/rain/annual/s/{run}.tif') as im:
            fields['annual_rain_0']=np.array(im,dtype=float)
        return fields,cmp['land'][0],masks
    def reference(self,field,s,width):
        if field=='skin_sst':return None
        if field=='annual_rain':
            p=ROOT/'refs/processed/climate/imerg_prec_annual.uwclim'
            with p.open('rb') as stream:stream.seek(12);w,h=struct.unpack('<II',stream.read(8))
            a=np.array(np.memmap(p,dtype='<f4',mode='r',offset=40,shape=(h,w)),dtype=float);a[a < -9000]=np.nan
            return remap(a,width)
        if field=='storage_depth':
            return remap(quarter(np.load(self.path/'mixed-layer-reference.npz')['metres'],s),width)
        result=[]
        for key in CACHE_KEYS[field]:
            source=self.refs[key]
            a=source[s*3] if field in ('surface_wind','850_wind','500_wind','pressure') else quarter(source,s)
            result.append(remap(a,width))
        if field=='evaporation':return result[0]*2.5e6/86400
        if field=='surface_radiation':return sum(result)
        return np.array(result) if FIELDS[field][2]=='vector' else result[0]
    def provenance(self):
        return read_json(self.path/'reference-provenance.json')
