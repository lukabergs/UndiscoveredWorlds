# /// script
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Couple an archived tectonic endpoint to climate and conservative surface erosion."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile
import numpy as np
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths
from tectonics_experiment import write_json, digest
from tectonics_report_data import load_bundle
from tectonics_batch_analysis import distribution
from map_rendering import colourize
from climate_palettes import anchors, DEFINITION as CLIMATE_PALETTE
from run_storage import ensure_archive_directory

ROOT=Path(__file__).resolve().parents[2]


def raster(table,field,shape,season=None):
    rows=table if season is None else table[table['season']==season]
    if 'x' in rows.dtype.names:
        result=np.full(shape,np.nan)
        result[rows['y'].astype(int),rows['x'].astype(int)]=rows[field]
    else:
        lats=np.unique(rows['latitude'])[::-1];lons=np.unique(rows['longitude'])
        result=np.full((len(lats),len(lons)),np.nan)
        ys=np.searchsorted(-lats,-rows['latitude']);xs=np.searchsorted(lons,rows['longitude'])
        result[ys,xs]=rows[field]
        # Cell-centred climate grids share the same world extent. Bilinear
        # interpolation supplies smooth forcing, never replaces terrain detail.
        if result.shape!=shape:
            result=np.asarray(Image.fromarray(result.astype('float32')).resize(shape[::-1],Image.Resampling.BILINEAR),dtype=float)
    if not np.isfinite(result).all(): raise ValueError(f'Incomplete climate field {field}')
    return result


def preview(path,values,stops):
    Image.fromarray(colourize(values,stops)).resize((1024,512),Image.Resampling.NEAREST).save(path)


def spherical_weights(shape):
    edges=np.linspace(np.pi/2,-np.pi/2,shape[0]+1)
    return np.broadcast_to((np.sin(edges[:-1])-np.sin(edges[1:]))[:,None],shape)


def climate_diagnostics(folder):
    result={}
    for name in ['climate_coupling','climate_energy_budget','climate_water_budget_area_weighted']:
        paths=list(folder.glob(f'validation/seed_*/{name}.csv'))
        if not paths:
            result[name]=dict(available=False)
            continue
        with paths[0].open() as stream:rows=list(csv.DictReader(stream))
        result[name]=dict(available=bool(rows),rows=rows)
    write_json(folder/'physics-diagnostics.json',result)
    return result


def climate(folder,z,exe,seed):
    folder.mkdir(parents=True)
    h,w=z.shape
    if np.max(abs(z))>23999: raise ValueError('Terrain exceeds existing app import datum; do not silently clip')
    land=np.where(z>0,np.maximum(1,np.rint(z)),0).astype('uint16')
    sea=np.where(z<=0,np.maximum(1,np.rint(-z)),0).astype('uint16')
    Image.fromarray(land).save(folder/'land.tif');Image.fromarray(sea).save(folder/'sea.tif')
    imported=land.astype(float)-sea.astype(float)
    env=os.environ.copy()
    for key in ['UW_CLIMATE_SURFACE_DRAG_INPUT','UW_CLIMATE_EXPERIMENT_INPUT']:env.pop(key,None)
    env.update(UW_CLIMATE_BENCHMARK_RUN_LOG=str(folder/'registry.json'),
        UW_CLIMATE_BENCHMARK_IMAGE_DIR=str(folder/'maps'),
        UW_CLIMATE_BENCHMARK_RUN_DIR=str(folder/'diagnostics'),
        UW_CLIMATE_VALIDATION_DIR=str(folder/'validation'),
        UW_REFERENCE_CLIMATE_DIR=str(folder/'no-observed-forcing'),
        UW_REFERENCE_PRECIP_GRID=str(folder/'no-observed-forcing/none.csv'))
    maps='koppen,temp,precip,precip-tif,'+','.join(f'{s}-{f}' for s in ['jan','apr','jul','oct'] for f in ['s-wind-speed','sst','ocean-current'])
    command=[str(exe),'generate','--climate-only','--seed',str(seed),'--resolution',str(w),
        '--import-land',str(folder/'land.tif'),'--import-sea',str(folder/'sea.tif'),
        '--climate-resolution','128','--climate-hydrology-resolution',str(w),
        '--climate-coupling-iterations','8','--climate-temperature-mode','calibrated',
        '--climate-mode','baseline','--climate-tropical-closure','legacy',
        '--climate-transport','ffsl','--climate-ocean-heat-workers','16',
        '--climate-ocean-circulation-solver','krylov','--climate-ocean-coastal-scheme','face-forcing',
        '--climate-atmosphere-solver','zonal',
        '--no-rivers','--no-lakes','--no-deltas','--no-profile-log','--map',maps,
        '--benchmark-info','Procedural tectonic world; climate-only feedback; no observed spatial forcing']
    start=time.monotonic()
    print(f'Climate solve: {folder.name}',flush=True)
    with (folder/'execution.log').open('x') as log:
        subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    elapsed=time.monotonic()-start
    cells=np.genfromtxt(folder/'diagnostics/1/climate_comparison_cells.csv',delimiter=',',names=True)
    thermal_path=next(folder.glob('validation/seed_*/circulation/run_*/climate_surface_thermal_fields.csv'))
    thermal=np.genfromtxt(thermal_path,delimiter=',',names=True)
    fields={key:[] for key in ['rain','temperature','wind','snow','ice']}
    for season in range(4):
        elevation=raster(cells,'elevation_m',z.shape,season)
        if not np.array_equal(elevation,imported):
            raise ValueError(f'Climate terrain differs from imported terrain by {np.max(abs(elevation-imported))} m')
        fields['rain'].append(12*raster(cells,'precipitation_mm_month',z.shape,season))
        fields['temperature'].append(raster(cells,'temperature_c',z.shape,season))
        u=raster(cells,'surface_u_mps',z.shape,season);v=raster(cells,'surface_v_north_mps',z.shape,season)
        fields['wind'].append(np.hypot(u,v))
        fields['snow'].append(raster(thermal,'land_snow_cover',z.shape,season))
        fields['ice'].append(raster(thermal,'permanent_land_ice_cover',z.shape,season))
    seasonal={key:np.asarray(value) for key,value in fields.items()}
    fields={key:value.mean(axis=0) for key,value in seasonal.items()}
    np.savez_compressed(folder/'seasonal-fields.npz',**seasonal)
    np.savez_compressed(folder/'annual-fields.npz',**fields,elevation=z,imported_elevation=imported)
    diagnostics=climate_diagnostics(folder)
    weights=spherical_weights(z.shape);land=z>0
    forcing_stats={key:dict(min=float(value.min()),max=float(value.max()),
        area_weighted_mean=float(np.sum(value*weights)/weights.sum()),
        area_weighted_land_mean=float(np.sum(value[land]*weights[land])/weights[land].sum()))
        for key,value in fields.items()}
    write_json(folder/'execution.json',dict(command=command,elapsed_seconds=elapsed,
        environment={k:v for k,v in env.items() if k.startswith('UW_')},terrain_preserved=True,
        maximum_import_rounding_m=float(np.max(abs(imported-z))),climate_reference_forcing=False,
        forcing_statistics=forcing_stats,
        coupling_converged=(diagnostics['climate_coupling']['rows'][-1]['converged']=='1'
            if diagnostics['climate_coupling'].get('rows') else None)))
    print(f'Climate complete: {folder.name}, {elapsed:.1f} s',flush=True)
    return fields,elapsed


def erosion(folder,z,fields,exe,scenario,duration):
    inputs=folder/'input';inputs.mkdir(parents=True)
    for key,value in dict(fields,elevation=z).items():np.asarray(value,dtype='<f8').tofile(inputs/f'{key}.f64')
    command=[str(exe),str(inputs),str(folder/'output'),str(z.shape[1]),str(z.shape[0]),str(duration),
        str(scenario['river-erosion']),str(scenario.get('wave-erosion-m-myr',0)),str(scenario['sediment-travel-km'])]
    start=time.monotonic()
    print(f'Surface erosion: {folder.name}, {duration} Myr',flush=True)
    with (folder/'execution.log').open('x') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    elapsed=time.monotonic()-start
    result=np.fromfile(folder/'output/elevation.f64',dtype='<f8').reshape(z.shape)
    write_json(folder/'execution.json',dict(command=command,elapsed_seconds=elapsed))
    return result,elapsed


def run(parent,run_id,cycles,duration):
    record=json.loads((ROOT/f'runs/manifests/tectonics/{parent}.json').read_text())
    assert record['status']=='complete'
    registry_path=ROOT/'runs/registry/tectonics.json';registry=json.loads(registry_path.read_text())
    assert max(item['id'] for item in registry['runs'])+1==run_id
    archive=ensure_archive_directory(ROOT/'runs/reports/tectonics')/str(run_id)
    if archive.exists():raise ValueError('Archive exists')
    archive.mkdir()
    snap=next(s for s in record['snapshots'] if s['label']=='final')
    parent_state,fields,_=load_bundle(ROOT/f'runs/checkpoints/tectonics/{parent}'/snap['directory'],record['field_sha256']['final'])
    z=fields['elevation'].copy();initial=z.copy();uniform=z.copy()
    manifest=dict(id=run_id,status='running',kind='tectonics-climate-feedback/v1',parent=parent,
        started_utc=datetime.now(timezone.utc).isoformat(),recipe=record['recipe'],
        tectonic_state_reused=True,parent_final_hashes=record['field_sha256']['final'],
        cycles=cycles,surface_interval_myr=duration,climate_reference_forcing=False,
        tectonic_endpoint_time_myr=parent_state['time_myr'],surface_experiment_duration_myr=cycles*duration,
        geometry='Parent torus raster interpreted as cell-centred longitude/latitude for climate and spherical surface erosion',
        limitations=['Frozen plates during surface experiment','No isostatic rebound','Regional glacier/wind closures; not a resolved ice-flow model'],stages=[])
    manifest_path=ROOT/f'runs/manifests/tectonics/{run_id}.json'
    registry['runs'].append(dict(id=run_id,status='running',kind=manifest['kind'],manifest=manifest_path.relative_to(ROOT).as_posix()))
    write_json(registry_path,registry);write_json(manifest_path,manifest);write_json(archive/'recipe.json',record['recipe'])
    binaries=archive/'binaries';binaries.mkdir()
    app=paths.executable('world_pipeline')
    for path in [app,*app.parent.glob('*.dll'),paths.executable('tectonic_climate_erosion')]:
        shutil.copy2(path,binaries/path.name)
    manifest['binary_sha256']={p.name:digest(p) for p in binaries.iterdir()}
    with zipfile.ZipFile(archive/'source.zip','w',zipfile.ZIP_DEFLATED) as package:
        for name in ['CMakeLists.txt','CMakePresets.json']:
            package.write(ROOT/name,name)
        for folder in ['src','tools','scripts/benchmarks','configs']:
            for path in sorted((ROOT/folder).rglob('*')):
                if path.is_file() and path.suffix in {'.cpp','.hpp','.h','.py','.json','.js','.html','.txt'}:
                    package.write(path,path.relative_to(ROOT).as_posix())
    manifest['source_sha256']=digest(archive/'source.zip')
    palette=json.loads((ROOT/'configs/tectonics-palettes.json').read_text())['fields']['elevation']['stops']
    preview(archive/'initial-elevation.png',z,palette)
    start=time.monotonic()
    try:
        for cycle in range(cycles+1):
            climate_fields,seconds=climate(archive/f'climate-{cycle}',z,binaries/app.name,record['recipe']['seed'])
            stage=dict(cycle=cycle,climate_seconds=seconds)
            if cycle<cycles:
                before=z.copy()
                z,seconds=erosion(archive/f'erosion-{cycle}',z,climate_fields,binaries/'tectonic_climate_erosion.exe',record['recipe']['scenario'],duration)
                stage['erosion_seconds']=seconds
                constant={key:np.full(z.shape,value) for key,value in dict(rain=1000,temperature=15,wind=8,snow=0,ice=0).items()}
                uniform,seconds=erosion(archive/f'uniform-control-{cycle}',uniform,constant,binaries/'tectonic_climate_erosion.exe',record['recipe']['scenario'],duration)
                stage['uniform_seconds']=seconds
                stage['change_rms_m']=float(np.sqrt(np.mean((z-before)**2)))
                preview(archive/f'elevation-{cycle+1}.png',z,palette)
            manifest['stages'].append(stage);write_json(manifest_path,manifest)
        np.savez_compressed(archive/'comparison-fields.npz',initial=initial,final=z,uniform_control=uniform,climate_effect=z-uniform,total_change=z-initial)
        preview(archive/'final-elevation.png',z,palette);preview(archive/'uniform-control.png',uniform,palette)
        stops=[(-500,[35,70,150]),(0,[238,237,220]),(500,[150,48,40])]
        preview(archive/'climate-effect.png',z-uniform,stops)
        for key,values,stops in [('rain',climate_fields['rain'],anchors(*CLIMATE_PALETTE['scales']['annual_rain'])),
                ('wind',climate_fields['wind'],anchors(*CLIMATE_PALETTE['scales']['surface_wind'])),
                ('temperature',climate_fields['temperature'],anchors(*CLIMATE_PALETTE['scales']['air_temperature'],signed=True)),
                ('snow',climate_fields['snow'],[(0,[54,72,61]),(1,[240,247,255])]),
                ('ice',climate_fields['ice'],[(0,[54,72,61]),(1,[195,236,255])])]:
            preview(archive/f'final-{key}.png',values,stops)
            manifest.setdefault('map_palettes',{})[key]=stops
        manifest.update(status='complete',elapsed_seconds=time.monotonic()-start,
            fraction_below={k:distribution(v) for k,v in [('initial',initial),('final',z),('uniform',uniform)]},
            fraction_below_spherical={k:distribution(v,spherical_weights(z.shape))
                for k,v in [('initial',initial),('final',z),('uniform',uniform)]},
            climate_effect_rms_m=float(np.sqrt(np.mean((z-uniform)**2))),maximum_total_change_m=float(np.max(abs(z-initial))),
            fields_sha256=digest(archive/'comparison-fields.npz'))
    except Exception as error:
        manifest.update(status='failed',error=str(error),elapsed_seconds=time.monotonic()-start);raise
    finally:
        write_json(manifest_path,manifest);registry['runs'][-1]['status']=manifest['status'];write_json(registry_path,registry)
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',type=int,required=True);parser.add_argument('--run',type=int,required=True)
    parser.add_argument('--cycles',type=int,default=3);parser.add_argument('--interval-myr',type=float,default=10)
    args=parser.parse_args();run(args.parent,args.run,args.cycles,args.interval_myr)
