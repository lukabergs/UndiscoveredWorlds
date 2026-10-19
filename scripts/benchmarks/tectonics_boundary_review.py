"""Boundary activity evidence from saved fields; never replays a simulation."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import grid, runs
from uwtools.tectonics import bundle_files, is_run_bundle, snapshot_manifest
from tectonics_report_data import elevation_m, read_field


def review(run):
    record=runs.tectonics(run)
    frames=[]
    for snap in record['snapshots']:
        path=snap['path']
        m=snapshot_manifest(path);files=bundle_files(path);native=is_run_bundle(path)
        if frames and frames[-1]['time']==m['time_myr']:continue
        owner=read_field(path/files['plate_id'],'plate_id',(m['height'],m['width']),native)
        stable=np.asarray(m['topology_diagnostics']['stable_ids'])[owner]
        scenario=record['recipe']['scenario'];value=lambda key:scenario[key] if key in scenario else scenario[key.replace('_','-')]  # records before the snake_case keys
        if m['plates'] and 'velocity_east_north_cm_yr' in m['plates'][0]:
            # Contract schema 8: each plate's centroid velocity (east, north) in cm/yr. On a cube
            # snapshot the rasters are the id-order view (N x 6N), so the rolls pair cells within a face row.
            v=np.asarray([p['velocity_east_north_cm_yr'] for p in m['plates']])[owner]
        else:
            width=next((scenario[k] for k in ('world_width_km','world-width-km') if k in scenario),None)
            cell_m=width*1000/m['width'] if width is not None else grid.cell_km(m.get('run_scenario'),m['width'])*1000
            scale=cell_m/value('delta_time_myr')/1e4
            v=np.asarray([p['linear_velocity'] for p in m['plates']])[owner]*scale
        a=np.stack([stable,stable]);b=np.stack([np.roll(stable,-1,1),np.roll(stable,-1,0)])
        closing=np.stack([v[...,0]-np.roll(v[...,0],-1,1),v[...,1]-np.roll(v[...,1],-1,0)])
        if not frames:
            z=elevation_m(read_field(path/files['heightmap'],'heightmap',owner.shape,native),m['sea_level_m'])
            land=np.stack([z>0,z>0])&np.stack([np.roll(z,-1,1)>0,np.roll(z,-1,0)>0])
            initial=dict(boundary_faces=int((a!=b).sum()),land_land_faces=int(((a!=b)&land).sum()),
                strong_converging_land_faces=int(((a!=b)&land&(closing>.5)).sum()))
        frames.append(dict(time=m['time_myr'],a=a,b=b,closing=closing))
    sampled=[]
    for threshold in (0,.1,.5):
        persistent=flips=triplets=roundtrips=0
        for old,new in zip(frames,frames[1:]):
            mask=(old['a']!=old['b'])&(old['a']==new['a'])&(old['b']==new['b'])
            mask&=(abs(old['closing'])>threshold)&(abs(new['closing'])>threshold)
            persistent+=int(mask.sum());flips+=int((mask&(old['closing']*new['closing']<0)).sum())
        for a,b,c in zip(frames,frames[1:],frames[2:]):
            mask=(a['a']!=a['b'])&(a['a']==b['a'])&(a['a']==c['a'])&(a['b']==b['b'])&(a['b']==c['b'])
            mask&=(abs(a['closing'])>threshold)&(abs(b['closing'])>threshold)&(abs(c['closing'])>threshold)
            triplets+=int(mask.sum());roundtrips+=int((mask&(a['closing']*b['closing']<0)&(b['closing']*c['closing']<0)).sum())
        sampled.append(dict(threshold_cm_yr=threshold,persistent_face_intervals=persistent,reversals=flips,
            reversal_fraction=flips/max(1,persistent),triplets=triplets,roundtrips=roundtrips))
    result=dict(run=run,initial=initial,sampled=sampled,step_diagnostics=m['motion_diagnostics'],
        method='Same fixed raster face with the same stable owner on each side; not material-boundary tracking. Sampled intervals are 10 Myr. Step counters, where available, use 1 Myr.')
    out=runs.output_directory('tectonics',run)/'boundary-review.json'
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('runs',type=int,nargs='+')
    for run in parser.parse_args().runs:review(run)
