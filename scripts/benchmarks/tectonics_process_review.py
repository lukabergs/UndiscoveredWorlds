# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Save initial terrain and early ocean-change diagnostics without producing HTML."""
import argparse
import json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from tectonics_experiment import ROOT, write_json, digest
from tectonics_report_data import load_bundle, diagnostics
from map_rendering import colourize


def expand(mask, steps=2):
    for _ in range(steps):
        mask = mask | np.roll(mask, 1, 0) | np.roll(mask, -1, 0) | np.roll(mask, 1, 1) | np.roll(mask, -1, 1)
    return mask


def stats(values):
    return None if not values.size else dict(cells=int(values.size), mean=float(values.mean()),
        p10=float(np.quantile(values,.1)), median=float(np.median(values)), p90=float(np.quantile(values,.9)))


def review(run):
    record=json.loads((ROOT/f'runs/manifests/tectonics/{run}.json').read_text())
    if record['status']!='complete':
        raise ValueError('Process review requires a completed run')
    output=ROOT/f'runs/diagnostics/tectonics/{run}'
    palette=json.loads((ROOT/'configs/tectonics-palettes.json').read_text())['fields']['elevation']['stops']
    result={'run':run,'source_sha256':digest(Path(__file__)),'snapshots':{}}
    panels=[]
    initial=None
    early=next(s['label'] for s in record['snapshots'] if s['update_count']>0)
    for snap in record['snapshots']:
        if snap['label'] not in ('initial',early,'final'):
            continue
        manifest,f,_=load_bundle(ROOT/f'runs/checkpoints/tectonics/{run}'/snap['directory'],record['field_sha256'][snap['label']])
        z=f['elevation'];ocean=z<0
        m=diagnostics(manifest,f)
        # The native boundary classification is a proxy, not the stored flexural load.
        near=expand(f['boundary_type']==1)
        neighbour=(np.roll(z,1,0)+np.roll(z,-1,0)+np.roll(z,1,1)+np.roll(z,-1,1))/4
        item={'time_myr':manifest['time_myr'],'land_fraction':m['land_fraction'],
            'ocean_depth_m':m['ocean_depth_m'],'land_elevation_m':m['land_elevation_m'],
            'elevation_area_bands':m['elevation_area_bands'],'terrain_detail':m['terrain_detail'],
            'ocean_geologic_age_myr':stats(f['crust_age_myr'][ocean]),
            'convergent_ocean_depression_m':stats((neighbour-z)[ocean & (f['boundary_type']==1)]),
            'deep_ocean_cells':int((z<-6000).sum()),
            'deep_ocean_near_convergence_cells':int(((z<-6000)&near).sum())}
        if initial is None:
            initial={k:v.copy() for k,v in f.items() if k in ('elevation','boundary_type')}
        elif snap['label']==early:
            remote=ocean & (initial['elevation']<0) & ~near & ~expand(initial['boundary_type']==1)
            delta=z-initial['elevation']
            item['early_remote_ocean_eulerian_delta_m']=stats(delta[remote])
            item['early_remote_ocean_lowered_fraction']=float((delta[remote]<-1).mean()) if remote.any() else None
            opening_land=(initial['elevation']>0)&expand(initial['boundary_type']==2)
            item['early_opening_land_eulerian_delta_m']=stats(delta[opening_land])
            item['early_opening_land_absolute_delta_m']=stats(abs(delta[opening_land]))
        result['snapshots'][snap['label']]=item
        image=Image.fromarray(colourize(z,palette)).resize((1024,512),Image.Resampling.NEAREST)
        panels.append((f'Run {run}: {snap["label"]}, {manifest["time_myr"]:g} Myr',image))
    preview=Image.new('RGB',(1024,552*len(panels)),'#101b28');draw=ImageDraw.Draw(preview)
    for i,(label,image) in enumerate(panels):
        draw.text((12,i*552+12),label,fill='white');preview.paste(image,(0,i*552+40))
    output.mkdir(parents=True,exist_ok=True)
    preview.save(output/'process-preview.png')
    result['interpretation']='Early changes are at fixed raster locations and include advection, cooling and surface processes; they do not isolate cooling. Boundary depression is relative to four neighbours, not a measured trench depth.'
    write_json(output/'process-review.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('runs',type=int,nargs='+')
    for run in parser.parse_args().runs:
        review(run)
