# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Render reference wind textures with the same kernels as standard comparisons."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/benchmarks'))
from climate_report_data import bundle
from climate_palettes import DEFINITION,VERSION
from climate_report_render import save_wind_flow,flow_path,flow_complete

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--reuse-maps',action='store_true');args=parser.parse_args()
    base=ROOT/'refs/processed/climate';out=base/'maps/standard-v2';records={};inputs={}
    with Image.open(ROOT/'refs/processed/land_ocean/maps/earth_land_ocean_2048x1024.png') as im:land=np.array(im.convert('L'))>127
    for field,layer,east_name,north_name in [('surface_wind','s','era5_u10m','era5_v10m'),('850_wind','850','era5_u850','era5_v850'),('500_wind','u','era5_u500','era5_v500')]:
        east,north=bundle(east_name),bundle(north_name)
        for name in (east_name,north_name):
            p=base/f'{name}_monthly.uwclim';inputs[p.relative_to(ROOT).as_posix()]=sha(p)
        for season in range(1,5):
            prefix=out/f'climate/wind/speed/{season}/{layer}/2048'
            value=np.array([east[(season-1)*3],north[(season-1)*3]])
            if not args.reuse_maps or not flow_complete(prefix):save_wind_flow(prefix,field,value,np.isfinite(value).all(axis=0),land)
            for style in ('lic','particles'):
                source=flow_path(prefix,style,'.tif');path=flow_path(prefix,style)
                key=f'{field}_{style}/{season-1}'
                records[key]=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path),numeric_texture=source.relative_to(ROOT).as_posix(),numeric_texture_sha256=sha(source))
        print(field,'reference textures rendered',flush=True)
    (out/'wind-flow-manifest.json').write_text(json.dumps(dict(palette_version=VERSION,inputs_sha256=inputs,products=records,
        seed=20260906,particle_days=4.5,method='Monthly ERA5 vectors; shared LIC and 4.5-day particle kernels and seed with standard model comparisons. Shared palettes and fixed layer speed limits.'),indent=2)+'\n')
if __name__=='__main__':main()
