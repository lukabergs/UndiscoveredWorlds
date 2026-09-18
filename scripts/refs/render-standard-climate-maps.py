# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Generate standalone standard climate reference PNGs without running a simulation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts/benchmarks'))
from climate_report_data import Archive,FIELDS,magnitude,SOURCES
from climate_report_render import full_reference,save_map,save_ocean_flow,flow_complete,flow_path
from climate_report import write_reference_index
from climate_palettes import VERSION

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reuse-maps',action='store_true');args=p.parse_args()
    out=ROOT/'refs/processed/climate/maps/standard-v2';parent=ROOT/'runs/reports';records={};coverage={};inputs={}
    # The annual-reference reader is independent of model/archive state.
    provider=Archive.__new__(Archive)
    for field,(_,_,_,_,semantic) in FIELDS.items():
        if field=='skin_sst':
            coverage[field]='Unavailable: no retained matching skin-temperature reference.';continue
        cache={}
        for s in ([0] if field=='annual_rain' else range(4)):
            value=full_reference(field,s,provider,cache)
            target=out/'climate'/semantic.format(s=s+1)/'2048.png'
            if not args.reuse_maps or not target.exists():save_map(target,field,value,np.isfinite(magnitude(value)))
            records[f'{field}/{s}']=Path(os.path.relpath(target,parent)).as_posix()
            if field=='current':
                prefix=target.with_suffix('')
                if not args.reuse_maps or not flow_complete(prefix):save_ocean_flow(prefix,value,np.isfinite(magnitude(value)))
                for style in ('lic','particles'):records[f'current_{style}/{s}']=Path(os.path.relpath(flow_path(prefix,style),parent)).as_posix()
        for name in SOURCES.get(field,()):
            source=ROOT/f'refs/processed/climate/{name}_monthly.uwclim';inputs[source.relative_to(ROOT).as_posix()]=sha(source)
        coverage[field]='GLORYS density-defined MLD; not a thermal storage observation.' if field=='storage_depth' else 'Available'
        print(field,flush=True)
    source=ROOT/'refs/processed/climate/imerg_prec_annual.uwclim';inputs[source.relative_to(ROOT).as_posix()]=sha(source)
    manifest=dict(palette_version=VERSION,input_sha256=inputs,reference_maps=records,coverage=coverage,
        method='Monthly representative wind/SLP; day-weighted quarters for other seasonal quantities. Annual IMERG. Shared palettes, fixed ranges, black nodata. Ocean particles: 45 days, seed 20260906.',
        images_sha256={str((parent/path).relative_to(ROOT)):sha(parent/path) for path in records.values()})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    write_reference_index(out,records,coverage,parent)
if __name__=='__main__':main()
