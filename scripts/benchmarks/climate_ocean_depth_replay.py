# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Frozen momentum sensitivity, not an additional coupled benchmark run."""
import argparse,json,hashlib,shutil,subprocess
from pathlib import Path
import numpy as np
from climate_report_data import Archive,ROOT,magnitude,weights
from climate_report_metrics import evaluate

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True);p.add_argument('--run',type=int,required=True)
    args=p.parse_args();archive=args.archive.absolute();out=archive/f'{args.run}-depth-replay';out.mkdir(exist_ok=True)
    assert not (out/'receipt.json').exists(), 'Completed replay evidence is immutable'
    exe=out/'climate_ocean_dynamics_tests.exe';shutil.copy2(ROOT/'out/build/x64-Debug/tests/Release'/exe.name,exe)
    a=Archive(archive,[args.run]);fields,land,masks=a.load(args.run)
    sources={}
    for path in (archive/str(args.run)/'ocean-capture').glob('ocean-*.txt'):
        with path.open() as f:s=int(f.readline().split()[-1])
        if s not in sources or int(path.stem.split('-')[-1])>int(sources[s].stem.split('-')[-1]):sources[s]=path
    results={};regions={'Global':(-90,90,-180,180),'Southern Ocean':(-65,-40,-180,180),
        'Kuroshio':(20,40,125,150),'Agulhas':(-45,-20,15,40),'Equatorial Pacific':(-5,5,-160,-90),
        'Arabian Sea':(5,20,50,75)}
    for s,source in sorted(sources.items()):
        lines=source.read_text().splitlines();parts=lines[0].split();w,h=int(parts[1]),int(parts[2])
        first=3 if parts[0]=='UW_OCEAN_V2' else 2
        depth_start=first+int(lines[first])+2
        assert int(lines[depth_start-1])==w*h
        original=np.asarray(lines[depth_start:depth_start+w*h],dtype=float)
        for cap in (300,150,75):
            modified=list(lines);modified[depth_start:depth_start+w*h]=[format(x,'.9g') for x in np.minimum(original,cap)]
            # Momentum-only replay: a shallower active layer cannot accept the
            # old deep thermal-storage state. Disable that state for every
            # case; the one-second heat result is discarded, never compared.
            if parts[0]=='UW_OCEAN_V2':
                storage=modified[2].split();storage[0]='0';modified[2]=' '.join(storage)
            fixture=out/f'{s}-{cap}.txt';fixture.write_text('\n'.join(modified)+'\n')
            csv=out/f'{s}-{cap}.csv';command=[str(exe),'--momentum-replay',str(fixture),str(csv),'2e-6','1.5e-6']
            run=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            if run.returncode:
                (out/f'{s}-{cap}-failure.json').write_text(json.dumps(dict(command=command,code=run.returncode,stdout=run.stdout,stderr=run.stderr),indent=2))
                raise RuntimeError(f'Frozen replay {s}/{cap} failed; preserved receipt')
            budget=json.loads(run.stdout);raw=np.genfromtxt(csv,delimiter=',',names=True)
            v=np.array([raw['east_mps'].reshape(h,w),raw['north_mps'].reshape(h,w)])
            assert np.isfinite(v).all() and budget['volume_divergence_mps']<1e-14 and abs(budget['heat_relative_residual'])<1e-8
            if cap==300:
                difference=float(np.max(np.abs(v-fields[f'current_{s}'])))
                assert difference<2e-6,difference
                budget['native_reproduction_max_error_mps']=difference
            ref=a.reference('current',s,w);lat=90-(np.arange(h)+.5)*180/h;lon=-180+(np.arange(w)+.5)*360/w
            metrics={}
            for name,(south,north,west,east) in regions.items():
                mask=masks[f'wet_{s}']&np.isfinite(magnitude(ref))&(lat[:,None]>=south)&(lat[:,None]<=north)&(lon[None,:]>=west)&(lon[None,:]<=east)
                metrics[name]=evaluate(v,ref,mask,'current')
            results[f'{s}/{cap}']=dict(budget=budget,metrics=metrics,command=command,source=str(source),source_sha256=sha(source),fixture_sha256=sha(fixture),output_sha256=sha(csv))
            print(s,cap,{name:round(data['loss'],4) for name,data in metrics.items()},flush=True)
    receipt=dict(method=__doc__,run=args.run,executable_sha256=sha(exe),
        limitations='Same frozen winds, temperature and coast mask; active dynamical depth is capped and variable thermal storage is disabled in all cases. Heat response is discarded. Shallower dynamical bathymetry would also cap thermal storage in a full run, so these are response diagnostics, not evidence of improved coupled climate.',results=results)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
