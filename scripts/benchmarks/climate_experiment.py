"""Run one archived climate experiment, analyze it, then require a decision before the next."""
import argparse
import datetime
import difflib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import zipfile
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths
from uwtools.hashing import sha256
from climate_report_data import Archive, FIELDS, ROOT, magnitude, table, weights
from climate_report_metrics import evaluate, rank
from climate_report import mask_for
from run_storage import ensure_archive_directory
from climate_benchmark import from_legacy_command, set_option

H = ROOT / 'runs/reports'
BUILD = paths.build_directory()
EXE = paths.executable('world_pipeline')

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')

def setup(archive, source_archive=None, baseline=None, count=10, skip_tests=False):
    selected = read(H/'climate-global-current-selection.json')
    baseline = selected['selected_run'] if baseline is None else baseline
    if source_archive is None:
        source_archive = (H/selected['source_manifest']).resolve().parent.parent
    ensure_archive_directory(archive)
    assert not (archive / 'batch.json').exists()
    source = Archive(source_archive, [baseline]); key = source.entries[baseline][0]
    for name in ('references.npz', 'reference-provenance.json', 'mixed-layer-reference.npz', 'mixed-layer-reference.json'):
        shutil.copy2(source_archive / name, archive / name)
    for suffix in ('fields.npz', 'comparison.npz', 'analysis.json'):
        shutil.copy2(source_archive / f'{key}-{suffix}', archive / f'{baseline}-{suffix}')
    (archive / str(baseline)).mkdir()
    for name in ('manifest.json', 'source.zip', 'decision.json', 'command.json'):
        shutil.copy2(source_archive / key / name, archive / str(baseline) / name)
    shutil.copy2(ROOT / 'runs/registry/climate.json', archive / 'registry-before.json')
    shutil.copy2(H / 'climate-global-current-selection.json', archive / 'selection-before.json')
    first = max(x['id'] for x in read(ROOT/'runs/registry/climate.json')['runs']) + 1
    write(archive / 'batch.json', dict(baseline=baseline, expected_runs=list(range(first, first+count)), skip_tests=skip_tests,
        pre_batch_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        intent=f'{count} serial adaptive experiments; fixed global ranking, distribution and regional numeric review. Visual assessment reserved for user. No imposed observed climate fields.',
        grids=dict(world=[1024,512], atmosphere=[128,64], ocean=[128,64], hydrology=[256,128], render=[2048,1024]),
        seed=1, coupling_updates=8, climate_reference_forcing=False))

def run(archive, run_id, plan):
    batch = read(archive/'batch.json'); expected = batch['expected_runs']; baseline = batch['baseline']
    assert run_id in expected
    first = expected[0]; parent = baseline if run_id == first else expected[expected.index(run_id)-1]
    if run_id != first:
        assert (archive / str(parent) / 'decision.json').is_file(), 'Analyze and decide before advancing'
    folder = archive / str(run_id)
    assert not folder.exists()
    registry = read(ROOT / 'runs/registry/climate.json')
    assert max(x['id'] for x in registry['runs']) + 1 == run_id
    with (archive / f'{run_id}-build.log').open('x') as stream:
        targets = ['world_pipeline'] if batch.get('skip_tests') else ['world_pipeline', 'world_tests']
        subprocess.run(['cmake','--build',str(BUILD),'--config','Release','--target',*targets,
            '--parallel','4'], cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
    if not batch.get('skip_tests'):
        with (archive / f'{run_id}-tests.log').open('x') as stream:
            subprocess.run(['ctest','--test-dir',str(BUILD),'-C','Release','--output-on-failure','-R',
                r'^Climate(Atmosphere|Hydrology|OceanDynamics|OceanHeat)\.'],
                cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
    folder.mkdir()
    paths = subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','src','tests','CMakeLists.txt','CMakePresets.json','configs','scripts/benchmarks'], cwd=ROOT, text=True).splitlines()
    paths.extend(p.relative_to(ROOT).as_posix() for p in (ROOT/'scripts/benchmarks').glob('climate_experiment*.py'))
    paths = sorted(set(p for p in paths if (ROOT / p).is_file()))
    with zipfile.ZipFile(folder / 'source.zip','w',zipfile.ZIP_DEFLATED) as z:
        for path in paths: z.write(ROOT / path, path)
    with zipfile.ZipFile(archive / str(parent) / 'source.zip') as previous:
        patches = []
        for path in paths:
            if not (path.startswith(('src/simulations/climate/', 'tests/simulations/climate/')) or path == 'src/wip/generation_tuning.hpp'): continue
            before = previous.read(path).decode('utf-8').replace('\r\n','\n') if path in previous.namelist() else ''
            after = (ROOT / path).read_text(encoding='utf-8')
            patches.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=f'{parent}/{path}',tofile=f'{run_id}/{path}'))
        (folder / 'changes-from-parent.patch').write_text(''.join(patches),encoding='utf-8')
    shutil.copy2(EXE, folder / EXE.name)
    for path in EXE.parent.glob('*.dll'): shutil.copy2(path,folder / path.name)
    inherited = read(archive / str(baseline) / 'manifest.json')
    tuning = (ROOT / 'src/wip/generation_tuning.hpp').read_text()
    env = os.environ.copy()
    env.pop('UW_CLIMATE_EXPERIMENT_INPUT',None)
    env['UW_CLIMATE_SURFACE_DRAG_INPUT'] = str(ROOT / 'refs/processed/climate/vegetation-drag-20260914/land-drag-contrast-1024.txt')
    inputs = inherited['inputs_sha256']
    stats = {path:dict(bytes=(ROOT/path).stat().st_size,mtime_ns=(ROOT/path).stat().st_mtime_ns) for path in inputs}
    verified = read(archive/str(first)/'manifest.json') if run_id != first else inherited
    previous_stats = verified.get('input_file_stats',{}) if verified else {}
    # Reuse baseline/first-case hashes only while both size and modification
    # time match the archived receipt; hash changed or unrecorded inputs.
    inputs = {path:(verified['inputs_sha256'][path] if verified and previous_stats.get(path)==stats[path]
                   else sha256(ROOT/path)) for path in inputs}
    manifest = dict(case=str(run_id),parent=str(parent),reason=plan['hypothesis'],utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={path:sha256(ROOT / path) for path in paths},source_zip_sha256=sha256(folder/'source.zip'),
        executable_sha256=sha256(folder/EXE.name),runtime_sha256={p.name:sha256(p) for p in folder.glob('*.dll')},
        inputs_sha256=inputs,input_file_stats=stats,input_hash_provenance='Reuse baseline/first-case hashes with matching size/mtime receipts; rehash changed or unrecorded inputs.',all_tuning_parameters=dict(re.findall(r'inline constexpr (?:float|double|bool|int) (\w+)\s*=\s*([^;]+);',tuning)),
        environment={'UW_CLIMATE_EXPERIMENT_INPUT':None,'UW_CLIMATE_SURFACE_DRAG_INPUT':env['UW_CLIMATE_SURFACE_DRAG_INPUT']},
        seed=batch['seed'],world=[1024,512],native_grids={'atmosphere':[128,64],'ocean':[128,64],'hydrology':[256,128]},
        coupling_updates=8,climate_reference_forcing=False,visual_review='User',tests_skipped=batch.get('skip_tests',False))
    write(folder/'manifest.json',manifest)
    write(archive/f'{run_id}-plan.json',dict(plan,case=str(run_id),parent=str(parent)))
    # Receipts from the retired PowerShell wrapper are converted to uw.py.
    command = from_legacy_command(read(archive/str(baseline)/'command.json'))
    if batch.get('skip_tests') and '--map' in command:
        # Flow textures are generated once by the shared report renderer.
        index = command.index('--map')+1
        command[index] = ','.join(x for x in command[index].split(',') if not x.endswith(('-lic','-part')))
    # Batches keep the seed recorded at setup; earlier batches retain 20260906.
    set_option(command,'--seed',batch['seed'])
    set_option(command,'--executable',folder/EXE.name)
    set_option(command,'--climate-atmosphere-capture',folder/'capture')
    set_option(command,'--climate-ocean-capture',folder/'ocean-capture')
    write(folder/'command.json',command)
    start=time.perf_counter()
    with (folder/'console.log').open('x') as stream:
        result=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
    ids=re.findall(r'Climate benchmark run ID: (\d+)',(folder/'console.log').read_text())
    receipt=dict(run_id=int(ids[-1]) if ids else None,parent=parent,exit_code=result.returncode,wrapper_seconds=time.perf_counter()-start)
    write(folder/'execution.json',receipt)
    assert result.returncode == 0 and receipt['run_id'] == run_id, receipt
    check_code = None
    if not batch.get('skip_tests'):
        with (folder/'checks.log').open('x') as stream:
            check=subprocess.run([sys.executable,'scripts/benchmarks/check_climate_baseline.py','--run-id',str(run_id)],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        check_code = check.returncode
    receipt['checks_exit']=check_code
    outputs=list((ROOT/f'runs/diagnostics/climate/{run_id}').rglob('*'))+[ROOT/f'runs/logs/climate/{run_id}.log']
    receipt['artifact_sha256']={p.relative_to(ROOT).as_posix():sha256(p) for p in outputs if p.is_file()}
    write(folder/'execution.json',receipt)
    assert check_code in (None,0), 'See preserved checks.log'
    analyze(archive,run_id)

def extract(archive,run_id,capture_directory=None):
    folder=archive/str(run_id)
    native=ROOT/f'runs/diagnostics/climate/{run_id}'
    raw=np.loadtxt(native/'climate_comparison_cells.csv',delimiter=',',skiprows=1,usecols=(5,12,13,14,15,16))
    low=raw.reshape(4,128,4,256,4,6).mean(axis=(2,4))
    cmp={name:low[...,i] for i,name in enumerate(('land','temperature','rain_month','pressure','u','v'))}
    upper=np.loadtxt(native/'climate_comparison_cells.csv',delimiter=',',skiprows=1,usecols=(17,18)).reshape(4,64,8,128,8,2).mean(axis=(2,4))
    process=table(native/'maps/climate_process_fields.csv')
    thermal=table(native/'maps/climate_surface_thermal_fields.csv')
    terrain=table(native/'maps/climate_terrain_fields.csv')
    ocean=table(native/'maps/ocean_fields.csv')
    fields={};captures=[]
    for s in range(4):
        capture=max((capture_directory or folder/'capture').glob(f'circulation-*-season-{s}.csv'),key=lambda p:int(p.name.split('-')[1]))
        captures.append(capture.relative_to(ROOT).as_posix())
        d=np.genfromtxt(capture,delimiter=',',names=True)
        c=lambda name:d[name].reshape(64,128)
        values=dict(surface_wind=np.array([cmp['u'][s],cmp['v'][s]]),
            **{'850_wind':np.array([c('low_level_850_east_mps'),-c('low_level_850_south_mps')]),
               '500_wind':np.array([upper[s,:,:,0],upper[s,:,:,1]])},
            pressure=c('solved_hpa'),rain=sum(terrain[n][s] for n in ('stratiform_mm_day','orographic_mm_day','convective_mm_day')),
            water=process['mean_column_water_mm'][s],
            moisture_flux=np.array([process['lower_flux_east_kg_m_s'][s]+process['upper_flux_east_kg_m_s'][s],-process['lower_flux_south_kg_m_s'][s]-process['upper_flux_south_kg_m_s'][s]]),
            convergence=process['transport_convergence_mm_day'][s],sst=ocean['mean_sst_c'][s],sst_oisst=ocean['mean_sst_c'][s],
            skin_sst=thermal['skin_temperature_c'][s],current=np.array([ocean['east_current_mps'][s],-ocean['south_current_mps'][s]]),
            evaporation=thermal['surface_latent_loss_wm2'][s],surface_radiation=thermal['surface_radiative_wm2'][s])
        fields.update({f'{name}_{s}':value for name,value in values.items()})
    return fields,cmp,captures

def analyze(archive,run_id):
    fields,cmp,captures=extract(archive,run_id)
    np.savez_compressed(archive/f'{run_id}-fields.npz',**fields)
    np.savez_compressed(archive/f'{run_id}-comparison.npz',**cmp)
    check=read(H/f'{run_id}-checks.json') if (H/f'{run_id}-checks.json').exists() else {'skipped':'User requested benchmark analysis without test suites.'}
    write(archive/f'{run_id}-analysis.json',dict(run=run_id,capture_paths=captures,checks=check,
        duration_note='Fixed eight coupling updates; convergence is diagnosed, not assumed.'))
    score(archive)

def score(archive):
    baseline=read(archive/'batch.json')['baseline']
    ids=sorted(read(p)['run'] for p in archive.glob('*-analysis.json'))
    source=Archive(archive,ids);inputs={i:source.load(i) for i in ids};metrics={str(i):{} for i in ids}
    for field,info in FIELDS.items():
        if not info[3]:continue
        for s in ([0] if field=='annual_rain' else range(4)):
            values={i:inputs[i][0][f'{field}_{s}'] for i in ids};shape=magnitude(values[baseline]).shape
            ref=source.reference(field,s,shape[-1]);mask=np.isfinite(magnitude(ref))
            for i in ids:mask &= mask_for(field,shape,*inputs[i][1:],s)
            for i in ids:
                assert np.isfinite(magnitude(values[i])[mask]).all(),(i,field,s)
                metrics[str(i)].setdefault(field,{})[str(s)]=evaluate(values[i],ref,mask,field)
    ranking=rank(metrics,baseline,ids)
    write(archive/'metrics.json',metrics);write(archive/'ranking.json',ranking)
    print(json.dumps(ranking,indent=2),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['setup','run','analyze','score','decision'])
    parser.add_argument('--archive',required=True,type=Path)
    parser.add_argument('--run',type=int)
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--source-archive',type=Path,help='Defaults to the selected baseline archive')
    parser.add_argument('--baseline',type=int,help='Defaults to the current selected baseline')
    parser.add_argument('--count',type=int,default=10)
    parser.add_argument('--skip-tests',action='store_true')
    args=parser.parse_args();archive=args.archive.absolute()
    if args.action=='setup':setup(archive,args.source_archive,args.baseline,args.count,args.skip_tests)
    elif args.action=='score':score(archive)
    elif args.action=='run':run(archive,args.run,read(args.plan))
    elif args.action=='analyze':analyze(archive,args.run)
    else:
        assert (archive/f'{args.run}-analysis.json').is_file()
        assert not (archive/str(args.run)/'decision.json').exists()
        write(archive/str(args.run)/'decision.json',read(args.plan))

if __name__=='__main__':main()
