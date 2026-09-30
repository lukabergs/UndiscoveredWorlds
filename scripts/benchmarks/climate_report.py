# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Build the standard visual gallery, numerical ranking and durable decision report.

uv run --offline scripts/benchmarks/climate_report.py --archive runs/reports/<archive> --baseline 353 --runs 353 354 355 --apply-selection
No simulations run. Existing sealed experiment evidence is read-only.
"""
import argparse
import datetime
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

from climate_report_data import Archive,FIELDS,REGIONS,ROOT,SOURCES,read_json,magnitude,weights
from climate_report_metrics import evaluate,bins,rank,quantile
from climate_report_render import save_map,save_ocean_flow,save_wind_flow,save_continent_outlines,nearest_path,full_reference,flow_path,flow_complete
from climate_palettes import anchors,DEFINITION,VERSION
from run_storage import ensure_archive_directory

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def dump(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def relative(path,parent):return Path(os.path.relpath(path,parent)).as_posix()
def reference_input_hashes():
    names=sorted({name for group in SOURCES.values() for name in group})
    paths=[ROOT/f'refs/processed/climate/{name}_monthly.uwclim' for name in names]
    paths.append(ROOT/'refs/processed/climate/imerg_prec_annual.uwclim')
    return {relative(path,ROOT):sha(path) for path in paths}
def mask_for(field,shape,land,masks,s):
    mask=np.ones(shape,bool)
    if field in ('sst','sst_oisst','current','storage_depth'):
        mask&=masks[f'wet_{s}']
    if field in ('skin_sst','evaporation','surface_radiation'):mask&=land<.1
    if field=='850_wind':mask&=masks[f'850_{s}']
    return mask
def clean_text(text,aliases):
    for token,run in sorted(aliases.items(),key=lambda x:-len(x[0])):
        text=re.sub(r'\b'+re.escape(token)+r'\b',str(run),text)
    return re.sub(r'\b(\d+)\s*/\s*\1\b',r'\1',text)
def profiles(a,ref,mask):
    a=magnitude(a);ref=magnitude(ref);lat=90-(np.arange(a.shape[0])+.5)*180/a.shape[0]
    lon=-180+(np.arange(a.shape[1])+.5)*360/a.shape[1]
    result={}
    for sector,(lo,hi) in {'Global':(-180,180),'Central / eastern Pacific':(-160,-90),'Western Pacific':(125,175)}.items():
        valid=mask&(lon[None,:]>=lo)&(lon[None,:]<=hi)
        rows=valid.any(axis=1)
        result[sector]=dict(latitude=lat[rows].tolist(),model=[float(a[y][valid[y]].mean()) for y in np.flatnonzero(rows)],reference=[float(ref[y][valid[y]].mean()) for y in np.flatnonzero(rows)])
    return result

def build(args):
    ids=sorted(set([args.baseline,*args.runs]));archive=Archive(ROOT/args.archive,ids)
    name=f'climate-comparison-{min(ids)}-{max(ids)}';parent=ROOT/'runs/reports'
    out=ensure_archive_directory(parent/name);out.mkdir(exist_ok=True)
    refroot=ROOT/'refs/processed/climate/maps/standard-v2';refroot.mkdir(parents=True,exist_ok=True)
    title=f'({args.baseline}) - {min(ids)}-{max(ids)}'
    inputs={run:archive.load(run) for run in ids};metrics={str(run):{} for run in ids}
    records={};refrecords={};nearest_records={};nearest_refs={};reference_audit={};fieldinfo={}
    wind_fields=('surface_wind','850_wind','500_wind')
    pool=ThreadPoolExecutor(max_workers=args.render_workers);flow_jobs=[]
    def flow_ready(prefix):
        return args.reuse_maps and flow_complete(prefix)
    for field,(label,unit,kind,weight,semantic) in FIELDS.items():
        print('Field',field,flush=True)
        lo,hi=DEFINITION['scales'][field];season_ids=[0] if field=='annual_rain' else list(range(4))
        fieldinfo[field]=dict(label=label,unit=unit,kind=kind,weight=weight,range=[lo,hi],
            palette=anchors(lo,hi,lo<0),bins=[None if not np.isfinite(x) else float(x) for x in bins(field)])
        full_cache={}
        for s in season_ids:
            models={run:inputs[run][0][f'{field}_{s}'] for run in ids}
            shape=magnitude(models[args.baseline]).shape
            ref=archive.reference(field,s,shape[-1]);mask=mask_for(field,shape,*inputs[args.baseline][1:],s)
            for run in ids:
                assert magnitude(models[run]).shape==shape,(run,field)
                # Availability differences must be explicit and common, never candidate-selected.
                mask&=mask_for(field,shape,*inputs[run][1:],s)
            if ref is not None:mask&=np.isfinite(magnitude(ref))
            for run in ids:
                if not np.isfinite(magnitude(models[run])[mask]).all():raise ValueError(f'Non-finite active field: {run}/{field}/{s}')
                metrics[str(run)].setdefault(field,{})
                if ref is not None:
                    item=evaluate(models[run],ref,mask,field);item['profiles']=profiles(models[run],ref,mask)
                    metrics[str(run)][field][str(s)]=item
                image=out/'maps'/semantic.format(s=s+1)/f'{run}.png'
                if not args.reuse_maps or args.rebuild_scalars or not image.exists():save_map(image,field,models[run],mask,args.width)
                records[f'{run}/{field}/{s}']=relative(image,parent)
                pixel_image=nearest_path(image)
                if not args.reuse_maps or args.rebuild_scalars or not pixel_image.exists():save_map(pixel_image,field,models[run],mask,args.width,False)
                nearest_records[f'{run}/{field}/{s}']=relative(pixel_image,parent)
                if field in wind_fields:
                    prefix=image.with_suffix('')
                    if not flow_ready(prefix):flow_jobs.append(pool.submit(save_wind_flow,prefix,field,models[run],mask,inputs[run][1],args.width))
                    for style in ('lic','particles'):records[f'{run}/{field}_{style}/{s}']=relative(flow_path(prefix,style),parent)
                if field=='current':
                    prefix=image.with_suffix('')
                    if not flow_ready(prefix):
                        flow_jobs.append(pool.submit(save_ocean_flow,prefix,models[run],mask,args.width))
                    for style in ('lic','particles'):records[f'{run}/current_{style}/{s}']=relative(flow_path(prefix,style),parent)
            if ref is None:
                reference_audit[field]='Unavailable: no retained ocean/ice skin-temperature reference. SST is not substituted.'
                continue
            reference_audit[field]=('Context only: GLORYS density-defined mixed layer is not model thermal storage depth.' if field=='storage_depth' else
                'GLORYS near-surface current versus model effective layer velocity; vertical definitions differ.' if field=='current' else
                'Available; source and period recorded in reference manifest.')
            full=full_reference(field,s,archive,full_cache)
            refimage=refroot/'climate'/semantic.format(s=s+1)/f'{args.width}.png'
            if not args.reuse_maps or args.rebuild_scalars or not refimage.exists():save_map(refimage,field,full,np.isfinite(magnitude(full)),args.width)
            refrecords[f'{field}/{s}']=relative(refimage,parent)
            pixel_ref=nearest_path(refimage)
            if not args.reuse_maps or args.rebuild_scalars or not pixel_ref.exists():save_map(pixel_ref,field,full,np.isfinite(magnitude(full)),args.width,False)
            nearest_refs[f'{field}/{s}']=relative(pixel_ref,parent)
            if field in wind_fields:
                prefix=refimage.with_suffix('')
                if not flow_ready(prefix):flow_jobs.append(pool.submit(save_wind_flow,prefix,field,full,np.isfinite(magnitude(full)),inputs[args.baseline][1],args.width))
                for style in ('lic','particles'):refrecords[f'{field}_{style}/{s}']=relative(flow_path(prefix,style),parent)
            if field=='current':
                prefix=refimage.with_suffix('')
                if not flow_ready(prefix):flow_jobs.append(pool.submit(save_ocean_flow,prefix,full,np.isfinite(magnitude(full)),args.width))
                for style in ('lic','particles'):refrecords[f'current_{style}/{s}']=relative(flow_path(prefix,style),parent)
        del full_cache
    for i,job in enumerate(as_completed(flow_jobs),1):
        job.result();print(f'Flow products {i}/{len(flow_jobs)} complete',flush=True)
    pool.shutdown()
    outline=out/'overlays/continents.png';save_continent_outlines(outline,args.width)
    wind_products={}
    for field in wind_fields:
        for s in range(4):
            prefix=(parent/refrecords[f'{field}/{s}']).resolve().with_suffix('')
            for style in ('lic','particles'):
                path=flow_path(prefix,style);texture=flow_path(prefix,style,'.tif')
                wind_products[f'{field}_{style}/{s}']=dict(path=relative(path,ROOT),sha256=sha(path),numeric_texture=relative(texture,ROOT),numeric_texture_sha256=sha(texture))
    dump(refroot/'wind-flow-manifest.json',dict(palette_version=VERSION,products=wind_products,seed=20260906,particle_days=4.5,method='Same spherical LIC and 4.5-day particle kernels and seed for model and reference. Wind layer limits 25/30/40 m/s.'))
    ranking=rank(metrics,args.baseline,ids)
    payload=dict(schema=1,title=title,baseline=args.baseline,runs=ids,fields=fieldinfo,regions=REGIONS,
        maps=records,references=refrecords,nearest_maps=nearest_records,nearest_references=nearest_refs,continent_outlines=relative(outline,parent),metrics=metrics,ranking=ranking,reference_audit=reference_audit,
        palette_version=VERSION,report=name+'.md')
    dump(out/'data.json',payload)
    dump(out/'ranking.json',ranking)
    dump(out/'reference-coverage.json',reference_audit)
    provenance=dict(palette_version=VERSION,palette_sha256=sha(ROOT/'configs/climate-palettes.json'),
        native_sources=archive.provenance(),input_sha256=reference_input_hashes(),reference_maps=refrecords,
        method='Wind and pressure: January/April/July/October monthly means. Other seasonal fields: day-weighted quarters beginning in those months. Native annual rainfall: IMERG. Full-reference images use retained bundles at 2048; numerical tables conservatively match each model grid. Ocean current LIC/particles use identical 45-day spherical kernels, seed 20260906, 0–1 m/s.',
        unavailable=reference_audit,images_sha256={relative(parent/p,ROOT):sha(parent/p) for p in refrecords.values()},
        ocean_texture_meaning='LIC TIFF: dimensionless texture intensity. Particle TIFF: dimensionless trail intensity. Neither is speed or concentration.',
        ocean_texture_sha256={relative(flow_path(refroot/'climate/ocean/current'/str(s)/'s'/str(args.width),style,'.tif'),ROOT):sha(flow_path(refroot/'climate/ocean/current'/str(s)/'s'/str(args.width),style,'.tif')) for s in range(1,5) for style in ('lic','particles')})
    dump(refroot/'manifest.json',provenance)
    write_reference_index(refroot,refrecords,reference_audit,parent)
    template=(ROOT/'scripts/benchmarks/climate_report_template.html').read_text(encoding='utf-8')
    encoded=json.dumps(payload,separators=(',',':'),allow_nan=False).replace('<','\\u003c')
    document=template.replace('@@TITLE@@',html.escape(title)).replace('@@DATA@@',encoded)
    (parent/f'{name}.html').write_text(document,encoding='utf-8')
    for run in ids:
        runpage=parent/f'{run}.html'
        if runpage.exists() and 'standard-report-pending' in runpage.read_text(encoding='utf-8'):
            target=f'{name}.html?run={run}'
            runpage.write_text(f'<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url={target}"><a href="{target}">{run}</a>',encoding='utf-8')
    write_report(parent/name,archive,ids,ranking,reference_audit)
    if args.apply_selection:apply_selection(out,archive,ranking,args.override_baseline,args.override_reason)
    (parent/'climate-latest-comparison.html').write_text(f'<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url={name}.html"><a href="{name}.html">{title}</a>',encoding='utf-8')
    dump(out/'generation.json',dict(command=sys.argv,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        simulations_executed=0,source_archive=str(archive.path),run_ids=ids,baseline=args.baseline,
        input_sha256={str(run):sha(archive.path/f'{archive.entries[run][0]}-fields.npz') for run in ids},
        scripts_sha256={p.name:sha(p) for p in (ROOT/'scripts/benchmarks').glob('climate_report*')},
        texture_kernel_sha256=sha(ROOT/'src/validation/climate/flow_texture_kernel.cpp'),
        texture_wrapper_sha256=sha(ROOT/'scripts/benchmarks/climate_flow_texture.py'),
        data_sha256=sha(out/'data.json'),html_sha256=sha(parent/f'{name}.html')))
    print(json.dumps(ranking,indent=2),flush=True)

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

def write_report(prefix,archive,ids,ranking,audit):
    baseline=ranking['baseline'];best=ranking['best_in_batch'];selected=ranking['recommended_baseline']
    aliases={key:run for run,(key,_) in archive.entries.items()}
    for run,(key,_) in archive.entries.items():
        plan=archive.path/f'{key}-plan.json'
        if plan.exists():
            info=read_json(plan)
            if info.get('case'):aliases[info['case']]=run
            match=re.fullmatch(r'([A-Za-z]+\d+)/(\d+)',info.get('parent',''))
            if match:aliases[match[1]]=int(match[2])
    lines=[f'# ({baseline}) - {min(ids)}-{max(ids)}','',f'[Visual comparison]({prefix.name}.html). No new simulations were run for this report.','',
        f'Numerical best in this batch: **{best}**. Recommended baseline: **{selected}**. Its aggregate cost changes by {ranking["improvement_percent"]:.3f}% relative to {baseline} (positive is improvement). This is a numerical recommendation that the user can override after visual inspection.','',
        '## Why this run','', '| Field | Best-run loss change vs baseline | Weight |','| --- | ---: | ---: |']
    for f,weight in ranking['field_weights'].items():
        change=100*(ranking['field_losses'][best][f]/ranking['field_losses'][baseline][f]-1)
        lines.append(f'| {FIELDS[f][0]} | {change:+.2f}% | {weight:.0%} |')
    changes={f:ranking['field_weights'][f]*(ranking['field_losses'][best][f]-ranking['field_losses'][baseline][f]) for f in ranking['field_weights']}
    good=sorted(changes,key=changes.get)[:3];bad=sorted(changes,key=changes.get,reverse=True)[:3]
    lines+=['', 'The largest contributions to the score improvement are '+(', '.join(FIELDS[f][0] for f in good if changes[f]<0) or 'none')+'. The main tradeoffs are '+(', '.join(FIELDS[f][0] for f in bad if changes[f]>0) or 'none')+'. The tables retain individual error and tail measurements. Global diagnostics can still miss regional pattern errors; the regional maps are retained for user visual review.','',
        '## Fixed selection rule','',ranking['formula'],'',
        'The score uses fixed physical display bins, unclipped numerical values including underflow/overflow bins, spherical area weights and common valid reference/terrain masks. Spatial errors keep a correct histogram with misplaced features from winning on distribution alone. The baseline changes only when the batch winner has a strictly lower cost; otherwise it is retained. Missing required fields are an error, not silently dropped. A visual override takes precedence until a later batch explicitly applies a new recommendation.','',
        'Mixed-layer depth has a different physical definition and skin temperature lacks an exact retained reference; neither enters the score. ERA5 SST is a secondary diagnostic; OISST supplies the SST score. All four seasons are equally weighted. Annual IMERG and seasonal ERA5 rainfall are distinct comparisons. The rank is within this fixed-resolution, bounded-duration Earth experiment, not evidence of convergence or transfer to procedural worlds.','',
        '## Report, log and evidence','',
        'This report is the readable development record: intended changes, findings and decisions. A log is raw program output: timing, solver messages, warnings and failures. JSON receipts hold machine-readable metrics, inputs and hashes. The HTML is the visual workspace and contains only selectors, maps, charts and compact measurement tables. These roles avoid repeating the full narrative in each place.','',
        f'[Ranking JSON]({prefix.name}/ranking.json) · [Reference coverage]({prefix.name}/reference-coverage.json) · [Generation receipt]({prefix.name}/generation.json).','',
        'The original sealed batch and its historical renderings remain unchanged. The new standard images use climate-colour-v2; native exporters use the same palette definition for subsequent runs. The 2048-wide images do not increase the native dynamical resolution.','', '## Run history','']
    for run in ids:
        key,analysis=archive.entries[run];folder=archive.path/key
        plan=archive.path/f'{key}-plan.json';decision=folder/'decision.json'
        lines +=[f'### {run}','']
        if plan.exists():lines +=[clean_text(read_json(plan)['hypothesis'],aliases),'']
        if decision.exists():lines +=[clean_text(read_json(decision)['insight'],aliases),'']
        links=[f'[Raw log](../logs/climate/{run}.log)',f'[Numerical analysis]({relative(archive.path/f"{key}-analysis.json",prefix.parent)})']
        if (folder/'changes-from-parent.patch').exists():links.append(f'[Source patch]({relative(folder/"changes-from-parent.patch",prefix.parent)})')
        lines+=[' · '.join(links),'']
    findings=archive.path/'batch-findings.json'
    if findings.exists():
        review=read_json(findings)
        lines+=['## Cross-run findings and experiment limits','']
        lines.extend(text+'\n' for text in review.get('findings',[]))
        for label,path in review.get('evidence',{}).items():
            lines.append(f'[{label}]({relative(archive.path/path,prefix.parent)})')
    prefix.with_suffix('.md').write_text('\n'.join(lines),encoding='utf-8')

def apply_selection(out,archive,ranking,override=None,reason=None):
    p=ROOT/'runs/reports/climate-global-current-selection.json';previous=read_json(p)
    before=out/'selection-before.json'
    if not before.exists():dump(before,previous)
    if previous.get('selection_mode')=='visual_override' and previous.get('batch')==out.name and override is None:
        dump(out/'selection-after.json',previous)
        return
    if override is not None and (override not in archive.ids or not reason):
        raise ValueError('Visual override requires a displayed run ID and a reason')
    run=override if override is not None else ranking['recommended_baseline'];key,analysis=archive.entries[run]
    folder=archive.path/key
    if not (folder/'manifest.json').exists():
        folder=(ROOT/analysis['capture_paths'][0]).parent.parent
    manifest=read_json(folder/'manifest.json')
    retained=archive.path/'retained-candidate.json'
    working=read_json(retained).get('run') if retained.exists() else None
    selection=dict(selected_run=run,reference_run_for_future_experiments=run,status='Numerically selected baseline; user visual override takes precedence',
        authority='User authorized automatic best-of-batch selection on 2026-09-18.',
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),numerical_recommendation=ranking['best_in_batch'],
        prior_baseline=ranking['baseline'],ranking_policy=ranking['version'],ranking=relative(out/'ranking.json',p.parent),
        comparison_policy='Use the selected run ID and observed references. No permanently protected historical wind run.',
        source_manifest=relative(folder/'manifest.json',p.parent),source_zip_sha256=manifest['source_zip_sha256'],
        executable_sha256=manifest['executable_sha256'],live_source_and_executable_match=False,
        working_physics_run=working,working_note='Working physics is recorded separately from the selected baseline. Use the selected archived source for the next batch unless a visual override changes the selection.',
        climate_reference_forcing=False,new_simulations_run_for_this_selection=0)
    selection.update(batch=out.name,selection_mode='visual_override' if override is not None else 'numerical',visual_override_reason=reason)
    dump(p,selection);dump(out/'selection-after.json',selection)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',required=True)
    p.add_argument('--baseline',type=int,required=True);p.add_argument('--runs',type=int,nargs='+',required=True)
    p.add_argument('--width',type=int,default=2048);p.add_argument('--reuse-maps',action='store_true')
    p.add_argument('--render-workers',type=int,default=4)
    p.add_argument('--rebuild-scalars',action='store_true',help='Re-export magnitude/scalar maps while reusing complete ocean textures')
    p.add_argument('--apply-selection',action='store_true')
    p.add_argument('--override-baseline',type=int)
    p.add_argument('--override-reason');args=p.parse_args()
    # Arguments intentionally describe a user decision; ordinary rebuilds preserve an override.
    if args.width!=2048:raise SystemExit('Standard report output is fixed to 2048x1024; native fields keep their own grids.')
    if args.override_baseline is not None and not args.apply_selection:raise SystemExit('An override requires --apply-selection')
    build(args)
if __name__=='__main__':main()
