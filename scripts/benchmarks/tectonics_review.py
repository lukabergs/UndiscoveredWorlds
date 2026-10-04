# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Review archived tectonic runs before selecting the next experiment."""
import argparse
import base64
import html
import json
from pathlib import Path
import sys
import numpy as np
from PIL import Image, ImageDraw
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import palettes
from uwtools.hashing import sha256
from tectonics_report_data import load_bundle, diagnostics
from tectonics_earth_reference import load_earth_reference, compare_to_earth
from tectonics_experiment import write_json
from map_rendering import colourize

ROOT = Path(__file__).resolve().parents[2]


def write_overview(results, output):
    """Small visual index; detailed comparisons and APNGs remain in the common report."""
    output = Path(output)
    if output.parent.resolve() != (ROOT / 'runs/reports/tectonics').resolve():
        raise ValueError('Overview must use the existing tectonics report directory')
    rows, cards = [], []
    for result in results:
        run, m = result['run'], result['final']
        note = ROOT / f'runs/reports/tectonics/{run}.md'
        links = [f'<a href="{run}.md">Findings</a>'] if note.exists() else []
        for other in [run-1] + [r['run'] for r in results if r['run'] < run]:
            report = ROOT / f'runs/reports/tectonics/{other}-{run}.html'
            link = f'<a href="{report.name}">{other} → {run}: maps + animations</a>'
            if other != run and report.is_file() and link not in links:
                links.append(link)
        d, islands = m['terrain_detail'], m['islands']
        rows.append(f'<tr><td><a href="#run{run}">{run}</a></td>'
            f'<td>{d["interior_ocean_relief_rms_m"]:.1f}</td>'
            f'<td>{d["interior_land_relief_rms_m"]:.1f}</td>'
            f'<td>{islands["islands_2_to_64_cells"]}</td><td>{islands["single_cell_islands"]}</td>'
            f'<td>{m["land_fraction"]*100:.2f}%</td><td>{m["ocean_depth_m"]["mean"]:.0f}</td>'
            f'<td>{m["earth_comparison"]["hypsography_mean_absolute_quantile_gap_m"]:.0f}</td></tr>')
        png = ROOT / f'runs/diagnostics/tectonics/{run}/iteration-preview.png'
        image = base64.b64encode(png.read_bytes()).decode('ascii')
        cards.append(f'<article id="run{run}"><h2>Run {run}</h2><p>{html.escape(result["recipe"])}</p>'
            f'<p>{" · ".join(links)}</p><img alt="Run {run}: final elevation and terrain detail" '
            f'src="data:image/png;base64,{image}"></article>')
    text = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Tectonics iteration review</title><style>
body{margin:24px;background:#101b28;color:#e3eaf0;font:16px system-ui;line-height:1.5}main{max-width:1500px;margin:auto}
a{color:#84d8ec}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{text-align:right;padding:9px;border-bottom:1px solid #415162}th:first-child,td:first-child{text-align:left}
.scroll{overflow:auto}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(400px,1fr));gap:20px;margin-top:24px}article{background:#192737;padding:15px;border-radius:10px}article img{width:100%;height:auto;image-rendering:pixelated}article p{overflow-wrap:anywhere}h1,h2{line-height:1.2}@media(max-width:480px){.cards{grid-template-columns:1fr}body{margin:12px}}
</style><main><h1>Tectonics iteration review</h1>
<p>Compare each experiment with the preceding run and consult its findings for the next decision.
All maps use the same elevation and detail palettes. The terrain-detail image is a diagnostic, not an added noise layer.</p>
'''
    note_name = output.stem.removesuffix('-overview') + '.md'
    text += f'<p><a href="{html.escape(note_name)}">Decisions and limitations</a></p>'
    text += '''
<div class="scroll"><table><thead><tr><th>Run</th><th>Ocean detail RMS (m)</th><th>Land detail RMS (m)</th><th>Islands, 2–64 cells</th><th>Single cells</th><th>Land area</th><th>Mean depth (m)</th><th>Earth curve gap (m)</th></tr></thead><tbody>'''
    text += ''.join(rows) + '</tbody></table></div><p>Detail samples require a full 5×5 window on the same surface and plate. Higher RMS can include artifacts. Islands use four-connected land with both periodic seams; one native cell is 156.25 km wide. Earth comparison is distribution context, not a geographical reconstruction.</p><section class="cards">'
    text += ''.join(cards) + '</section></main></html>'
    output.write_text(text, encoding='utf-8')
    print(f'Wrote {output}')


def review(run_id, earth):
    record = json.loads((ROOT / f'runs/manifests/tectonics/{run_id}.json').read_text())
    if record['status'] != 'complete':
        raise ValueError('Review requires a completed run; preserve pending replay evidence')
    states, previous = [], None
    for snap in record['snapshots']:
        manifest, fields, _ = load_bundle(ROOT / f'runs/checkpoints/tectonics/{run_id}' / snap['directory'],
                                          record['field_sha256'][snap['label']])
        m = diagnostics(manifest, fields)
        temporal = None
        if previous is not None and manifest['time_myr'] > previous[0]:
            delta = fields['elevation'] - previous[1]
            temporal = {'interval_myr': manifest['time_myr']-previous[0],
                        'elevation_change_rms_m': float(np.sqrt(np.mean(delta**2))),
                        'elevation_change_p99_m': float(np.quantile(abs(delta), .99))}
        previous = (manifest['time_myr'], fields['elevation'].copy())
        states.append({'label': snap['label'], 'time_myr': manifest['time_myr'],
                       'invariants_pass': m['invariants_pass'], 'detail': m['terrain_detail'],
                       'islands': m['islands'], 'temporal': temporal,
                       'land_fraction': m['land_fraction'], 'peak_m': m['land_elevation_m']['max'],
                       'deepest_m': m['ocean_depth_m']['max'],
                       'coastal_p95_m': m['coastal_terrain']['coastal_land_m']['p95'],
                       'coastal_jump_p95_m': m['coastal_terrain']['coastal_face_jump_m']['p95'],
                       'coastal_highland_fraction': m['coastal_terrain']['highland_fraction_on_coast'],
                       'motion': m['motion_budget'], 'events': m['plate_events'],
                       'mass_residual': m.get('column_budget', {}).get('relative_mass_residual')})
    m['earth_comparison'] = compare_to_earth(m, earth)
    result = {'run': run_id, 'seed': record['recipe']['seed'], 'recipe': record['recipe']['name'],
              'exact_replay': record.get('replay', {}).get('exact_match'),
              'all_checkpoints_pass': all(s['invariants_pass'] for s in states),
              'states': states, 'final': m,
              'review_source_sha256': {p.name: sha256(p) for p in [Path(__file__),
                  ROOT / 'scripts/benchmarks/tectonics_report_data.py', palettes.TECTONICS_PATH]}}
    output = ROOT / f'runs/diagnostics/tectonics/{run_id}'
    write_json(output / 'iteration-review.json', result)
    field_palettes = palettes.tectonics()['fields']
    panels = []
    for field in ['elevation', 'local_relief']:
        pixels = colourize(fields[field], field_palettes[field]['stops'])
        im = Image.fromarray(pixels).resize((1024,512), Image.Resampling.NEAREST)
        panels.append(im)
    preview = Image.new('RGB',(1024,1084),'#101b28')
    draw = ImageDraw.Draw(preview)
    draw.text((12,8),f'Run {run_id} final elevation; fixed shared palette',fill='white')
    preview.paste(panels[0],(0,30))
    draw.text((12,550),'Terrain detail: departure from a 5 x 5 cell mean (not added relief)',fill='white')
    preview.paste(panels[1],(0,572))
    preview.save(output / 'iteration-preview.png')
    summary = {'run': run_id, 'replay': result['exact_replay'], 'all_states_pass': result['all_checkpoints_pass'],
        'initial_ocean_detail': states[0]['detail']['interior_ocean_relief_rms_m'],
        'final_detail': m['terrain_detail'], 'initial_islands': states[0]['islands'], 'final_islands': m['islands'],
        'land':m['land_fraction'], 'peak':m['land_elevation_m']['max'], 'mean_depth':m['ocean_depth_m']['mean'],
        'coastal_highland_fraction':m['coastal_terrain']['highland_fraction_on_coast'],
        'earth_gap':m['earth_comparison']['hypsography_mean_absolute_quantile_gap_m'],
        'mass_residual':m.get('column_budget',{}).get('relative_mass_residual')}
    print(json.dumps(summary,indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',type=int,nargs='+')
    parser.add_argument('--overview',type=Path)
    args = parser.parse_args()
    earth = load_earth_reference(ROOT)
    results = [review(run,earth) for run in args.runs]
    if args.overview:
        write_overview(results,args.overview)
