"""Compact, cached distribution and runtime evidence for tectonic parameter batches."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from PIL import Image
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import grid, palettes, runs
from uwtools.hashing import sha256
from uwtools.tectonics import snapshot_manifest
from tectonics_report_data import load_bundle, write_json
from tectonics_review import review
from uwtools.tectonics_earth_reference import load_earth_reference
from map_rendering import colourize

ROOT = Path(__file__).resolve().parents[2]
LEVELS = [-6000, -4000, -2000, -1000, -258, -200, -160, 0, 1000, 2000, 4000]


def distribution(z, weights=None):
    if weights is None:
        weights = np.ones_like(z)
    total = weights.sum()
    return {str(level): float(weights[z < level].sum()/total) for level in LEVELS}


def earth_thresholds():
    cache = ROOT/'refs/processed/metadata/2048/tectonics-shelf-thresholds.json'
    source = ROOT/'refs/processed/heightmap/2048.tif'
    source_hash = sha256(source)
    if cache.exists():
        data = json.loads(cache.read_text())
        if data['input_sha256'] != source_hash:
            raise ValueError('Earth threshold source changed')
        return data
    z = np.array(Image.open(source))
    weights = grid.area_weights(*z.shape)
    data = dict(schema='tectonics-shelf-thresholds/v1', input_sha256=source_hash,
        weighting='spherical cell area', fraction_below=distribution(z,weights),
        shelf_proxy='Signed elevation -200 to 0 m; bathymetric proxy, not a geomorphological shelf mask',
        shelf_break_study='https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL071130',
        shelf_break_study_median_m=160, shelf_break_study_mean_m=258)
    write_json(cache,data)
    return data


def analyze(run_id):
    record = runs.tectonics(run_id)
    output = runs.output_directory('tectonics', run_id)
    cached = output/'iteration-review.json'
    result = json.loads(cached.read_text()) if cached.exists() else review(run_id,load_earth_reference(ROOT))
    m = result['final']
    thresholds = {}
    final_bundle = None
    for label in ['initial','final']:
        snap = next(s for s in record['snapshots'] if s['label']==label)
        bundle, fields, _ = load_bundle(snap['path'], record.get('field_sha256', {}).get(label))
        thresholds[label] = distribution(fields['elevation'])
        if label == 'final': final_bundle = bundle
    fracture_history = []
    for snap in record['snapshots']:
        state = snapshot_manifest(snap['path'])
        diagnostics = state.get('topology_diagnostics', {})
        fracture_history.append(dict(time_myr=state['time_myr'],
            **{key: diagnostics.get(key, 0) for key in
               ['max_rift_damage', 'max_tension_n_m', 'damaged_rift_cells']}))
    peak_index = max(range(len(fracture_history)), key=lambda i: fracture_history[i]['damaged_rift_cells'])
    if fracture_history[peak_index]['damaged_rift_cells'] > 0:
        snap = record['snapshots'][peak_index]
        _, damaged, _ = load_bundle(snap['path'])
        palette = palettes.tectonics()
        for name, field in [('rift-damage-peak', 'rift_damage'), ('rift-terrain-peak', 'elevation')]:
            stops = palette['fields'][field].get('stops', palette['sequential'])
            Image.fromarray(colourize(damaged[field], stops)).resize(
                (1024, 512), Image.Resampling.NEAREST).save(output/f'{name}.png')
    c = m['coastal_terrain']
    summary = dict(run=run_id, seed=record['recipe']['seed'], elapsed_seconds=record.get('elapsed_seconds'),
        land_percent=100*m['land_fraction'], peak_m=m['elevation_m']['max'],
        mean_depth_m=m['ocean_depth_m']['mean'], coast_p95_m=c['coastal_land_m']['p95'],
        coast_jump_p95_m=c['coastal_face_jump_m']['p95'],
        earth_gap_m=m['earth_comparison']['hypsography_mean_absolute_quantile_gap_m'],
        ocean_detail_m=m['terrain_detail']['interior_ocean_relief_rms_m'],
        land_detail_m=m['terrain_detail']['interior_land_relief_rms_m'],
        events=m['plate_events'], checks_pass=result['all_checkpoints_pass'],
        fraction_below=thresholds, earth_fraction_below=earth_thresholds()['fraction_below'],
        final_column_budget=final_bundle.get('column_budget', {}),
        fracture_history=fracture_history)
    write_json(output/'batch-summary.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in
        ['fraction_below','earth_fraction_below','events','final_column_budget','fracture_history']}))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',type=int,nargs='+')
    args=parser.parse_args()
    for run_id in args.runs: analyze(run_id)
