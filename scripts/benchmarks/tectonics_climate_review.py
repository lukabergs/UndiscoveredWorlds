# /// script
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Summarize archived climate-feedback fields and budgets without creating HTML."""
import argparse
import json
from pathlib import Path
import numpy as np
from tectonics_experiment import write_json, digest
from tectonics_earth_reference import load_earth_reference, weighted_distribution, QUANTILES
from tectonics_climate_feedback import spherical_weights
from tectonics_batch_analysis import distribution
from tectonics_report_data import load_bundle, local_relief, terrain_detail, coastal_terrain

ROOT = Path(__file__).resolve().parents[2]


def review(run):
    archive = ROOT / f'runs/reports/tectonics/{run}'
    manifest = json.loads((ROOT / f'runs/manifests/tectonics/{run}.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Feedback experiment is incomplete')
    if digest(archive / 'comparison-fields.npz') != manifest['fields_sha256']:
        raise ValueError('Feedback field hash mismatch')
    data = dict(np.load(archive / 'comparison-fields.npz'))
    matched = None
    if (archive / 'matched-control.json').exists():
        matched = json.loads((archive / 'matched-control.json').read_text())
        if matched['status'] == 'complete':
            if digest(archive / 'matched-comparison-fields.npz') != matched['fields_sha256']:
                raise ValueError('Matched comparison hash mismatch')
            data['matched_control'] = np.load(archive / 'matched-comparison-fields.npz')['matched_control']
    parent = json.loads((ROOT / f'runs/manifests/tectonics/{manifest["parent"]}.json').read_text())
    snapshot = next(s for s in parent['snapshots'] if s['label'] == 'final')
    native, parent_fields, _ = load_bundle(ROOT / f'runs/checkpoints/tectonics/{manifest["parent"]}' / snapshot['directory'])
    reference = load_earth_reference(ROOT)
    earth = np.asarray(reference['metrics']['hypsography']['elevation_m'])
    results = {}
    for name in ['initial', 'final', 'uniform_control', *(['matched_control'] if 'matched_control' in data else [])]:
        z = data[name]
        if not np.isfinite(z).all():
            raise ValueError('Nonfinite terrain')
        fields = dict(parent_fields, elevation=z, local_relief=local_relief(z))
        fields.pop('coast_distance_km', None)  # Recompute after shoreline changes.
        results[name] = dict(native_detail=terrain_detail(fields), native_coast=coastal_terrain(native, fields))
        for geometry, weights in [('native_equal_area', np.ones_like(z)), ('spherical', spherical_weights(z.shape))]:
            stats, curve = weighted_distribution(z, weights)
            # Native reports use interpolated unweighted quantiles.
            if geometry == 'native_equal_area':
                curve['elevation_m'] = np.quantile(z, QUANTILES).tolist()
            results[name][geometry] = dict(elevation_m=stats, hypsography=curve,
                fraction_below=distribution(z, weights),
                land_fraction=float(weights[z > 0].sum() / weights.sum()),
                earth_gap_m=float(np.trapezoid(abs(np.asarray(curve['elevation_m']) - earth), QUANTILES)))
    climates = []
    for folder in sorted(archive.glob('climate-*')):
        if not folder.is_dir():
            continue
        execution = json.loads((folder / 'execution.json').read_text())
        physics = json.loads((folder / 'physics-diagnostics.json').read_text())
        climates.append(dict(stage=folder.name, **execution, physics=physics))
    erosion = []
    for pattern in ['erosion-*', 'uniform-control-*', 'matched-control-*']:
        for folder in sorted(archive.glob(pattern)):
            if not folder.is_dir():
                continue
            budget = json.loads((folder / 'output/budget.json').read_text())
            if abs(budget['relative_volume_residual']) > 1e-9:
                raise ValueError('Erosion volume budget failed')
            execution = json.loads((folder / 'execution.json').read_text())
            erosion.append(dict(stage=folder.name, **budget, **execution))
    change = data['final'] - data['initial']
    effect = data['final'] - data['uniform_control']
    summary = dict(run=run, parent=manifest['parent'], elapsed_seconds=manifest['elapsed_seconds'],
        review_source_sha256=digest(Path(__file__)), earth_reference_cache=reference['cache'],
        distributions=results, climate_stages=climates, erosion_stages=erosion,
        change_rms_m=float(np.sqrt(np.mean(change**2))), climate_effect_rms_m=float(np.sqrt(np.mean(effect**2))),
        maximum_lowering_m=float(-change.min()), maximum_deposition_m=float(change.max()),
        matched_control=matched, limitations=manifest['limitations'], geometry=manifest['geometry'])
    write_json(archive / 'review.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ['distributions', 'climate_stages', 'erosion_stages']}))
    for name, metrics in results.items():
        m = metrics['native_equal_area']
        print(name, json.dumps(dict(land_percent=100*m['land_fraction'], peak_m=m['elevation_m']['max'],
            earth_gap_m=m['earth_gap_m'], above_2000_percent=100*(1-m['fraction_below']['2000']))))
    print('climate convergence:', [c['coupling_converged'] for c in climates])
    print('maximum volume residual:', max(abs(e['relative_volume_residual']) for e in erosion))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=int)
    review(parser.parse_args().run)
