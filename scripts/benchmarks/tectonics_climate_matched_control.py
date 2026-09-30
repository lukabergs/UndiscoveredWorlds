# /// script
# dependencies = ["numpy==2.5.2", "pillow==12.3.0"]
# ///
"""Compare spatial climate forcing with matching land-mean input fields."""
import argparse
import json
from pathlib import Path
import numpy as np
from tectonics_climate_feedback import erosion, spherical_weights, preview
from tectonics_experiment import digest, write_json

ROOT = Path(__file__).resolve().parents[2]


def run(run_id):
    archive = ROOT / f'runs/reports/tectonics/{run_id}'
    record = json.loads((ROOT / f'runs/manifests/tectonics/{run_id}.json').read_text())
    exe = archive / 'binaries/tectonic_climate_erosion.exe'
    if digest(exe) != record['binary_sha256'][exe.name]:
        raise ValueError('Archived erosion executable changed')
    initial = np.load(archive / 'climate-0/annual-fields.npz')['elevation']
    z = initial.copy()
    stages = []
    for cycle in range(record['cycles']):
        source = archive / f'climate-{cycle}/annual-fields.npz'
        receipt = source.parent / 'execution.json'
        if not receipt.exists():
            break  # Can be resumed as the main experiment finishes each climate stage.
        data = np.load(source)
        weights = spherical_weights(z.shape)
        land = data['elevation'] > 0
        means = {key: float(np.sum(data[key][land]*weights[land])/weights[land].sum())
                 for key in ['rain', 'temperature', 'wind', 'snow', 'ice']}
        folder = archive / f'matched-control-{cycle}'
        forcing_hash = digest(source)
        if folder.exists():
            stage = json.loads((folder / 'control.json').read_text())
            if stage['forcing_sha256'] != forcing_hash or stage['means'] != means:
                raise ValueError('Matched control forcing changed')
            if digest(folder / 'output/elevation.f64') != stage['elevation_sha256']:
                raise ValueError('Matched control terrain changed')
            z = np.fromfile(folder / 'output/elevation.f64', dtype='<f8').reshape(z.shape)
        else:
            fields = {key: np.full(z.shape, value) for key, value in means.items()}
            z, seconds = erosion(folder, z, fields, exe, record['recipe']['scenario'], record['surface_interval_myr'])
            stage = dict(cycle=cycle, means=means, forcing_sha256=forcing_hash,
                         elapsed_seconds=seconds, elevation_sha256=digest(folder / 'output/elevation.f64'))
            write_json(folder / 'control.json', stage)
        stages.append(stage)
    result = dict(run=run_id, completed_cycles=len(stages), stages=stages,
        status='complete' if len(stages)==record['cycles'] else 'partial',
        definition='Spatially uniform inputs match spherical land-area means of each corresponding climate field. Nonlinear erosion rates and evolving land masks are not forced to match.',
        source_sha256=digest(Path(__file__)), elapsed_seconds=sum(s['elapsed_seconds'] for s in stages))
    actual_path = archive / f'erosion-{record["cycles"]-1}/output/elevation.f64'
    if result['status']=='complete' and not (actual_path.parent.parent / 'execution.json').exists():
        result['status']='awaiting_surface'
    if result['status']=='complete':
        actual = np.fromfile(actual_path, dtype='<f8').reshape(z.shape)
        output = archive / 'matched-comparison-fields.npz'
        if not output.exists():
            np.savez_compressed(output, initial=initial, matched_control=z, final=actual, spatial_effect=actual-z)
        result['fields_sha256']=digest(output)
        result['spatial_effect_rms_m']=float(np.sqrt(np.mean((actual-z)**2)))
        palette=json.loads((ROOT/'configs/tectonics-palettes.json').read_text())['fields']['elevation']['stops']
        preview(archive/'matched-control.png',z,palette)
        preview(archive/'spatial-climate-effect.png',actual-z,[(-500,[35,70,150]),(0,[238,237,220]),(500,[150,48,40])])
    write_json(archive/'matched-control.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=int)
    run(parser.parse_args().run)
