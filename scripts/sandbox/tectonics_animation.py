"""Tectonics animations for the terrain sandbox: the plates (identities tracked between frames) and the
surface elevation, one frame per checkpoint bundle, with the evaluation's palettes and frame rendering
(scripts/evaluation/maps.py, the tectonics profile's map specs).

  uv run python scripts/sandbox/tectonics_animation.py [--run out/scratch/sandbox/tectonics] [--maps plates,elevation]

The run directory holds bundle/ (the run's end) and checkpoints/ (world_pipeline --checkpoint-every); the
APNGs go to <run>/maps/ beside an index.html.
"""
import argparse
import json
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'evaluation'))
from uwtools.bundle import Bundle  # noqa: E402
from maps import animations  # noqa: E402

PROFILE = ROOT / 'configs' / 'evaluation' / 'tectonics.json'


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run', type=Path, default=ROOT / 'out' / 'scratch' / 'sandbox' / 'tectonics')
    parser.add_argument('--maps', default='plates,elevation')
    arguments = parser.parse_args()

    wanted = arguments.maps.split(',')
    specs = [spec for spec in json.loads(PROFILE.read_text(encoding='utf-8'))['maps'] if spec['name'] in wanted]
    checkpoints = [Bundle.open(path) for path in sorted((arguments.run / 'checkpoints').iterdir())
                   if (path / 'manifest.json').exists()]
    # The run's end bundle; while the run is still going, the last checkpoint stands in for it.
    bundle = Bundle.open(arguments.run / 'bundle') if (arguments.run / 'bundle' / 'manifest.json').exists() else checkpoints[-1]
    written, _ = animations(bundle, specs, arguments.run, checkpoints)
    figures = ''.join(f'<figure><img src="{entry["path"]}"><figcaption>{entry["title"]} ({entry["frame_count"]} frames)'
                      f'</figcaption></figure>' for entry in written)
    (arguments.run / 'index.html').write_text(
        '<!doctype html><html><head><meta charset="utf-8"><title>Tectonics animation</title><style>'
        'body{background:#15181d;color:#ddd;font:14px system-ui;margin:16px}img{width:100%;max-width:1024px}'
        'figure{margin:0 0 24px}</style></head><body><h1>Tectonics, seed 1, 256 columns</h1>'
        f'{figures}</body></html>', encoding='utf-8')
    for entry in written:
        print(arguments.run / entry['path'], entry['frame_count'], 'frames')


if __name__ == '__main__':
    main()
