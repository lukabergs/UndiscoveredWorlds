"""Command-line checks of world_pipeline's world grid (--grid, cube by default since T8 W3.2); seed 1.

    uv run --offline python tests/scripts/check_grid_cli.py out/build/headless/bin/Release/world_pipeline.exe

Here: a run from a bundle without --grid keeps the bundle's grid (a lat-lon bundle stays lat-lon,
a cube bundle stays cube), a --grid that names the other kind is refused (exit 1) and the matching
one accepted; and the imported lakes of an earth-hydrology run on the cube (the default) against
the same run on lat-lon (--import-lakes refs/processed/lakes/2048.png at 128 columns): the share of
the sphere's area under imported lake within 15 % of each other. The recipes run with the main
tree as working directory (refs/ lives there; a worktree has none).
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from uwtools import paths  # noqa: E402
from uwtools.bundle import Bundle  # noqa: E402

LAKES = 'refs/processed/lakes/2048.png'
CLIMATE = ['--climate-resolution', '32', '--climate-coupling-iterations', '2']


def main(executable):
    executable = str(Path(executable).absolute())
    data_root = paths.main_root()

    def run(extra):
        return subprocess.run([executable, 'run', '--seed', '1'] + extra, cwd=data_root, capture_output=True, text=True)

    def kind(bundle):
        return json.loads((bundle / 'manifest.json').read_text(encoding='utf-8'))['grid']['kind']

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        tectonics = ['--recipe', 'tectonics', '--resolution', '32', '--evolution-duration-myr', '2']
        for name, grid in (('latlon', ['--grid', 'latlon']), ('cube', [])):
            result = run(tectonics + grid + ['--bundle', str(root / name)])
            assert result.returncode == 0, f'tectonics run on {name} failed:\n{result.stderr}'
            assert kind(root / name) == name, f'--grid {name or "default"} gave a {kind(root / name)} bundle'
        resume = ['--recipe', 'tectonics', '--only', 'geology.plate_tectonics', '--evolution-duration-myr', '2']
        for name, other in (('latlon', 'cube'), ('cube', 'latlon')):
            # Without --grid the bundle's grid applies.
            result = run(resume + ['--from-bundle', str(root / name), '--bundle', str(root / f'{name}-resumed')])
            assert result.returncode == 0, f'--from-bundle of a {name} bundle failed:\n{result.stderr}'
            assert kind(root / f'{name}-resumed') == name, f'a run from a {name} bundle changed its grid'
            # A --grid naming the other kind is refused; the bundle's own kind is accepted.
            result = run(resume + ['--from-bundle', str(root / name), '--grid', other, '--bundle', str(root / f'{name}-{other}')])
            assert result.returncode == 1, f'--grid {other} from a {name} bundle: exit {result.returncode}, not 1'
            result = run(resume + ['--from-bundle', str(root / name), '--grid', name, '--bundle', str(root / f'{name}-{name}')])
            assert result.returncode == 0, f'--grid {name} from a {name} bundle failed:\n{result.stderr}'

        shares = {}
        for name, grid in (('cube', []), ('latlon', ['--grid', 'latlon'])):
            bundle = root / f'lakes-{name}'
            result = run(['--recipe', 'earth-hydrology', '--resolution', '128', '--import-lakes', LAKES] + CLIMATE + grid +
                         ['--bundle', str(bundle)])
            assert result.returncode == 0, f'earth-hydrology with imported lakes on {name} failed:\n{result.stderr}'
            opened = Bundle.open(bundle)
            assert opened.grid_kind == name
            if name == 'cube':
                lake = opened.load_native('hydrosphere.imported_lake').reshape(-1) != 0
                area = opened.cells()['area_m2'].astype(float)
            else:
                lake = opened.load('hydrosphere.imported_lake') != 0
                rows, columns = lake.shape
                edges = np.radians(90.0 - 180.0 * np.arange(rows + 1) / rows)
                area = np.repeat((np.sin(edges[:-1]) - np.sin(edges[1:]))[:, None], columns, axis=1).ravel()
                lake = lake.ravel()
            shares[name] = (float(area[lake].sum() / area.sum()), int(lake.sum()))
        (cube, cubecells), (latlon, latloncells) = shares['cube'], shares['latlon']
        print(f'Imported lake: cube {cubecells} cells, {cube:.5f} of the sphere; lat-lon {latloncells} cells, {latlon:.5f}')
        assert latlon > 0 and abs(cube / latlon - 1) < 0.15, 'imported lake area on the cube differs from lat-lon by 15 % or more'
    print('Grid CLI passed: --from-bundle keeps the bundle grid, mismatched --grid refused, imported lakes on the cube')


if __name__ == '__main__':
    main(sys.argv[1])
