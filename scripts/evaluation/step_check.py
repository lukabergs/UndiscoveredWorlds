"""Check of one refactoring step at the minimal configuration.

    uv run --offline scripts/evaluation/step_check.py LABEL [--reference LABEL] [--desktop] [--social]
                                                       [--grid cube|latlon] [--directory out/scratch/phase5]

Builds the headless preset (and the desktop preset with --desktop), makes the three seed-1
bundles at 128 export columns (tectonics 50 Myr; climate on 32 columns with 2 coupling iterations;
--audit-access) as <directory>/<LABEL>-<recipe>, on the world_pipeline default grid, the cube of
N 32 (--grid latlon: on 128 x 64 lat-lon, for lat-lon parity), compares each with the bundle of the
reference label, then runs parameters.py --check, check_layout.py, check_structure.py, check_records.py and
process_access.py --check on the stored report sets. --social adds a procedural world with society (<LABEL>-social).
Build and run output goes to <directory>/<LABEL>-*.log. Exit status: 0 when every build,
run and check passed and every comparison was identical.

From a git worktree (.claude/worktrees/<name>, no refs/ or runs/) the worktree is built and
checked as above, but the recipes run with the main tree as working directory (so refs/
resolve), using the worktree's world_pipeline.exe, and process_access.py --check reads the
main tree's runs/reports. Pass --directory as an absolute path in the main tree's out/scratch.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths

ROOT = paths.ROOT
# The main tree: ROOT itself, or the main checkout when ROOT is a worktree; refs/ and runs/ live there.
DATA_ROOT = paths.main_root()
COMMON = ['--seed', '1', '--resolution', '128', '--audit-access']
CLIMATE = ['--climate-resolution', '32', '--climate-coupling-iterations', '2']
RECIPES = {
    'procedural-world': ['--recipe', 'procedural-world', '--duration', '50'] + CLIMATE,
    'earth-climate': ['--recipe', 'earth-climate'] + CLIMATE,
    'tectonics': ['--recipe', 'tectonics', '--duration', '50'],
}
SOCIAL = {'social': RECIPES['procedural-world'] + ['--social']}
ACCESS_SETS = ('phase4-access', 'phase4b-access', 'phase4c-access', 'erosion-a-access', 'erosion-b-access',
               'erosion-c-access', 'erosion-d-access', 'erosion-g-access',
               'multiscale-s3-access', 'multiscale-s5-access', 't412-access', 't412-sea-access',
               't8-w2t2-access', 't8-e3-access', 't8-w3-access')


def logged(command, log, cwd=ROOT):
    with open(log, 'w', encoding='utf-8', errors='replace') as stream:
        return subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT).returncode


def build(preset, log):
    status = logged(['cmake', '--build', f'out/build/{preset}', '--config', paths.CONFIGURATION, '--parallel', '8'], log)
    print(f'{preset} build exit {status}')
    return status == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('label')
    parser.add_argument('--reference')
    parser.add_argument('--desktop', action='store_true')
    parser.add_argument('--social', action='store_true')
    parser.add_argument('--grid', choices=('cube', 'latlon'), default='cube',
                        help='the world grid (default cube, the world_pipeline default; latlon for lat-lon parity)')
    parser.add_argument('--directory', default='out/scratch/phase5')
    arguments = parser.parse_args()
    directory = ROOT / arguments.directory
    directory.mkdir(parents=True, exist_ok=True)
    label = arguments.label
    if any(directory.glob(f'{label}-*')):
        print(f'{arguments.directory}/{label}-* exists: a label is used once, so that bundles are never overwritten')
        return 1

    if not build('headless', directory / f'{label}-build-headless.log'):
        return 1
    passed = True
    if arguments.desktop:
        passed = build('desktop', directory / f'{label}-build-desktop.log')

    executable = paths.executable('world_pipeline')
    runs = dict(RECIPES)
    grid = ['--grid', 'latlon'] if arguments.grid == 'latlon' else []
    if arguments.social:
        runs.update(SOCIAL)
    for name, options in runs.items():
        bundle = directory / f'{label}-{name}'
        start = time.time()
        status = logged([str(executable), 'run'] + options + COMMON + grid + ['--bundle', str(bundle)], directory / f'{label}-{name}.log',
                        cwd=DATA_ROOT)
        print(f'{name} run exit {status} in {time.time() - start:.0f} s')
        passed = passed and status == 0
        if arguments.reference:
            log = directory / f'{label}-{name}-compare.log'
            status = logged([sys.executable, 'scripts/evaluation/compare_bundles.py',
                             str(directory / f'{arguments.reference}-{name}'), str(bundle)], log)
            lines = log.read_text(encoding='utf-8', errors='replace').strip().splitlines()
            print(f'{name} compare exit {status}: {lines[-1] if lines else ""}')
            passed = passed and status == 0

    status = subprocess.run([sys.executable, 'scripts/parameters.py', '--check'], cwd=ROOT).returncode
    passed = passed and status == 0
    status = subprocess.run([sys.executable, 'scripts/check_layout.py'], cwd=ROOT).returncode
    print(f'check_layout exit {status}')
    passed = passed and status == 0
    status = subprocess.run([sys.executable, 'scripts/check_structure.py'], cwd=ROOT).returncode
    print(f'check_structure exit {status}')
    passed = passed and status == 0
    status = subprocess.run([sys.executable, 'scripts/check_records.py'], cwd=ROOT).returncode
    print(f'check_records exit {status}')
    passed = passed and status == 0
    reports = [str(path) for name in ACCESS_SETS for path in sorted((DATA_ROOT / 'runs/reports' / name).glob('*.json'))]
    status = subprocess.run([sys.executable, 'scripts/evaluation/process_access.py', '--check'] + reports, cwd=ROOT).returncode
    print(f'process_access exit {status}')
    return 0 if passed and status == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
