"""Check of one refactoring step at the minimal configuration.

    uv run --offline scripts/evaluation/step_check.py TAG [--reference REF] [--desktop] [--social]
                                                       [--grid cube|latlon]

Builds the headless preset (and the desktop preset with --desktop), makes the three seed-1
bundles at 128 export columns (tectonics 50 Myr; climate on 32 columns with 2 coupling iterations;
--audit-access) as <directory>/<recipe>, on the world_pipeline default grid, the cube of
N 32 (--grid latlon: on 128 x 64 lat-lon, for lat-lon parity), compares each with the reference's
bundle, then runs parameters.py --check, check_layout.py, check_structure.py, check_records.py and
process_access.py --check on the stored report sets. --social adds a procedural world with society (<directory>/social).
The directory is the main tree's out/scratch/step-check/<id>-<TAG> (scripts/uwtools/naming.py: a new id
each time, recorded in runs/registry/scratch.json; TAG names the topic, without a version suffix). Build
and run output goes to <directory>/*.log; <directory>/result.json records passed, grid, social and the
reference. Exit status: 0 when every build, run and check passed and every comparison was identical.

REF is an id, an <id>-<tag> name, a label of the earlier layout (flat <label>-<recipe> bundles in
out/scratch/step-check, out/scratch/sphere/step-check or out/scratch/phase5) or an absolute path.
Without --reference the newest step check with the same grid and --social setting that passed is the
reference; without one, nothing is compared.

From a git worktree (.claude/worktrees/<name>, no refs/ or runs/) the worktree is built and
checked as above, but the recipes run with the main tree as working directory (so refs/
resolve), using the worktree's world_pipeline.exe, and process_access.py --check reads the
main tree's runs/reports; the directory is the main tree's in any case.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import naming, paths

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


def reference_bundle(reference, name):
    """The reference's bundle of one recipe: <directory>/<recipe>, or <prefix>-<recipe> for a flat step check of
    the earlier layout (naming.find returns its prefix)."""
    return reference / name if reference.is_dir() else reference.with_name(f'{reference.name}-{name}')


def default_reference(grid, social):
    """The newest step check with the same grid and --social setting whose result.json says it passed, or None."""
    root = naming.main_root()
    for entry in reversed(naming.records('step-check')):
        result_path = root / entry['run'] / 'result.json'
        if not result_path.is_file():
            continue
        result = json.loads(result_path.read_text(encoding='utf-8'))
        if result.get('passed') and result.get('grid') == grid and bool(result.get('social')) == social:
            return root / entry['run']
    return None


def shown(path):
    """A path relative to the main tree when it lies there, else as given."""
    root = naming.main_root()
    path = Path(os.path.abspath(path))
    return path.relative_to(root).as_posix() if path.is_relative_to(root) else str(path)


def checks(arguments, directory, reference, result):
    """Builds, runs, comparisons and the repository checks; True when every one passed."""
    if not build('headless', directory / 'build-headless.log'):
        return False
    passed = True
    if arguments.desktop:
        passed = build('desktop', directory / 'build-desktop.log')

    executable = paths.executable('world_pipeline')
    runs = dict(RECIPES)
    grid = ['--grid', 'latlon'] if arguments.grid == 'latlon' else []
    if arguments.social:
        runs.update(SOCIAL)
    for name, options in runs.items():
        bundle = directory / name
        start = time.time()
        status = logged([str(executable), 'run'] + options + COMMON + grid + ['--bundle', str(bundle)], directory / f'{name}.log',
                        cwd=DATA_ROOT)
        print(f'{name} run exit {status} in {time.time() - start:.0f} s')
        result['runs'][name] = status
        passed = passed and status == 0
        if reference is not None:
            log = directory / f'{name}-compare.log'
            status = logged([sys.executable, 'scripts/evaluation/compare_bundles.py',
                             str(reference_bundle(reference, name)), str(bundle)], log)
            lines = log.read_text(encoding='utf-8', errors='replace').strip().splitlines()
            print(f'{name} compare exit {status}: {lines[-1] if lines else ""}')
            result['comparisons'][name] = status
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
    return passed and status == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tag', help='the topic of the check (naming.py tag rule); the directory is <id>-<tag>')
    parser.add_argument('--reference', metavar='REF',
                        help='an id, <id>-<tag>, a label of the earlier layout or an absolute path '
                             '(default: the newest passing step check with the same --grid and --social)')
    parser.add_argument('--desktop', action='store_true')
    parser.add_argument('--social', action='store_true')
    parser.add_argument('--grid', choices=('cube', 'latlon'), default='cube',
                        help='the world grid (default cube, the world_pipeline default; latlon for lat-lon parity)')
    arguments = parser.parse_args()
    try:
        naming.check_tag(arguments.tag)
        reference = naming.find('step-check', arguments.reference) if arguments.reference else \
            default_reference(arguments.grid, arguments.social)
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    named = naming.new_directory('step-check', arguments.tag,
                                 label=f"{arguments.grid}{', social' if arguments.social else ''}")
    directory = named.path
    print(f'Step check {named.name}: {directory}')
    print(f'Reference: {shown(reference) if reference is not None else "none (nothing compared)"}')
    result = {'schema': 'uw-step-check/v1', 'id': named.id, 'tag': named.tag, 'name': named.name,
              'grid': arguments.grid, 'social': arguments.social, 'desktop': arguments.desktop,
              'reference': shown(reference) if reference is not None else None,
              'passed': False, 'runs': {}, 'comparisons': {}}
    try:
        result['passed'] = checks(arguments, directory, reference, result)
    finally:
        (directory / 'result.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
