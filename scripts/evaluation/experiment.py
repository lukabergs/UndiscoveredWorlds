"""Run a module experiment: write a run bundle with a headless tool, then evaluate it.

    uv run --offline python scripts/evaluation/experiment.py <module> \\
        [--name NAME] [--archive runs/reports/<experiment>] [--seed 1] [--resolution W] \\
        [--baseline BUNDLE] [--repeat] [--executable EXE] [-- native options]

The module's profile (configs/evaluation/<module>.json) names the tool and its default
arguments. Runs go to out/evaluation/<module>/<name>/ (scratch) unless --archive names a
runs/reports/<experiment> archive, which is routed to the data drive first. A run
directory holds run.json (command, seed, timing, code state), run.log, outputs/ (every
file the tool would otherwise write under runs/), bundle/ and evaluation/. Modules with
their own registry (climate, tectonics) register runs through uw.py, not here.
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths
from uwtools.bundle import Bundle
from evaluate import code_state, evaluate, load_profile
import report

ROOT = Path(__file__).resolve().parents[2]
# Environment overrides that keep a tool's own outputs inside the run directory.
OUTPUT_VARIABLES = {
    'UW_CLIMATE_VALIDATION_DIR': 'validation', 'UW_CLIMATE_BENCHMARK_IMAGE_DIR': 'maps',
    'UW_CLIMATE_BENCHMARK_RUN_DIR': 'diagnostics', 'UW_CLIMATE_BENCHMARK_RUN_LOG': 'registry.json',
    'UW_PROFILING_WORKBOOK': 'profiling.xlsx', 'UW_CLIMATE_WORKBOOK': 'climate.xlsx',
    'UW_DEFAULT_IMAGE_DIR': 'images', 'UW_DEFAULT_WORLD_DIR': 'worlds', 'UW_DEFAULT_UWS_DIR': 'worlds'}


def expand(arguments, resolution, directory):
    """Profile arguments with {width}, {height} and {directory} (the run directory) filled in."""
    values = {'width': resolution, 'height': resolution // 2, 'directory': directory.as_posix()}
    return [str(argument).format(**values) for argument in arguments]


def run_tool(command, directory, label):
    outputs = directory / 'outputs'
    environment = dict(os.environ, **{name: str(outputs / path) for name, path in OUTPUT_VARIABLES.items()})
    started = time.time()
    with (directory / f'{label}.log').open('x', encoding='utf-8') as log:
        code = subprocess.run(command, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT).returncode
    return {'command': command, 'exit_code': code, 'seconds': round(time.time() - started, 3), 'log': f'{label}.log'}


def run_directory(module, name, archive):
    if archive is None:
        return ROOT / 'out/evaluation' / module / name
    sys.path.insert(0, str(ROOT / 'scripts/benchmarks'))
    from run_storage import ensure_archive_directory
    archive = ensure_archive_directory(ROOT / archive)
    archive.mkdir(parents=True, exist_ok=True)
    return archive / name


def record_in_archive(archive, module, name, record):
    path = ROOT / archive / 'experiment.json'
    document = json.loads(path.read_text(encoding='utf-8')) if path.is_file() else \
        {'schema': 'uw-experiment/v1', 'module': module, 'runs': []}
    document['runs'].append({'name': name, 'recorded': record['finished'], 'verdict': record.get('verdict')})
    path.write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('module', help='Profile name in configs/evaluation/')
    parser.add_argument('--name', help='Run directory name (default: UTC time and seed)')
    parser.add_argument('--archive', help='runs/reports/<experiment> archive for kept evidence')
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--resolution', type=int, help='Grid width W (default from the profile)')
    parser.add_argument('--baseline', type=Path, help='Bundle compared with the run')
    parser.add_argument('--repeat', action='store_true', help='Run twice and check the bundles are identical')
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--label', help='Change summary stored with the run')
    args, native = parser.parse_known_args()
    native = [argument for argument in native if argument != '--']
    profile, profile_path = load_profile(args.module)
    if 'run' not in profile:
        raise SystemExit(f'Profile {args.module} has no run section')
    spec = profile['run']
    resolution = args.resolution or spec['resolution']
    if resolution < 2 or resolution % 2:
        raise SystemExit('--resolution must be an even width')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    name = args.name or f'{stamp}-seed{args.seed}-{resolution}'
    directory = run_directory(args.module, name, args.archive)
    directory.mkdir(parents=True, exist_ok=False)
    executable = args.executable or paths.executable(spec['tool'])
    base = [str(executable), *expand(spec['arguments'], resolution, directory), '--seed', str(args.seed), *native]

    record = {'schema': 'uw-experiment-run/v1', 'module': args.module, 'name': name, 'label': args.label,
              'seed': args.seed, 'resolution': [resolution, resolution // 2], 'tool': spec['tool'],
              'executable': str(executable), 'configuration': paths.CONFIGURATION, 'code': code_state(),
              'started': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), 'runs': []}
    labels = ['run', 'repeat'] if args.repeat else ['run']
    for label in labels:
        bundle = directory / ('bundle' if label == 'run' else 'bundle-repeat')
        outcome = run_tool([*base, '--bundle', str(bundle)], directory, label)
        record['runs'].append(outcome)
        print(f"{label}: exit {outcome['exit_code']} in {outcome['seconds']} s", flush=True)
        if outcome['exit_code']:
            break
    record['finished'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    succeeded = all(outcome['exit_code'] == 0 for outcome in record['runs'])
    if succeeded:
        repeat = Bundle.open(directory / 'bundle-repeat') if args.repeat else None
        baseline = Bundle.open(args.baseline) if args.baseline else None
        document = evaluate(profile, profile_path, Bundle.open(directory / 'bundle'), directory / 'evaluation',
                            baseline=baseline, repeat=repeat)
        record['verdict'] = document['verdict']
        print(json.dumps(document['verdict'], indent=2))
    (directory / 'run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    if succeeded:
        report.run_page(directory / 'evaluation')  # again, now with run.json
    if args.archive:
        record_in_archive(args.archive, args.module, name, record)
    print(f'Run directory: {directory}')
    if not succeeded:
        raise SystemExit(f"Tool failed; see {directory / record['runs'][-1]['log']}")
    sys.exit(1 if record['verdict']['checks'] == 'fail' else 0)


if __name__ == '__main__':
    main()
