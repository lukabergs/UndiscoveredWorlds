"""Run a module experiment: write a run bundle with a headless tool, then evaluate it.

    uv run --offline python scripts/evaluation/experiment.py <module> [--register] \\
        [--name NAME] [--archive runs/reports/<experiment>] [--seed 1] [--resolution W] \\
        [--baseline BUNDLE] [--repeat] [--explicit-recipe] [--executable EXE] [--label TEXT] [--topic TEXT] \\
        [-- native options]

The module's profile (configs/evaluation/<module>.json, "run" section) holds the tool,
its arguments and its declared inputs (inputs.py); "registered" adds the settings of a
registered run (more arguments, a resolution, a recipe file whose scenario becomes flags).
Native options come last, so they override. --repeat runs the same command again
(bundle-repeat, check `repeat`); --explicit-recipe makes the recipe-equality pair
(bundle-explicit, check `recipe-equality`): a scratch run is run again with every value of
the recipe file as an option, a registered run (which carries those options) again without
them. Both second runs write no checkpoints.

A profile whose run has "seeded": false (the `regional` module: `world_pipeline region`)
takes its seed and grid from its inputs: no --seed is passed unless given (the tool then
checks it), --resolution is refused, the options of "required_options" (groups, one of
each present) must be among the arguments, every bundle directory named in the command is
declared as an input by its manifest's SHA-256, and run.json records the seed, grid and
window of the bundle written.

Without --register the run is scratch, in out/evaluation/<module>/<name>/, or in a
runs/reports/<experiment> archive (routed to the storage drive) with --archive. With
--register (`uw.py <module> run`) the run takes the next id of runs/registry/<module>.json,
runs in runs/<module>/<id>/ and appends a registry record naming that directory when it
ends, with status `failed` (and no score) when the tool or the evaluation failed or the
run was interrupted. A run directory holds run.json (command, executable and inputs with SHA-256, code
state, timing, steps, exit codes, verdict), log.txt, bundle/, checkpoints/ when the
settings ask for them, evaluation/ (evaluation.json, review request, run page, maps) and
ledger-entry.json, a draft to complete and append with `uw.py <module> ledger --add`.
Files the tool writes on its own (validation tables, images) are scratch and go to
out/evaluation/outputs/<module>/<name>/, not into the run directory.
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
from uwtools import ledger, paths, runs, storage
from uwtools.bundle import Bundle
from uwtools.hashing import sha256
from evaluate import code_state, evaluate, load_profile
import inputs
import report

ROOT = Path(__file__).resolve().parents[2]
# Environment overrides that send a tool's own outputs to scratch (OUTPUTS/<module>/<name>).
OUTPUT_VARIABLES = {
    'UW_CLIMATE_DIAGNOSTICS_DIR': 'diagnostics',
    'UW_DEFAULT_IMAGE_DIR': 'images', 'UW_DEFAULT_WORLD_DIR': 'worlds', 'UW_DEFAULT_UWS_DIR': 'worlds'}


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def relative(path):
    path = Path(path)
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path)


def recipe_options(recipe):
    """The scenario of a recipe file as --key value options (a key trench_depth is the option --trench-depth)."""
    return [part for key, value in recipe['scenario'].items() for part in ('--' + key.replace('_', '-'), str(value))]


def settings(spec, register, resolution):
    """Profile arguments, recipe options and registered arguments (with placeholders), grid
    width, placeholder values, recipe and recipe file of a run; the command is them in that order."""
    arguments, scenario, recipe, recipe_path = list(spec['arguments']), [], None, None
    values = {}
    registered = spec.get('registered', {}) if register else {}
    if registered.get('recipe'):
        recipe_path = ROOT / registered['recipe']
        recipe = json.loads(recipe_path.read_text(encoding='utf-8'))
        if recipe['height'] * 2 != recipe['width']:
            raise SystemExit('The world grid is W x W/2')
        scenario = recipe_options(recipe)
        interval = recipe['checkpoint_every_updates'] * float(recipe['scenario'].get('delta_time_myr', 1.0))
        values['checkpoint_every'] = f'{interval:g}'
    if not spec.get('seeded', True):
        # The tool takes its seed and grid from its inputs (the region command: the world bundle).
        if resolution is not None:
            raise SystemExit('This module takes its grid from its inputs; --resolution is not accepted')
        return arguments, scenario, list(registered.get('arguments', [])), None, values, recipe, recipe_path
    width = resolution or (recipe['width'] if recipe else None) or registered.get('resolution') or spec['resolution']
    if width < 2 or width % 2:
        raise SystemExit('--resolution must be an even width')
    return arguments, scenario, list(registered.get('arguments', [])), width, values, recipe, recipe_path


def expand(arguments, resolution, directory, values):
    """Profile arguments with {width}, {height}, {directory} (the run directory) and recipe values filled in."""
    values = dict(values, directory=directory.as_posix())
    if resolution is not None:
        values.update(width=resolution, height=resolution // 2)
    return [str(argument).format(**values) for argument in arguments]


def run_tool(command, directory, label, environment, outputs):
    environment = dict(environment, **{name: str(outputs / path) for name, path in OUTPUT_VARIABLES.items()})
    log = 'log.txt' if label == 'run' else f'log-{label}.txt'
    started = time.time()
    with (directory / log).open('x', encoding='utf-8') as stream:
        code = subprocess.run(command, cwd=ROOT, env=environment, stdout=stream, stderr=subprocess.STDOUT).returncode
    return {'label': label, 'command': command, 'exit_code': code, 'seconds': round(time.time() - started, 3), 'log': log}


def scratch_directory(module, name, archive):
    if archive is None:
        return ROOT / 'out/evaluation' / module / name
    return storage.route(ROOT / archive) / name


def record_in_archive(archive, module, name, record):
    path = ROOT / archive / 'experiment.json'
    document = json.loads(path.read_text(encoding='utf-8')) if path.is_file() else \
        {'schema': 'uw-experiment/v1', 'module': module, 'runs': []}
    document['runs'].append({'name': name, 'recorded': record['finished'], 'verdict': record.get('verdict')})
    path.write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')


def without_checkpoints(command):
    result, skip = [], False
    for part in command:
        if skip or part in ('--checkpoint-at', '--checkpoint-every', '--checkpoint-dir'):
            skip = not skip
            continue
        result.append(part)
    return result


def bundle_steps(bundle):
    manifest = json.loads((bundle / 'manifest.json').read_text(encoding='utf-8'))
    return [{'process': step['process'], 'label': step.get('label'), 'seconds': step.get('seconds')}
            for step in manifest.get('steps', [])]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('module', help='Profile name in configs/evaluation/')
    parser.add_argument('--register', action='store_true', help='Registered run in runs/<module>/<id>/')
    parser.add_argument('--name', help='Scratch run directory name (default: UTC time and seed)')
    parser.add_argument('--archive', help='runs/reports/<experiment> archive for kept scratch evidence')
    parser.add_argument('--seed', type=int, help='World seed (default 1; a profile with "seeded": false passes it only when given)')
    parser.add_argument('--resolution', type=int, help='Grid width W (default from the profile)')
    parser.add_argument('--baseline', type=Path, help='Bundle compared with the run')
    parser.add_argument('--repeat', action='store_true', help='Run twice and check the bundles are identical')
    parser.add_argument('--explicit-recipe', action='store_true',
                        help="Run again with the profile's recipe file as options (registered: without them) "
                             "and check the bundles are identical")
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--label', help='Change summary stored with the run')
    parser.add_argument('--topic', help='Print the ledger entries naming this text before the run (default: the last five)')
    args, native = parser.parse_known_args()
    native = [argument for argument in native if argument != '--']
    if args.register and (args.name or args.archive):
        raise SystemExit('A registered run lives in runs/<module>/<id>/; --name and --archive are for scratch runs')
    profile, profile_path = load_profile(args.module)
    if 'run' not in profile:
        raise SystemExit(f'Profile {args.module} has no run section')
    spec = profile['run']
    arguments, scenario, extra, resolution, values, recipe, recipe_path = settings(spec, args.register, args.resolution)
    seeded = spec.get('seeded', True)
    given = arguments + scenario + extra + native
    for group in spec.get('required_options', []):
        if not any(option in given for option in group):
            raise SystemExit(f"{args.module} runs require {' or '.join(group)} (a native option or a registered argument)")
    # What was already tried is shown before anything runs, so it does not depend on an instruction to read the ledger.
    previous = ledger.relevant(args.module, args.topic)
    if previous:
        print(f'Ledger ({ledger.path(args.module).relative_to(ROOT).as_posix()}):', *map(ledger.brief, previous), sep='\n', flush=True)
    # Every argument is checked before the run directory exists, so an error leaves nothing behind.
    executable = (args.executable or paths.executable(spec['tool'])).absolute()
    if not executable.is_file():
        raise SystemExit(f'Build {spec["tool"]} first: {executable}')
    explicit = []
    if args.explicit_recipe and not scenario:
        named = spec.get('registered', {}).get('recipe')
        if not named or not (ROOT / named).is_file():
            raise SystemExit(f'Profile {args.module} names no recipe file for --explicit-recipe')
        recipe_path = ROOT / named
        explicit = recipe_options(json.loads(recipe_path.read_text(encoding='utf-8')))
    if args.archive:
        runs.refuse_archive(ROOT / args.archive)
    baseline = Bundle.open(args.baseline) if args.baseline else None
    environment = dict(os.environ)
    declared = []
    for item in spec.get('inputs', []):
        variable, path = inputs.resolve(item, resolution)
        environment[variable] = str(path)
        declared.append({'path': relative(path), 'sha256': sha256(path), 'environment': variable})
    if args.register:
        identifier = runs.next_id(args.module)
        name = str(identifier)
        directory = storage.run_directory(args.module, identifier)
    else:
        identifier = None
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        name = args.name or (f'{stamp}-seed{1 if args.seed is None else args.seed}-{resolution}' if seeded else f'{stamp}-{args.module}')
        directory = scratch_directory(args.module, name, args.archive)
        directory.mkdir(parents=True, exist_ok=False)
    tool_command = lambda options: [str(executable), *expand(options, resolution, directory, values)]
    if seeded:
        seed = ['--seed', str(1 if args.seed is None else args.seed)]
    else:
        seed = [] if args.seed is None else ['--seed', str(args.seed)]
    base = [*tool_command(arguments + scenario + extra), *seed, *native]
    files = [part for part in base[1:] if not part.startswith('-') and (ROOT / part).is_file()]
    commands = {'run': base, 'repeat': without_checkpoints(base)}
    if args.explicit_recipe:
        # The recipe-equality pair: a run given every value of the recipe file as an option and a
        # run given none (the compiled defaults). A registered run carries the options itself, so
        # its second run is the one without them; a scratch run's second run adds them.
        commands['explicit'] = without_checkpoints([*tool_command(arguments + extra), *explicit, *seed, *native])
    files += [relative(recipe_path)] if recipe_path else []
    declared += [{'path': relative(ROOT / part), 'sha256': sha256(ROOT / part)} for part in files]
    if not seeded:
        # A bundle the tool reads (region --from-bundle) is declared by its manifest.
        manifests = [ROOT / part / 'manifest.json' for part in base[1:]
                     if not part.startswith('-') and (ROOT / part / 'manifest.json').is_file()]
        declared += [{'path': relative(path), 'sha256': sha256(path)} for path in manifests]

    record = {'schema': 'uw-run/v1', 'module': args.module, 'id': identifier, 'name': name, 'label': args.label,
              'seed': (1 if args.seed is None else args.seed) if seeded else args.seed,
              'grid': [resolution, resolution // 2] if seeded else None, 'tool': spec['tool'],
              'command': [*base, '--bundle', str(directory / 'bundle')],
              'executable': {'path': relative(executable), 'sha256': sha256(executable)},
              'configuration': paths.CONFIGURATION, 'code': code_state(), 'inputs': declared,
              'profile': {'path': relative(profile_path), 'sha256': sha256(profile_path)},
              'started': now(), 'runs': []}
    if recipe is not None:
        record['recipe'] = recipe
    print(f"{'Run ' + args.module + ' ' + name if args.register else 'Scratch run ' + name}: {' '.join(record['command'])}", flush=True)
    labels = ['run'] + ['repeat'] * args.repeat + ['explicit'] * args.explicit_recipe
    outputs = ROOT / 'out/evaluation/outputs' / args.module / name
    start = time.time()
    document, failure, succeeded = None, None, False
    try:
        for label in labels:
            bundle = directory / ('bundle' if label == 'run' else f'bundle-{label}')
            # The second runs check the final state only; checkpoints are not written twice.
            outcome = run_tool([*commands[label], '--bundle', str(bundle)], directory, label, environment, outputs)
            record['runs'].append(outcome)
            print(f"{label}: exit {outcome['exit_code']} in {outcome['seconds']} s", flush=True)
            if outcome['exit_code']:
                break
        record['finished'] = now()
        record['seconds'] = round(time.time() - start, 3)
        record['exit_codes'] = [outcome['exit_code'] for outcome in record['runs']]
        succeeded = all(code == 0 for code in record['exit_codes'])
        if succeeded:
            try:
                record['steps'] = bundle_steps(directory / 'bundle')
                if not seeded:
                    manifest = json.loads((directory / 'bundle/manifest.json').read_text(encoding='utf-8'))
                    record['seed'] = manifest['seed']
                    opened = Bundle.open(directory / 'bundle')
                    record['grid'] = [opened.columns, opened.rows]  # the lat-lon (export) grid in every kind
                    if opened.grid_kind != 'latlon':
                        record['grid_kind'] = opened.grid_kind
                    record['window'] = manifest.get('window')
                repeat = Bundle.open(directory / 'bundle-repeat') if args.repeat else None
                explicit_bundle = Bundle.open(directory / 'bundle-explicit') if args.explicit_recipe else None
                document = evaluate(profile, profile_path, Bundle.open(directory / 'bundle'), directory / 'evaluation',
                                    baseline=baseline, repeat=repeat, explicit=explicit_bundle)
            except Exception as error:  # the run stays recorded; evaluate again into another directory
                failure = error
                record['evaluation_error'] = f'{type(error).__name__}: {error}'
        if document:
            record['verdict'] = document['verdict']
            for check_name, key, wanted in (('repeat', 'repeat', args.repeat),
                                            ('recipe-equality', 'recipe_equality', args.explicit_recipe)):
                check = next((item for item in document['checks'] if item['check'] == check_name), None) if wanted else None
                if check:
                    record[key] = {'identical': check['status'] == 'pass', 'summary': check['summary']}
            print(json.dumps(document['verdict'], indent=2))
        (directory / 'run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
        if document:
            report.run_page(directory / 'evaluation')  # again, now with run.json
        # The entry to complete (hypothesis, expected, assessment, decision) and append with `uw.py <module> ledger --add`.
        draft = ledger.template(args.module, change={'code': record['code'], 'label': args.label, 'native_options': native},
                                observed={'run': relative(directory), 'seed': record['seed'], 'resolution': record['grid'],
                                          'verdict': record.get('verdict')})
        (directory / 'ledger-entry.json').write_text(json.dumps(draft, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
        if args.archive:
            record_in_archive(args.archive, args.module, name, record)
    except BaseException as error:
        # An interrupted or failed run still explains its directory.
        record.setdefault('finished', now())
        record['error'] = f'{type(error).__name__}: {error}'
        (directory / 'run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
        raise
    finally:
        if args.register:
            # Every registered id gets its record, so it is never reused; only a run evaluated to the end is complete.
            complete = succeeded and document is not None and 'error' not in record
            entry = {'id': identifier, 'status': 'complete' if complete else 'failed', 'run': relative(directory),
                     'date': record['started'], 'label': args.label or ''}
            score = document['validation']['scores'].get('run', {}).get('score') if complete else None
            if score is not None:
                entry['score'] = round(score, 4)
            runs.append_record(args.module, entry)
            print(f"Registered {args.module} run {identifier} ({entry['status']}) in "
                  f"{runs.registry_path(args.module).relative_to(ROOT).as_posix()}")
    print(f'Run directory: {directory}')
    if not succeeded:
        raise SystemExit(f"Tool failed; see {directory / record['runs'][-1]['log']}")
    if failure:
        raise failure
    sys.exit(1 if record['verdict']['checks'] == 'fail' else 0)


if __name__ == '__main__':
    main()
