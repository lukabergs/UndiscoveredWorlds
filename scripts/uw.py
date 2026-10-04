"""Common entry point for module benchmarks.

    uv run --offline scripts/uw.py <module> <command> [shared options] [native options]

Modules and commands:
    tectonics run      registered run (runs/registry/tectonics.json): world_pipeline recipe with checkpoints,
                       archived under runs/reports|checkpoints/tectonics/<id>
    tectonics inspect  unregistered final-state run bundle in out/bench/
    tectonics report   HTML comparison of registered runs (--runs A B ...)
    tectonics experiment|evaluate  run bundle + evaluation with configs/evaluation/tectonics.json
    geology   run      tectonics bundles at each width, scored against Earth references
    climate   run      registered Earth climate benchmark (runs/registry/climate.json) with gallery
    climate   check    numerical checks of a run (--runs ID [--compare ID])
    climate   report   gallery of a run against another (--runs ID [--compare ID])
    world     run      full procedural world (terrain, climate, optional society) in out/bench/
    climate|world experiment  run bundle + evaluation with configs/evaluation/<module>.json
                       (scripts/evaluation/experiment.py; out/evaluation/ unless --archive)
    climate|world evaluate    evaluate an existing bundle (--bundle DIR [--baseline DIR] [--checkpoints DIR])
    climate|tectonics|world report --evaluation DIR [DIR ...] [--output DIR]
                       pages from evaluation.json (scripts/evaluation/report.py): one DIR rewrites its run
                       page (index.html, report.md); several give a comparison page, the first DIR being the
                       baseline, in --output (default out/evaluation/<module>/comparison)
    <module>  ledger   experiment ledger (runs/ledger/<module>.jsonl; also module `world`):
                       ledger [--topic TEXT]   entries naming the topic, or the last five
                       ledger --draft          a new entry with every key, printed for filling in
                       ledger --add FILE.json  validate and append one entry; LEDGER.md is regenerated
    climate|tectonics|world dashboard [--output DIR]
                       registry runs and ledger entries (default out/evaluation/<module>/dashboard)

Shared options: --seed (default 1), --resolution W (grid W x W/2), --executable, --output, --label.
Any other --flag value pairs go to the module's native tool unchanged: world_pipeline
flags for every module (tectonics: the scenario options of `run --recipe tectonics`).
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / 'benchmarks'))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools import ledger, paths
import climate_benchmark

# `uv run` exports its own path; subprocesses may not find uv on PATH.
UV = os.environ.get('UV', 'uv')
ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts/benchmarks'
WORLD_EXE = paths.executable('world_pipeline')
RECIPE = ROOT / 'configs/tectonics-baseline.json'
SCRATCH = ROOT / 'out/bench'
EVALUATION = ROOT / 'scripts/evaluation'
COMMANDS = {'tectonics': ('run', 'inspect', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger'), 'geology': ('run',),
            'climate': ('run', 'check', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger'),
            'world': ('run', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger')}


def call(command):
    print('>', ' '.join(map(str, command)), flush=True)
    code = subprocess.run([str(part) for part in command], cwd=ROOT).returncode
    if code:
        raise SystemExit(code)


def python_script(name, *arguments):
    call([UV, 'run', '--offline', SCRIPTS / name, *arguments])


def evaluation(args, extra):
    """Profile-driven runs and evaluations (scripts/evaluation/, configs/evaluation/<module>.json)."""
    if args.command == 'evaluate':
        if not args.bundle:
            raise SystemExit(f'{args.module} evaluate requires --bundle DIR')
        output = args.output or args.bundle.parent / 'evaluation'
        call([UV, 'run', '--offline', EVALUATION / 'evaluate.py', '--profile', args.module,
              '--bundle', args.bundle, '--output', output, *extra])
        return
    options = ['--seed', args.seed]
    options += ['--resolution', args.resolution] if args.resolution is not None else []
    options += ['--executable', args.executable] if args.executable else []
    options += ['--label', args.label] if args.label else []
    options += ['--topic', args.topic] if args.topic else []
    call([UV, 'run', '--offline', EVALUATION / 'experiment.py', args.module, *options, *extra])


def ledger_command(args, extra):
    """Show, draft or append entries of runs/ledger/<module>.jsonl."""
    if extra:
        raise SystemExit(f"Unexpected arguments: {' '.join(extra)}")
    if args.add:
        entry = ledger.append(json.loads(args.add.read_text(encoding='utf-8')))
        print(f"Appended {entry['id']} to {ledger.path(entry['module']).relative_to(ROOT)}")
    elif args.draft:
        print(json.dumps(ledger.template(args.module), indent=1, ensure_ascii=False))
    else:
        found = ledger.relevant(args.module, args.topic)
        print('\n'.join(map(ledger.brief, found)) or f'No {args.module} ledger entry' + (f' names {args.topic!r}' if args.topic else ''))


def pages(args, extra):
    """Run, comparison and dashboard pages rendered from JSON by scripts/evaluation/report.py."""
    if extra:
        raise SystemExit(f"Unexpected arguments: {' '.join(extra)}")
    script = EVALUATION / 'report.py'
    if args.command == 'dashboard':
        output = args.output or ROOT / 'out/evaluation' / args.module / 'dashboard'
        call([UV, 'run', '--offline', script, 'dashboard', args.module, '--output', output])
    elif len(args.evaluation) == 1:
        call([UV, 'run', '--offline', script, 'run', args.evaluation[0]])
    else:
        output = args.output or ROOT / 'out/evaluation' / args.module / 'comparison'
        call([UV, 'run', '--offline', script, 'compare', *args.evaluation, '--output', output])


def native_pairs(extra):
    """Parse pass-through `--flag value` pairs; bare flags get value None."""
    pairs, index = [], 0
    while index < len(extra):
        flag = extra[index]
        if not flag.startswith('--'):
            raise SystemExit(f'Unexpected argument: {flag}')
        value = extra[index + 1] if index + 1 < len(extra) and not extra[index + 1].startswith('--') else None
        pairs.append((flag, value))
        index += 1 if value is None else 2
    return pairs


def flatten(pairs):
    return [part for flag, value in pairs for part in ((flag,) if value is None else (flag, value))]


def tectonics_bundle(exe, seed, resolution, extra, output):
    """Unregistered final state of the tectonics recipe as a run bundle."""
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f'Output exists; bundles are never overwritten: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    grid = [] if resolution is None else ['--resolution', resolution]
    call([exe, 'run', '--recipe', 'tectonics', '--seed', seed, *grid, *extra, '--bundle', output])


def tectonics(args, extra):
    exe = args.executable or WORLD_EXE
    if args.command == 'report':
        if not args.runs:
            raise SystemExit('tectonics report requires --runs')
        python_script('tectonics_report.py', '--runs', *map(str, args.runs), '--external-assets',
                      *(['--output', args.output] if args.output else []), *extra)
        compare = args.runs[-2] if len(args.runs) > 1 else None
        python_script('run_report.py', 'tectonics', '--run', str(args.runs[-1]), *(['--compare', str(compare)] if compare else []))
        return
    if args.command == 'inspect':
        output = Path(args.output or SCRATCH / f'tectonics-seed{args.seed}-{args.resolution or 256}')
        tectonics_bundle(exe, args.seed, args.resolution, extra, output)
        return
    recipe = json.loads(RECIPE.read_text(encoding='utf-8'))
    recipe['seed'] = args.seed
    if args.resolution is not None:
        recipe['width'], recipe['height'] = args.resolution, args.resolution // 2
    for flag, value in native_pairs(extra):
        if value is None:
            raise SystemExit(f'Recipe flags need a value: {flag}')
        key = flag[2:]
        if key in ('seed', 'width', 'height', 'checkpoint_every_updates'):
            recipe[key] = int(value)
        else:
            recipe['scenario'][key] = value
    arguments = ['--executable', exe]
    if recipe != json.loads(RECIPE.read_text(encoding='utf-8')):
        path = SCRATCH / f'tectonics-recipe-{time.strftime("%Y%m%d-%H%M%S")}.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(recipe, indent=2) + '\n', encoding='utf-8')
        arguments += ['--recipe', path]
    if args.seed != 1:
        arguments.append('--allow-seed-comparison')
    python_script('tectonics_experiment.py', *arguments)
    # Same outputs as a climate run: maps, animations, comparison report and run page.
    registry = json.loads((ROOT / 'runs/registry/tectonics.json').read_text(encoding='utf-8'))['runs']
    complete = [r['id'] for r in registry if r.get('status') == 'complete']
    run_id = max(complete)
    earlier = [run for run in complete if run < run_id]
    runs = [max(earlier), run_id] if earlier else [run_id]
    python_script('tectonics_report.py', '--runs', *map(str, runs), '--external-assets')
    python_script('run_report.py', 'tectonics', '--run', str(run_id), *(['--compare', str(runs[0])] if earlier else []))


def geology(args, extra):
    exe = args.executable or WORLD_EXE
    output = Path(args.output or SCRATCH / 'geology')
    widths = [args.resolution] if args.resolution else [128, 256]
    bundles = []
    for width in widths:
        bundle = output / f'seed-{args.seed}-{width}'
        tectonics_bundle(exe, args.seed, width, extra, bundle)
        bundles.append(bundle)
    python_script('benchmark_geology.py', '--bundle', *bundles, '--width', str(args.reference_width),
                  '--output', output / 'report.json')


def climate(args, extra):
    if args.command in ('check', 'report'):
        if not args.runs or len(args.runs) != 1:
            raise SystemExit(f'climate {args.command} requires --runs ID')
        compare = ['--compare-run-id', str(args.compare)] if args.compare else []
        if args.command == 'check':
            python_script('check_climate_baseline.py', '--run-id', str(args.runs[0]), *compare, *extra)
            return
        climate_pages(args.runs[0], args.compare)
        return
    arguments = ['--seed', str(args.seed), *extra]
    if args.resolution is not None:
        arguments += ['--resolution', str(args.resolution)]
    run_id = climate_benchmark.run(arguments, args.executable or climate_benchmark.EXECUTABLE, args.label or '',
                                   args.surface_drag)
    climate_pages(run_id, args.compare)


def climate_pages(run_id, compare=None):
    """Standard run page plus the side-by-side gallery against observed references."""
    python_script('run_report.py', 'climate', '--run', str(run_id), *(['--compare', str(compare)] if compare else []))
    gallery_compare = compare or run_id - 1
    python_script('write_benchmark_report.py', '--run-id', str(run_id), '--legacy-layout',
                  '--compare-run-id', str(gallery_compare))
    # The gallery refreshes the record; rewrite the run page so its links include the gallery.
    python_script('run_report.py', 'climate', '--run', str(run_id), *(['--compare', str(compare)] if compare else []))


def world(args, extra):
    output = Path(args.output or SCRATCH / f'world-seed{args.seed}-{args.resolution or "default"}.uww')
    output.parent.mkdir(parents=True, exist_ok=True)
    resolution = ['--resolution', str(args.resolution)] if args.resolution is not None else []
    call([args.executable or WORLD_EXE, 'generate', '--seed', args.seed, *resolution, *extra, '--save', output])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     allow_abbrev=False)
    parser.add_argument('module', choices=COMMANDS)
    parser.add_argument('command')
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--resolution', type=int, help='Grid width W; height is W/2')
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--label', help='Change summary stored with the run')
    parser.add_argument('--runs', type=int, nargs='+')
    parser.add_argument('--compare', type=int, help='Run ID to compare against (climate check/report)')
    parser.add_argument('--bundle', type=Path, help='Run bundle to evaluate (evaluate)')
    parser.add_argument('--evaluation', type=Path, nargs='+',
                        help='Evaluation directories for report pages (report); the first is the baseline')
    parser.add_argument('--topic', help='Ledger entries naming this text (ledger; printed before an experiment)')
    parser.add_argument('--draft', action='store_true', help='ledger: print a new entry to fill in')
    parser.add_argument('--add', type=Path, help='ledger: JSON file of one entry to append')
    parser.add_argument('--reference-width', type=int, default=512, help='Earth reference width (geology)')
    parser.add_argument('--surface-drag', default='earth',
                        help="climate: 'earth' vegetation drag contrast (default), 'none', or a UWROUGH1 file")
    args, extra = parser.parse_known_args()
    if args.command not in COMMANDS[args.module]:
        parser.error(f'{args.module} commands: {", ".join(COMMANDS[args.module])}')
    if args.resolution is not None and (args.resolution < 2 or args.resolution % 2):
        parser.error('--resolution must be an even width')
    if args.command == 'ledger':
        ledger_command(args, extra)
        return
    if args.command in ('experiment', 'evaluate'):
        evaluation(args, extra)
        return
    if args.command == 'dashboard' or args.evaluation:
        if args.command not in ('report', 'dashboard'):
            parser.error('--evaluation belongs to report')
        pages(args, extra)
        return
    if args.module == 'world' and args.command == 'report':
        parser.error('world report requires --evaluation DIR [DIR ...]')
    {'tectonics': tectonics, 'geology': geology, 'climate': climate, 'world': world}[args.module](args, extra)


if __name__ == '__main__':
    main()
