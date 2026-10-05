"""Common entry point for module benchmarks.

    uv run --offline scripts/uw.py <module> <command> [shared options] [native options]

Modules and commands:
    climate|tectonics|world run  registered run: the next id of runs/registry/<module>.json, run in runs/<module>/<id>/
                       with the registered settings of configs/evaluation/<module>.json (run.json, log.txt,
                       bundle/, checkpoints/, evaluation/, ledger-entry.json; scripts/evaluation/experiment.py)
    climate|tectonics|world experiment  unregistered scratch run + evaluation with the profile's settings
                       (out/evaluation/<module>/ unless --archive runs/reports/<name>); --repeat runs twice
                       and checks that the bundles are identical
    climate|tectonics|world evaluate    evaluate an existing bundle (--bundle DIR [--baseline DIR] [--checkpoints DIR])
    tectonics report   comparison page of the evaluations of registered runs (--runs A B ...), in runs/reports/tectonics/
    geology   run      tectonics bundles at each width, scored against Earth references (out/bench/geology)
    climate|tectonics|world report --evaluation DIR [DIR ...] [--output DIR]
                       pages from evaluation.json (scripts/evaluation/report.py): one DIR rewrites its run
                       page (index.html, report.md); several give a comparison page, the first DIR being the
                       baseline, in --output (default out/evaluation/<module>/comparison)
    <module>  review   review.json of an evaluation: review --evaluation DIR [--verdict approved|rejected|revise
                       --notes TEXT --reviewer user|agent --answers FILE.json]; the user's verdict overrides
    <module>  ledger   experiment ledger (runs/ledger/<module>.jsonl; also module `world`):
                       ledger [--topic TEXT]   entries naming the topic, or the last five
                       ledger --draft          a new entry with every key, printed for filling in
                       ledger --add FILE.json  validate and append one entry; LEDGER.md is regenerated
    climate|tectonics|world dashboard [--output DIR]
                       registry runs and ledger entries (default out/evaluation/<module>/dashboard)

Shared options: --seed (default 1), --resolution W (grid W x W/2), --executable, --output, --label.
Any other --flag value pairs go to the module's native tool unchanged and after the profile's
arguments, so they override them: world_pipeline flags for every module (tectonics: the
scenario options of `run --recipe tectonics`, and --duration).
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools import ledger, paths, runs

# `uv run` exports its own path; subprocesses may not find uv on PATH.
UV = os.environ.get('UV', 'uv')
ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts/benchmarks'
WORLD_EXE = paths.executable('world_pipeline')
SCRATCH = ROOT / 'out/bench'
EVALUATION = ROOT / 'scripts/evaluation'
COMMANDS = {'tectonics': ('run', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger', 'review'), 'geology': ('run',),
            'climate': ('run', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger', 'review'),
            'world': ('run', 'report', 'experiment', 'evaluate', 'dashboard', 'ledger', 'review')}


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
        output = runs.refuse_archive(args.output or args.bundle.parent / 'evaluation')
        call([UV, 'run', '--offline', EVALUATION / 'evaluate.py', '--profile', args.module,
              '--bundle', args.bundle, '--output', output, *extra])
        return
    options = ['--register'] if args.command == 'run' else []
    options += ['--seed', args.seed]
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
        entry = json.loads(args.add.read_text(encoding='utf-8'))
        if entry.get('module') != args.module:
            raise SystemExit(f"{args.add} is a {entry.get('module')!r} ledger entry, not {args.module!r}")
        entry = ledger.append(entry)
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


def tectonics_bundle(exe, seed, resolution, extra, output):
    """Unregistered final state of the tectonics recipe as a run bundle."""
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f'Output exists; bundles are never overwritten: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    grid = [] if resolution is None else ['--resolution', resolution]
    call([exe, 'run', '--recipe', 'tectonics', '--seed', seed, *grid, *extra, '--bundle', output])


def tectonics(args, extra):
    """Comparison page of registered tectonics runs, from their evaluations (the first run is the baseline).

    A sealed run is evaluated once into runs/reports/tectonics/<id>/evaluation when its final state is a
    run bundle; runs up to 75 are raw snapshots and cannot be evaluated (run 77 holds the fields of 75)."""
    if not args.runs:
        raise SystemExit('tectonics report requires --runs')
    if extra:
        raise SystemExit(f"Unexpected arguments: {' '.join(extra)}")
    directories = []
    for identifier in args.runs:
        run = runs.directory('tectonics', identifier)
        directory = run / 'evaluation' if run else ROOT / 'runs/reports/tectonics' / str(identifier) / 'evaluation'
        if not (directory / 'evaluation.json').is_file():
            if run:
                raise SystemExit(f'Run {identifier} has no evaluation: {directory}')
            final = runs.tectonics(identifier)['snapshots'][-1]['path']
            if not (final / 'manifest.json').is_file() or 'uw-run-bundle' not in (final / 'manifest.json').read_text(encoding='utf-8')[:400]:
                raise SystemExit(f'Run {identifier} is stored as raw snapshots, not as a run bundle, and cannot be evaluated; '
                                 'compare a run of the bundle format instead (run 77 holds the fields of runs 75 and 76)')
            call([UV, 'run', '--offline', EVALUATION / 'evaluate.py', '--profile', 'tectonics', '--bundle', final,
                  '--output', directory, '--checkpoints', final.parent])
        directories.append(directory)
    if len(directories) == 1:
        call([UV, 'run', '--offline', EVALUATION / 'report.py', 'run', directories[0]])
        return
    output = args.output or ROOT / 'runs/reports/tectonics' / ('-'.join(map(str, args.runs)) + '-comparison')
    call([UV, 'run', '--offline', EVALUATION / 'report.py', 'compare', *directories, '--output', output])


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
    parser.add_argument('--bundle', type=Path, help='Run bundle to evaluate (evaluate)')
    parser.add_argument('--evaluation', type=Path, nargs='+',
                        help='Evaluation directories for report pages (report); the first is the baseline')
    parser.add_argument('--topic', help='Ledger entries naming this text (ledger; printed before an experiment)')
    parser.add_argument('--draft', action='store_true', help='ledger: print a new entry to fill in')
    parser.add_argument('--add', type=Path, help='ledger: JSON file of one entry to append')
    parser.add_argument('--reference-width', type=int, default=512, help='Earth reference width (geology)')
    args, extra = parser.parse_known_args()
    if args.command not in COMMANDS[args.module]:
        parser.error(f'{args.module} commands: {", ".join(COMMANDS[args.module])}')
    if args.resolution is not None and (args.resolution < 2 or args.resolution % 2):
        parser.error('--resolution must be an even width')
    if args.output:
        runs.refuse_archive(args.output)
    if args.command == 'ledger':
        ledger_command(args, extra)
        return
    if args.command in ('experiment', 'evaluate') or (args.command == 'run' and args.module in ('climate', 'tectonics', 'world')):
        evaluation(args, extra)
        return
    if args.command == 'review':
        if not args.evaluation or len(args.evaluation) != 1:
            parser.error('review requires --evaluation DIR')
        call([UV, 'run', '--offline', EVALUATION / 'review.py', args.evaluation[0], *extra])
        return
    if args.command == 'dashboard' or args.evaluation:
        if args.command not in ('report', 'dashboard'):
            parser.error('--evaluation belongs to report')
        pages(args, extra)
        return
    if args.module in ('world', 'climate') and args.command == 'report':
        parser.error(f'{args.module} report requires --evaluation DIR [DIR ...]')
    {'tectonics': tectonics, 'geology': geology}[args.module](args, extra)


if __name__ == '__main__':
    main()
