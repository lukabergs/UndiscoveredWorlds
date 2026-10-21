"""Evaluate a run bundle with a module profile: checks, validation and review maps.

    uv run --offline python scripts/evaluation/evaluate.py \\
        --profile climate --bundle <run>/bundle --output <run>/evaluation \\
        [--baseline <bundle>] [--repeat <bundle>] [--explicit <bundle>] [--input <bundle>] [--checkpoints DIR] [--no-maps]
        [--earth NAME] [--coarser <bundle>]

Profiles are configs/evaluation/<module>.json. Writes evaluation.json and maps/ to a new
or empty output directory, then the run page (index.html) and report.md from it with
report.py; evaluations are evidence and are never overwritten. Checkpoint bundles
(--checkpoints, default <bundle>/../checkpoints when present) give map animations. A profile with an
"earth" block draws an Earth reference window beside a regional map (earth_windows.py; --earth NAME
chooses the catalogue entry, otherwise the nearest-latitude rule). A regional refinement hierarchy
bundle carries its per-level diagnostics into evaluation.json ("levels", shown on the run page);
--coarser names a coarser hierarchy bundle of the same world for level-consistency and
resolution-invariance.
Exit code 1 when a check fails.
"""
import argparse
import datetime
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import SCHEMAS as BUNDLE_SCHEMAS, Bundle, read_manifest
from uwtools.hashing import sha256_bytes
from uwtools import runs
from checks import hierarchy_levels, run_checks
from references import References
from validation import run_validation
import features
import report
import review

ROOT = Path(__file__).resolve().parents[2]
PROFILE_DIRECTORY = ROOT / 'configs/evaluation'
PROFILE_SCHEMA = 'uw-evaluation-profile/v1'
SCHEMA = 'uw-evaluation/v1'
FAILED = 'evaluation-error.json'


def load_profile(name):
    path = Path(name)
    if not path.suffix:
        path = PROFILE_DIRECTORY / f'{name}.json'
    profile = json.loads(path.read_text(encoding='utf-8'))
    if profile.get('schema') != PROFILE_SCHEMA:
        raise ValueError(f'{path}: expected schema {PROFILE_SCHEMA}')
    return profile, path


def code_state():
    def git(*arguments):
        completed = subprocess.run(['git', *arguments], cwd=ROOT, capture_output=True)
        return completed.stdout if completed.returncode == 0 else None
    revision = git('rev-parse', 'HEAD')
    # The run registries change with every registered or named run; they are not code (naming.py).
    diff = git('diff', 'HEAD', '--binary', '--', '.', ':!runs/registry')
    return {'revision': revision.decode().strip() if revision else None,
            'dirty': bool(diff), 'diff_sha256': sha256_bytes(diff) if diff else None}


def describe(bundle):
    if bundle is None:
        return None
    manifest = (bundle.path / 'manifest.json').read_bytes()
    return {'bundle': bundle.path.as_posix(), 'manifest_sha256': sha256_bytes(manifest), 'seed': bundle.seed,
            'recipe': bundle.recipe, 'grid': [bundle.columns, bundle.rows], 'created': bundle.manifest.get('created'),
            'command': bundle.manifest.get('command', [])}


def score_value(entry):
    return entry.get('score', entry.get('partial_score'))


def verdict(profile, checks, validation, review_entry):
    failed = [item['check'] for item in checks if item['status'] == 'fail']
    result = {'checks': 'fail' if failed else 'pass', 'failed_checks': failed}
    scores = validation.get('scores', {})
    if 'baseline' in scores:
        run, base = score_value(scores['run']), score_value(scores['baseline'])
        if run is None or base is None:  # no weighted field was compared: nothing to rank
            result['validation'] = 'scored'
        else:
            result['validation'] = 'improved' if run < base - 1e-9 else 'worse' if run > base + 1e-9 else 'equal'
    elif scores:
        result['validation'] = 'scored'
    elif validation.get('distributions'):
        result['validation'] = 'reported'
    else:
        result['validation'] = 'none'
    required = profile['evidence'].get('review', 'optional')
    result['review'] = review_entry['verdict'] if review_entry else ('pending' if required != 'optional' else 'optional')
    return result


def checkpoint_bundles(directory):
    """Run bundles in the subdirectories of directory, in name order (the checkpoint order)."""
    if directory is None or not Path(directory).is_dir():
        return []
    return [Bundle.open(path) for path in sorted(Path(directory).iterdir())
            if (path / 'manifest.json').is_file() and read_manifest(path).get('schema') in BUNDLE_SCHEMAS]


def evaluate(profile, profile_path, bundle, output, baseline=None, repeat=None, input_bundle=None, maps=True,
             checkpoints=None, explicit=None, earth=None, coarser=None):
    output = runs.refuse_archive(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        failed = f' (it holds a failed evaluation: {FAILED})' if (output / FAILED).is_file() else ''
        raise SystemExit(f'{output} is not empty{failed}; evaluations are never overwritten: use another --output')
    try:
        return evaluate_into(profile, profile_path, bundle, output, baseline, repeat, input_bundle, maps, checkpoints, explicit,
                             earth, coarser)
    except BaseException as error:
        # What an interrupted or failed evaluation left behind stays as evidence, marked as incomplete.
        (output / FAILED).write_text(json.dumps({
            'schema': 'uw-evaluation-error/v1', 'error': f'{type(error).__name__}: {error}',
            'recorded': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
            'note': 'This evaluation did not complete; files beside this one are partial.'}, indent=2) + '\n', encoding='utf-8')
        raise


def evaluate_into(profile, profile_path, bundle, output, baseline, repeat, input_bundle, maps, checkpoints, explicit, earth=None,
                  coarser=None):
    context = {'input': input_bundle, 'repeat': repeat, 'explicit': explicit, 'coarser': coarser}
    checks = run_checks(bundle, profile.get('checks', []), context)
    references = References()
    validation = run_validation(bundle, profile.get('validation', []), references, baseline)
    if validation['maps']:
        validation['features'] = features.features(bundle, profile['validation'], references, baseline)
        if baseline is not None:
            validation['footprint'] = features.footprint(bundle, baseline, profile['validation'])
    rendered = []
    if maps and profile.get('maps'):
        from maps import render
        if checkpoints is None:
            checkpoints = bundle.path.parent / 'checkpoints'
        rendered = render(bundle, profile['maps'], output, profile.get('validation', []), validation['maps'], references,
                          checkpoint_bundles(checkpoints), baseline, validation['distributions'])
    earth_window = None
    if maps and profile.get('earth'):
        import earth_windows
        written, earth_window = earth_windows.render_maps(bundle, profile, output, earth)
        rendered += written
    document = {
        'schema': SCHEMA,
        'evaluated': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'profile': {'module': profile['module'], 'path': profile_path.relative_to(ROOT).as_posix() if profile_path.is_relative_to(ROOT) else str(profile_path),
                    'sha256': sha256_bytes(profile_path.read_bytes())},
        'evidence': profile['evidence'], 'run': describe(bundle), 'baseline': describe(baseline),
        'repeat': describe(repeat), 'explicit': describe(explicit), 'input': describe(input_bundle),
        'coarser': describe(coarser), 'code': code_state(),
        'checks': checks, 'validation': validation, 'references': references.provenance(), 'maps': rendered,
        'earth': earth_window, 'levels': hierarchy_levels(bundle)}
    document['verdict'] = verdict(profile, checks, validation, review.latest(output))
    (output / 'evaluation.json').write_text(json.dumps(document, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if profile.get('review'):
        (output / 'review-request.json').write_text(json.dumps(review.request(profile, document), indent=2) + '\n', encoding='utf-8')
    report.run_page(output)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--profile', required=True, help='Module name in configs/evaluation/ or a profile path')
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline', type=Path, help='Bundle compared with the run on common masks')
    parser.add_argument('--repeat', type=Path, help='Second run of the same seed and code (determinism check)')
    parser.add_argument('--explicit', type=Path, help='Run given the recipe file explicitly (recipe-equality check)')
    parser.add_argument('--input', type=Path, help='Bundle holding the inputs of frozen components')
    parser.add_argument('--coarser', type=Path, help='Coarser refinement hierarchy bundle of the same world '
                        '(regional level-consistency and resolution-invariance checks)')
    parser.add_argument('--checkpoints', type=Path, help='Directory of checkpoint bundles (default <bundle>/../checkpoints)')
    parser.add_argument('--no-maps', action='store_true')
    parser.add_argument('--earth', metavar='NAME', help='Earth reference window of the catalogue beside a regional map '
                        '(profiles with an "earth" block; default: the nearest-latitude entry holding the map\'s metres per pixel)')
    args = parser.parse_args()
    profile, path = load_profile(args.profile)
    if args.earth is not None:
        if not profile.get('earth'):
            raise SystemExit(f"--earth: the {profile['module']} profile draws no Earth windows")
        import earth_windows
        try:
            earth_windows.choose_entry(earth_windows.load_catalogue(ROOT / profile['earth'].get('catalogue', earth_windows.CATALOGUE)),
                                       0.0, 0.0, args.earth)
        except (ValueError, OSError) as error:
            raise SystemExit(f'--earth: {error}')
    open_bundle = lambda p: Bundle.open(p) if p else None
    document = evaluate(profile, path, Bundle.open(args.bundle), args.output.absolute(), open_bundle(args.baseline),
                        open_bundle(args.repeat), open_bundle(args.input), not args.no_maps, args.checkpoints,
                        open_bundle(args.explicit), args.earth, open_bundle(args.coarser))
    print(json.dumps(document['verdict'], indent=2))
    for key, entry in document['validation']['scores'].items():
        if 'score' in entry:
            print(f"{key}: {entry['version']} on {entry['inputs']} = {entry['score']:.4f}")
        else:
            print(f"{key}: partial {entry['version']} = {report.fmt(entry['partial_score'])} over weight {entry['covered_weight']:.2f}")
    print(f"Report: {args.output / 'report.md'}")
    sys.exit(1 if document['verdict']['checks'] == 'fail' else 0)


if __name__ == '__main__':
    main()
