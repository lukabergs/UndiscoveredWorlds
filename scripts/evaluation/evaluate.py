"""Evaluate a run bundle with a module profile: checks, validation and review maps.

    uv run --offline python scripts/evaluation/evaluate.py \\
        --profile climate --bundle <run>/bundle --output <run>/evaluation \\
        [--baseline <bundle>] [--repeat <bundle>] [--input <bundle>] [--no-maps]

Profiles are configs/evaluation/<module>.json. Writes evaluation.json, report.md,
index.html and maps/ to a new or empty output directory; evaluations are evidence
and are never overwritten. Exit code 1 when a check fails.
"""
import argparse
import datetime
import html
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import SEASONS, Bundle
from uwtools.hashing import sha256_bytes
from checks import run_checks
from references import References
from validation import run_validation
import review

ROOT = Path(__file__).resolve().parents[2]
PROFILE_DIRECTORY = ROOT / 'configs/evaluation'
PROFILE_SCHEMA = 'uw-evaluation-profile/v1'
SCHEMA = 'uw-evaluation/v1'


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
    diff = git('diff', 'HEAD', '--binary')
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


def fmt(value, digits=4):
    if value is None:
        return '–'
    if isinstance(value, float):
        return f'{value:.{digits}g}'
    return str(value)


def validation_rows(validation):
    """(field, season label, run loss, baseline loss or None, key metrics) per map comparison."""
    rows = []
    for field, seasons in validation['maps'].items():
        for season, entries in seasons.items():
            label = 'annual' if len(seasons) == 1 else SEASONS[int(season)]
            metrics = {name: value for name, value, _ in entries['run']['rows']}
            rows.append((field, label, entries['run']['loss'], entries.get('baseline', {}).get('loss'), metrics))
    return rows


def write_markdown(path, document):
    lines = [f"# Evaluation · {document['profile']['module']}", '',
             f"- Run bundle: `{document['run']['bundle']}` (seed {document['run']['seed']}, recipe "
             f"{document['run']['recipe']}, {document['run']['grid'][0]} x {document['run']['grid'][1]})"]
    if document['baseline']:
        lines.append(f"- Baseline bundle: `{document['baseline']['bundle']}`")
    v = document['verdict']
    lines += [f"- Checks: {v['checks']}" + (f" ({', '.join(v['failed_checks'])})" if v['failed_checks'] else ''),
              f"- Validation: {v['validation']}", f"- Review: {v['review']}", '',
              f"Evidence weighting: {json.dumps(document['evidence'])}", '', '## Checks', '',
              '| Check | Status | Summary |', '|---|---|---|']
    lines += [f"| {item['check']} | {item['status']} | {item['summary']} |" for item in document['checks']]
    for item in document['checks']:
        if item['details'] and item['status'] == 'fail':
            lines += ['', f"{item['check']}:", *[f'- {detail}' for detail in item['details'][:20]]]
    validation = document['validation']
    if validation['maps'] or validation['distributions']:
        lines += ['', '## Validation', '']
        for key, entry in validation['scores'].items():
            value = entry.get('score')
            if value is not None:
                lines.append(f"- {key}: {entry['version']} on {entry['inputs']} = {value:.4f} (lower is better). {entry['note']}")
            else:
                lines.append(f"- {key}: partial score {fmt(entry['partial_score'])} over weight "
                             f"{entry['covered_weight']:.2f}; {entry['note']}")
        if validation['maps']:
            lines += ['', '| Field | Season | Loss | Baseline loss | RMSE | Bias | Correlation |', '|---|---|---|---|---|---|---|']
        for field, label, loss, base, metrics in validation_rows(validation):
            rmse = metrics.get('RMSE', metrics.get('Vector RMSE'))
            bias = metrics.get('Bias', metrics.get('Speed bias'))
            correlation = metrics.get('Spatial correlation', metrics.get('Speed spatial correlation'))
            lines.append(f'| {field} | {label} | {fmt(loss)} | {fmt(base)} | {fmt(rmse)} | {fmt(bias)} | {fmt(correlation)} |')
        for name, entry in validation['distributions'].items():
            lines += ['', f'### {name} ({entry["reference"]})', '']
            for key in ('run', 'baseline'):
                if key in entry:
                    lines.append(f"- {key}: " + ', '.join(f'{metric} {fmt(value)}' for metric, value in entry[key]['earth_comparison'].items()))
        for item in validation['skipped']:
            lines.append(f"- skipped {item['id']}: {item['reason']}")
    if document['maps']:
        lines += ['', '## Maps', '', *[f"- [{item['title']}]({item['path']})" for item in document['maps']]]
    lines += ['', '## Provenance', '', f"- Profile: `{document['profile']['path']}` ({document['profile']['sha256'][:16]})",
              f"- Code: {document['code']['revision']} (dirty: {document['code']['dirty']})",
              *[f'- Reference `{name}`: {digest[:16]}' for name, digest in document['references'].items()]]
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def write_html(path, document):
    e = html.escape
    v = document['verdict']
    checks = ''.join(f"<tr class='{item['status']}'><td>{e(item['check'])}</td><td>{item['status']}</td><td>{e(item['summary'])}</td></tr>"
                     for item in document['checks'])
    rows = ''.join(f'<tr><td>{e(field)}</td><td>{label}</td><td>{fmt(loss)}</td><td>{fmt(base)}</td></tr>'
                   for field, label, loss, base, _ in validation_rows(document['validation']))
    scores = ''.join(f"<li>{e(key)}: {e(entry['version'])} {fmt(entry.get('score'))} " +
                     (f"(partial {fmt(entry['partial_score'])}, weight {entry['covered_weight']:.2f})" if 'partial_score' in entry else '') + '</li>'
                     for key, entry in document['validation']['scores'].items())
    maps = ''.join(f"<figure><img src='{e(item['path'])}' alt='{e(item['title'])}'><figcaption>{e(item['title'])}</figcaption></figure>"
                   for item in document['maps'])
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Evaluation {e(document['profile']['module'])}</title><style>
:root{{--bg:#fbfaf7;--fg:#1f2328;--muted:#5d6670;--line:#d9d6cf;--fail:#b3261e;--pass:#1e6b3a}}
@media (prefers-color-scheme: dark){{:root{{--bg:#16181c;--fg:#e8e6e1;--muted:#a3a9b1;--line:#363a40;--fail:#ff8a80;--pass:#7fd49b}}}}
body{{margin:0 auto;max-width:1200px;padding:16px;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif}}
table{{border-collapse:collapse;width:100%;margin:8px 0 24px}}td,th{{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left}}
tr.fail td:nth-child(2){{color:var(--fail);font-weight:600}}tr.pass td:nth-child(2){{color:var(--pass)}}
.maps{{display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:12px}}
figure{{margin:0}}img{{width:100%;image-rendering:pixelated;border:1px solid var(--line)}}figcaption,.muted{{color:var(--muted)}}
</style></head><body>
<h1>Evaluation · {e(document['profile']['module'])}</h1>
<p class="muted">{e(document['run']['bundle'])} · seed {document['run']['seed']} · {e(document['run']['recipe'])} · {document['run']['grid'][0]} × {document['run']['grid'][1]}</p>
<p>Checks: <b>{v['checks']}</b> · Validation: <b>{v['validation']}</b> · Review: <b>{v['review']}</b></p>
<h2>Checks</h2><table><tr><th>Check</th><th>Status</th><th>Summary</th></tr>{checks}</table>
<h2>Validation</h2><ul>{scores}</ul><table><tr><th>Field</th><th>Season</th><th>Loss</th><th>Baseline loss</th></tr>{rows}</table>
<h2>Maps</h2><div class="maps">{maps}</div></body></html>"""
    path.write_text(page, encoding='utf-8')


def evaluate(profile, profile_path, bundle, output, baseline=None, repeat=None, input_bundle=None, maps=True):
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise SystemExit(f'{output} is not empty; evaluations are never overwritten')
    context = {'input': input_bundle, 'repeat': repeat}
    checks = run_checks(bundle, profile.get('checks', []), context)
    references = References()
    validation = run_validation(bundle, profile.get('validation', []), references, baseline)
    rendered = []
    if maps and profile.get('maps'):
        from maps import render
        rendered = render(bundle, profile['maps'], output)
    document = {
        'schema': SCHEMA,
        'evaluated': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'profile': {'module': profile['module'], 'path': profile_path.relative_to(ROOT).as_posix() if profile_path.is_relative_to(ROOT) else str(profile_path),
                    'sha256': sha256_bytes(profile_path.read_bytes())},
        'evidence': profile['evidence'], 'run': describe(bundle), 'baseline': describe(baseline),
        'repeat': describe(repeat), 'input': describe(input_bundle), 'code': code_state(),
        'checks': checks, 'validation': validation, 'references': references.provenance(), 'maps': rendered}
    document['verdict'] = verdict(profile, checks, validation, review.latest(output))
    (output / 'evaluation.json').write_text(json.dumps(document, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    write_markdown(output / 'report.md', document)
    write_html(output / 'index.html', document)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--profile', required=True, help='Module name in configs/evaluation/ or a profile path')
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--baseline', type=Path, help='Bundle compared with the run on common masks')
    parser.add_argument('--repeat', type=Path, help='Second run of the same seed and code (determinism check)')
    parser.add_argument('--input', type=Path, help='Bundle holding the inputs of frozen components')
    parser.add_argument('--no-maps', action='store_true')
    args = parser.parse_args()
    profile, path = load_profile(args.profile)
    open_bundle = lambda p: Bundle.open(p) if p else None
    document = evaluate(profile, path, Bundle.open(args.bundle), args.output.absolute(), open_bundle(args.baseline),
                        open_bundle(args.repeat), open_bundle(args.input), not args.no_maps)
    print(json.dumps(document['verdict'], indent=2))
    for key, entry in document['validation']['scores'].items():
        if 'score' in entry:
            print(f"{key}: {entry['version']} on {entry['inputs']} = {entry['score']:.4f}")
        else:
            print(f"{key}: partial {entry['version']} = {fmt(entry['partial_score'])} over weight {entry['covered_weight']:.2f}")
    print(f"Report: {args.output / 'report.md'}")
    sys.exit(1 if document['verdict']['checks'] == 'fail' else 0)


if __name__ == '__main__':
    main()
