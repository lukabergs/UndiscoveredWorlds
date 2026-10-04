"""Record and show review verdicts: the user's judgement of an evaluated run.

    uv run --offline python scripts/evaluation/review.py <evaluation-dir> --verdict approved --notes "..."
    uv run --offline python scripts/evaluation/review.py <evaluation-dir>

review-request.json in the evaluation directory holds the module's checklist (the "review"
entries of the profile) with the maps of each entry and the metric table; the verifier agent
(.claude/agents/verifier.md) answers it as --reviewer agent with --answers FILE.json. The
user's verdict overrides an agent's. Verdicts are appended to <evaluation-dir>/review.json with the evaluated bundle's
manifest hash, so a verdict stays attached to the exact fields that were reviewed.
The run page (index.html, report.md) is rendered again with the new entry.
"""
import argparse
import datetime
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools.bundle import SEASONS

VERDICTS = ('approved', 'rejected', 'revise')


def load(directory):
    path = Path(directory) / 'review.json'
    if not path.is_file():
        return {'schema': 'uw-review/v1', 'entries': []}
    return json.loads(path.read_text(encoding='utf-8'))


def latest(directory):
    """The verdict that counts: the user's latest entry, else the latest of any reviewer."""
    entries = load(directory)['entries']
    own = [entry for entry in entries if entry['reviewer'] == 'user']
    return (own or entries)[-1] if entries else None


def blind(description):
    """A bundle description without the command and recipe, which would tell a reviewer what was changed."""
    return description and {key: description[key] for key in ('seed', 'grid', 'manifest_sha256') if key in description}


def request(profile, document):
    """What a reviewer needs: the profile's checklist with the maps of each entry, and the metric table."""
    checklist = []
    for entry in profile.get('review', []):
        names = set(entry.get('maps', []))
        maps = [item for item in document['maps'] if item.get('map') in names]
        if entry.get('seasons', 'all') != 'all':
            maps = [item for item in maps if item.get('season') in (None, *entry['seasons'])]
        checklist.append(dict(entry, images=[{key: item.get(key) for key in ('title', 'path', 'kind', 'season', 'palette', 'half_range')} for item in maps]))
    metrics = []
    for field, seasons in document['validation']['maps'].items():
        for season, entries in seasons.items():
            loss, base = entries['run']['loss'], entries.get('baseline', {}).get('loss')
            metrics.append({'field': field, 'season': 'annual' if len(seasons) == 1 else SEASONS[int(season)], 'loss': loss, 'parent_loss': base,
                            'delta': None if base is None else loss - base,
                            'metrics': {name: value for name, value, _ in entries['run']['rows']}})
    return {'schema': 'uw-review-request/v1', 'module': profile['module'], 'run': blind(document['run']), 'parent': blind(document['baseline']),
            'palettes': 'configs/climate-palettes.json and configs/tectonics-palettes.json hold the colour scales; a difference map spans -half_range to +half_range',
            'checks_not_passed': [item for item in document['checks'] if item['status'] in ('fail', 'reported')],
            'scores': document['validation']['scores'], 'metrics': metrics, 'checklist': checklist}


def record(directory, verdict, notes, reviewer, against, answers=None):
    directory = Path(directory)
    evaluation = json.loads((directory / 'evaluation.json').read_text(encoding='utf-8'))
    document = load(directory)
    document['entries'].append({
        'recorded': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'verdict': verdict, 'reviewer': reviewer, 'notes': notes, 'against': against, 'answers': answers or [],
        'bundle_manifest_sha256': evaluation['run']['manifest_sha256']})
    (directory / 'review.json').write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    import report  # report imports this module
    report.run_page(directory)
    return document['entries'][-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('evaluation', type=Path, help='Directory written by evaluate.py')
    parser.add_argument('--verdict', choices=VERDICTS)
    parser.add_argument('--notes', default='')
    parser.add_argument('--reviewer', default='user')
    parser.add_argument('--against', help='Run or evaluation the judgement compares with')
    parser.add_argument('--answers', type=Path, help='JSON list of {id, answer, evidence} per checklist entry')
    args = parser.parse_args()
    if args.verdict:
        answers = json.loads(args.answers.read_text(encoding='utf-8')) if args.answers else None
        print(json.dumps(record(args.evaluation, args.verdict, args.notes, args.reviewer, args.against, answers), indent=2))
    else:
        print(json.dumps(load(args.evaluation), indent=2))


if __name__ == '__main__':
    main()
