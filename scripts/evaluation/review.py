# /// script
# dependencies = ["numpy", "pillow"]
# ///
"""Record and show review verdicts: the user's judgement of an evaluated run.

    uv run --offline python scripts/evaluation/review.py <evaluation-dir> --verdict approved --notes "..."
    uv run --offline python scripts/evaluation/review.py <evaluation-dir>

Verdicts are appended to <evaluation-dir>/review.json with the evaluated bundle's
manifest hash, so a verdict stays attached to the exact fields that were reviewed.
"""
import argparse
import datetime
import json
from pathlib import Path

VERDICTS = ('approved', 'rejected', 'revise')


def load(directory):
    path = Path(directory) / 'review.json'
    if not path.is_file():
        return {'schema': 'uw-review/v1', 'entries': []}
    return json.loads(path.read_text(encoding='utf-8'))


def latest(directory):
    entries = load(directory)['entries']
    return entries[-1] if entries else None


def record(directory, verdict, notes, reviewer, against):
    directory = Path(directory)
    evaluation = json.loads((directory / 'evaluation.json').read_text(encoding='utf-8'))
    document = load(directory)
    document['entries'].append({
        'recorded': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'verdict': verdict, 'reviewer': reviewer, 'notes': notes, 'against': against,
        'bundle_manifest_sha256': evaluation['run']['manifest_sha256']})
    (directory / 'review.json').write_text(json.dumps(document, indent=2) + '\n', encoding='utf-8')
    return document['entries'][-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('evaluation', type=Path, help='Directory written by evaluate.py')
    parser.add_argument('--verdict', choices=VERDICTS)
    parser.add_argument('--notes', default='')
    parser.add_argument('--reviewer', default='user')
    parser.add_argument('--against', help='Run or evaluation the judgement compares with')
    args = parser.parse_args()
    if args.verdict:
        print(json.dumps(record(args.evaluation, args.verdict, args.notes, args.reviewer, args.against), indent=2))
    else:
        print(json.dumps(load(args.evaluation), indent=2))


if __name__ == '__main__':
    main()
