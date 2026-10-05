"""Experiment ledger: runs/ledger/<module>.jsonl, one JSON line per experiment, append-only.

An entry records what was tried and what came of it (docs/architecture/consolidation-plan.md,
section 3). LEDGER.md beside the ledgers is generated from them and starts with what is
settled and what must not be retried.
"""
import datetime
import json
import re

from .paths import ROOT

DIRECTORY = ROOT / 'runs/ledger'
SUMMARY = DIRECTORY / 'LEDGER.md'
REQUIRED = ('id', 'date', 'module', 'parent', 'hypothesis', 'change', 'expected', 'observed', 'assessment',
            'verdict', 'decision', 'lessons', 'do_not_retry')
DECISIONS = ('keep', 'revert', 'park')
FILLED = ('hypothesis', 'expected', 'assessment', 'verdict')


def path(module):
    return DIRECTORY / f'{module}.jsonl'


def modules():
    return sorted(file.stem for file in DIRECTORY.glob('*.jsonl'))


def entries(module):
    file = path(module)
    if not file.is_file():
        return []
    return [json.loads(line) for line in file.read_text(encoding='utf-8').splitlines() if line.strip()]


def next_id(module):
    numbers = [int(match.group(1)) for entry in entries(module) if (match := re.fullmatch(rf'{module}-(\d+)', entry['id']))]
    return f'{module}-{max(numbers, default=0) + 1:04d}'


def validate(entry, existing):
    missing = [key for key in REQUIRED if key not in entry]
    if missing:
        raise ValueError(f"ledger entry lacks {', '.join(missing)}")
    identifiers = {item['id'] for item in existing}
    if entry['id'] in identifiers:
        raise ValueError(f"{entry['id']} is already in the ledger; entries are never rewritten")
    if not re.fullmatch(rf"{entry['module']}-\d{{4}}", entry['id']):
        raise ValueError(f"id {entry['id']!r} is not {entry['module']}-NNNN")
    if entry['parent'] is not None and entry['parent'] not in identifiers:
        raise ValueError(f"parent {entry['parent']} is not in the {entry['module']} ledger")
    # What was tried and what came of it must be filled in: a draft's empty placeholders are refused.
    empty = [key for key in FILLED if not (entry[key].strip() if isinstance(entry[key], str) else entry[key])]
    if empty:
        raise ValueError(f"ledger entry leaves {', '.join(empty)} empty")
    if entry['decision'] not in DECISIONS:
        raise ValueError(f"decision must be one of {', '.join(DECISIONS)}")
    for key in ('lessons', 'do_not_retry'):
        if not isinstance(entry[key], list):
            raise ValueError(f'{key} must be a list')


def template(module, **values):
    """An entry with every required key, for the caller to fill in."""
    existing = entries(module)
    entry = {'id': next_id(module), 'date': datetime.date.today().isoformat(), 'module': module,
             'parent': existing[-1]['id'] if existing else None, 'hypothesis': '', 'change': {}, 'expected': '',
             'observed': {}, 'not_run': [], 'assessment': '', 'verdict': '', 'decision': 'keep', 'lessons': [],
             'do_not_retry': []}
    entry.update(values)
    return entry


def append(entry):
    """Validate and append one entry; rewrite LEDGER.md."""
    existing = entries(entry['module'])
    validate(entry, existing)
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    with path(entry['module']).open('a', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(entry, ensure_ascii=False) + '\n')
    write_summary()
    return entry


def matches(entry, topic):
    return topic.lower() in json.dumps(entry, ensure_ascii=False).lower()


def relevant(module, topic=None, last=5):
    """Entries naming the topic anywhere, or the last few without a topic."""
    found = entries(module)
    return [entry for entry in found if matches(entry, topic)] if topic else found[-last:]


def brief(entry):
    lines = [f"{entry['id']} ({entry['date']}, {entry['verdict']}, {entry['decision']}): {entry.get('phase') or entry['hypothesis']}"]
    lines += [f'  lesson: {text}' for text in entry['lessons']]
    lines += [f'  DO NOT RETRY: {text}' for text in entry['do_not_retry']]
    return '\n'.join(lines)


def summary():
    lines = ['# Experiment ledger', '',
             'Generated from `runs/ledger/*.jsonl` by `scripts/uwtools/ledger.py`; do not edit. '
             'Show entries with `uv run scripts/uw.py <module> ledger [--topic TEXT]`.', '', "## Settled / don't retry", '']
    blocked = [(entry['id'], text) for module in modules() for entry in entries(module) for text in entry['do_not_retry']]
    lines += [f'- `{identifier}`: {text}' for identifier, text in blocked] or ['Nothing recorded.']
    for module in modules():
        lines += ['', f'## {module}', '', '| Id | Date | Verdict | Decision | What |', '|---|---|---|---|---|']
        for entry in entries(module):
            what = (entry.get('phase') or entry['hypothesis']).replace('|', '/')
            lines.append(f"| `{entry['id']}` | {entry['date']} | {entry['verdict']} | {entry['decision']} | {what} |")
        lessons = [(entry['id'], text) for entry in entries(module) for text in entry['lessons']]
        if lessons:
            lines += ['', 'Lessons:', '', *[f'- `{identifier}`: {text}' for identifier, text in lessons]]
    return '\n'.join(lines) + '\n'


def write_summary():
    SUMMARY.write_text(summary(), encoding='utf-8', newline='\n')
