"""Check the machine-readable records: the experiment ledger, the master plan's status rows and agents' reports.

    python scripts/check_records.py                              # the working tree's ledgers and plan
    python scripts/check_records.py --staged                     # the staged ledgers and plan (pre-commit hook)
    python scripts/check_records.py --report FILE --kind worker  # a worker's or mechanic's report (or reviewer)

Ledger (runs/ledger/*.jsonl): every line parses and ids are unique; entries above the cutoffs of
uwtools.ledger.GRANDFATHERED pass ledger.problems, their duration starts from a commit, and a master-plan step
named in their phase is a row of the plan's section 3. Plan (docs/architecture/master-plan.md, section 4): every
status row but the grandfathered ones names a section-3 step once, a status 'done|part|stopped: ...', commit
hashes or '(this commit)', a duration with '<n> min|h' and 'from <commit>', and only ledger ids that exist.
Reports: the headings of docs/reference-guides/agent-workflow.md ('The report', 'Review') in order and non-empty.
Field lists: docs/reference-guides/agent-workflow.md. --staged checks only when a ledger or the plan is staged.
Exit status: 0 clean, 1 findings (path:line: step or id: problem), 2 cannot check. Standard library only.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools import ledger  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PLAN = 'docs/architecture/master-plan.md'
LEDGERS = 'runs/ledger'
# Status rows written before the record fields were fixed (master plan 2.8).
GRANDFATHERED_ROWS = {'1.1', '1.6', 'protection (new)', '1.2', '1.3', '1.4', '2.1', '3.1', '3.2 (diagnosis)', '5.1',
                      'safeguards (new)'}
STEP_ID = r'\d+\.\d+[a-z]?'
ROW_STEP = re.compile(rf'^({STEP_ID})( \([a-z0-9 ]+\))?$')
ROW_STATUS = re.compile(r'^(done|part|stopped): \S')
HASH = r'[0-9a-f]{7,40}'
ROW_COMMIT = re.compile(rf'^{HASH}(, {HASH})*$')
ROW_SPAN = re.compile(r'\d+ (min|h)')
ROW_FROM = re.compile(rf'\bfrom ({HASH})\b')
THIS_COMMIT = ('this commit', '(this commit)')
PHASE_STEP = re.compile(rf'^master plan ({STEP_ID})\b')
REPORTS = {
    # (heading, required)
    'worker': [('Files changed', True), ('Checks', True), ('Not done', True), ('Not run', True), ('Questions', True),
               ('Triage', False), ('Commit', False), ('Duration', True)],
    'reviewer': [('Verdict', True), ('Findings', True), ('Claims', True), ('Not exercised', True), ('Scope', True)],
}
VERDICTS = ('pass with notes', 'pass', 'rework')


def section(lines, start, end):
    """(line number, text) of the lines after the heading starting with `start` up to the one starting with `end`."""
    inside, found = False, []
    for number, line in enumerate(lines, 1):
        if line.startswith(start):
            inside = True
        elif inside and line.startswith(end):
            break
        elif inside:
            found.append((number, line))
    return found


def plan_steps(text):
    """The step ids of the plan's section 3 tables."""
    return {match.group(1) for _, line in section(text.splitlines(), '## 3.', '## 4.')
            if (match := re.match(rf'^\| ({STEP_ID}) \|', line))}


def clean(cell):
    """A cell's text without Markdown emphasis or code marks."""
    return cell.replace('`', '').replace('*', '').strip()


def ledger_id_pattern(modules):
    """Ledger ids ('<module>-NNNN') of the given modules (the ledger file stems)."""
    names = '|'.join(re.escape(module) for module in sorted(modules)) or '(?!)'
    return re.compile(rf'\b((?:{names})-\d{{4}})\b')


def status_rows(text):
    """(line number, cleaned cells) of the section 4 table rows, header and separator excluded."""
    rows = []
    for number, line in section(text.splitlines(), '## 4. Status', '## 5.'):
        if not line.startswith('|') or line.startswith('|---') or line.startswith('| Step |'):
            continue
        rows.append((number, [clean(cell) for cell in line.strip().strip('|').split('|')]))
    return rows


def row_hashes(cells):
    if len(cells) != 4:
        return []
    commits = ROW_COMMIT.fullmatch(cells[2]) and re.findall(HASH, cells[2]) or []
    return commits + ROW_FROM.findall(cells[3])


def check_plan(text, known, ledger_ids, modules, path=PLAN):
    """Findings (path, line, step, problem) of the status rows; modules are the ledger file stems."""
    found, steps, seen, ledger_id = [], plan_steps(text), set(), ledger_id_pattern(modules)
    for number, cells in status_rows(text):
        if len(cells) != 4:
            found.append((path, number, cells[0] if cells else '?', f'row has {len(cells)} cells, not 4'))
            continue
        step, status, commit, duration = cells
        if step in GRANDFATHERED_ROWS:
            continue
        add = lambda problem: found.append((path, number, step, problem))  # noqa: E731
        match = ROW_STEP.fullmatch(step)
        if not match:
            add("step is not '<phase>.<step>' with an optional '(words)'")
        elif match.group(1) not in steps:
            add(f'step {match.group(1)} is not a row of section 3')
        if step in seen:
            add('step has two status rows')
        seen.add(step)
        if not ROW_STATUS.match(status):
            add("status does not start with 'done: ', 'part: ' or 'stopped: '")
        if commit not in THIS_COMMIT:
            if not ROW_COMMIT.fullmatch(commit):
                add("commit is not 'hash, hash...' or '(this commit)'")
            else:
                for value in re.findall(HASH, commit):
                    if not known(value):
                        add(f'commit {value} is not a commit')
        if not ROW_SPAN.search(duration):
            add("duration names no '<n> min' or '<n> h'")
        start = ROW_FROM.search(duration)
        if not start:
            add("duration names no 'from <commit>'")
        elif not known(start.group(1)):
            add(f'duration starts from {start.group(1)}, which is not a commit')
        for identifier in ledger_id.findall(' | '.join(cells)):
            if identifier not in ledger_ids:
                add(f'ledger id {identifier} is not in the ledger')
    return found


def parse_ledger(text, path):
    """(entries with their line numbers, findings) of one ledger file's text."""
    entries, found = [], []
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        try:
            entry = json.loads(line)
        except json.JSONDecodeError as error:
            found.append((path, number, '?', f'line does not parse: {error}'))
            continue
        if not isinstance(entry, dict):
            found.append((path, number, '?', 'line is not a JSON object'))
            continue
        entries.append((number, entry))
    return entries, found


def ledger_hashes(files):
    """The start commits of the strict entries of {path: text}."""
    hashes = []
    for path, text in files.items():
        for _, entry in parse_ledger(text, path)[0]:
            if ledger.strict(entry) and (value := ledger.start_commit(entry)):
                hashes.append(value)
    return hashes


def check_ledgers(files, known, steps):
    """Findings of the ledgers {path: text}, and the set of ids they hold."""
    found, identifiers = [], set()
    for path, text in sorted(files.items()):
        entries, problems = parse_ledger(text, path)
        found += problems
        earlier = []
        for number, entry in entries:
            identifier = str(entry.get('id', '?'))
            if identifier in identifiers:
                found.append((path, number, identifier, 'id is already used by an earlier entry'))
            identifiers.add(identifier)
            if ledger.strict(entry):
                others = [item for item in earlier if item.get('id') != entry.get('id')]
                found += [(path, number, identifier, problem) for problem in ledger.problems(entry, others)]
                start = ledger.start_commit(entry)
                if start and not known(start):
                    found.append((path, number, identifier, f'duration starts from {start}, which is not a commit'))
                phase = entry.get('phase')
                match = PHASE_STEP.match(phase) if isinstance(phase, str) else None
                if match and match.group(1) not in steps:
                    found.append((path, number, identifier, f'master plan step {match.group(1)} is not in the plan'))
            earlier.append(entry)
    return found, identifiers


def check_report(text, kind):
    """Problems (line, problem) of a report against the headings of its kind."""
    expected = REPORTS[kind]
    names = [name for name, _ in expected]
    lines = text.splitlines()
    headings, fenced = [], False
    for number, line in enumerate(lines, 1):
        if line.lstrip().startswith(('```', '~~~')):
            fenced = not fenced
        elif not fenced and line.startswith('## '):
            headings.append((number, line[3:].strip()))
    found, positions = [], {}
    for index, (number, name) in enumerate(headings):
        if name not in names:
            continue
        if name in positions:
            found.append((number, f"heading '{name}' appears twice"))
            continue
        end = headings[index + 1][0] - 1 if index + 1 < len(headings) else len(lines)
        body = '\n'.join(lines[number:end]).strip()
        positions[name] = (number, body)
        if not body:
            found.append((number, f"'{name}' is empty (write 'none')"))
    for name, required in expected:
        if required and name not in positions:
            found.append((0, f"heading '## {name}' is missing"))
    order = [name for name in names if name in positions]
    present = sorted(order, key=lambda name: positions[name][0])
    if present != order:
        found.append((0, f"headings are out of order: {', '.join(present)}; expected {', '.join(order)}"))
    if kind == 'reviewer' and 'Verdict' in positions and positions['Verdict'][1]:
        number, body = positions['Verdict']
        first = clean(body.splitlines()[0]).lower()
        if not re.match(rf"^({'|'.join(VERDICTS)})\b", first):
            found.append((number, f"verdict '{first}' is not one of {', '.join(VERDICTS)}"))
    return found


def git(*arguments, input=None):
    return subprocess.run(['git', *arguments], cwd=ROOT, input=input, capture_output=True, text=True, check=True,
                          encoding='utf-8')


def resolver(hashes):
    """known(hash): whether each hash names a commit, resolved in one git cat-file call."""
    unique = sorted(set(hashes))
    commits = set()
    if unique:
        output = git('cat-file', '--batch-check', input=''.join(f'{value}^{{commit}}\n' for value in unique)).stdout
        for value, line in zip(unique, output.splitlines()):
            parts = line.split()
            if len(parts) == 3 and parts[1] == 'commit':
                commits.add(value)
    return lambda value: value in commits


def working_files():
    plan = (ROOT / PLAN).read_text(encoding='utf-8')
    ledgers = {f'{LEDGERS}/{path.name}': path.read_text(encoding='utf-8')
               for path in sorted((ROOT / LEDGERS).glob('*.jsonl'))}
    return plan, ledgers


def staged_files():
    """The staged plan and ledgers, or None when neither is staged."""
    staged = git('diff', '--cached', '--name-only', '--diff-filter=ACMR').stdout.splitlines()
    if not any(path == PLAN or re.fullmatch(rf'{LEDGERS}/[^/]+\.jsonl', path) for path in staged):
        return None
    plan = git('show', f':{PLAN}').stdout
    names = [path for path in git('ls-files', '--cached', '--', f'{LEDGERS}/*.jsonl').stdout.splitlines()
             if re.fullmatch(rf'{LEDGERS}/[^/]+\.jsonl', path)]
    return plan, {path: git('show', f':{path}').stdout for path in names}


def check_records(plan, ledgers):
    known = resolver(ledger_hashes(ledgers) + [value for _, cells in status_rows(plan) for value in row_hashes(cells)])
    found, identifiers = check_ledgers(ledgers, known, plan_steps(plan))
    return found + check_plan(plan, known, identifiers, {Path(path).stem for path in ledgers})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--staged', action='store_true', help='check the staged ledgers and plan')
    parser.add_argument('--report', type=Path, help="an agent's report to check")
    parser.add_argument('--kind', choices=sorted(REPORTS), default='worker', help='the report kind')
    arguments = parser.parse_args(argv)
    try:
        if arguments.report:
            problems = check_report(arguments.report.read_text(encoding='utf-8'), arguments.kind)
            for number, problem in problems:
                print(f'{arguments.report}:{number}: report: {problem}')
            return 1 if problems else 0
        files = staged_files() if arguments.staged else working_files()
        if files is None:
            return 0
        found = check_records(*files)
    except (OSError, subprocess.CalledProcessError, UnicodeDecodeError) as error:
        print(f'check_records: cannot check: {error}', file=sys.stderr)
        return 2
    for path, number, label, problem in found:
        print(f'{path}:{number}: {label}: {problem}')
    return 1 if found else 0


if __name__ == '__main__':
    sys.exit(main())
