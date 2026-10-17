"""Agent worktrees in .claude/worktrees: audit, preserve to storage, retire.

    uv run python scripts/uwtools/worktrees.py audit [--json FILE]
    uv run python scripts/uwtools/worktrees.py preserve NAME... | --all
    uv run python scripts/uwtools/worktrees.py retire NAME... --yes

audit (read-only): per worktree its links (never followed), whether its commits or uncommitted diff reached
main (the share of added lines found in some version of the file in main's history since the worktree's
base), its untracked files and sizes. preserve (copies only): writes <storage root>/worktree-archive/<name>/
with state.json (head, branch, status), uncommitted.patch (git diff --binary HEAD) and every untracked or
ignored file except .venv/, out/build/ and caches, each copy SHA-256-verified; a worktree holding a link is
refused. retire: deletes a preserved worktree's directory with storage.remove_tree (no delete when any link
is inside), only if nothing changed since the archive was written, then `git worktree prune`; the branch is
kept, so its commits stay reachable. Never `git worktree remove --force` (2026-10-07: it followed a junction
into refs/).
"""
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import storage  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402
from uwtools.paths import ROOT  # noqa: E402

WORKTREES = ROOT / '.claude' / 'worktrees'
# Regenerable: not archived, and changes inside them do not block a retire.
SKIPPED = ('.venv/', 'out/build/', '__pycache__/', '.pytest_cache/', '.ruff_cache/')


def git(*args, cwd=ROOT, binary=False):
    result = subprocess.run(['git', *args], cwd=cwd, capture_output=True)
    if binary:
        return result.returncode, result.stdout
    return result.returncode, result.stdout.decode('utf-8', 'replace').strip()


def skipped(relative):
    return any(relative.startswith(prefix) or f'/{prefix}' in relative for prefix in SKIPPED)


def archive_root():
    root = storage.storage_root()
    if root is None:
        raise OSError('No storage configured (runs/storage.json)')
    return root / 'worktree-archive'


def landed(worktree, head):
    """(added lines, found in main's history, the first missing ones, description of the change)."""
    in_main = git('merge-base', '--is-ancestor', head, 'main')[0] == 0
    base = head if in_main else git('merge-base', head, 'main')[1]
    diff = '' if in_main else git('diff', base, head)[1]
    describe = 'uncommitted diff' if in_main else 'commits ' + '; '.join(git('log', '--format=%h %s', f'{base}..{head}')[1].splitlines())
    uncommitted = git('diff', 'HEAD', cwd=worktree)[1]
    if uncommitted and not in_main:
        describe += ' + uncommitted diff'
    diff += '\n' + uncommitted
    added, path = {}, None
    for line in diff.splitlines():
        if line.startswith('+++ '):
            path = line[6:] if line.startswith('+++ b/') else None
        elif path and line.startswith('+') and len(line[1:].strip()) > 3:
            added.setdefault(path, []).append(line[1:].strip())
    total = found = 0
    missing = []
    for path, lines in added.items():
        versions = set()
        for commit in git('rev-list', f'{base}..main', '--', path)[1].split() + ['main']:
            versions.update(text.strip() for text in git('show', f'{commit}:{path}')[1].splitlines())
        for text in lines:
            total += 1
            if text in versions:
                found += 1
            else:
                missing.append(f'{path}: {text[:100]}')
    return total, found, missing, describe


def sizes(worktree):
    """Bytes per top-level entry (per out/<child> under out), links not followed."""
    totals = {}
    stack = [str(worktree)]
    while stack:
        try:
            entries = list(os.scandir(stack.pop()))
        except OSError:
            continue
        for entry in entries:
            if storage.is_link(entry):
                continue
            if entry.is_dir(follow_symlinks=False):
                stack.append(entry.path)
                continue
            parts = Path(entry.path).relative_to(worktree).parts
            key = '/'.join(parts[:2]) if parts[0] == 'out' and len(parts) > 2 else parts[0]
            totals[key] = totals.get(key, 0) + entry.stat(follow_symlinks=False).st_size
    return totals


def audit_one(worktree):
    head = git('rev-parse', 'HEAD', cwd=worktree)[1]
    total, found, missing, describe = landed(worktree, head)
    status = git('status', '--porcelain', '--untracked-files=all', cwd=worktree)[1].splitlines()
    totals = sizes(worktree)
    receipt = archive_root() / worktree.name / 'receipt.json'
    return {
        'name': worktree.name, 'head': head, 'branch': git('rev-parse', '--abbrev-ref', 'HEAD', cwd=worktree)[1],
        'date': git('log', '-1', '--format=%ci', 'HEAD', cwd=worktree)[1], 'change': describe,
        'added_lines': total, 'found_in_main': found, 'missing': missing[:20],
        'modified': [line[3:] for line in status if not line.startswith('??')],
        'untracked': [line[3:] for line in status if line.startswith('??')],
        'links': storage.find_links(worktree),
        'gb': round(sum(totals.values()) / 1e9, 2),
        'gb_by_part': {key: round(value / 1e9, 2) for key, value in sorted(totals.items(), key=lambda kv: -kv[1]) if value > 2e7},
        'archive': json.loads(receipt.read_text(encoding='utf-8'))['status'] if receipt.is_file() else None,
    }


def preserve(worktree):
    links = storage.find_links(worktree)
    if links:
        raise OSError(f'{worktree.name}: links inside, preserve and retire by hand after review: {links[:3]}')
    destination = archive_root() / worktree.name
    receipt_path = destination / 'receipt.json'
    if receipt_path.is_file():
        raise FileExistsError(f'{worktree.name}: already preserved ({receipt_path})')
    destination.mkdir(parents=True, exist_ok=True)
    state = audit_one(worktree)
    (destination / 'state.json').write_text(json.dumps(state, indent=1) + '\n', encoding='utf-8')
    patch = git('diff', '--binary', 'HEAD', cwd=worktree, binary=True)[1]
    (destination / 'uncommitted.patch').write_bytes(patch)
    listed = git('ls-files', '--others', '-z', cwd=worktree, binary=True)[1].decode('utf-8', 'replace').split('\0')
    digests = {}
    for relative in sorted(item for item in listed if item and not skipped(item)):
        source = worktree / relative
        if not source.is_file() or storage.is_link(source):
            continue
        target = destination / 'files' / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        digests[relative] = sha256(source)
        if sha256(target) != digests[relative]:
            raise OSError(f'Copy verification failed: {source}')
    receipt = {'worktree': str(worktree), 'head': state['head'], 'branch': state['branch'], 'status': 'preserved',
               'preserved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
               'patch_sha256': sha256(destination / 'uncommitted.patch'), 'sha256': digests}
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    return len(digests), sum((worktree / item).stat().st_size for item in digests) / 1e9


def retire(worktree):
    receipt_path = archive_root() / worktree.name / 'receipt.json'
    if not receipt_path.is_file():
        raise OSError(f'{worktree.name}: not preserved; run preserve first')
    receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
    if receipt['status'] != 'preserved':
        raise OSError(f'{worktree.name}: archive status {receipt["status"]}')
    stamp = receipt_path.stat().st_mtime
    if git('rev-parse', 'HEAD', cwd=worktree)[1] != receipt['head']:
        raise OSError(f'{worktree.name}: HEAD moved since the archive')
    if sha256_bytes(git('diff', '--binary', 'HEAD', cwd=worktree, binary=True)[1]) != receipt['patch_sha256']:
        raise OSError(f'{worktree.name}: the uncommitted diff changed since the archive')
    newer = [item for item in worktree.rglob('*') if item.is_file() and item.stat().st_mtime > stamp
             and not skipped(item.relative_to(worktree).as_posix()) and item.name != '.git']
    if newer:
        raise OSError(f'{worktree.name}: {len(newer)} files newer than the archive, e.g. {newer[0]}')
    storage.remove_tree(worktree)
    git('worktree', 'prune')
    receipt['status'] = 'retired'
    receipt['retired_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')


def sha256_bytes(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest='command', required=True)
    audit = commands.add_parser('audit', help='read-only report of every worktree')
    audit.add_argument('--json', help='also write the full report here')
    for name in ('preserve', 'retire'):
        command = commands.add_parser(name)
        command.add_argument('names', nargs='*')
        command.add_argument('--all', action='store_true')
        if name == 'retire':
            command.add_argument('--yes', action='store_true', required=True, help='confirm the delete')
    args = parser.parse_args()
    every = sorted(item for item in WORKTREES.iterdir() if item.is_dir()) if WORKTREES.is_dir() else []
    if args.command == 'audit':
        report = []
        for worktree in every:
            item = audit_one(worktree)
            report.append(item)
            print(f"{item['name']}  {item['date'][:16]}  {item['gb']:5.2f} GB  links {len(item['links'])}  "
                  f"landed {item['found_in_main']}/{item['added_lines']}  modified {len(item['modified'])}  "
                  f"untracked {len(item['untracked'])}  archive {item['archive']}  {item['change'][:80]}")
        print(f'total {sum(item["gb"] for item in report):.1f} GB in {len(report)} worktrees')
        if args.json:
            Path(args.json).write_text(json.dumps(report, indent=1) + '\n', encoding='utf-8')
        return
    chosen = every if args.all else [WORKTREES / name for name in args.names]
    for worktree in chosen:
        if args.command == 'preserve':
            files, gb = preserve(worktree)
            print(f'{worktree.name}: preserved {files} files ({gb:.2f} GB) and the uncommitted patch')
        else:
            retire(worktree)
            print(f'{worktree.name}: retired (branch kept)')


if __name__ == '__main__':
    main()
