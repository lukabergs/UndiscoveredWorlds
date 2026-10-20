"""The agents' trash: untracked data they no longer need is moved here, never deleted.

    uv run python scripts/uwtools/trash.py move PATH... [--reason TEXT]
    uv run python scripts/uwtools/trash.py list
    uv run python scripts/uwtools/trash.py restore ENTRY

The user's rule (2026-10-20): agents may get rid of files as long as code makes an accidental loss impossible.
A delete of untracked data is refused by delete_guard.py; this moves it instead to
<storage root>/trash/<entry>/ (the entry is a UTC time stamp and a number) with a receipt (original path, reason,
files, bytes, SHA-256 of every file). On the same drive the move is a rename, so nothing is copied or deleted; from
another drive every file is copied, verified by digest, and only then removed from its source. Only the user
empties the trash (guard_catastrophe.py protects the folder from agent commands).

Refuses: links (the path itself or anything inside it: a recursive move must never follow a junction,
2026-10-07), tracked files (use git rm), anything outside the repository and the storage root, and the protected
places: the repository root, .git, refs/, runs/registry, runs/ledger, runs/archive, the storage root, its
worktree-archive, receipts and trash, and the junction points themselves.
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


def trash_root():
    root = storage.storage_root()
    if root is None:
        raise OSError('No storage configured (runs/storage.json)')
    return root / 'trash'


def protected_places():
    storage_root = storage.storage_root()
    return [ROOT, ROOT / '.git', ROOT / 'refs', ROOT / 'runs' / 'registry', ROOT / 'runs' / 'ledger',
            ROOT / 'runs' / 'archive', storage_root, storage_root / 'worktree-archive', storage_root / 'receipts',
            storage_root / 'trash']


def inside(path, parent):
    try:
        Path(path).relative_to(parent)
        return True
    except ValueError:
        return False


def check(path):
    """(the resolved source, its label inside the trash) or an OSError saying why it may not move."""
    given = Path(os.path.abspath(path))
    if not os.path.lexists(given):
        raise OSError(f'{given}: does not exist')
    if storage.is_link(given):
        raise OSError(f'{given}: is a link or junction; links are never moved')
    real = Path(os.path.realpath(given))
    storage_root = storage.storage_root()
    if inside(real, ROOT):
        label = Path('repo') / real.relative_to(ROOT)
    elif inside(real, storage_root):
        label = Path('storage') / real.relative_to(storage_root)
    else:
        raise OSError(f'{real}: outside the repository and the storage root')
    for place in protected_places():
        place = Path(os.path.realpath(place))
        if inside(place, real) or inside(real, place) and place.name in ('.git', 'refs', 'registry', 'ledger',
                                                                          'archive', 'worktree-archive', 'receipts',
                                                                          'trash'):
            raise OSError(f'{real}: is, holds or lies in the protected place {place}')
    if real.is_dir():
        links = storage.find_links(real)
        if links:
            raise OSError(f'{real}: {len(links)} links or unreadable folders inside, e.g. {links[0]}')
    if inside(given, ROOT):
        relative = given.relative_to(ROOT).as_posix()
        tracked = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '--', relative],
                                 capture_output=True, text=True).stdout.strip()
        if tracked:
            raise OSError(f'{given}: holds tracked files (e.g. {tracked.splitlines()[0]}); use git rm')
    return real, label


def files_of(real):
    if real.is_file():
        return [real]
    return sorted(item for item in storage.long_path(real).rglob('*') if item.is_file())


def move(paths, reason):
    checked = [check(path) for path in paths]  # every path is checked before anything moves
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = trash_root()
    number = 0
    while (root / f'{stamp}-{number}').exists():
        number += 1
    entry = root / f'{stamp}-{number}'
    entry.mkdir(parents=True)
    receipt = {'entry': entry.name, 'reason': reason, 'utc': stamp, 'items': []}
    for real, label in checked:
        target = entry / label
        target.parent.mkdir(parents=True, exist_ok=True)
        files = files_of(real)
        base = storage.long_path(real)
        digests = {(item.relative_to(base).as_posix() if real.is_dir() else real.name): sha256(item) for item in files}
        size = sum(item.stat().st_size for item in files)
        same_drive = os.path.splitdrive(str(real))[0].lower() == os.path.splitdrive(str(entry))[0].lower()
        if same_drive:
            os.rename(storage.long_path(real), storage.long_path(target))
        else:
            if real.is_dir():
                shutil.copytree(storage.long_path(real), storage.long_path(target), symlinks=True)
            else:
                shutil.copy2(storage.long_path(real), storage.long_path(target))
            copied = storage.long_path(target)
            for key, digest in digests.items():
                if sha256(copied / key if real.is_dir() else copied) != digest:
                    raise OSError(f'Copy verification failed for {key}; the source {real} is untouched')
            if real.is_dir():
                storage.remove_tree(real)
            else:
                os.unlink(storage.long_path(real))
        receipt['items'].append({'original': str(real), 'label': label.as_posix(), 'files': len(files), 'bytes': size,
                                 'method': 'rename' if same_drive else 'copy, verify, remove', 'sha256': digests})
    (entry / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    return entry, receipt


def entries():
    root = trash_root()
    return sorted(item for item in root.iterdir() if (item / 'receipt.json').is_file()) if root.is_dir() else []


def restore(name):
    entry = trash_root() / name
    receipt = json.loads((entry / 'receipt.json').read_text(encoding='utf-8'))
    for item in receipt['items']:
        if os.path.lexists(item['original']):
            raise OSError(f'{item["original"]} exists again; nothing restored')
    for item in receipt['items']:
        original = Path(item['original'])
        original.parent.mkdir(parents=True, exist_ok=True)
        source = entry / item['label']
        if os.path.splitdrive(str(original))[0].lower() != os.path.splitdrive(str(entry))[0].lower():
            raise OSError(f'{original}: on another drive; copy it back by hand from {source}')
        os.rename(storage.long_path(source), storage.long_path(original))
    (entry / 'receipt.json').rename(entry / 'restored.json')
    return len(receipt['items'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest='command', required=True)
    moving = commands.add_parser('move', help='move untracked files or folders into the trash')
    moving.add_argument('paths', nargs='+')
    moving.add_argument('--reason', default='')
    commands.add_parser('list', help='entries with their sizes')
    back = commands.add_parser('restore', help='move an entry back where it was')
    back.add_argument('entry')
    args = parser.parse_args()
    if args.command == 'move':
        entry, receipt = move(args.paths, args.reason)
        print(f'{len(receipt["items"])} items, {sum(i["files"] for i in receipt["items"])} files, '
              f'{sum(i["bytes"] for i in receipt["items"]) / 1e6:.1f} MB moved to {entry}')
    elif args.command == 'list':
        total = 0
        for entry in entries():
            receipt = json.loads((entry / 'receipt.json').read_text(encoding='utf-8'))
            size = sum(item['bytes'] for item in receipt['items'])
            total += size
            print(f'{entry.name}  {size / 1e6:9.1f} MB  {receipt["reason"][:60]}  '
                  + ', '.join(item['original'] for item in receipt['items'])[:120])
        print(f'total {total / 1e9:.2f} GB')
    else:
        print(f'{restore(args.entry)} items restored')


if __name__ == '__main__':
    main()
