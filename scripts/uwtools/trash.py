"""The agents' trash: untracked data they no longer need is moved here, never deleted.

    uv run python scripts/uwtools/trash.py move PATH... [--reason TEXT]
    uv run python scripts/uwtools/trash.py list
    uv run python scripts/uwtools/trash.py restore ENTRY

The user's rule (2026-10-20): agents may get rid of files as long as code makes an accidental loss impossible.
A delete of untracked data is refused by delete_guard.py; this moves it instead to
<storage root>/trash/<entry>/ (the entry is a UTC time stamp and a number) with a receipt (original path, reason,
files, bytes, SHA-256 of every file). On the same drive the move is a rename, so nothing is copied or deleted; from
another drive every file is copied, verified by digest, and only then removed from its source. Each moved item's
permissions are reset to inherit the trash's. Only the user empties the trash; `restore` copies an entry back and
leaves it in place.
A read-only file or folder is made writable and its removal tried once more. A move that stopped half-way (its
item still 'moving' in an entry's receipt) is finished by naming the same path again, alone: see resume.

Refuses: links (the path itself or anything inside it: a recursive move must never follow a junction,
2026-10-07), tracked files (use git rm; a path on storage is checked at its logical path in the main tree), anything
outside the repository and the storage root, and the protected places: the repository root, .git, refs/, meta/,
runs/registry, runs/ledger, runs/archive, out/, runs/, the storage root, its worktree-archive, receipts and trash, and
the storage-side roots tools open directly (STORAGE_PLACES and every folder of <storage root>/runs, runs/reports and
runs/work; no junctions since 2026-11-01). A relative PATH is logical (storage.data_path: out/scratch/x is on the
storage drive); an absolute PATH is taken as given.
"""
import argparse
import datetime
import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import storage  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402
from uwtools.paths import ROOT, main_root  # noqa: E402


def trash_root():
    root = storage.storage_root()
    if root is None:
        raise OSError('No storage configured (runs/storage.json)')
    return root / 'trash'


# Places nothing may be moved out of (their whole subtree), by name.
SEALED = ('.git', 'refs', 'meta', 'registry', 'ledger', 'archive', 'worktree-archive', 'receipts', 'trash')


# The storage-side roots that are never removed whole (no junctions since 2026-11-01: tools open these folders
# directly), below the storage root; every folder of <storage root>/runs (one per module, reports, work, archive) too.
STORAGE_PLACES = ('worktree-archive', 'receipts', 'trash', 'refs', 'refs/processed', 'refs/source',
                  'refs/processed/physical', 'refs/source/physical', 'out', 'out/scratch', 'out/evaluation', 'runs',
                  'runs/reports', 'runs/work', 'runs/archive')


# The folders whose every direct child is a root tools open directly too (runs/reports/<name> and runs/work/<name>
# were junction targets before 2026-11-01).
STORAGE_PARENTS = ('runs/reports', 'runs/work')


def storage_places(storage_root):
    """The storage root, its STORAGE_PLACES, every folder of its runs/ and every folder of STORAGE_PARENTS (none
    when storage_root is None)."""
    if storage_root is None:
        return []
    storage_root = Path(storage_root)
    places = [storage_root] + [storage_root / place for place in STORAGE_PLACES]
    for parent in [storage_root / 'runs'] + [storage_root / place for place in STORAGE_PARENTS]:
        places += [item for item in sorted(parent.iterdir()) if item.is_dir()] if parent.is_dir() else []
    return places


def protected_places():
    """Places that may not be moved, nor anything holding them; the SEALED ones not even partly."""
    places = [ROOT, ROOT / '.git', ROOT / 'refs', ROOT / 'meta', ROOT / 'runs' / 'registry', ROOT / 'runs' / 'ledger',
              ROOT / 'runs' / 'archive', ROOT / 'out', ROOT / 'runs']
    return places + storage_places(storage.storage_root())


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
        if inside(place, real) or (inside(real, place) and place.name in SEALED):
            raise OSError(f'{real}: is, holds or lies in the protected place {place}')
    if real.is_dir():
        links = storage.find_links(real)
        if links:
            raise OSError(f'{real}: {len(links)} links or unreadable folders inside, e.g. {links[0]}')
    # Tracked files at the path in this checkout, or, for a path on storage, at its logical path in the main tree.
    # git lists the top folder and the names are filtered here, so git never has to follow a link in the path.
    if inside(given, ROOT):
        tree, relative = ROOT, given.relative_to(ROOT).as_posix()
    else:
        tree, relative = main_root(), storage.logical(real)
        if relative is None:
            raise OSError(f'{real}: no logical repository path; nothing moved')
    result = subprocess.run(['git', '-C', str(tree), 'ls-files', '-z', '--', relative.split('/')[0]],
                            capture_output=True, text=True, encoding='utf-8', errors='replace')
    if result.returncode != 0:
        raise OSError(f'{given}: git ls-files failed; nothing moved')
    folded = relative.lower().rstrip('/')
    tracked = [name for name in result.stdout.split('\0')
               if name and (name.lower() == folded or name.lower().startswith(folded + '/'))]
    if tracked:
        raise OSError(f'{given}: holds tracked files (e.g. {tracked[0]}); use git rm')
    return real, label


def files_of(real):
    if real.is_file():
        return [real]
    return sorted(item for item in storage.long_path(real).rglob('*') if item.is_file())


def writable_retry(action, path):
    """action(path); on a PermissionError clear the read-only attribute and try once more (never on a link)."""
    try:
        action(path)
    except PermissionError:
        if storage.is_link(Path(path)):
            raise
        os.chmod(path, stat.S_IWRITE, follow_symlinks=False)
        action(path)


def remove_verified(real, verified):
    """Remove only the files whose copies were verified and that still match them (re-hashed just before), then
    the folders bottom-up; fail if anything is left."""
    base = storage.long_path(real)
    if not real.is_dir():
        if sha256(base) != next(iter(verified.values())):
            raise OSError(f'{real} changed after its copy was verified; the source is kept')
        writable_retry(os.unlink, base)
        return
    for key, digest in verified.items():
        if sha256(base / key) != digest:
            raise OSError(f'{base / key} changed after its copy was verified; it and the rest are kept')
        writable_retry(os.unlink, base / key)
    for folder in sorted((item for item in base.rglob('*') if item.is_dir()), key=lambda item: -len(item.parts)):
        writable_retry(os.rmdir, folder)  # fails on a folder that still holds a file nobody verified
    writable_retry(os.rmdir, base)


def pending(real):
    """(entry, receipt, item) of an interrupted move of real (its item still 'moving'), or None."""
    found = []
    for entry in entries():
        receipt = json.loads((entry / 'receipt.json').read_text(encoding='utf-8'))
        for item in receipt['items']:
            if item.get('status') == 'moving' and os.path.normcase(item['original']) == os.path.normcase(str(real)):
                found.append((entry, receipt, item))
    if len(found) > 1:
        raise OSError(f'{real}: {len(found)} interrupted moves in the trash ({[entry.name for entry, _, _ in found]}); '
                      'sort them out by hand')
    return found[0] if found else None


def finish(target, item, receipt, receipt_path):
    # A rename keeps the item's old permissions; reset them to the trash's inherited ones.
    reset = subprocess.run(['icacls', str(target), '/reset', '/T', '/C', '/Q'], capture_output=True, text=True,
                           encoding='mbcs', errors='replace')
    item['status'] = 'moved' if reset.returncode == 0 else 'moved; permissions not reset: ' + reset.stdout[-200:]
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')


def resume(real, label, entry, receipt, item):
    """Finish an interrupted move of real into its entry, by the receipt's digests (taken before the move began).

    A rename is atomic: done now if the trash copy does not exist yet, refused if both exist. A copy: every file
    left in the source must be in the receipt and still match its digest; every file of the receipt must have a
    copy in the entry matching its digest, copied again from the source where it is missing or wrong (refused when
    the source no longer has it). Only then are the files left in the source removed (remove_verified), and the
    item is marked moved. Anything else is refused with nothing more removed.
    """
    if item['label'] != label.as_posix():
        raise OSError(f'{real}: its interrupted move in {entry.name} has the label {item["label"]}, not {label.as_posix()}')
    receipt_path = entry / 'receipt.json'
    target = entry / label
    copied = storage.long_path(target)
    base = storage.long_path(real)
    digests = item['sha256']
    if item['method'] == 'rename':
        if os.path.lexists(copied):
            raise OSError(f'{real}: both it and {target} exist after an interrupted rename; sort it out by hand')
        os.rename(base, copied)
    else:
        is_dir = real.is_dir()
        left = {(path.relative_to(base).as_posix() if is_dir else real.name): path for path in files_of(real)}
        unknown = [key for key in left if key not in digests]
        if unknown:
            raise OSError(f'{real}: {len(unknown)} files not in the receipt of {entry.name}, e.g. {unknown[0]}; '
                          'nothing more removed')
        changed = [key for key, path in left.items() if sha256(path) != digests[key]]
        if changed:
            raise OSError(f'{real}: {len(changed)} files changed since the move began, e.g. {changed[0]}; '
                          'nothing more removed')
        for key, digest in digests.items():
            destination = copied / key if is_dir else copied
            if destination.is_file() and sha256(destination) == digest:
                continue
            if key not in left:
                raise OSError(f'{key}: gone from the source while its copy in {entry.name} is missing or wrong; '
                              'nothing more removed')
            destination.parent.mkdir(parents=True, exist_ok=True)
            writable_retry(lambda path, source=left[key]: shutil.copy2(source, path), destination)
            if sha256(destination) != digest:
                raise OSError(f'Copy verification failed for {key}; nothing more removed')
        remove_verified(real, {key: digests[key] for key in left})
    item['resumed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    finish(target, item, receipt, receipt_path)
    return entry, receipt


def move(paths, reason):
    if len(paths) == 1 and not os.path.lexists(os.path.abspath(paths[0])):
        # A move whose source is already gone but whose receipt was not updated: a rename that went through, or a
        # copy whose source was removed. Only the permission reset and the status are left; a copy is finished only
        # when every file of the receipt is in the trash and matches its digest. Anything else for a path that does
        # not exist is refused by check below.
        found = pending(Path(os.path.realpath(os.path.abspath(paths[0]))))
        if found is not None:
            entry, receipt, item = found
            target = entry / item['label']
            copied = storage.long_path(target)
            complete = False
            if item['method'] == 'rename':
                complete = os.path.lexists(copied)
            elif copied.is_dir():
                complete = all((copied / key).is_file() and sha256(copied / key) == digest
                               for key, digest in item['sha256'].items())
            elif copied.is_file():
                complete = len(item['sha256']) == 1 and sha256(copied) == next(iter(item['sha256'].values()))
            if complete:
                item['resumed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
                finish(target, item, receipt, entry / 'receipt.json')
                return entry, receipt
            raise OSError(f'{paths[0]}: gone, and its interrupted move in {entry.name} is not complete in the trash '
                          '(a file missing or not matching its digest); nothing done, sort it out by hand')
    checked =[check(path) for path in paths]  # every path is checked before anything moves
    for first, _ in checked:
        for second, _ in checked:
            if first != second and inside(second, first):
                raise OSError(f'{second} lies inside {first}: name each place once')
    interrupted = [(real, label, found) for real, label in checked if (found := pending(real)) is not None]
    if interrupted:
        if len(checked) != 1:
            raise OSError(f'{interrupted[0][0]}: an interrupted move; name it alone to finish it')
        real, label, (entry, receipt, item) = interrupted[0]
        return resume(real, label, entry, receipt, item)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = trash_root()
    number = 0
    while (root / f'{stamp}-{number}').exists():
        number += 1
    entry = root / f'{stamp}-{number}'
    entry.mkdir(parents=True)
    receipt_path = entry / 'receipt.json'
    receipt = {'entry': entry.name, 'reason': reason, 'utc': stamp, 'items': []}
    for real, label in checked:
        target = entry / label
        target.parent.mkdir(parents=True, exist_ok=True)
        files = files_of(real)
        base = storage.long_path(real)
        digests = {(item.relative_to(base).as_posix() if real.is_dir() else real.name): sha256(item) for item in files}
        size = sum(item.stat().st_size for item in files)
        same_drive = os.path.splitdrive(str(real))[0].lower() == os.path.splitdrive(str(entry))[0].lower()
        item = {'original': str(real), 'label': label.as_posix(), 'files': len(files), 'bytes': size,
                'method': 'rename' if same_drive else 'copy, verify, remove', 'status': 'moving', 'sha256': digests}
        receipt['items'].append(item)
        # The receipt is written before each move, so an entry interrupted half-way still says what it holds.
        receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
        if same_drive:
            os.rename(base, storage.long_path(target))
        else:
            if real.is_dir():
                shutil.copytree(base, storage.long_path(target), symlinks=True)
            else:
                shutil.copy2(base, storage.long_path(target))
            copied = storage.long_path(target)
            for key, digest in digests.items():
                if sha256(copied / key if real.is_dir() else copied) != digest:
                    raise OSError(f'Copy verification failed for {key}; the source {real} is untouched')
            remove_verified(real, digests)
        finish(target, item, receipt, receipt_path)
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
        if not item.get('status', '').startswith('moved'):
            raise OSError(f'{item["original"]}: its move did not finish (status {item.get("status")}); '
                          f'copy it back by hand from {entry / item["label"]}')
    for item in receipt['items']:
        # A copy back: the entry stays until the user empties the trash.
        original = Path(item['original'])
        original.parent.mkdir(parents=True, exist_ok=True)
        source = storage.long_path(entry / item['label'])
        if source.is_dir():
            shutil.copytree(source, storage.long_path(original), symlinks=True)
        else:
            shutil.copy2(source, storage.long_path(original))
        copied = storage.long_path(original)
        for key, digest in item['sha256'].items():
            if sha256(copied / key if source.is_dir() else copied) != digest:
                raise OSError(f'restored copy of {key} does not match the receipt')
    receipt['restored_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    (entry / 'receipt.json').write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    return len(receipt['items'])


def argument(text):
    """A path of the command line: a relative one is logical (storage.data_path: out/scratch/x and runs/... on the
    storage drive, out/build/x in this checkout); an absolute one is taken as given, so a leftover copy on C is
    never mistaken for its storage counterpart."""
    return text if Path(text).is_absolute() else str(storage.data_path(text))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    moving = commands.add_parser('move', help='move untracked files or folders into the trash', allow_abbrev=False)
    moving.add_argument('paths', nargs='+', type=argument)
    moving.add_argument('--reason', default='')
    commands.add_parser('list', help='entries with their sizes')
    back = commands.add_parser('restore', help='copy an entry back where it was (the entry stays)')
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
