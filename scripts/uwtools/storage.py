"""Bulk run data on the storage drive: runs/storage.json and directory junctions.

Run directories keep their repository-relative paths; a routed directory is a junction
to the same relative path under the storage root. Without runs/storage.json nothing is
routed and directories are ordinary ones.
"""
import datetime
import json
import os
import re
import shutil
from pathlib import Path

if __package__:
    from .hashing import sha256
    from .paths import ROOT
else:  # run as a script: uv run python scripts/uwtools/storage.py status
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from uwtools.hashing import sha256
    from uwtools.paths import ROOT

CONFIG = ROOT / 'runs/storage.json'
# Only generated data may leave the workspace: one directory level under runs/, or a
# direct child of the two parents that hold one directory per experiment or module.
# runs/archive is sealed (already routed, never written), so it is not routable.
# Run data under runs/, and the bulk scratch and evaluation trees under out/ (routed 2026-10-16 when C filled up).
ROUTABLE = re.compile(r'runs/(reports/[^/]+|work/[^/]+|[a-z][a-z0-9_-]*)|out/(scratch|evaluation)')
RESERVED = {'archive', 'registry', 'ledger', 'reports', 'work'}


def storage_root():
    """The storage root of runs/storage.json, or None when storage is not configured."""
    if not CONFIG.is_file():
        return None
    root = Path(json.loads(CONFIG.read_text(encoding='utf-8-sig'))['storage_root']).resolve()
    if root == ROOT or ROOT in root.parents:
        raise ValueError('Storage must be outside the code workspace')
    if not Path(root.anchor).exists():
        raise OSError(f'Storage drive is unavailable: {root}')
    return root


def relative_path(path):
    relative = Path(os.path.abspath(path)).relative_to(ROOT).as_posix()
    if not ROUTABLE.fullmatch(relative) or relative.split('/')[1] in RESERVED and relative.count('/') == 1:
        raise ValueError(f'Only run data directories are routed to storage: {relative}')
    return relative


def junction_target(path):
    """Resolved target of a junction, or None for anything else."""
    if not os.path.isjunction(path):
        return None
    # readlink returns the extended-length form (\\?\D:\...).
    return Path(os.readlink(path).removeprefix('\\\\?\\')).resolve()


def long_path(path):
    """The extended-length form (\\\\?\\C:\\...): files below it stay reachable past 260 characters.

    Without it a deep file is invisible to exists(), open() and rglob's is_file() (2026-10-17: the clean-up
    of out/scratch stopped on one, and a scan without it would skip such files silently).
    """
    text = os.path.abspath(path)
    return Path(text if text.startswith('\\\\?\\') else '\\\\?\\' + text)


def is_link(entry):
    """A symbolic link, a junction or any other reparse point (os.DirEntry or Path)."""
    if isinstance(entry, Path):
        if entry.is_symlink() or os.path.isjunction(entry):
            return True
        return bool(os.lstat(entry).st_file_attributes & 0x400)
    return entry.is_symlink() or entry.is_junction() or bool(entry.stat(follow_symlinks=False).st_file_attributes & 0x400)


def find_links(root):
    """Every link below root, never followed (an unreadable directory counts as one: it cannot be checked)."""
    found, stack = [], [str(long_path(root))]
    while stack:
        directory = stack.pop()
        try:
            entries = list(os.scandir(directory))
        except OSError:
            found.append(directory)
            continue
        for entry in entries:
            if is_link(entry):
                found.append(entry.path)
            elif entry.is_dir(follow_symlinks=False):
                stack.append(entry.path)
    return found


def remove_tree(root):
    """Delete a directory tree that holds no link (2026-10-07: a recursive delete followed a junction into refs/).

    Scans the whole tree first and deletes nothing if it finds a link; while deleting it never descends
    into a link and stops at one that appeared since the scan.
    """
    root = long_path(root)
    if is_link(root):
        raise OSError(f'Refusing to delete through a link: {root}')
    links = find_links(root)
    if links:
        raise OSError(f'Refusing to delete {root}: {len(links)} links or unreadable directories inside, e.g. {links[0]}')

    def remove(directory):
        for entry in os.scandir(directory):
            if is_link(entry):
                raise OSError(f'A link appeared during the delete: {entry.path}')
            if entry.is_dir(follow_symlinks=False):
                remove(entry.path)
            else:
                try:
                    os.unlink(entry.path)
                except PermissionError:
                    os.chmod(entry.path, 0o666)
                    os.unlink(entry.path)
        os.rmdir(directory)

    remove(root)


def make_junction(path, target):
    import _winapi
    _winapi.CreateJunction(str(target), str(path))
    if junction_target(path) != Path(target).resolve():
        raise OSError(f'Junction verification failed: {path}')


def copy_verified(source, destination):
    """Copy a directory tree, checking every file's SHA-256; returns {relative path: digest}."""
    digests = {}
    source, destination = long_path(source), long_path(destination)
    items = list(source.rglob('*'))
    links = [item for item in items if os.path.isjunction(item) or item.is_symlink()]
    if links:
        raise ValueError(f'Nested links require separate review: {links[0]}')
    for file in sorted(item for item in items if item.is_file()):
        relative = file.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(file, target)
        digests[relative.as_posix()] = sha256(file)
        if sha256(target) != digests[relative.as_posix()]:
            raise OSError(f'Copy verification failed: {file}')
    for directory in (item for item in items if item.is_dir()):
        (destination / directory.relative_to(source)).mkdir(parents=True, exist_ok=True)
    return digests


def route(path):
    """Make path a junction into the storage root, moving existing contents there.

    Returns the logical path (never the resolved one: helpers must keep seeing the
    workspace path). A directory that is already routed is left alone.
    """
    path = Path(os.path.abspath(path))
    root = storage_root()
    if root is None:
        path.mkdir(parents=True, exist_ok=True)
        return path
    relative = relative_path(path)
    destination = root / relative
    current = junction_target(path)
    if current is not None:
        if current != destination.resolve():
            raise OSError(f'Unexpected existing link: {path} -> {current}')
        return path
    if path.exists() and not path.is_dir():
        raise NotADirectoryError(path)
    receipt_path = root / 'receipts' / (relative.replace('/', '__') + '.json')
    if destination.exists() and any(destination.iterdir()) and not receipt_path.is_file():
        raise FileExistsError(f'Destination already exists without a receipt: {destination}')
    if receipt_path.is_file() and path.is_dir():
        # A verified copy whose swap was refused (an open handle, 2026-10-17): finish it without copying
        # again, unless something on the source side is newer than the receipt.
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if receipt.get('status') == 'verified':
            stamp = receipt_path.stat().st_mtime
            newer = [item for item in path.rglob('*') if item.is_file() and item.stat().st_mtime > stamp]
            if newer:
                raise OSError(f'{len(newer)} files newer than the verified copy, e.g. {newer[0]}: route again after removing them or re-copy')
            backup = path.with_name(path.name + '.storage-migration-backup')
            path.rename(backup)
            make_junction(path, destination)
            remove_tree(backup)
            receipt['status'] = 'complete'
            receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
            receipt_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
            return path
    destination.mkdir(parents=True, exist_ok=True)
    receipt = {'source': str(path), 'destination': str(destination), 'status': 'copying', 'sha256': {}}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    write = lambda: receipt_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    write()
    backup = path.with_name(path.name + '.storage-migration-backup')
    if backup.exists():
        raise FileExistsError(f'Unexpected backup exists: {backup}')
    if path.exists():
        receipt['sha256'] = copy_verified(path, destination)
        receipt['status'] = 'verified'
        write()
        path.rename(backup)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        make_junction(path, destination)
    except Exception:
        if backup.exists() and not path.exists():
            backup.rename(path)
        raise
    if backup.exists():
        remove_tree(backup)
    receipt['status'] = 'complete'
    receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    write()
    return path


def run_directory(module, identifier):
    """runs/<module>/<id>: the module directory is routed, the run directory is created inside it."""
    directory = route(ROOT / 'runs' / module) / str(identifier)
    directory.mkdir(parents=True, exist_ok=False)
    return directory


def check_backup(path):
    """Compare a leftover <path>.storage-migration-backup with its receipt and the storage copy.

    Returns counts: ok (the file is unchanged since the copy and the storage copy has the same digest),
    newer_on_storage (the storage copy was rewritten since, by a later run), and the blocking kinds
    not_in_receipt, changed_since_copy, missing_on_storage, links.
    """
    path = Path(os.path.abspath(path))
    relative = relative_path(path)
    backup = long_path(path.with_name(path.name + '.storage-migration-backup'))
    receipt = json.loads((storage_root() / 'receipts' / (relative.replace('/', '__') + '.json')).read_text(encoding='utf-8'))
    destination, digests = long_path(receipt['destination']), receipt['sha256']
    counts = dict.fromkeys(('files', 'ok', 'newer_on_storage', 'not_in_receipt', 'changed_since_copy',
                            'missing_on_storage', 'links'), 0)
    examples = {}
    counts['links'] = len(find_links(backup))
    for file in (item for item in backup.rglob('*') if item.is_file() and not is_link(item)):
        key = file.relative_to(backup).as_posix()
        counts['files'] += 1
        target = destination / key
        kind = ('not_in_receipt' if key not in digests else 'changed_since_copy' if sha256(file) != digests[key]
                else 'missing_on_storage' if not target.is_file() else 'ok' if sha256(target) == digests[key]
                else 'newer_on_storage')
        counts[kind] += 1
        examples.setdefault(kind, key)
    return backup, receipt, counts, examples


def finish(path, accept_newer=False):
    """Delete the leftover backup of a routed directory once every file in it is verified on storage.

    A swap whose clean-up stopped (an open handle, 2026-10-17) leaves the junction in place and the
    old tree beside it. Refuses unless every file is unchanged since the receipt's copy and the storage
    copy matches its digest; a storage copy rewritten since by a later run passes only with accept_newer.
    A backup already gone (removed by hand) only marks the receipt complete, once the junction is verified.
    """
    path = Path(os.path.abspath(path))
    receipt_path = storage_root() / 'receipts' / (relative_path(path).replace('/', '__') + '.json')
    if not path.with_name(path.name + '.storage-migration-backup').exists():
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if junction_target(path) != Path(receipt['destination']).resolve():
            raise OSError(f'{path} is not the junction to {receipt["destination"]}')
        receipt['status'] = 'complete'
        receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
        receipt['backup_check'] = 'backup already absent; nothing removed'
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
        return receipt['backup_check']
    backup, receipt, counts, examples = check_backup(path)
    blocking = [kind for kind in ('not_in_receipt', 'changed_since_copy', 'missing_on_storage', 'links') if counts[kind]]
    if not accept_newer and counts['newer_on_storage']:
        blocking.append('newer_on_storage')
    if junction_target(path) != Path(receipt['destination']).resolve():
        raise OSError(f'{path} is not the junction to {receipt["destination"]}')
    if blocking:
        raise OSError(f'Backup not removed: {counts}; first of each: {examples}')
    remove_tree(backup)
    receipt['status'] = 'complete'
    receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    receipt['backup_check'] = counts
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return counts


# A session starts with a warning below this much free space on the workspace drive.
LOW_SPACE_GB = 60


def status():
    """Lines describing the workspace drive: free space, leftover backups, unrouted bulk directories, worktrees."""
    free = shutil.disk_usage(ROOT).free / 1e9
    lines = [f'C free {free:.0f} GB' + (f' (below {LOW_SPACE_GB} GB)' if free < LOW_SPACE_GB else '')]
    for parent in (ROOT / 'out', ROOT / 'runs', ROOT / 'runs' / 'reports'):
        for item in parent.glob('*.storage-migration-backup'):
            lines.append(f'leftover backup {item.relative_to(ROOT).as_posix()}: '
                         f'storage.py finish {item.relative_to(ROOT).as_posix().removesuffix(".storage-migration-backup")}')
    for name in ('scratch', 'evaluation'):
        if junction_target(ROOT / 'out' / name) is None and (ROOT / 'out' / name).exists():
            lines.append(f'out/{name} is a directory on C, not the junction: storage.py route out/{name}')
    worktrees = ROOT / '.claude' / 'worktrees'
    count = sum(1 for item in worktrees.iterdir() if item.is_dir()) if worktrees.is_dir() else 0
    if count:
        lines.append(f'{count} agent worktrees in .claude/worktrees (about 1.5-3 GB each): scripts/uwtools/worktrees.py audit')
    return lines


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Bulk data on the storage drive (runs/storage.json).', allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('status', help='free space and pending storage work')
    for name, text in (('route', 'make a directory a junction into storage, moving its contents'),
                       ('check', 'compare a leftover backup with its receipt and the storage copy (read-only)'),
                       ('finish', 'delete a leftover backup once every file in it is verified on storage')):
        command = commands.add_parser(name, help=text, allow_abbrev=False)
        command.add_argument('path')
        if name == 'finish':
            command.add_argument('--accept-newer', action='store_true',
                                 help='also when a later run rewrote the storage copy of some files')
    args = parser.parse_args()
    if args.command == 'status':
        print('\n'.join(status()))
    elif args.command == 'route':
        print(route(ROOT / args.path))
    elif args.command == 'check':
        print(check_backup(ROOT / args.path)[2:])
    else:
        print(finish(ROOT / args.path, args.accept_newer))


if __name__ == '__main__':
    main()
