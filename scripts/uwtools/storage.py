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

from .hashing import sha256
from .paths import ROOT

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


def make_junction(path, target):
    import _winapi
    _winapi.CreateJunction(str(target), str(path))
    if junction_target(path) != Path(target).resolve():
        raise OSError(f'Junction verification failed: {path}')


def copy_verified(source, destination):
    """Copy a directory tree, checking every file's SHA-256; returns {relative path: digest}."""
    digests = {}
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
            shutil.rmtree(backup)
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
        shutil.rmtree(backup)
    receipt['status'] = 'complete'
    receipt['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    write()
    return path


def run_directory(module, identifier):
    """runs/<module>/<id>: the module directory is routed, the run directory is created inside it."""
    directory = route(ROOT / 'runs' / module) / str(identifier)
    directory.mkdir(parents=True, exist_ok=False)
    return directory
