"""Bulk data on the storage drive (runs/storage.json): where a logical repository path is on disk.

Since 2026-11-01 (the user) nothing is linked: generated data and the large reference trees live under the storage
root (D) at the same relative path as their logical repository path, and tools open them there through data_path.
Records (runs/registry/*.json) keep logical repository-relative paths; logical() maps a real path back. Without
runs/storage.json nothing is on storage and data_path refuses the storage-side paths.

    uv run --offline python scripts/uwtools/storage.py status
"""
import functools
import json
import os
import re
import shutil
import stat
from pathlib import Path

if __package__:
    from .paths import ROOT, main_root
else:  # run as a script: uv run python scripts/uwtools/storage.py status
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from uwtools.paths import ROOT, main_root

CONFIG = ROOT / 'runs/storage.json'


def storage_root(config=None):
    """The storage root of a runs/storage.json (default this checkout's CONFIG; naming.py passes the main
    tree's), or None when storage is not configured."""
    config = Path(CONFIG if config is None else config)
    if not config.is_file():
        return None
    root = Path(json.loads(config.read_text(encoding='utf-8-sig'))['storage_root']).resolve()
    for workspace in (ROOT, config.resolve().parents[1]):
        if root == workspace or workspace in root.parents:
            raise ValueError('Storage must be outside the code workspace')
    if not Path(root.anchor).exists():
        raise OSError(f'Storage drive is unavailable: {root}')
    return root


# Since 2026-11-01 (the user) there are no junctions: generated data and the large reference trees live under the
# storage root at the same relative path as their logical repository path, and tools open them there. On storage:
# runs/<anything but registry, ledger, storage.json>, out/scratch, out/evaluation, refs/processed/physical,
# refs/source/physical (and everything below them).
ON_STORAGE = re.compile(r'(runs/(?!(registry|ledger|storage\.json)(/|$))[^/]+|out/(scratch|evaluation)|'
                        r'refs/(processed|source)/physical)(/.*)?')


@functools.lru_cache(maxsize=None)
def _main():
    return main_root()


def on_storage(relative):
    """True when a logical repository path (posix, relative) lives under the storage root."""
    return bool(ON_STORAGE.fullmatch(Path(relative).as_posix().strip('/')))


def data_path(relative):
    """Where a logical repository path (runs/erosion/9, out/scratch/x, refs/processed/physical/y) is on disk:
    under the storage root when on_storage, else in the main tree. An absolute path comes back unchanged."""
    if Path(relative).is_absolute():
        return Path(relative)
    relative = Path(relative).as_posix().strip('/')
    if on_storage(relative):
        storage = storage_root(_main() / 'runs/storage.json')
        if storage is None:
            raise OSError('runs/storage.json is missing: data on storage cannot be found')
        return storage / relative
    return _main() / relative


def logical(path):
    """The logical repository path (posix) of a path under the storage root or the main tree, or None."""
    path = Path(os.path.abspath(path))
    storage = storage_root(_main() / 'runs/storage.json')
    for base in ([storage] if storage else []) + [_main(), ROOT]:
        try:
            return path.relative_to(base).as_posix()
        except ValueError:
            continue
    return None


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
    into a link and stops at one that appeared since the scan. A read-only file or folder (git packs, a
    read-only directory) is made writable and its delete tried once more; a link is never made writable.
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
                    os.chmod(entry.path, stat.S_IWRITE, follow_symlinks=False)
                    os.unlink(entry.path)
        try:
            os.rmdir(directory)
        except PermissionError:
            if is_link(Path(directory)):
                raise
            os.chmod(directory, stat.S_IWRITE, follow_symlinks=False)
            os.rmdir(directory)

    remove(root)


# A session starts with a warning below this much free space on the workspace drive.
LOW_SPACE_GB = 60
# status lists the folders under out/ other than these that hold more than OUT_LISTED_MB.
OUT_KEPT = ('build', 'scratch', 'evaluation', 'inputs')  # out/inputs: derived inputs (scripts/evaluation/inputs.py)
OUT_LISTED_MB = 10
COMMAND = 'uv run --offline python scripts/uwtools/'
# The logical places outside runs/ whose data belongs on storage (ON_STORAGE); every child of runs/ but registry,
# ledger and storage.json too.
STORAGE_PLACES = ('out/scratch', 'out/evaluation', 'refs/processed/physical', 'refs/source/physical')


def folder_bytes(path):
    """The bytes of the files below path, links never followed and unreadable entries skipped."""
    total, stack = 0, [str(long_path(path))]
    while stack:
        try:
            entries = list(os.scandir(stack.pop()))
        except OSError:
            continue
        for entry in entries:
            try:
                if is_link(entry):
                    continue
                if entry.is_dir(follow_symlinks=False):
                    stack.append(entry.path)
                else:
                    total += entry.stat(follow_symlinks=False).st_size
            except OSError:
                continue
    return total


def plain_folders(parent):
    """The folders of parent that are not links."""
    if not parent.is_dir():
        return []
    return [item for item in sorted(parent.iterdir()) if item.is_dir() and not is_link(item)]


def misplaced():
    """Lines for data on C that belongs on storage: a path of this checkout at a storage place (a child of runs/ other
    than registry, ledger and storage.json, out/scratch, out/evaluation, refs/processed/physical,
    refs/source/physical), a link there (no junctions since 2026-11-01) or a folder holding links."""
    places = [f'runs/{item.name}' for item in sorted((ROOT / 'runs').iterdir())] if (ROOT / 'runs').is_dir() else []
    places = [place for place in places + list(STORAGE_PLACES) if on_storage(place) and os.path.lexists(ROOT / place)]
    lines = []
    for relative in places:
        path = ROOT / relative
        try:
            target = data_path(relative)
        except (OSError, ValueError) as error:
            target = f'the storage root ({error})'
        if is_link(path):
            lines.append(f'{relative} is a link on C (no junctions since 2026-11-01): tools read {target} directly')
        elif path.is_dir():
            links = find_links(path)
            inside = f', {len(links)} links inside, e.g. {links[0]}' if links else ''
            lines.append(f'{relative} is a folder on C ({folder_bytes(path) / 1e6:.0f} MB{inside}): its data belongs at '
                         f'{target}')
        else:
            lines.append(f'{relative} is a file on C: it belongs at {target}')
    return lines


def primary_clutter():
    """Lines for the bulky folders on C that are no storage place: folders under out/ other than OUT_KEPT over
    OUT_LISTED_MB, and out/build folders that are neither a CMakePresets.json build nor the tracking build."""
    lines = []
    for item in plain_folders(ROOT / 'out'):
        if item.name in OUT_KEPT:
            continue
        size = folder_bytes(item)
        if size > OUT_LISTED_MB * 1e6:
            lines.append(f'out/{item.name} on C ({size / 1e6:.0f} MB): when no longer needed, '
                         f'{COMMAND}trash.py move out/{item.name} --reason TEXT')
    presets = preset_builds()
    for item in plain_folders(ROOT / 'out' / 'build') if presets else []:  # no presets read: nothing to judge by
        if item.name in presets or tracking_build(item):
            continue
        lines.append(f'out/build/{item.name} is not a build of CMakePresets.json ({folder_bytes(item) / 1e6:.0f} MB): '
                     f'when no longer needed, {COMMAND}trash.py move out/build/{item.name} --reason TEXT')
    return lines


def preset_builds():
    """The names of the out/build/<name> folders CMakePresets.json's configure presets build in (binaryDir), or None
    when the file cannot be read."""
    try:
        document = json.loads((ROOT / 'CMakePresets.json').read_text(encoding='utf-8-sig'))
    except (OSError, ValueError):
        return None
    names = set()
    for preset in document.get('configurePresets', []) if isinstance(document, dict) else []:
        text = str(preset.get('binaryDir', '')).replace('\\', '/').rstrip('/') if isinstance(preset, dict) else ''
        if '/out/build/' in f'/{text}':
            names.add(text.rsplit('/', 1)[-1])
    return names


def tracking_build(folder):
    """True for the field-access tracking build (its CMakeCache.txt sets UW_TRACK_FIELD_READS on; AGENTS.md: measure
    access only when a change alters which fields a process reads or writes), which is kept beside the presets."""
    try:
        return 'UW_TRACK_FIELD_READS:BOOL=ON' in (folder / 'CMakeCache.txt').read_text(encoding='utf-8', errors='replace')
    except OSError:
        return False


def status():
    """Lines describing the workspace drive: free space, data on C that belongs on storage, bulky folders, worktrees."""
    free = shutil.disk_usage(ROOT).free / 1e9
    lines = [f'C free {free:.0f} GB' + (f' (below {LOW_SPACE_GB} GB)' if free < LOW_SPACE_GB else '')]
    lines += misplaced()
    lines += primary_clutter()
    worktrees = ROOT / '.claude' / 'worktrees'
    count = sum(1 for item in worktrees.iterdir() if item.is_dir()) if worktrees.is_dir() else 0
    if count:
        lines.append(f'{count} agent worktrees in .claude/worktrees (about 1.5-3 GB each): scripts/uwtools/worktrees.py audit')
    return lines


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Bulk data on the storage drive (runs/storage.json).', allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('status', help='free space and data on C that belongs on storage')
    parser.parse_args()
    print('\n'.join(status()))


if __name__ == '__main__':
    main()
