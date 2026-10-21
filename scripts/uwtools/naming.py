"""One naming scheme for experiments (master plan 2.2): a numeric id plus a readable tag.

    uv run --offline python scripts/uwtools/naming.py new KIND TAG [--module M] [--label TEXT]

Registered module runs keep runs/<module>/<id> and the ids of runs/registry/<module>.json (kind
`run`; experiment.py appends their record, which now carries the tag). Every other directory a tool
creates is `<parent>/<id>-<tag>`, its id drawn from one sequence shared by all those kinds and
recorded in runs/registry/scratch.json:

    sandbox      out/scratch/sandbox             terrain_sandbox.py
    step-check   out/scratch/step-check          step_check.py
    geology      out/scratch/geology             uw.py geology run
    experiment   out/evaluation/<module>         experiment.py without --register or --archive
    page         out/evaluation                  acceptance.py without --output

The registry and the directories are always the main tree's (paths.main_root()), never a worktree's
copy, so parallel workers share one sequence; an OS byte-range lock on out/scratch/naming.lock (a file
that stays; the system releases the lock when its process ends, so a lock is never stale) serialises the
appends. A new directory must resolve under the storage
root of runs/storage.json (D). The next id is one more than the registry's highest and every numeric
prefix under the kinds' parent folders. A tag is lower case, `^[a-z][a-z0-9]*(-[a-z0-9]+)*$`, at most
40 characters, and does not end in `-<digits>` or `-v<digits>` (the id is the version). Labels written
before the scheme start with a letter and are not renamed; find() still resolves them. The topic of a
directory is (kind, module, tag); previous() gives the topic's newest earlier record.

`new` prints the created path, for hand-run tools (world_pipeline --bundle <path>/bundle).
"""
import argparse
import datetime
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

if __package__:
    from . import paths, runs, storage
else:  # run as a script: uv run python scripts/uwtools/naming.py new ...
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from uwtools import paths, runs, storage

REGISTRY = 'scratch'
PARENTS = {'sandbox': 'out/scratch/sandbox', 'step-check': 'out/scratch/step-check', 'geology': 'out/scratch/geology',
           'experiment': 'out/evaluation/{module}', 'page': 'out/evaluation'}
KINDS = ('run',) + tuple(PARENTS)
# Step checks written before the scheme: flat <parent>/<label>-<recipe> folders.
LEGACY_STEP_CHECKS = ('out/scratch/step-check', 'out/scratch/sphere/step-check', 'out/scratch/phase5')
LEGACY_STEP_CHECK_RECIPES = ('procedural-world', 'earth-climate', 'tectonics', 'social')
TAG = re.compile(r'[a-z][a-z0-9]*(-[a-z0-9]+)*')
VERSION_SUFFIX = re.compile(r'-v?[0-9]+')
TAG_LENGTH = 40
NUMBERED = re.compile(r'([0-9]+)(?:-|$)')
LOCK = 'out/scratch/naming.lock'
LOCK_TIMEOUT_SECONDS = 600.0
_main_root = None


@dataclass(frozen=True)
class Named:
    id: int
    kind: str
    tag: str
    module: str
    name: str
    path: Path


def main_root():
    """The main tree (looked up once): the registry and every named directory live there."""
    global _main_root
    if _main_root is None:
        _main_root = paths.main_root()
    return _main_root


def check_tag(tag):
    """tag itself; ValueError when it breaks the tag rule."""
    if not isinstance(tag, str) or not TAG.fullmatch(tag) or len(tag) > TAG_LENGTH:
        raise ValueError(f'Tag {tag!r}: lower-case words of letters and digits joined by "-", starting with a letter, '
                         f'at most {TAG_LENGTH} characters')
    match = VERSION_SUFFIX.search(tag)
    if match and match.end() == len(tag):
        raise ValueError(f'Tag {tag!r} ends in a version ({match.group(0)}); the id is the version, drop the suffix')
    return tag


def _check_kind(kind, module):
    if kind not in KINDS:
        raise ValueError(f'Kind {kind!r}: one of {", ".join(KINDS)}')
    if kind in ('run', 'experiment') and not module:
        raise ValueError(f'Kind {kind} needs a module')


def _parent(kind, module, root):
    return root / PARENTS[kind].format(module=module)


def _under(path, store):
    real = Path(os.path.realpath(path))
    store = Path(os.path.realpath(store))
    return real == store or store in real.parents


def _storage():
    store = storage.storage_root()
    if store is None:
        raise SystemExit('runs/storage.json names no storage root: named directories are written to the storage drive')
    return store


def _storage_parent(parent, root, store):
    """parent created inside the routed out/scratch or out/evaluation of the main tree."""
    top = root.joinpath(*parent.relative_to(root).parts[:2])
    if not _under(top, store):
        if Path(os.path.abspath(storage.ROOT)) != Path(os.path.abspath(root)):
            raise SystemExit(f'{top} is not routed to the storage drive; route it from the main tree '
                             '(scripts/uwtools/storage.py route)')
        storage.route(top)
    parent.mkdir(parents=True, exist_ok=True)
    if not _under(parent, store):
        raise SystemExit(f'{parent} does not resolve under the storage root {store}')
    return parent


def _code_state():
    def git(*arguments):
        return subprocess.run(['git', *arguments], cwd=paths.ROOT, capture_output=True, text=True)
    revision = git('rev-parse', 'HEAD')
    diff = git('diff', 'HEAD', '--quiet')
    return {'revision': revision.stdout.strip() if revision.returncode == 0 else None, 'dirty': diff.returncode != 0}


def _tool():
    return Path(sys.argv[0]).stem if sys.argv and sys.argv[0] else 'python'


def _numbered(directory):
    if not directory.is_dir():
        return []
    return [int(match.group(1)) for entry in directory.iterdir()
            if entry.is_dir() and (match := NUMBERED.match(entry.name))]


def next_id(root):
    """One more than the scratch registry's highest id and every numeric prefix under the kinds' parents."""
    used = [entry['id'] for entry in runs.registry(REGISTRY, root=root)['runs']]
    for kind in ('sandbox', 'step-check', 'geology', 'page'):
        used += _numbered(_parent(kind, None, root))
    evaluation = root / 'out/evaluation'
    if evaluation.is_dir():
        for module in evaluation.iterdir():
            if module.is_dir() and not NUMBERED.match(module.name):
                used += _numbered(module)
    return max(used, default=0) + 1


class _Lock:
    """A byte-range lock (msvcrt.locking) on the first byte of out/scratch/naming.lock. The file is created
    when missing and never removed; the system releases the lock when the handle closes or the process
    ends. Waits up to LOCK_TIMEOUT_SECONDS, then raises TimeoutError."""

    def __init__(self, root):
        self.path = root / LOCK
        self.handle = None

    def __enter__(self):
        import msvcrt
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = os.open(self.path, os.O_CREAT | os.O_RDWR)
        deadline = time.monotonic() + LOCK_TIMEOUT_SECONDS
        while True:
            os.lseek(self.handle, 0, os.SEEK_SET)
            try:
                msvcrt.locking(self.handle, msvcrt.LK_NBLCK, 1)
                return self
            except OSError:
                if time.monotonic() >= deadline:
                    os.close(self.handle)
                    self.handle = None
                    raise TimeoutError(f'{self.path} stayed locked for {LOCK_TIMEOUT_SECONDS:g} s')
                time.sleep(0.1)

    def __exit__(self, *exception):
        import msvcrt
        try:
            os.lseek(self.handle, 0, os.SEEK_SET)
            msvcrt.locking(self.handle, msvcrt.LK_UNLCK, 1)
        finally:
            os.close(self.handle)
            self.handle = None


def new_directory(kind, tag, *, module=None, label=''):
    """A new, empty directory of a kind, named by the next id and the tag (see the module docstring)."""
    _check_kind(kind, module)
    check_tag(tag)
    store = _storage()
    if kind == 'run':
        identifier = runs.next_id(module)
        path = storage.run_directory(module, identifier)
        if not _under(path, store):
            raise SystemExit(f'{path} does not resolve under the storage root {store}')
        return Named(identifier, kind, tag, module, str(identifier), path)
    root = main_root()
    parent = _storage_parent(_parent(kind, module, root), root, store)
    with _Lock(root):
        identifier = next_id(root)
        name = f'{identifier}-{tag}'
        path = parent / name
        path.mkdir(exist_ok=False)
        entry = {'id': identifier, 'kind': kind, 'tag': tag}
        if module:
            entry['module'] = module
        entry.update(run=path.relative_to(root).as_posix(),
                     date=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
                     label=label or '', tool=_tool(), code=_code_state())
        runs.append_record(REGISTRY, entry, root=root)
    return Named(identifier, kind, tag, module, name, path)


def records(kind=None, module=None, root=None):
    """The scratch registry's records, oldest first, of a kind (and module) when given."""
    root = root or main_root()
    return [entry for entry in runs.registry(REGISTRY, root=root)['runs']
            if (kind is None or entry.get('kind') == kind) and (module is None or entry.get('module') == module)]


def _legacy_step_check(ref, root):
    for parent in LEGACY_STEP_CHECKS:
        folder = root / parent
        if (folder / ref).is_dir():
            return folder / ref
        if any((folder / f'{ref}-{recipe}').is_dir() for recipe in LEGACY_STEP_CHECK_RECIPES):
            return folder / ref  # the prefix: its bundles are <prefix>-<recipe>
    return None


def find(kind, ref, module=None):
    """The directory of ref: an absolute path as given; digits, the id; `<id>-<tag>`, the new name; else a
    label written before the scheme in the kind's parent. A flat step check (<parent>/<label>-<recipe>, in
    out/scratch/step-check, out/scratch/sphere/step-check or out/scratch/phase5) is returned as the prefix
    <parent>/<label>, a path that is not itself a directory. FileNotFoundError when nothing matches."""
    ref = str(ref)
    if Path(ref).is_absolute():
        return Path(ref)
    _check_kind(kind, module)
    if kind == 'run':
        if not ref.isdigit():
            raise FileNotFoundError(f'A registered {module} run is named by its id, not {ref!r}')
        path = runs.directory(module, int(ref))
        if path is None:
            raise FileNotFoundError(f'{module} run {ref} is sealed under runs/archive')
        return path
    root = main_root()
    parent = _parent(kind, module, root)
    if ref.isdigit():
        for entry in records(kind, module if kind == 'experiment' else None, root):
            if entry['id'] == int(ref):
                return root / entry['run']
        found = [entry for entry in parent.glob(f'{ref}-*') if entry.is_dir()] if parent.is_dir() else []
        if len(found) == 1:
            return found[0]
        raise FileNotFoundError(f'No {kind} with id {ref} in {REGISTRY}.json or {parent}')
    if (parent / ref).is_dir():
        return parent / ref
    if kind == 'step-check':
        legacy = _legacy_step_check(ref, root)
        if legacy is not None:
            return legacy
    raise FileNotFoundError(f'No {kind} named {ref!r} in {parent}')


def previous(named):
    """The newest earlier record of named's topic (kind, module, tag), or None."""
    if named.kind == 'run':
        earlier = [entry for entry in runs.registry(named.module)['runs']
                   if entry.get('tag') == named.tag and entry['id'] < named.id]
    else:
        earlier = [entry for entry in records(named.kind)
                   if entry.get('module') == named.module and entry.get('tag') == named.tag and entry['id'] < named.id]
    return max(earlier, key=lambda entry: entry['id'], default=None)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    new = commands.add_parser('new', help='create <parent>/<id>-<tag>, record it and print its path')
    new.add_argument('kind', choices=tuple(PARENTS))
    new.add_argument('tag')
    new.add_argument('--module', help='the module (required for experiment)')
    new.add_argument('--label', default='', help='a description stored with the record')
    arguments = parser.parse_args(argv)
    try:
        named = new_directory(arguments.kind, arguments.tag, module=arguments.module, label=arguments.label)
    except ValueError as error:
        parser.error(str(error))
    print(named.path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
