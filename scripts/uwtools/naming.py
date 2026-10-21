"""One naming scheme for experiments (master plan 2.2): a numeric id plus a readable tag.

    uv run --offline python scripts/uwtools/naming.py new KIND TAG [--module M] [--label TEXT]

Registered module runs keep runs/<module>/<id> and the ids of runs/registry/<module>.json (kind
`run`; experiment.py appends their record, which now carries the tag, to the main tree's registry,
also from a worktree). Every other directory a tool
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

Automatic moves (master plan 2.4). finish(named, passed, inputs) runs at the end of every tool that names its
directory (experiment.py, step_check.py, terrain_sandbox.py, nested_window.py, acceptance.py, compare_runs.py,
hypsometry.py, terrain_metrics.py, uw.py geology run). When the new run passed, every earlier run of its topic is
archived (scripts/uwtools/archive.py: its re-renderable map PNGs and pages dropped, with a receipt), so a failed run
in between leaves nothing behind; a folder archive.py does not take as a run (a sandbox, step-check, geology or
page folder: nothing in it is re-renderable) gets `superseded_by` and `superseded_utc` on its scratch.json record
instead. Nothing happens to an earlier run that is pinned: the registry's `baseline` of a module, an evaluation
awaiting its review (verdict.review 'pending' or a review-request.json, and no verdict in review.json), a path
named in the main tree's (or this checkout's) configs/sandbox/*.json or configs/evaluation/*.json (also as the
value of a key="path" setting), or an input of the new run (its declared inputs, --baseline, the step check's
reference, --parent-run, the sandbox's base bundle); nor to one already archived (archive-receipt.json) or
superseded; nor to one without its completion marker (COMPLETE: ledger-entry.json, result.json, config.json, report.json or
index.html by kind; a run still in progress); nor to a step check of another grid or --social setting. A failed run archives nothing, and an earlier run that
cannot be archived is reported and left as it is: the finished run never fails because of it.

tool_environment(directory) is the environment of a headless tool whose own outputs (climate diagnostics, images,
worlds; the UW_*_DIR variables) go to directory on the storage drive instead of the pipeline's defaults under out/
on C. page_directory(named) is out/evaluation/<name> of the main tree, where a scratch run's HTML page goes.
"""
import argparse
import copy
import datetime
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

if __package__:
    from . import archive, paths, runs, storage
else:  # run as a script: uv run python scripts/uwtools/naming.py new ...
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from uwtools import archive, paths, runs, storage

REGISTRY = 'scratch'
PARENTS = {'sandbox': 'out/scratch/sandbox', 'step-check': 'out/scratch/step-check', 'geology': 'out/scratch/geology',
           'experiment': 'out/evaluation/{module}', 'page': 'out/evaluation'}
KINDS = ('run',) + tuple(PARENTS)
# Step checks written before the scheme: flat <parent>/<label>-<recipe> folders.
LEGACY_STEP_CHECKS = ('out/scratch/step-check', 'out/scratch/sphere/step-check', 'out/scratch/phase5')
LEGACY_STEP_CHECK_RECIPES = ('procedural-world', 'earth-climate', 'tectonics', 'social')
TAG = re.compile(r'[a-z][a-z0-9]*(-[a-z0-9]+)*')
VERSION_SUFFIX = re.compile(r'-v?[0-9]+$')
TAG_LENGTH = 40
NUMBERED = re.compile(r'([0-9]+)(?:-|$)')
LOCK = 'out/scratch/naming.lock'
LOCK_TIMEOUT_SECONDS = 600.0
# Where a headless tool writes its own outputs (src: app.env, command_line_generation.cpp), below tool_environment's
# directory.
TOOL_OUTPUTS = {'UW_CLIMATE_DIAGNOSTICS_DIR': 'diagnostics', 'UW_DEFAULT_IMAGE_DIR': 'images',
                'UW_DEFAULT_WORLD_DIR': 'worlds', 'UW_DEFAULT_UWS_DIR': 'worlds'}
# Tracked configurations whose paths pin a run (finish never archives it).
PINNING_CONFIGS = ('configs/sandbox', 'configs/evaluation')
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
    if match:
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
    """The storage root named by the main tree's runs/storage.json (a worktree's copy is not read)."""
    store = storage.storage_root(main_root() / 'runs/storage.json')
    if store is None:
        raise SystemExit('runs/storage.json names no storage root: named directories are written to the storage drive')
    return store


def _storage_parent(parent, root, store):
    """parent created inside the routed out/scratch, out/evaluation or runs/<module> of the main tree."""
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
    # The registries change with every named run; they do not make the code dirty.
    diff = git('diff', 'HEAD', '--quiet', '--', '.', ':!runs/registry')
    return {'revision': revision.stdout.strip() if revision.returncode == 0 else None, 'dirty': diff.returncode != 0}


def _tool():
    return Path(sys.argv[0]).stem if sys.argv and sys.argv[0] else 'python'


def _numbered(directory):
    if not directory.is_dir():
        return []
    return [int(match.group(1)) for entry in directory.iterdir()
            if entry.is_dir() and (match := NUMBERED.match(entry.name))]


def next_id(root):
    """One more than the scratch registry's highest id and every numeric prefix under the kinds' parents
    (under out/evaluation: the folders of the modules configs/evaluation/*.json names)."""
    used = [entry['id'] for entry in runs.registry(REGISTRY, root=root)['runs']]
    for kind in ('sandbox', 'step-check', 'geology', 'page'):
        used += _numbered(_parent(kind, None, root))
    for profile in sorted((root / 'configs/evaluation').glob('*.json')):
        used += _numbered(_parent('experiment', profile.stem, root))
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
    root = main_root()
    if kind == 'run':
        # runs/<module>/<id> of the main tree, numbered by its registry (storage.run_directory, from any checkout).
        parent = _storage_parent(root / 'runs' / module, root, store)
        identifier = runs.next_id(module, root=root)
        path = parent / str(identifier)
        path.mkdir(exist_ok=False)
        return Named(identifier, kind, tag, module, str(identifier), path)
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
    root = main_root()
    if kind == 'run':
        if not ref.isdigit():
            raise FileNotFoundError(f'A registered {module} run is named by its id, not {ref!r}')
        entry = next((entry for entry in runs.registry(module, root=root)['runs'] if entry['id'] == int(ref)), None)
        if entry is None:
            raise FileNotFoundError(f'{module} run {ref} is not in the main tree\'s runs/registry/{module}.json')
        if not entry.get('run'):
            raise FileNotFoundError(f'{module} run {ref} is sealed under runs/archive')
        return root / entry['run']
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


def earlier(named):
    """The earlier records of named's topic (kind, module, tag), newest first."""
    if named.kind == 'run':
        found = [entry for entry in runs.registry(named.module, root=main_root())['runs']
                 if entry.get('tag') == named.tag and entry['id'] < named.id]
    else:
        found = [entry for entry in records(named.kind)
                 if entry.get('module') == named.module and entry.get('tag') == named.tag and entry['id'] < named.id]
    return sorted(found, key=lambda entry: entry['id'], reverse=True)


def previous(named):
    """The newest earlier record of named's topic (kind, module, tag), or None."""
    return next(iter(earlier(named)), None)


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def tool_environment(directory, base=None):
    """base (default os.environ) with the UW_*_DIR variables pointing below directory (TOOL_OUTPUTS), so a headless
    tool's own outputs land there; SystemExit when directory does not resolve under the storage root."""
    directory = Path(os.path.abspath(directory))
    store = _storage()
    if not _under(directory, store):
        raise SystemExit(f'{directory} does not resolve under the storage root {store}; a tool writes its own '
                         'outputs to the storage drive')
    environment = dict(os.environ if base is None else base)
    environment.update({name: str(directory / folder) for name, folder in TOOL_OUTPUTS.items()})
    return environment


def page_directory(named):
    """out/evaluation/<name> of the main tree (created when missing): the page folder of a scratch run, named like
    the run (same id and tag; no record of its own). Its pages link the run's images by relative path."""
    root = main_root()
    parent = _storage_parent(root / PARENTS['page'], root, _storage())
    path = parent / named.name
    path.mkdir(exist_ok=True)
    return path


def _forms(path):
    """A path as given (absolute) and resolved through links, case-folded."""
    path = Path(path)
    forms = {Path(os.path.normcase(os.path.abspath(path)))}
    try:
        forms.add(Path(os.path.normcase(os.path.realpath(path))))
    except (OSError, ValueError):
        pass
    return forms


def _inside(path, folder):
    """True when path is folder or lies inside it, compared as given and as resolved."""
    return any(inner == outer or outer in inner.parents for inner in _forms(path) for outer in _forms(folder))


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)


# A `--set`-style override string, key=value (configs/sandbox/presets.json: base_bundle="out/scratch/sandbox/base-v2").
SETTING = re.compile(r'[A-Za-z_][\w.\-]*=(.+)')


def _path_texts(text):
    """The path-like texts a configuration string names: the string itself, and the value of a key=value setting
    (a JSON string value unquoted)."""
    texts = [text]
    match = SETTING.fullmatch(text.strip())
    if match:
        value = match.group(1).strip()
        try:
            parsed = json.loads(value)
        except ValueError:
            parsed = value.strip('"\'')
        if isinstance(parsed, str):
            texts.append(parsed)
    return [item for item in texts if ('/' in item or '\\' in item) and '\n' not in item and len(item) < 1000]


def _config_paths(root):
    """(config file, path) of every path-like string in the main tree's and this checkout's PINNING_CONFIGS (also the
    value of a key="path" setting); a relative path is read against the main tree, where the runs are."""
    found = []
    for tree in dict.fromkeys([root, paths.ROOT]):
        for folder in PINNING_CONFIGS:
            for file in sorted((tree / folder).glob('*.json')):
                try:
                    document = json.loads(file.read_text(encoding='utf-8-sig'))
                except (OSError, ValueError):
                    continue
                for text in (item for string in _strings(document) for item in _path_texts(string)):
                    try:
                        found.append((file.relative_to(tree).as_posix(), Path(text) if Path(text).is_absolute()
                                      else root / text))
                    except (OSError, ValueError):
                        continue
    return found


def _review_open(folder):
    """True when an evaluation folder awaits its review: its evaluation.json's verdict.review is 'pending' or it holds
    a review-request.json, and its review.json holds no entry with a verdict. An unreadable file counts as open."""
    evaluation, request, review = folder / 'evaluation.json', folder / 'review-request.json', folder / 'review.json'
    asked = request.is_file()
    if evaluation.is_file():
        try:
            document = json.loads(evaluation.read_text(encoding='utf-8'))
        except (OSError, ValueError):
            return True
        verdict = document.get('verdict') if isinstance(document, dict) else None
        asked = asked or (isinstance(verdict, dict) and verdict.get('review') == 'pending')
    if not asked:
        return False
    if not review.is_file():
        return True
    try:
        entries = json.loads(review.read_text(encoding='utf-8')).get('entries') or []
    except (OSError, ValueError, AttributeError):
        return True
    return not any(isinstance(entry, dict) and entry.get('verdict') for entry in entries)


def awaiting_review(path):
    """The first evaluation folder of a run (its root or an evaluation* child) that awaits its review, or None."""
    folders = [path] + (sorted(item for item in path.glob('evaluation*') if item.is_dir()) if path.is_dir() else [])
    return next((folder for folder in folders if _review_open(folder)), None)


def pinned(named, record, path, inputs=(), config_paths=None):
    """Why an earlier run (record, at path) of named's topic must stay as it is, or None: the module registry's
    baseline, an evaluation awaiting its review, a path named in configs/sandbox/*.json or configs/evaluation/*.json
    (config_paths: _config_paths, read when not given), or an input of the new run."""
    if named.kind == 'run':
        baseline = runs.registry(named.module, root=main_root()).get('baseline')
        if baseline is not None and str(baseline) == str(record['id']):
            return f'the baseline of runs/registry/{named.module}.json'
    waiting = awaiting_review(path)
    if waiting is not None:
        return f'its evaluation awaits a review ({waiting.name})'
    for file, named_path in _config_paths(main_root()) if config_paths is None else config_paths:
        if _inside(named_path, path):
            return f'named in {file}'
    for item in inputs:
        if item and _inside(Path(item), path):
            return f'an input of {named.name} ({item})'
    return None


def _mark_superseded(record, named):
    """superseded_by (the new run's name) and superseded_utc on record in the main tree's scratch.json, under the
    naming lock; every other byte of meaning stays (the document is checked after the write)."""
    root = main_root()
    path = runs.registry_path(REGISTRY, root)
    with _Lock(root):
        data = path.read_bytes()
        document = json.loads(data.decode('utf-8'))
        before = copy.deepcopy(document)
        indices = [index for index, entry in enumerate(document['runs'])
                   if entry.get('id') == record['id'] and entry.get('run') == record.get('run')]
        if len(indices) != 1:
            raise RuntimeError(f'{path}: {len(indices)} records of id {record["id"]} at {record.get("run")}')
        entry = document['runs'][indices[0]]
        if entry.get('superseded_by'):
            return entry
        entry.update(superseded_by=named.name, superseded_utc=_now())
        newline = '\r\n' if b'\r\n' in data else '\n'
        text = json.dumps(document, indent=2, ensure_ascii=False) + '\n'
        path.write_bytes(text.replace('\n', newline).encode('utf-8'))
        written = json.loads(path.read_bytes().decode('utf-8'))
        restored = copy.deepcopy(written)
        for key in ('superseded_by', 'superseded_utc'):
            restored['runs'][indices[0]].pop(key, None)
        if written != document or restored != before:
            raise RuntimeError(f'{path}: the record of {record["run"]} did not read back as written; check the registry')
    return entry


# The file a tool writes when a run of a kind has ended (written last, or in a finally block): finish leaves an
# earlier run without it alone, so a run still in progress is never touched.
COMPLETE = {'run': 'ledger-entry.json', 'experiment': 'ledger-entry.json', 'step-check': 'result.json', 'sandbox': 'config.json',
            'geology': 'report.json', 'page': 'index.html'}
# Outcomes finish prints one line each for; the others (repeated at every run of a topic) only as counts.
SHOWN = ('archived', 'superseded', 'pinned', 'refused', 'error')


def finish(named, passed, inputs=()):
    """The automatic move after a run of named's topic ends (module docstring, "Automatic moves"): every earlier run
    of the topic that is not pinned, archived or superseded yet is archived or marked, newest first, so a failed run
    in between leaves nothing behind.

    Returns None when the run did not pass, else a list (newest first) of {'previous', 'id', 'action' (archived,
    superseded, pinned, already-archived, already-superseded, missing, not-a-run, refused or error), 'reason'}.
    Never raises because of an earlier run: a refusal or an error is printed and that run stays as it is."""
    if not passed:
        return None
    try:
        found = [record for record in earlier(named) if record.get('run')]
        config_paths = _config_paths(main_root()) if found else []
    except Exception as error:  # an unreadable registry or configuration: the finished run stands
        print(f'finish {named.name}: earlier runs unknown ({type(error).__name__}: {error})', flush=True)
        return [{'previous': None, 'id': None, 'action': 'error', 'reason': f'{type(error).__name__}: {error}'}]
    outcomes = [_finish_one(named, record, inputs, config_paths) for record in found]
    for outcome in outcomes:
        if outcome['action'] in SHOWN:
            print(f'finish {named.name}: earlier {outcome["previous"]} {outcome["action"]}'
                  + (f' ({outcome["reason"]})' if outcome['reason'] else ''), flush=True)
    others = {}
    for outcome in outcomes:
        if outcome['action'] not in SHOWN:
            others[outcome['action']] = others.get(outcome['action'], 0) + 1
    if others:
        print(f'finish {named.name}: ' + ', '.join(f'{count} {action}' for action, count in sorted(others.items())),
              flush=True)
    return outcomes


def _step_check_setting(path):
    """(grid, social) of a step check's result.json, or None when it cannot be read."""
    try:
        result = json.loads((path / 'result.json').read_text(encoding='utf-8'))
        return result.get('grid'), bool(result.get('social'))
    except (OSError, ValueError, AttributeError):
        return None


def _other_setting(named, path):
    """Why an earlier step check is another topic in effect (its result.json names another grid or --social setting
    than the new check's, or either cannot be read), or ''; '' for the other kinds."""
    if named.kind != 'step-check':
        return ''
    new, old = _step_check_setting(named.path), _step_check_setting(path)
    if new is None or old is None:
        return 'a result.json cannot be read'
    if new != old:
        return f'grid {old[0]}{", social" if old[1] else ""}, not {new[0]}{", social" if new[1] else ""}'
    return ''


def _finish_one(named, record, inputs, config_paths):
    """finish for one earlier record: its outcome; never raises."""
    path = main_root() / record['run']
    outcome = {'previous': record['run'], 'id': record['id'], 'action': None, 'reason': ''}
    try:
        if record.get('superseded_by'):
            outcome.update(action='already-superseded', reason=f'by {record["superseded_by"]}')
        elif not path.is_dir():
            outcome.update(action='missing', reason='the folder is gone')
        elif not (path / COMPLETE[named.kind]).is_file():
            outcome.update(action='incomplete', reason=f'no {COMPLETE[named.kind]}: running, or stopped before its end')
        elif os.path.lexists(storage.long_path(path) / archive.RECEIPT):
            outcome.update(action='already-archived', reason=archive.RECEIPT)
        elif (setting := _other_setting(named, path)):
            outcome.update(action='other-setting', reason=setting)
        elif (reason := pinned(named, record, path, inputs, config_paths)):
            outcome.update(action='pinned', reason=reason)
        else:
            try:
                receipt = archive.archive(path, f'superseded by {named.name}', apply=True)
                outcome.update(action='archived', reason=f'{len(receipt["dropped"])} files, '
                                                         f'{receipt["dropped_bytes"] / 1e6:.1f} MB dropped')
            except archive.NotARun:
                if named.kind == 'run':  # a module registry is append-only: nothing to mark
                    outcome.update(action='not-a-run', reason='no evaluation folder; left as it is')
                else:
                    _mark_superseded(record, named)
                    outcome.update(action='superseded', reason=f'superseded_by {named.name} in {REGISTRY}.json')
            except archive.Refused as error:
                outcome.update(action='refused', reason=str(error))
    except Exception as error:  # the finished run stands; an earlier run's problem is only reported
        outcome.update(action='error', reason=f'{type(error).__name__}: {error}')
    return outcome


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
