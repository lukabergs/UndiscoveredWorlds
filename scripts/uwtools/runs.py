"""Where run N of module M is: one resolver for every reader of registered runs.

A run registered since 2026-10-05 lives in the run directory its registry record names
(`"run": "runs/<module>/<id>"`: run.json, log.txt, bundle/, checkpoints/, evaluation/).
Runs registered before are sealed under runs/archive/ (runs/archive/README.md): a recorded
path `runs/<kind>/...` is now `runs/archive/<kind>/...`. Read sealed runs, never write
there: derived outputs of a run go to output_directory(), the run's evaluation/analysis/
or, for a sealed run, runs/reports/<module>/<id>/analysis/.
"""
import json
import os
from pathlib import Path

from .paths import ROOT
from .tectonics import run_snapshots

ARCHIVE = ROOT / 'runs/archive'
SEALED_KINDS = ('checkpoints', 'diagnostics', 'fields', 'logs', 'maps', 'manifests', 'parameters', 'metrics')


def registry_path(module, root=None):
    """runs/registry/<module>.json of root (default this checkout; naming.py passes the main tree)."""
    return (root or ROOT) / 'runs/registry' / f'{module}.json'


def registry(module, root=None):
    path = registry_path(module, root)
    return json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {'runs': []}


def record(module, identifier):
    for entry in registry(module)['runs']:
        if entry['id'] == int(identifier):
            return entry
    raise KeyError(f'{module} run {identifier} is not in {registry_path(module).relative_to(ROOT).as_posix()}')


def directory(module, identifier):
    """runs/<module>/<id> of a run in the current layout; None for a sealed run."""
    path = record(module, identifier).get('run')
    return ROOT / path if path else None


def archived(path):
    """A repository path recorded before the seal, as it is now."""
    parts = Path(path).as_posix().split('/')
    if len(parts) > 1 and parts[0] == 'runs' and parts[1] in SEALED_KINDS:
        return ARCHIVE.joinpath(*parts[1:])
    return ROOT / path


def output_directory(module, identifier):
    """Where scripts write what they derive from one run (not created here)."""
    run = directory(module, identifier)
    return run / 'evaluation/analysis' if run else ROOT / 'runs/reports' / module / str(identifier) / 'analysis'


def derived(module, identifier, name):
    """An existing derived file of a run: its output directory first, then a sealed run's diagnostics."""
    candidates = [output_directory(module, identifier) / name]
    if directory(module, identifier) is None:
        candidates.append(ARCHIVE / 'diagnostics' / module / str(identifier) / name)
    return next((path for path in candidates if path.is_file()), None)


def next_id(module):
    """One more than every registered id and every run directory of the module."""
    used = [entry['id'] for entry in registry(module)['runs']]
    parent = ROOT / 'runs' / module
    used += [int(path.name) for path in parent.iterdir() if path.name.isdigit()] if parent.is_dir() else []
    return max(used, default=0) + 1


def refuse_archive(path):
    """path as an absolute path; SystemExit when it lies in the sealed runs/archive (by its
    workspace path or through the junction's target), which nothing may write into."""
    path = Path(os.path.abspath(path))
    for candidate in (path, path.resolve()):
        for sealed in (ARCHIVE, ARCHIVE.resolve()):
            if candidate == sealed or sealed in candidate.parents:
                raise SystemExit(f'{path} is in the sealed runs/archive; write derived outputs elsewhere '
                                 '(runs/reports/<module>/<id>/analysis or out/)')
    return path


def append_record(module, entry, root=None):
    """Append one record to runs/registry/<module>.json (of root, default this checkout) as text: every
    existing byte stays as it is.

    A missing registry is created as {"schema": "<module>-run-registry/v1", "runs": []}. One
    session registers runs at a time: concurrent appends are not supported (naming.py holds a lock
    around its appends to scratch.json), and a registry changed between reading and writing, or a
    record missing after the write, raises."""
    path = registry_path(module, root)
    if not path.is_file():
        path.write_bytes(json.dumps({'schema': f'{module}-run-registry/v1', 'runs': []}, indent=2).encode('utf-8') + b'\n')
    data = path.read_bytes()
    newline = b'\r\n' if b'\r\n' in data else b'\n'
    end = data.rstrip().rfind(b']')
    if end < 0:
        raise ValueError(f'{path}: no runs list')
    body = data[:end].rstrip()
    text = json.dumps(entry, indent=2, ensure_ascii=False).encode('utf-8')
    lines = [b'    ' + line for line in text.split(b'\n')]
    separator = b'' if body.endswith(b'[') else b','
    updated = body + separator + newline + newline.join(lines) + newline + b'  ' + data[end:]
    document = json.loads(updated.decode('utf-8'))
    if document['runs'][:-1] != json.loads(data.decode('utf-8'))['runs'] or document['runs'][-1] != entry:
        raise ValueError(f'{path}: appending would change existing records')
    if path.read_bytes() != data:
        raise RuntimeError(f'{path} changed while appending run {entry["id"]}; register one run at a time')
    path.write_bytes(updated)
    if path.read_bytes() != updated or registry(module, root)['runs'][-1] != entry:
        raise RuntimeError(f'{path}: run {entry["id"]} is not the last record after writing; check the registry')


def tectonics(identifier):
    """A tectonics run described as its former manifest: id, status, recipe, elapsed_seconds,
    snapshots (label, update_count, time_myr, directory; path = the bundle), field_sha256 (sealed
    runs), replay, and for feedback runs their own keys. Snapshot paths are resolved for either layout."""
    identifier = int(identifier)
    run = directory('tectonics', identifier)
    if run is None:
        manifest = json.loads((ARCHIVE / f'manifests/tectonics/{identifier}.json').read_text(encoding='utf-8'))
        for snapshot in manifest.get('snapshots', []):
            snapshot['path'] = ARCHIVE / f'checkpoints/tectonics/{identifier}' / snapshot['directory']
        return manifest
    document = json.loads((run / 'run.json').read_text(encoding='utf-8'))
    snapshots = run_snapshots(run)
    for snapshot in snapshots:
        snapshot['path'] = run / snapshot['directory']
    repeat = document.get('repeat')
    # The recipe file's scenario as run, with the run's own seed and grid.
    recipe = dict(document['recipe'], seed=document['seed'], width=document['grid'][0], height=document['grid'][1])
    return {'id': identifier, 'status': record('tectonics', identifier).get('status'), 'recipe': recipe,
            'elapsed_seconds': document.get('seconds'), 'snapshots': snapshots, 'run': run,
            'replay': {'exact_match': repeat['identical']} if repeat else {'status': 'skipped'}}
