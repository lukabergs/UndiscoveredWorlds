"""Archive a finished run: drop the files the scripts re-render, keep everything else (master plan 2.3).

    uv run python scripts/uwtools/archive.py plan PATH...                 (dry run: what would be dropped)
    uv run python scripts/uwtools/archive.py apply PATH... --reason TEXT
    uv run python scripts/uwtools/archive.py list                         (the archive receipts)

A path names one run-level folder: one with an evaluation.json at its root or in an evaluation* child folder (else
NotARun, a Refused: a sandbox, step-check or geology folder), and no evaluation.json anywhere else except nested in
an evaluation* folder (a comparison, kept whole); another would be several runs at once. Its evaluation folders are
the root itself and
its evaluation* children, each holding a uw-evaluation/v1 evaluation.json; in them two classes are dropped, nothing
else:
  R1  map PNGs under the evaluation folder's maps/ that its evaluation.json lists in maps[].path, only while the
      evaluated bundle's manifest.json exists and has the digest evaluation.json recorded (run.manifest_sha256);
      a kind 'parent' map only while the baseline bundle's manifest matches baseline.manifest_sha256
      (regenerate: uw.py <module> evaluate --bundle B --output <new folder>). The bundle may lie outside the
      archived folder (run.bundle, accepted 2026-10-21): trashing it later leaves those PNGs without a source,
      which the trash's restore covers. No PNG is dropped from a folder with a kind 'animation' map (its
      checkpoints are not recorded) or with a review*.json holding a verdict (a judged run's images are evidence);
  R2  an evaluation folder's index.html and report.md beside its evaluation.json (regenerate: report.py run DIR).
Each dropped file's receipt entry holds its regenerate command and the code {revision, dirty} evaluation.json
recorded.
Kept: .npy, bundles, checkpoints, run.json, evaluation.json, review and metrics JSON, maps/palette.json, CSV, logs,
any PNG evaluation.json does not list (sandbox images, analysis maps) and every PNG of a run without its bundle.

Files are deleted in place, never moved to the trash (the trash is deny-delete: moving there frees nothing).
Refused, with nothing dropped: a path that is or holds a link, lies in runs/archive, is or holds a protected place
or a junction target of the main tree (or lies in a sealed place), holds tracked files (checked in the main tree,
the path mapped back through its junctions), is not one run-level folder, or does not resolve under the storage
root (nothing on C is dropped). The receipt <dir>/archive-receipt.json (status 'dropping', every file with its bytes
and SHA-256) is written before the first file goes, with a copy in <storage root>/receipts/archive/; naming the
same path again finishes a stopped job. Every file is re-hashed just before its removal; a changed file stops the
job. With several paths, every path is checked (receipt states included) before the first file goes.

archive(path, reason, apply) -> receipt is the interface naming.finish (master plan 2.4) calls.
"""
import argparse
import datetime
import json
import os
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import runs, storage, trash  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402
from uwtools.paths import ROOT, main_root  # noqa: E402

SCHEMA = 'uw-archive-receipt/v1'
RULES = 'uw-archive-rules/v2'
RECEIPT = 'archive-receipt.json'
EVALUATION = re.compile(r'evaluation.*')
EVALUATION_SCHEMA = 'uw-evaluation/v1'
PAGES = ('index.html', 'report.md')


class Refused(OSError):
    """The path may not be archived (or the job must stop); nothing more is dropped."""


class NotARun(Refused):
    """The folder holds no evaluation folder (a sandbox, step-check or geology folder): nothing is regenerable."""


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def real_path(path):
    return Path(os.path.realpath(path))


def storage_root():
    """The resolved storage root, or Refused (none configured, an unreadable runs/storage.json, the drive missing)."""
    try:
        root = storage.storage_root()
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise Refused(f'The storage root is unavailable ({error}); nothing on C is archived') from None
    if root is None:
        raise Refused('No storage configured (runs/storage.json); nothing on C is archived')
    return real_path(root)


def main_tree():
    try:
        return real_path(main_root())
    except (OSError, subprocess.CalledProcessError) as error:
        raise Refused(f'The main tree is unknown ({error})') from None


def junctions(tree):
    """(junction, resolved target) of every routed folder of a tree."""
    found = []
    for parent in (tree / 'runs', tree / 'runs' / 'reports', tree / 'runs' / 'work', tree / 'out'):
        for item in parent.iterdir() if parent.is_dir() else []:
            target = storage.junction_target(item)
            if target is not None:
                found.append((item, real_path(target)))
    return found


def protected_places(trees, root):
    """trash.protected_places, built from the main tree (and the tree this runs in), plus every junction target."""
    places = [root, root / 'worktree-archive', root / 'receipts', root / 'trash', root / 'refs', root / 'out',
              root / 'out' / 'scratch', root / 'out' / 'evaluation', root / 'runs']
    for tree in trees:
        places += [tree, tree / '.git', tree / 'refs', tree / 'meta', tree / 'runs' / 'registry',
                   tree / 'runs' / 'ledger', tree / 'runs' / 'archive', tree / 'out', tree / 'runs']
        places += [target for _, target in junctions(tree)]
    return [real_path(place) for place in places]


def tracked_files(tree, relative):
    """Tracked files of tree at or below relative (posix); lists the top folder and filters, so git never has to
    follow a junction in the path."""
    top = relative.split('/')[0]
    result = subprocess.run(['git', '-C', str(tree), 'ls-files', '-z', '--', top], capture_output=True, text=True,
                            encoding='utf-8', errors='replace')
    if result.returncode != 0:
        raise Refused(f'{tree / relative}: git ls-files failed; nothing dropped')
    folded = relative.lower().rstrip('/')
    return [name for name in result.stdout.split('\0')
            if name and (name.lower() == folded or name.lower().startswith(folded + '/'))]


def logical_paths(given, real, main, root):
    """(tree, path inside it) of every repository path that may name real: the given path when it lies in this tree,
    and real mapped back through the main tree's junctions (else by its path below the storage root)."""
    found = []
    if trash.inside(given, ROOT):
        found.append((ROOT, given))
    mapped = [(main, junction / real.relative_to(target)) for junction, target in junctions(main)
              if trash.inside(real, target)]
    found += mapped or [(main, main / real.relative_to(root))]
    return found


def evaluation_folders(base):
    """The evaluation folders of a run-level folder (long paths), or Refused when it is not one."""
    folders = [base] if (base / 'evaluation.json').is_file() else []
    folders += sorted(item for item in base.iterdir()
                      if item.is_dir() and EVALUATION.fullmatch(item.name) and (item / 'evaluation.json').is_file())
    if not folders:
        raise NotARun(f'{base}: no evaluation.json at its root or in an evaluation* folder; not a run-level folder')
    # An evaluation.json nested in an evaluation* folder (a comparison, runs/erosion/9/evaluation/comparison-8) is
    # part of the run; such folders are kept whole (the rules look only at the evaluation folders' own files).
    named = [folder for folder in folders if EVALUATION.fullmatch(folder.name)]
    allowed = {folder / 'evaluation.json' for folder in folders}
    stray = [item for item in base.rglob('evaluation.json')
             if item not in allowed and not any(folder in item.parents for folder in named)]
    if stray:
        raise Refused(f'{base}: {len(stray)} evaluation.json in subfolders (e.g. {stray[0]}): several runs at once; '
                      'name each run-level folder')
    return folders


def check(path):
    """(the given absolute path, its resolved path, the storage root), or Refused saying why it may not be archived."""
    given = Path(os.path.abspath(path))
    if not os.path.lexists(given):
        raise Refused(f'{given}: does not exist')
    if storage.is_link(given):
        raise Refused(f'{given}: is a link or junction; links are never archived')
    try:
        runs.refuse_archive(given)
    except SystemExit as error:
        raise Refused(str(error)) from None
    real = real_path(given)
    if not real.is_dir():
        raise Refused(f'{real}: not a folder')
    root = storage_root()
    if not trash.inside(real, root) or real == root:
        raise Refused(f'{real}: does not resolve under the storage root {root}; nothing on C is archived')
    main = main_tree()
    for place in protected_places(dict.fromkeys([main, real_path(ROOT)]), root):
        if trash.inside(place, real) or (trash.inside(real, place) and place.name in trash.SEALED):
            raise Refused(f'{real}: is, holds or lies in the protected place {place}')
    links = storage.find_links(real)
    if links:
        raise Refused(f'{real}: {len(links)} links or unreadable folders inside, e.g. {links[0]}')
    evaluation_folders(storage.long_path(real))
    for tree, logical in logical_paths(given, real, main, root):
        tracked = tracked_files(tree, logical.relative_to(tree).as_posix())
        if tracked:
            raise Refused(f'{logical}: holds tracked files (e.g. {tracked[0]})')
    return given, real, root


def shown(path):
    """A path as a command line names it: repository-relative when inside the repository."""
    path = Path(path)
    return path.relative_to(ROOT).as_posix() if trash.inside(path, ROOT) else path.as_posix()


def resolve_bundle(text):
    """A bundle path as evaluation.json records it: absolute, or relative to the repository (the main tree's, when
    run from a worktree that lacks runs/)."""
    path = Path(text)
    if path.is_absolute() or (ROOT / path).exists():
        return ROOT / path
    try:
        return main_root() / path
    except (OSError, subprocess.CalledProcessError):
        return ROOT / path


def bundle_matches(record):
    """True when record names a bundle whose manifest.json exists and has the digest record holds."""
    if not isinstance(record, dict) or not record.get('bundle') or not record.get('manifest_sha256'):
        return False
    manifest = resolve_bundle(str(record['bundle'])) / 'manifest.json'
    return manifest.is_file() and sha256(manifest) == record['manifest_sha256']


def judged(folder):
    """True when a review*.json of the folder holds a verdict (or cannot be read: then it may)."""
    for file in folder.glob('review*.json'):
        try:
            document = json.loads(file.read_text(encoding='utf-8'))
        except (OSError, ValueError):
            return True
        if not isinstance(document, dict):
            continue
        entries = document.get('entries') if isinstance(document.get('entries'), list) else []
        if document.get('verdict') or any(isinstance(entry, dict) and entry.get('verdict') for entry in entries):
            return True
    return False


def code_of(document):
    """The code version evaluation.json recorded: {revision, dirty} (None where it recorded nothing)."""
    code = document.get('code') if isinstance(document.get('code'), dict) else {}
    return {'revision': code.get('revision'), 'dirty': code.get('dirty')}


def evaluate_command(document, output_parent):
    module = (document.get('profile') or {}).get('module', '<module>')
    command = (f"uv run --offline python scripts/uw.py {module} evaluate --bundle {document['run']['bundle']} "
               f"--output {output_parent}/<new folder>")
    for key, flag in (('baseline', '--baseline'), ('repeat', '--repeat'), ('input', '--input')):
        if isinstance(document.get(key), dict) and document[key].get('bundle'):
            command += f" {flag} {document[key]['bundle']}"
    if isinstance(document.get('earth'), dict) and document['earth'].get('entry'):
        command += f" --earth {document['earth']['entry']}"
    return command


def map_paths(document, folder):
    """{relative path inside folder: kind} of the PNG maps the evaluation's R1 rule drops."""
    entries = [entry for entry in document.get('maps') or [] if isinstance(entry, dict)]
    if any(entry.get('kind') == 'animation' for entry in entries) or judged(folder):
        return {}
    if not bundle_matches(document.get('run')):
        return {}
    parent = bundle_matches(document.get('baseline'))
    found = {}
    for entry in entries:
        text = entry.get('path')
        if not isinstance(text, str) or (entry.get('kind') == 'parent' and not parent):
            continue
        parts = PurePosixPath(text.replace('\\', '/')).parts
        if (len(parts) < 2 or parts[0] != 'maps' or '..' in parts or ':' in text
                or not parts[-1].lower().endswith('.png')):
            continue
        file = folder.joinpath(*parts)
        if file.is_file() and not storage.is_link(file):
            found['/'.join(parts)] = entry.get('kind')
    return found


def regenerable(given, real):
    """{path relative to real (posix): (class, regenerate command, code)} of the files the two rules drop."""
    found = {}
    base = storage.long_path(real)
    for folder in evaluation_folders(base):
        try:
            document = json.loads((folder / 'evaluation.json').read_text(encoding='utf-8'))
        except (OSError, ValueError):
            continue  # unreadable: nothing in it is known to be regenerable
        if not isinstance(document, dict) or document.get('schema') != EVALUATION_SCHEMA:
            continue
        relative = folder.relative_to(base)
        logical = given / relative
        code = code_of(document)
        page = f'uv run --offline python scripts/evaluation/report.py run {shown(logical)}'
        for name in PAGES:
            if (folder / name).is_file():
                found[(relative / name).as_posix()] = ('R2', page, code)
        maps = map_paths(document, folder)
        if maps:
            command = evaluate_command(document, shown(logical.parent))
            for text in maps:
                found[(relative / text).as_posix()] = ('R1', command, code)
    return found


def plan(given, real):
    """The receipt the archive would write, with nothing written or dropped."""
    base = storage.long_path(real)
    rules = regenerable(given, real)
    dropped, kept_files, kept_bytes = [], 0, 0
    for file in sorted(item for item in base.rglob('*') if item.is_file()):
        key = file.relative_to(base).as_posix()
        size = file.stat().st_size
        if key in rules:
            dropped.append({'path': key, 'bytes': size, 'sha256': sha256(file), 'class': rules[key][0],
                            'regenerate': rules[key][1], 'code': rules[key][2]})
        elif key != RECEIPT:
            kept_files += 1
            kept_bytes += size
    return {'schema': SCHEMA, 'path': str(given), 'real': str(real), 'utc': now(), 'reason': '', 'rules': RULES,
            'status': 'plan', 'dropped': dropped, 'kept_files': kept_files, 'kept_bytes': kept_bytes,
            'dropped_bytes': sum(item['bytes'] for item in dropped)}


def copy_path(real, root):
    relative = real.relative_to(root).as_posix()
    return root / 'receipts' / 'archive' / (relative.replace('/', '__') + '.json')


def write_receipt(receipt, real, root):
    """The receipt in the folder and its copy under the storage root, both overwritten in place."""
    text = json.dumps(receipt, indent=1) + '\n'
    (storage.long_path(real) / RECEIPT).write_text(text, encoding='utf-8')
    copy = copy_path(real, root)
    copy.parent.mkdir(parents=True, exist_ok=True)
    copy.write_text(text, encoding='utf-8')


def unlink(path):
    trash.writable_retry(os.unlink, path)


def drop(receipt, real, root):
    """Remove the receipt's files that are still there, each re-hashed first; then mark the receipt complete."""
    base = storage.long_path(real)
    present = [(base / item['path'], item) for item in receipt['dropped'] if os.path.lexists(base / item['path'])]
    for file, item in present:  # every file checked before any goes
        if storage.is_link(file) or sha256(file) != item['sha256']:
            raise Refused(f'{file}: changed since the receipt was written; nothing more dropped')
    for file, item in present:
        if storage.is_link(file) or sha256(file) != item['sha256']:
            raise Refused(f'{file}: changed since the receipt was written; nothing more dropped')
        unlink(file)
    receipt['status'] = 'complete'
    receipt['completed_utc'] = now()
    write_receipt(receipt, real, root)
    return receipt


def read_receipt(real, planned):
    """The folder's existing receipt checked against the plan: None, or the receipt; Refused when it is corrupt, of an
    unknown status, or does not fit what is on disk now."""
    receipt_file = storage.long_path(real) / RECEIPT
    if not os.path.lexists(receipt_file):
        return None
    try:
        existing = json.loads(receipt_file.read_text(encoding='utf-8'))
        status = existing['status']
        items = {item['path']: item['sha256'] for item in existing['dropped']}
        if not all(isinstance(key, str) and isinstance(value, str) for key, value in items.items()):
            raise TypeError('a dropped entry without a text path and digest')
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise Refused(f'{receipt_file}: unreadable or corrupt ({error}); sort it out by hand') from None
    rules = {item['path'] for item in planned['dropped']}
    new = sorted(rules - set(items))
    if status == 'dropping':
        if new:
            raise Refused(f'{real}: a stopped archive job, and {len(new)} regenerable files not in its receipt '
                          f'(e.g. {new[0]}); sort it out by hand')
        base = storage.long_path(real)
        unruled = [key for key in items if os.path.lexists(base / key) and key not in rules]
        if unruled:
            raise Refused(f'{real}: {len(unruled)} files of the stopped job no longer fall under the rules '
                          f'(e.g. {unruled[0]}); sort it out by hand')
    elif status == 'complete':
        if new:
            raise Refused(f'{real}: archived before ({existing.get("completed_utc")}), and {len(new)} regenerable '
                          f'files since (e.g. {new[0]})')
    else:
        raise Refused(f'{receipt_file}: unknown status {status!r}')
    return existing


def prepare(path, apply):
    """Every check of one path, receipt state included, with nothing written: (given, real, root, plan, existing)."""
    given, real, root = check(path)
    planned = plan(given, real)
    return given, real, root, planned, read_receipt(real, planned) if apply else None


def execute(prepared, reason):
    given, real, root, planned, existing = prepared
    if existing is None:
        receipt = dict(planned, reason=reason, status='dropping')
        write_receipt(receipt, real, root)  # before the first file goes
        return drop(receipt, real, root)
    if existing['status'] == 'complete':
        return existing
    existing['resumed_utc'] = now()
    write_receipt(existing, real, root)
    return drop(existing, real, root)


def archive(path: Path, reason: str, apply: bool) -> dict:
    """The receipt of archiving path: a plan (status 'plan', nothing written) unless apply.

    With apply: a stopped job (receipt status 'dropping') is finished; an archive already complete is returned as
    it is when nothing new is regenerable, and refused otherwise.
    """
    prepared = prepare(path, apply)
    return execute(prepared, reason) if apply else prepared[3]


def archive_all(paths, reason):
    """Archive several paths: every path is checked, receipt states included, before the first file goes."""
    prepared = [prepare(path, True) for path in paths]
    reals = [item[1] for item in prepared]
    for first in reals:
        for second in reals:
            if first != second and trash.inside(second, first):
                raise Refused(f'{second} lies inside {first}: name each place once')
    if len(set(reals)) != len(reals):
        raise Refused('a folder is named twice')
    return [execute(item, reason) for item in prepared]


def summary(receipt):
    classes = {}
    for item in receipt['dropped']:
        count = classes.setdefault(item['class'], [0, 0])
        count[0] += 1
        count[1] += item['bytes']
    lines = [f"{receipt['path']}  ({receipt['status']})"]
    for name in ('R1', 'R2'):
        files, size = classes.get(name, (0, 0))
        lines.append(f'  {name} {"map PNGs" if name == "R1" else "pages"}: {files} files, {size / 1e6:.1f} MB')
    lines.append(f"  kept: {receipt['kept_files']} files, {receipt['kept_bytes'] / 1e6:.1f} MB")
    return '\n'.join(lines)


def receipts():
    folder = storage_root() / 'receipts' / 'archive'
    return sorted(folder.glob('*.json')) if folder.is_dir() else []


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    planning = commands.add_parser('plan', help='dry run: per class, the files and bytes that would be dropped')
    planning.add_argument('paths', nargs='+')
    applying = commands.add_parser('apply', help='drop the regenerable files, with a receipt', allow_abbrev=False)
    applying.add_argument('paths', nargs='+')
    applying.add_argument('--reason', required=True)
    commands.add_parser('list', help='the archive receipts')
    args = parser.parse_args()
    try:
        if args.command == 'list':
            total = 0
            for path in receipts():
                receipt = json.loads(path.read_text(encoding='utf-8'))
                total += receipt['dropped_bytes']
                print(f"{receipt['status']:9} {receipt['dropped_bytes'] / 1e6:9.1f} MB  {receipt['path']}  "
                      f"{receipt['reason'][:60]}")
            print(f'total {total / 1e9:.2f} GB dropped')
            return 0
        if args.command == 'plan':
            failed = 0
            for path in args.paths:
                try:
                    print(summary(archive(Path(path), '', apply=False)))
                except Refused as error:
                    failed += 1
                    print(f'{path}: refused: {error}')
            return 1 if failed else 0
        for receipt in archive_all([Path(path) for path in args.paths], args.reason):
            print(summary(receipt))
    except Refused as error:
        print(f'refused: {error}')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
