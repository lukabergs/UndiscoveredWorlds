"""Reject local-only data and large blobs from the index or reachable history."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
LOCAL_ROOTS = {'meta', 'docs', 'refs', 'runs', 'extra', 'references', 'data',
               'out', 'build', 'bin', 'obj', 'saved_worlds', '0'}
LOCAL_FILES = {'climate_benchmark_runs.json', 'climate_reference_data.json', 'UWinstaller.zip'}
DATA_SUFFIXES = {'.zip', '.7z', '.tar', '.gz', '.bz2', '.xz', '.uwclim', '.nc',
                 '.hdf', '.h5', '.tif', '.tiff', '.npz', '.npy', '.xlsx', '.xls',
                 '.exe', '.dll', '.lib', '.so', '.a', '.dylib', '.pdb'}
MAX_BYTES = 5 * 1024 * 1024


def path_issue(path):
    parts = path.replace('\\', '/').split('/')
    if parts[0].lower() in LOCAL_ROOTS or path in LOCAL_FILES:
        return 'local-only directory or artifact'
    if path.startswith('tools/ripgrep/'):
        return 'downloaded tool distribution'
    fixture = parts[0] == 'tests' and 'fixtures' in parts
    if not fixture and Path(path).suffix.lower() in DATA_SUFFIXES:
        return 'downloaded binary or generated data'
    return None


def git(root, *args, input=None):
    return subprocess.check_output(['git', *args], cwd=root, input=input)


def inspect(root, revisions=None):
    issues = []
    if revisions:
        names = git(root, 'log', '--full-history', '-m', '--root', '--format=',
                    '--name-only', '-z', *revisions).decode('utf-8').split('\0')
        names = {name.strip('\n') for name in names if name.strip('\n')}
        objects = git(root, 'rev-list', '--objects', '--no-object-names', *revisions).splitlines()
        labels = {}
    else:
        entries = git(root, 'ls-files', '--stage', '-z').split(b'\0')
        labels = {}
        names = set()
        for entry in entries:
            if not entry:
                continue
            meta, name = entry.split(b'\t', 1)
            mode, oid, stage = meta.split()
            name = name.decode('utf-8'); names.add(name)
            if stage != b'0':
                issues.append(f'{name}: unresolved index stage {stage.decode()}')
            if mode != b'160000':
                labels[oid] = name
        objects = list(labels)
    for name in sorted(names):
        reason = path_issue(name)
        if reason:
            issues.append(f'{name}: {reason}')
    if objects:
        rows = git(root, 'cat-file', '--batch-check=%(objectname) %(objecttype) %(objectsize)',
                   input=b'\n'.join(objects) + b'\n').splitlines()
        for row in rows:
            oid, kind, size = row.split()
            if kind == b'blob' and int(size) > MAX_BYTES:
                issues.append(f'{labels.get(oid, oid.decode())}: {int(size):,} bytes exceeds {MAX_BYTES:,}')
    return issues, len(names)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--staged', action='store_true', help='Inspect the entire index')
    group.add_argument('--history', nargs='+', metavar='REV', help='Inspect all reachable history')
    group.add_argument('--pre-push', action='store_true', help='Read Git pre-push ref records from stdin')
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args(argv)
    revisions = args.history
    if args.pre_push:
        revisions = sorted({line.split()[1] for line in sys.stdin if line.strip()
                            and set(line.split()[1]) != {'0'}})
        if not revisions:
            return 0
    issues, count = inspect(args.root, revisions)
    if issues:
        print('Repository scope check failed:\n' + '\n'.join(issues[:40]), file=sys.stderr)
        if len(issues) > 40:
            print(f'... and {len(issues)-40} further violations', file=sys.stderr)
        return 1
    print(f'Repository scope passed: {count} paths; no local-only data or blobs over 5 MiB.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
