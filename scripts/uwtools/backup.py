"""Copy the untracked data listed in configs/backup.json to the backup destination.

    uv run --offline python scripts/uwtools/backup.py [--only refs runs] [--dry-run]

Additive by design (erosion plan Phase D, 2026-10-07, after most of refs/processed
was deleted through a worktree junction): a file is copied when the backup lacks it
or its size or modification time differs; nothing is deleted in the destination;
a backup file that would be overwritten is first moved to
<destination>/.versions/<timestamp>/<source>/<path>. Junctions inside a source
(runs/<module>, refs/processed/physical) are read through, with a guard against
cycles. Each run writes <destination>/receipts/<timestamp>.json.
"""
import argparse
import datetime
import json
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'configs' / 'backup.json'


def walk(source):
    """Yield (path, relative) for every file below source, reading through junctions once."""
    seen = set()
    stack = [source]
    while stack:
        directory = stack.pop()
        real = os.path.realpath(directory)
        if real in seen:
            continue
        seen.add(real)
        try:
            entries = list(os.scandir(directory))
        except OSError as error:
            print(f'  skipped unreadable {directory}: {error}', file=sys.stderr)
            continue
        for entry in entries:
            try:
                if entry.is_dir():  # follows junctions and links
                    stack.append(entry.path)
                elif entry.is_file():
                    yield entry.path, os.path.relpath(entry.path, source)
            except OSError as error:
                print(f'  skipped {entry.path}: {error}', file=sys.stderr)


def same(source, target):
    a, b = os.stat(source), os.stat(target)
    return a.st_size == b.st_size and abs(a.st_mtime - b.st_mtime) < 2.0


def back_up(name, source, destination, stamp, dry_run, counts):
    target_root = destination / name
    for path, relative in walk(source):
        target = target_root / relative
        counts['files'] += 1
        if target.exists():
            if same(path, target):
                continue
            kept = destination / '.versions' / stamp / name / relative
            counts['versioned'] += 1
            if not dry_run:
                kept.parent.mkdir(parents=True, exist_ok=True)
                os.replace(target, kept)
        counts['copied'] += 1
        counts['bytes'] += os.path.getsize(path)
        if not dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name(target.name + '.partial')
            shutil.copy2(path, temporary)
            os.replace(temporary, target)


def main():
    config = json.loads(CONFIG.read_text(encoding='utf-8'))
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--only', nargs='+', choices=sorted(config['sources']), help='back up these sources only')
    parser.add_argument('--dry-run', action='store_true', help='count what would be copied, copy nothing')
    args = parser.parse_args()
    destination = Path(config['destination'])
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    receipt = {'schema': config['schema'], 'started': stamp, 'dry_run': args.dry_run, 'sources': {}}
    for name, location in config['sources'].items():
        if args.only and name not in args.only:
            continue
        source = Path(os.path.expanduser(location))
        source = source if source.is_absolute() else ROOT / source
        if not source.exists():
            print(f'{name}: {source} does not exist, skipped')
            receipt['sources'][name] = {'path': str(source), 'missing': True}
            continue
        counts = {'files': 0, 'copied': 0, 'versioned': 0, 'bytes': 0}
        back_up(name, source, destination, stamp, args.dry_run, counts)
        receipt['sources'][name] = {'path': str(source), **counts}
        print(f"{name}: {counts['files']} files, {counts['copied']} copied ({counts['bytes'] / 1e9:.2f} GB), "
              f"{counts['versioned']} earlier copies kept in .versions")
    receipt['finished'] = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    if not args.dry_run:
        (destination / 'receipts').mkdir(parents=True, exist_ok=True)
        (destination / 'receipts' / f'{stamp}.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
        print(f"receipt: {destination / 'receipts' / (stamp + '.json')}")


if __name__ == '__main__':
    main()
