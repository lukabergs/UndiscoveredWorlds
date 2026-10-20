"""Refuse files in places the repository layout does not name (configs/structure.json).

    python scripts/check_structure.py              # every tracked and untracked, not ignored file
    python scripts/check_structure.py --path P...  # judge paths as new files
    python scripts/check_structure.py --hook       # PreToolUse hook (Write, Edit, ...): JSON on stdin

A file is allowed when its folder is a named place and its kind (suffix) or name is one the place takes; per-run
HTML is allowed only where scripts/evaluation/ renders pages, scripts never under runs/, and runs/archive/ is
sealed. The hook judges only paths inside this repository, and only files that do not exist yet (except the
sealed folders): exit 0 allows, exit 2 blocks with the reason on stderr. Standard library only, so that the hook
runs with the plain interpreter. Beside check_layout.py (root items, includes, run directories) and
repository_scope.py (local-only data, blob size).
"""
import argparse
import fnmatch
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'configs/structure.json'
ADVICE = ('The named places are in configs/structure.json (AGENTS.md, "Where things go"): scratch goes to '
          'out/scratch/<topic>/ or the session scratchpad, pages are rendered by scripts/evaluation/. Use a named '
          'place; if a new place is needed, report the path to the orchestrator instead of working around this.')


def load(path=CONFIG):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def paths_of(entry):
    path = entry['path']
    return [path] if isinstance(path, str) else list(path)


def split(path):
    return [part for part in path.replace('\\', '/').split('/') if part not in ('', '.')]


def segment_matches(name, pattern, folder_name, free):
    if pattern == '*':
        return free or folder_name.match(name) is not None
    if pattern == '#':
        return name.isdigit()
    if any(char in pattern for char in '*?['):
        return fnmatch.fnmatchcase(name, pattern)
    return name == pattern


def folder_matches(folders, patterns, folder_name, free=True):
    """Whether the folder segments match the pattern segments ('**' is any depth)."""
    if not patterns:
        return not folders
    head, rest = patterns[0], patterns[1:]
    if head == '**':
        for depth in range(len(folders) + 1):
            if all(free or folder_name.match(name) for name in folders[:depth]) and \
                    folder_matches(folders[depth:], rest, folder_name, free):
                return True
        return False
    return bool(folders) and segment_matches(folders[0], head, folder_name, free) and \
        folder_matches(folders[1:], rest, folder_name, free)


def matches_any(folders, patterns, folder_name):
    return any(folder_matches(folders, split(pattern), folder_name) for pattern in patterns)


def judge(path, config, existing=False):
    """The reason a repository-relative file path is refused, or None when it is allowed.

    existing: the file is already there (the hook lets it be written again unless its folder is sealed).
    """
    folder_name = re.compile(config['folder_name'])
    parts = split(path)
    if not parts:
        return None
    for checkout in config.get('checkouts', ()):
        depth = len(split(checkout))
        if len(parts) > depth and folder_matches(parts[:depth], split(checkout), folder_name):
            parts = parts[depth:]
            break
    folders, name = parts[:-1], parts[-1]
    for sealed in config.get('sealed', ()):
        if folder_matches(folders, split(sealed['path']), folder_name):
            return sealed['reason']
    if existing:
        return None
    kind = os.path.splitext(name)[1].lower()
    relative = '/'.join(parts)
    for restriction in config.get('restrictions', ()):
        if kind not in restriction['kinds']:
            continue
        if 'only_in' in restriction and relative not in restriction.get('files', ()) and \
                not matches_any(folders, restriction['only_in'], folder_name):
            return f'{kind} files only in {", ".join(restriction["only_in"])}: {restriction["reason"]}'
        if 'not_in' in restriction and matches_any(folders, restriction['not_in'], folder_name):
            return f'no {kind} files in {", ".join(restriction["not_in"])}: {restriction["reason"]}'
    named = False
    for place in config['places']:
        free = place.get('free', False)
        if not any(folder_matches(folders, split(pattern), folder_name, free) for pattern in paths_of(place)):
            continue
        named = True
        kinds = place.get('kinds', ())
        if '*' in kinds or kind in kinds or name in place.get('names', ()):
            return None
    where = '/'.join(folders) + '/' if folders else 'the repository root'
    if named:
        return f'{name!r} is not a file kind or name that {where} takes'
    return f'the folder {where} is not a place the layout names'


def repository_path(path, root=ROOT, cwd=None):
    """The path relative to root with forward slashes, or None when it lies outside root (junctions not followed)."""
    absolute = os.path.abspath(os.path.join(cwd or os.getcwd(), path))
    base = os.path.abspath(str(root))
    if os.path.normcase(absolute) == os.path.normcase(base):
        return None
    prefix = os.path.normcase(base.rstrip('\\/')) + os.sep
    if not os.path.normcase(absolute).startswith(prefix):
        return None
    return absolute[len(prefix):].replace('\\', '/')


def hook_decision(event, config, root=ROOT):
    """The block message for a PreToolUse event, or None to allow it."""
    tool_input = event.get('tool_input') or {}
    path = tool_input.get('file_path') or tool_input.get('notebook_path')
    if not path:
        return None
    relative = repository_path(path, root, event.get('cwd'))
    if relative is None:
        return None
    existing = os.path.exists(os.path.join(str(root), relative))
    reason = judge(relative, config, existing)
    return f'BLOCKED: {relative}: {reason}. {ADVICE}' if reason else None


def repository_files(root=ROOT):
    listed = subprocess.check_output(['git', 'ls-files', '-co', '--exclude-standard', '-z'], cwd=root)
    return sorted({name for name in listed.decode('utf-8').split('\0') if name})


def full_issues(files, config):
    issues = []
    for name in files:
        reason = judge(name, config)
        if reason:
            issues.append(f'{name}: {reason}')
    return issues


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--hook', action='store_true', help='read a PreToolUse event on stdin; exit 2 blocks')
    group.add_argument('--path', nargs='+', help='judge these paths as new files')
    arguments = parser.parse_args(argv)
    config = load()
    if arguments.hook:
        message = hook_decision(json.loads(sys.stdin.read() or '{}'), config)
        if message:
            sys.stderr.write(message + '\n')
            return 2
        return 0
    if arguments.path:
        issues = []
        for path in arguments.path:
            relative = repository_path(path)
            reason = judge(relative, config) if relative is not None else None
            if reason:
                issues.append(f'{relative}: {reason}')
        print('\n'.join(issues) if issues else 'Allowed.')
        return 1 if issues else 0
    files = repository_files()
    issues = full_issues(files, config)
    if issues:
        print('Structure check failed:\n' + '\n'.join(issues))
        return 1
    print(f'Structure check passed: {len(files)} files in places configs/structure.json names.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
