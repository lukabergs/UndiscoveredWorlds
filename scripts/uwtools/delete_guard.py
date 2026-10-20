"""PreToolUse guard (Bash and PowerShell): refuse deleting anything git cannot give back.

The user's rule (2026-10-18): agents work without asking; a delete is allowed when git can restore it, so the
guard refuses, before the command runs:
- deleting an untracked or ignored file or folder (refs/, runs/, out/scratch, D: data, anything outside a
  repository), or a tracked file with uncommitted changes;
- git commands that discard untracked files or uncommitted work: `clean` (unless a dry run), `reset --hard`,
  `checkout --`/`restore` on a dirty tree, `stash drop`/`clear`, `worktree remove`, `branch -D`;
- the repository's own deleting commands only in their user-only forms: `storage.py finish --accept-newer`,
  `worktrees.py retire --allow-unlanded` (since 2026-10-20 the user lets agents run `worktrees.py retire`,
  `worktrees.py drop-branches` and `storage.py finish`: each verifies its own preconditions in code first);
- a deleter with no target on its command line (fed by a pipeline);
- deletes it cannot resolve to paths (variables, command substitution, xargs, find -delete, inline code with
  rmtree/remove/unlink, .NET Delete), and untracked or modified scripts whose new code deletes files.
Allowed without a check: the system temp directory (the session scratchpads), `out/build` trees and
`__pycache__` (all regenerable). Tracked, unmodified files and folders may be deleted: git restores them.

Untracked data an agent no longer needs goes to the trash instead (`scripts/uwtools/trash.py move PATH`: a move,
never a delete; only the user empties the trash).

Reads the hook JSON on stdin; exit 2 blocks the call and the message tells the agent what to do instead. Fails
closed: an error inside the guard, or a check running past DEADLINE, blocks (only exit 2 blocks in Claude Code,
so an uncaught exception or the hook's own timeout would let the command pass). guard_catastrophe.py runs before
it in every project and refuses drive-level and protected-tree deletes whatever git says.
The guard sees the command line only: a tracked, unmodified program that deletes internally passes (the
repository's are listed above).
"""
import glob
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time

DELETERS = {'rm', 'rmdir', 'unlink', 'shred', 'del', 'erase', 'rd', 'ri', 'remove-item', 'rmdir.exe', 'rm.exe'}
CODE_DELETE = re.compile(r'rmtree|os\.remove\b|os\.unlink|\.unlink\(|os\.rmdir|\.rmdir\(|removedirs|send2trash|'
                         r'\[(system\.)?io\.(file|directory)\]::delete|remove-item|\brm\s+-|\bdel\s+/|\brd\s+/s',
                         re.IGNORECASE)
SCRIPT_SUFFIXES = ('.py', '.ps1', '.sh', '.cmd', '.bat')
INTERPRETERS = {'python', 'python.exe', 'python3', 'py', 'uv', 'uv.exe', 'pwsh', 'powershell', 'bash', 'sh', 'source', '.'}
TEMP = os.path.normcase(os.path.realpath(tempfile.gettempdir()))
SEPARATORS = re.compile(r'&&|\|\||[;|\n]')
DEADLINE = time.monotonic() + 40  # the hook's timeout is 60 s; past it Claude Code would run the command


def block(reason):
    sys.stderr.write(
        'BLOCKED by scripts/uwtools/delete_guard.py: ' + reason + '\n'
        'Deleting untracked or ignored data, or uncommitted changes, is for the user. Do not retry or work around '
        'it: move untracked data you no longer need with scripts/uwtools/trash.py move PATH, or note the path and '
        'why in the master plan\'s "For the user" list, and continue without it.\n')
    sys.exit(2)


def git(top, *args):
    if time.monotonic() > DEADLINE:
        block('the check ran out of time (too many targets); name fewer paths per command.')
    try:
        result = subprocess.run(['git', '-C', top, *args], capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError) as error:
        block(f'git could not check the target ({error}).')
    return result.returncode, result.stdout


def to_native(token, cwd):
    token = token.strip('"\'')
    if token.startswith('~'):
        token = os.path.expanduser(token)
    match = re.match(r'^/([a-zA-Z])(/.*)?$', token)
    if match:
        token = match.group(1).upper() + ':' + (match.group(2) or '/')
    elif token == '/tmp' or token.startswith('/tmp/'):
        token = tempfile.gettempdir() + token[4:]
    return os.path.normpath(token if os.path.isabs(token) else os.path.join(cwd, token))


def repository(path):
    probe = path
    while not os.path.isdir(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            return None
        probe = parent
    code, out = git(probe, 'rev-parse', '--show-toplevel')
    return os.path.normpath(out.strip()) if code == 0 and out.strip() else None


def regenerable(path):
    norm = os.path.normcase(os.path.realpath(path))
    parts = norm.replace('\\', '/').split('/')
    if norm == TEMP or norm.startswith(TEMP + os.sep):
        return True
    if '__pycache__' in parts or '.pytest_cache' in parts:
        return True
    return any(parts[i] == 'out' and parts[i + 1] == 'build' for i in range(len(parts) - 1))


def check_target(path):
    if not os.path.lexists(path) or regenerable(path):
        return
    top = repository(path)
    if top is None:
        block(f'{path} is outside any git repository (untracked by definition).')
    relative = os.path.relpath(path, top)
    if relative == '.' or relative.split(os.sep)[0] == '.git':
        block(f'{path} is a repository root or its .git.')
    code, others = git(top, 'ls-files', '--others', '--directory', '-z', '--', relative)
    if code != 0 or others.strip('\0'):
        block(f'{path} is or contains untracked or ignored files ({others.split(chr(0))[0] or "git failed"}).')
    code, tracked = git(top, 'ls-files', '-z', '--', relative)
    if code != 0 or not tracked.strip('\0'):
        block(f'{path} is not tracked by git.')
    code, dirty = git(top, 'status', '--porcelain', '-uno', '--', relative)
    if code != 0 or dirty.strip():
        block(f'{path} has uncommitted changes ({dirty.strip().splitlines()[0] if dirty.strip() else "git failed"}).')


def expand(token, cwd):
    if re.search(r'[$`%{}]', token) or token.startswith('@('):
        block(f'cannot resolve the delete target {token!r} to a path (variables or expansion): name paths literally.')
    native = to_native(token, cwd)
    if any(ch in token for ch in '*?['):
        return glob.glob(native, recursive=True)
    return [native]


def dirty_tree(cwd):
    top = repository(cwd)
    if top is None:
        return False
    code, out = git(top, 'status', '--porcelain', '-uno')
    return code != 0 or bool(out.strip())


def check_git(words, cwd):
    args = words[1:]
    while args and args[0].startswith('-'):
        if args[0] == '-C' and len(args) > 1:
            cwd = to_native(args[1], cwd)
            args = args[2:]
        else:
            args = args[1:]
    if not args:
        return
    sub, rest = args[0].lower(), [a.lower() for a in args[1:]]
    if sub == 'clean' and not ({'-n', '--dry-run'} & set(rest)):
        block('git clean deletes untracked files.')
    if sub == 'reset' and {'--hard', '--merge', '--keep'} & set(rest) and dirty_tree(cwd):
        block('git reset --hard on a tree with uncommitted changes.')
    if sub == 'checkout' and ('--' in rest or '.' in rest or '-f' in rest or '--force' in rest) and dirty_tree(cwd):
        block('git checkout over uncommitted changes.')
    staged_only = ('--staged' in rest or '-s' in rest) and not ({'--worktree', '-w'} & set(rest))
    if sub == 'restore' and not staged_only and dirty_tree(cwd):
        block('git restore over uncommitted changes.')
    if sub == 'stash' and rest and rest[0] in ('drop', 'clear'):
        block('git stash drop/clear discards stashed work.')
    if sub == 'worktree' and rest and rest[0] == 'remove':
        block('git worktree remove deletes the worktree\'s untracked files (scripts/uwtools/worktrees.py preserve, '
              'then the user retires it).')
    if sub == 'branch' and ('-D' in args[1:] or ('--delete' in rest and ('--force' in rest or '-f' in rest))):
        block('git branch -D can drop unmerged commits: use scripts/uwtools/worktrees.py drop-branches, which '
              'bundles and verifies them first.')


def check_script(path):
    if not os.path.isfile(path):
        return
    top = repository(path)
    if top is None:
        text = open(path, encoding='utf-8', errors='replace').read()
    else:
        relative = os.path.relpath(path, top)
        code, tracked = git(top, 'ls-files', '--', relative)
        if code == 0 and tracked.strip():
            code, diff = git(top, 'diff', 'HEAD', '--unified=0', '--', relative)
            text = '\n'.join(line for line in diff.splitlines() if line.startswith('+') and not line.startswith('+++'))
        else:
            text = open(path, encoding='utf-8', errors='replace').read()
    if CODE_DELETE.search(text):
        block(f'{path} is untracked or modified and its new code deletes files; deletes in scripts cannot be checked.')


def check_segment(segment, cwd):
    try:
        words = shlex.split(segment, posix=False)
    except ValueError:
        words = segment.split()
    words = [w.strip('"\'') for w in words if w.strip('"\'')]
    while words and words[0].lower() in ('sudo', 'command', 'cmd', 'cmd.exe', '/c', '&', 'call', 'time', 'nohup'):
        words = words[1:]
    if not words:
        return cwd
    head = os.path.basename(words[0]).lower()
    if head in ('cd', 'set-location', 'sl', 'pushd', 'chdir') and len(words) > 1:
        target = words[1] if not words[1].startswith('-') else (words[2] if len(words) > 2 else '')
        return to_native(target, cwd) if target and not re.search(r'[$`%]', target) else cwd
    lowered = segment.lower()
    if 'storage.py' in lowered and re.search(r'\bfinish\b', lowered) and '--accept-newer' in lowered:
        block('storage.py finish --accept-newer deletes old versions the storage copy no longer has; the user runs it.')
    if 'worktrees.py' in lowered and re.search(r'\bretire\b', lowered) and '--allow-unlanded' in lowered:
        block('worktrees.py retire --allow-unlanded deletes a worktree whose work is not in main; the user runs it.')
    if head in ('git', 'git.exe'):
        check_git(words, cwd)
        return cwd
    if head == 'xargs' and any(os.path.basename(w).lower() in DELETERS for w in words[1:]):
        block('xargs feeds delete targets that cannot be checked.')
    if head == 'find' and ('-delete' in words or any(os.path.basename(w).lower() in DELETERS for w in words)):
        block('find -delete / -exec rm targets cannot be checked.')
    if head in ('mv', 'move-item', 'mi', 'move') and words[-1].lower() in ('/dev/null', '$null', 'nul'):
        block('moving to the null device deletes.')
    if head in DELETERS:
        if not [w for w in words[1:] if not w.startswith('-') and not (w.startswith('/') and len(w) <= 3)]:
            block(f'{head} has no target on its command line (fed by a pipeline): name paths literally.')
        for token in words[1:]:
            if token.startswith('-') or (token.startswith('/') and len(token) <= 3 and head in ('del', 'erase', 'rd', 'rmdir')):
                continue
            for path in expand(token, cwd):
                check_target(path)
    if head in INTERPRETERS or head.endswith(SCRIPT_SUFFIXES):
        for token in words:
            if token.lower().endswith(SCRIPT_SUFFIXES) and not re.search(r'[$`%]', token):
                check_script(to_native(token, cwd))
    return cwd


def main():
    try:
        check(json.load(sys.stdin))
    except SystemExit:
        raise
    except BaseException as error:  # fail closed: a guard that cannot decide blocks
        block(f'the guard could not check this command ({type(error).__name__}: {error}).')


def check(payload):
    command = (payload.get('tool_input') or {}).get('command') or ''
    cwd = payload.get('cwd') or os.getcwd()
    lowered = ' ' + command.lower() + ' '
    if re.search(r'\s-c\s|\s-command\s|\s-e\s|<<', lowered) and CODE_DELETE.search(command):
        block('inline code or a heredoc deletes files; its targets cannot be checked.')
    if re.search(r'\[(system\.)?io\.(file|directory)\]::delete|\.delete\(\)', lowered):
        block('.NET Delete targets cannot be checked.')
    for segment in SEPARATORS.split(command):
        if segment.strip():
            cwd = check_segment(segment.strip(), cwd)


if __name__ == '__main__':
    main()
