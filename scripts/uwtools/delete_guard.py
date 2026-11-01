"""PreToolUse guard (Bash and PowerShell): refuse git commands that discard work git cannot give back.

Since 2026-11-01 (the user: too many false positives) this guard no longer judges deleted paths: untracked data in
runs/, out/scratch and the like may be cleaned up, and guard_catastrophe.py (global) still refuses drive-level and
protected-tree deletes. It keeps one job, because several sessions share this working tree and one agent's git
command can wipe another's uncommitted work. Refused, before the command runs:
- `clean` (unless a dry run), `reset --hard/--merge/--keep` on a dirty tree, `checkout`/`restore`/`switch --force`,
  `rm -f`, `checkout-index -f`, `read-tree -u` and `merge`/`rebase`/`am`/`cherry-pick`/`revert --abort` over
  uncommitted changes, `stash drop`/`clear`, `worktree remove`, `branch -D` (use
  `scripts/uwtools/worktrees.py drop-branches`, which bundles first) and `branch -M`/`-C` (they overwrite a branch);
- `worktrees.py retire --allow-unlanded` (the user's only).

Reads the hook JSON on stdin; exit 2 blocks the call. Fails closed: an error inside the guard, or a git check past
DEADLINE, blocks (hook_runner.py turns any other failure into a block too).
"""
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guard_text  # noqa: E402  (beside this file; hook_runner.py turns an import failure into a block)

DEADLINE = time.monotonic() + 40  # the hook's timeout is 60 s; past it Claude Code would run the command
WRAPPERS = {'sudo', 'command', 'builtin', 'exec', 'cmd', 'cmd.exe', '/c', '//c', '/k', '//k', '&', 'call', 'time',
            'nohup', 'env', 'do', 'then', 'else', 'elif', '!', '(', '{', 'npx', 'nice', 'timeout'}
SHELLS = {'bash', 'sh', 'zsh', 'cmd', 'pwsh', 'powershell', 'iex', 'invoke-expression', 'start-process'}


def block(reason):
    sys.stderr.write(
        'BLOCKED by scripts/uwtools/delete_guard.py: ' + reason + '\n'
        'Other sessions\' uncommitted work shares this tree. Do not retry or work around it: commit or leave the '
        'changes, or note what you wanted and why in the master plan\'s "For the user" list, and continue.\n')
    sys.exit(2)


def git(top, *args):
    if time.monotonic() > DEADLINE:
        block('the check ran out of time.')
    try:
        result = subprocess.run(['git', '-C', top, *args], capture_output=True, text=True, encoding='utf-8',
                                errors='replace', timeout=max(0.5, min(20.0, DEADLINE + 10 - time.monotonic())))
    except (OSError, subprocess.SubprocessError) as error:
        block(f'git could not check the tree ({error}).')
    return result.returncode, result.stdout


def to_native(token, cwd):
    token = re.sub(r'["\']', '', token)
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


def dirty_tree(cwd):
    top = repository(cwd)
    if top is None:
        return False
    code, out = git(top, 'status', '--porcelain', '-uno')
    return code != 0 or bool(out.strip())


def dirty_paths(cwd, paths):
    """True when any of the named paths (relative to cwd) has uncommitted changes."""
    top = repository(cwd)
    if top is None:
        return False
    named = [p for p in paths if not p.startswith('-') and os.path.lexists(to_native(p, cwd))]
    if not named:
        return False
    code, out = git(cwd, 'status', '--porcelain', '-uno', '--', *named)
    return code != 0 or bool(out.strip())


def check_git(words, cwd):
    args = words[1:]
    while args and args[0].startswith('-'):  # global options; -C and -c take a value
        if args[0] == '-C' and len(args) > 1:
            cwd = to_native(args[1], cwd)
            args = args[2:]
        elif args[0] in ('-c', '--git-dir', '--work-tree', '--namespace') and len(args) > 1:
            args = args[2:]
        else:
            args = args[1:]
    if not args:
        return
    sub, rest = args[0].lower(), [a.lower() for a in args[1:]]
    flags = {a for a in rest if a.startswith('-')}
    if sub == 'clean' and not ({'-n', '--dry-run'} & set(rest)):
        block('git clean deletes untracked files.')
    if sub == 'rm' and any(re.fullmatch(r'-[a-z]*f[a-z]*|--force', a) for a in rest) and '--cached' not in rest \
            and dirty_paths(cwd, [a for a in args[1:] if not a.startswith('-')]):
        block('git rm -f on a file with uncommitted changes.')
    if sub in ('merge', 'rebase', 'am', 'cherry-pick', 'revert') and '--abort' in rest and dirty_tree(cwd):
        block(f'git {sub} --abort on a tree with uncommitted changes (it resets them).')
    if sub == 'checkout-index' and ({'-f', '--force'} & set(rest)) and dirty_tree(cwd):
        block('git checkout-index -f over uncommitted changes.')
    if sub == 'read-tree' and '-u' in rest and dirty_tree(cwd):
        block('git read-tree -u over uncommitted changes.')
    if sub == 'reset' and {'--hard', '--merge', '--keep'} & set(rest) and dirty_tree(cwd):
        block('git reset --hard on a tree with uncommitted changes.')
    if sub == 'checkout' and (('--' in rest or '.' in rest or '-f' in rest or '--force' in rest) and dirty_tree(cwd)
                              or dirty_paths(cwd, args[1:])):
        block('git checkout over uncommitted changes (a path argument discards its edits).')
    if sub == 'switch' and ({'-f', '--force', '--discard-changes'} & flags) and dirty_tree(cwd):
        block('git switch --force / --discard-changes over uncommitted changes.')
    staged_only = ('--staged' in args or '-S' in args) and not ({'--worktree', '-W'} & set(args))  # -s is --source
    if sub == 'restore' and not staged_only:
        specs = [a for a in args[1:] if not a.startswith('-')]
        unusual = any(not os.path.lexists(to_native(a, cwd)) or a.startswith(':') or re.search(r'[*?\[]', a)
                      for a in specs)  # pathspecs, globs, sources: fall back to the whole tree
        if (unusual and dirty_tree(cwd)) or dirty_paths(cwd, args[1:]) or '.' in rest:
            block('git restore over uncommitted changes.')
    if sub == 'stash' and rest and rest[0] in ('drop', 'clear'):
        block('git stash drop/clear discards stashed work.')
    if sub == 'worktree' and rest and rest[0] == 'remove':
        block('git worktree remove deletes the worktree\'s untracked files (scripts/uwtools/worktrees.py preserve, '
              'then retire).')
    if sub == 'branch' and any(re.fullmatch(r'-[a-zA-Z]*[MC][a-zA-Z]*', a) for a in args[1:]):
        block('git branch -M / -C overwrites an existing branch, whose commits may be unmerged.')
    if sub == 'branch' and (any(re.fullmatch(r'-[a-zA-Z]*D[a-zA-Z]*', a) for a in args[1:])
                            or any(re.fullmatch(r'-[a-z]*d[a-z]*f[a-z]*|-[a-z]*f[a-z]*d[a-z]*', a) for a in rest)
                            or (({'--delete', '-d'} & set(rest)) and ({'--force', '-f'} & set(rest)))):
        block('git branch -D can drop unmerged commits: use scripts/uwtools/worktrees.py drop-branches, which '
              'bundles and verifies them first.')


def check_segment(segment, cwd, level=0):
    """Check one command (and the quoted commands a shell runs); return the working folder after it."""
    try:
        words = shlex.split(segment, posix=False)
    except ValueError:
        words = segment.split()
    words = [w.strip('"\'') for w in words if w.strip('"\'')]
    while words:
        first = words[0].lower()
        if first in WRAPPERS or re.fullmatch(r'\w+=\S*', words[0]) or re.fullmatch(r'-?\d+[smh]?|-n', first):
            words = words[1:]
        elif first[:1] in ('(', '{'):
            words[0] = words[0][1:]
        else:
            break
    if not words:
        return cwd
    head = os.path.basename(words[0]).lower().removesuffix('.exe')
    if head in ('cd', 'set-location', 'sl', 'pushd', 'chdir', 'push-location'):
        rest = [w for w in words[1:] if not w.startswith('-') and w.lower() != '/d']
        return to_native(rest[0], cwd) if rest and not re.search(r'[$`%(]', rest[0]) else cwd
    lowered = segment.lower()
    # The parser refuses abbreviated options (allow_abbrev=False); the prefix catches them anyway.
    if 'worktrees' in lowered and re.search(r'\bretire\b', lowered) and '--allow' in lowered:
        block('worktrees.py retire --allow-unlanded deletes a worktree whose work is not in main; the user runs it.')
    if head == 'git':
        check_git(words, cwd)
    if head in SHELLS or re.search(r'\s', words[0]):
        if level >= 3:
            block('commands nested more than three deep.')
        for word in (w for w in words if re.search(r'\s', w)):  # bash -c "...", cmd /c "...": their commands
            inner = cwd
            for part, _, _ in guard_text.split_command(word):
                inner = check_segment(part, inner, level + 1)
    return cwd


def main():
    try:
        payload = json.loads(sys.stdin.buffer.read().decode('utf-8', errors='replace'))
        command = (payload.get('tool_input') or {}).get('command') or ''
        cwd = to_native(payload.get('cwd') or os.getcwd(), os.getcwd())  # also a Git Bash /c/... form
        if re.search(r'\bgit\b|worktrees', command, re.IGNORECASE):
            # Text (a heredoc or here-string commit message) is not commands (guard_text.readable).
            command, _ = guard_text.readable(command, payload.get('tool_name') == 'PowerShell')
            for segment, _, _ in guard_text.split_command(command):
                cwd = check_segment(segment, cwd)
    except SystemExit:
        raise
    except BaseException as error:  # fail closed: a guard that cannot decide blocks
        block(f'the guard could not check this command ({type(error).__name__}: {error}).')
    sys.exit(0)  # hook_runner.py treats a run without a decision as a block


if __name__ == '__main__':
    main()
