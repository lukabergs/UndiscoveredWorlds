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
                         r'\[(system\.)?io\.(file|directory)\]::delete|remove-item|\brm\s+-|\bdel\s+/|\brd\s+/s|'
                         r'remove_tree|accept_newer\s*=\s*true|allow_unlanded\s*=\s*true|[\'"]rm[\'"]|deletefolder|'
                         r'\brm\s*\(|rmsync|unlinksync|icacls|cacls|takeown|set-acl|setaccesscontrol|'
                         r'setsecuritydescriptor', re.IGNORECASE)
INLINE = re.compile(r'(?<![\w.-])(python3?|py|node|perl|ruby|bash|sh|zsh|pwsh|powershell)(\.exe)?\b[^|;&\n]*?'
                    r'\s-(c|e|command)\b|<<', re.IGNORECASE)
SCRIPT_SUFFIXES = ('.py', '.ps1', '.sh', '.cmd', '.bat', '.js', '.mjs', '.cjs')
INTERPRETERS = {'python', 'python.exe', 'python3', 'py', 'uv', 'uv.exe', 'pwsh', 'powershell', 'bash', 'sh', 'source', '.',
                'node', 'node.exe'}
READ_ONLY = {'rg', 'grep', 'egrep', 'fgrep', 'findstr', 'select-string', 'sls', 'cat', 'type', 'get-content', 'gc',
             'less', 'more', 'head', 'tail', 'ls', 'dir', 'get-childitem', 'gci', 'wc', 'echo', 'write-output',
             'printf', 'test-path', 'sed', 'awk'}  # not git: git rm -f and git -c aliases delete
CREATORS = re.compile(r'(?<![\w.-])(mv|move|move-item|mi|ren|rename|rename-item|rni|cp|copy|copy-item|cpi|xcopy|'
                      r'robocopy|mkdir|md|new-item|ni|ln|mklink|tar|unzip|expand-archive)(?![\w-])', re.IGNORECASE)
CREATES = False  # set per command: a missing delete target may be one the same command creates
TEMP = os.path.normcase(os.path.realpath(tempfile.gettempdir()))
# Captured: a deleter after a single | is fed by a pipeline. A lone & runs the next command too (cmd, bash); 2>&1,
# &> and >& are redirections.
SEPARATORS = re.compile(r'(&&|\|\||[;\n]|\||(?<![<>])&(?!>))')
DELETE_TEXT = re.compile(
    r'(?<![\w.-])(rm|rmdir|rd|del|erase|unlink|shred|rimraf|remove-item|ri|clear-recyclebin|forfiles|clear-content|'
    r'clc)(?![\w-])|rmtree|rmsync|rmdirsync|unlinksync|remove_directory|removedirs|os\.remove|\.unlink\(|\.rmdir\(|'
    r'deletedirectory|deletefile|deletefolder|::delete|\.delete\s*\(|\brm\s*\(|send2trash|--remove-files|--delete|'
    r'(?<![\w-])-delete\b|-e\s+remove\b|'
    r'(?<![\w-])[/-](mir|purge|mov|move)(?![\w-])', re.IGNORECASE)
TOKEN_SPLIT = re.compile(r'''[\s;|&<>(){}\[\],="'`]+''')
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
        result = subprocess.run(['git', '-C', top, *args], capture_output=True, text=True, encoding='utf-8',
                                errors='replace', timeout=max(0.5, min(20.0, DEADLINE + 10 - time.monotonic())))
    except (OSError, subprocess.SubprocessError) as error:
        block(f'git could not check the target ({error}).')
    return result.returncode, result.stdout


def to_native(token, cwd):
    token = re.sub(r'["\']', '', token)  # quotes anywhere: "/c/dev/meta"/ref is /c/dev/meta/ref
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


def check_target(path, target=False):
    """target: the path is a deleter's own argument (not a word the backstop scanned)."""
    if not os.path.lexists(path):
        if target and CREATES and not regenerable(path):
            block(f'{path} does not exist yet but this command also moves, copies or creates something: '
                  'delete in a separate command.')
        return
    if regenerable(path):
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
    if cwd is None and not os.path.isabs(to_native(token, '.')) and not re.match(r'^/[a-zA-Z](/|$)', token):
        block(f'cannot resolve the delete target {token!r}: an earlier cd could not be followed; use full paths.')
    native = to_native(token, cwd or '.')
    if '**' in token:
        block(f'{token!r}: recursive wildcards in a delete cannot be checked in time; name folders literally.')
    if any(ch in token for ch in '*?['):
        found = []
        for path in glob.iglob(native, recursive=True, include_hidden=True):
            if time.monotonic() > DEADLINE:
                block('the check ran out of time (a wildcard over too many files); name fewer paths per command.')
            found.append(path)
        return found
    return [native]


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
              'then the user retires it).')
    if sub == 'branch' and (any(re.fullmatch(r'-[a-zA-Z]*D[a-zA-Z]*', a) for a in args[1:])
                            or any(re.fullmatch(r'-[a-z]*d[a-z]*f[a-z]*|-[a-z]*f[a-z]*d[a-z]*', a) for a in rest)
                            or (({'--delete', '-d'} & set(rest)) and ({'--force', '-f'} & set(rest)))):
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


WRAPPERS = {'sudo', 'command', 'builtin', 'exec', 'cmd', 'cmd.exe', '/c', '//c', '/k', '//k', '&', 'call', 'time',
            'nohup', 'env', 'do', 'then', 'else', 'elif', '!', '(', '{', 'npx', 'nice', 'timeout'}


def check_segment(segment, cwd, piped=False):
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
    head = os.path.basename(words[0]).lower()
    if head in ('cd', 'chdir') and len(words) == 1:
        return os.path.expanduser('~')
    if head in ('cd', 'set-location', 'sl', 'pushd', 'chdir', 'push-location') and len(words) > 1:
        rest = [w for w in words[1:] if not w.startswith('-') and w.lower() != '/d']
        target = rest[0] if rest else ''
        # A folder the guard cannot follow makes every later relative target unknown (None): deleters block.
        return to_native(target, cwd) if target and cwd and not re.search(r'[$`%(]', target) else None
    if head in ('popd', 'pop-location'):
        return None
    lowered = segment.lower()
    # Both parsers refuse abbreviated options (allow_abbrev=False); the prefixes catch them anyway.
    if 'storage' in lowered and re.search(r'\bfinish\b', lowered) and '--acc' in lowered:
        block('storage.py finish --accept-newer deletes old versions the storage copy no longer has; the user runs it.')
    if 'worktrees' in lowered and re.search(r'\bretire\b', lowered) and '--allow' in lowered:
        block('worktrees.py retire --allow-unlanded deletes a worktree whose work is not in main; the user runs it.')
    if os.path.basename(words[0]).lower().removesuffix('.exe') == 'robocopy':
        switches = {w.lower().lstrip('/-') for w in words[1:] if w.startswith(('/', '-'))}
        places = [w for w in words[1:] if not w.startswith(('/', '-'))][:2]
        if switches & {'mir', 'purge', 'move', 'mov'}:
            if len(places) < 2:
                block('robocopy /MIR, /PURGE or /MOVE without a literal source and destination.')
            for path in expand(places[1], cwd):
                check_target(path)  # the mirror's destination (it may not exist yet: robocopy creates it)
            for path in expand(places[0], cwd) if switches & {'move', 'mov'} else []:
                check_target(path, target=True)
    if head in ('git', 'git.exe'):
        check_git(words, cwd)
        return cwd
    if head == 'xargs' and any(os.path.basename(w).lower() in DELETERS for w in words[1:]):
        block('xargs feeds delete targets that cannot be checked.')
    if head == 'find' and ('-delete' in words or any(os.path.basename(w).lower() in DELETERS for w in words)):
        block('find -delete / -exec rm targets cannot be checked.')
    if head in ('mv', 'move-item', 'mi', 'move') and words[-1].lower() in ('/dev/null', '$null', 'nul'):
        block('moving to the null device deletes.')
    if head in ('mv', 'move-item', 'mi', 'move'):
        # A move into the temp folder or out/build is a delete in two steps (both are exempt): check the sources.
        places = [w for w in words[1:] if not w.startswith('-') and not re.fullmatch(r'//?[a-zA-Z]{1,2}', w)]
        if len(places) >= 2 and not re.search(r'[$`%]', places[-1]) and regenerable(to_native(places[-1], cwd or '.')):
            for source in places[:-1]:
                for path in expand(source, cwd):
                    if not regenerable(path):
                        check_target(path)
    if head in DELETERS:
        if piped:
            block(f'{head} fed by a pipeline: its targets cannot be checked; name paths literally.')
        if head in ('del', 'erase') and {'/s', '//s'} & {w.lower() for w in words} and any(
                ch in w for w in words[1:] for ch in '*?'):
            block('del /s with a wildcard deletes in every subfolder; the guard checks only the first.')
        if not [w for w in words[1:] if not w.startswith('-') and not (w.startswith('/') and len(w) <= 3)]:
            block(f'{head} has no target on its command line (fed by a pipeline): name paths literally.')
        for token in words[1:]:
            if token.startswith('-') or (token.startswith('/') and len(token) <= 3 and head in ('del', 'erase', 'rd', 'rmdir')):
                continue
            for path in expand(token, cwd):
                check_target(path, target=True)
    if head in INTERPRETERS or head.endswith(SCRIPT_SUFFIXES) or re.search(r'\s', words[0]):
        for token in (part for word in words for part in word.split()):  # also inside bash -c "python x.py"
            if token.lower().endswith(SCRIPT_SUFFIXES) and not re.search(r'[$`%]', token):
                if cwd is None and not os.path.isabs(to_native(token, '.')):
                    block(f'cannot find the script {token!r}: an earlier cd could not be followed; use a full path.')
                check_script(to_native(token, cwd or '.'))
    return cwd


def main():
    try:
        check(json.loads(sys.stdin.buffer.read().decode('utf-8', errors='replace')))
    except SystemExit:
        raise
    except BaseException as error:  # fail closed: a guard that cannot decide blocks
        block(f'the guard could not check this command ({type(error).__name__}: {error}).')
    sys.exit(0)  # hook_runner.py treats a run without a decision as a block


def backstop(command, cwd):
    """Independent of the parser: a command holding a delete word anywhere may not name, anywhere, an existing path
    git cannot restore (quoted inner commands, comma lists, cd and Push-Location targets, cmake -E, rsync, ...)."""
    if re.search(r'(^|\s)//?s(\s|$)', command, re.IGNORECASE) and re.search(r'[*?]', command):
        block('a recursive delete (/s) with a wildcard: it acts in every subfolder; name folders literally.')
    if re.search(r'\\ (?=[^\s-])', command):
        block('a delete with an escaped space (run\\ 7): quote the path instead.')
    if re.search(r'[\\`]\r?\n', command):
        block('a delete in a command continued over lines (a trailing \\ or backtick joins words); one line each.')
    previous = ''
    for raw in TOKEN_SPLIT.split(command):
        after_cd, previous = previous in ('cd', 'set-location', 'sl', 'pushd', 'push-location', 'chdir', '-c', '/d'), \
            raw.lower()
        if not raw or raw.startswith('-') or re.fullmatch(r'//?[a-zA-Z]{1,3}', raw) or re.search(r'[$%]', raw):
            continue  # options, cmd switches and variables (the parser refuses deleters with variable targets)
        if raw in ('/', '\\', '.') or raw.lower() in ('nul', '/dev/null') or raw.startswith('/dev/'):
            continue
        for path in expand(raw, cwd):
            top = repository(path) if after_cd and os.path.isdir(path) else None
            if top is not None and os.path.normcase(os.path.normpath(path)) == os.path.normcase(top) \
                    and os.path.isdir(os.path.join(top, '.git')):
                continue  # a main checkout's root right after cd or git -C; the parser checks the targets
            check_target(path)


def heads(command):
    """The command word of each segment, after the common wrappers."""
    found = []
    for part in SEPARATORS.split(command)[::2]:
        words = [w.strip('"\'') for w in part.split() if w.strip('"\'')]
        while words and (words[0].lower() in WRAPPERS or re.fullmatch(r'\w+=\S*', words[0])):
            words = words[1:]
        if words:
            found.append(os.path.basename(words[0]).lower())
    return found


def check(payload):
    global CREATES
    command = (payload.get('tool_input') or {}).get('command') or ''
    cwd = to_native(payload.get('cwd') or os.getcwd(), os.getcwd())  # also a Git Bash /c/... form
    CREATES = bool(CREATORS.search(command))
    # The backstop runs unless every command is a read-only tool with no unquoted ( (PowerShell runs a grouped
    # command inside echo (...)).
    unquoted = re.sub(r'"[^"]*"|\'[^\']*\'', '', command)
    if DELETE_TEXT.search(command) and not (all(head in READ_ONLY for head in heads(command)) and '(' not in unquoted):
        backstop(command, cwd)
    lowered = ' ' + command.lower() + ' '
    searching = all(head in READ_ONLY for head in heads(command))
    if INLINE.search(command) and CODE_DELETE.search(command) and not searching:
        block('inline code or a heredoc deletes files; its targets cannot be checked.')
    if re.search(r'\[(system\.)?io\.(file|directory)\]::delete|\.delete\s*\(|deletefolder', lowered) and not searching:
        block('.NET or COM Delete targets cannot be checked.')
    if re.search(r'(?<![\w-])git(\.exe)?\s.*\s-c\s', command) and DELETE_TEXT.search(command):  # -c, not -C
        block('git -c with a delete (an alias or a hook command can run it).')
    parts = SEPARATORS.split(command)
    for position in range(0, len(parts), 2):
        if parts[position].strip():
            cwd = check_segment(parts[position].strip(), cwd, position > 0 and parts[position - 1] == '|')


if __name__ == '__main__':
    main()
