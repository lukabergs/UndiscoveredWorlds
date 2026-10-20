"""PreToolUse guard (Bash and PowerShell, every project): refuse commands that could wipe a drive or a protected tree.

Added 2026-10-20 at the user's request ("people have complained about their AI agents trying to delete whole
drives"). Registered in ~/.claude/settings.json through hook_runner.py, which turns any failure (a syntax error,
an exception, a run without a decision) into a block: Claude Code blocks only on exit code 2. It needs no git and
is independent of delete_guard.py (which decides what git can restore). The user asked that it never be tested
by attempting deletes: change it by reading, and have a reviewer read the change.

The rule is that a delete must be simple enough to check:
1. Disk, partition, volume, shadow-copy and registry-delete tools are refused (diskpart, Clear-Disk, format X:, ...).
2. A command that contains a delete word (rm, rd, del, erase, Remove-Item, rmtree, rmSync, unlink, ...) must be a
   recognised deleter (rm, rmdir, rd, del, erase, unlink, shred, rimraf, Remove-Item, cmake -E rm, find -delete,
   robocopy /MIR, rsync --delete, git clean) whose targets are literal paths on its own command line: not fed by a
   pipeline, no variables, braces, command substitution or line continuation, wildcards only in the last path
   component, and no relative path after a `cd` the guard could not follow. Deletes inside inline code (python -c,
   node -e, heredoc lines), loops and script blocks are refused: write a literal command or a script file.
3. Every delete target, move source and mirror place is checked against the protected places: a drive root or any
   path within two levels of one, the home folder, its direct children and its ancestors, Windows, Program Files,
   ProgramData, any .git, PROTECTED_TREES (whole subtrees) and ROOT_ONLY folders (the folder itself; its children
   may go). Literal destinations of copies, moves, writers and redirections may not lie in DESTINATION_TREES.
Commands are split on ; && || | outside quotes and on every newline (inside quotes too: more segments is the safe
side). Read-only commands (rg, grep, cat, ls, echo, git log, ...) are not checked for delete words unless they
hold command substitution or an unbalanced quote; a git commit message given with -m "..." is text.
"""
import json
import os
import re
import shlex
import sys
import tempfile

HOME = os.path.normcase(os.path.normpath(os.path.expanduser('~')))
# Whole subtrees no agent command may delete, move or overwrite (the user's data, backups and the guards).
PROTECTED_TREES = [
    r'C:\dev\undiscovered-worlds\refs',
    r'C:\dev\undiscovered-worlds\.git',
    r'C:\dev\meta',
    r'D:\dev\bak',
    r'D:\dev\undiscovered-worlds\worktree-archive',
    r'D:\dev\undiscovered-worlds\receipts',
    r'D:\dev\undiscovered-worlds\trash',
    os.path.join(HOME, '.claude', 'hooks'),
    os.path.join(HOME, '.claude', 'settings.json'),
    os.path.join(HOME, '.codex'),
]
# Folders that may not be deleted or moved whole (nor any folder above them); their children may be.
ROOT_ONLY = [
    r'C:\dev\undiscovered-worlds\out', r'C:\dev\undiscovered-worlds\runs', r'C:\dev\undiscovered-worlds\scripts',
    r'C:\dev\undiscovered-worlds\src', r'D:\dev\undiscovered-worlds\out', r'D:\dev\undiscovered-worlds\out\scratch',
    r'D:\dev\undiscovered-worlds\out\evaluation', r'D:\dev\undiscovered-worlds\runs',
    os.path.join(HOME, '.claude', 'projects'),
]
SYSTEM_TREES = [r'C:\Windows', r'C:\Program Files', r'C:\Program Files (x86)', r'C:\ProgramData', r'C:\Users\Public',
                r'C:\Users\Default']
# Destinations no copy, move, writer or redirection may write into (refs/ is left out: data scripts write there).
DESTINATION_TREES = SYSTEM_TREES + PROTECTED_TREES[1:]

DISK_TOOLS = re.compile(
    r'(?<![\w-])(diskpart|clear-disk|format-volume|remove-partition|initialize-disk|clear-volume|remove-physicaldisk|'
    r'reset-physicaldisk|set-partition|resize-partition|mkfs(\.\w+)?|wipefs|sfdisk|cfdisk|sgdisk|sdelete(64)?|bcdedit|'
    r'manage-bde|remove-itemproperty)(?![\w-])|\b(vssadmin|wbadmin)\b.*\bdelete\b|\bwmic\b.*\bshadowcopy\b.*\bdelete\b|'
    r'\breg(\.exe)?\s+delete\b', re.IGNORECASE)
DELETE_WORDS = re.compile(
    r'(?<![\w.-])(rm|rmdir|rd|del|erase|unlink|shred|rimraf|remove-item|ri|clear-recyclebin|forfiles)(?![\w-])|'
    r'rmtree|rmsync|rmdirsync|unlinksync|remove_directory|removedirs|os\.remove|\.unlink\(|\.rmdir\(|'
    r'deletedirectory|deletefile|::delete|\.delete\(\)|send2trash', re.IGNORECASE)
SUBSTITUTION = re.compile(r'\$\(|`|<\(|>\(')
DELETERS = {'rm', 'rmdir', 'rd', 'del', 'erase', 'unlink', 'shred', 'rimraf', 'remove-item', 'ri'}
CMD_STYLE = {'rd', 'del', 'erase', 'rmdir'}  # their /s /q /f are switches; rm's /c is a path
MOVERS = {'mv', 'move', 'move-item', 'mi', 'ren', 'rename', 'rename-item', 'rni'}
COPIERS = {'cp', 'copy', 'copy-item', 'cpi', 'xcopy'}
WRITERS = {'set-content', 'out-file', 'add-content', 'ac', 'tee', 'tee-object'}
READ_ONLY = {'rg', 'grep', 'egrep', 'fgrep', 'findstr', 'select-string', 'sls', 'cat', 'type', 'get-content', 'gc',
             'less', 'more', 'head', 'tail', 'ls', 'dir', 'get-childitem', 'gci', 'wc', 'echo', 'write-output',
             'write-host', 'printf', 'test-path', 'which', 'where', 'get-command', 'git-read'}
WRAPPERS = {'sudo', 'command', 'builtin', 'exec', 'nohup', 'time', 'env', 'do', 'then', 'else', 'elif', '!', '&',
            'call', 'start', 'start-process', 'iex', 'invoke-expression', 'cmd', '/c', '//c', '/k', '//k', '/s',
            '//s', '/q', '//q', 'powershell', 'pwsh', 'bash', 'sh', 'zsh', '-c', '-lc', '-command', '-noprofile',
            '-nologo', '-noninteractive', 'npx'}
VALUE_PARAMETERS = {'-erroraction', '-ea', '-warningaction', '-wa', '-informationaction', '-ia', '-errorvariable',
                    '-ev', '-outvariable', '-ov', '-filter', '-include', '-exclude', '-credential', '-stream',
                    '-executionpolicy', '-encoding', '-value'}
PATH_PARAMETERS = {'-path', '-literalpath', '-lp', '-pspath', '-filepath'}
DESTINATION_PARAMETERS = {'-destination', '-t', '--target-directory'}
MESSAGE = re.compile(r'''(\s(?:-[a-z]*m|--message)[ =]?\s*)("[^"]*"|'[^']*')''', re.DOTALL)
REDIRECT = re.compile(r'''(?<![<>&\d=-])\d?>>?\s*("[^"]+"|'[^']+'|[^\s;|&<>()]+)''')


def block(reason):
    sys.stderr.write('BLOCKED by scripts/uwtools/guard_catastrophe.py: ' + reason + '\n'
                     'This guard protects drives, system folders, the home folder, .git and the user\'s data trees, '
                     'and accepts only simple deletes with literal paths. Do not retry, rephrase or work around it. '
                     'Delete with one literal command (rm / Remove-Item and full paths, no pipeline or variables), '
                     'move untracked data with scripts/uwtools/trash.py where it exists, or tell the user what you '
                     'wanted to do and why.\n')
    sys.exit(2)


def split_command(command):
    """[(segment, separator before it, quote left open)] on ; && || | outside quotes and on every newline.

    Backslash and backtick are not taken as escapes (PowerShell's "C:\\" ends a string, bash's \\" does not): a
    quote may close early, which only splits more. A newline always ends a segment and resets the quote.
    """
    segments, current, before, quote, index = [], [], '', None, 0

    def close(separator):
        nonlocal current, before
        segments.append((''.join(current).strip(), before, quote is not None))
        current, before = [], separator

    while index < len(command):
        char = command[index]
        if char in '\r\n':
            close('\n')
            quote = None
        elif quote:
            if char == quote:
                quote = None
            current.append(char)
        elif char in '"\'':
            quote = char
            current.append(char)
        elif command[index:index + 2] in ('&&', '||'):
            close(command[index:index + 2])
            index += 1
        elif char in ';|':
            close(char)
        else:
            current.append(char)
        index += 1
    close('')
    return [item for item in segments if item[0]]


def clean_components(path):
    """Windows ignores trailing dots and spaces in a name (C:\\dev\\meta. is C:\\dev\\meta)."""
    drive, rest = os.path.splitdrive(path)
    parts = re.split(r'([\\/])', rest)
    return drive + ''.join(part if part in ('.', '..', '\\', '/', '') else part.rstrip(' .') or part for part in parts)


def expand(token, cwd):
    """A native, normalised path for a literal token, or None when it is not a literal path."""
    token = re.sub(r'["\']', '', token)
    if token[:4] in ('\\\\?\\', '//?/'):  # extended-length prefix: \\?\C:\x is C:\x
        token = token[4:]
    if token.startswith('~'):
        token = os.path.expanduser(token)
    if not token or re.search(r'[$`%{}()]|::|@', token) or re.match(r'^[A-Za-z]{2,}:', token):
        return None  # variables, braces, subexpressions, provider paths (HKCU:, Env:, FileSystem::)
    if token in ('\\', '`'):
        return None  # a line continuation
    match = re.match(r'^/([a-zA-Z])(/.*)?$', token)
    if match:  # Git Bash drive paths: /c/dev -> C:/dev
        token = match.group(1).upper() + ':' + (match.group(2) or '/')
    elif token == '/tmp' or token.startswith('/tmp/'):  # Git Bash mounts /tmp on the user's temp folder
        token = tempfile.gettempdir() + token[4:]
    elif token.startswith('/') and not token.startswith('//'):
        token = 'C:/Program Files/Git' + token  # Git Bash's own root (/, /usr, /etc): protected
    if re.match(r'^[A-Za-z]:$', token):
        token += '\\'
    if cwd is None and not os.path.isabs(token):
        return None
    path = token if os.path.isabs(token) else os.path.join(cwd, token)
    return os.path.normcase(clean_components(os.path.normpath(path)))


def depth(path):
    drive, rest = os.path.splitdrive(path)
    return len([part for part in re.split(r'[\\/]', rest) if part])


def within(path, tree):
    tree = os.path.normcase(os.path.normpath(tree))
    return path == tree or path.startswith(tree.rstrip('\\/') + os.sep)


def forms(path):
    """The path, its resolved form (junctions, 8.3 names), each cut back to the folder above any wildcard."""
    candidates = {path}
    try:
        candidates.add(os.path.normcase(os.path.realpath(path)))
    except (OSError, ValueError):
        pass
    return {os.path.dirname(re.split(r'[*?]', item)[0] + 'x') if re.search(r'[*?]', item) else item
            for item in candidates}


def protected(path):
    """The reason a path may not be deleted or moved away, or None."""
    for item in forms(path):
        if item.startswith('\\\\') and depth(item) <= 3:
            return f'{item} is a network share or device root'
        if depth(item) <= 2:
            return f'{item} is a drive root or within two levels of one'
        if item == HOME or os.path.dirname(item) == HOME or within(HOME, item):
            return f'{item} is the home folder, one of its direct children or above it'
        if '.git' in re.split(r'[\\/]', item):
            return f'{item} is or is inside a .git directory'
        for tree in SYSTEM_TREES + PROTECTED_TREES:
            normal = os.path.normcase(os.path.normpath(tree))
            if within(item, normal) or within(normal, item):
                return f'{item} is or holds the protected tree {tree}'
        for folder in ROOT_ONLY:
            if within(os.path.normcase(os.path.normpath(folder)), item):
                return f'{item} is or holds the protected folder {folder} (its children may be deleted)'
    return None


def protected_destination(path):
    """The reason a path may not be written or copied into, or None."""
    for item in forms(path):
        if depth(item) == 0:
            return f'{item} is a drive root'
        for tree in DESTINATION_TREES:
            if within(item, os.path.normcase(os.path.normpath(tree))):
                return f'{item} lies in the protected tree {tree}'
    return None


def literal(token, cwd, what, rule, strict=True):
    """Check one target (commas split PowerShell lists). Strict: block a target that is not a literal path."""
    for part in (piece for piece in token.split(',') if piece):
        if strict and re.search(r'[*?]', os.path.dirname(part.replace('\\', '/'))):
            block(f'{what}: wildcards only in the last path component ({part!r}).')
        path = expand(part, cwd)
        if path is None:
            if strict:
                block(f'{what}: {part!r} is not a literal path the guard can resolve (variables, braces, a provider '
                      'path, or a relative path after a cd it could not follow); name full paths literally.')
            continue
        reason = rule(path)
        if reason:
            block(f'{what} of a protected location: {reason}.')


def arguments(words, head):
    """(positional, path parameter values, destination parameter values) of a command's words."""
    positional, paths, destinations = [], [], []
    index = 1
    while index < len(words):
        word, lowered = words[index], words[index].lower()
        if lowered in PATH_PARAMETERS | DESTINATION_PARAMETERS | VALUE_PARAMETERS:
            if index + 1 < len(words):
                if lowered in PATH_PARAMETERS:
                    paths.append(words[index + 1])
                elif lowered in DESTINATION_PARAMETERS:
                    destinations.append(words[index + 1])
            index += 2
            continue
        if lowered.startswith('--target-directory='):
            destinations.append(word.split('=', 1)[1])
        elif word.startswith('-') or (head in CMD_STYLE | MOVERS | COPIERS and re.fullmatch(r'//?[a-zA-Z]{1,2}', word)):
            pass
        else:
            positional.append(word)
        index += 1
    return positional, paths, destinations


SHELLS = {'cmd', 'powershell', 'pwsh', 'bash', 'sh', 'zsh', 'start-process', 'start', 'iex', 'invoke-expression'}


def name_of(word):
    """A command word's name: the file name of a path, lower case, without .exe or .com."""
    return os.path.basename(word.replace('\\', '/')).lower().removesuffix('.exe').removesuffix('.com')


def words_of(segment):
    """(words without leading wrappers, whether a shell wrapper was stripped: its quoted words are commands)."""
    try:
        words = shlex.split(segment, posix=False)
    except ValueError:
        words = segment.split()
    words = [w for w in (re.sub(r'^["\']|["\']$', '', w) for w in words) if w]
    shell = False
    while words:
        first = words[0].lower()
        readable = not re.search(r'\s', first) or re.search(r'\.(exe|com)$', first)  # a quoted program path
        if readable and (first in WRAPPERS or name_of(first) in WRAPPERS):
            shell = shell or name_of(first) in SHELLS
            words = words[1:]
        elif re.fullmatch(r'\w+=\S*', words[0]):
            words = words[1:]
        elif first in ('timeout', 'nice'):
            words = words[1:]
            while words and re.fullmatch(r'-\w+|\d+[smh]?', words[0]):
                words = words[1:]
        elif first in ('(', '{', ')', '}'):
            words = words[1:]
        elif first[:1] in ('(', '{'):
            words[0] = words[0][1:]
        else:
            break
    return words, shell


def check_segment(segment, cwd, piped, unbalanced=False, level=0):
    """Check one command; return the working folder after it (None when a cd could not be followed)."""
    words, shell = words_of(segment)
    if not words:
        return cwd
    head = name_of(words[0])
    lowered = [w.lower() for w in words]
    if head == 'git' and len(lowered) > 1 and lowered[1] in ('log', 'show', 'diff', 'grep', 'status', 'blame'):
        head = 'git-read'
    for match in REDIRECT.finditer(segment):
        target = match.group(1)
        if target.lower() not in ('/dev/null', '$null', 'nul') and not target.startswith('&'):
            literal(target, cwd, 'redirection', protected_destination, strict=False)
    deleting = bool(DELETE_WORDS.search(segment))
    if deleting and SUBSTITUTION.search(segment):
        block('a delete together with command substitution or backticks: the guard cannot see what runs.')
    if head in READ_ONLY and not unbalanced:
        return cwd
    if DISK_TOOLS.search(segment):
        block('disk, partition, volume, shadow-copy and registry-delete tools are never run by an agent.')
    if head in ('fdisk', 'parted') or (head == 'format' and any(re.fullmatch(r'[a-z]:\\?|//?fs:\S+', w) for w in lowered[1:])):
        block(f'{head} is never run by an agent.')
    if head == 'cipher' and any(re.fullmatch(r'//?w(:.*)?', w) for w in lowered):
        block('cipher /w wipes free space.')
    if head == 'dd' and any(w.startswith('of=') for w in lowered):
        block('dd with of= writes raw devices or files.')
    if deleting and segment.rstrip().endswith(('\\', '`')):
        block('a delete continued on the next line; write it on one line.')
    if shell:  # after cmd /c, bash -c, pwsh -c, Start-Process ...: every quoted word is a command of its own
        inner_words = [w for w in words if re.search(r'\s', w)]
        if inner_words and level >= 3:
            block('commands nested more than three deep.')
        for word in inner_words:
            inner_cwd = cwd
            for inner, separator, open_quote in split_command(word):
                inner_cwd = check_segment(inner, inner_cwd, separator == '|', open_quote, level + 1)
        if re.search(r'\s', words[0]):
            if deleting and any(not re.search(r'\s', w) and not w.startswith(('-', '/')) for w in words[1:]):
                block('a delete next to a quoted inner command: the guard cannot tell which runs.')
            return cwd
    elif re.search(r'\s', words[0]) and not re.search(r'\.(exe|com|bat|cmd|ps1)$', words[0], re.IGNORECASE):
        if deleting:  # quoted text that is not a program path (heredoc lines, JSON): readable only without deletes
            block(f'a delete in a line the guard cannot read as a command ({words[0][:60]!r}).')
        return cwd
    if head in ('cd', 'set-location', 'sl', 'pushd', 'chdir'):
        positional, paths, _ = arguments([w for w in words if w.lower() != '/d'], head)
        return expand((paths or positional or ['~'])[0], cwd)
    if head == 'popd':
        return None
    positional, paths, destinations = arguments(words, head)
    if head in DELETERS:
        if piped:
            block(f'{head} fed by a pipeline: the guard cannot see its targets; name paths literally.')
        targets = positional + paths
        if not targets:
            block(f'{head} without a target on its command line; name paths literally.')
        for target in targets:
            literal(target, cwd, 'delete', protected)
        return cwd
    if head in MOVERS or head in COPIERS or head in WRITERS:
        if piped and head in MOVERS:
            block(f'{head} fed by a pipeline: the guard cannot see its sources; name paths literally.')
        if head in WRITERS:
            sources = []
            targets = destinations + paths + positional if head.startswith('tee') else (destinations + paths
                                                                                      or positional[:1])
        else:
            named = paths + positional
            sources = named if destinations else named[:-1]
            targets = destinations or named[-1:]
        if head in MOVERS:
            if not sources:
                block(f'{head} without a literal source; name paths literally.')
            for source in sources:
                literal(source, cwd, 'move', protected)
        for target in targets:
            literal(target, cwd, f'{head} destination', protected_destination, strict=False)
        if deleting and head not in MOVERS:
            block(f'a delete inside {head!r}: the guard cannot see its targets.')
        return cwd
    if head == 'robocopy':
        switches = {w.lstrip('/-') for w in lowered[1:] if w.startswith(('/', '-'))}
        if switches & {'mir', 'purge', 'move', 'mov'}:
            places = [w for w in words[1:] if not w.startswith(('/', '-'))][:2]
            if len(places) < 2:
                block('robocopy with /MIR, /PURGE or /MOVE needs a literal source and destination.')
            for place in places:
                literal(place, cwd, 'robocopy mirror or move', protected)
        return cwd
    if head == 'rsync':
        if any(w.startswith(('--delete', '--remove-source')) for w in lowered):
            places = [w for w in words[1:] if not w.startswith('-')]
            if not places:
                block('rsync with deletes needs literal paths.')
            for place in places:
                literal(place, cwd, 'rsync with deletes', protected)
        return cwd
    if head == 'find':
        if '-delete' in lowered or deleting:
            starts = [w for w in words[1:] if not w.startswith('-')][:1] or ['.']
            literal(starts[0], cwd, 'find with a delete', protected)
        return cwd
    if head == 'cmake' and len(lowered) > 2 and lowered[1] == '-e' and lowered[2] in ('rm', 'remove', 'remove_directory'):
        targets = [w for w in words[3:] if not w.startswith('-')]
        if not targets:
            block('cmake -E rm without a target.')
        for target in targets:
            literal(target, cwd, 'cmake -E rm', protected)
        return cwd
    if head == 'git':
        if 'clean' in lowered and not {'-n', '--dry-run'} & set(lowered):
            if any(re.fullmatch(r'-[a-z]*x[a-z]*', w, re.IGNORECASE) for w in words[1:lowered.index('clean') + 6]):
                block('git clean -x / -X deletes ignored files, which often hold data; the user runs it.')
            where = cwd
            if '-c' in lowered[1:lowered.index('clean')]:
                index = lowered.index('-c')
                where = expand(words[index + 1], cwd) if index + 1 < len(words) else None
            if where is None:
                block('git clean in a folder the guard could not follow.')
            reason = protected_destination(where)
            if reason is None and (depth(where) <= 1 or where == HOME or within(HOME, where)):
                reason = f'{where} is a drive root, a top folder, the home folder or above it'
            if reason:
                block(f'git clean in a protected place: {reason}.')
        return cwd
    if deleting:
        block(f'a delete inside {head!r} (inline code, a loop or script block, xargs, or an unrecognised command): '
              'the guard cannot see its targets. Write one literal rm / Remove-Item command, or a script file.')
    return cwd


def check(command, cwd):
    command = MESSAGE.sub(r'\1"(message)"', command)  # a commit message is text, not a command
    for segment, separator, open_quote in split_command(command):
        cwd = check_segment(segment, cwd, separator == '|', open_quote)


def main():
    try:
        payload = json.loads(sys.stdin.buffer.read().decode('utf-8', errors='replace'))
        command = (payload.get('tool_input') or {}).get('command') or ''
        start = payload.get('cwd') or os.getcwd()
        check(command, os.path.normcase(os.path.normpath(start)))
    except SystemExit:
        raise
    except BaseException as error:  # fail closed: a guard that cannot decide blocks
        block(f'the guard could not check this command ({type(error).__name__}: {error}).')
    sys.exit(0)


if __name__ == '__main__':
    main()
