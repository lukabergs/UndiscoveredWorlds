"""PreToolUse guard (Bash and PowerShell, every project): refuse commands that could wipe a drive or a protected tree.

Added 2026-10-20 at the user's request ("people have complained about their AI agents trying to delete whole
drives"). Registered in ~/.claude/settings.json before delete_guard.py, so it runs in every project. It is
independent of delete_guard.py (which decides what git can restore), needs no git, and fails closed: any error
inside the guard blocks the command (Claude Code treats only exit code 2 as a block; an uncaught exception would
let it pass). The user asked that it never be tested by attempting deletes: change it by reading, and have a
reviewer read the change.

Blocks, before the command runs:
1. disk, partition, volume and shadow-copy tools: format, diskpart, Clear-Disk, Format-Volume, Remove-Partition,
   Initialize-Disk, mkfs, wipefs, fdisk, parted, dd of=, cipher /w, sdelete, vssadmin/wbadmin/wmic deletes, bcdedit;
2. a delete, a mirror or a move whose target is protected (`protected()`): a drive root or any path within two levels
   of one (C:\\, C:\\dev, C:\\dev\\project, C:\\Users\\name), the home folder and its direct children, Windows,
   Program Files, ProgramData, any `.git` directory, and the trees in PROTECTED_TREES (whole subtree);
3. a delete whose targets cannot be read from the command line: a deleter fed by a pipeline or xargs, a target
   with variables the guard cannot expand.
Segments are split on ; && || | and newlines, so inline code (python -c, heredocs, pwsh -Command) is checked line
by line: a segment that names a delete API (rmtree, os.remove, Remove-Item, rd /s, ...) has every absolute path
literal in it checked against the protected set.
Nothing here deletes or runs anything; it reads the hook JSON on stdin and exits 0 (allow) or 2 (block).
"""
import json
import os
import re
import shlex
import sys
import tempfile

HOME = os.path.normcase(os.path.normpath(os.path.expanduser('~')))
# Whole subtrees no agent command may delete, move or mirror onto (the user's data and backups).
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

DISK_TOOLS = re.compile(
    r'\b(diskpart|clear-disk|format-volume|remove-partition|initialize-disk|clear-volume|remove-physicaldisk|'
    r'reset-physicaldisk|set-partition|resize-partition|mkfs(\.\w+)?|wipefs|fdisk|sfdisk|cfdisk|parted|sgdisk|'
    r'sdelete(64)?|bcdedit|manage-bde)\b', re.IGNORECASE)
SHADOW_DELETE = re.compile(r'\b(vssadmin|wbadmin)\b.*\bdelete\b|\bwmic\b.*\bshadowcopy\b.*\bdelete\b', re.IGNORECASE)
DELETERS = {'rm', 'rm.exe', 'rmdir', 'rmdir.exe', 'rd', 'del', 'erase', 'unlink', 'shred', 'remove-item', 'ri',
            'rimraf', 'trash', 'clear-recyclebin'}
MOVERS = {'mv', 'mv.exe', 'move', 'move-item', 'mi', 'ren', 'rename', 'rename-item', 'rni'}
DELETE_API = re.compile(
    r'rmtree|os\.remove\b|os\.unlink|\.unlink\(|os\.rmdir|\.rmdir\(|removedirs|send2trash|::delete\(|\.delete\(\)|'
    r'remove-item|\brm\s+-|\brd\s+/s|\brmdir\s+/s|\bdel\s+/|\berase\s+/|deletedirectory|deletefile',
    re.IGNORECASE)
SEPARATORS = re.compile(r'&&|\|\||[;|\n]')
# Absolute path literals inside a segment that names a delete API: drive paths (C:\x, C:/x, C:), Git Bash drive
# paths with at least one component (/c/x; a bare /s is a cmd flag), home paths (~/x) and UNC paths (\\host\x).
ABSOLUTE = re.compile(r'''(?:^|(?<=[\s'"(=,\[]))(?:[A-Za-z]:(?:[\\/][^'"\s;|&<>)\]]*)?|/[a-zA-Z]/[^'"\s;|&<>)\]]*|'''
                      r'''~(?:[\\/][^'"\s;|&<>)\]]*)?|\\\\[^'"\s;|&<>)\]]+)(?=$|[\s'"),;\]])''')
CMD_STYLE = {'rd', 'del', 'erase', 'rmdir', 'rmdir.exe'}  # their /s /q /f are flags; rm's /c is a path
WRAPPERS = {'sudo', 'command', 'cmd', 'cmd.exe', '/c', '/k', '&', 'call', 'time', 'nohup', 'exec', 'start',
            'start-process', 'invoke-expression', 'iex', 'powershell', 'powershell.exe', 'pwsh', 'pwsh.exe',
            '-command', '-c', '-noprofile', '-nologo', '-noninteractive'}


def block(reason):
    sys.stderr.write('BLOCKED by ~/.claude/hooks/guard_catastrophe.py: ' + reason + '\n'
                     'This guard protects drives, system folders, the home folder, .git and the user\'s data trees. '
                     'Do not retry, rephrase or work around it; tell the user what you wanted to do and why.\n')
    sys.exit(2)


def expand(token, cwd):
    """A native, normalised path for a token, or None when it cannot be resolved (variables left)."""
    token = token.strip('"\'')
    if token[:4] in ('\\\\?\\', '//?/'):  # extended-length prefix: \\?\C:\x is C:\x
        token = token[4:]
    token = re.sub(r'\$env:(\w+)', lambda m: os.environ.get(m.group(1), '$UNRESOLVED'), token, flags=re.IGNORECASE)
    token = re.sub(r'\$\{?(HOME|USERPROFILE|TEMP|TMP|LOCALAPPDATA|APPDATA)\}?',
                   lambda m: os.environ.get(m.group(1), os.path.expanduser('~') if m.group(1) == 'HOME' else '$X'), token)
    token = os.path.expandvars(token)
    if token.startswith('~'):
        token = os.path.expanduser(token)
    if re.search(r'[$`%]|\$\(|@\(', token):
        return None
    match = re.match(r'^/([a-zA-Z])(/.*)?$', token)
    if match:  # Git Bash drive paths: /c/dev -> C:/dev
        token = match.group(1).upper() + ':' + (match.group(2) or '/')
    elif token == '/tmp' or token.startswith('/tmp/'):  # Git Bash mounts /tmp on the user's temp folder
        token = tempfile.gettempdir() + token[4:]
    elif token.startswith('/') and not token.startswith('//'):
        token = 'C:/Program Files/Git' + token  # Git Bash's own root (/, /usr, /etc): protected
    if re.match(r'^[A-Za-z]:$', token):
        token += '\\'
    path = token if os.path.isabs(token) else os.path.join(cwd, token)
    return os.path.normcase(os.path.normpath(path))


def depth(path):
    drive, rest = os.path.splitdrive(path)
    return len([part for part in re.split(r'[\\/]', rest) if part])


def within(path, tree):
    tree = os.path.normcase(os.path.normpath(tree))
    return path == tree or path.startswith(tree.rstrip('\\/') + os.sep)


def protected(path):
    """The reason a path may not be deleted, moved or mirrored onto, or None."""
    candidates = {path}
    try:
        candidates.add(os.path.normcase(os.path.realpath(path)))
    except (OSError, ValueError):
        pass
    for item in candidates:
        if '*' in item or '?' in item:
            item = os.path.dirname(item.split('*')[0].split('?')[0] + 'x')
        if item.startswith('\\\\') and depth(item) <= 3:
            return f'{item} is a network share or device root'
        if depth(item) <= 2:
            return f'{item} is a drive root or within two levels of one'
        if item == HOME or os.path.dirname(item) == HOME or within(HOME, item):
            return f'{item} is the home folder or one of its direct children'
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


def words_of(segment):
    try:
        words = shlex.split(segment, posix=False)
    except ValueError:
        words = segment.split()
    words = [w.strip('"\'') for w in words if w.strip('"\'')]
    while words and words[0].lower() in WRAPPERS:
        words = words[1:]
    return words


def check_targets(tokens, cwd, what, rule=None):
    rule = rule or protected
    for token in (part for item in tokens for part in item.split(',') if part):  # PowerShell lists: a,b
        path = expand(token, cwd)
        if path is None:
            block(f'cannot resolve the {what} target {token!r} to a path; name paths literally.')
        reason = rule(path)
        if reason:
            block(f'{what} of a protected location: {reason}.')


def protected_tree(path):
    """The narrower rule for commands that act inside a folder (git clean): system, home and data trees only."""
    if path == HOME or within(HOME, path):
        return f'{path} is the home folder or above it'
    for tree in SYSTEM_TREES + PROTECTED_TREES:
        normal = os.path.normcase(os.path.normpath(tree))
        if within(path, normal) or within(normal, path):
            return f'{path} is or holds the protected tree {tree}'
    return None


def check_segment(segment, cwd, level=0):
    if DISK_TOOLS.search(segment) or SHADOW_DELETE.search(segment):
        block('disk, partition, volume or shadow-copy tools are never run by an agent.')
    words = words_of(segment)
    if not words:
        return cwd
    if level < 3:  # a quoted inner command (cmd /c "rd /s x", pwsh -c "...") is checked as its own command
        for word in words:
            if re.search(r'\s', word):
                for inner in SEPARATORS.split(word):
                    if inner.strip():
                        check_segment(inner.strip(), cwd, level + 1)
    head = os.path.basename(words[0].replace('\\', '/')).lower()
    lowered = [w.lower() for w in words]
    if head in ('cd', 'set-location', 'sl', 'pushd', 'chdir') and len(words) > 1:
        target = next((w for w in words[1:] if not w.startswith('-')), '')
        resolved = expand(target, cwd) if target else None
        return resolved or cwd
    if head in ('format', 'format.com'):
        block('format is never run by an agent.')
    if head in ('cipher', 'cipher.exe') and any(w.startswith('/w') for w in lowered):
        block('cipher /w wipes free space.')
    if head == 'dd' and any(w.startswith('of=') for w in lowered):
        block('dd with of= writes raw devices or files.')
    if head == 'xargs' and any(os.path.basename(w).lower() in DELETERS for w in words[1:]):
        block('xargs feeds delete targets the guard cannot read.')
    if head in DELETERS:
        targets = [w for w in words[1:] if not w.startswith('-')
                   and not (head in CMD_STYLE and re.fullmatch(r'/[a-zA-Z]{1,2}', w))]
        if not targets:
            block(f'{head} without a target on the command line (fed by a pipeline); name paths literally.')
        check_targets(targets, cwd, 'delete')
    if head in MOVERS:
        places = [w for w in words[1:] if not w.startswith('-') and not re.fullmatch(r'/[a-zA-Z]{1,2}', w)]
        if not places:
            block(f'{head} without a source on the command line (fed by a pipeline); name paths literally.')
        check_targets(places[:-1] or places, cwd, 'move')
    if head in ('git', 'git.exe') and 'clean' in lowered:
        where = cwd
        if '-c' in lowered[:lowered.index('clean')]:
            index = lowered.index('-c')
            where = expand(words[index + 1], cwd) if index + 1 < len(words) else None
        dry = {'-n', '--dry-run'} & set(lowered)
        if not dry and any(re.fullmatch(r'-[a-z]*x[a-z]*', w, re.IGNORECASE) for w in words[1:]):
            block('git clean -x / -X deletes ignored files, which often hold data; the user runs it.')
        if not dry:
            check_targets([where or '$UNRESOLVED'], cwd, 'git clean', protected_tree)
    if head in ('robocopy', 'robocopy.exe'):
        flags = {w.lower() for w in words if w.startswith('/')}
        if flags & {'/mir', '/purge', '/move', '/mov'}:
            places = [w for w in words[1:] if not w.startswith('/')][:2]
            if len(places) < 2:
                block('robocopy with /MIR, /PURGE or /MOVE needs a literal source and destination.')
            check_targets(places, cwd, 'robocopy mirror or move')
    if head in ('rsync', 'rsync.exe') and any(w.startswith('--delete') or w.startswith('--remove-source') for w in lowered):
        places = [w for w in words[1:] if not w.startswith('-')]
        check_targets(places[-1:] if places else ['$UNRESOLVED'], cwd, 'rsync --delete')
    if head == 'find' and ('-delete' in lowered or any(os.path.basename(w) in DELETERS for w in lowered[1:])):
        starts = [w for w in words[1:] if not w.startswith('-')][:1] or ['.']
        check_targets(starts, cwd, 'find -delete')
    if DELETE_API.search(segment):
        for match in ABSOLUTE.finditer(segment):
            path = expand(match.group(0), cwd)
            reason = protected(path) if path else None
            if reason:
                block(f'a delete, move or mirror in this command names a protected location: {reason}.')
    return cwd


def main():
    try:
        payload = json.load(sys.stdin)
        command = (payload.get('tool_input') or {}).get('command') or ''
        cwd = os.path.normcase(os.path.normpath(payload.get('cwd') or os.getcwd()))
        for segment in SEPARATORS.split(command):
            if segment.strip():
                cwd = check_segment(segment.strip(), cwd)
    except SystemExit:
        raise
    except BaseException as error:  # fail closed: a guard that cannot decide blocks
        block(f'the guard could not check this command ({type(error).__name__}: {error}).')
    sys.exit(0)


if __name__ == '__main__':
    main()
