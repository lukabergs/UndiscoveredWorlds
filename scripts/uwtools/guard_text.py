"""What in a command runs and what is only text, for the delete guards.

Used by guard_catastrophe.py and delete_guard.py (and ~/.claude/hooks/guard_recursive_delete.py). The user asked on
2026-11-01 that words such as rm, del or erase in a commit message, an echo, a search pattern, a comment or a heredoc
of file content no longer count as deletes, and that a < or > inside text no longer counts as a redirection.

A command is inert (its quoted strings, # comments, heredoc bodies and PowerShell here-string bodies are text) only
when every command in it is one that cannot run text: a head in TEXT_HEADS (git outside RUNNING_GIT, gh, echo, rg,
cat, the file commands, cmake but not -E env/time/chdir, ...) or an interpreter running a script file (python x.py,
uv run python x.py; not -c, -Ic, -e, -, <, | python). A script's own options are not looked into (x.py --exec "..."
is text). Besides, it is not inert when it holds:
- a shell, eval or launcher anywhere (bash, sh, cmd, pwsh, eval, xargs, Invoke-Expression, Start-Process, ...), a call
  of a string or variable (& "...", & $x, . $x, .Invoke(), [scriptblock]), git -c, or a bash $x command word;
- outside quotes: $(, and in bash a backtick, <( or >(, in PowerShell any ( (what they run could hide in a blanked
  string);
- a heredoc or here-string whose body expands (unquoted delimiter, @"...") and holds $( or a backtick;
- a heredoc or here-string together with an interpreter (the body may be the script it runs);
- quotes the scanner cannot pair: an escaped quote, a quote left open, curly quotes, a heredoc without its end line.
A non-inert command is read whole by the guards, as before. In an inert command a double-quoted string holding $(
(in bash also a backtick) still counts, and so does a quoted command word.
Keep this module small and readable: both guards depend on it, and hook_runner.py blocks when it fails.
"""
import os
import re

UNPAIRABLE = re.compile('[\\\\`]["\']|[\u201c\u201d\u201e\u2018\u2019\u201a]')  # escaped and curly quotes
HEREDOC = re.compile(r'(?<!<)<<(?!<)-?\s*(["\']?)([A-Za-z_][\w.-]*)\1')
HERESTRING = re.compile(r'@(["\'])\s*$')
RUNS_TEXT = re.compile(
    r'(?<![\w.-])(bash|sh|zsh|dash|ksh|fish|cmd|pwsh|powershell|wsl|ssh|su|eval|source|xargs|parallel|watch|iex|'
    r'invoke-expression|invoke-command|icm|start-process|saps|start|start-job|sajb|start-threadjob)(\.exe)?(?![\w.-])|'
    r'\[scriptblock\]|newscriptblock|invokecommand|\.invoke\w*\s*\(|addscript|(^|[\s;|&({])[&.]\s*["\'$]|'
    r'(?-i:(?<![\w-])git(\.exe)?\s(.*\s)?-c\s)', re.IGNORECASE)
BASH_VARIABLE_COMMAND = re.compile(r'(^|[\n;|&(`]|\b(then|do|else)\b)\s*\$\{?\w')
INTERPRETERS = {'python', 'python3', 'pythonw', 'py', 'uv', 'node'}
INTERPRETER = re.compile(r'(\|\s*)?(?<![\w.-])(python[\d.]*|pythonw|py|node|perl|ruby|deno|bun|php|lua|rscript|'
                         r'osascript)(\.exe)?(?![\w.-])([^;&|\n]*)', re.IGNORECASE)
CODE_INPUT = re.compile(r'\s-[A-Za-z]*[ceEpr]\b|\s--?(command|eval|print)\b|\s-(\s|$)|<', re.IGNORECASE)
TEXT_HEADS = {
    'git', 'gh', 'echo', 'printf', 'write-output', 'write-host', 'write-error', 'write-warning', 'write-verbose',
    'rg', 'grep', 'egrep', 'fgrep', 'findstr', 'select-string', 'sls', 'cat', 'type', 'tee', 'tee-object',
    'set-content', 'add-content', 'out-file', 'get-content', 'gc', 'ls', 'dir', 'get-childitem', 'gci', 'get-item',
    'head', 'tail', 'wc', 'sort', 'uniq', 'cut', 'diff', 'jq', 'du', 'df', 'pwd', 'date', 'which', 'where',
    'get-command', 'test-path', 'cd', 'set-location', 'sl', 'pushd', 'popd', 'push-location', 'pop-location', 'mkdir',
    'md', 'new-item', 'ni', 'touch', 'cp', 'copy', 'copy-item', 'mv', 'move', 'move-item', 'rm', 'rmdir', 'rd', 'del',
    'remove-item', 'ri', 'cmake', 'ctest', 'ninja', 'true', 'false', 'exit', 'sleep'}
RUNNING_GIT = {'config', 'rebase', 'submodule', 'bisect', 'filter-branch', 'filter-repo', 'difftool', 'mergetool'}
RUNNING_GH = {'alias', 'extension', 'ext'}


def split_command(command):
    """[(segment, separator before it, quote left open)] on ; && || | & outside quotes and on every newline.

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
        elif char == '&' and command[index - 1:index] not in ('>', '<') and command[index + 1:index + 2] != '>':
            close('&')  # cmd's and bash's lone & runs the next command too (2>&1, &> and >& are redirections)
        else:
            current.append(char)
        index += 1
    close('')
    return [item for item in segments if item[0]]


def scan(text):
    """([(start, end, kind)], start of a quote left open or None): quoted strings (kind " or ') and # comments
    (kind #, a # that starts a word, to the end of its line)."""
    spans, index = [], 0
    while index < len(text):
        char = text[index]
        if char in '"\'':
            end = text.find(char, index + 1)
            if end < 0:
                return spans, index
            spans.append((index, end + 1, char))
            index = end + 1
        elif char == '#' and (index == 0 or text[index - 1] in ' \t\n;|&('):
            end = text.find('\n', index)
            end = len(text) if end < 0 else end
            spans.append((index, end, '#'))
            index = end
        else:
            index += 1
    return spans, None


def runs(body, powershell):
    """True when a quoted string runs part of itself: "$(...)", or a backtick in bash."""
    return body.startswith('"') and ('$(' in body or (not powershell and '`' in body))


def blank(text, powershell=False, redirects=False, everything=False):
    """The text with its comments and the quoted strings that are text emptied (a string that runs, a quoted
    command word and, with redirects, a redirection's quoted target stay; with everything, none stays). Unchanged
    when the quotes cannot be paired."""
    if UNPAIRABLE.search(text):
        return text
    spans, open_at = scan(text)
    if open_at is not None:
        return text
    out, last = [], 0
    for start, end, kind in spans:
        out.append(text[last:start])
        prefix = text[:start].rstrip(' \t')
        if kind == '#':
            out.append(' ')
        elif not everything and (runs(text[start:end], powershell) or not prefix or prefix[-1] in '\n;|&(){}`' or (
                redirects and prefix[-1] == '>')):
            out.append(text[start:end])
        else:
            out.append(kind * 2)
        last = end
    out.append(text[last:])
    return ''.join(out)


def heredocs(before, line):
    """[(end word, quoted)] of the heredocs a bash line opens (not those inside quoted text or a comment); before is
    the command's text up to the line."""
    found = []
    for match in HEREDOC.finditer(line):
        text = before + line[:match.start()]
        spans, open_at = scan(text)
        if open_at is not None and not (text[open_at] == '"' and '$(' in text[open_at:]):
            continue  # inside quoted text (a heredoc inside "$(...)" is real)
        if spans and spans[-1][2] == '#' and spans[-1][1] >= len(text):
            continue  # inside a comment
        found.append((match.group(2), bool(match.group(1))))
    return found


def head_of(segment):
    """(command name, its words): after leading assignments, without brackets, quotes, a path or .exe."""
    words = segment.split()
    while words and (re.fullmatch(r'\w+=\S*', words[0]) or re.fullmatch(r'[(){}]+', words[0])):
        words = words[1:]  # assignments and the brackets of ( ... ) and { ...; } groups
    first = words[0].strip('(){}"\'@') if words else ''
    return os.path.basename(first.replace('\\', '/')).lower().removesuffix('.exe'), words


def text_head(segment, ends):
    """True when a segment cannot run text: its command is in TEXT_HEADS (git and gh outside their running
    subcommands) or an interpreter running a script file; heredoc end lines and closing brackets count too."""
    head, words = head_of(segment)
    if not head or segment.strip() in ends:
        return True
    if head in ('git', 'gh'):
        index, sub = 1, ''
        while index < len(words):
            if words[index] in ('-C', '--git-dir', '--work-tree', '--namespace'):
                index += 2
            elif words[index].startswith('-'):
                index += 1
            else:
                sub = words[index].lower()
                break
        return sub not in (RUNNING_GIT if head == 'git' else RUNNING_GH)
    if head == 'uv':  # uv run runs any program: only python or a .py script counts
        program = next((word for word in words[2:] if not word.startswith('-')), '') if words[1:2] == ['run'] else ''
        if not (re.fullmatch(r'python[\d.]*(\.exe)?', program.lower()) or program.lower().endswith('.py')):
            return False
    if head in INTERPRETERS:
        return not CODE_INPUT.search(' ' + ' '.join(words[1:]))
    if head == 'cmake' and len(words) > 2 and words[1] == '-E' and words[2].lower() in ('env', 'time', 'chdir'):
        return False  # they run the command that follows
    return head in TEXT_HEADS


def readable(command, powershell=False):
    """(command, inert). inert is True when nothing in the command can run text (see the module's rules): the
    command then comes back without the bodies of its heredocs (bash) and here-strings (PowerShell), with the line
    breaks inside its text strings made spaces (a quoted string is one argument), and blank() may be applied to it.
    Otherwise it comes back unchanged."""
    lines = command.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    kept, ends, index, dropped = [], set(), 0, False
    while index < len(lines):
        line = lines[index]
        if powershell:
            opener = HERESTRING.search(line)
            found = [(opener.group(1) + '@', opener.group(1) == "'")] if opener else []
        else:
            found = heredocs(''.join(item + '\n' for item in kept), line) if '<<' in line else []
        kept.append(line)
        index += 1
        for name, literal in found:
            closes = (lambda text: text.startswith(name)) if powershell else (lambda text: text.strip() == name)
            end = next((i for i in range(index, len(lines)) if closes(lines[i])), None)
            if end is None:
                return command, False
            body = '\n'.join(lines[index:end])
            if not literal and ('$(' in body or (not powershell and '`' in body)):
                return command, False  # an expanding body runs its substitutions
            kept.append(lines[end])
            ends.add(lines[end].strip())
            index, dropped = end + 1, True
    reduced = '\n'.join(kept)
    spans, open_at = scan(reduced)
    if open_at is not None or UNPAIRABLE.search(reduced):
        return command, False
    probe = blank(reduced, powershell)
    bare = blank(reduced, powershell, everything=True)
    if RUNS_TEXT.search(probe) or '$(' in bare or (powershell and '(' in bare) or (not powershell and (
            '`' in bare or '<(' in bare or '>(' in bare or BASH_VARIABLE_COMMAND.search(
            probe))) or any(match.group(1) or CODE_INPUT.search(match.group(4)) for match in INTERPRETER.finditer(probe)):
        return command, False
    out, last = [], 0
    for start, end, kind in spans:
        if kind != '#' and not runs(reduced[start:end], powershell):
            out.append(reduced[last:start] + reduced[start:end].replace('\n', ' '))
            last = end
    out.append(reduced[last:])
    joined = ''.join(out)
    segments = [segment for segment, _, _ in split_command(blank(joined, powershell))]
    if not all(text_head(segment, ends) for segment in segments):
        return command, False
    if dropped and any(head_of(segment)[0] in INTERPRETERS for segment in segments):
        return command, False  # a heredoc body may be the script an interpreter runs
    return joined, True
