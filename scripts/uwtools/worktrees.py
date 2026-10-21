"""Agent worktrees in .claude/worktrees: audit, preserve to storage, retire.

    uv run python scripts/uwtools/worktrees.py audit [--json FILE]
    uv run python scripts/uwtools/worktrees.py preserve NAME... | --all
    uv run python scripts/uwtools/worktrees.py retire NAME... | --all --yes
    uv run python scripts/uwtools/worktrees.py drop-branches NAME... | --retired

audit (read-only): per worktree its links (never followed), whether its commits or uncommitted diff reached
main (the share of added lines found in some version of the file in main's history since the worktree's
base), its untracked files and sizes. preserve (copies only): writes <storage root>/worktree-archive/<name>/
with state.json (head, branch, status), uncommitted.patch (git diff --binary HEAD) and every untracked or
ignored file except .venv/, out/build/ and caches, each copy SHA-256-verified; a worktree holding a link is
refused. retire: deletes a preserved worktree's directory with storage.remove_tree (no delete when any link
is inside), only if nothing changed since the archive was written, then `git worktree prune`; the branch is
kept, so its commits stay reachable. Never `git worktree remove --force` (2026-10-07: it followed a junction
into refs/); remove_tree makes read-only files and folders writable. A retire that stopped after the delete
began (the worktree's .git file already gone, 2026-10-20) is finished by running retire again: see
resume_retire for what it checks first.
"""
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import storage  # noqa: E402
from uwtools.hashing import sha256  # noqa: E402
from uwtools.paths import ROOT  # noqa: E402

WORKTREES = ROOT / '.claude' / 'worktrees'
# Regenerable: not archived, and changes inside them do not block a retire.
SKIPPED = ('.venv/', 'out/build/', '__pycache__/', '.pytest_cache/', '.ruff_cache/')


def git(*args, cwd=ROOT, binary=False):
    result = subprocess.run(['git', *args], cwd=cwd, capture_output=True)
    if binary:
        return result.returncode, result.stdout
    return result.returncode, result.stdout.decode('utf-8', 'replace').strip()


def skipped(relative):
    return any(relative.startswith(prefix) or f'/{prefix}' in relative for prefix in SKIPPED)


def archive_root():
    root = storage.storage_root()
    if root is None:
        raise OSError('No storage configured (runs/storage.json)')
    return root / 'worktree-archive'


def landed(worktree, head, uncommitted=None):
    """(added lines, found in main's history, the first missing ones, description of the change).

    The uncommitted diff is read from the worktree unless given (a worktree without its .git file: git run
    there would read the main repository instead).
    """
    in_main = git('merge-base', '--is-ancestor', head, 'main')[0] == 0
    base = head if in_main else git('merge-base', head, 'main')[1]
    diff = '' if in_main else git('diff', base, head)[1]
    describe = 'uncommitted diff' if in_main else 'commits ' + '; '.join(git('log', '--format=%h %s', f'{base}..{head}')[1].splitlines())
    if uncommitted is None:
        uncommitted = git('diff', 'HEAD', cwd=worktree)[1]
    if uncommitted and not in_main:
        describe += ' + uncommitted diff'
    diff += '\n' + uncommitted
    added, path = {}, None
    for line in diff.splitlines():
        if line.startswith('+++ '):
            path = line[6:] if line.startswith('+++ b/') else None
        elif path and line.startswith('+') and len(line[1:].strip()) > 3:
            added.setdefault(path, []).append(line[1:].strip())
    total = found = 0
    missing = []
    for path, lines in added.items():
        versions = set()
        for commit in git('rev-list', f'{base}..main', '--', path)[1].split() + ['main']:
            versions.update(text.strip() for text in git('show', f'{commit}:{path}')[1].splitlines())
        for text in lines:
            total += 1
            if text in versions:
                found += 1
            else:
                missing.append(f'{path}: {text[:100]}')
    return total, found, missing, describe


def sizes(worktree):
    """Bytes per top-level entry (per out/<child> under out), links not followed."""
    totals = {}
    root = storage.long_path(worktree)
    stack = [str(root)]
    while stack:
        try:
            entries = list(os.scandir(stack.pop()))
        except OSError:
            continue
        for entry in entries:
            if storage.is_link(entry):
                continue
            if entry.is_dir(follow_symlinks=False):
                stack.append(entry.path)
                continue
            parts = Path(entry.path).relative_to(root).parts
            key = '/'.join(parts[:2]) if parts[0] == 'out' and len(parts) > 2 else parts[0]
            totals[key] = totals.get(key, 0) + entry.stat(follow_symlinks=False).st_size
    return totals


def half_retired(worktree):
    """True when the folder has no .git, or a .git file naming a gitdir that no longer exists: git run there would
    read the main repository (the folder lies inside it), so it is never run there."""
    dot_git = storage.long_path(worktree) / '.git'
    if not os.path.lexists(dot_git):
        return True
    if storage.is_link(dot_git) or not dot_git.is_file():
        return False
    text = dot_git.read_text(encoding='utf-8', errors='replace').strip()
    if not text.startswith('gitdir:'):
        return False
    gitdir = Path(text[len('gitdir:'):].strip())
    return not (gitdir if gitdir.is_absolute() else worktree / gitdir).is_dir()


def audit_one(worktree):
    receipt = archive_root() / worktree.name / 'receipt.json'
    if not worktree.is_dir() or half_retired(worktree):
        # The folder itself, or its .retiring quarantine when only that is left.
        folder = worktree if worktree.is_dir() else worktree.with_name(worktree.name + '.retiring')
        totals = sizes(folder)
        archived = json.loads(receipt.read_text(encoding='utf-8')) if receipt.is_file() else {}
        return {
            'name': worktree.name, 'head': archived.get('head'), 'branch': archived.get('branch'),
            'date': datetime.datetime.fromtimestamp(folder.stat().st_mtime).strftime('%Y-%m-%d %H:%M'),
            'change': f'no .git (half-retired): resume with retire {worktree.name} --yes', 'half_retired': True,
            'added_lines': 0, 'found_in_main': 0, 'missing': [], 'modified': [], 'untracked': [],
            'links': storage.find_links(folder),
            'gb': round(sum(totals.values()) / 1e9, 2),
            'gb_by_part': {key: round(value / 1e9, 2) for key, value in sorted(totals.items(), key=lambda kv: -kv[1]) if value > 2e7},
            'archive': archived.get('status'),
        }
    head = git('rev-parse', 'HEAD', cwd=worktree)[1]
    total, found, missing, describe = landed(worktree, head)
    status = git('status', '--porcelain', '--untracked-files=all', cwd=worktree)[1].splitlines()
    totals = sizes(worktree)
    return {
        'name': worktree.name, 'head': head, 'branch': git('rev-parse', '--abbrev-ref', 'HEAD', cwd=worktree)[1],
        'date': git('log', '-1', '--format=%ci', 'HEAD', cwd=worktree)[1], 'change': describe,
        'added_lines': total, 'found_in_main': found, 'missing': missing[:20],
        'modified': [line[3:] for line in status if not line.startswith('??')],
        'untracked': [line[3:] for line in status if line.startswith('??')],
        'links': storage.find_links(worktree),
        'gb': round(sum(totals.values()) / 1e9, 2),
        'gb_by_part': {key: round(value / 1e9, 2) for key, value in sorted(totals.items(), key=lambda kv: -kv[1]) if value > 2e7},
        'archive': json.loads(receipt.read_text(encoding='utf-8'))['status'] if receipt.is_file() else None,
    }


def preserve(worktree):
    if half_retired(worktree):
        raise OSError(f'{worktree.name}: no .git (half-retired); git there would read the main repository. '
                      f'Resume with retire {worktree.name} --yes')
    links = storage.find_links(worktree)
    if links:
        raise OSError(f'{worktree.name}: links inside, preserve and retire by hand after review: {links[:3]}')
    destination = archive_root() / worktree.name
    receipt_path = destination / 'receipt.json'
    if receipt_path.is_file():
        # Run again on an archive: add the files it lacks (2026-10-17: the first archives were written
        # without extended-length paths, so a file past 260 characters could have been skipped).
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if receipt['status'] != 'preserved':
            raise OSError(f'{worktree.name}: archive status {receipt["status"]}')
        digests = receipt['sha256']
    else:
        destination.mkdir(parents=True, exist_ok=True)
        state = audit_one(worktree)
        (destination / 'state.json').write_text(json.dumps(state, indent=1) + '\n', encoding='utf-8')
        patch = git('diff', '--binary', 'HEAD', cwd=worktree, binary=True)[1]
        (destination / 'uncommitted.patch').write_bytes(patch)
        receipt = {'worktree': str(worktree), 'head': state['head'], 'branch': state['branch'], 'status': 'preserved',
                   'patch_sha256': sha256(destination / 'uncommitted.patch')}
        digests = {}
    listed = git('ls-files', '--others', '-z', cwd=worktree, binary=True)[1].decode('utf-8', 'replace').split('\0')
    source_root, target_root = storage.long_path(worktree), storage.long_path(destination / 'files')
    added = 0
    for relative in sorted(item for item in listed if item and not skipped(item) and item not in digests):
        source = source_root / relative
        if not source.is_file() or storage.is_link(source):
            if not source.exists():
                raise OSError(f'{worktree.name}: listed by git but unreachable: {relative}')
            continue
        added += 1
        target = target_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        digests[relative] = sha256(source)
        if sha256(target) != digests[relative]:
            raise OSError(f'Copy verification failed: {source}')
    receipt['sha256'] = digests
    receipt['preserved_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    return len(digests), added, sum((source_root / item).stat().st_size for item in digests) / 1e9


def retire(worktree, allow_unlanded=False):
    receipt_path = archive_root() / worktree.name / 'receipt.json'
    if not receipt_path.is_file():
        raise OSError(f'{worktree.name}: not preserved; run preserve first')
    receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
    if receipt['status'] != 'preserved':
        raise OSError(f'{worktree.name}: archive status {receipt["status"]}')
    if not worktree.is_dir() or half_retired(worktree):
        return resume_retire(worktree, receipt, receipt_path, allow_unlanded)
    stamp = receipt_path.stat().st_mtime
    if git('rev-parse', 'HEAD', cwd=worktree)[1] != receipt['head']:
        raise OSError(f'{worktree.name}: HEAD moved since the archive')
    # Agents may retire on their own since 2026-10-20 (the user), but only a worktree whose work is all in main:
    # every added line of its commits and uncommitted diff found in main's history. The user can override.
    total, found, missing, _ = landed(worktree, receipt['head'])
    if found < total and not allow_unlanded:
        raise OSError(f'{worktree.name}: {total - found} of {total} added lines are not in main, e.g. {missing[:2]}; '
                      'merge them first (or the user retires it with --allow-unlanded)')
    if sha256_bytes(git('diff', '--binary', 'HEAD', cwd=worktree, binary=True)[1]) != receipt['patch_sha256']:
        raise OSError(f'{worktree.name}: the uncommitted diff changed since the archive')
    code, listing = git('ls-files', '--others', '-z', cwd=worktree, binary=True)
    if code != 0:
        raise OSError(f'{worktree.name}: git ls-files failed; nothing deleted')
    listed = [item for item in listing.decode('utf-8', 'replace').split('\0') if item and not skipped(item)]
    unarchived = [item for item in listed if item not in receipt['sha256']]
    if unarchived:
        raise OSError(f'{worktree.name}: {len(unarchived)} untracked files not in the archive, e.g. {unarchived[0]}; '
                      'run preserve again')
    root_path = storage.long_path(worktree)
    changed = [item for item in listed if sha256(root_path / item) != receipt['sha256'][item]]
    if changed:
        raise OSError(f'{worktree.name}: {len(changed)} untracked files changed since the archive, e.g. {changed[0]}; '
                      'run preserve again')
    root = storage.long_path(worktree)
    newer = [item for item in root.rglob('*') if item.is_file() and item.stat().st_mtime > stamp
             and not skipped(item.relative_to(root).as_posix()) and item.name != '.git']
    if newer:
        raise OSError(f'{worktree.name}: {len(newer)} files newer than the archive, e.g. {newer[0]}')
    # The head must stay reachable after the folder goes: on a branch or in main (a detached HEAD's commits
    # would otherwise be lost to git's garbage collection).
    in_main = git('merge-base', '--is-ancestor', receipt['head'], 'main')[0] == 0
    if not in_main and not git('branch', '--contains', receipt['head'])[1].strip() and not allow_unlanded:
        raise OSError(f'{worktree.name}: its HEAD {receipt["head"][:9]} is on no branch and not in main')
    if git('update-ref', f'refs/archive/{worktree.name}', receipt['head'])[0] != 0:  # keeps the commits reachable
        raise OSError(f'{worktree.name}: could not write refs/archive/{worktree.name}; nothing deleted')
    links = storage.find_links(worktree)
    if links:
        raise OSError(f'{worktree.name}: links inside ({links[0]}); nothing deleted')
    delete_folder(worktree, worktree)
    git('worktree', 'prune')
    receipt['status'] = 'retired'
    receipt['retired_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')


def delete_folder(folder, worktree):
    """Delete a worktree's folder (or its .retiring quarantine) through storage.remove_tree, never following a link.

    Renamed first (unless already the quarantine): the rename fails while a process still holds a file inside (a
    worker still running). remove_tree makes read-only files and folders writable itself. On failure the folder
    is renamed back to the worktree's name, for a later attempt (which resumes: resume_retire).
    """
    quarantine = worktree.with_name(worktree.name + '.retiring')
    if folder != quarantine:
        os.rename(storage.long_path(folder), storage.long_path(quarantine))
    try:
        storage.remove_tree(quarantine)
    except OSError:
        if quarantine.exists() and not worktree.exists():
            os.rename(storage.long_path(quarantine), storage.long_path(worktree))  # back, for a later attempt
        raise


def listed_worktree(path):
    """The lines of `git worktree list --porcelain` for path, or None when git does not list it."""
    wanted = os.path.normcase(os.path.abspath(path))
    for block in git('worktree', 'list', '--porcelain')[1].split('\n\n'):
        lines = block.splitlines()
        if lines and lines[0].startswith('worktree ') and os.path.normcase(os.path.abspath(lines[0][9:])) == wanted:
            return lines
    return None


def resume_retire(worktree, receipt, receipt_path, allow_unlanded=False):
    """Finish a retire that stopped after its delete began: the .git file is gone or names a missing gitdir
    (half_retired), or the whole folder is gone.

    git cannot be asked about such a folder (run there it reads the main repository), so everything is checked
    against the archive instead, and nothing is deleted unless all of it holds: the folder is a direct child of
    .claude/worktrees and the one the receipt names; the archive's patch and files match their digests; git lists
    no worktree there, or lists it as prunable; the archived HEAD exists, the receipt's branch (if it still
    exists) and any HEAD git still records for the worktree point at it, and its added lines (commits and the
    archived uncommitted patch) are all in main's history (or the user passes --allow-unlanded) and stay
    reachable; no link is inside; every file left outside the regenerable folders is an archived untracked file
    still matching its digest or a file tracked at the archived HEAD (whose content, hashed by git, is its blob
    there unless the archived patch touches it), and none is newer than the archive. Then
    `git worktree prune`, the rest of the folder is deleted (none left: nothing to delete) and the receipt marked
    retired.
    """
    name = worktree.name
    quarantine = worktree.with_name(name + '.retiring')
    if name in ('', '.', '..') or os.path.normcase(os.path.abspath(worktree.parent)) != os.path.normcase(os.path.abspath(WORKTREES)):
        raise OSError(f'{worktree}: not a folder directly in {WORKTREES}; nothing deleted')
    if os.path.normcase(os.path.abspath(receipt['worktree'])) != os.path.normcase(os.path.abspath(worktree)):
        raise OSError(f'{name}: the archive receipt names {receipt["worktree"]}; nothing deleted')
    archive = receipt_path.parent
    patch_path = archive / 'uncommitted.patch'
    if not patch_path.is_file() or sha256(patch_path) != receipt['patch_sha256']:
        raise OSError(f'{name}: the archived uncommitted.patch is missing or does not match its digest; nothing deleted')
    files_root = storage.long_path(archive / 'files')
    damaged = [item for item, digest in receipt['sha256'].items()
               if not (files_root / item).is_file() or sha256(files_root / item) != digest]
    if damaged:
        raise OSError(f'{name}: {len(damaged)} archived files missing or not matching, e.g. {damaged[0]}; nothing deleted')
    listed = listed_worktree(worktree)
    if listed is not None and not any(line.startswith('prunable') for line in listed):
        raise OSError(f'{name}: its .git file is gone but git still lists the worktree (not prunable: {listed[1:]}); '
                      'nothing deleted')
    head = receipt['head']
    if git('cat-file', '-e', f'{head}^{{commit}}')[0] != 0:
        raise OSError(f'{name}: the archived HEAD {head[:9]} is not in the repository; nothing deleted')
    # Work committed after the archive: the branch moved on, or git still records another HEAD for the worktree.
    code, branch_head = git('rev-parse', '--verify', '--quiet', f'refs/heads/{receipt["branch"]}')
    if code == 0 and branch_head and branch_head != head:
        raise OSError(f'{name}: branch {receipt["branch"]} is at {branch_head[:9]}, not the archived {head[:9]}; '
                      'nothing deleted')
    recorded = [line[5:] for line in listed or [] if line.startswith('HEAD ')]
    if recorded and recorded[0] != head:
        raise OSError(f'{name}: git records its HEAD at {recorded[0][:9]}, not the archived {head[:9]}; nothing deleted')
    patch = patch_path.read_bytes().decode('utf-8', 'replace').strip()
    total, found, missing, _ = landed(None, head, patch)
    if found < total and not allow_unlanded:
        raise OSError(f'{name}: {total - found} of {total} added lines are not in main, e.g. {missing[:2]}; '
                      'merge them first (or the user retires it with --allow-unlanded)')
    in_main = git('merge-base', '--is-ancestor', head, 'main')[0] == 0
    if not in_main and not git('branch', '--contains', head)[1].strip() and not allow_unlanded:
        raise OSError(f'{name}: its HEAD {head[:9]} is on no branch and not in main')
    if git('update-ref', f'refs/archive/{name}', head)[0] != 0:  # keeps the commits reachable
        raise OSError(f'{name}: could not write refs/archive/{name}; nothing deleted')
    folder = worktree if worktree.is_dir() else quarantine if quarantine.is_dir() else None
    if folder is not None:
        if storage.is_link(folder):
            raise OSError(f'{folder}: is a link; nothing deleted')
        links = storage.find_links(folder)
        if links:
            raise OSError(f'{name}: links inside ({links[0]}); nothing deleted')
        code, listing = git('ls-tree', '-r', '-z', head, binary=True)
        if code != 0:
            raise OSError(f'{name}: git ls-tree failed; nothing deleted')
        tracked = {}  # path -> blob id at the archived HEAD
        for record in listing.decode('utf-8', 'replace').split('\0'):
            meta, _, path = record.partition('\t')
            if path and meta.split(' ')[1:2] == ['blob']:
                tracked[path] = meta.split(' ')[2]
        touched = set()  # tracked paths the archived patch changes: their content is in the patch, not at HEAD
        if patch_path.stat().st_size:
            # --numstat only reads the patch (nothing is applied); -z gives unquoted paths, a rename or copy as
            # "added<TAB>deleted<TAB>" followed by its old and new path.
            code, numstat = git('apply', '--numstat', '-z', str(patch_path), binary=True)
            if code != 0:
                raise OSError(f'{name}: git apply --numstat could not read the archived patch; nothing deleted')
            fields = iter(numstat.decode('utf-8', 'replace').split('\0'))
            for field in fields:
                if not field:
                    continue
                path = field.split('\t', 2)[2] if field.count('\t') >= 2 else ''
                if path:
                    touched.add(path)
                else:
                    touched.update((next(fields, ''), next(fields, '')))
            touched.discard('')
        unchanged = []  # tracked, untouched by the patch: must hash to their blob at HEAD
        stamp = receipt_path.stat().st_mtime
        root = storage.long_path(folder)
        stack = [str(root)]
        while stack:
            for entry in os.scandir(stack.pop()):
                if storage.is_link(entry):
                    raise OSError(f'{name}: a link appeared ({entry.path}); nothing deleted')
                relative = Path(entry.path).relative_to(root).as_posix()
                if entry.is_dir(follow_symlinks=False):
                    if not skipped(relative + '/'):
                        stack.append(entry.path)
                    continue
                if relative == '.git' or skipped(relative):
                    continue
                if relative in receipt['sha256']:
                    if sha256(entry.path) != receipt['sha256'][relative]:
                        raise OSError(f'{name}: {relative} changed since the archive; run preserve again')
                elif relative not in tracked:
                    raise OSError(f'{name}: {relative} is neither in the archive nor tracked at {head[:9]}; nothing deleted')
                elif relative not in touched:
                    unchanged.append(relative)
                if entry.stat(follow_symlinks=False).st_mtime > stamp:
                    raise OSError(f'{name}: {relative} is newer than the archive; nothing deleted')
        # Hashed by the main repository's git as the worktree's own checkout (its .gitattributes, the eol settings):
        # an explicit --git-dir, so git never looks for a repository in or above the folder. (--path, which would
        # name the attribute path, cannot be combined with --stdin-paths.)
        for start in range(0, len(unchanged), 500):
            chunk = unchanged[start:start + 500]
            result = subprocess.run(['git', f'--git-dir={ROOT / ".git"}', f'--work-tree={folder}', 'hash-object',
                                     '--stdin-paths'], cwd=folder, input='\n'.join(chunk) + '\n', capture_output=True,
                                    text=True, encoding='utf-8')
            blobs = result.stdout.split()
            if result.returncode != 0 or len(blobs) != len(chunk):
                raise OSError(f'{name}: git hash-object failed ({result.stderr.strip()[:200]}); nothing deleted')
            differing = [relative for relative, blob in zip(chunk, blobs) if blob != tracked[relative]]
            if differing:
                raise OSError(f'{name}: {len(differing)} tracked files differ from {head[:9]} and the archived patch, '
                              f'e.g. {differing[0]}; nothing deleted')
    if listed is not None:
        git('worktree', 'prune')
        if listed_worktree(worktree) is not None:
            raise OSError(f'{name}: git worktree prune left it listed; nothing deleted')
    if folder is not None:
        delete_folder(folder, worktree)
    receipt['status'] = 'retired'
    receipt['retired_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    receipt['resumed'] = True
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')


def drop_branches(names):
    """Bundle the named local branches to storage, verify the bundle and its heads, then delete the branches.

    Refuses main, the current branch and any branch checked out in a worktree. Deletes nothing unless the bundle
    verifies and lists every branch at its current commit; the bundle and a receipt stay in worktree-archive.
    """
    if not names:
        raise OSError('No branches named')
    current = git('rev-parse', '--abbrev-ref', 'HEAD')[1]
    checked_out = {line.removeprefix('branch refs/heads/') for line in git('worktree', 'list', '--porcelain')[1].splitlines()
                   if line.startswith('branch refs/heads/')}
    heads = {}
    for name in names:
        if name in ('main', current) or name in checked_out:
            raise OSError(f'{name}: main, the current branch or checked out in a worktree')
        code, sha = git('rev-parse', '--verify', '--quiet', f'refs/heads/{name}')
        if code != 0 or not sha:
            raise OSError(f'{name}: no such local branch')
        heads[name] = sha
    stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    number = 0
    while (archive_root() / f'branches-{stamp}-{number}.bundle').exists():
        number += 1
    bundle = archive_root() / f'branches-{stamp}-{number}.bundle'  # never overwrites an earlier bundle
    bundle.parent.mkdir(parents=True, exist_ok=True)
    # git writes a bundle through a lock file it then renames; the archive is delete-protected (protect.py), so
    # the bundle is made in the temp folder and copied in.
    with tempfile.TemporaryDirectory() as scratch:
        made = Path(scratch) / bundle.name
        code, out = git('bundle', 'create', str(made), *(f'refs/heads/{name}' for name in heads))
        if code != 0 or not made.is_file():
            raise OSError(f'git bundle create failed: {out}')
        shutil.copy2(made, bundle)
        if sha256(bundle) != sha256(made):
            raise OSError(f'{bundle}: copy verification failed; no branch deleted')
    if git('bundle', 'verify', str(bundle))[0] != 0:
        raise OSError(f'{bundle}: git bundle verify failed; no branch deleted')
    listed = {}
    for line in git('bundle', 'list-heads', str(bundle))[1].splitlines():
        sha, _, ref = line.partition(' ')
        listed[ref.removeprefix('refs/heads/')] = sha
    wrong = [name for name, sha in heads.items() if listed.get(name) != sha]
    if wrong:
        raise OSError(f'{bundle}: heads missing or different for {wrong}; no branch deleted')
    receipt = {'bundle': str(bundle), 'bundle_sha256': sha256(bundle), 'heads': heads, 'status': 'verified'}
    receipt_path = bundle.with_suffix('.json')
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    for name, sha in heads.items():
        if git('rev-parse', f'refs/heads/{name}')[1] != sha:
            raise OSError(f'{name}: moved since the bundle; stopped (branches before it were deleted)')
        code, out = git('branch', '-D', name)
        if code != 0:
            raise OSError(f'{name}: git branch -D failed: {out}')
    receipt['status'] = 'deleted'
    receipt_path.write_text(json.dumps(receipt, indent=1) + '\n', encoding='utf-8')
    return bundle, len(heads)


def sha256_bytes(data):
    import hashlib
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    audit = commands.add_parser('audit', help='read-only report of every worktree')
    audit.add_argument('--json', help='also write the full report here')
    for name in ('preserve', 'retire'):
        command = commands.add_parser(name, allow_abbrev=False)
        command.add_argument('names', nargs='*')
        command.add_argument('--all', action='store_true')
        if name == 'retire':
            command.add_argument('--yes', action='store_true', required=True, help='confirm the delete')
            command.add_argument('--allow-unlanded', action='store_true',
                                 help='the user only: retire although some added lines are not in main')
    drop = commands.add_parser('drop-branches', help='bundle, verify, then delete local branches', allow_abbrev=False)
    drop.add_argument('names', nargs='*')
    drop.add_argument('--retired', action='store_true', help='every branch of a retired archive that still exists')
    args = parser.parse_args()
    if args.command == 'drop-branches':
        names = list(args.names)
        if args.retired:
            for receipt_path in sorted(archive_root().glob('*/receipt.json')):
                receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
                if receipt.get('status') == 'retired' and git('rev-parse', '--verify', '--quiet',
                                                               f'refs/heads/{receipt["branch"]}')[0] == 0:
                    names.append(receipt['branch'])
        bundle, count = drop_branches(sorted(set(names)))
        print(f'{count} branches bundled to {bundle} (verified, heads matched) and deleted')
        return
    # A <name>.retiring folder is a retire stopped after its rename: listed (half-retired) under <name>.
    every = sorted({WORKTREES / item.name.removesuffix('.retiring') for item in WORKTREES.iterdir() if item.is_dir()}) \
        if WORKTREES.is_dir() else []
    if args.command == 'audit':
        report = []
        for worktree in every:
            item = audit_one(worktree)
            report.append(item)
            print(f"{item['name']}  {item['date'][:16]}  {item['gb']:5.2f} GB  links {len(item['links'])}  "
                  f"landed {item['found_in_main']}/{item['added_lines']}  modified {len(item['modified'])}  "
                  f"untracked {len(item['untracked'])}  archive {item['archive']}  "
                  f"{item['change'] if item.get('half_retired') else item['change'][:80]}")
        print(f'total {sum(item["gb"] for item in report):.1f} GB in {len(report)} worktrees')
        if args.json:
            Path(args.json).write_text(json.dumps(report, indent=1) + '\n', encoding='utf-8')
        return
    chosen = every if args.all else [WORKTREES / name for name in args.names]
    for worktree in chosen:
        if args.command == 'preserve':
            if args.all and (not worktree.is_dir() or half_retired(worktree)):
                print(f'{worktree.name}: half-retired, not preserved (resume with retire {worktree.name} --yes)')
                continue
            files, added, gb = preserve(worktree)
            print(f'{worktree.name}: {files} files in the archive ({gb:.2f} GB), {added} added now, and the uncommitted patch')
        else:
            retire(worktree, args.allow_unlanded)
            print(f'{worktree.name}: retired (branch kept; drop-branches --retired bundles and deletes it)')


if __name__ == '__main__':
    main()
