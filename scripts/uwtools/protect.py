"""Physical delete protection for the data that cannot be regenerated: an NTFS deny-delete permission.

    uv run python scripts/uwtools/protect.py status
    uv run python scripts/uwtools/protect.py lock            (agents may run it)
    uv run python scripts/uwtools/protect.py harden          (the user, from an elevated terminal: owner = Administrators)
    uv run python scripts/uwtools/protect.py unlock PLACE... (the user only; the guards refuse it for agents)

The user's request (2026-10-20): a safeguard that makes accidental deletes physically impossible. `lock` adds one
explicit, inheritable ACE to each place in PLACES: Everyone (S-1-1-0) denied DELETE and DELETE_CHILD, for the folder,
its subfolders and files ((OI)(CI)(DE,DC)). A deny ACE wins over every allow, for every process of every account,
scripts and Explorer included: files can still be created, read and overwritten in place, but nothing inside can be
deleted, renamed or moved away, and replacing a file by rename (atomic writes) fails. Writers in these places copy
instead of renaming (backup.py, trash.py, worktrees.py).

While the user owns a place, any of the user's processes may change its permissions, so an agent could remove the
ACE; the guards refuse icacls, cacls, takeown, Set-Acl and `protect.py unlock` for agents. `harden` closes that too:
it makes BUILTIN\\Administrators the owner, so that changing the permissions needs an elevated process (a UAC prompt,
which agents never get). Removing the protection, for the user: `unlock PLACE` from an elevated terminal once
hardened (from a normal one before).

Since 2026-11-01 (the user could not delete D:/dev/bak or empty the trash by hand) only refs/ is locked: D:/dev/bak,
worktree-archive and the trash are RELEASED, `lock` no longer takes them and `unlock` lifts what is left on them.
The guards still refuse agents' deletes there (guard_catastrophe.py PROTECTED_TREES).
"""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

PLACES = [
    Path(r'C:\dev\undiscovered-worlds\refs'),
    Path(r'D:\dev\undiscovered-worlds\refs'),  # the storage side of refs/processed/physical and refs/source/physical
]
RELEASED = [  # locked until 2026-11-01; `unlock` lifts them, `status` shows them
    Path(r'D:\dev\bak'),
    Path(r'D:\dev\undiscovered-worlds\worktree-archive'),
    Path(r'D:\dev\undiscovered-worlds\trash'),
]
EVERYONE = '*S-1-1-0'
ADMINISTRATORS = '*S-1-5-32-544'
ACE = f'{EVERYONE}:(OI)(CI)(DE,DC)'


def icacls(*args):
    result = subprocess.run(['icacls', *args], capture_output=True, text=True, encoding='mbcs', errors='replace')
    return result.returncode, result.stdout + result.stderr


def powershell(command):
    # Windows PowerShell cannot load Get-Acl with the module path PowerShell 7 hands down (status read every place
    # as unlocked when run from pwsh): drop it.
    env = {key: value for key, value in os.environ.items() if key.upper() != 'PSMODULEPATH'}
    result = subprocess.run(['powershell', '-NoProfile', '-Command', command], capture_output=True, text=True, env=env)
    return result.stdout.strip()


def sddl(place):
    return powershell(f"(Get-Acl -LiteralPath '{place}').Sddl")


def denies(text, inherited):
    """True when an SDDL string holds a deny ACE for Everyone (WD) with DELETE (SD) and DELETE_CHILD (DT); the
    folder itself must carry it explicitly with OI and CI, a child as inherited (ID). Read in the language-neutral
    SDDL form: icacls prints localised names."""
    for flags, rights in re.findall(r'\(D;([A-Z]*);([A-Za-z0-9]*);;;WD\)', text):
        if 'SD' not in rights:
            continue
        if inherited and 'ID' in flags:
            return True
        if not inherited and 'ID' not in flags and 'OI' in flags and 'CI' in flags and 'DT' in rights:
            return True
    return False


def locked(place):
    """True when the place itself carries the explicit deny-delete ACE (not inherited)."""
    return denies(sddl(place), inherited=False)


def covered(path):
    """True when a file or folder inside a place carries the inherited deny-delete ACE."""
    return denies(sddl(path), inherited=True)


def owner(place):
    return powershell(f"(Get-Acl -LiteralPath '{place}').Owner")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('status', allow_abbrev=False)
    lock = commands.add_parser('lock', allow_abbrev=False)
    lock.add_argument('places', nargs='*', help='some of PLACES (default: all)')
    commands.add_parser('harden', allow_abbrev=False)
    unlock = commands.add_parser('unlock', allow_abbrev=False)
    unlock.add_argument('places', nargs='+')
    args = parser.parse_args()
    if args.command == 'status':
        for place in PLACES + RELEASED:
            state = 'missing' if not place.exists() else 'locked' if locked(place) else 'NOT locked'
            # One file deep inside shows whether the protection reached the existing tree.
            sample = next((item for item in place.rglob('*') if item.is_file()), None) if state == 'locked' else None
            reach = '' if sample is None else ('  (reaches ' if covered(sample) else '  (NOT on ') + str(sample) + ')'
            if place.exists() and not sddl(place):
                state = 'UNKNOWN'  # the permissions could not be read (do not report it as unlocked)
            elif place in RELEASED and state == 'NOT locked':
                state = 'released'
            note = '  (released: run unlock)' if place in RELEASED and state == 'locked' else ''
            print(f'{state:10}  owner {owner(place) if place.exists() else "-":36}  {place}{reach}{note}')
        return
    if args.command == 'lock':
        chosen = [Path(item).resolve() for item in args.places] or PLACES
        unknown = [str(item) for item in chosen if item not in PLACES]
        if unknown:
            sys.exit(f'not protected places: {unknown}')
        for place in chosen:
            if not place.exists():
                print(f'missing     {place}')
                continue
            if locked(place):
                print(f'locked      {place} (already)')
                continue
            code, out = icacls(str(place), '/deny', ACE)  # icacls propagates the inheritable ACE to the tree
            if code != 0 or not locked(place):
                sys.exit(f'lock failed for {place}: {out.strip()}')
            print(f'locked      {place}')
        return
    if args.command == 'harden':
        for place in PLACES:
            if place.exists():
                # /T: every existing file too. Files created later belong to their creator (the user), whose
                # processes could change their permissions: the guards refuse the tools for that.
                code, out = icacls(str(place), '/setowner', ADMINISTRATORS, '/T', '/C', '/Q')
                print(('owner set   ' if code == 0 else 'FAILED      ') + str(place) +
                      ('' if code == 0 else f' ({out.strip()}; run from an elevated terminal)'))
        return
    chosen = [Path(item).resolve() for item in args.places]
    unknown = [str(place) for place in chosen if place not in PLACES + RELEASED]
    if unknown:
        sys.exit(f'not protected or released places: {unknown}')
    failed = False
    for place in chosen:
        code, out = icacls(str(place), '/remove:d', EVERYONE)
        print(out.strip())
        failed = failed or code != 0 or not sddl(place) or locked(place)  # unreadable is not verified
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
