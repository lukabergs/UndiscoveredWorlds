"""Physical delete protection for the data that cannot be regenerated: an NTFS deny-delete permission.

    uv run python scripts/uwtools/protect.py status
    uv run python scripts/uwtools/protect.py lock            (agents may run it)
    uv run python scripts/uwtools/protect.py harden          (the user, from an elevated terminal: owner = Administrators)
    uv run python scripts/uwtools/protect.py unlock PLACE    (the user only; the guards refuse it for agents)

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
"""
import argparse
import subprocess
import sys
from pathlib import Path

PLACES = [
    Path(r'C:\dev\undiscovered-worlds\refs'),
    Path(r'D:\dev\undiscovered-worlds\refs'),  # the storage side of refs/processed/physical and refs/source/physical
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


def locked(place):
    """True when the place itself carries the explicit deny-delete ACE (not inherited)."""
    code, out = icacls(str(place))
    lines = [line for line in out.splitlines() if 'DENY' in line.upper() and '(I)' not in line]
    return code == 0 and any('(DE,DC)' in line and '(OI)(CI)' in line for line in lines)


def owner(place):
    result = subprocess.run(['powershell', '-NoProfile', '-Command', f"(Get-Acl -LiteralPath '{place}').Owner"],
                            capture_output=True, text=True)
    return result.stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('status', allow_abbrev=False)
    lock = commands.add_parser('lock', allow_abbrev=False)
    lock.add_argument('places', nargs='*', help='some of PLACES (default: all)')
    commands.add_parser('harden', allow_abbrev=False)
    unlock = commands.add_parser('unlock', allow_abbrev=False)
    unlock.add_argument('place')
    args = parser.parse_args()
    if args.command == 'status':
        for place in PLACES:
            state = 'missing' if not place.exists() else 'locked' if locked(place) else 'NOT locked'
            print(f'{state:10}  owner {owner(place) if place.exists() else "-":36}  {place}')
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
    place = Path(args.place).resolve()
    if place not in PLACES:
        sys.exit(f'{place} is not one of the protected places')
    code, out = icacls(str(place), '/remove:d', EVERYONE)
    print(out.strip())
    sys.exit(0 if code == 0 and not locked(place) else 1)


if __name__ == '__main__':
    main()
