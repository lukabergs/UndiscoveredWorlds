"""Common header for raw run logs, so every module records the same provenance first."""
import datetime
import hashlib
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def header(module, executable, command=None):
    """Return the standard first lines of a run log."""
    executable = Path(executable)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=normal'],
                                         cwd=ROOT, text=True).strip())
    lines = [f'Module: {module}',
             f'Started UTC: {datetime.datetime.now(datetime.timezone.utc).isoformat()}',
             f'Executable: {executable}',
             f'Executable SHA256: {hashlib.sha256(executable.read_bytes()).hexdigest().upper()}',
             f'Source HEAD: {head}; working tree modified: {dirty}']
    if command:
        lines.append('Command: ' + subprocess.list2cmdline([str(part) for part in command]))
    return '\n'.join(lines) + '\n'
