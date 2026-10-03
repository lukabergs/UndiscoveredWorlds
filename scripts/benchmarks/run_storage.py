"""Route new local experiment archives through the configured bulk-data drive."""
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
# PowerShell 7 when installed; the storage script also runs on Windows PowerShell 5.1.
POWERSHELL=shutil.which('pwsh') or 'powershell'

def ensure_archive_directory(path):
    # Keep the logical C: path: resolving the archive's junction would make
    # historical helpers mistake the data drive for the source workspace.
    path=Path(path).absolute()
    if not (ROOT/'runs/storage.json').is_file():return path
    relative=path.relative_to(ROOT)
    if relative.parts[:2]!=('runs','reports') or len(relative.parts)!=3:
        raise ValueError('Expected a direct runs/reports experiment archive')
    subprocess.run([POWERSHELL,'-NoProfile','-ExecutionPolicy','Bypass','-File',str(ROOT/'scripts/benchmarks/configure-run-storage.ps1'),
                    '-RelativePaths',relative.as_posix(),'-EnsureArchive'],check=True,cwd=ROOT)
    return path
