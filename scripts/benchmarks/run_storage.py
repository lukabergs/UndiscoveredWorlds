"""Route new local experiment archives through the configured bulk-data drive."""
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]

def ensure_archive_directory(path):
    # Keep the logical C: path: resolving the archive's junction would make
    # historical helpers mistake the data drive for the source workspace.
    path=Path(path).absolute()
    if not (ROOT/'runs/storage.json').is_file():return path
    relative=path.relative_to(ROOT)
    if relative.parts[:2]!=('runs','reports') or len(relative.parts)!=3:
        raise ValueError('Expected a direct runs/reports experiment archive')
    subprocess.run(['pwsh','-NoProfile','-File',str(ROOT/'scripts/benchmarks/configure-run-storage.ps1'),
                    '-RelativePaths',relative.as_posix(),'-EnsureArchive'],check=True,cwd=ROOT)
    return path
