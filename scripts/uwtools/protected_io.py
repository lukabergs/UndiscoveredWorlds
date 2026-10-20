"""Writes into refs/ that work under its deny-delete permission (scripts/uwtools/protect.py).

Inside the protected trees a file can be created or overwritten in place, but nothing can
be deleted, renamed or moved away: a temporary file renamed onto its target fails, and so
do libraries that delete a file before writing it again (GDAL/rasterio 'w', NetCDF writers,
download clients). Such writes go to a fresh folder in the system temp, outside the
protected trees, and are copied onto the target, which is overwritten in place:

    with staged(target) as temporary:
        write(temporary)        # any writer; the target is untouched if this raises

The copy is checked by SHA-256 and the temp folder goes on exit, also after an error.
Copying is not atomic: a failure during the copy itself can leave the target partly
written, where a rename could not.
"""
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path

from .hashing import sha256


@contextmanager
def staging():
    """A fresh folder in the system temp; it and its contents go on exit."""
    with tempfile.TemporaryDirectory(prefix='uw-staged-', ignore_cleanup_errors=True) as folder:
        yield Path(folder)


def install(source, target):
    """Copy source onto target (parent folders created, an existing file overwritten in place)
    and check the copy by SHA-256; returns the target."""
    source, target = Path(source), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha256(target) != sha256(source):
        raise OSError(f'The copy of {source} onto {target} does not match it')
    return target


@contextmanager
def staged(target, companions=False):
    """A temporary path with the target's file name; installed onto the target on normal exit.

    companions also installs the other files the writer left beside it (GDAL's
    <name>.aux.xml) next to the target. On an exception the target is untouched."""
    target = Path(target)
    with staging() as folder:
        temporary = folder / target.name
        yield temporary
        install(temporary, target)
        if companions:
            for path in sorted(folder.iterdir()):
                if path.is_file() and path != temporary:
                    install(path, target.with_name(path.name))


def write_text_into(target, text, encoding='utf-8'):
    """Path.write_text through staged; returns the target."""
    with staged(target) as temporary:
        temporary.write_text(text, encoding=encoding)
    return Path(target)


def write_bytes_into(target, data):
    """Path.write_bytes through staged; returns the target."""
    with staged(target) as temporary:
        temporary.write_bytes(data)
    return Path(target)
