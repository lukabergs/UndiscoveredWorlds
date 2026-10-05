"""Declared input files of a run (profile "run" section, "inputs").

An input is a file the tool reads besides its arguments, handed over through an
environment variable; its path and SHA-256 go into run.json. "derive" names how the file
for the run's grid width is made from "source"; derived files are written once to
out/inputs/ (never under refs/).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DERIVED = ROOT / 'out/inputs'


def drag_contrast(source, width):
    """The vegetation surface-drag contrast (UWROUGH1 text) for a `width` world.

    The 1024-wide ERA5 field is used as it is at its own width and area-averaged once to
    any other width (stress factors are remapped by area, as in the 1024 preparation)."""
    values = source.read_text(encoding='utf-8').split()
    if values[0] != 'UWROUGH1':
        raise SystemExit(f'Unexpected drag input header: {source}')
    columns, rows = int(values[1]), int(values[2])
    if columns == width:
        return source
    path = DERIVED / source.parent.parent.name / f'{source.stem.rsplit("-", 1)[0]}-{width}.txt'
    if path.exists():
        return path
    sys.path.insert(0, str(ROOT / 'scripts/refs'))
    from cell_grid import conservative_remap
    field = np.clip(conservative_remap(np.asarray(values[3:], dtype=np.float64).reshape(rows, columns), width), 0.25, 4.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    with temporary.open('w', encoding='utf-8') as stream:
        stream.write(f'UWROUGH1 {width} {width // 2}\n')
        for row in field:
            stream.write(' '.join(f'{value:.6g}' for value in row) + '\n')
    temporary.replace(path)
    print(f'Derived {path.relative_to(ROOT).as_posix()} from {source.name} by conservative area remap.', flush=True)
    return path


DERIVATIONS = {'drag-contrast': drag_contrast}


def resolve(spec, width):
    """(environment variable, file) of one declared input for a grid of the given width."""
    source = ROOT / spec['source']
    path = DERIVATIONS[spec['derive']](source, width) if spec.get('derive') else source
    return spec['environment'], path
