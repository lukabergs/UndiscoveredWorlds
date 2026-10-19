"""Processed observed reference fields: the only code that knows their layout and format.

    refs/processed/<dataset>/<width>/manifest.json
    refs/processed/<dataset>/<width>/fields/<name>.npy

A field is a version 1.0 .npy of little-endian float32, C order, shape (layers, rows,
columns): 12 monthly layers or 1 annual layer on a cell-centred W x W/2 grid, north row
first. Missing cells are NaN or the fill -9999.9 (below MISSING_BELOW) depending on the
writer. Callers name a field logically as `<dataset>_<name>` (era5_u10m_monthly):
path('era5_u10m') adds the monthly suffix, path('imerg_prec_annual', annual=True) does not.

manifest.json (schema uw-reference/v1) follows the run-bundle manifest where it applies:
dataset, grid, missing_below and one entry per field (id <dataset>.<name>, variable,
layout monthly|annual, dtype, shape, units, file, sha256).
"""
import hashlib
import json
import os
from pathlib import Path
from typing import NamedTuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
from .paths import DATA_ROOT  # noqa: E402
PROCESSED = DATA_ROOT / 'refs/processed'
SCHEMA = 'uw-reference/v1'
DTYPE = '<f4'
# Every processed field is prepared on this grid width (rows = width / 2).
WIDTH = 2048
# Values below this are the -9999.9 fill some writers use for missing cells.
MISSING_BELOW = -9000
# Dataset prefixes of logical names, longest first so era5land is not read as era5.
DATASETS = ('era5land', 'worldclim', 'glorys', 'ceres', 'imerg', 'oisst', 'era5')
GRID = {'columns': None, 'rows': None, 'registration': 'cell-centred', 'layout': 'row-major, north row first'}


class Header(NamedTuple):
    columns: int
    rows: int
    layers: int
    variable: str

    @property
    def shape(self):
        return self.layers, self.rows, self.columns


def split(stem):
    """('era5', 'u10m_monthly') for the logical stem 'era5_u10m_monthly'."""
    for dataset in DATASETS:
        if stem.startswith(dataset + '_'):
            return dataset, stem[len(dataset) + 1:]
    raise ValueError(f'Unknown reference dataset in {stem!r}; known: {", ".join(DATASETS)}')


def directory(dataset, width=WIDTH, root=PROCESSED):
    """The dataset directory that holds manifest.json, fields/ and the preparation receipts."""
    return Path(root) / dataset / str(width)


def field_path(dataset, name, width=WIDTH, root=PROCESSED):
    return directory(dataset, width, root) / 'fields' / f'{name}.npy'


def path(name, annual=False, width=None, root=PROCESSED):
    """File of a logical reference: 'era5_u10m' is era5_u10m_monthly, 'imerg_prec_annual' with annual=True.

    Without width, the one width the field exists at, else WIDTH."""
    dataset, field = split(name if annual else f'{name}_monthly')
    if width is None:
        found = sorted(Path(root, dataset).glob(f'*/fields/{field}.npy'))
        if len(found) == 1:
            return found[0]
        width = WIDTH
    return field_path(dataset, field, width, root)


def manifest_path(field_file):
    return Path(field_file).parent.parent / 'manifest.json'


def read_manifest(folder):
    path = Path(folder) / 'manifest.json'
    manifest = json.loads(path.read_text(encoding='utf-8'))
    if manifest.get('schema') != SCHEMA:
        raise ValueError(f'{path}: not a {SCHEMA} manifest')
    return manifest


def manifests(root=PROCESSED):
    """(directory, manifest) for every processed reference dataset under root."""
    for path in sorted(Path(root).glob('*/*/manifest.json')):
        manifest = json.loads(path.read_text(encoding='utf-8'))
        if manifest.get('schema') == SCHEMA:
            yield path.parent, manifest


def entry(field_file):
    """The manifest entry of a field file, or None when no manifest lists it."""
    field_file = Path(field_file)
    manifest = manifest_path(field_file)
    if not manifest.is_file():
        return None
    relative = f'fields/{field_file.name}'
    return next((item for item in read_manifest(manifest.parent)['fields'] if item['file'] == relative), None)


def _open_npy(path):
    """(shape, payload offset) of a version 1.0 little-endian float32 C-order .npy."""
    with Path(path).open('rb') as stream:
        version = np.lib.format.read_magic(stream)
        if version != (1, 0):
            raise ValueError(f'{path}: .npy version {version}, expected 1.0')
        shape, fortran, dtype = np.lib.format.read_array_header_1_0(stream)
        offset = stream.tell()
    if dtype != np.dtype(DTYPE) or fortran or len(shape) != 3:
        raise ValueError(f'{path}: expected a C-order {DTYPE} array of (layers, rows, columns), found {dtype} {shape}')
    return shape, offset


def read_header(path):
    """Header of a field; the variable comes from the manifest (the file stem when none lists it)."""
    shape, _ = _open_npy(path)
    item = entry(path)
    return Header(shape[2], shape[1], shape[0], item['variable'] if item else Path(path).stem)


def payload_offset(path):
    """Bytes before the values: the .npy magic, version and header."""
    return _open_npy(path)[1]


def check(path):
    """Header of a field after checking that the file holds exactly the declared values."""
    path = Path(path)
    shape, offset = _open_npy(path)
    expected = offset + int(np.prod(shape)) * 4
    if path.stat().st_size != expected:
        raise ValueError(f'{path}: {path.stat().st_size} bytes, the header describes {expected}')
    return read_header(path)


def read(path, missing=None, dtype=None, memmap=False):
    """(header, values[layers, rows, columns]); values are float32 as stored unless:

    dtype converts them (float for the validation and report readers); missing replaces
    values below MISSING_BELOW (np.nan where the fill must not enter statistics); memmap
    returns a read-only map of the file instead of reading it (not with missing alone).
    """
    path = Path(path)
    header = check(path)
    shape, offset = header.shape, payload_offset(path)
    if memmap:
        values = np.memmap(path, dtype=DTYPE, mode='r', offset=offset, shape=shape)
    else:
        values = np.fromfile(path, dtype=DTYPE, offset=offset).reshape(shape)
    if dtype is not None:
        values = values.astype(dtype)
    if missing is not None:
        values[values < MISSING_BELOW] = missing
    return header, values


def _header_bytes(shape):
    import io
    stream = io.BytesIO()
    np.lib.format.write_array_header_1_0(stream, {'descr': DTYPE, 'fortran_order': False, 'shape': tuple(shape)})
    return stream.getvalue()


def write_npy(path, layers, shape=None, mode='wb'):
    """Write a field file; returns its sha256. layers is a [layers, rows, columns] array or a
    sequence of 2-D arrays; an iterator, written layer by layer as it yields, needs shape =
    (layers, rows, columns). mode 'xb' refuses an existing file."""
    if shape is None:
        shape = (len(layers), *np.shape(layers[0]))
    shape = tuple(int(extent) for extent in shape)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    written = 0
    with path.open(mode) as stream:
        header = _header_bytes(shape)
        stream.write(header)
        digest.update(header)
        for values in layers:
            data = np.ascontiguousarray(np.asarray(values).astype(DTYPE, copy=False)).tobytes(order='C')
            stream.write(data)
            digest.update(data)
            written += len(data)
    if written != shape[0] * shape[1] * shape[2] * 4:
        raise ValueError(f'{path}: wrote {written} payload bytes for shape {shape}')
    return digest.hexdigest()


def update_manifest(dataset, width, entry_value, root=PROCESSED):
    """Add or replace one field entry in <root>/<dataset>/<width>/manifest.json."""
    folder = directory(dataset, width, root)
    path = folder / 'manifest.json'
    if path.is_file():
        manifest = read_manifest(folder)
    else:
        manifest = {'schema': SCHEMA, 'dataset': dataset, 'grid': dict(GRID, columns=width, rows=width // 2),
                    'missing_below': MISSING_BELOW, 'fields': []}
    fields = [item for item in manifest['fields'] if item['id'] != entry_value['id']] + [entry_value]
    manifest['fields'] = sorted(fields, key=lambda item: item['id'])
    lines = ['{']
    for key in ('schema', 'dataset', 'grid', 'missing_below'):
        lines.append(f'  {json.dumps(key)}: {json.dumps(manifest[key])},')
    lines.append('  "fields": [')
    lines.append(',\n'.join(f'    {json.dumps(item)}' for item in manifest['fields']))
    lines += ['  ]', '}', '']
    temporary = path.with_suffix('.json.part')
    temporary.write_text('\n'.join(lines), encoding='utf-8')
    os.replace(temporary, path)
    return path


def write(dataset, name, variable, layers, shape=None, units='', root=PROCESSED, mode='wb'):
    """Write <root>/<dataset>/<width>/fields/<name>.npy and its manifest entry; returns the path.

    name is the logical stem without the dataset prefix (u10m_monthly, prec_annual); a
    one-layer field is annual, twelve layers monthly. See write_npy for layers and shape."""
    if shape is None:
        shape = (len(layers), *np.shape(layers[0]))
    count, rows, columns = (int(extent) for extent in shape)
    if count not in (1, 12) or columns != 2 * rows:
        raise ValueError(f'{dataset}.{name}: shape {shape} is not 1 or 12 layers of W x W/2')
    if not name.startswith(variable):
        raise ValueError(f'{dataset}.{name}: the field name must begin with its variable {variable!r}')
    target = field_path(dataset, name, columns, root)
    temporary = target.with_name(target.name + '.part')
    digest = write_npy(temporary, layers, (count, rows, columns), mode='wb')
    if mode == 'xb' and target.exists():
        temporary.unlink()
        raise FileExistsError(target)
    os.replace(temporary, target)
    update_manifest(dataset, columns, {
        'id': f'{dataset}.{name}', 'variable': variable, 'layout': 'annual' if count == 1 else 'monthly',
        'dtype': DTYPE, 'shape': [count, rows, columns], 'units': units, 'file': f'fields/{name}.npy',
        'sha256': digest}, root)
    return target
