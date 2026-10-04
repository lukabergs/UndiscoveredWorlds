"""Processed reference fields (.uwclim): the only code that knows the file format.

A .uwclim file is a 40-byte header (magic UWCLIM1\\0; little-endian uint32 version 1,
columns, rows, layers; 16-byte NUL-padded ASCII variable name) followed by
layers x rows x columns little-endian float32 values, cell-centred, north row first.
Missing values are NaN or the fill -9999.9 depending on the writer.
"""
import struct
from pathlib import Path
from typing import NamedTuple

import numpy as np

MAGIC = b'UWCLIM1\0'
HEADER_BYTES = 40
# Values below this are the -9999.9 fill some writers use for missing cells.
MISSING_BELOW = -9000


class Header(NamedTuple):
    columns: int
    rows: int
    layers: int
    variable: str

    @property
    def shape(self):
        return self.layers, self.rows, self.columns


def file_size(layers, rows, columns):
    """Bytes of a complete file with this payload shape."""
    return HEADER_BYTES + layers * rows * columns * 4


def read_header(path):
    path = Path(path)
    with path.open('rb') as stream:
        header = stream.read(HEADER_BYTES)
    if len(header) != HEADER_BYTES or header[:8] != MAGIC:
        raise ValueError(f'{path}: not a UWCLIM1 file')
    version, columns, rows, layers = struct.unpack('<IIII', header[8:24])
    if version != 1:
        raise ValueError(f'{path}: unsupported UWCLIM version {version}')
    return Header(columns, rows, layers, header[24:].rstrip(b'\0').decode('ascii'))


def read(path, missing=None, dtype=None, memmap=False):
    """(header, values[layers, rows, columns]); values are float32 as stored unless:

    dtype converts them (float for the validation and report readers); missing replaces
    values below MISSING_BELOW (np.nan where the fill must not enter statistics); memmap
    returns a read-only map of the file instead of reading it (not with missing alone).
    """
    path = Path(path)
    header = read_header(path)
    size = path.stat().st_size
    if size != file_size(*header.shape):
        raise ValueError(f'{path}: {size} bytes, the header describes {file_size(*header.shape)}')
    if memmap:
        values = np.memmap(path, dtype='<f4', mode='r', offset=HEADER_BYTES, shape=header.shape)
    else:
        values = np.fromfile(path, dtype='<f4', offset=HEADER_BYTES).reshape(header.shape)
    if dtype is not None:
        values = values.astype(dtype)
    if missing is not None:
        values[values < MISSING_BELOW] = missing
    return header, values


def write(path, variable, layers, shape=None, mode='wb'):
    """Write a .uwclim file. layers is a [layers, rows, columns] array or a sequence of 2-D arrays;
    an iterator, written layer by layer as it yields, needs shape = (layers, rows, columns).
    mode 'xb' refuses an existing file."""
    if shape is None:
        shape = (len(layers), *np.shape(layers[0]))
    count, rows, columns = shape
    name = variable.encode('ascii')
    if len(name) > 16:
        raise ValueError(f'UWCLIM variable exceeds sixteen bytes: {variable}')
    with Path(path).open(mode) as stream:
        stream.write(MAGIC + struct.pack('<IIII', 1, columns, rows, count) + name.ljust(16, b'\0'))
        for values in layers:
            stream.write(np.asarray(values).astype('<f4', copy=False).tobytes(order='C'))
