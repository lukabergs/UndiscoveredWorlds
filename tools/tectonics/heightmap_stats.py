#!/usr/bin/env python3
"""Print min/max/mean/stddev for a heightmap file.

Supported inputs:
- PNG grayscale files with 8-bit or 16-bit samples.
- PNG grayscale+alpha files. Alpha is ignored.
- Uncompressed grayscale TIFF files with 8-bit or 16-bit unsigned samples.

Samples are reported in stored sample units. For the TIFF16 terrain pipeline,
that means metres above the exported zero baseline.
"""

from __future__ import annotations

import argparse
import binascii
import math
import pathlib
import struct
import sys
import zlib


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
ADAM7_PASSES = (
    (0, 0, 8, 8),
    (4, 0, 8, 8),
    (0, 4, 4, 8),
    (2, 0, 4, 4),
    (0, 2, 2, 4),
    (1, 0, 2, 2),
    (0, 1, 1, 2),
)
TIFF_TYPE_SIZES = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 8,
    6: 1,
    7: 1,
    8: 2,
    9: 4,
    10: 8,
    11: 4,
    12: 8,
}
TIFF_TAG_IMAGE_WIDTH = 256
TIFF_TAG_IMAGE_LENGTH = 257
TIFF_TAG_BITS_PER_SAMPLE = 258
TIFF_TAG_COMPRESSION = 259
TIFF_TAG_PHOTOMETRIC_INTERPRETATION = 262
TIFF_TAG_STRIP_OFFSETS = 273
TIFF_TAG_SAMPLES_PER_PIXEL = 277
TIFF_TAG_ROWS_PER_STRIP = 278
TIFF_TAG_STRIP_BYTE_COUNTS = 279
TIFF_TAG_PLANAR_CONFIGURATION = 284
TIFF_TAG_SAMPLE_FORMAT = 339
TIFF_COMPRESSION_NONE = 1
TIFF_PHOTOMETRIC_MINISBLACK = 1
TIFF_PLANAR_CONFIGURATION_CONTIG = 1
TIFF_SAMPLE_FORMAT_UINT = 1


class StatsAccumulator:
    def __init__(self) -> None:
        self.count = 0
        self.minimum = math.inf
        self.maximum = -math.inf
        self.mean = 0.0
        self.m2 = 0.0

    def add(self, value: float) -> None:
        if not math.isfinite(value):
            raise ValueError("Encountered a non-finite sample value.")

        self.count += 1
        if value < self.minimum:
            self.minimum = value
        if value > self.maximum:
            self.maximum = value

        delta = value - self.mean
        self.mean += delta / self.count
        delta2 = value - self.mean
        self.m2 += delta * delta2

    @property
    def stddev(self) -> float:
        if self.count == 0:
            return float("nan")
        variance = self.m2 / self.count
        if variance < 0.0:
            variance = 0.0
        return math.sqrt(variance)


def format_number(value: float) -> str:
    return format(value, ".17g")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filepath", help="Path to a PNG or TIFF heightmap.")
    return parser.parse_args()


def pass_length(length: int, start: int, step: int) -> int:
    if start >= length:
        return 0
    return 1 + (length - 1 - start) // step


def paeth_predictor(a: int, b: int, c: int) -> int:
    prediction = a + b - c
    distance_a = abs(prediction - a)
    distance_b = abs(prediction - b)
    distance_c = abs(prediction - c)
    if distance_a <= distance_b and distance_a <= distance_c:
        return a
    if distance_b <= distance_c:
        return b
    return c


def unfilter_row(filter_type: int, row: bytes, previous_row: bytes, bytes_per_pixel: int) -> bytearray:
    result = bytearray(row)
    if filter_type == 0:
        return result

    for index in range(len(result)):
        left = result[index - bytes_per_pixel] if index >= bytes_per_pixel else 0
        up = previous_row[index] if previous_row else 0
        up_left = previous_row[index - bytes_per_pixel] if previous_row and index >= bytes_per_pixel else 0

        if filter_type == 1:
            predictor = left
        elif filter_type == 2:
            predictor = up
        elif filter_type == 3:
            predictor = (left + up) // 2
        elif filter_type == 4:
            predictor = paeth_predictor(left, up, up_left)
        else:
            raise ValueError(f"Unsupported PNG filter type: {filter_type}")

        result[index] = (result[index] + predictor) & 0xFF

    return result


def read_png_chunks(path: pathlib.Path) -> tuple[int, int, int, int, int, bytes]:
    with path.open("rb") as handle:
        signature = handle.read(len(PNG_SIGNATURE))
        if signature != PNG_SIGNATURE:
            raise ValueError("Not a PNG file.")

        width = None
        height = None
        bit_depth = None
        color_type = None
        interlace_method = None
        idat_parts: list[bytes] = []

        while True:
            length_bytes = handle.read(4)
            if len(length_bytes) != 4:
                raise ValueError("Unexpected end of file while reading PNG chunks.")

            length = struct.unpack(">I", length_bytes)[0]
            chunk_type = handle.read(4)
            if len(chunk_type) != 4:
                raise ValueError("Unexpected end of file while reading PNG chunk type.")

            chunk_data = handle.read(length)
            if len(chunk_data) != length:
                raise ValueError("Unexpected end of file while reading PNG chunk data.")

            crc_bytes = handle.read(4)
            if len(crc_bytes) != 4:
                raise ValueError("Unexpected end of file while reading PNG chunk CRC.")

            expected_crc = struct.unpack(">I", crc_bytes)[0]
            actual_crc = binascii.crc32(chunk_type)
            actual_crc = binascii.crc32(chunk_data, actual_crc) & 0xFFFFFFFF
            if actual_crc != expected_crc:
                raise ValueError(f"CRC mismatch in PNG chunk {chunk_type.decode('ascii', errors='replace')}.")

            if chunk_type == b"IHDR":
                if length != 13:
                    raise ValueError("Invalid IHDR chunk length.")
                width, height, bit_depth, color_type, compression, filtering, interlace_method = struct.unpack(
                    ">IIBBBBB", chunk_data
                )
                if compression != 0 or filtering != 0:
                    raise ValueError("Unsupported PNG compression or filter method.")
            elif chunk_type == b"IDAT":
                idat_parts.append(chunk_data)
            elif chunk_type == b"IEND":
                break

        if width is None or height is None or bit_depth is None or color_type is None or interlace_method is None:
            raise ValueError("PNG file is missing an IHDR chunk.")
        if not idat_parts:
            raise ValueError("PNG file has no IDAT data.")

        return width, height, bit_depth, color_type, interlace_method, b"".join(idat_parts)


def describe_png(bit_depth: int, color_type: int) -> str:
    if color_type == 0:
        return f"png grayscale {bit_depth}-bit"
    if color_type == 4:
        return f"png grayscale+alpha {bit_depth * 2}-bit-per-pixel"
    raise ValueError("Unsupported PNG color type.")


def analyze_png(path: pathlib.Path) -> tuple[str, int, int, StatsAccumulator]:
    width, height, bit_depth, color_type, interlace_method, idat_data = read_png_chunks(path)

    if color_type not in (0, 4):
        raise ValueError("Unsupported PNG color type. Expected grayscale or grayscale+alpha.")
    if bit_depth not in (8, 16):
        raise ValueError("Unsupported PNG bit depth. Expected 8-bit or 16-bit grayscale samples.")
    if interlace_method not in (0, 1):
        raise ValueError("Unsupported PNG interlace method.")

    channels = 1 if color_type == 0 else 2
    bytes_per_sample = 1 if bit_depth == 8 else 2
    bytes_per_pixel = channels * bytes_per_sample

    try:
        raw_data = zlib.decompress(idat_data)
    except zlib.error as exc:
        raise ValueError("Failed to decompress PNG image data.") from exc

    offset = 0
    stats = StatsAccumulator()
    passes = ((0, 0, 1, 1),) if interlace_method == 0 else ADAM7_PASSES

    for start_x, start_y, step_x, step_y in passes:
        pass_width = pass_length(width, start_x, step_x)
        pass_height = pass_length(height, start_y, step_y)
        if pass_width == 0 or pass_height == 0:
            continue

        row_bytes = pass_width * bytes_per_pixel
        previous_row = b""
        for _ in range(pass_height):
            if offset >= len(raw_data):
                raise ValueError("PNG image data ended before all rows were decoded.")

            filter_type = raw_data[offset]
            offset += 1
            row = raw_data[offset : offset + row_bytes]
            if len(row) != row_bytes:
                raise ValueError("PNG image data ended in the middle of a scanline.")
            offset += row_bytes

            decoded_row = unfilter_row(filter_type, row, previous_row, bytes_per_pixel)
            previous_row = bytes(decoded_row)

            if bit_depth == 8:
                for index in range(0, row_bytes, channels):
                    stats.add(float(decoded_row[index]))
            else:
                step = channels * 2
                for index in range(0, row_bytes, step):
                    sample = struct.unpack_from(">H", decoded_row, index)[0]
                    stats.add(float(sample))

    if offset != len(raw_data):
        raise ValueError("PNG image data contains trailing bytes after the final scanline.")

    return describe_png(bit_depth, color_type), width, height, stats


def read_exact(handle, byte_count: int, description: str) -> bytes:
    data = handle.read(byte_count)
    if len(data) != byte_count:
        raise ValueError(f"Unexpected end of file while reading {description}.")
    return data


def read_tiff_header(handle) -> tuple[str, str, int]:
    byte_order = read_exact(handle, 2, "TIFF byte order")
    if byte_order == b"II":
        endian = "<"
        endian_name = "little"
    elif byte_order == b"MM":
        endian = ">"
        endian_name = "big"
    else:
        raise ValueError("Not a TIFF file.")

    version = struct.unpack(endian + "H", read_exact(handle, 2, "TIFF version"))[0]
    if version == 43:
        raise ValueError("BigTIFF is not supported.")
    if version != 42:
        raise ValueError(f"Unsupported TIFF version: {version}")

    ifd_offset = struct.unpack(endian + "I", read_exact(handle, 4, "TIFF first IFD offset"))[0]
    if ifd_offset == 0:
        raise ValueError("TIFF file has no image file directory.")
    return endian, endian_name, ifd_offset


def read_tiff_ifd(handle, endian: str, ifd_offset: int) -> dict[int, tuple[int, int, bytes]]:
    handle.seek(ifd_offset)
    entry_count = struct.unpack(endian + "H", read_exact(handle, 2, "TIFF IFD entry count"))[0]
    entries: dict[int, tuple[int, int, bytes]] = {}
    for _ in range(entry_count):
        entry = read_exact(handle, 12, "TIFF IFD entry")
        tag, type_id, count = struct.unpack(endian + "HHI", entry[:8])
        entries[tag] = (type_id, count, entry[8:12])
    read_exact(handle, 4, "TIFF next IFD offset")
    return entries


def read_tiff_values(handle, endian: str, entry: tuple[int, int, bytes]) -> tuple[int | float, ...]:
    type_id, count, value_or_offset = entry
    type_size = TIFF_TYPE_SIZES.get(type_id)
    if type_size is None:
        raise ValueError(f"Unsupported TIFF field type: {type_id}")

    total_size = type_size * count
    if total_size <= 4:
        raw = value_or_offset[:total_size]
    else:
        data_offset = struct.unpack(endian + "I", value_or_offset)[0]
        current_offset = handle.tell()
        handle.seek(data_offset)
        raw = read_exact(handle, total_size, "TIFF field data")
        handle.seek(current_offset)

    if type_id in (1, 7):
        return tuple(raw[index] for index in range(count))
    if type_id == 3:
        return struct.unpack(endian + f"{count}H", raw)
    if type_id == 4:
        return struct.unpack(endian + f"{count}I", raw)
    if type_id == 6:
        return struct.unpack(endian + f"{count}b", raw)
    if type_id == 8:
        return struct.unpack(endian + f"{count}h", raw)
    if type_id == 9:
        return struct.unpack(endian + f"{count}i", raw)
    if type_id == 11:
        return struct.unpack(endian + f"{count}f", raw)
    if type_id == 12:
        return struct.unpack(endian + f"{count}d", raw)
    raise ValueError(f"Unsupported TIFF field type: {type_id}")


def get_tiff_scalar(handle, endian: str, entries: dict[int, tuple[int, int, bytes]], tag: int, default=None):
    entry = entries.get(tag)
    if entry is None:
        if default is not None:
            return default
        raise ValueError(f"Missing required TIFF tag: {tag}")

    values = read_tiff_values(handle, endian, entry)
    if len(values) != 1:
        raise ValueError(f"Expected one value for TIFF tag {tag}, found {len(values)}.")
    return values[0]


def get_tiff_values(handle, endian: str, entries: dict[int, tuple[int, int, bytes]], tag: int) -> tuple[int | float, ...]:
    entry = entries.get(tag)
    if entry is None:
        raise ValueError(f"Missing required TIFF tag: {tag}")
    return read_tiff_values(handle, endian, entry)


def describe_tiff(bits_per_sample: int, sample_format: int) -> str:
    if sample_format == TIFF_SAMPLE_FORMAT_UINT:
        return f"tiff grayscale {bits_per_sample}-bit"
    raise ValueError("Unsupported TIFF sample format.")


def analyze_tiff(path: pathlib.Path) -> tuple[str, int, int, StatsAccumulator]:
    with path.open("rb") as handle:
        endian, endian_name, ifd_offset = read_tiff_header(handle)
        entries = read_tiff_ifd(handle, endian, ifd_offset)

        width = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_IMAGE_WIDTH))
        height = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_IMAGE_LENGTH))
        if width <= 0 or height <= 0:
            raise ValueError("TIFF dimensions must be positive.")

        samples_per_pixel = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_SAMPLES_PER_PIXEL, 1))
        if samples_per_pixel != 1:
            raise ValueError("Unsupported TIFF layout. Expected a single grayscale sample per pixel.")

        bits_values = [int(value) for value in get_tiff_values(handle, endian, entries, TIFF_TAG_BITS_PER_SAMPLE)]
        if len(bits_values) != 1:
            raise ValueError("Unsupported TIFF BitsPerSample layout.")
        bits_per_sample = bits_values[0]

        compression = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_COMPRESSION, TIFF_COMPRESSION_NONE))
        if compression != TIFF_COMPRESSION_NONE:
            raise ValueError("Unsupported TIFF compression. Expected uncompressed strips.")

        photometric = int(
            get_tiff_scalar(handle, endian, entries, TIFF_TAG_PHOTOMETRIC_INTERPRETATION)
        )
        if photometric != TIFF_PHOTOMETRIC_MINISBLACK:
            raise ValueError("Unsupported TIFF photometric interpretation. Expected grayscale minisblack.")

        planar_configuration = int(
            get_tiff_scalar(
                handle,
                endian,
                entries,
                TIFF_TAG_PLANAR_CONFIGURATION,
                TIFF_PLANAR_CONFIGURATION_CONTIG,
            )
        )
        if planar_configuration != TIFF_PLANAR_CONFIGURATION_CONTIG:
            raise ValueError("Unsupported TIFF planar configuration.")

        sample_format = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_SAMPLE_FORMAT, TIFF_SAMPLE_FORMAT_UINT))
        if bits_per_sample not in (8, 16):
            raise ValueError("Unsupported TIFF sample size. Expected 8 or 16 bits.")
        if bits_per_sample in (8, 16) and sample_format != TIFF_SAMPLE_FORMAT_UINT:
            raise ValueError("Unsupported TIFF sample format. Expected unsigned integer grayscale samples.")

        bytes_per_sample, remainder = divmod(bits_per_sample, 8)
        if remainder != 0:
            raise ValueError("Unsupported TIFF sample size. Expected whole-byte samples.")

        strip_offsets = [int(value) for value in get_tiff_values(handle, endian, entries, TIFF_TAG_STRIP_OFFSETS)]
        strip_byte_counts = [
            int(value) for value in get_tiff_values(handle, endian, entries, TIFF_TAG_STRIP_BYTE_COUNTS)
        ]
        if len(strip_offsets) != len(strip_byte_counts):
            raise ValueError("TIFF StripOffsets and StripByteCounts length mismatch.")

        rows_per_strip = int(get_tiff_scalar(handle, endian, entries, TIFF_TAG_ROWS_PER_STRIP, height))
        if rows_per_strip <= 0:
            raise ValueError("Invalid TIFF RowsPerStrip value.")

        expected_strip_count = (height + rows_per_strip - 1) // rows_per_strip
        if len(strip_offsets) != expected_strip_count:
            raise ValueError("Unexpected TIFF strip count for the declared image size.")

        stats = StatsAccumulator()
        row_start = 0
        for strip_index, (strip_offset, strip_byte_count) in enumerate(zip(strip_offsets, strip_byte_counts)):
            rows_in_strip = min(rows_per_strip, height - row_start)
            sample_count = rows_in_strip * width
            expected_bytes = sample_count * bytes_per_sample
            if strip_byte_count < expected_bytes:
                raise ValueError(f"TIFF strip {strip_index} is shorter than expected.")

            handle.seek(strip_offset)
            strip_data = read_exact(handle, strip_byte_count, f"TIFF strip {strip_index}")
            sample_data = strip_data[:expected_bytes]

            if bits_per_sample == 8:
                for value in sample_data:
                    stats.add(float(value))
            elif bits_per_sample == 16:
                for (value,) in struct.iter_unpack(endian + "H", sample_data):
                    stats.add(float(value))
            row_start += rows_in_strip

        if row_start != height:
            raise ValueError("TIFF strip decoding did not cover the full image.")

    return f"{describe_tiff(bits_per_sample, sample_format)} ({endian_name}-endian)", width, height, stats


def analyze_heightmap(path: pathlib.Path) -> tuple[str, int | None, int | None, StatsAccumulator]:
    suffix = path.suffix.lower()
    if suffix == ".png":
        return analyze_png(path)
    if suffix in (".tif", ".tiff"):
        return analyze_tiff(path)
    raise ValueError(f"Unsupported file type: {path.suffix}")


def main() -> int:
    args = parse_args()
    path = pathlib.Path(args.filepath)

    if not path.is_file():
        print(f"File not found: {path}", file=sys.stderr)
        return 1

    try:
        format_name, width, height, stats = analyze_heightmap(path)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"file: {path}")
    print(f"format: {format_name}")
    if width is not None and height is not None:
        print(f"width: {width}")
        print(f"height: {height}")
    print(f"samples: {stats.count}")
    print(f"min: {format_number(stats.minimum)}")
    print(f"max: {format_number(stats.maximum)}")
    print(f"mean: {format_number(stats.mean)}")
    print(f"stddev: {format_number(stats.stddev)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
