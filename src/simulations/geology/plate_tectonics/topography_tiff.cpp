#include "topography_tiff.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <vector>

namespace platec::io {
namespace {

constexpr uint16_t kTiffVersion = 42;
constexpr uint16_t kTypeShort = 3;
constexpr uint16_t kTypeLong = 4;

constexpr uint16_t kTagImageWidth = 256;
constexpr uint16_t kTagImageLength = 257;
constexpr uint16_t kTagBitsPerSample = 258;
constexpr uint16_t kTagCompression = 259;
constexpr uint16_t kTagPhotometricInterpretation = 262;
constexpr uint16_t kTagStripOffsets = 273;
constexpr uint16_t kTagSamplesPerPixel = 277;
constexpr uint16_t kTagRowsPerStrip = 278;
constexpr uint16_t kTagStripByteCounts = 279;
constexpr uint16_t kTagPlanarConfiguration = 284;
constexpr uint16_t kTagSampleFormat = 339;

constexpr uint16_t kCompressionNone = 1;
constexpr uint16_t kPhotometricMinIsBlack = 1;
constexpr uint16_t kPlanarConfigurationContig = 1;
constexpr uint16_t kSampleFormatUInt = 1;

struct TiffEntry {
    uint16_t tag = 0;
    uint16_t type = 0;
    uint32_t count = 0;
    std::array<unsigned char, 4> value_or_offset{};
};

void writeU16Le(std::ostream& output, uint16_t value)
{
    const unsigned char bytes[2] = {
        static_cast<unsigned char>(value & 0xFFu),
        static_cast<unsigned char>((value >> 8) & 0xFFu),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writeU32Le(std::ostream& output, uint32_t value)
{
    const unsigned char bytes[4] = {
        static_cast<unsigned char>(value & 0xFFu),
        static_cast<unsigned char>((value >> 8) & 0xFFu),
        static_cast<unsigned char>((value >> 16) & 0xFFu),
        static_cast<unsigned char>((value >> 24) & 0xFFu),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writeEntryLong(std::ostream& output, uint16_t tag, uint32_t value)
{
    writeU16Le(output, tag);
    writeU16Le(output, kTypeLong);
    writeU32Le(output, 1U);
    writeU32Le(output, value);
}

void writeEntryShort(std::ostream& output, uint16_t tag, uint16_t value)
{
    writeU16Le(output, tag);
    writeU16Le(output, kTypeShort);
    writeU32Le(output, 1U);
    writeU16Le(output, value);
    writeU16Le(output, 0U);
}

bool readExact(std::ifstream& input, unsigned char* buffer, size_t size)
{
    input.read(reinterpret_cast<char*>(buffer), static_cast<std::streamsize>(size));
    return static_cast<size_t>(input.gcount()) == size;
}

uint16_t readU16(const unsigned char* bytes, bool little_endian)
{
    if (little_endian) {
        return static_cast<uint16_t>(bytes[0]) |
               static_cast<uint16_t>(static_cast<uint16_t>(bytes[1]) << 8);
    }
    return static_cast<uint16_t>(bytes[1]) |
           static_cast<uint16_t>(static_cast<uint16_t>(bytes[0]) << 8);
}

uint32_t readU32(const unsigned char* bytes, bool little_endian)
{
    if (little_endian) {
        return static_cast<uint32_t>(bytes[0]) |
               (static_cast<uint32_t>(bytes[1]) << 8) |
               (static_cast<uint32_t>(bytes[2]) << 16) |
               (static_cast<uint32_t>(bytes[3]) << 24);
    }
    return static_cast<uint32_t>(bytes[3]) |
           (static_cast<uint32_t>(bytes[2]) << 8) |
           (static_cast<uint32_t>(bytes[1]) << 16) |
           (static_cast<uint32_t>(bytes[0]) << 24);
}

size_t typeSize(uint16_t type)
{
    switch (type) {
    case kTypeShort:
        return sizeof(uint16_t);
    case kTypeLong:
        return sizeof(uint32_t);
    default:
        return 0;
    }
}

bool readEntryPayload(std::ifstream& input, const TiffEntry& entry, bool little_endian,
                      std::vector<unsigned char>& payload)
{
    const size_t element_size = typeSize(entry.type);
    if (element_size == 0U) {
        return false;
    }

    const uint64_t total_size =
        static_cast<uint64_t>(entry.count) * static_cast<uint64_t>(element_size);
    if (total_size > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        return false;
    }

    payload.resize(static_cast<size_t>(total_size));
    if (total_size <= entry.value_or_offset.size()) {
        std::copy_n(entry.value_or_offset.begin(), payload.size(), payload.begin());
        return true;
    }

    const uint32_t offset = readU32(entry.value_or_offset.data(), little_endian);
    input.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    return readExact(input, payload.data(), payload.size());
}

bool readEntryValues(std::ifstream& input, const TiffEntry& entry, bool little_endian,
                     std::vector<uint32_t>& values)
{
    std::vector<unsigned char> payload;
    if (!readEntryPayload(input, entry, little_endian, payload)) {
        return false;
    }

    values.clear();
    values.reserve(entry.count);
    switch (entry.type) {
    case kTypeShort:
        for (uint32_t i = 0; i < entry.count; ++i) {
            values.push_back(readU16(payload.data() + static_cast<size_t>(i) * 2U,
                                     little_endian));
        }
        return true;
    case kTypeLong:
        for (uint32_t i = 0; i < entry.count; ++i) {
            values.push_back(readU32(payload.data() + static_cast<size_t>(i) * 4U,
                                     little_endian));
        }
        return true;
    default:
        return false;
    }
}

const TiffEntry* findEntry(const std::vector<TiffEntry>& entries, uint16_t tag)
{
    for (const TiffEntry& entry : entries) {
        if (entry.tag == tag) {
            return &entry;
        }
    }
    return nullptr;
}

bool readScalar(std::ifstream& input, const std::vector<TiffEntry>& entries, uint16_t tag,
                bool little_endian, uint32_t& value, bool has_default = false,
                uint32_t default_value = 0U)
{
    const TiffEntry* entry = findEntry(entries, tag);
    if (entry == nullptr) {
        if (!has_default) {
            return false;
        }
        value = default_value;
        return true;
    }

    std::vector<uint32_t> values;
    if (!readEntryValues(input, *entry, little_endian, values) || values.size() != 1U) {
        return false;
    }
    value = values[0];
    return true;
}

bool readValues(std::ifstream& input, const std::vector<TiffEntry>& entries, uint16_t tag,
                bool little_endian, std::vector<uint32_t>& values)
{
    const TiffEntry* entry = findEntry(entries, tag);
    if (entry == nullptr) {
        return false;
    }
    return readEntryValues(input, *entry, little_endian, values);
}

} // namespace

bool writeGrayscaleTiff16(const char* filename, uint32_t width, uint32_t height,
                          const uint16_t* samples)
{
    if (filename == nullptr || samples == nullptr || width == 0U || height == 0U) {
        return false;
    }

    const uint64_t sample_count = static_cast<uint64_t>(width) * static_cast<uint64_t>(height);
    const uint64_t byte_count = sample_count * sizeof(uint16_t);
    if (byte_count > static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
        return false;
    }

    constexpr uint16_t entry_count = 11U;
    constexpr uint32_t ifd_offset = 8U;
    constexpr uint32_t pixel_offset = ifd_offset + 2U + entry_count * 12U + 4U;

    std::ofstream output(filename, std::ios::binary);
    if (!output) {
        return false;
    }

    output.write("II", 2);
    writeU16Le(output, kTiffVersion);
    writeU32Le(output, ifd_offset);

    writeU16Le(output, entry_count);
    writeEntryLong(output, kTagImageWidth, width);
    writeEntryLong(output, kTagImageLength, height);
    writeEntryShort(output, kTagBitsPerSample, 16U);
    writeEntryShort(output, kTagCompression, kCompressionNone);
    writeEntryShort(output, kTagPhotometricInterpretation, kPhotometricMinIsBlack);
    writeEntryLong(output, kTagStripOffsets, pixel_offset);
    writeEntryShort(output, kTagSamplesPerPixel, 1U);
    writeEntryLong(output, kTagRowsPerStrip, height);
    writeEntryLong(output, kTagStripByteCounts, static_cast<uint32_t>(byte_count));
    writeEntryShort(output, kTagPlanarConfiguration, kPlanarConfigurationContig);
    writeEntryShort(output, kTagSampleFormat, kSampleFormatUInt);
    writeU32Le(output, 0U);

    for (uint64_t index = 0; index < sample_count; ++index) {
        writeU16Le(output, samples[index]);
    }

    return static_cast<bool>(output);
}

bool readGrayscaleTiff16(const char* filename, std::vector<uint16_t>& samples,
                         uint32_t& width, uint32_t& height)
{
    samples.clear();
    width = 0U;
    height = 0U;

    if (filename == nullptr) {
        return false;
    }

    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        return false;
    }

    std::array<unsigned char, 8> header{};
    if (!readExact(input, header.data(), header.size())) {
        return false;
    }

    bool little_endian = false;
    if (header[0] == 'I' && header[1] == 'I') {
        little_endian = true;
    } else if (header[0] == 'M' && header[1] == 'M') {
        little_endian = false;
    } else {
        return false;
    }

    if (readU16(header.data() + 2U, little_endian) != kTiffVersion) {
        return false;
    }

    const uint32_t ifd_offset = readU32(header.data() + 4U, little_endian);
    if (ifd_offset == 0U) {
        return false;
    }

    input.seekg(static_cast<std::streamoff>(ifd_offset), std::ios::beg);
    std::array<unsigned char, 2> entry_count_bytes{};
    if (!readExact(input, entry_count_bytes.data(), entry_count_bytes.size())) {
        return false;
    }

    const uint16_t entry_count = readU16(entry_count_bytes.data(), little_endian);
    std::vector<TiffEntry> entries;
    entries.reserve(entry_count);
    for (uint16_t i = 0; i < entry_count; ++i) {
        std::array<unsigned char, 12> raw_entry{};
        if (!readExact(input, raw_entry.data(), raw_entry.size())) {
            return false;
        }

        TiffEntry entry;
        entry.tag = readU16(raw_entry.data(), little_endian);
        entry.type = readU16(raw_entry.data() + 2U, little_endian);
        entry.count = readU32(raw_entry.data() + 4U, little_endian);
        std::copy_n(raw_entry.begin() + 8, 4, entry.value_or_offset.begin());
        entries.push_back(entry);
    }

    std::array<unsigned char, 4> next_ifd_bytes{};
    if (!readExact(input, next_ifd_bytes.data(), next_ifd_bytes.size())) {
        return false;
    }

    uint32_t samples_per_pixel = 0U;
    uint32_t bits_per_sample = 0U;
    uint32_t compression = 0U;
    uint32_t photometric = 0U;
    uint32_t rows_per_strip = 0U;
    uint32_t planar_configuration = 0U;
    uint32_t sample_format = kSampleFormatUInt;
    if (!readScalar(input, entries, kTagImageWidth, little_endian, width) ||
        !readScalar(input, entries, kTagImageLength, little_endian, height) ||
        !readScalar(input, entries, kTagSamplesPerPixel, little_endian, samples_per_pixel, true,
                    1U) ||
        !readScalar(input, entries, kTagBitsPerSample, little_endian, bits_per_sample) ||
        !readScalar(input, entries, kTagCompression, little_endian, compression, true,
                    kCompressionNone) ||
        !readScalar(input, entries, kTagPhotometricInterpretation, little_endian, photometric) ||
        !readScalar(input, entries, kTagRowsPerStrip, little_endian, rows_per_strip, true,
                    height) ||
        !readScalar(input, entries, kTagPlanarConfiguration, little_endian, planar_configuration,
                    true, kPlanarConfigurationContig) ||
        !readScalar(input, entries, kTagSampleFormat, little_endian, sample_format, true,
                    kSampleFormatUInt)) {
        return false;
    }

    if (width == 0U || height == 0U || samples_per_pixel != 1U || bits_per_sample != 16U ||
        compression != kCompressionNone || photometric != kPhotometricMinIsBlack ||
        planar_configuration != kPlanarConfigurationContig ||
        sample_format != kSampleFormatUInt || rows_per_strip == 0U) {
        return false;
    }

    std::vector<uint32_t> strip_offsets;
    std::vector<uint32_t> strip_byte_counts;
    if (!readValues(input, entries, kTagStripOffsets, little_endian, strip_offsets) ||
        !readValues(input, entries, kTagStripByteCounts, little_endian, strip_byte_counts) ||
        strip_offsets.size() != strip_byte_counts.size()) {
        return false;
    }

    const uint32_t expected_strip_count = (height + rows_per_strip - 1U) / rows_per_strip;
    if (strip_offsets.size() != expected_strip_count) {
        return false;
    }

    const uint64_t sample_count = static_cast<uint64_t>(width) * static_cast<uint64_t>(height);
    if (sample_count > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        return false;
    }
    samples.resize(static_cast<size_t>(sample_count));

    uint32_t row_start = 0U;
    for (size_t strip_index = 0; strip_index < strip_offsets.size(); ++strip_index) {
        const uint32_t rows_in_strip = std::min(rows_per_strip, height - row_start);
        const uint64_t strip_sample_count =
            static_cast<uint64_t>(rows_in_strip) * static_cast<uint64_t>(width);
        const uint64_t expected_bytes = strip_sample_count * sizeof(uint16_t);
        if (strip_byte_counts[strip_index] < expected_bytes ||
            expected_bytes > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
            return false;
        }

        std::vector<unsigned char> strip_data(static_cast<size_t>(expected_bytes));
        input.seekg(static_cast<std::streamoff>(strip_offsets[strip_index]), std::ios::beg);
        if (!readExact(input, strip_data.data(), strip_data.size())) {
            return false;
        }

        const size_t destination_offset =
            static_cast<size_t>(row_start) * static_cast<size_t>(width);
        for (size_t i = 0; i < static_cast<size_t>(strip_sample_count); ++i) {
            samples[destination_offset + i] =
                readU16(strip_data.data() + i * sizeof(uint16_t), little_endian);
        }

        row_start += rows_in_strip;
    }

    return row_start == height;
}

} // namespace platec::io
