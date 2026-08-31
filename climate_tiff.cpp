#include "climate_tiff.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <limits>

namespace climateio
{
namespace
{
constexpr std::uint16_t tiffversion = 42;
constexpr std::uint16_t typeascii = 2;
constexpr std::uint16_t typeshort = 3;
constexpr std::uint16_t typelong = 4;
constexpr std::uint16_t typedouble = 12;

void writeu16(std::ostream& output, std::uint16_t value)
{
    const unsigned char bytes[] = {
        static_cast<unsigned char>(value & 0xffu),
        static_cast<unsigned char>((value >> 8) & 0xffu),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writeu32(std::ostream& output, std::uint32_t value)
{
    const unsigned char bytes[] = {
        static_cast<unsigned char>(value & 0xffu),
        static_cast<unsigned char>((value >> 8) & 0xffu),
        static_cast<unsigned char>((value >> 16) & 0xffu),
        static_cast<unsigned char>((value >> 24) & 0xffu),
    };
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void writeu64(std::ostream& output, std::uint64_t value)
{
    for (int shift = 0; shift < 64; shift += 8)
        output.put(static_cast<char>((value >> shift) & 0xffu));
}

void writefloat(std::ostream& output, float value)
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "float32 GeoTIFF requires 32-bit float");
    std::memcpy(&bits, &value, sizeof(bits));
    writeu32(output, bits);
}

void writedouble(std::ostream& output, double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "GeoTIFF requires 64-bit double");
    std::memcpy(&bits, &value, sizeof(bits));
    writeu64(output, bits);
}

void writeentry(std::ostream& output, std::uint16_t tag, std::uint16_t type,
    std::uint32_t count, std::uint32_t valueoroffset)
{
    writeu16(output, tag);
    writeu16(output, type);
    writeu32(output, count);
    writeu32(output, valueoroffset);
}

void writeentryshort(std::ostream& output, std::uint16_t tag, std::uint16_t value)
{
    writeu16(output, tag);
    writeu16(output, typeshort);
    writeu32(output, 1);
    writeu16(output, value);
    writeu16(output, 0);
}
}

bool writefloat32geotiff(const char* filename, std::uint32_t width,
    std::uint32_t height, const float* samples, float nodatavalue)
{
    if (filename == nullptr || samples == nullptr || width == 0 || height == 0 ||
        !std::isfinite(nodatavalue))
    {
        return false;
    }

    const std::uint64_t samplecount =
        static_cast<std::uint64_t>(width) * static_cast<std::uint64_t>(height);
    const std::uint64_t bytecount = samplecount * sizeof(float);

    if (bytecount > std::numeric_limits<std::uint32_t>::max())
        return false;

    constexpr std::uint16_t entrycount = 15;
    constexpr std::uint32_t ifdoffset = 8;
    constexpr std::uint32_t ifdend = ifdoffset + 2 + entrycount * 12 + 4;
    constexpr std::uint32_t pixelscaleoffset = ifdend;
    constexpr std::uint32_t tiepointoffset = pixelscaleoffset + 3 * sizeof(double);
    constexpr std::uint32_t geokeyoffset = tiepointoffset + 6 * sizeof(double);
    constexpr std::uint32_t nodataoffset = geokeyoffset + 16 * sizeof(std::uint16_t);
    constexpr std::uint32_t pixeloffset = nodataoffset + 8 + 2;

    std::ofstream output(filename, std::ios::binary);

    if (!output)
        return false;

    output.write("II", 2);
    writeu16(output, tiffversion);
    writeu32(output, ifdoffset);

    writeu16(output, entrycount);
    writeentry(output, 256, typelong, 1, width);
    writeentry(output, 257, typelong, 1, height);
    writeentryshort(output, 258, 32);
    writeentryshort(output, 259, 1);
    writeentryshort(output, 262, 1);
    writeentry(output, 273, typelong, 1, pixeloffset);
    writeentryshort(output, 277, 1);
    writeentry(output, 278, typelong, 1, height);
    writeentry(output, 279, typelong, 1, static_cast<std::uint32_t>(bytecount));
    writeentryshort(output, 284, 1);
    writeentryshort(output, 339, 3);
    writeentry(output, 33550, typedouble, 3, pixelscaleoffset);
    writeentry(output, 33922, typedouble, 6, tiepointoffset);
    writeentry(output, 34735, typeshort, 16, geokeyoffset);
    writeentry(output, 42113, typeascii, 8, nodataoffset);
    writeu32(output, 0);

    writedouble(output, 360.0 / static_cast<double>(width));
    writedouble(output, 180.0 / static_cast<double>(height));
    writedouble(output, 0.0);

    for (double value : {0.0, 0.0, 0.0, -180.0, 90.0, 0.0})
        writedouble(output, value);

    for (std::uint16_t value : {
             1, 1, 0, 3,
             1024, 0, 1, 2,
             1025, 0, 1, 1,
             2048, 0, 1, 4326,
         })
    {
        writeu16(output, value);
    }

    output.write("-9999.9", 8);
    writeu16(output, 0);

    for (std::uint64_t index = 0; index < samplecount; ++index)
    {
        const float value = std::isfinite(samples[index]) ? samples[index] : nodatavalue;
        writefloat(output, value);
    }

    return static_cast<bool>(output);
}
}
