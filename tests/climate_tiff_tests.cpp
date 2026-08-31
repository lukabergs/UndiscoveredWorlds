#include "climate_tiff.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>

namespace
{
struct tiffentry
{
    std::uint16_t type = 0;
    std::uint32_t count = 0;
    std::uint32_t valueoroffset = 0;
};

std::uint16_t readu16(std::istream& input)
{
    unsigned char bytes[2] = {};
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    return static_cast<std::uint16_t>(bytes[0]) |
        static_cast<std::uint16_t>(static_cast<std::uint16_t>(bytes[1]) << 8);
}

std::uint32_t readu32(std::istream& input)
{
    unsigned char bytes[4] = {};
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    return static_cast<std::uint32_t>(bytes[0]) |
        (static_cast<std::uint32_t>(bytes[1]) << 8) |
        (static_cast<std::uint32_t>(bytes[2]) << 16) |
        (static_cast<std::uint32_t>(bytes[3]) << 24);
}

float readfloat(std::istream& input)
{
    const std::uint32_t bits = readu32(input);
    float value = 0.0f;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

bool expect(bool condition, const std::string& message)
{
    if (!condition)
        std::cerr << message << '\n';
    return condition;
}
}

int main()
{
    const char* path = "climate_tiff_test.tif";
    const float samples[] = {
        1.25f,
        2.5f,
        std::numeric_limits<float>::quiet_NaN(),
        4.0f,
    };

    if (!expect(climateio::writefloat32geotiff(path, 2, 2, samples),
            "GeoTIFF writer returned failure"))
    {
        return 1;
    }

    std::ifstream input(path, std::ios::binary);
    char byteorder[2] = {};
    input.read(byteorder, sizeof(byteorder));
    const std::uint16_t version = readu16(input);
    const std::uint32_t ifdoffset = readu32(input);
    input.seekg(ifdoffset, std::ios::beg);
    const std::uint16_t entrycount = readu16(input);
    std::map<std::uint16_t, tiffentry> entries;

    for (std::uint16_t index = 0; index < entrycount; ++index)
    {
        const std::uint16_t tag = readu16(input);
        entries[tag] = {readu16(input), readu32(input), readu32(input)};
    }

    bool passed = true;
    passed &= expect(byteorder[0] == 'I' && byteorder[1] == 'I', "GeoTIFF is not little-endian");
    passed &= expect(version == 42, "GeoTIFF version is not 42");
    passed &= expect(entries[256].valueoroffset == 2, "GeoTIFF width is incorrect");
    passed &= expect(entries[257].valueoroffset == 2, "GeoTIFF height is incorrect");
    passed &= expect((entries[258].valueoroffset & 0xffffu) == 32, "GeoTIFF is not 32-bit");
    passed &= expect((entries[339].valueoroffset & 0xffffu) == 3, "GeoTIFF samples are not IEEE float");
    passed &= expect(entries.count(33550) == 1 && entries.count(33922) == 1,
        "GeoTIFF georeferencing tags are missing");
    passed &= expect(entries.count(34735) == 1, "GeoTIFF WGS84 key is missing");
    passed &= expect(entries.count(42113) == 1, "GeoTIFF no-data tag is missing");

    input.seekg(entries[273].valueoroffset, std::ios::beg);
    passed &= expect(std::fabs(readfloat(input) - 1.25f) < 0.0001f, "First sample changed");
    passed &= expect(std::fabs(readfloat(input) - 2.5f) < 0.0001f, "Second sample changed");
    passed &= expect(std::fabs(readfloat(input) - (-9999.9f)) < 0.01f, "NaN was not written as no-data");
    passed &= expect(std::fabs(readfloat(input) - 4.0f) < 0.0001f, "Fourth sample changed");

    input.close();
    std::remove(path);
    return passed ? 0 : 1;
}
