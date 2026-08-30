#include "climate_reference.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}

template <typename T>
void writeValue(std::ofstream& output, const T& value)
{
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void writeFixture(const std::filesystem::path& path, const std::string& variable, bool truncateValues)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    const std::array<char, 8> magic = { 'U', 'W', 'C', 'L', 'I', 'M', '1', '\0' };
    const uint32_t version = 1;
    const uint32_t width = 2;
    const uint32_t height = 1;
    const uint32_t monthCount = 12;
    std::array<char, 16> variableField{};
    std::copy_n(variable.c_str(), std::min(variable.size(), variableField.size() - 1), variableField.data());

    output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    writeValue(output, version);
    writeValue(output, width);
    writeValue(output, height);
    writeValue(output, monthCount);
    output.write(variableField.data(), static_cast<std::streamsize>(variableField.size()));

    const int valueCount = truncateValues ? 4 : 24;

    for (int index = 0; index < valueCount; index++)
    {
        const float value = static_cast<float>(index) + 0.25f;
        writeValue(output, value);
    }
}
}

int main()
{
    const std::filesystem::path fixture = std::filesystem::current_path() / "climate-reference-test.uwclim";
    climatereference::MonthlyGrid grid;
    std::string error;

    writeFixture(fixture, "tavg", false);
    expect(climatereference::loadMonthlyGrid(fixture, "tavg", grid, &error), "valid monthly grid must load");
    expect(grid.width == 2 && grid.height == 1 && grid.monthCount == 12, "monthly grid dimensions must match its header");
    expect(grid.value(11, 1, 0) == 23.25f, "monthly grid values must preserve month-major ordering");
    expect(!climatereference::loadMonthlyGrid(fixture, "prec", grid, &error), "variable mismatch must be rejected");

    writeFixture(fixture, "tavg", true);
    expect(!climatereference::loadMonthlyGrid(fixture, "tavg", grid, &error), "truncated grid data must be rejected");

    std::error_code ignored;
    std::filesystem::remove(fixture, ignored);
    return failures == 0 ? 0 : 1;
}
