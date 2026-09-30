#include "climate_reference.hpp"
#include "climate_experiment_input.hpp"

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

void writeFixture(
    const std::filesystem::path& path,
    const std::string& variable,
    bool truncateValues,
    uint32_t monthCount = 12,
    uint32_t width = 2,
    uint32_t height = 1)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    const std::array<char, 8> magic = { 'U', 'W', 'C', 'L', 'I', 'M', '1', '\0' };
    const uint32_t version = 1;
    std::array<char, 16> variableField{};
    std::copy_n(variable.c_str(), std::min(variable.size(), variableField.size() - 1), variableField.data());

    output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    writeValue(output, version);
    writeValue(output, width);
    writeValue(output, height);
    writeValue(output, monthCount);
    output.write(variableField.data(), static_cast<std::streamsize>(variableField.size()));

    const int valueCount = truncateValues ? 4 :
        static_cast<int>(width * height * monthCount);

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

    writeFixture(fixture, "prec_annual", false, 1);
    expect(climatereference::loadMonthlyGrid(fixture, "prec_annual", grid, &error), "annual grid must load");
    expect(grid.monthCount == 1 && grid.value(0, 1, 0) == 1.25f, "annual grid values must be readable");

    writeFixture(fixture, "tavg", true);
    expect(!climatereference::loadMonthlyGrid(fixture, "tavg", grid, &error), "truncated grid data must be rejected");

    writeFixture(fixture, "tavg", false, 12, 2, 2);
    expect(!climatereference::loadMonthlyGrid(fixture, "tavg", grid, &error),
        "a reference with an extra latitude row must be rejected");

    climatereference::MonthlyGrid seasonal;
    seasonal.width = 2; seasonal.height = 1; seasonal.monthCount = 12;
    seasonal.values.resize(24);
    for (int month = 0; month < 12; ++month)
        for (int x = 0; x < 2; ++x) seasonal.values[month * 2 + x] = static_cast<float>(month * 10 + x);
    expect(climateio::experiment::sample(&seasonal, 0, 1, 0.5f, 1, 0, -1.0f) == 16.0f,
        "reference interpolation must use January and April, not January and February");
    expect(climateio::experiment::sample(&seasonal, 3, 0, 0.5f, 0, 0, -1.0f) == 45.0f,
        "reference interpolation must wrap October to January");
    seasonal.values[6] = -9999.9f;
    expect(climateio::experiment::sample(&seasonal, 0, 1, 0.5f, 0, 0, -7.0f) == -7.0f,
        "missing reference coverage must retain the native field, not become zero");
    expect(climateio::experiment::sample(&seasonal, 0, 1, 0.0f, 0, 0, -7.0f) == 0.0f,
        "unused missing endpoint must not invalidate an exact seasonal sample");
    expect(climateio::experiment::sample(nullptr, 0, 1, 0.5f, 0, 0, -7.0f) == -7.0f,
        "disabled reference input must preserve the native value");
    expect(climateio::experiment::Inputs{}.field(2, "qa") == nullptr,
        "default input must not read reference data");

    const auto directory = std::filesystem::current_path() / "reference-experiment-test";
    std::filesystem::create_directories(directory / "2");
    const auto input = directory / "2" / "qa.uwclim";
    writeFixture(input, "qa", false);
    expect(climateio::experiment::Inputs(directory).field(2, "qa") != nullptr,
        "explicit matching reference input must load");
    writeFixture(input, "qa", false, 1);
    bool rejected = false;
    try { climateio::experiment::Inputs invalid(directory); }
    catch (const std::runtime_error&) { rejected = true; }
    expect(rejected, "annual data must not be silently used as seasonal forcing");

    std::error_code ignored;
    std::filesystem::remove(fixture, ignored);
    std::filesystem::remove(input, ignored);
    std::filesystem::remove(directory / "2", ignored);
    std::filesystem::remove(directory, ignored);
    return failures == 0 ? 0 : 1;
}
