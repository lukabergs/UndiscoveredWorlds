#include "climate_reference.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <limits>

namespace climatereference
{
namespace
{
constexpr std::array<char, 8> expectedMagic = { 'U', 'W', 'C', 'L', 'I', 'M', '1', '\0' };

template <typename T>
bool readValue(std::ifstream& input, T& value)
{
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    return input.good();
}

void setError(std::string* errorMessage, const std::string& message)
{
    if (errorMessage != nullptr)
        *errorMessage = message;
}
}

float MonthlyGrid::value(int month, int x, int y) const
{
    if (month < 0 || month >= monthCount || x < 0 || x >= width || y < 0 || y >= height)
        return std::numeric_limits<float>::quiet_NaN();

    const size_t cellCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    const size_t index = static_cast<size_t>(month) * cellCount +
        static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
    return values[index];
}

bool loadMonthlyGrid(
    const std::filesystem::path& path,
    const std::string& expectedVariable,
    MonthlyGrid& grid,
    std::string* errorMessage)
{
    grid = {};
    std::ifstream input(path, std::ios::binary);

    if (!input.is_open())
    {
        setError(errorMessage, "file not found: " + path.string());
        return false;
    }

    std::array<char, 8> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));

    uint32_t version = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t monthCount = 0;
    std::array<char, 16> variable{};

    if (magic != expectedMagic ||
        !readValue(input, version) ||
        !readValue(input, width) ||
        !readValue(input, height) ||
        !readValue(input, monthCount))
    {
        setError(errorMessage, "invalid or truncated climate reference header: " + path.string());
        return false;
    }

    input.read(variable.data(), static_cast<std::streamsize>(variable.size()));

    if (!input.good() || version != 1 || width == 0 || height == 0 || monthCount != 12)
    {
        setError(errorMessage, "unsupported climate reference header: " + path.string());
        return false;
    }

    variable.back() = '\0';
    const std::string variableName(variable.data());

    if (variableName != expectedVariable)
    {
        setError(errorMessage, "unexpected variable in climate reference: " + variableName);
        return false;
    }

    const uint64_t valueCount = static_cast<uint64_t>(width) * static_cast<uint64_t>(height) * static_cast<uint64_t>(monthCount);

    if (valueCount > static_cast<uint64_t>(std::numeric_limits<size_t>::max() / sizeof(float)))
    {
        setError(errorMessage, "climate reference is too large: " + path.string());
        return false;
    }

    grid.width = static_cast<int>(width);
    grid.height = static_cast<int>(height);
    grid.monthCount = static_cast<int>(monthCount);
    grid.variable = variableName;
    grid.values.resize(static_cast<size_t>(valueCount));
    input.read(
        reinterpret_cast<char*>(grid.values.data()),
        static_cast<std::streamsize>(grid.values.size() * sizeof(float)));

    if (!input.good())
    {
        grid = {};
        setError(errorMessage, "truncated climate reference data: " + path.string());
        return false;
    }

    return true;
}
}
