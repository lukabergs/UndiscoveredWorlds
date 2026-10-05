#include "io/climate_reference.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>

namespace climatereference
{
namespace
{
void setError(std::string* errorMessage, const std::string& message)
{
    if (errorMessage != nullptr)
        *errorMessage = message;
}

// Reads the (layers, rows, columns) shape of a version 1.0 little-endian float32 C-order .npy header.
bool readNpyShape(std::ifstream& input, std::array<uint64_t, 3>& shape)
{
    std::array<unsigned char, 10> prefix{};
    input.read(reinterpret_cast<char*>(prefix.data()), static_cast<std::streamsize>(prefix.size()));

    if (!input.good() || std::string(reinterpret_cast<const char*>(prefix.data()), 6) != "\x93NUMPY" ||
        prefix[6] != 1 || prefix[7] != 0)
        return false;

    std::string header(static_cast<size_t>(prefix[8]) | (static_cast<size_t>(prefix[9]) << 8), ' ');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));

    if (!input.good() ||
        header.find("'descr': '<f4'") == std::string::npos ||
        header.find("'fortran_order': False") == std::string::npos)
        return false;

    const size_t open = header.find("'shape': (");
    const size_t close = open == std::string::npos ? open : header.find(')', open);

    if (close == std::string::npos)
        return false;

    std::string extents = header.substr(open + 10, close - open - 10);

    for (char& character : extents)
        if (character == ',')
            character = ' ';

    std::istringstream values(extents);

    for (uint64_t& extent : shape)
        if (!(values >> extent))
            return false;

    uint64_t extra = 0;
    return !(values >> extra);
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

    // Field names begin with their variable: u10m_monthly, prec_annual_2001_2022 (refs/README.md).
    const std::string name = path.stem().string();

    if (name != expectedVariable && name.rfind(expectedVariable + "_", 0) != 0)
    {
        setError(errorMessage, "unexpected variable in climate reference: " + name);
        return false;
    }

    std::ifstream input(path, std::ios::binary);

    if (!input.is_open())
    {
        setError(errorMessage, "file not found: " + path.string());
        return false;
    }

    std::array<uint64_t, 3> shape{};

    if (!readNpyShape(input, shape))
    {
        setError(errorMessage, "invalid or truncated climate reference header: " + path.string());
        return false;
    }

    const uint64_t monthCount = shape[0];
    const uint64_t height = shape[1];
    const uint64_t width = shape[2];

    if (width < 2 || width % 2 != 0 || width > static_cast<uint64_t>(std::numeric_limits<int>::max()) ||
        height != width / 2 ||
        (monthCount != 1 && monthCount != 12))
    {
        setError(errorMessage, "unsupported cell-centred climate reference shape: " + path.string());
        return false;
    }

    const uint64_t valueCount = width * height * monthCount;

    if (valueCount > static_cast<uint64_t>(std::numeric_limits<size_t>::max() / sizeof(float)))
    {
        setError(errorMessage, "climate reference is too large: " + path.string());
        return false;
    }

    grid.width = static_cast<int>(width);
    grid.height = static_cast<int>(height);
    grid.monthCount = static_cast<int>(monthCount);
    grid.variable = expectedVariable;
    grid.values.resize(static_cast<size_t>(valueCount));
    input.read(
        reinterpret_cast<char*>(grid.values.data()),
        static_cast<std::streamsize>(grid.values.size() * sizeof(float)));

    if (!input.good() || input.peek() != std::ifstream::traits_type::eof())
    {
        grid = {};
        setError(errorMessage, "climate reference data do not match the header shape: " + path.string());
        return false;
    }

    return true;
}
}
