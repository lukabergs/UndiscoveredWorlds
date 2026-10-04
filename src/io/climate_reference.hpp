#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace climatereference
{
struct MonthlyGrid
{
    int width = 0;
    int height = 0;
    int monthCount = 0;
    std::string variable;
    std::vector<float> values;

    float value(int month, int x, int y) const;
};

// Reads a processed reference field: a version 1.0 .npy of little-endian float32,
// shape (12 or 1, width / 2, width), cell-centred, north row first. The file name
// must be the variable or begin with it and an underscore (u10m_monthly).
bool loadMonthlyGrid(
    const std::filesystem::path& path,
    const std::string& expectedVariable,
    MonthlyGrid& grid,
    std::string* errorMessage = nullptr);
}
