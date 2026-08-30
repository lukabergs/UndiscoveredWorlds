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

bool loadMonthlyGrid(
    const std::filesystem::path& path,
    const std::string& expectedVariable,
    MonthlyGrid& grid,
    std::string* errorMessage = nullptr);
}
