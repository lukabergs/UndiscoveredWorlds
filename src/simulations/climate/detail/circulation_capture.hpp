#pragma once
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace climateatmosphere::detail
{
// Optional diagnostic only: unrounded inputs and responses, in cell-centred
// row-major order. Existing UW_STATIONARY_V1 replay/save formats are untouched.
inline void captureCirculation(const std::string& directory, int columns, int rows,
    int season, const std::vector<std::pair<std::string, std::vector<float>>>& fields)
{
    if (directory.empty()) return;
    for (const auto& field : fields)
        if (field.second.size() != static_cast<std::size_t>(columns * rows))
            throw std::runtime_error("Invalid circulation diagnostic dimensions");
    static std::atomic<unsigned> sequence{0};
    std::filesystem::create_directories(directory);
    const auto path = std::filesystem::path(directory) /
        ("circulation-" + std::to_string(sequence++) + "-season-" + std::to_string(season) + ".csv");
    if (std::filesystem::exists(path)) throw std::runtime_error("Circulation diagnostic already exists");
    std::ofstream out(path);
    out << std::setprecision(std::numeric_limits<float>::max_digits10) << "x,y";
    for (const auto& field : fields) out << ',' << field.first;
    out << '\n';
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            out << x << ',' << y;
            for (const auto& field : fields) out << ',' << field.second[y * columns + x];
            out << '\n';
        }
    if (!out) throw std::runtime_error("Circulation diagnostic write failed");
}
}
