#pragma once
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace climateatmosphere::detail
{
// Optional diagnostic fixture, independent of world/save formats. Decimal
// max_digits10 preserves every float input exactly on text roundtrip.
struct StationaryReplay
{
    int columns = 0, rows = 0;
    float depth = 0, damping = 0, density = 0, radius = 0, rotation = 0, direction = 0;
    bool rowProjection = true;
    int limit = 0;
    float tolerance = 0;
    int restart = 0;
    std::vector<float> forcing, drag;
    void write(const std::string& directory) const
    {
        if (directory.empty()) return;
        static std::atomic<unsigned> sequence{0};
        std::filesystem::create_directories(directory);
        const auto path = std::filesystem::path(directory) / ("stationary-" + std::to_string(sequence++) + ".txt");
        if (std::filesystem::exists(path)) throw std::runtime_error("Stationary capture already exists: " + path.string());
        std::ofstream out(path);
        out << std::setprecision(std::numeric_limits<float>::max_digits10)
            << "UW_STATIONARY_V1 " << columns << ' ' << rows << ' ' << depth << ' ' << damping << ' '
            << density << ' ' << radius << ' ' << rotation << ' ' << direction << ' ' << rowProjection << ' '
            << limit << ' ' << tolerance << ' ' << restart << '\n';
        for (std::size_t i = 0; i < forcing.size(); ++i) out << forcing[i] << ' ' << drag[i] << '\n';
        if (!out) throw std::runtime_error("Stationary capture write failed: " + path.string());
    }
    static StationaryReplay read(const std::string& path)
    {
        StationaryReplay f;
        std::ifstream in(path);
        std::string version;
        in >> version >> f.columns >> f.rows >> f.depth >> f.damping >> f.density >> f.radius
            >> f.rotation >> f.direction >> f.rowProjection >> f.limit >> f.tolerance >> f.restart;
        if (!in || version != "UW_STATIONARY_V1" || f.columns < 8 || f.columns > 2048 || f.rows != f.columns / 2)
            throw std::runtime_error("Invalid stationary replay header");
        f.forcing.resize(f.columns * f.rows); f.drag.resize(f.forcing.size());
        for (std::size_t i = 0; i < f.forcing.size(); ++i) in >> f.forcing[i] >> f.drag[i];
        if (!in) throw std::runtime_error("Truncated stationary replay");
        return f;
    }
};
}
