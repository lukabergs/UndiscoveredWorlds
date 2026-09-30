#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tectonics {
// Reinterpret native torus coordinates on the application's cell-centred world
// raster. This is a sampling boundary, not a spherical tectonic simulation.
struct WorldGrid {
    std::uint32_t native_width, native_height, world_width, world_height;

    WorldGrid(std::uint32_t nw, std::uint32_t nh, std::uint32_t ww, std::uint32_t wh)
        : native_width(nw), native_height(nh), world_width(ww), world_height(wh) {
        if (!nw || !nh || !ww || !wh) throw std::invalid_argument("Empty tectonic grid");
        if (static_cast<std::uint64_t>(nw) * wh != static_cast<std::uint64_t>(nh) * ww)
            throw std::invalid_argument("Tectonic and world rasters must have the same aspect ratio");
    }
    double scale() const { return static_cast<double>(world_width) / native_width; }
    std::size_t nativeIndex(std::uint32_t x, std::uint32_t y) const {
        const auto sx = std::min(native_width - 1, static_cast<std::uint32_t>((x + .5) / scale()));
        const auto sy = std::min(native_height - 1, static_cast<std::uint32_t>((y + .5) / scale()));
        return static_cast<std::size_t>(sy) * native_width + sx;
    }
    float height(const float* native, std::uint32_t x, std::uint32_t y) const {
        const double sx = (x + .5) / scale() - .5;
        const double sy = std::clamp((y + .5) / scale() - .5, 0.0, static_cast<double>(native_height - 1));
        const int ix = static_cast<int>(std::floor(sx));
        const auto iy = static_cast<std::uint32_t>(std::floor(sy));
        const auto x0 = static_cast<std::uint32_t>((ix + static_cast<int>(native_width)) % native_width);
        const auto x1 = (x0 + 1) % native_width;
        const auto y1 = std::min(iy + 1, native_height - 1);
        const double fx = sx - ix, fy = sy - iy;
        const auto row = [&](std::uint32_t row_index) {
            const auto start = static_cast<std::size_t>(row_index) * native_width;
            return (1 - fx) * native[start + x0] + fx * native[start + x1];
        };
        return static_cast<float>((1 - fy) * row(iy) + fy * row(y1));
    }
    float worldCoordinate(float native_coordinate) const {
        return static_cast<float>((native_coordinate + .5) * scale() - .5);
    }
};
}
