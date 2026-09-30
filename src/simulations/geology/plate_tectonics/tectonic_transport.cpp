#include "tectonic_transport.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tectonics::transport {
namespace {
struct Profile { double linear, quadratic; };
Profile profile(double ll, double l, double q, double r, double rr) {
    double lo = std::clamp((7 * (l + q) - ll - r) / 12, std::min(l, q), std::max(l, q));
    double hi = std::clamp((7 * (q + r) - l - rr) / 12, std::min(q, r), std::max(q, r));
    const double curvature = l - 2 * q + r;
    const bool smooth_extremum = curvature * (ll - 2 * l + q) > 0
        && curvature * (q - 2 * r + rr) > 0;
    // Retain resolved smooth extrema. Flattening every peak at each remap
    // would itself diffuse mountains even with a higher-order reconstruction.
    if (!smooth_extremum) {
        const double delta = hi - lo, six = 6 * q - 3 * (lo + hi);
        if ((hi - q) * (q - lo) <= 0) lo = hi = q;
        else if (delta * six > delta * delta) lo = 3 * q - 2 * hi;
        else if (delta * six < -delta * delta) hi = 3 * q - 2 * lo;
    }
    return {hi - lo, 3 * (lo + hi - 2 * q)};
}
double minimum(const Profile& p) {
    double result = -.5 * std::abs(p.linear) + p.quadratic / 6;
    if (p.quadratic > 0) {
        const double x = -.5 * p.linear / p.quadratic;
        if (std::abs(x) < .5) result = std::min(result, p.linear * x + p.quadratic * (x * x - 1. / 12));
    }
    return result;
}
std::vector<double> remap_impl(uint32_t width, uint32_t height,
    const std::vector<double>& values, const std::vector<Displacement>& displacement,
    const std::vector<double>* moment, std::vector<double>* moved_moment) {
    if (width < 2 || height < 2 || values.size() != static_cast<size_t>(width) * height
        || displacement.size() != values.size()) throw std::runtime_error("Invalid transport grid");
    if (moment && moment->size() != values.size()) throw std::runtime_error("Invalid material moment grid");
    if (moved_moment) moved_moment->assign(values.size(), 0);
    for (size_t i = 0; i < values.size(); ++i)
        if (!std::isfinite(values[i]) || values[i] < 0 || !std::isfinite(displacement[i].x)
            || !std::isfinite(displacement[i].y)
            || std::max(std::abs(displacement[i].x), std::abs(displacement[i].y)) > .95)
            throw std::runtime_error("Invalid transport value or displacement");
    if (moment) for (size_t i = 0; i < values.size(); ++i)
        if (!std::isfinite((*moment)[i]) || (*moment)[i] < 0
            || (values[i] == 0 && (*moment)[i] != 0)) throw std::runtime_error("Invalid material moment");
    const auto index = [&](int x, int y) {
        return static_cast<size_t>((y + static_cast<int>(height)) % static_cast<int>(height)) * width
            + (x + static_cast<int>(width)) % static_cast<int>(width);
    };
    std::vector<double> result(values.size(), 0);
    for (uint32_t y = 0; y < height; ++y) for (uint32_t x = 0; x < width; ++x) {
        const auto i = index(x, y);
        auto px = profile(values[index(static_cast<int>(x) - 2, y)], values[index(static_cast<int>(x) - 1, y)],
            values[i], values[index(x + 1, y)], values[index(x + 2, y)]);
        auto py = profile(values[index(x, static_cast<int>(y) - 2)], values[index(x, static_cast<int>(y) - 1)],
            values[i], values[index(x, y + 1)], values[index(x, y + 2)]);
        // The sum of the two polynomial minima bounds the entire cell, not
        // just its corners. Scaling zero-mean terms preserves the integral.
        const double extent = -minimum(px) - minimum(py);
        if (extent > values[i]) {
            const double factor = values[i] / extent;
            px.linear *= factor; px.quadratic *= factor;
            py.linear *= factor; py.quadratic *= factor;
        }
        const double tx = x + displacement[i].x, ty = y + displacement[i].y;
        const int ix = static_cast<int>(std::floor(tx)), iy = static_cast<int>(std::floor(ty));
        const double fx = tx - ix, fy = ty - iy;
        for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) {
            const double weight = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy);
            const double cx = dx ? .5 * (1 - fx) : -.5 * fx;
            const double cy = dy ? .5 * (1 - fy) : -.5 * fy;
            const double qx = dx ? (1 - fx) * (.5 - fx) / 3 : fx * (fx - .5) / 3;
            const double qy = dy ? (1 - fy) * (.5 - fy) / 3 : fy * (fy - .5) / 3;
            const double average = values[i] + px.linear * cx + px.quadratic * qx
                + py.linear * cy + py.quadratic * qy;
            const auto destination = index(ix + dx, iy + dy);
            const double transferred = weight * std::max(0.0, average); // Roundoff at zero only.
            result[destination] += transferred;
            if (moment && values[i] > 0)
                (*moved_moment)[destination] += transferred * ((*moment)[i] / values[i]);
        }
    }
    return result;
}
}
std::vector<double> remap(uint32_t width, uint32_t height,
    const std::vector<double>& values, const std::vector<Displacement>& displacement) {
    return remap_impl(width, height, values, displacement, nullptr, nullptr);
}
std::vector<double> remap_with_moment(uint32_t width, uint32_t height,
    const std::vector<double>& values, const std::vector<double>& moment,
    const std::vector<Displacement>& displacement, std::vector<double>& moved_moment) {
    return remap_impl(width, height, values, displacement, &moment, &moved_moment);
}
}
