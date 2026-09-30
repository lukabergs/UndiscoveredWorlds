#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

namespace climateatmosphere::detail
{
// Invert the original centred pressure/divergence stencil after averaging its
// momentum mobilities in longitude. Fourier modes then have five latitude
// diagonals, including Coriolis cross terms and the original polar closures.
// This changes the preconditioner only; acceptance uses the original operator.
class StationaryPreconditioner
{
    using Complex = std::complex<double>;
    using Band = std::array<Complex, 5>;
    int columns_, rows_;
    bool valid_ = true;
    std::vector<Band> factors_;
    std::vector<Complex> field_, line_, scratch_;

    void transform(bool inverse)
    {
        constexpr double tau = 6.283185307179586476925286766559;
        const int n = columns_;
        if ((n & (n - 1)) == 0)
        {
            for (int i = 1, j = 0; i < n; ++i)
            {
                int bit = n >> 1;
                for (; j & bit; bit >>= 1) j ^= bit;
                j ^= bit;
                if (i < j) std::swap(line_[i], line_[j]);
            }
            for (int length = 2; length <= n; length *= 2)
            {
                const auto step = std::polar(1.0, (inverse ? tau : -tau) / length);
                for (int i = 0; i < n; i += length)
                {
                    Complex weight = 1.0;
                    for (int j = 0; j < length / 2; ++j)
                    {
                        const auto a = line_[i + j], b = weight * line_[i + j + length / 2];
                        line_[i + j] = a + b;
                        line_[i + j + length / 2] = a - b;
                        weight *= step;
                    }
                }
            }
        }
        else
        {
            // Preserve support for every valid grid width. Power-of-two grids
            // use FFT above; unusual widths use this deterministic DFT fallback.
            std::fill(scratch_.begin(), scratch_.end(), Complex{});
            for (int k = 0; k < n; ++k)
            {
                const auto step = std::polar(1.0, (inverse ? tau : -tau) * k / n);
                Complex weight = 1.0;
                for (int x = 0; x < n; ++x) { scratch_[k] += line_[x] * weight; weight *= step; }
            }
            line_.swap(scratch_);
        }
        if (inverse) for (auto& value : line_) value /= n;
    }

public:
    StationaryPreconditioner(int columns, int rows, const std::vector<float>& drag,
        double coupling, double density, double radius, double rotation, double direction, double pi)
        : columns_(columns), rows_(rows), factors_(columns * rows),
          field_(columns * rows), line_(columns), scratch_(columns)
    {
        std::vector<double> direct(rows), cross(rows), dx(rows), cosine(rows), cap(rows), span(rows);
        const double dy = pi * radius / rows;
        for (int y = 0; y < rows; ++y)
        {
            const double latitude = (90.0 - 180.0 * (y + 0.5) / rows) * pi / 180.0;
            cosine[y] = std::max(0.0, std::cos(latitude));
            cap[y] = std::max(0.02, std::abs(std::cos(latitude)));
            dx[y] = 2.0 * pi * radius * cap[y] / columns;
            span[y] = (std::min(rows - 1, y + 1) - std::max(0, y - 1)) * dy;
            const double f = 2.0 * rotation * direction * std::sin(latitude);
            for (int x = 0; x < columns; ++x)
                if (drag[y * columns + x] > 0.0f)
                {
                    const double r = 1.0 / drag[y * columns + x];
                    direct[y] += 100.0 / density * r / (r * r + f * f) / columns;
                    cross[y] += 100.0 / density * f / (r * r + f * f) / columns;
                }
        }
        constexpr double tau = 6.283185307179586476925286766559;
        for (int k = 0; k < columns; ++k)
        {
            auto* bands = factors_.data() + k * rows;
            const double wave = std::sin(tau * k / columns);
            for (int y = 0; y < rows; ++y)
            {
                const int n = std::max(0, y - 1), s = std::min(rows - 1, y + 1);
                const Complex derivative(0.0, wave / dx[y]);
                const auto add = [&](int target, Complex value) { bands[y][target - y + 2] += value; };
                add(y, 1.0 - coupling * direct[y] * derivative * derivative);
                add(n, -coupling * derivative * cross[y] / span[y]);
                add(s, coupling * derivative * cross[y] / span[y]);
                for (int side = 0; side < 2; ++side)
                {
                    const int v = side ? s : n;
                    const double weight = (side ? 1.0 : -1.0) * coupling * cosine[v] / (span[y] * cap[y]);
                    add(std::max(0, v - 1), weight * direct[v] / span[v]);
                    add(std::min(rows - 1, v + 1), -weight * direct[v] / span[v]);
                    add(v, -weight * cross[v] * Complex(0.0, wave / dx[v]));
                }
            }
            for (int i = 0; i < rows; ++i)
            {
                if (!std::isfinite(std::abs(bands[i][2])) || std::abs(bands[i][2]) < 1.0e-20)
                { valid_ = false; return; }
                for (int j = i + 1; j <= std::min(rows - 1, i + 2); ++j)
                {
                    auto& factor = bands[j][i - j + 2];
                    factor /= bands[i][2];
                    for (int c = i + 1; c <= std::min(rows - 1, i + 2); ++c)
                        bands[j][c - j + 2] -= factor * bands[i][c - i + 2];
                }
            }
        }
    }

    bool valid() const { return valid_; }
    void solve(std::vector<double>& values)
    {
        for (int y = 0; y < rows_; ++y)
        {
            for (int x = 0; x < columns_; ++x) line_[x] = values[y * columns_ + x];
            transform(false);
            for (int k = 0; k < columns_; ++k) field_[k * rows_ + y] = line_[k];
        }
        for (int k = 0; k < columns_; ++k)
        {
            auto* b = field_.data() + k * rows_;
            const auto* a = factors_.data() + k * rows_;
            for (int y = 0; y < rows_; ++y)
                for (int x = std::max(0, y - 2); x < y; ++x) b[y] -= a[y][x - y + 2] * b[x];
            for (int y = rows_; y-- > 0;)
            {
                for (int x = y + 1; x <= std::min(rows_ - 1, y + 2); ++x) b[y] -= a[y][x - y + 2] * b[x];
                b[y] /= a[y][2];
            }
        }
        for (int y = 0; y < rows_; ++y)
        {
            for (int k = 0; k < columns_; ++k) line_[k] = field_[k * rows_ + y];
            transform(true);
            for (int x = 0; x < columns_; ++x) values[y * columns_ + x] = line_[x].real();
        }
    }
};
}
