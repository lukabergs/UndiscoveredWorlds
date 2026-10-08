// Tests of landscape/upsampling.hpp (erosion plan Phase E, step E2a), suite Landscape.

#include <gtest/gtest.h>

#include "core/deterministic_random.hpp"
#include "landscape_fixture.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
// The window over global columns c0 .. c0 + columns - 1 (modulo W) and rows
// r0 .. r0 + rows - 1, from a box a quarter cell inside them.
landscape::RegionalWindow cellwindow(int width, int scale, int halo, int c0, int columns, int r0, int rows)
{
    const int height = width / 2;
    double west = -180.0 + 360.0 * (c0 + 0.25) / width;
    double east = -180.0 + 360.0 * (c0 + columns - 0.25) / width;
    if (west > 180.0)
        west -= 360.0;
    if (east > 180.0)
        east -= 360.0;
    const double north = 90.0 - 180.0 * (r0 + 0.25) / height;
    const double south = 90.0 - 180.0 * (r0 + rows - 0.25) / height;
    landscape::RegionalWindow window = landscape::regionalwindow(width, scale, halo, west, south, east, north);
    EXPECT_EQ(window.column0, c0 % width);
    EXPECT_EQ(window.columns, columns);
    EXPECT_EQ(window.row0, r0);
    EXPECT_EQ(window.rows, rows);
    return window;
}

std::vector<float> globalfield(int width, int height, double (*value)(int, int))
{
    std::vector<float> field(static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
            field[static_cast<std::size_t>(j) * static_cast<std::size_t>(width) + static_cast<std::size_t>(i)] =
                static_cast<float>(value(i, j));
    return field;
}

double smooth(int i, int j) { return 400.0 * std::sin(0.5 * i) * std::cos(0.45 * j) + 30.0 * j; }

// Per-fine-cell noise of +-10 m, a function of the absolute fine cell only.
float finenoise(std::int64_t column, std::int64_t row)
{
    const std::uint64_t hash = deterministiccombine(std::uint64_t{ 99 }, column, row);
    return static_cast<float>(20.0 * (static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0) - 0.5));
}

std::vector<float> noisyupsample(const std::vector<float>& field, int width, int height,
                                 const landscape::RegionalWindow& window)
{
    std::vector<float> fine = landscape::upsample(field, width, height, window);
    const int finewidth = window.domainfinewidth();
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x < finewidth; x++)
            fine[static_cast<std::size_t>(y) * static_cast<std::size_t>(finewidth) + static_cast<std::size_t>(x)] +=
                finenoise(window.absolutefinecolumn(x), window.absolutefinerow(y));
    return fine;
}

// The largest |block mean - global value| over the domain cells at least `margin`
// global cells from the domain's edges.
double largestblockerror(const std::vector<float>& field, int width, const landscape::RegionalWindow& window,
                         const std::vector<float>& fine, int margin)
{
    const int k = window.scale, finewidth = window.domainfinewidth();
    double largest = 0.0;
    for (int dr = margin; dr < window.domainrows - margin; dr++)
        for (int dc = margin; dc < window.domaincolumns - margin; dc++)
        {
            double sum = 0.0;
            for (int y = dr * k; y < (dr + 1) * k; y++)
                for (int x = dc * k; x < (dc + 1) * k; x++)
                    sum += fine[static_cast<std::size_t>(y) * static_cast<std::size_t>(finewidth) + static_cast<std::size_t>(x)];
            const int i = (window.domaincolumn0 + dc) % width, j = window.domainrow0 + dr;
            const double global = field[static_cast<std::size_t>(j) * static_cast<std::size_t>(width) + static_cast<std::size_t>(i)];
            largest = std::max(largest, std::fabs(sum / (k * k) - global));
        }
    return largest;
}

// Global cells between a global cell and the nearest domain edge that is not a pole.
int insidedistance(const landscape::RegionalWindow& window, int i, int j)
{
    const int relative = ((i - window.domaincolumn0) % window.worldwidth + window.worldwidth) % window.worldwidth;
    int distance = std::min(relative, window.domaincolumns - 1 - relative);
    if (window.domainrow0 > 0)
        distance = std::min(distance, j - window.domainrow0);
    if (window.domainrow0 + window.domainrows < window.worldheight)
        distance = std::min(distance, window.domainrow0 + window.domainrows - 1 - j);
    return distance;
}
}

TEST_F(Landscape, UpsamplingReproducesAConstantExactly)
{
    const std::vector<float> field(16 * 8, 123.25f);
    const landscape::RegionalWindow windows[] = {
        cellwindow(16, 4, 2, 3, 3, 2, 2),   // interior
        cellwindow(16, 4, 2, 14, 4, 3, 2),  // across the date line
        cellwindow(16, 3, 2, 5, 2, 0, 1),   // at the north pole
        cellwindow(16, 5, 1, 9, 2, 7, 1),   // at the south pole
    };
    for (const landscape::RegionalWindow& window : windows)
    {
        std::vector<float> fine = landscape::upsample(field, 16, 8, window);
        ASSERT_EQ(fine.size(), static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
        for (float value : fine)
            ASSERT_EQ(value, 123.25f);
        landscape::correctblockmeans(field, 16, 8, window, 2, fine);
        for (float value : fine)
            ASSERT_EQ(value, 123.25f);
    }
}

TEST_F(Landscape, UpsamplingReproducesALinearField)
{
    // 3 i + 2 j + 10 on a 16 x 8 world; the domain is columns 4 .. 11 and rows 2 .. 5,
    // so the stencil reaches neither the wrap nor a clamped row.
    constexpr int k = 4;
    const std::vector<float> field = globalfield(16, 8, [](int i, int j) { return 3.0 * i + 2.0 * j + 10.0; });
    const landscape::RegionalWindow window = cellwindow(16, k, 1, 5, 6, 3, 2);
    ASSERT_EQ(window.domaincolumn0, 4);
    ASSERT_EQ(window.domainrow0, 2);
    const std::vector<float> fine = landscape::upsample(field, 16, 8, window);
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x < window.domainfinewidth(); x++)
        {
            const double gx = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / k;
            const double gy = (static_cast<double>(window.absolutefinerow(y)) + 0.5) / k;
            const double expected = 3.0 * (gx - 0.5) + 2.0 * (gy - 0.5) + 10.0;
            EXPECT_NEAR(landscape::bicubic(field, 16, 8, gx, gy), expected, 1e-9) << x << ", " << y;
            EXPECT_NEAR(fine[static_cast<std::size_t>(y) * window.domainfinewidth() + x], expected, 1e-5) << x << ", " << y;
        }

    // Linear in the column only: every row, the polar ones included, away from the wrap.
    const std::vector<float> columns = globalfield(16, 8, [](int i, int) { return 3.0 * i + 10.0; });
    for (double gy : { 0.1, 0.5, 3.3, 7.6, 7.95 })
        for (double gx : { 2.0, 4.3, 8.5, 12.9 })
            EXPECT_NEAR(landscape::bicubic(columns, 16, 8, gx, gy), 3.0 * (gx - 0.5) + 10.0, 1e-9) << gx << ", " << gy;
}

TEST_F(Landscape, UpsamplingWrapsColumnsAndClampsRows)
{
    const std::vector<float> field = globalfield(16, 8, smooth);
    // The column index wraps: gx and gx + 16 are the same point (up to the rounding of gx + 16).
    for (double gx : { 0.1, 0.5, 15.7 })
        EXPECT_NEAR(landscape::bicubic(field, 16, 8, gx, 3.2), landscape::bicubic(field, 16, 8, gx + 16.0, 3.2), 1e-9);
    // At a cell centre the global value is returned exactly (the spline interpolates).
    EXPECT_EQ(landscape::bicubic(field, 16, 8, 3.5, 2.5), static_cast<double>(field[2 * 16 + 3]));
    EXPECT_EQ(landscape::bicubic(field, 16, 8, 0.5, 0.5), static_cast<double>(field[0]));
    EXPECT_EQ(landscape::bicubic(field, 16, 8, 15.5, 7.5), static_cast<double>(field[7 * 16 + 15]));
    // Rows beyond the poles are the polar rows: a field whose two rows nearest a pole
    // are equal is flat between the pole and the first row's centre.
    std::vector<float> polar = field;
    for (int i = 0; i < 16; i++)
        polar[16 + i] = polar[i];
    EXPECT_NEAR(landscape::bicubic(polar, 16, 8, 5.5, 0.1), field[5], 1e-4);
}

TEST_F(Landscape, UpsamplingBlockMeansKeepTheGlobalValues)
{
    // A smooth 32 x 16 field, upsampled at k = 4 with +-10 m of per-fine-cell noise,
    // on a domain away from the poles (rows 1 .. 12, columns 2 .. 17).
    constexpr int width = 32, height = 16;
    const std::vector<float> field = globalfield(width, height, smooth);
    const landscape::RegionalWindow window = cellwindow(width, 4, 4, 6, 8, 5, 4);
    ASSERT_EQ(window.domainrow0, 1);
    ASSERT_EQ(window.domaincolumn0, 2);
    std::vector<float> fine = noisyupsample(field, width, height, window);
    const double before = largestblockerror(field, width, window, fine, 3);
    landscape::correctblockmeans(field, width, height, window, 2, fine);
    const double after = largestblockerror(field, width, window, fine, 3);
    // Measured (E2a review): 8.47 m before, 0.443 m after 2 passes; the tolerance is 0.5 m.
    std::cout << "block-mean error over the interior: " << before << " m before, " << after << " m after 2 passes\n";
    EXPECT_GT(before, 5.0);
    EXPECT_LT(after, 0.5);
}

TEST_F(Landscape, UpsamplingIsIdenticalInOverlappingWindows)
{
    // On a 64 x 32 world, window a: columns 56 .. 3 with a halo of 6 (domain columns
    // 50 .. 9, rows 2 .. 19); window b: columns 60 .. 5 with a halo of 5 (domain
    // columns 55 .. 10, rows 5 .. 20). Both domains cross the date line and start
    // at different non-zero rows.
    constexpr int width = 64, height = 32, k = 4, passes = 2;
    const std::vector<float> field = globalfield(width, height, smooth);
    const landscape::RegionalWindow a = cellwindow(width, k, 6, 56, 12, 8, 6);
    const landscape::RegionalWindow b = cellwindow(width, k, 5, 60, 10, 10, 6);
    ASSERT_EQ(a.domaincolumn0, 50);
    ASSERT_EQ(a.domainrow0, 2);
    ASSERT_EQ(a.domainrows, 18);
    ASSERT_EQ(b.domaincolumn0, 55);
    ASSERT_EQ(b.domainrow0, 5);
    ASSERT_EQ(b.domainrows, 16);

    const std::vector<float> upa = landscape::upsample(field, width, height, a);
    const std::vector<float> upb = landscape::upsample(field, width, height, b);
    std::vector<float> fa = noisyupsample(field, width, height, a);
    std::vector<float> fb = noisyupsample(field, width, height, b);
    landscape::correctblockmeans(field, width, height, a, passes, fa);
    landscape::correctblockmeans(field, width, height, b, passes, fb);

    // The reach of the correction is 2 x passes global cells from a non-polar domain
    // edge: cells at least that far inside both domains are compared bit for bit.
    constexpr int reach = 2 * passes;
    const std::int64_t finewidth = static_cast<std::int64_t>(width) * k;
    int shared = 0, inside = 0, inside0 = 0, nearer = 0, nearerdiffer = 0;
    for (int y = 0; y < a.domainfineheight(); y++)
        for (int x = 0; x < a.domainfinewidth(); x++)
        {
            const std::int64_t column = a.absolutefinecolumn(x), row = a.absolutefinerow(y);
            const std::int64_t xb = ((column - static_cast<std::int64_t>(b.domaincolumn0) * k) % finewidth + finewidth) % finewidth;
            const std::int64_t yb = row - static_cast<std::int64_t>(b.domainrow0) * k;
            if (xb >= b.domainfinewidth() || yb < 0 || yb >= b.domainfineheight())
                continue;
            const std::size_t ia = static_cast<std::size_t>(y) * a.domainfinewidth() + x;
            const std::size_t ib = static_cast<std::size_t>(yb) * b.domainfinewidth() + static_cast<std::size_t>(xb);
            shared++;
            ASSERT_EQ(upa[ia], upb[ib]) << "upsample at fine cell " << column << ", " << row;
            const int i = a.globalcolumn(x), j = a.globalrow(y);
            const int distance = std::min(insidedistance(a, i, j), insidedistance(b, i, j));
            if (distance < reach)
            {
                nearer++;
                nearerdiffer += fa[ia] != fb[ib] ? 1 : 0;
                continue;
            }
            inside++;
            inside0 += i == 0 ? 1 : 0;
            ASSERT_EQ(fa[ia], fb[ib]) << "corrected at fine cell " << column << ", " << row;
        }
    std::cout << "shared fine cells " << shared << ", compared " << inside << ", nearer than the reach " << nearer
              << " of which " << nearerdiffer << " differ\n";
    // Shared: columns 55 .. 9 and rows 5 .. 19; inside both by 4: columns 59 .. 5
    // (across the date line) and rows 9 .. 15.
    EXPECT_EQ(shared, 19 * 15 * k * k);
    EXPECT_EQ(inside, 11 * 7 * k * k);
    EXPECT_EQ(inside0, 7 * k * k);
    // The reach is tight: some cells nearer an edge than 2 x passes do differ.
    EXPECT_GT(nearerdiffer, 0);
}

TEST_F(Landscape, UpsamplingThrowsOnBadInput)
{
    const std::vector<float> field(16 * 8, 1.0f);
    const landscape::RegionalWindow window = cellwindow(16, 2, 1, 3, 2, 2, 2);
    EXPECT_THROW(landscape::bicubic(field, 16, 7, 1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(landscape::upsample(field, 8, 16, window), std::invalid_argument);
    std::vector<float> fine = landscape::upsample(field, 16, 8, window);
    EXPECT_THROW(landscape::correctblockmeans(field, 16, 8, window, -1, fine), std::invalid_argument);
    fine.pop_back();
    EXPECT_THROW(landscape::correctblockmeans(field, 16, 8, window, 1, fine), std::invalid_argument);
}
