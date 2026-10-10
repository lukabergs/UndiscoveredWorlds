// Tests of landscape/refinement_projection.hpp (docs/architecture/multiscale-plan.md,
// section 3, "Tests"), suite Landscape: Landscape.RefinementProjection*. S0.4 shows that
// a child and a parent window built by hand (refinement_test_windows.hpp) find each other's
// cells, which is all the projection needs from refinement_grid.hpp; S2.2 adds exactness,
// the weights, the border step, the bit identity in two windows and the resample.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"
#include "refinement_test_windows.hpp"
#include "regional_test_world.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_projection.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
using landscape::LevelWindow;

std::uint64_t mixkey(std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// The test's "detail": a deterministic perturbation of +-40 m, a function of the absolute
// key and the level only.
double detail(std::uint64_t key, int level)
{
    const double unit = static_cast<double>(mixkey(key ^ (0x51ULL * static_cast<std::uint64_t>(level + 1))) >> 11) * 0x1.0p-53;
    return 40.0 * (2.0 * unit - 1.0);
}

// A channel cell by hash of the absolute key (about one cell in seven).
std::vector<std::uint8_t> channelmask(const LevelWindow& window)
{
    std::vector<std::uint8_t> mask(window.domainsize(), 0);
    for (int y = 0; y < window.domainrows; y++)
        for (int x = 0; x < window.domaincolumns; x++)
            mask[window.index(x, y)] = mixkey(window.absolutekey(x, y) + 7u) % 7u == 0u ? 1 : 0;
    return mask;
}

std::vector<std::uint8_t> seamask(const std::vector<float>& field)
{
    std::vector<std::uint8_t> mask(field.size(), 0);
    for (std::size_t cell = 0; cell < field.size(); cell++)
        mask[cell] = field[cell] <= 0.0f ? 1 : 0;
    return mask;
}

// The level-0 window's domain values of a global W x W/2 field.
std::vector<float> globalvalues(const LevelWindow& window, const std::vector<float>& global)
{
    std::vector<float> field(window.domainsize());
    for (int y = 0; y < window.domainrows; y++)
        for (int x = 0; x < window.domaincolumns; x++)
            field[window.index(x, y)] = global[static_cast<std::size_t>(window.absoluterow(y)) * static_cast<std::size_t>(window.worldwidth) +
                                               static_cast<std::size_t>(window.absolutecolumn(x))];
    return field;
}

// Bilinear prolongation in the parent's domain (clamped at its edges) plus the detail.
std::vector<float> bilinear(const LevelWindow& child, const LevelWindow& parent, const std::vector<float>& parentfield)
{
    const int f = child.factor;
    const int column0 = parent.localcolumn(child.parentcolumn(0));
    const int row0 = static_cast<int>(parent.localrow(child.parentrow(0)));
    auto value = [&](int i, int j) {
        i = std::clamp(i, 0, parent.domaincolumns - 1);
        j = std::clamp(j, 0, parent.domainrows - 1);
        return static_cast<double>(parentfield[parent.index(i, j)]);
    };
    std::vector<float> field(child.domainsize());
    for (int y = 0; y < child.domainrows; y++)
        for (int x = 0; x < child.domaincolumns; x++)
        {
            const double u = column0 + (x + 0.5) / f - 0.5, v = row0 + (y + 0.5) / f - 0.5;
            const int i = static_cast<int>(std::floor(u)), j = static_cast<int>(std::floor(v));
            const double fu = u - i, fv = v - j;
            const double e = (1.0 - fu) * (1.0 - fv) * value(i, j) + fu * (1.0 - fv) * value(i + 1, j) +
                             (1.0 - fu) * fv * value(i, j + 1) + fu * fv * value(i + 1, j + 1);
            field[child.index(x, y)] = static_cast<float>(e + detail(child.absolutekey(x, y), child.level));
        }
    return field;
}

// Nearest prolongation (each child takes its parent's value) plus the detail: a function of
// the absolute position and the parent value only, so two windows agree on it.
std::vector<float> nearest(const LevelWindow& child, const LevelWindow& parent, const std::vector<float>& parentfield)
{
    std::vector<float> field(child.domainsize());
    for (int y = 0; y < child.domainrows; y++)
        for (int x = 0; x < child.domaincolumns; x++)
        {
            const int i = parent.localcolumn(child.parentcolumn(x));
            const int j = static_cast<int>(parent.localrow(child.parentrow(y)));
            field[child.index(x, y)] = static_cast<float>(static_cast<double>(parentfield[parent.index(i, j)]) +
                                                          detail(child.absolutekey(x, y), child.level));
        }
    return field;
}

void printreport(const char* label, int factor, int level, const landscape::ProjectionReport& report)
{
    std::printf("[projection] %s, factor %d, level %d: %d parents, parent-mean error %.3g m, %d all-pinned, "
                "lambda max %.3f m, lambda step max %.3f m, border step max %.3f m, 99th percentile %.3f m\n",
                label, factor, level, report.parents, report.maxparentmeanerror_m, report.allpinnedparents, report.maxlambda_m,
                report.maxlambdastep_m, report.maxborderstep_m, report.borderstep99_m);
}

// The border step recomputed from the projected field (d = F - F*, cells with a positive
// effective weight), for a cross-check of the report: the largest and the nearest-rank 99th
// percentile.
std::pair<double, double> borderstepof(const LevelWindow& level, const std::vector<float>& before, const std::vector<float>& after,
                                       const std::vector<float>& weights)
{
    const int f = level.factor;
    std::vector<std::uint8_t> moved(level.domainsize(), 0);
    for (int q = 0; q < level.domainparentrows(); q++)
        for (int p = 0; p < level.domainparentcolumns(); p++)
        {
            double sum = 0.0;
            for (int y = q * f; y < (q + 1) * f; y++)
                for (int x = p * f; x < (p + 1) * f; x++)
                    sum += weights[level.index(x, y)];
            for (int y = q * f; y < (q + 1) * f; y++)
                for (int x = p * f; x < (p + 1) * f; x++)
                    moved[level.index(x, y)] = sum == 0.0 || weights[level.index(x, y)] > 0.0f ? 1 : 0;
        }
    std::vector<double> steps;
    const int x0 = level.keptcolumnoffset(), y0 = level.keptrowoffset();
    for (int y = y0; y < y0 + level.rows; y++)
        for (int x = x0; x < x0 + level.columns; x++)
        {
            const std::size_t a = level.index(x, y);
            const double da = static_cast<double>(after[a]) - before[a];
            const int dx[2] = {1, 0}, dy[2] = {0, 1};
            for (int k = 0; k < 2; k++)
            {
                const int nx = x + dx[k], ny = y + dy[k];
                if (nx >= x0 + level.columns || ny >= y0 + level.rows)
                    continue;
                if ((dx[k] == 1 && (nx - x0) % f != 0) || (dy[k] == 1 && (ny - y0) % f != 0))
                    continue;
                const std::size_t b = level.index(nx, ny);
                if (moved[a] == 0 || moved[b] == 0)
                    continue;
                steps.push_back(std::fabs(da - (static_cast<double>(after[b]) - before[b])));
            }
        }
    if (steps.empty())
        return {0.0, 0.0};
    std::sort(steps.begin(), steps.end());
    const std::size_t n = steps.size();
    const std::size_t rank = (99 * n + 99) / 100 - 1;  // ceil(0.99 n) - 1
    return {steps.back(), steps[rank]};
}

// Levels 0 .. 2 on the 32-column synthetic world (the level-0 kept box spans -1676 ..
// 2799 m), the continent across the date line: level 0 keeps global columns 26 .. 37 (mod
// 32) by rows 2 .. 13; level 1 keeps 28 .. 35 by 4 .. 11 with a halo of 1 parent cell;
// level 2 keeps 29 .. 34 by 5 .. 10 with a halo of 2 parent cells. Each level is the
// bilinear prolongation of the projected level above plus the detail, projected with
// weights 0 on the sea and the hashed channels.
struct Chain
{
    std::vector<LevelWindow> windows;
    std::vector<std::vector<float>> fields;
    std::vector<std::vector<float>> candidates;  // F* of levels 1 and 2 (index = level)
    std::vector<std::vector<float>> weights;
    std::vector<landscape::ProjectionReport> reports;  // index = level, [0] unused
    std::vector<float> global;
    double before[3] = {0.0, 0.0, 0.0};  // parent-mean errors before the projection
};

Chain buildchain(int f)
{
    Chain chain;
    chain.global = syntheticworld(32).elevation_m;
    chain.windows.push_back(handlevelwindow(32, f, 0, 0, 26, 12, 2, 12));
    chain.windows.push_back(handlevelwindow(32, f, 1, 1, 28 * f, 8 * f, 4 * f, 8 * f));
    chain.windows.push_back(handlevelwindow(32, f, 2, 2, 29 * f * f, 6 * f * f, 5 * f * f, 6 * f * f));
    chain.fields.push_back(globalvalues(chain.windows[0], chain.global));
    chain.candidates.emplace_back();
    chain.weights.emplace_back();
    chain.reports.emplace_back();
    for (int level = 1; level <= 2; level++)
    {
        const LevelWindow& window = chain.windows[static_cast<std::size_t>(level)];
        const LevelWindow& parent = chain.windows[static_cast<std::size_t>(level - 1)];
        const std::vector<float>& parentfield = chain.fields[static_cast<std::size_t>(level - 1)];
        std::vector<float> field = bilinear(window, parent, parentfield);
        chain.before[level] = landscape::levelparentmeanerror(window, field, parent, parentfield);
        chain.candidates.push_back(field);
        chain.weights.push_back(landscape::projectionweights(window, channelmask(window), seamask(field), {}, 0.0));
        chain.reports.push_back(landscape::projectlevel(window, field, chain.weights.back(), parent, parentfield));
        chain.fields.push_back(std::move(field));
    }
    return chain;
}

void checkexactness(int f)
{
    const Chain chain = buildchain(f);
    float lowest = std::numeric_limits<float>::max(), highest = std::numeric_limits<float>::lowest();
    for (const float value : chain.fields[0])
    {
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
    }
    std::printf("[projection] factor %d: level-0 kept box %.0f .. %.0f m\n", f, lowest, highest);
    EXPECT_GT(highest - lowest, 3000.0f);
    for (int level = 1; level <= 2; level++)
    {
        const landscape::ProjectionReport& report = chain.reports[static_cast<std::size_t>(level)];
        const LevelWindow& window = chain.windows[static_cast<std::size_t>(level)];
        printreport("synthetic world", f, level, report);
        std::printf("[projection] factor %d, level %d: parent-mean error before the projection %.3f m\n", f, level,
                    chain.before[level]);
        EXPECT_GT(chain.before[level], 0.1) << "the candidates must miss the parent means before the projection";
        EXPECT_LT(report.maxparentmeanerror_m, 1.0e-3);
        EXPECT_EQ(report.parents, (window.columns / f) * (window.rows / f));
        EXPECT_EQ(report.maxparentmeanerror_m,
                  landscape::levelparentmeanerror(window, chain.fields[static_cast<std::size_t>(level)],
                                                  chain.windows[static_cast<std::size_t>(level - 1)],
                                                  chain.fields[static_cast<std::size_t>(level - 1)]));
        EXPECT_GT(report.maxlambda_m, 0.0);
        EXPECT_LE(report.borderstep99_m, report.maxborderstep_m);
        const std::pair<double, double> step =
            borderstepof(window, chain.candidates[static_cast<std::size_t>(level)], chain.fields[static_cast<std::size_t>(level)],
                         chain.weights[static_cast<std::size_t>(level)]);
        EXPECT_NEAR(report.maxborderstep_m, step.first, 1.0e-2);
        EXPECT_NEAR(report.borderstep99_m, step.second, 1.0e-2);
    }

    // The level-2 field averaged over each global cell of its kept box is the global value.
    const LevelWindow& fine = chain.windows[2];
    const std::vector<float>& field = chain.fields[2];
    const int s = f * f;
    double largest = 0.0;
    for (int gj = 0; gj < fine.rows / s; gj++)
        for (int gi = 0; gi < fine.columns / s; gi++)
        {
            double sum = 0.0;
            for (int y = 0; y < s; y++)
                for (int x = 0; x < s; x++)
                    sum += field[fine.index(fine.keptcolumnoffset() + gi * s + x, fine.keptrowoffset() + gj * s + y)];
            const int column = fine.globalcolumn(fine.keptcolumnoffset() + gi * s);
            const int row = fine.globalrow(fine.keptrowoffset() + gj * s);
            const double global = chain.global[static_cast<std::size_t>(row) * 32u + static_cast<std::size_t>(column)];
            largest = std::max(largest, std::fabs(sum / (static_cast<double>(s) * s) - global));
        }
    std::printf("[projection] factor %d: level 2 box-averaged to level 0, largest error %.3g m\n", f, largest);
    EXPECT_LT(largest, 1.0e-3);
}

bool samebits(float a, float b) { return std::memcmp(&a, &b, sizeof(float)) == 0; }
}

// A 16-column world at factor 2: the child at level 2 (64 x 32 cells) crosses the date
// line; the parent at level 1 (32 x 16 cells) keeps every parent cell of the child's
// domain. The parent's field holds each cell's absolute key, so the value a child reads
// names the cell it reached.
TEST_F(Landscape, RefinementProjectionHandBuiltWindowsFindTheirParents)
{
    const landscape::LevelWindow child = handlevelwindow(16, 2, 2, 1, 60, 8, 8, 4);
    const landscape::LevelWindow parent = handlevelwindow(16, 2, 1, 1, 28, 8, 2, 6);
    ASSERT_EQ(child.domaincolumn0, 58);
    ASSERT_EQ(child.domaincolumns, 12);
    ASSERT_EQ(child.domainrow0, 6);
    ASSERT_EQ(child.domainrows, 8);
    ASSERT_EQ(parent.domaincolumn0, 26);
    ASSERT_EQ(parent.domainrows, 10);

    std::vector<float> parentfield(parent.domainsize());
    for (int y = 0; y < parent.domainrows; y++)
        for (int x = 0; x < parent.domaincolumns; x++)
            parentfield[parent.index(x, y)] = static_cast<float>(parent.absolutekey(x, y));

    std::map<std::pair<std::int64_t, std::int64_t>, int> children;
    for (int y = 0; y < child.domainrows; y++)
        for (int x = 0; x < child.domaincolumns; x++)
        {
            const std::int64_t column = child.parentcolumn(x), row = child.parentrow(y);
            // The domain's p-th parent column is the parent of its cells factor p ..
            // factor p + factor - 1.
            EXPECT_EQ(column, landscape::wraplevel(child.domaincolumn0 / child.factor + x / child.factor,
                                                   parent.levelwidth()));
            EXPECT_EQ(row, child.domainrow0 / child.factor + y / child.factor);
            // Every parent cell of the child's domain lies in the parent's kept box.
            ASSERT_TRUE(parent.containscolumn(column));
            ASSERT_TRUE(parent.containsrow(row));
            const int px = parent.localcolumn(column);
            const int py = static_cast<int>(parent.localrow(row));
            EXPECT_GE(px, parent.keptcolumnoffset());
            EXPECT_LT(px, parent.keptcolumnoffset() + parent.columns);
            EXPECT_GE(py, parent.keptrowoffset());
            EXPECT_LT(py, parent.keptrowoffset() + parent.rows);
            EXPECT_EQ(parentfield[parent.index(px, py)], static_cast<float>(landscape::levelkey(column, row, 32)));
            EXPECT_EQ(child.globalcolumn(x), parent.globalcolumn(px));
            EXPECT_EQ(child.globalrow(y), parent.globalrow(py));
            children[{column, row}]++;
        }
    // Six parent columns (29 .. 31 and 0 .. 2) by four rows, each with factor^2 children.
    EXPECT_EQ(children.size(), static_cast<std::size_t>(child.domainparentcolumns() * child.domainparentrows()));
    EXPECT_EQ(children.size(), 24u);
    for (const auto& entry : children)
        EXPECT_EQ(entry.second, 4);
    EXPECT_EQ(children.count({29, 3}), 1u);
    EXPECT_EQ(children.count({2, 6}), 1u);
}

// Exactness on the synthetic world: level 1 projected on level 0, level 2 on level 1, the
// largest parent-mean error below 1e-3 m at each, and the level-2 field box-averaged to
// level 0 equal to the global values to 1e-3 m.
TEST_F(Landscape, RefinementProjectionExactFactor2) { checkexactness(2); }

TEST_F(Landscape, RefinementProjectionExactFactor4) { checkexactness(4); }

// The weights: 0 on sea, lake beds and (at weight 0) channel beds, which the projection
// leaves unchanged bit for bit; an all-sea parent takes the uniform fallback and is counted.
TEST_F(Landscape, RefinementProjectionWeights)
{
    const int f = 4;
    const LevelWindow parent = handlevelwindow(16, f, 0, 0, 4, 2, 2, 1);
    const LevelWindow level = handlevelwindow(16, f, 1, 0, 4 * f, 2 * f, 2 * f, f);
    ASSERT_EQ(level.domainsize(), 32u);
    std::vector<float> parentfield = {100.0f, -50.0f};
    std::vector<float> field(level.domainsize());
    std::vector<std::uint8_t> channel(level.domainsize(), 0), sea(level.domainsize(), 0), lake(level.domainsize(), 0);
    for (int y = 0; y < f; y++)
        for (int x = 0; x < 2 * f; x++)
        {
            const std::size_t cell = level.index(x, y);
            field[cell] = x < f ? static_cast<float>(80.0 + 3.0 * x - 2.0 * y) : static_cast<float>(-70.0 + x + y);
            sea[cell] = x >= f ? 1 : 0;
        }
    channel[level.index(1, 0)] = channel[level.index(1, 1)] = channel[level.index(2, 2)] = 1;
    sea[level.index(3, 3)] = 1;
    channel[level.index(3, 3)] = 1;  // sea and channel: 0
    channel[level.index(5, 1)] = 1;  // sea and channel: 0
    lake[level.index(0, 3)] = 1;

    const std::vector<float> zero = landscape::projectionweights(level, channel, sea, lake, 0.0);
    const std::vector<float> half = landscape::projectionweights(level, channel, sea, lake, 0.5);
    for (std::size_t cell = 0; cell < level.domainsize(); cell++)
    {
        const bool pinned = sea[cell] != 0 || lake[cell] != 0;
        EXPECT_EQ(zero[cell], pinned || channel[cell] != 0 ? 0.0f : 1.0f) << cell;
        EXPECT_EQ(half[cell], pinned ? 0.0f : (channel[cell] != 0 ? 0.5f : 1.0f)) << cell;
    }
    EXPECT_EQ(landscape::projectionweights(level, {}, {}, {}, 0.0), std::vector<float>(level.domainsize(), 1.0f));
    EXPECT_THROW(landscape::projectionweights(level, std::vector<std::uint8_t>(31, 0), {}, {}, 0.0), std::invalid_argument);
    EXPECT_THROW(landscape::projectionweights(level, {}, std::vector<std::uint8_t>(33, 0), {}, 0.0), std::invalid_argument);
    EXPECT_THROW(landscape::projectionweights(level, {}, {}, std::vector<std::uint8_t>(1, 0), 0.0), std::invalid_argument);
    EXPECT_THROW(landscape::projectionweights(level, {}, {}, {}, -0.1), std::invalid_argument);
    EXPECT_THROW(landscape::projectionweights(level, {}, {}, {}, std::nan("")), std::invalid_argument);
    EXPECT_THROW(landscape::projectionweights(level, {}, {}, {}, std::numeric_limits<double>::infinity()),
                 std::invalid_argument);

    const std::vector<float> before = field;
    const landscape::ProjectionReport report = landscape::projectlevel(level, field, zero, parent, parentfield);
    printreport("weights", f, 1, report);
    EXPECT_EQ(report.parents, 2);
    EXPECT_EQ(report.allpinnedparents, 1);
    EXPECT_LT(report.maxparentmeanerror_m, 1.0e-4);
    int unchanged = 0, changed = 0;
    for (int y = 0; y < f; y++)
        for (int x = 0; x < f; x++)
        {
            const std::size_t cell = level.index(x, y);
            if (zero[cell] == 0.0f)
            {
                EXPECT_TRUE(samebits(field[cell], before[cell])) << "pinned cell " << x << ", " << y;
                unchanged++;
            }
            else
            {
                EXPECT_FALSE(samebits(field[cell], before[cell])) << "land cell " << x << ", " << y;
                changed++;
            }
        }
    EXPECT_EQ(unchanged, 5);
    EXPECT_EQ(changed, 11);
    // The all-sea parent: every child lifted by the same lambda.
    double lift = static_cast<double>(field[level.index(f, 0)]) - before[level.index(f, 0)];
    for (int y = 0; y < f; y++)
        for (int x = f; x < 2 * f; x++)
            EXPECT_NEAR(static_cast<double>(field[level.index(x, y)]) - before[level.index(x, y)], lift, 1.0e-4);
}

TEST_F(Landscape, RefinementProjectionRejectsBadInput)
{
    const int f = 2;
    const LevelWindow parent = handlevelwindow(16, f, 0, 0, 2, 5, 1, 5);
    const LevelWindow level = handlevelwindow(16, f, 1, 1, 3 * f, 3 * f, 2 * f, 3 * f);
    std::vector<float> field(level.domainsize(), 1.0f);
    const std::vector<float> weights(level.domainsize(), 1.0f), parentfield(parent.domainsize(), 2.0f);
    EXPECT_NO_THROW(landscape::projectlevel(level, field, weights, parent, parentfield));

    std::vector<float> bad = weights;
    bad[3] = -1.0f;
    EXPECT_THROW(landscape::projectlevel(level, field, bad, parent, parentfield), std::invalid_argument);
    bad[3] = std::numeric_limits<float>::quiet_NaN();
    EXPECT_THROW(landscape::projectlevel(level, field, bad, parent, parentfield), std::invalid_argument);
    bad[3] = std::numeric_limits<float>::infinity();
    EXPECT_THROW(landscape::projectlevel(level, field, bad, parent, parentfield), std::invalid_argument);
    EXPECT_THROW(landscape::projectlevel(level, field, std::vector<float>(3, 1.0f), parent, parentfield), std::invalid_argument);
    std::vector<float> shortfield(level.domainsize() - 1, 1.0f);
    EXPECT_THROW(landscape::projectlevel(level, shortfield, weights, parent, parentfield), std::invalid_argument);
    EXPECT_THROW(landscape::levelparentmeanerror(level, shortfield, parent, parentfield), std::invalid_argument);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, parent, std::vector<float>(2, 0.0f)), std::invalid_argument);
    // Level 0 has no parent; the parent must be one level coarser, of the same factor and world.
    EXPECT_THROW(landscape::levelparentmeanerror(parent, parentfield, parent, parentfield), std::invalid_argument);
    const LevelWindow levelparent = handlevelwindow(16, f, 1, 0, 2 * f, 5 * f, 1 * f, 5 * f);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, levelparent, std::vector<float>(levelparent.domainsize(), 0.0f)),
                 std::invalid_argument);
    const LevelWindow otherfactor = handlevelwindow(16, 4, 0, 0, 2, 5, 1, 5);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, otherfactor, parentfield), std::invalid_argument);
    const LevelWindow otherworld = handlevelwindow(32, f, 0, 0, 2, 5, 1, 5);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, otherworld, parentfield), std::invalid_argument);
    // A parent kept box that misses a parent cell of the level's domain (east, then south).
    const LevelWindow narrow = handlevelwindow(16, f, 0, 0, 2, 4, 1, 5);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, narrow, std::vector<float>(narrow.domainsize(), 0.0f)),
                 std::invalid_argument);
    const LevelWindow low = handlevelwindow(16, f, 0, 0, 2, 5, 1, 4);
    EXPECT_THROW(landscape::levelparentmeanerror(level, field, low, std::vector<float>(low.domainsize(), 0.0f)),
                 std::invalid_argument);
    // A parent with a halo: its domain is wider than its kept box, which is what counts.
    const LevelWindow haloed = handlevelwindow(16, f, 1, 1, 2 * f, 4 * f, 1 * f, 5 * f);
    const LevelWindow deeper = handlevelwindow(16, f, 2, 0, 3 * f * f, 2 * f * f, 2 * f * f, 3 * f * f);
    EXPECT_NO_THROW(landscape::levelparentmeanerror(deeper, std::vector<float>(deeper.domainsize(), 0.0f), haloed,
                                                    std::vector<float>(haloed.domainsize(), 0.0f)));
    const LevelWindow wide = handlevelwindow(16, f, 2, 3, 3 * f * f, 2 * f * f, 2 * f * f, 3 * f * f);
    EXPECT_THROW(landscape::levelparentmeanerror(wide, std::vector<float>(wide.domainsize(), 0.0f), haloed,
                                                 std::vector<float>(haloed.domainsize(), 0.0f)),
                 std::invalid_argument);
}

// The border step against section 2's sizes (I = 300 m): a one-cell channel line of
// incision I across a parent of factor^2 children lifts its land by I / (factor - 1); a
// single incised cell by I / (factor^2 - 1); a neighbouring parent with a lift of its own
// sees the difference. The level window keeps 3 x 3 parent cells (global columns 3 .. 5,
// rows 2 .. 4) with a halo of one parent cell whose parents carry 1000 m, which no
// statistic may see.
TEST_F(Landscape, RefinementProjectionBorderStep)
{
    const double incision = 300.0;
    for (const int f : {2, 4})
    {
        const LevelWindow parent = handlevelwindow(16, f, 0, 0, 2, 5, 1, 5);
        const LevelWindow level = handlevelwindow(16, f, 1, 1, 3 * f, 3 * f, 2 * f, 3 * f);
        ASSERT_EQ(level.domainparentcolumns(), 5);
        ASSERT_EQ(level.domainparentrows(), 5);
        // Pairs across the kept box's inner parent edges: 2 x 3f east-west and 2 x 3f
        // north-south.
        const int pairs = 12 * f;

        auto parentvalues = [&](double kept) {
            std::vector<float> values(parent.domainsize(), 1000.0f);
            for (int y = 1; y <= 3; y++)
                for (int x = 1; x <= 3; x++)
                    values[parent.index(x, y)] = static_cast<float>(kept);
            return values;
        };
        auto expectedpercentile = [&](int n, int nonzero, double step) {
            const int rank = (99 * n + 99) / 100 - 1;
            return rank >= n - nonzero ? step : 0.0;
        };

        // The line: child column 0 of the centre parent (domain parent 2, 2), every row.
        for (const double delta : {0.0, 20.0})
        {
            const std::vector<float> parentfield = parentvalues(delta);
            std::vector<float> field(level.domainsize(), 0.0f);
            std::vector<std::uint8_t> channel(level.domainsize(), 0);
            for (int y = 2 * f; y < 3 * f; y++)
            {
                field[level.index(2 * f, y)] = static_cast<float>(-incision);
                channel[level.index(2 * f, y)] = 1;
            }
            const std::vector<float> weights = landscape::projectionweights(level, channel, {}, {}, 0.0);
            const std::vector<float> before = field;
            const landscape::ProjectionReport report = landscape::projectlevel(level, field, weights, parent, parentfield);
            printreport(delta == 0.0 ? "channel line" : "channel line beside lifted parents", f, 1, report);
            const double lambda = delta + incision / f;
            const double lift = lambda * f / (f - 1.0);
            if (delta == 0.0)
                EXPECT_NEAR(lift, incision / (f - 1.0), 1.0e-9);
            EXPECT_NEAR(report.maxlambda_m, lambda, 1.0e-4);
            EXPECT_NEAR(report.maxlambdastep_m, incision / f, 1.0e-4);
            EXPECT_NEAR(report.maxborderstep_m, lift - delta, 1.0e-4);
            // The pairs touching the centre parent: f east, f - 1 north and south (the line's
            // cell is pinned), none west.
            EXPECT_NEAR(report.borderstep99_m, expectedpercentile(pairs - f, 3 * f - 2, lift - delta), 1.0e-4);
            EXPECT_LT(report.maxparentmeanerror_m, 1.0e-4);
            EXPECT_EQ(report.allpinnedparents, 0);
            for (int y = 2 * f; y < 3 * f; y++)
            {
                EXPECT_TRUE(samebits(field[level.index(2 * f, y)], before[level.index(2 * f, y)]));
                EXPECT_NEAR(field[level.index(2 * f + 1, y)], lift, 1.0e-4);
            }
        }

        // A single incised cell inside the centre parent.
        {
            const std::vector<float> parentfield = parentvalues(0.0);
            std::vector<float> field(level.domainsize(), 0.0f);
            std::vector<std::uint8_t> channel(level.domainsize(), 0);
            field[level.index(2 * f + 1, 2 * f + 1)] = static_cast<float>(-incision);
            channel[level.index(2 * f + 1, 2 * f + 1)] = 1;
            const std::vector<float> weights = landscape::projectionweights(level, channel, {}, {}, 0.0);
            const landscape::ProjectionReport report = landscape::projectlevel(level, field, weights, parent, parentfield);
            printreport("single incised cell", f, 1, report);
            const double lift = incision / (f * f - 1.0);
            EXPECT_NEAR(report.maxlambda_m, incision / (f * f), 1.0e-4);
            EXPECT_NEAR(report.maxborderstep_m, lift, 1.0e-4);
            // At factor 2 the cell is the parent's south-east corner: its east and south
            // pairs are pinned.
            const int pinned = f == 2 ? 2 : 0;
            EXPECT_NEAR(report.borderstep99_m, expectedpercentile(pairs - pinned, 4 * f - pinned, lift), 1.0e-4);
        }
    }
}

// Two windows of different column and row origin, window A across the date line, project
// levels 1 and 2 at factor 2 (nearest prolongation plus the detail, both functions of the
// absolute position): every cell the two domains share is bit-identical at both levels.
TEST_F(Landscape, RefinementProjectionBitIdenticalInTwoWindows)
{
    const int f = 2;
    const std::vector<float> global = syntheticworld(32).elevation_m;
    struct Windows
    {
        LevelWindow w0, w1, w2;
        std::vector<float> f0, f1, f2;
    };
    auto run = [&](const LevelWindow& w0, const LevelWindow& w1, const LevelWindow& w2) {
        Windows windows{w0, w1, w2, globalvalues(w0, global), {}, {}};
        windows.f1 = nearest(w1, w0, windows.f0);
        landscape::projectlevel(w1, windows.f1, landscape::projectionweights(w1, channelmask(w1), seamask(windows.f1), {}, 0.0),
                                w0, windows.f0);
        windows.f2 = nearest(w2, w1, windows.f1);
        const landscape::ProjectionReport report = landscape::projectlevel(
            w2, windows.f2, landscape::projectionweights(w2, channelmask(w2), seamask(windows.f2), {}, 0.0), w1, windows.f1);
        EXPECT_LT(report.maxparentmeanerror_m, 1.0e-3);
        return windows;
    };
    // A keeps global columns 28 .. 35 (mod 32) by rows 3 .. 9 at level 2; B keeps 0 .. 7 by 5 .. 11.
    const Windows a = run(handlevelwindow(32, f, 0, 0, 26, 12, 1, 11), handlevelwindow(32, f, 1, 1, 54, 20, 4, 18),
                          handlevelwindow(32, f, 2, 2, 112, 32, 12, 28));
    const Windows b = run(handlevelwindow(32, f, 0, 0, 30, 12, 3, 11), handlevelwindow(32, f, 1, 1, 62, 20, 8, 18),
                          handlevelwindow(32, f, 2, 2, 0, 32, 20, 28));

    auto compare = [&](const LevelWindow& wa, const std::vector<float>& fa, const LevelWindow& wb, const std::vector<float>& fb,
                       int level) {
        int shared = 0, differing = 0;
        for (int y = 0; y < wa.domainrows; y++)
            for (int x = 0; x < wa.domaincolumns; x++)
            {
                const std::int64_t column = wa.absolutecolumn(x), row = wa.absoluterow(y);
                if (!wb.containscolumn(column) || !wb.containsrow(row))
                    continue;
                shared++;
                const float va = fa[wa.index(x, y)];
                const float vb = fb[wb.index(wb.localcolumn(column), static_cast<int>(wb.localrow(row)))];
                differing += samebits(va, vb) ? 0 : 1;
            }
        std::printf("[projection] two windows, level %d: %d shared cells, %d differ\n", level, shared, differing);
        EXPECT_GT(shared, 100);
        EXPECT_EQ(differing, 0) << "level " << level;
    };
    compare(a.w1, a.f1, b.w1, b.f1, 1);
    compare(a.w2, a.f2, b.w2, b.f2, 2);
}

// The resample, aligned: at factor 4 a pixel of 2 level-2 cells divides factor^2, so each
// pixel is the exact mean of its 2 x 2 cells and a global cell's 8 x 8 pixels average to
// the global value; a pixel of one level-1 cell is chosen at level 1 and returns its cells.
TEST_F(Landscape, RefinementProjectionResampleAligned)
{
    const int f = 4;
    const Chain chain = buildchain(f);
    const double cell2 = pi * radius_m / static_cast<double>(chain.windows[2].levelheight());
    const LevelWindow& w2 = chain.windows[2];

    const landscape::ResampledField two = landscape::resampleto(chain.windows, chain.fields, radius_m, 2.0 * cell2);
    EXPECT_EQ(two.report.sourcelevel, 2);
    EXPECT_TRUE(two.report.aligned);
    EXPECT_FALSE(two.report.resampled());
    EXPECT_NEAR(two.report.pixelcells, 2.0, 1.0e-12);
    EXPECT_NEAR(two.report.sourcecell_m, cell2, 1.0e-9);
    EXPECT_EQ(two.report.metresperpixel, 2.0 * cell2);
    ASSERT_EQ(two.width, w2.columns / 2);
    ASSERT_EQ(two.height, w2.rows / 2);
    int differing = 0;
    for (int q = 0; q < two.height; q++)
        for (int p = 0; p < two.width; p++)
        {
            double sum = 0.0;
            for (int y = 2 * q; y < 2 * q + 2; y++)
                for (int x = 2 * p; x < 2 * p + 2; x++)
                    sum += chain.fields[2][w2.index(w2.keptcolumnoffset() + x, w2.keptrowoffset() + y)];
            differing += samebits(two.values[static_cast<std::size_t>(q * two.width + p)], static_cast<float>(sum / 4.0)) ? 0 : 1;
        }
    EXPECT_EQ(differing, 0);
    double largest = 0.0;
    for (int gj = 0; gj < 6; gj++)
        for (int gi = 0; gi < 6; gi++)
        {
            double sum = 0.0;
            for (int q = 8 * gj; q < 8 * gj + 8; q++)
                for (int p = 8 * gi; p < 8 * gi + 8; p++)
                    sum += two.values[static_cast<std::size_t>(q * two.width + p)];
            const double global = chain.global[static_cast<std::size_t>(5 + gj) * 32u + static_cast<std::size_t>((29 + gi) % 32)];
            largest = std::max(largest, std::fabs(sum / 64.0 - global));
        }
    std::printf("[projection] resample at 2 level-2 cells: %d x %d pixels, global means within %.3g m\n", two.width,
                two.height, largest);
    EXPECT_LT(largest, 1.0e-3);

    const landscape::ResampledField one = landscape::resampleto(chain.windows, chain.fields, radius_m, 4.0 * cell2);
    EXPECT_EQ(one.report.sourcelevel, 1);
    EXPECT_TRUE(one.report.aligned);
    const LevelWindow& w1 = chain.windows[1];
    ASSERT_EQ(one.width, 6 * f);
    ASSERT_EQ(one.height, 6 * f);
    // The target's kept box at level 1 starts one global cell into level 1's kept box.
    for (int q = 0; q < one.height; q++)
        for (int p = 0; p < one.width; p++)
            EXPECT_TRUE(samebits(one.values[static_cast<std::size_t>(q * one.width + p)],
                                 chain.fields[1][w1.index(w1.keptcolumnoffset() + f + p, w1.keptrowoffset() + f + q)]));
}

// The resample, not aligned: a whole n that does not divide factor^level (3 at factor 4),
// an unaligned n (2.5) with the overlap weights, and a whole n of 2 at level 0, exact but
// reported as resampled.
TEST_F(Landscape, RefinementProjectionResampleUnaligned)
{
    const int f = 4;
    const Chain chain = buildchain(f);
    const LevelWindow& w2 = chain.windows[2];
    const double cell2 = pi * radius_m / static_cast<double>(w2.levelheight());
    auto source = [&](int i, int j) {
        return static_cast<double>(chain.fields[2][w2.index(w2.keptcolumnoffset() + i, w2.keptrowoffset() + j)]);
    };

    const landscape::ResampledField three = landscape::resampleto(chain.windows, chain.fields, radius_m, 3.0 * cell2);
    EXPECT_EQ(three.report.sourcelevel, 2);
    EXPECT_FALSE(three.report.aligned);
    EXPECT_TRUE(three.report.resampled());
    ASSERT_EQ(three.width, w2.columns / 3);
    ASSERT_EQ(three.height, w2.rows / 3);
    double largest = 0.0;
    for (int q = 0; q < three.height; q++)
        for (int p = 0; p < three.width; p++)
        {
            double sum = 0.0;
            for (int y = 3 * q; y < 3 * q + 3; y++)
                for (int x = 3 * p; x < 3 * p + 3; x++)
                    sum += source(x, y);
            largest = std::max(largest, std::fabs(three.values[static_cast<std::size_t>(q * three.width + p)] - sum / 9.0));
        }
    EXPECT_LT(largest, 1.0e-3);

    // n = 2.5: each source cell split into 2 x 2 half cells, a pixel is the mean of its
    // 5 x 5 half cells, computed independently of the overlap lengths.
    const landscape::ResampledField unaligned = landscape::resampleto(chain.windows, chain.fields, radius_m, 2.5 * cell2);
    EXPECT_EQ(unaligned.report.sourcelevel, 2);
    EXPECT_FALSE(unaligned.report.aligned);
    EXPECT_NEAR(unaligned.report.pixelcells, 2.5, 1.0e-12);
    ASSERT_EQ(unaligned.width, static_cast<int>(w2.columns / 2.5));
    ASSERT_EQ(unaligned.height, static_cast<int>(w2.rows / 2.5));
    largest = 0.0;
    for (int q = 0; q < unaligned.height; q++)
        for (int p = 0; p < unaligned.width; p++)
        {
            double sum = 0.0;
            for (int b = 5 * q; b < 5 * q + 5; b++)
                for (int a = 5 * p; a < 5 * p + 5; a++)
                    sum += source(a / 2, b / 2);
            largest = std::max(largest, std::fabs(unaligned.values[static_cast<std::size_t>(q * unaligned.width + p)] - sum / 25.0));
        }
    std::printf("[projection] resample at 2.5 level-2 cells: %d x %d pixels, largest difference from the half-cell means %.3g m\n",
                unaligned.width, unaligned.height, largest);
    EXPECT_LT(largest, 1.0e-3);

    const double cell0 = pi * radius_m / 16.0;
    const landscape::ResampledField coarse = landscape::resampleto(chain.windows, chain.fields, radius_m, 2.0 * cell0);
    EXPECT_EQ(coarse.report.sourcelevel, 0);
    EXPECT_FALSE(coarse.report.aligned);
    ASSERT_EQ(coarse.width, 3);
    ASSERT_EQ(coarse.height, 3);
    for (int q = 0; q < 3; q++)
        for (int p = 0; p < 3; p++)
        {
            double sum = 0.0;
            for (int y = 0; y < 2; y++)
                for (int x = 0; x < 2; x++)
                    sum += chain.global[static_cast<std::size_t>(5 + 2 * q + y) * 32u + static_cast<std::size_t>((29 + 2 * p + x) % 32)];
            EXPECT_NEAR(coarse.values[static_cast<std::size_t>(q * 3 + p)], sum / 4.0, 1.0e-3);
        }
}

// The source level is the coarsest whose cell is at most the pixel (to 1e-9); the
// exceptions.
TEST_F(Landscape, RefinementProjectionResampleSourceLevelAndErrors)
{
    const int f = 4;
    const Chain chain = buildchain(f);
    double cell[3];
    for (int l = 0; l < 3; l++)
        cell[l] = pi * radius_m / static_cast<double>(chain.windows[static_cast<std::size_t>(l)].levelheight());
    auto sourcelevel = [&](double metres) {
        return landscape::resampleto(chain.windows, chain.fields, radius_m, metres).report.sourcelevel;
    };
    EXPECT_EQ(sourcelevel(1.5 * cell[0]), 0);
    EXPECT_EQ(sourcelevel(cell[0]), 0);
    EXPECT_EQ(sourcelevel(cell[0] * (1.0 - 1.0e-12)), 0);
    EXPECT_EQ(sourcelevel(cell[0] * 0.999), 1);
    EXPECT_EQ(sourcelevel(cell[1]), 1);
    EXPECT_EQ(sourcelevel(cell[1] * 0.999), 2);
    EXPECT_EQ(sourcelevel(cell[2]), 2);
    const landscape::ResampledField finest = landscape::resampleto(chain.windows, chain.fields, radius_m, cell[2]);
    EXPECT_TRUE(finest.report.aligned);
    EXPECT_EQ(finest.width, chain.windows[2].columns);
    const landscape::ResampledField near = landscape::resampleto(chain.windows, chain.fields, radius_m, cell[1] * 0.999);
    EXPECT_FALSE(near.report.aligned);
    EXPECT_NEAR(near.report.pixelcells, 3.996, 1.0e-9);

    const std::vector<landscape::LevelWindow>& windows = chain.windows;
    const std::vector<std::vector<float>>& fields = chain.fields;
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, 0.0), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, -cell[0]), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, std::nan("")), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, std::numeric_limits<double>::infinity()), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, 0.0, cell[0]), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, std::nan(""), cell[0]), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto({}, {}, radius_m, cell[0]), std::invalid_argument);
    // No level fine enough, no pixel (6 global cells, a pixel of 7).
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, 0.5 * cell[2]), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, 7.0 * cell[0]), std::invalid_argument);
    // An absurdly large pixel: no pixel, rejected before n is taken as an integer.
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, 1.0e300), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, fields, radius_m, std::numeric_limits<double>::max()), std::invalid_argument);
    // A target kept box that is not whole cells of the source level: level 2 shifted by a
    // quarter of a global cell reads at level 2 but not at level 0.
    {
        std::vector<landscape::LevelWindow> shifted = windows;
        std::vector<std::vector<float>> shiftedfields = fields;
        shifted[2] = handlevelwindow(32, f, 2, 2, 29 * f * f + f, 6 * f * f, 5 * f * f, 6 * f * f);
        shiftedfields[2].assign(shifted[2].domainsize(), 0.0f);
        EXPECT_NO_THROW(landscape::resampleto(shifted, shiftedfields, radius_m, cell[2]));
        EXPECT_THROW(landscape::resampleto(shifted, shiftedfields, radius_m, cell[0]), std::invalid_argument);
    }
    // A source level whose kept box misses the target's: level 1 kept one global cell east
    // (30 .. 37 instead of 28 .. 35), so the target's column 29 is outside it.
    {
        std::vector<landscape::LevelWindow> moved = windows;
        std::vector<std::vector<float>> movedfields = fields;
        moved[1] = handlevelwindow(32, f, 1, 1, 30 * f, 8 * f, 4 * f, 8 * f);
        movedfields[1].assign(moved[1].domainsize(), 0.0f);
        EXPECT_NO_THROW(landscape::resampleto(moved, movedfields, radius_m, cell[2]));
        EXPECT_THROW(landscape::resampleto(moved, movedfields, radius_m, cell[1]), std::invalid_argument);
    }
    // Malformed: levels out of order, a field of the wrong size, fields missing.
    std::vector<landscape::LevelWindow> swapped = windows;
    std::swap(swapped[0], swapped[1]);
    std::vector<std::vector<float>> swappedfields = fields;
    std::swap(swappedfields[0], swappedfields[1]);
    EXPECT_THROW(landscape::resampleto(swapped, swappedfields, radius_m, cell[0]), std::invalid_argument);
    std::vector<std::vector<float>> shortfields = fields;
    shortfields[1].pop_back();
    EXPECT_THROW(landscape::resampleto(windows, shortfields, radius_m, cell[0]), std::invalid_argument);
    EXPECT_THROW(landscape::resampleto(windows, std::vector<std::vector<float>>(fields.begin(), fields.end() - 1), radius_m, cell[0]),
                 std::invalid_argument);
}
