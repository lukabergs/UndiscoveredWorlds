#include "simulations/landscape/refinement_projection.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/upsampling.hpp"

// The projection, its error measure and the resample (step S2.2,
// docs/architecture/multiscale-plan.md; contracts in the header, written in S0.4; the
// bounds that replaced S2.2's weights are S2.3's rework, plan section 3, "S2.3's projection
// weights and detail levels"; the smoothing before the projection is S2.7's, "S2.7: smooth
// before exact"). A
// parent's children are visited in row-major order of their absolute position: the
// domain of a level window is whole parent cells and its columns never split a parent
// across the date line (the level width is a multiple of factor), so the local order
// inside a parent is the absolute order every window shares.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;

// The parent windows' agreement and the sizes of projectlevel and levelparentmeanerror.
void checkpair(const LevelWindow& level, std::size_t fieldsize, const LevelWindow& parent, std::size_t parentsize,
               const char* name)
{
    const std::string prefix = std::string(name) + ": ";
    if (level.level < 1 || parent.level != level.level - 1)
        throw std::invalid_argument(prefix + "the level must be at least 1 and the parent one level coarser");
    if (level.factor < 2 || parent.factor != level.factor || parent.worldwidth != level.worldwidth ||
        parent.worldheight != level.worldheight || level.worldwidth < 4)
        throw std::invalid_argument(prefix + "the level and its parent must share the world grid and the factor");
    if (level.domaincolumn0 % level.factor != 0 || level.domaincolumns % level.factor != 0 ||
        level.domainrow0 % level.factor != 0 || level.domainrows % level.factor != 0 ||
        level.columns % level.factor != 0 || level.rows % level.factor != 0 ||
        level.keptcolumnoffset() % level.factor != 0 || level.keptrowoffset() % level.factor != 0)
        throw std::invalid_argument(prefix + "the level window must be whole parent cells");
    if (fieldsize != level.domainsize() || parentsize != parent.domainsize())
        throw std::invalid_argument(prefix + "the fields must have their windows' domain sizes");
    // Every parent cell of the level's domain inside the parent's kept box: the domain's
    // parent columns are contiguous, so its first and last parent column and row decide.
    const int firstcolumn = parent.localcolumn(level.parentcolumn(0));
    const int lastcolumn = parent.localcolumn(level.parentcolumn(level.domaincolumns - 1));
    const std::int64_t firstrow = parent.localrow(level.parentrow(0));
    const std::int64_t lastrow = parent.localrow(level.parentrow(level.domainrows - 1));
    const int keptcolumn = parent.keptcolumnoffset(), keptrow = parent.keptrowoffset();
    if (firstcolumn < keptcolumn || firstcolumn + level.domainparentcolumns() - 1 != lastcolumn ||
        lastcolumn >= keptcolumn + parent.columns || firstrow < keptrow || lastrow >= keptrow + parent.rows)
        throw std::invalid_argument(prefix + "every parent cell of the level's domain must lie in the parent's kept box");
}

// v_c of the domain's parent cell (p, q).
double parentvalue(const LevelWindow& level, const LevelWindow& parent, const std::vector<float>& parentfield, int p, int q)
{
    const int x = parent.localcolumn(level.parentcolumn(p * level.factor));
    const int y = static_cast<int>(parent.localrow(level.parentrow(q * level.factor)));
    return static_cast<double>(parentfield[parent.index(x, y)]);
}

// The unweighted mean of the children of the domain's parent cell (p, q), in double, the
// children in row-major order.
double childmean(const LevelWindow& level, const std::vector<float>& field, int p, int q)
{
    const int f = level.factor;
    double sum = 0.0;
    for (int y = q * f; y < (q + 1) * f; y++)
        for (int x = p * f; x < (p + 1) * f; x++)
            sum += static_cast<double>(field[level.index(x, y)]);
    return sum / (static_cast<double>(f) * f);
}

// The kept box's parent cells in the domain's parent indices: p in first .. first + count - 1.
struct KeptParents
{
    int column0 = 0, columns = 0, row0 = 0, rows = 0;
};

KeptParents keptparents(const LevelWindow& level)
{
    KeptParents kept;
    kept.column0 = level.keptcolumnoffset() / level.factor;
    kept.columns = level.columns / level.factor;
    kept.row0 = level.keptrowoffset() / level.factor;
    kept.rows = level.rows / level.factor;
    return kept;
}

// A headroom in metres stored as a float no larger than the double value (so a child moved
// by its whole headroom never passes its bound after the store), at least 0.
float headroom(double metres)
{
    if (!(metres > 0.0))
        return 0.0f;
    float value = static_cast<float>(metres);
    if (static_cast<double>(value) > metres)
        value = std::nextafter(value, 0.0f);
    return value;
}

// The headroom of a child at `value` whose stored result must stay at or above `floor`
// (down) or at or below `ceiling` (up): the difference as a float, reduced while the stored
// float(value -+ h) would pass the bound (the difference rounds when the bound is far smaller
// in magnitude than the value, e.g. the next float above a sea level of 0); 0 if a few steps
// do not suffice.
float roomdown(double value, float floor)
{
    float room = headroom(value - static_cast<double>(floor));
    for (int step = 0; step < 4 && room > 0.0f; step++)
    {
        if (static_cast<float>(value - static_cast<double>(room)) >= floor)
            return room;
        room = std::nextafter(room, 0.0f);
    }
    return room > 0.0f && static_cast<float>(value - static_cast<double>(room)) >= floor ? room : 0.0f;
}

float roomup(double value, float ceiling)
{
    float room = headroom(static_cast<double>(ceiling) - value);
    for (int step = 0; step < 4 && room > 0.0f; step++)
    {
        if (static_cast<float>(value + static_cast<double>(room)) <= ceiling)
            return room;
        room = std::nextafter(room, 0.0f);
    }
    return room > 0.0f && static_cast<float>(value + static_cast<double>(room)) <= ceiling ? room : 0.0f;
}

// The bounds' sizes and values, as projectlevel and smoothparentmeans require them.
void checkbounds(const LevelWindow& level, const ProjectionBounds& bounds, const char* name)
{
    const std::size_t size = level.domainsize();
    if (bounds.rate.size() != size || bounds.up_m.size() != size || bounds.down_m.size() != size ||
        (!bounds.bed.empty() && bounds.bed.size() != size))
        throw std::invalid_argument(std::string(name) + ": the bounds must have the level's domain size");
    for (std::size_t cell = 0; cell < size; cell++)
        if (!std::isfinite(bounds.rate[cell]) || bounds.rate[cell] < 0.0f || std::isnan(bounds.up_m[cell]) ||
            std::isnan(bounds.down_m[cell]) || bounds.up_m[cell] < 0.0f || bounds.down_m[cell] < 0.0f)
            throw std::invalid_argument(std::string(name) +
                                        ": a rate is negative or not finite, or a headroom negative or NaN");
}

double largestmeanerror(const LevelWindow& level, const std::vector<float>& field, const LevelWindow& parent,
                        const std::vector<float>& parentfield)
{
    const KeptParents kept = keptparents(level);
    double largest = 0.0;
    for (int q = kept.row0; q < kept.row0 + kept.rows; q++)
        for (int p = kept.column0; p < kept.column0 + kept.columns; p++)
            largest = std::max(largest, std::fabs(childmean(level, field, p, q) - parentvalue(level, parent, parentfield, p, q)));
    return largest;
}

// The largest |lambda_c| and |lambda_c - lambda_c'| (kept parents adjacent east-west or
// north-south) over the kept box's parents, lambda_c = v_c - mean_C(c)(field).
void lambdastatistics(const LevelWindow& level, const std::vector<float>& field, const LevelWindow& parent,
                      const std::vector<float>& parentfield, double& largest, double& largeststep)
{
    const KeptParents kept = keptparents(level);
    std::vector<double> lambda(static_cast<std::size_t>(kept.columns) * static_cast<std::size_t>(kept.rows));
    auto lambdaat = [&](int i, int j) -> double& {
        return lambda[static_cast<std::size_t>(j) * static_cast<std::size_t>(kept.columns) + static_cast<std::size_t>(i)];
    };
    for (int j = 0; j < kept.rows; j++)
        for (int i = 0; i < kept.columns; i++)
            lambdaat(i, j) = parentvalue(level, parent, parentfield, kept.column0 + i, kept.row0 + j) -
                             childmean(level, field, kept.column0 + i, kept.row0 + j);
    largest = 0.0;
    largeststep = 0.0;
    for (int j = 0; j < kept.rows; j++)
        for (int i = 0; i < kept.columns; i++)
        {
            largest = std::max(largest, std::fabs(lambdaat(i, j)));
            if (i + 1 < kept.columns)
                largeststep = std::max(largeststep, std::fabs(lambdaat(i, j) - lambdaat(i + 1, j)));
            if (j + 1 < kept.rows)
                largeststep = std::max(largeststep, std::fabs(lambdaat(i, j) - lambdaat(i, j + 1)));
        }
}
}

ProjectionBounds projectionbounds(const LevelWindow& window, const std::vector<std::uint8_t>& channel,
                                  const std::vector<std::uint8_t>& sea, const std::vector<float>& lakesurface_m,
                                  const std::vector<float>& elevation_m, float sealevel, double channelrate,
                                  const std::vector<float>& ownlake_m, const std::vector<float>& bed_m)
{
    const std::size_t size = window.domainsize();
    if ((!channel.empty() && channel.size() != size) || (!sea.empty() && sea.size() != size) ||
        (!lakesurface_m.empty() && lakesurface_m.size() != size) || elevation_m.size() != size ||
        (!ownlake_m.empty() && ownlake_m.size() != size) || (!bed_m.empty() && bed_m.size() != size))
        throw std::invalid_argument("projectionbounds: every field must cover the window's domain");
    if (!std::isfinite(channelrate) || channelrate < 0.0)
        throw std::invalid_argument("projectionbounds: the channel rate must be finite and not negative");
    const float infinity = std::numeric_limits<float>::infinity();
    const float seabelow = std::nextafter(sealevel, -infinity);
    const float seaabove = std::nextafter(sealevel, infinity);
    auto wet = [&](std::size_t cell) {
        return !lakesurface_m.empty() && lakesurface_m[cell] > 0.0f && (sea.empty() || sea[cell] == 0);
    };
    ProjectionBounds bounds;
    bounds.rate.assign(size, 1.0f);
    bounds.up_m.assign(size, infinity);
    bounds.down_m.assign(size, infinity);
    if (!bed_m.empty())
        bounds.bed.assign(size, 0);
    for (int y = 0; y < window.domainrows; y++)
        for (int x = 0; x < window.domaincolumns; x++)
        {
            const std::size_t cell = window.index(x, y);
            const double value = static_cast<double>(elevation_m[cell]);
            if (!sea.empty() && sea[cell] != 0)
                bounds.up_m[cell] = roomup(value, seabelow);
            else if (wet(cell))
                bounds.up_m[cell] = roomup(value, std::nextafter(lakesurface_m[cell], -infinity));
            else
            {
                // A land channel cell: the channel rate on its bed, or (S3.2, with the beds given) rate 1
                // on its cell mean with its bed as a floor.
                const bool channelcell = !channel.empty() && channel[cell] != 0;
                if (channelcell && bed_m.empty())
                    bounds.rate[cell] = static_cast<float>(channelrate);
                if (channelcell && !bed_m.empty())
                    bounds.bed[cell] = 1;
                // The floor: above the sea level and above the water beside the cell.
                float level = -infinity;
                if (!ownlake_m.empty() && ownlake_m[cell] > 0.0f)
                    level = ownlake_m[cell];
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        const int nx = x + dx, ny = y + dy;
                        if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= window.domaincolumns || ny >= window.domainrows)
                            continue;
                        if (wet(window.index(nx, ny)))
                            level = std::max(level, lakesurface_m[window.index(nx, ny)]);
                    }
                float floor = level > -infinity ? std::max(seaabove, std::nextafter(level, infinity)) : seaabove;
                if (channelcell && !bed_m.empty())
                    floor = std::max(floor, bed_m[cell]);
                bounds.down_m[cell] = roomdown(value, floor);
            }
        }
    return bounds;
}

ProjectionReport projectlevel(const LevelWindow& level, std::vector<float>& field, const ProjectionBounds& bounds,
                              const LevelWindow& parent, const std::vector<float>& parentfield)
{
    checkpair(level, field.size(), parent, parentfield.size(), "projectlevel");
    const std::size_t size = level.domainsize();
    checkbounds(level, bounds, "projectlevel");

    const int f = level.factor;
    const int count = f * f;
    const double children = static_cast<double>(f) * f;
    const int pcolumns = level.domainparentcolumns(), prows = level.domainparentrows();
    const KeptParents kept = keptparents(level);
    auto inkept = [&](int p, int q) {
        return p >= kept.column0 && p < kept.column0 + kept.columns && q >= kept.row0 && q < kept.row0 + kept.rows;
    };

    ProjectionReport report;
    report.parents = kept.columns * kept.rows;
    std::vector<double> lambda(static_cast<std::size_t>(pcolumns) * static_cast<std::size_t>(prows), 0.0);
    // The shift of every domain cell and whether it moved (rated, or shifted by a fallback).
    std::vector<double> correction(size, 0.0);
    std::vector<std::uint8_t> moved(size, 0);
    std::vector<std::size_t> cells(static_cast<std::size_t>(count));
    std::vector<double> shift(static_cast<std::size_t>(count));
    std::vector<int> active, next;
    active.reserve(static_cast<std::size_t>(count));
    next.reserve(static_cast<std::size_t>(count));
    for (int q = 0; q < prows; q++)
        for (int p = 0; p < pcolumns; p++)
        {
            // The children in row-major absolute order.
            int n = 0;
            for (int y = q * f; y < (q + 1) * f; y++)
                for (int x = p * f; x < (p + 1) * f; x++)
                    cells[static_cast<std::size_t>(n++)] = level.index(x, y);
            const double lambdac = parentvalue(level, parent, parentfield, p, q) - childmean(level, field, p, q);
            std::fill(shift.begin(), shift.end(), 0.0);
            if (lambdac != 0.0)
            {
                const bool up = lambdac > 0.0;
                auto room = [&](int i) {
                    const std::size_t cell = cells[static_cast<std::size_t>(i)];
                    return static_cast<double>(up ? bounds.up_m[cell] : bounds.down_m[cell]);
                };
                auto rate = [&](int i) { return static_cast<double>(bounds.rate[cells[static_cast<std::size_t>(i)]]); };
                // The beds of the fallback: the unrated children, or those bounds.bed marks.
                auto isbed = [&](int i) {
                    return bounds.bed.empty() ? !(rate(i) > 0.0) : bounds.bed[cells[static_cast<std::size_t>(i)]] != 0;
                };
                double remainder = lambdac * children;
                active.clear();
                for (int i = 0; i < count; i++)
                    if (rate(i) > 0.0 && room(i) > 0.0)
                        active.push_back(i);
                // The common shift, capped: each round the children whose headroom is below
                // their share take it and leave; a round that caps none ends it.
                while (!active.empty())
                {
                    double rates = 0.0;
                    for (const int i : active)
                        rates += rate(i);
                    const double delta = remainder / rates;
                    next.clear();
                    for (const int i : active)
                    {
                        const double want = rate(i) * delta;
                        if (std::fabs(want) >= room(i))
                        {
                            shift[static_cast<std::size_t>(i)] = up ? room(i) : -room(i);
                            remainder -= shift[static_cast<std::size_t>(i)];
                        }
                        else
                            next.push_back(i);
                    }
                    if (next.size() == active.size())
                    {
                        for (const int i : active)
                            shift[static_cast<std::size_t>(i)] = rate(i) * delta;
                        remainder = 0.0;
                        active.clear();
                    }
                    else
                        active.swap(next);
                }
                // The headrooms exhausted: the rest on the beds (the unrated children), or
                // with none on every child, so the mean is exact.
                if (remainder != 0.0)
                {
                    const bool counted = inkept(p, q) && std::fabs(remainder) > 1.0e-6 * children;
                    int beds = 0;
                    for (int i = 0; i < count; i++)
                        beds += isbed(i) ? 1 : 0;
                    if (beds > 0)
                    {
                        for (int i = 0; i < count; i++)
                            if (isbed(i))
                                shift[static_cast<std::size_t>(i)] += remainder / beds;
                        report.bedsmovedparents += counted ? 1 : 0;
                    }
                    else
                    {
                        for (int i = 0; i < count; i++)
                            shift[static_cast<std::size_t>(i)] += remainder / children;
                        if (counted)
                        {
                            report.forcedparents++;
                            if (std::fabs(remainder / children) > report.maxforcedmove_m)
                            {
                                report.maxforcedmove_m = std::fabs(remainder / children);
                                report.maxforcedcolumn = wraplevel(level.domaincolumn0 / f + p, level.levelwidth() / f);
                                report.maxforcedrow = level.domainrow0 / f + q;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < count; i++)
            {
                const std::size_t cell = cells[static_cast<std::size_t>(i)];
                const double d = shift[static_cast<std::size_t>(i)];
                field[cell] = static_cast<float>(static_cast<double>(field[cell]) + d);
                correction[cell] = d;
                moved[cell] = bounds.rate[cell] > 0.0f || d != 0.0 ? 1 : 0;
            }
            lambda[static_cast<std::size_t>(q) * static_cast<std::size_t>(pcolumns) + static_cast<std::size_t>(p)] = lambdac;
        }

    // The statistics over the kept box's parents.
    auto lambdaat = [&](int p, int q) {
        return lambda[static_cast<std::size_t>(q) * static_cast<std::size_t>(pcolumns) + static_cast<std::size_t>(p)];
    };
    for (int q = kept.row0; q < kept.row0 + kept.rows; q++)
        for (int p = kept.column0; p < kept.column0 + kept.columns; p++)
        {
            report.maxlambda_m = std::max(report.maxlambda_m, std::fabs(lambdaat(p, q)));
            if (p + 1 < kept.column0 + kept.columns)
                report.maxlambdastep_m = std::max(report.maxlambdastep_m, std::fabs(lambdaat(p, q) - lambdaat(p + 1, q)));
            if (q + 1 < kept.row0 + kept.rows)
                report.maxlambdastep_m = std::max(report.maxlambdastep_m, std::fabs(lambdaat(p, q) - lambdaat(p, q + 1)));
        }
    report.maxparentmeanerror_m = largestmeanerror(level, field, parent, parentfield);

    // The border step: edge-adjacent kept cells in different parents, both moved. Inside the
    // kept box such a pair straddles a parent edge: x + 1 a multiple of factor (east-west)
    // or y + 1 (north-south), in the domain's local indices.
    const int x0 = level.keptcolumnoffset(), y0 = level.keptrowoffset();
    const int x1 = x0 + level.columns, y1 = y0 + level.rows;
    std::vector<double> steps;
    double largeststep = -1.0;
    auto addpair = [&](int x, int y, std::size_t a, std::size_t b) {
        if (moved[a] == 0 || moved[b] == 0)
            return;
        const double step = std::fabs(correction[a] - correction[b]);
        steps.push_back(step);
        if (step > largeststep)
        {
            largeststep = step;
            report.maxborderstepcolumn = wraplevel(level.domaincolumn0 / f + x / f, level.levelwidth() / f);
            report.maxbordersteprow = level.domainrow0 / f + y / f;
        }
    };
    for (int y = y0; y < y1; y++)
        for (int x = x0 + f - 1; x + 1 < x1; x += f)
            addpair(x, y, level.index(x, y), level.index(x + 1, y));
    for (int y = y0 + f - 1; y + 1 < y1; y += f)
        for (int x = x0; x < x1; x++)
            addpair(x, y, level.index(x, y), level.index(x, y + 1));
    if (!steps.empty())
    {
        const std::size_t n = steps.size();
        const std::size_t rank = (99 * n + 99) / 100 - 1;  // ceil(0.99 n) - 1
        std::nth_element(steps.begin(), steps.begin() + static_cast<std::ptrdiff_t>(rank), steps.end());
        report.borderstep99_m = steps[rank];
        report.maxborderstep_m = *std::max_element(steps.begin() + static_cast<std::ptrdiff_t>(rank), steps.end());
    }
    return report;
}

SmoothingReport smoothparentmeans(const LevelWindow& level, std::vector<float>& field, const ProjectionBounds& bounds,
                                  const LevelWindow& parent, const std::vector<float>& parentfield, int passes,
                                  double relaxation, int basis)
{
    checkpair(level, field.size(), parent, parentfield.size(), "smoothparentmeans");
    checkbounds(level, bounds, "smoothparentmeans");
    if (passes < 0 || !std::isfinite(relaxation) || !(relaxation > 0.0))
        throw std::invalid_argument("smoothparentmeans: passes must not be negative and the relaxation positive and finite");
    if (basis != 0 && basis != 1)
        throw std::invalid_argument("smoothparentmeans: the basis must be 0 (Catmull-Rom) or 1 (the cubic B-spline)");
    // The interpolation of mu between p1 and p2 at the fraction t: Catmull-Rom (C1, through the
    // values) or the uniform cubic B-spline (C2, approximating; S3.2's fallback).
    const auto spline = [basis](double p0, double p1, double p2, double p3, double t) {
        if (basis == 0)
            return catmullrom(p0, p1, p2, p3, t);
        const double u = 1.0 - t, t2 = t * t, t3 = t2 * t;
        return (u * u * u * p0 + (3.0 * t3 - 6.0 * t2 + 4.0) * p1 + (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) * p2 + t3 * p3) /
               6.0;
    };
    SmoothingReport report;
    report.passes = passes;
    lambdastatistics(level, field, parent, parentfield, report.maxlambdabefore_m, report.maxlambdastepbefore_m);
    if (passes == 0)
    {
        report.maxlambdaafter_m = report.maxlambdabefore_m;
        report.maxlambdastepafter_m = report.maxlambdastepbefore_m;
        return report;
    }

    const int f = level.factor;
    const double children = static_cast<double>(f) * f;
    const std::size_t size = level.domainsize();
    const int pcolumns = level.domainparentcolumns(), prows = level.domainparentrows();
    const bool wholeturn = static_cast<std::int64_t>(level.domaincolumns) >= level.levelwidth();
    // The interpolation position of a child at offset o inside its parent (correctparentmeans's):
    // t = (o + 0.5) / f - 0.5 parent cells from the parent's centre, from the parent itself
    // (lower 0) or its predecessor (lower -1, t < 0).
    std::vector<int> lower(static_cast<std::size_t>(f));
    std::vector<double> fraction(static_cast<std::size_t>(f));
    for (int o = 0; o < f; o++)
    {
        const double t = (o + 0.5) / f - 0.5;
        lower[static_cast<std::size_t>(o)] = t < 0.0 ? -1 : 0;
        fraction[static_cast<std::size_t>(o)] = t < 0.0 ? t + 1.0 : t;
    }
    const auto column = [&](int c) {
        return wholeturn ? ((c % pcolumns) + pcolumns) % pcolumns : std::clamp(c, 0, pcolumns - 1);
    };

    std::vector<double> origin(size), total(size, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
        origin[cell] = static_cast<double>(field[cell]);
    std::vector<double> mu(static_cast<std::size_t>(pcolumns) * static_cast<std::size_t>(prows));
    for (int pass = 0; pass < passes; pass++)
    {
        // mu_c = lambda_c f^2 / the rates of the children with room in lambda_c's direction.
        for (int q = 0; q < prows; q++)
            for (int p = 0; p < pcolumns; p++)
            {
                const double lambdac = parentvalue(level, parent, parentfield, p, q) - childmean(level, field, p, q);
                double rates = 0.0;
                if (lambdac != 0.0)
                    for (int y = q * f; y < (q + 1) * f; y++)
                        for (int x = p * f; x < (p + 1) * f; x++)
                        {
                            const std::size_t cell = level.index(x, y);
                            const double room = lambdac > 0.0 ? static_cast<double>(bounds.up_m[cell]) - total[cell]
                                                              : static_cast<double>(bounds.down_m[cell]) + total[cell];
                            if (bounds.rate[cell] > 0.0f && room > 0.0)
                                rates += static_cast<double>(bounds.rate[cell]);
                        }
                mu[static_cast<std::size_t>(q) * static_cast<std::size_t>(pcolumns) + static_cast<std::size_t>(p)] =
                    rates > 0.0 ? lambdac * children / rates : 0.0;
            }
        // Lambda at the child centres; each rated child moves rate x Lambda, its total capped.
        for (int y = 0; y < level.domainrows; y++)
        {
            const std::size_t oy = static_cast<std::size_t>(y % f);
            const int rowbase = y / f + lower[oy];
            const double ty = fraction[oy];
            for (int x = 0; x < level.domaincolumns; x++)
            {
                const std::size_t cell = level.index(x, y);
                const double rate = static_cast<double>(bounds.rate[cell]);
                if (!(rate > 0.0))
                    continue;
                const std::size_t ox = static_cast<std::size_t>(x % f);
                const int columnbase = x / f + lower[ox];
                const double tx = fraction[ox];
                double interpolated[4];
                for (int b = 0; b < 4; b++)
                {
                    const std::size_t base = static_cast<std::size_t>(std::clamp(rowbase - 1 + b, 0, prows - 1)) *
                                             static_cast<std::size_t>(pcolumns);
                    const auto value = [&](int a) { return mu[base + static_cast<std::size_t>(column(columnbase - 1 + a))]; };
                    interpolated[b] = spline(value(0), value(1), value(2), value(3), tx);
                }
                const double move =
                    relaxation * rate * spline(interpolated[0], interpolated[1], interpolated[2], interpolated[3], ty);
                const double moved = std::clamp(total[cell] + move, -static_cast<double>(bounds.down_m[cell]),
                                                static_cast<double>(bounds.up_m[cell]));
                if (moved == total[cell])
                    continue;
                total[cell] = moved;
                field[cell] = static_cast<float>(origin[cell] + moved);
            }
        }
    }
    lambdastatistics(level, field, parent, parentfield, report.maxlambdaafter_m, report.maxlambdastepafter_m);
    const int x0 = level.keptcolumnoffset(), y0 = level.keptrowoffset();
    for (int y = y0; y < y0 + level.rows; y++)
        for (int x = x0; x < x0 + level.columns; x++)
            report.maxmove_m = std::max(report.maxmove_m, std::fabs(total[level.index(x, y)]));
    return report;
}

double levelparentmeanerror(const LevelWindow& level, const std::vector<float>& field, const LevelWindow& parent,
                            const std::vector<float>& parentfield)
{
    checkpair(level, field.size(), parent, parentfield.size(), "levelparentmeanerror");
    return largestmeanerror(level, field, parent, parentfield);
}

ResampledField resampleto(const std::vector<LevelWindow>& windows, const std::vector<std::vector<float>>& fields,
                          double radius_m, double metresperpixel)
{
    if (!std::isfinite(metresperpixel) || metresperpixel <= 0.0 || !std::isfinite(radius_m) || radius_m <= 0.0)
        throw std::invalid_argument("resampleto: metres per pixel and the radius must be positive and finite");
    if (windows.empty())
        throw std::invalid_argument("resampleto: the chain is empty");
    if (fields.size() != windows.size())
        throw std::invalid_argument("resampleto: one field per level of the chain");
    const LevelWindow& target = windows.back();
    const int finest = target.level;
    for (std::size_t l = 0; l < windows.size(); l++)
    {
        const LevelWindow& window = windows[l];
        if (window.level != static_cast<int>(l) || fields[l].size() != window.domainsize() ||
            window.worldwidth != target.worldwidth || window.worldheight != target.worldheight ||
            window.factor != target.factor || (finest >= 1 && window.factor < 2))
            throw std::invalid_argument("resampleto: the chain is malformed");
    }

    // The source level: the coarsest whose cell is at most the pixel.
    int source = -1;
    double cell_m = 0.0;
    for (int l = 0; l <= finest; l++)
    {
        const double size = pi * radius_m / static_cast<double>(windows[static_cast<std::size_t>(l)].levelheight());
        if (size <= metresperpixel * (1.0 + 1.0e-9))
        {
            source = l;
            cell_m = size;
            break;
        }
    }
    if (source < 0)
        throw std::invalid_argument("resampleto: no level of the chain is as fine as the requested pixel");

    ResampledField result;
    ResampleReport& report = result.report;
    report.sourcelevel = source;
    report.metresperpixel = metresperpixel;
    report.sourcecell_m = cell_m;
    report.pixelcells = metresperpixel / cell_m;
    const double n = report.pixelcells;
    const double whole = std::round(n);
    const bool iswhole = std::fabs(n - whole) <= 1.0e-9 * n;
    const LevelWindow& window = windows[static_cast<std::size_t>(source)];

    // The target's kept box at the source level.
    const std::int64_t down = levelscale(target.factor, finest - source);
    if (target.column0 % down != 0 || target.columns % down != 0 || target.row0 % down != 0 || target.rows % down != 0)
        throw std::invalid_argument("resampleto: the target's kept box is not whole cells of the source level");
    const std::int64_t columns = target.columns / down, rows = target.rows / down;
    const int xstart = window.localcolumn(target.column0 / down);
    const std::int64_t ystart = window.localrow(target.row0 / down);
    if (xstart < window.keptcolumnoffset() || xstart + columns > window.keptcolumnoffset() + window.columns ||
        ystart < window.keptrowoffset() || ystart + rows > window.keptrowoffset() + window.rows)
        throw std::invalid_argument("resampleto: the source level's kept box does not cover the target's");

    // A whole n takes the exact whole value for the pixel geometry.
    const double edge = iswhole ? whole : n;
    result.width = static_cast<int>(std::floor(static_cast<double>(columns) / edge));
    result.height = static_cast<int>(std::floor(static_cast<double>(rows) / edge));
    if (result.width < 1 || result.height < 1)
        throw std::invalid_argument("resampleto: the output would have no pixel");
    // At least one pixel, so a whole n is at most the target's extent and fits an integer.
    report.aligned = iswhole && window.scale() % static_cast<std::int64_t>(whole) == 0;
    result.values.assign(static_cast<std::size_t>(result.width) * static_cast<std::size_t>(result.height), 0.0f);
    const std::vector<float>& field = fields[static_cast<std::size_t>(source)];
    auto value = [&](std::int64_t i, std::int64_t j) {
        return static_cast<double>(field[window.index(xstart + static_cast<int>(i), static_cast<int>(ystart + j))]);
    };

    for (int q = 0; q < result.height; q++)
        for (int p = 0; p < result.width; p++)
        {
            double mean = 0.0;
            if (report.aligned)
            {
                const std::int64_t m = static_cast<std::int64_t>(whole);
                double sum = 0.0;
                for (std::int64_t j = q * m; j < (q + 1) * m; j++)
                    for (std::int64_t i = p * m; i < (p + 1) * m; i++)
                        sum += value(i, j);
                mean = sum / (static_cast<double>(m) * static_cast<double>(m));
            }
            else
            {
                const double west = p * edge, east = (p + 1) * edge, north = q * edge, south = (q + 1) * edge;
                const std::int64_t ifirst = static_cast<std::int64_t>(std::floor(west));
                const std::int64_t ilast = std::min<std::int64_t>(columns - 1, static_cast<std::int64_t>(std::ceil(east)) - 1);
                const std::int64_t jfirst = static_cast<std::int64_t>(std::floor(north));
                const std::int64_t jlast = std::min<std::int64_t>(rows - 1, static_cast<std::int64_t>(std::ceil(south)) - 1);
                double sum = 0.0;
                for (std::int64_t j = jfirst; j <= jlast; j++)
                {
                    const double oy = std::min(static_cast<double>(j + 1), south) - std::max(static_cast<double>(j), north);
                    if (oy <= 0.0)
                        continue;
                    for (std::int64_t i = ifirst; i <= ilast; i++)
                    {
                        const double ox = std::min(static_cast<double>(i + 1), east) - std::max(static_cast<double>(i), west);
                        if (ox > 0.0)
                            sum += value(i, j) * ox * oy;
                    }
                }
                mean = sum / (edge * edge);
            }
            result.values[static_cast<std::size_t>(q) * static_cast<std::size_t>(result.width) + static_cast<std::size_t>(p)] =
                static_cast<float>(mean);
        }
    return result;
}

void smoothinglevel(const LevelContext& context, LevelState& state)
{
    // The refinement loop's smoothing before the projection (S2.7), moved here from
    // refinement_levels.cpp unchanged in S3.0: the projection by bounds (sea and wet cells may
    // sink but not rise past their water, dry land may rise but not sink to the sea level or the
    // water beside it, its own global lake cell's surface included, land channel beds take the
    // channel rate), the smoothing with the same bounds, then the bounds rebuilt on the smoothed
    // field (the same floors and ceilings).
    if (context.refinement == nullptr)
        return;
    const RegionalGlobalInputs& inputs = context.inputs;
    const RefinementParameters& refinement = *context.refinement;
    const LevelWindow& level = *context.levelwindow;
    const std::size_t count = level.domainsize();
    std::vector<float>& elevation = state.elevation_m;
    const std::vector<float>& lakesurface = state.water.lakesurface_m;
    std::vector<float>& ownlake = state.ownlake_m;
    ownlake.assign(count, 0.0f);
    for (int y = 0; y < level.domainrows; y++)
        for (int x = 0; x < level.domaincolumns; x++)
        {
            const std::size_t global = static_cast<std::size_t>(level.globalrow(y)) * static_cast<std::size_t>(inputs.width) +
                                       static_cast<std::size_t>(level.globalcolumn(x));
            if (inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                ownlake[level.index(x, y)] = inputs.lakesurface_m[global];
        }
    // S3.2: with the channel cells at their cell means, the beds are their floors (rate 1).
    static const std::vector<float> nobeds;
    const std::vector<float>& beds = refinement.channelcellslope > 0.0 ? state.bed_m : nobeds;
    state.bounds = projectionbounds(level, state.channels.channel, state.sea, lakesurface, elevation, inputs.sealevel,
                                    refinement.projectionweightchannel, ownlake, beds);
    state.lift_m = elevation;  // the elevation entering the smoothing (lift_m's reference)
    const SmoothingReport smoothing = smoothparentmeans(level, elevation, state.bounds, *context.projectionparent,
                                                        *context.projectionfield, refinement.smoothpasses,
                                                        refinement.smoothrelaxation, refinement.smoothbasis);
    if (state.diagnostics != nullptr)
    {
        state.diagnostics->smoothing = smoothing;
        // The 99th percentile (nearest rank) of the kept cells' |move| by the smoothing (S3.2's check).
        std::vector<double> moves;
        moves.reserve(level.keptsize());
        const int x0 = level.keptcolumnoffset(), y0 = level.keptrowoffset();
        for (int y = y0; y < y0 + level.rows; y++)
            for (int x = x0; x < x0 + level.columns; x++)
            {
                const std::size_t cell = level.index(x, y);
                moves.push_back(std::fabs(static_cast<double>(elevation[cell]) - static_cast<double>(state.lift_m[cell])));
            }
        if (!moves.empty())
        {
            const std::size_t rank = static_cast<std::size_t>(std::ceil(0.99 * static_cast<double>(moves.size()))) - 1;
            std::nth_element(moves.begin(), moves.begin() + static_cast<std::ptrdiff_t>(rank), moves.end());
            state.diagnostics->smoothmove99_m = moves[rank];
        }
    }
    if (refinement.smoothpasses > 0)
    {
        recordlevelstep(state, LevelStep::smoothing);
        state.bounds = projectionbounds(level, state.channels.channel, state.sea, lakesurface, elevation, inputs.sealevel,
                                        refinement.projectionweightchannel, ownlake, beds);
    }
}

void projectionlevel(const LevelContext& context, LevelState& state)
{
    // The refinement loop's exact projection and the drying after it (S2.3), moved here from
    // refinement_levels.cpp unchanged in S3.0; then lift_m.
    if (context.refinement == nullptr)
        return;
    const LevelWindow& level = *context.levelwindow;
    const std::size_t count = level.domainsize();
    std::vector<float>& elevation = state.elevation_m;
    std::vector<float>& lakesurface = state.water.lakesurface_m;
    std::vector<uint8_t> wet(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        wet[cell] = lakesurface[cell] > 0.0f ? 1 : 0;
    const ProjectionReport projection =
        projectlevel(level, elevation, state.bounds, *context.projectionparent, *context.projectionfield);
    if (state.diagnostics != nullptr)
        state.diagnostics->projection = projection;
    recordlevelstep(state, LevelStep::projection);

    // After the projection: wet non-channel cells at or above their surface dried (counted in
    // the kept box).
    const int keptx0 = level.keptcolumnoffset(), kepty0 = level.keptrowoffset();
    const auto kept = [&](int x, int y) {
        return x >= keptx0 && x < keptx0 + level.columns && y >= kepty0 && y < kepty0 + level.rows;
    };
    for (int y = 0; y < level.domainrows; y++)
        for (int x = 0; x < level.domaincolumns; x++)
        {
            const std::size_t cell = level.index(x, y);
            if (wet[cell] == 0 || state.channels.channel[cell] != 0 || elevation[cell] < lakesurface[cell])
                continue;
            lakesurface[cell] = 0.0f;
            if (kept(x, y) && state.diagnostics != nullptr)
                state.diagnostics->driedcells++;
        }

    // The level's lift: the final elevation less the elevation entering the smoothing (S3.2
    // prolongs it to the next level).
    for (std::size_t cell = 0; cell < count; cell++)
        state.lift_m[cell] =
            static_cast<float>(static_cast<double>(elevation[cell]) - static_cast<double>(state.lift_m[cell]));
}
}
