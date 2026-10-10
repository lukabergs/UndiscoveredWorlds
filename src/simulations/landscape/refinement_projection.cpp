#include "simulations/landscape/refinement_projection.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

// The projection, its error measure and the resample (step S2.2,
// docs/architecture/multiscale-plan.md; contracts in the header, written in S0.4; the
// bounds that replaced S2.2's weights are S2.3's rework, plan section 3, "S2.3's projection
// weights and detail levels"). A
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
}

ProjectionBounds projectionbounds(const LevelWindow& window, const std::vector<std::uint8_t>& channel,
                                  const std::vector<std::uint8_t>& sea, const std::vector<float>& lakesurface_m,
                                  const std::vector<float>& elevation_m, float sealevel, double channelrate,
                                  const std::vector<float>& ownlake_m)
{
    const std::size_t size = window.domainsize();
    if ((!channel.empty() && channel.size() != size) || (!sea.empty() && sea.size() != size) ||
        (!lakesurface_m.empty() && lakesurface_m.size() != size) || elevation_m.size() != size ||
        (!ownlake_m.empty() && ownlake_m.size() != size))
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
                if (!channel.empty() && channel[cell] != 0)
                    bounds.rate[cell] = static_cast<float>(channelrate);
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
                const float floor = level > -infinity ? std::max(seaabove, std::nextafter(level, infinity)) : seaabove;
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
    if (bounds.rate.size() != size || bounds.up_m.size() != size || bounds.down_m.size() != size)
        throw std::invalid_argument("projectlevel: the bounds must have the level's domain size");
    for (std::size_t cell = 0; cell < size; cell++)
        if (!std::isfinite(bounds.rate[cell]) || bounds.rate[cell] < 0.0f || std::isnan(bounds.up_m[cell]) ||
            std::isnan(bounds.down_m[cell]) || bounds.up_m[cell] < 0.0f || bounds.down_m[cell] < 0.0f)
            throw std::invalid_argument("projectlevel: a rate is negative or not finite, or a headroom negative or NaN");

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
                        beds += rate(i) > 0.0 ? 0 : 1;
                    if (beds > 0)
                    {
                        for (int i = 0; i < count; i++)
                            if (!(rate(i) > 0.0))
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
                            report.maxforcedmove_m = std::max(report.maxforcedmove_m, std::fabs(remainder / children));
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
    auto addpair = [&](std::size_t a, std::size_t b) {
        if (moved[a] != 0 && moved[b] != 0)
            steps.push_back(std::fabs(correction[a] - correction[b]));
    };
    for (int y = y0; y < y1; y++)
        for (int x = x0 + f - 1; x + 1 < x1; x += f)
            addpair(level.index(x, y), level.index(x + 1, y));
    for (int y = y0 + f - 1; y + 1 < y1; y += f)
        for (int x = x0; x < x1; x++)
            addpair(level.index(x, y), level.index(x, y + 1));
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
}
