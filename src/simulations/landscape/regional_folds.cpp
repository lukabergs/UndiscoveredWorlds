#include "simulations/landscape/regional_folds.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <queue>
#include <tuple>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/upsampling.hpp"

// Multiscale plan S5.2, S5.3b and S5.3c (regional_folds.hpp): the fold fabric, A sin of one phase field,
// the least-squares potential of the across-strike wavenumber on the global grid, evaluated at a fine
// cell by bicubic Hermite interpolation.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Salt of the orientation flood's hashed ties (deterministiccombine(seed, salt, global cell)).
constexpr std::int64_t foldsalt = 0x464f4c44;  // "FOLD"
// The solve, numerical settings, not a look: successive over-relaxed Gauss-Seidel sweeps (factor solverelaxation)
// until the largest change of a sweep is below solvetolerance_rad, at most solvemaxsweeps sweeps.
constexpr double solverelaxation = 1.97;
constexpr double solvetolerance_rad = 1.0e-2;
constexpr int solvemaxsweeps = 20000;
// The rings of global cells around the mask that take an extrapolated phase: every 2 x 2 corner of a fine
// cell whose bicubic strain stencil (4 x 4 global cells) holds a mask cell lies within two cells of it.
constexpr int extensionrings = 2;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapglobal(int column, int width) { return ((column % width) + width) % width; }

double smoothstep(double value, double lower, double upper)
{
    if (!(upper > lower))
        return value >= lower ? 1.0 : 0.0;
    const double t = std::clamp((value - lower) / (upper - lower), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

double strainweight(const FoldParameters& folds, double strain)
{
    return smoothstep(strain, folds.strainmin, 2.0 * folds.strainmin);
}

// The nearest-rank percentile of `values` (sorted in place; 0 when empty).
double nearestrank(std::vector<double>& values, double fraction)
{
    if (values.empty())
        return 0.0;
    const std::size_t rank = static_cast<std::size_t>(std::ceil(fraction * static_cast<double>(values.size())));
    const std::size_t index = std::min(values.size() - 1, rank > 0 ? rank - 1 : 0);
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
    return values[index];
}

// The orientation flood (S5.3c item 2): a Dijkstra search over the mask's 8-neighbours by ground distance from
// the cells already oriented; a cell takes the sign that keeps its axis within 90 degrees of the oriented axis of
// the cell it is reached from (on equal distances the predecessor of the smaller hash). Fills `orientation` of
// every mask cell the sources reach; returns how many it oriented.
int orientationflood(const FoldPhaseField& field, const std::vector<uint8_t>& mask, const std::vector<double>& axiseast,
                     const std::vector<double>& axissouth, std::uint64_t seed, const std::vector<std::size_t>& sources,
                     std::vector<int8_t>& orientation, std::vector<uint8_t>& done)
{
    const int width = field.width, height = field.height;
    const std::size_t size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    using Entry = std::tuple<double, std::uint64_t, std::size_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    std::vector<double> distance(size, std::numeric_limits<double>::infinity());
    std::vector<std::int64_t> from(size, -1);
    const auto hash = [seed](std::size_t cell) { return deterministiccombine(seed, foldsalt, static_cast<std::int64_t>(cell)); };
    for (const std::size_t cell : sources)
    {
        distance[cell] = 0.0;
        queue.emplace(0.0, hash(cell), cell);
    }
    int oriented = 0;
    while (!queue.empty())
    {
        const auto [d, unused, cell] = queue.top();
        (void)unused;
        queue.pop();
        if (done[cell] != 0 || d > distance[cell])
            continue;
        done[cell] = 1;
        if (orientation[cell] == 0)
        {
            const std::size_t previous = static_cast<std::size_t>(from[cell]);
            const double dot = axiseast[cell] * axiseast[previous] + axissouth[cell] * axissouth[previous];
            orientation[cell] = static_cast<int8_t>(dot * static_cast<double>(orientation[previous]) >= 0.0 ? 1 : -1);
            oriented++;
        }
        const int x = static_cast<int>(cell % static_cast<std::size_t>(width));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(width));
        for (int dy = -1; dy <= 1; dy++)
        {
            const int ny = y + dy;
            if (ny < 0 || ny >= height)
                continue;
            for (int dx = -1; dx <= 1; dx++)
            {
                if (dx == 0 && dy == 0)
                    continue;
                const std::size_t next = at(wrapglobal(x + dx, width), ny, width);
                if (mask[next] == 0 || done[next] != 0)
                    continue;
                const double east = dx == 0 ? 0.0 : 0.5 * (field.rowdx_km[static_cast<std::size_t>(y)] + field.rowdx_km[static_cast<std::size_t>(ny)]);
                const double south = dy == 0 ? 0.0 : field.dy_km;
                const double candidate = d + std::sqrt(east * east + south * south);
                if (candidate < distance[next] ||
                    (candidate == distance[next] && from[next] >= 0 && hash(cell) < hash(static_cast<std::size_t>(from[next]))))
                {
                    const bool improved = candidate < distance[next];
                    distance[next] = candidate;
                    from[next] = static_cast<std::int64_t>(cell);
                    if (improved)
                        queue.emplace(candidate, hash(next), next);
                }
            }
        }
    }
    return oriented;
}

// One corner's Hermite data: the value, the slopes along the two cell parameters and the cross term.
struct Corner
{
    double value = 0.0, slopes = 0.0, slopet = 0.0, cross = 0.0;
};

Corner corner(const FoldPhaseField& field, int column, int row)
{
    const std::size_t cell = at(column, row, field.width);
    Corner c;
    c.value = field.phase_rad[cell];
    c.slopes = field.wavenumbereast[cell] * field.rowdx_km[static_cast<std::size_t>(row)];
    c.slopet = field.wavenumbersouth[cell] * field.dy_km;
    c.cross = field.cross_rad[cell];
    return c;
}

// The cubic Hermite basis on [0, 1]: h00, h01 (the values at 0 and 1), h10, h11 (the slopes at 0 and 1).
void hermite(double s, double& h00, double& h01, double& h10, double& h11)
{
    const double s2 = s * s, s3 = s2 * s;
    h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    h01 = -2.0 * s3 + 3.0 * s2;
    h10 = s3 - 2.0 * s2 + s;
    h11 = s3 - s2;
}

// The basis' derivatives on [0, 1].
void hermitederivative(double s, double& h00, double& h01, double& h10, double& h11)
{
    const double s2 = s * s;
    h00 = 6.0 * s2 - 6.0 * s;
    h01 = -6.0 * s2 + 6.0 * s;
    h10 = 3.0 * s2 - 4.0 * s + 1.0;
    h11 = 3.0 * s2 - 2.0 * s;
}
}

bool foldfabricon(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    return parameters.folds.amplitudemax_m > 0.0 && size > 0 && inputs.foldstrain.size() == size;
}

double foldwavelength_km(const FoldParameters& folds, double sediment_m, bool withsediment)
{
    const double lower = std::min(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double upper = std::max(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double cover_km = withsediment ? std::max(0.0, sediment_m) / 1000.0 : 0.0;
    const double raw = folds.wavelengthbase_km + folds.wavelengthpercover_km * cover_km;
    return std::clamp(std::isfinite(raw) ? raw : 0.0, lower, upper);
}

double foldamplitude_m(const FoldParameters& folds, double strain, double wavelength_km, double heightabovesea_m)
{
    if (!(folds.amplitudemax_m > 0.0) || !(strain > 0.0) || !(wavelength_km > 0.0))
        return 0.0;
    const double structural = folds.relieffraction * wavelength_km * 1000.0 / pi * std::sqrt(strain);
    const double elevation = folds.reliefelevationfraction * std::max(0.0, heightabovesea_m);
    const double amplitude = std::min(std::min(structural, folds.amplitudemax_m), elevation);
    return std::isfinite(amplitude) ? std::max(0.0, amplitude) * strainweight(folds, strain) : 0.0;
}

FoldPhaseField foldphasefield(const RegionalGlobalInputs& inputs, const FoldParameters& folds)
{
    const auto start = std::chrono::steady_clock::now();
    FoldPhaseField field;
    const int width = inputs.width, height = inputs.height;
    const std::size_t size = static_cast<std::size_t>(std::max(0, width)) * static_cast<std::size_t>(std::max(0, height));
    if (size == 0 || inputs.foldstrain.size() != size || !(inputs.radius_m > 0.0))
        return field;
    field.width = width;
    field.height = height;
    field.phase_rad.assign(size, 0.0);
    field.wavenumbereast.assign(size, 0.0);
    field.wavenumbersouth.assign(size, 0.0);
    field.cross_rad.assign(size, 0.0);
    field.state.assign(size, 0);
    field.orientation.assign(size, 0);

    // The sphere's metric: the north-south cell, the east-west cell of each row at its centre and at the edge
    // between rows j and j + 1 (km).
    field.dy_km = pi * inputs.radius_m / static_cast<double>(height) / 1000.0;
    const double column_km = 2.0 * pi * inputs.radius_m / static_cast<double>(width) / 1000.0;
    field.column_km = column_km;
    field.rowdx_km.resize(static_cast<std::size_t>(height));
    std::vector<double> edgedx_km(static_cast<std::size_t>(height), 0.0);
    for (int j = 0; j < height; j++)
    {
        const double centre = 0.5 * pi - pi * (static_cast<double>(j) + 0.5) / static_cast<double>(height);
        const double edge = 0.5 * pi - pi * (static_cast<double>(j) + 1.0) / static_cast<double>(height);
        field.rowdx_km[static_cast<std::size_t>(j)] = column_km * std::cos(centre);
        edgedx_km[static_cast<std::size_t>(j)] = column_km * std::max(0.0, std::cos(edge));
    }
    const auto rowdx = [&field](int j) { return field.rowdx_km[static_cast<std::size_t>(j)]; };

    // The mask (w(eps) > 0), the axis (east, south) and the wavenumber's size 2 pi / L per global cell.
    const bool withsediment = inputs.sediment_m.size() == size;
    const bool withdirection = inputs.folddirection_rad.size() == size;
    std::vector<uint8_t> mask(size, 0);
    std::vector<double> axiseast(size, 1.0), axissouth(size, 0.0), wavenumber(size, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const double strain = std::max(0.0, static_cast<double>(inputs.foldstrain[cell]));
        mask[cell] = strainweight(folds, strain) > 0.0 ? 1 : 0;
        field.maskcells += mask[cell];
        const double theta = withdirection && std::isfinite(inputs.folddirection_rad[cell]) ? static_cast<double>(inputs.folddirection_rad[cell]) : 0.0;
        axiseast[cell] = std::cos(theta);
        axissouth[cell] = std::sin(theta);
        const double wavelength = foldwavelength_km(folds, withsediment ? static_cast<double>(inputs.sediment_m[cell]) : 0.0, withsediment);
        wavenumber[cell] = wavelength > 0.0 ? 2.0 * pi / wavelength : 0.0;
    }
    if (field.maskcells == 0)
        return field;

    // The axis smoothed (round two): the doubled-angle vectors w(eps) (cos 2 theta, sin 2 theta) of the mask cells
    // within a (2 r + 1) x (2 r + 1) block, r = ceil(axissmoothing_cells), summed with the Gaussian weight
    // exp(-d^2 / (2 sigma^2)) of the ground distance d in global north-south cells (columns by the mean of the two
    // rows' east-west cells), then the angle halved; a cell whose sum vanishes keeps its own axis. sigma <= 0: none.
    if (folds.axissmoothing_cells > 0.0)
    {
        const double sigma = folds.axissmoothing_cells;
        const int reach = static_cast<int>(std::ceil(sigma));
        std::vector<double> doubledcos(size, 0.0), doubledsin(size, 0.0);
        for (std::size_t cell = 0; cell < size; cell++)
            if (mask[cell] != 0)
            {
                const double weight = strainweight(folds, std::max(0.0, static_cast<double>(inputs.foldstrain[cell])));
                doubledcos[cell] = weight * (axiseast[cell] * axiseast[cell] - axissouth[cell] * axissouth[cell]);
                doubledsin[cell] = weight * 2.0 * axiseast[cell] * axissouth[cell];
            }
        std::vector<double> smoothedeast = axiseast, smoothedsouth = axissouth;
        for (int j = 0; j < height; j++)
            for (int i = 0; i < width; i++)
            {
                const std::size_t cell = at(i, j, width);
                if (mask[cell] == 0)
                    continue;
                double c = 0.0, s = 0.0;
                for (int dj = -reach; dj <= reach; dj++)
                {
                    const int nj = j + dj;
                    if (nj < 0 || nj >= height)
                        continue;
                    const double eastcell = 0.5 * (rowdx(j) + rowdx(nj)) / field.dy_km;
                    for (int di = -reach; di <= reach; di++)
                    {
                        const std::size_t next = at(wrapglobal(i + di, width), nj, width);
                        if (mask[next] == 0)
                            continue;
                        const double east = static_cast<double>(di) * eastcell, south = static_cast<double>(dj);
                        const double gauss = std::exp(-(east * east + south * south) / (2.0 * sigma * sigma));
                        c += gauss * doubledcos[next];
                        s += gauss * doubledsin[next];
                    }
                }
                if (!(std::hypot(c, s) > 0.0))
                    continue;
                const double theta = 0.5 * std::atan2(s, c);
                smoothedeast[cell] = std::cos(theta);
                smoothedsouth[cell] = std::sin(theta);
            }
        axiseast = std::move(smoothedeast);
        axissouth = std::move(smoothedsouth);
    }

    // The orientation (item 2, round two: by continuity alone): each component's first cell (row-major) points north
    // (east when its axis is east-west), the flood orients the rest from it.
    std::vector<uint8_t> done(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (mask[cell] == 0 || done[cell] != 0)
            continue;
        field.orientation[cell] = static_cast<int8_t>(axissouth[cell] < 0.0 || (axissouth[cell] == 0.0 && axiseast[cell] >= 0.0) ? 1 : -1);
        field.orientedbynorth++;
        field.orientedbyflood += orientationflood(field, mask, axiseast, axissouth, inputs.seed, {cell}, field.orientation, done);
    }
    for (std::size_t cell = 0; cell < size; cell++)
        if (mask[cell] != 0)
        {
            const double sign = static_cast<double>(field.orientation[cell]);
            field.wavenumbereast[cell] = sign * wavenumber[cell] * axiseast[cell];
            field.wavenumbersouth[cell] = sign * wavenumber[cell] * axissouth[cell];
            field.state[cell] = 1;
        }

    // The least-squares potential (item 1): the edges between 4-neighbouring mask cells (periodic east-west, none
    // across the poles), each asking Phi_n - Phi_c = the mean wavenumber of its two cells times the step on the
    // ground, weighted by the finite-volume Laplacian's metric (dy / dx across a row, dx at the edge / dy between
    // rows); natural (Neumann) conditions at the mask's edge.
    struct Link
    {
        std::size_t cell = 0;
        double weight = 0.0, target = 0.0;
    };
    std::vector<std::size_t> cells;
    std::vector<std::size_t> first;
    std::vector<Link> links;
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            if (mask[cell] == 0)
                continue;
            cells.push_back(cell);
            first.push_back(links.size());
            const auto link = [&](int ni, int nj, double weight, double step_km, bool eastward) {
                const std::size_t next = at(wrapglobal(ni, width), nj, width);
                if (mask[next] == 0)
                    return;
                const double k = eastward ? 0.5 * (field.wavenumbereast[cell] + field.wavenumbereast[next])
                                          : 0.5 * (field.wavenumbersouth[cell] + field.wavenumbersouth[next]);
                links.push_back({next, weight, k * step_km});
            };
            if (j > 0)
                link(i, j - 1, edgedx_km[static_cast<std::size_t>(j - 1)] / field.dy_km, -field.dy_km, false);
            if (width > 1)
            {
                link(i - 1, j, field.dy_km / rowdx(j), -rowdx(j), true);
                link(i + 1, j, field.dy_km / rowdx(j), rowdx(j), true);
            }
            if (j + 1 < height)
                link(i, j + 1, edgedx_km[static_cast<std::size_t>(j)] / field.dy_km, field.dy_km, false);
        }
    first.push_back(links.size());

    // The start: the wavenumber integrated along a breadth-first tree from each component's first cell; the
    // components numbered in that order.
    std::vector<int> component(size, -1);
    std::vector<std::size_t> slot(size, 0);
    for (std::size_t index = 0; index < cells.size(); index++)
        slot[cells[index]] = index;
    for (const std::size_t root : cells)
    {
        if (component[root] >= 0)
            continue;
        const int number = field.components++;
        component[root] = number;
        field.phase_rad[root] = 0.0;
        std::deque<std::size_t> queue{root};
        while (!queue.empty())
        {
            const std::size_t cell = queue.front();
            queue.pop_front();
            const std::size_t index = slot[cell];
            for (std::size_t l = first[index]; l < first[index + 1]; l++)
            {
                const Link& link = links[l];
                if (component[link.cell] >= 0)
                    continue;
                component[link.cell] = number;
                field.phase_rad[link.cell] = field.phase_rad[cell] + link.target;
                queue.push_back(link.cell);
            }
        }
    }

    // Over-relaxed Gauss-Seidel in double, the cells in row-major order.
    for (field.sweeps = 0; field.sweeps < solvemaxsweeps;)
    {
        double largest = 0.0;
        for (std::size_t index = 0; index < cells.size(); index++)
        {
            double numerator = 0.0, denominator = 0.0;
            for (std::size_t l = first[index]; l < first[index + 1]; l++)
            {
                numerator += links[l].weight * (field.phase_rad[links[l].cell] - links[l].target);
                denominator += links[l].weight;
            }
            if (!(denominator > 0.0))
                continue;
            double& value = field.phase_rad[cells[index]];
            const double change = solverelaxation * (numerator / denominator - value);
            largest = std::max(largest, std::fabs(change));
            value += change;
        }
        field.sweeps++;
        field.lastchange_rad = largest;
        if (largest < solvetolerance_rad)
            break;
    }

    // The gauge: each component's mean phase 0.
    {
        std::vector<double> sums(static_cast<std::size_t>(field.components), 0.0);
        std::vector<double> counts(static_cast<std::size_t>(field.components), 0.0);
        for (const std::size_t cell : cells)
        {
            sums[static_cast<std::size_t>(component[cell])] += field.phase_rad[cell];
            counts[static_cast<std::size_t>(component[cell])] += 1.0;
        }
        for (const std::size_t cell : cells)
            field.phase_rad[cell] -= sums[static_cast<std::size_t>(component[cell])] / counts[static_cast<std::size_t>(component[cell])];
    }

    // The misfit: the weighted rms of the edges' Phi_n - Phi_c - target relative to their targets' (each edge once).
    {
        double misfit = 0.0, scale = 0.0;
        for (std::size_t index = 0; index < cells.size(); index++)
            for (std::size_t l = first[index]; l < first[index + 1]; l++)
            {
                const Link& link = links[l];
                if (link.cell < cells[index])
                    continue;
                if (field.wavenumbereast[link.cell] * field.wavenumbereast[cells[index]] +
                        field.wavenumbersouth[link.cell] * field.wavenumbersouth[cells[index]] <
                    0.0)
                    field.reversaledges++;
                const double error = field.phase_rad[link.cell] - field.phase_rad[cells[index]] - link.target;
                misfit += link.weight * error * error;
                scale += link.weight * link.target * link.target;
            }
        field.misfit = scale > 0.0 ? std::sqrt(misfit / scale) : 0.0;
    }

    // The extension rings: a cell next to the defined ones (8-neighbours) takes the mean of their phases carried
    // over by their wavenumbers, and the mean of their wavenumbers; ring by ring, each from the cells defined before.
    for (int ring = 0; ring < extensionrings; ring++)
    {
        const std::vector<uint8_t> defined = field.state;
        for (int j = 0; j < height; j++)
            for (int i = 0; i < width; i++)
            {
                const std::size_t cell = at(i, j, width);
                if (defined[cell] != 0)
                    continue;
                double phase = 0.0, east = 0.0, south = 0.0;
                int count = 0;
                for (int dj = -1; dj <= 1; dj++)
                {
                    const int nj = j + dj;
                    if (nj < 0 || nj >= height)
                        continue;
                    for (int di = -1; di <= 1; di++)
                    {
                        const std::size_t next = at(wrapglobal(i + di, width), nj, width);
                        if ((di == 0 && dj == 0) || defined[next] == 0)
                            continue;
                        const double stepeast = -static_cast<double>(di) * 0.5 * (rowdx(j) + rowdx(nj));
                        const double stepsouth = -static_cast<double>(dj) * field.dy_km;
                        phase += field.phase_rad[next] + field.wavenumbereast[next] * stepeast + field.wavenumbersouth[next] * stepsouth;
                        east += field.wavenumbereast[next];
                        south += field.wavenumbersouth[next];
                        count++;
                    }
                }
                if (count == 0)
                    continue;
                field.phase_rad[cell] = phase / count;
                field.wavenumbereast[cell] = east / count;
                field.wavenumbersouth[cell] = south / count;
                field.state[cell] = 2;
            }
    }

    // The Hermite cross term per defined cell (rad per cell squared): the mean of d/dt (k_east dx) and d/ds (k_south
    // dy) by central differences, one-sided beside an undefined cell, 0 without a defined neighbour.
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            if (field.state[cell] == 0)
                continue;
            const auto difference = [&](bool before, bool after, double previous, double here, double next) {
                if (before && after)
                    return 0.5 * (next - previous);
                if (after)
                    return next - here;
                if (before)
                    return here - previous;
                return 0.0;
            };
            const bool north = j > 0 && field.state[at(i, j - 1, width)] != 0;
            const bool south = j + 1 < height && field.state[at(i, j + 1, width)] != 0;
            const double alongt = difference(north, south, north ? field.wavenumbereast[at(i, j - 1, width)] * rowdx(j - 1) : 0.0,
                                             field.wavenumbereast[cell] * rowdx(j),
                                             south ? field.wavenumbereast[at(i, j + 1, width)] * rowdx(j + 1) : 0.0);
            const std::size_t westcell = at(wrapglobal(i - 1, width), j, width), eastcell = at(wrapglobal(i + 1, width), j, width);
            const bool west = field.state[westcell] != 0, east = field.state[eastcell] != 0;
            const double alongs = difference(west, east, field.wavenumbersouth[westcell] * field.dy_km, field.wavenumbersouth[cell] * field.dy_km,
                                             field.wavenumbersouth[eastcell] * field.dy_km);
            field.cross_rad[cell] = 0.5 * (alongt + alongs);
        }
    field.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return field;
}

bool foldphase_rad(const FoldPhaseField& field, double u, double v, double& phase, double* misfit)
{
    const int width = field.width, height = field.height;
    if (width <= 0 || height <= 0 || !std::isfinite(u) || !std::isfinite(v))
        return false;
    const double x = u - 0.5, y = v - 0.5;
    const double fx = std::floor(x), fy = std::floor(y);
    const double s = x - fx, t = y - fy;
    const int i0 = wrapglobal(static_cast<int>(static_cast<std::int64_t>(fx) % width), width), i1 = wrapglobal(i0 + 1, width);
    const int j0 = static_cast<int>(fy), j1 = j0 + 1;
    if (j1 < 0 || j0 >= height)
        return false;
    double hs0 = 0.0, hs1 = 0.0, gs0 = 0.0, gs1 = 0.0;
    hermite(s, hs0, hs1, gs0, gs1);
    double ds0 = 0.0, ds1 = 0.0, es0 = 0.0, es1 = 0.0;
    hermitederivative(s, ds0, ds1, es0, es1);
    double alongs = 0.0, alongt = 0.0, keast = 0.0, ksouth = 0.0;
    if (j0 < 0 || j1 >= height)
    {
        // The polar half rows: the edge row's Hermite curve carried north or south by its slope along t.
        const int row = j0 < 0 ? j1 : j0;
        if (field.state[at(i0, row, width)] == 0 || field.state[at(i1, row, width)] == 0)
            return false;
        const Corner a = corner(field, i0, row), b = corner(field, i1, row);
        const double offset = j0 < 0 ? t - 1.0 : t;
        phase = hs0 * a.value + hs1 * b.value + gs0 * a.slopes + gs1 * b.slopes +
                offset * (hs0 * a.slopet + hs1 * b.slopet + gs0 * a.cross + gs1 * b.cross);
        if (misfit == nullptr)
            return true;
        alongs = ds0 * a.value + ds1 * b.value + es0 * a.slopes + es1 * b.slopes +
                 offset * (ds0 * a.slopet + ds1 * b.slopet + es0 * a.cross + es1 * b.cross);
        alongt = hs0 * a.slopet + hs1 * b.slopet + gs0 * a.cross + gs1 * b.cross;
        const std::size_t ca = at(i0, row, width), cb = at(i1, row, width);
        keast = (1.0 - s) * field.wavenumbereast[ca] + s * field.wavenumbereast[cb];
        ksouth = (1.0 - s) * field.wavenumbersouth[ca] + s * field.wavenumbersouth[cb];
    }
    else
    {
        if (field.state[at(i0, j0, width)] == 0 || field.state[at(i1, j0, width)] == 0 || field.state[at(i0, j1, width)] == 0 ||
            field.state[at(i1, j1, width)] == 0)
            return false;
        double ht0 = 0.0, ht1 = 0.0, gt0 = 0.0, gt1 = 0.0;
        hermite(t, ht0, ht1, gt0, gt1);
        const Corner c00 = corner(field, i0, j0), c10 = corner(field, i1, j0), c01 = corner(field, i0, j1), c11 = corner(field, i1, j1);
        const auto term = [](const Corner& c, double hs, double gs, double ht, double gt) {
            return hs * ht * c.value + gs * ht * c.slopes + hs * gt * c.slopet + gs * gt * c.cross;
        };
        phase = term(c00, hs0, gs0, ht0, gt0) + term(c10, hs1, gs1, ht0, gt0) + term(c01, hs0, gs0, ht1, gt1) + term(c11, hs1, gs1, ht1, gt1);
        if (misfit == nullptr)
            return true;
        double dt0 = 0.0, dt1 = 0.0, et0 = 0.0, et1 = 0.0;
        hermitederivative(t, dt0, dt1, et0, et1);
        alongs = term(c00, ds0, es0, ht0, gt0) + term(c10, ds1, es1, ht0, gt0) + term(c01, ds0, es0, ht1, gt1) + term(c11, ds1, es1, ht1, gt1);
        alongt = term(c00, hs0, gs0, dt0, et0) + term(c10, hs1, gs1, dt0, et0) + term(c01, hs0, gs0, dt1, et1) + term(c11, hs1, gs1, dt1, et1);
        const std::size_t a00 = at(i0, j0, width), a10 = at(i1, j0, width), a01 = at(i0, j1, width), a11 = at(i1, j1, width);
        keast = (1.0 - t) * ((1.0 - s) * field.wavenumbereast[a00] + s * field.wavenumbereast[a10]) +
                t * ((1.0 - s) * field.wavenumbereast[a01] + s * field.wavenumbereast[a11]);
        ksouth = (1.0 - t) * ((1.0 - s) * field.wavenumbersouth[a00] + s * field.wavenumbersouth[a10]) +
                 t * ((1.0 - s) * field.wavenumbersouth[a01] + s * field.wavenumbersouth[a11]);
    }
    // The gradient on the ground: d/ds over the column's east-west width at v, d/dt over the north-south cell.
    const double eastwidth = field.column_km * std::cos(0.5 * pi - pi * v / static_cast<double>(height));
    const double gradeast = eastwidth > 0.0 ? alongs / eastwidth : keast, gradsouth = alongt / field.dy_km;
    const double size = std::hypot(keast, ksouth);
    *misfit = size > 0.0 ? std::hypot(gradeast - keast, gradsouth - ksouth) / size : 0.0;
    return true;
}

void addfoldfabric(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                   const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                   int parentscale, std::vector<uint8_t>* floorclamped)
{
    if (!foldfabricon(inputs, parameters))
        return;
    const FoldParameters& folds = parameters.folds;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || sea.size() != count || elevation_m.size() != count)
        return;
    if (floorclamped != nullptr)
        floorclamped->assign(count, 0);
    const int globalwidth = inputs.width, globalheight = inputs.height;
    const std::size_t globalsize = static_cast<std::size_t>(globalwidth) * static_cast<std::size_t>(globalheight);
    const int k = window.scale;
    const std::shared_ptr<const mesh::CellGrid> grid = regionalwindowgrid(window, inputs.radius_m);
    const double cell_km = northsouthspacing_m(*grid) / 1000.0;
    // The global north-south cell on the ground and the octave band of the level (addlevelnoise's).
    const double globalcell_km = pi * inputs.radius_m / static_cast<double>(globalheight) / 1000.0;
    if (!(cell_km > 0.0) || !(globalcell_km > 0.0))
        return;
    // The level that draws a cell (round two): the coarsest whose cell is at most L / 4.
    const double parentcell_km = parentscale > 0 ? globalcell_km / static_cast<double>(parentscale) : 0.0;

    // The phase field on the global grid (S5.3c), a pure function of the global fields.
    const FoldPhaseField phasefield = foldphasefield(inputs, folds);
    if (phasefield.maskcells == 0)
        return;

    // The global fields: the wavelength and the height above the sea level per global cell.
    const bool withsediment = inputs.sediment_m.size() == globalsize;
    std::vector<float> wavelengthfield(globalsize), heightfield(globalsize);
    const double sealevel = static_cast<double>(inputs.sealevel);
    for (std::size_t cell = 0; cell < globalsize; cell++)
    {
        wavelengthfield[cell] = static_cast<float>(
            foldwavelength_km(folds, withsediment ? static_cast<double>(inputs.sediment_m[cell]) : 0.0, withsediment));
        heightfield[cell] = static_cast<float>(std::max(0.0, static_cast<double>(inputs.elevation_m[cell]) - sealevel));
    }

    // The cell fields, interpolated as the noise amplitude's relief is.
    const std::vector<float> strain = upsample(inputs.foldstrain, globalwidth, globalheight, window);
    const std::vector<float> cellwavelength = upsample(wavelengthfield, globalwidth, globalheight, window);
    const std::vector<float> above = upsample(heightfield, globalwidth, globalheight, window);
    const bool withice = inputs.glaciation.size() == globalsize;
    const std::vector<float> glaciation = withice ? upsample(inputs.glaciation, globalwidth, globalheight, window)
                                                  : std::vector<float>{};
    const bool withlakes = inputs.lakesurface_m.size() == globalsize && inputs.network.sea.size() == globalsize;
    const double wavelengthlower = std::min(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const double wavelengthupper = std::max(folds.wavelengthmin_km, folds.wavelengthmax_km);
    const float floorlevel = std::max(std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity()),
                                      static_cast<float>(sealevel + std::max(0.0, parameters.coast.shorefloor_m)));

    // Every cell reads the elevation before the fabric only at itself, so the order does not matter.
    for (int y = 0; y < height; y++)
    {
        const double ratio = eastwestratio(*grid, at(0, y, width));  // the row's first cell
        const double v = (static_cast<double>(window.absolutefinerow(y)) + 0.5) / static_cast<double>(k);
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = at(x, y, width);
            if (sea[cell] != 0 || channels.channel[cell] != 0)
                continue;
            // Dry land at or below the shore floor (a filled sea pocket, a basin below the sea level)
            // takes no fabric: it could only rise, and the projection must then sink the parent's
            // other children (the first judging run forced parents there).
            if (!(elevation_m[cell] > floorlevel))
                continue;
            // The gates (S5.3b): under the global ice cover, and in a global lake cell's block.
            if (withice && std::clamp(static_cast<double>(glaciation[cell]), 0.0, 1.0) > folds.icemax)
                continue;
            if (withlakes)
            {
                const std::size_t global = window.globalcell(x, y);
                if (inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                    continue;
            }
            const double epsilon = std::max(0.0, static_cast<double>(strain[cell]));
            if (!(strainweight(folds, epsilon) > 0.0))
                continue;
            const double wavelengthcell = std::clamp(static_cast<double>(cellwavelength[cell]), wavelengthlower, wavelengthupper);
            if (!(cell_km <= 0.25 * wavelengthcell) || (parentscale > 0 && !(parentcell_km > 0.25 * wavelengthcell)))
                continue;  // drawn at another level of the refinement (or at none: fewer than 4 cells per wavelength)
            const double amplitude = foldamplitude_m(folds, epsilon, wavelengthcell, static_cast<double>(above[cell]));
            if (!(amplitude > 0.0))
                continue;

            // The water gap: the largest falloff over the channel cells within gapfraction x L.
            double gap = 1.0;
            const double gapcells = folds.gapfraction * wavelengthcell / cell_km;
            if (gapcells > 0.0)
            {
                const int rows = std::min(k, static_cast<int>(std::floor(gapcells)));
                const int columns = columnextent(gapcells, ratio, k);
                double strongest = 0.0;
                for (int dy = -rows; dy <= rows; dy++)
                {
                    const int ny = y + dy;
                    if (ny < 0 || ny >= height)
                        continue;
                    for (int dx = -columns; dx <= columns; dx++)
                    {
                        const int nx = x + dx;
                        if (nx < 0 || nx >= width || channels.channel[at(nx, ny, width)] == 0)
                            continue;
                        const double east = static_cast<double>(dx) * ratio;
                        const double distance = std::sqrt(east * east + static_cast<double>(dy * dy));
                        if (distance >= gapcells)
                            continue;
                        strongest = std::max(strongest, 1.0 - smoothstep(distance / gapcells, 0.0, 1.0));
                    }
                }
                gap = 1.0 - strongest;
            }
            if (!(gap > 0.0))
                continue;

            // The phase at the cell's centre (S5.3c item 3).
            const double u = (static_cast<double>(window.absolutefinecolumn(x)) + 0.5) / static_cast<double>(k);
            double phase = 0.0, misfit = 0.0;
            if (!foldphase_rad(phasefield, u, v, phase, &misfit))
                continue;
            // Where the phase cannot follow k, the folds end (round two): the misfit fade.
            const double fade = 1.0 - smoothstep(misfit, folds.misfitfadestart, folds.misfitfadeend);
            const double fabric = amplitude * gap * fade * std::sin(phase);
            if (fabric == 0.0 || !std::isfinite(fabric))
                continue;
            const float before = elevation_m[cell];
            const float value = static_cast<float>(static_cast<double>(before) + fabric);
            // Never lowered below the shore floor above the sea level.
            if (value < floorlevel && floorclamped != nullptr)
                (*floorclamped)[cell] = 1;
            elevation_m[cell] = std::max(value, floorlevel);
        }
    }
}

void foldfabriclevel(const LevelContext& context, LevelState& state)
{
    if (!foldfabricon(context.inputs, context.parameters))
        return;
    const RegionalWindow& window = context.window;
    const int parentscale = context.parent != nullptr ? static_cast<int>(levelscale(context.factor, context.level - 1)) : 0;
    const double cap = channelcellcap_m(context);
    const bool record = state.diagnostics != nullptr;
    std::vector<float> before = state.elevation_m;
    std::vector<uint8_t> floorclamped;
    addfoldfabric(context.inputs, window, context.parameters, state.channels, state.sea, state.elevation_m, parentscale,
                  record ? &floorclamped : nullptr);
    if (record)
    {
        const int domainwidth = window.domainfinewidth();
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        std::vector<double> changes;
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (!floorclamped.empty() && floorclamped[cell] != 0)
                    state.diagnostics->foldfloorclamped++;
                if (before[cell] != state.elevation_m[cell])
                    changes.push_back(std::fabs(static_cast<double>(state.elevation_m[cell]) - static_cast<double>(before[cell])));
            }
        state.diagnostics->foldcells = static_cast<int>(changes.size());
        state.diagnostics->foldamplitude_p50_m = nearestrank(changes, 0.5);
        state.diagnostics->foldamplitude_p99_m = nearestrank(changes, 0.99);
    }
    channelcellmeans(window, state.channels, state.sea, cap, before, state.elevation_m);
    // S5.3b: the surfaces before and after the fabric, for the lakes step (unfoldforlakes).
    state.folded_m = state.elevation_m;
    state.unfolded_m = std::move(before);
    recordlevelstep(state, LevelStep::foldfabric);
}

void unfoldforlakes(const LevelContext& context, LevelState& state)
{
    (void)context;
    std::vector<float>& elevation = state.elevation_m;
    const std::size_t count = elevation.size();
    if (state.folded_m.size() != count || state.unfolded_m.size() != count || state.sea.size() != count ||
        state.channels.channel.size() != count)
    {
        state.folded_m.clear();
        state.unfolded_m.clear();
        return;
    }
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const float folded = state.folded_m[cell], prior = state.unfolded_m[cell];
        if (folded == prior)
            continue;
        if (state.sea[cell] != 0 || state.channels.channel[cell] != 0)
        {
            state.unfolded_m[cell] = folded;  // a channel cell's mean: left, and not refolded
            continue;
        }
        const float present = elevation[cell];
        const float unfolded = present == folded ? prior : std::min(present, prior);
        state.folded_m[cell] = present;
        state.unfolded_m[cell] = unfolded;
        elevation[cell] = unfolded;
    }
}

void refoldafterlakes(const LevelContext& context, LevelState& state)
{
    std::vector<float>& elevation = state.elevation_m;
    const std::size_t count = elevation.size();
    const int width = context.window.domainfinewidth(), height = context.window.domainfineheight();
    if (state.folded_m.size() == count && state.unfolded_m.size() == count &&
        count == static_cast<std::size_t>(width) * static_cast<std::size_t>(height))
    {
        const std::vector<float>& surface = state.water.lakesurface_m;
        const bool withwater = surface.size() == count;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                const std::size_t cell = at(x, y, width);
                const float folded = state.folded_m[cell], unfolded = state.unfolded_m[cell];
                if (folded == unfolded || elevation[cell] != unfolded)
                    continue;  // not folded, or changed by the lakes (a lake bed, a shore raise, a rim)
                float level = 0.0f;
                if (withwater)
                {
                    if (surface[cell] > 0.0f)
                        continue;  // under the water
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int nx = x + dx, ny = y + dy;
                            if ((dx != 0 || dy != 0) && nx >= 0 && nx < width && ny >= 0 && ny < height)
                                level = std::max(level, surface[at(nx, ny, width)]);
                        }
                }
                if (level > 0.0f && !(unfolded > level))
                    continue;  // already at or below the water beside it: no fabric on the shore
                elevation[cell] =
                    level > 0.0f && !(folded > level) ? std::nextafter(level, std::numeric_limits<float>::infinity()) : folded;
            }
    }
    state.folded_m.clear();
    state.folded_m.shrink_to_fit();
    state.unfolded_m.clear();
    state.unfolded_m.shrink_to_fit();
}
}
