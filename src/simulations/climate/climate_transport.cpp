#include "climate_transport.hpp"
#include "climate_grid.hpp"
#include "parallel_rows.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace climatehydrology
{
namespace
{
// Integrate a bounded linear reconstruction in cumulative cell-area space.
// Prefix integrals make cost independent of the number of crossed cells.
struct LineIntegral
{
    std::vector<double> edge, mass, value, slope;
    void prepare(const std::vector<double>& q, const std::vector<double>& volume, bool periodic)
    {
        const int n = static_cast<int>(q.size());
        value = q; edge.resize(n + 1); mass.resize(n + 1); slope.resize(n);
        edge[0] = mass[0] = 0.0;
        for (int i = 0; i < n; ++i)
        {
            edge[i + 1] = edge[i] + volume[i];
            mass[i + 1] = mass[i] + volume[i] * q[i];
            const int l = periodic ? (i + n - 1) % n : std::max(0, i - 1);
            const int r = periodic ? (i + 1) % n : std::min(n - 1, i + 1);
            const double a = q[i] - q[l], b = q[r] - q[i];
            const double centred = (q[r] - q[l]) * volume[i] /
                (0.5 * volume[l] + volume[i] + 0.5 * volume[r]);
            slope[i] = a * b <= 0.0 ? 0.0 : std::copysign(
                std::min({2.0 * std::abs(a), 2.0 * std::abs(b), std::abs(centred)}), a);
        }
    }
    double integral(double position, bool periodic) const
    {
        double whole = 0.0;
        if (periodic)
        {
            const double turns = std::floor(position / edge.back());
            whole = turns * mass.back();
            position -= turns * edge.back();
        }
        else
            position = std::clamp(position, 0.0, edge.back());
        if (position >= edge.back()) return whole + mass.back();
        const auto i = static_cast<std::size_t>(std::upper_bound(edge.begin(), edge.end(), position) - edge.begin() - 1);
        const double volume = edge[i + 1] - edge[i];
        const double f = (position - edge[i]) / volume;
        return whole + mass[i] + volume * (value[i] * f + 0.5 * slope[i] * (f * f - f));
    }
};

struct LineWorkspace
{
    std::vector<double> q, volume;
    LineIntegral integral;
};

struct TransportWorkspace
{
    climategrid::SphericalGrid grid;
    std::vector<double> volume, east0, east1, south0, south1;
    std::vector<double> q, east, south, fx, fy, gx, gy, qx, qy, vx, vy;
    void resize(std::size_t cells)
    {
        for (auto* v : {&volume, &east0, &east1, &south0, &south1, &q, &east,
            &south, &fx, &fy, &gx, &gy, &qx, &qy, &vx, &vy}) v->resize(cells);
    }
};

void fluxes(int columns, int rows, bool east, const std::vector<double>& q,
    const std::vector<double>& volume, const std::vector<double>& swept, std::vector<double>& flux, bool reuse)
{
    const int length = east ? columns : rows, lines = east ? rows : columns;
    flux.resize(q.size());
    parallelforrows(0, lines - 1, [&](int first, int last)
    {
    // Each executing thread owns its scratch, including concurrent callers of
    // the shared row executor. No tracer values are reused across calls.
    thread_local LineWorkspace retained;
    LineWorkspace fresh;
    auto& work = reuse ? retained : fresh;
    auto& lineQ = work.q; auto& lineVolume = work.volume; auto& integral = work.integral;
    lineQ.resize(length); lineVolume.resize(length);
    for (int line = first; line <= last; ++line)
    {
        const auto index = [&](int k) { return east ? line * columns + k : k * columns + line; };
        for (int k = 0; k < length; ++k)
        {
            lineQ[k] = q[index(k)]; lineVolume[k] = volume[index(k)];
        }
        integral.prepare(lineQ, lineVolume, east);
        for (int k = 0; k < length; ++k)
        {
            const int i = index(k);
            flux[i] = (!east && k == length - 1) ? 0.0 :
                integral.mass[k + 1] - integral.integral(integral.edge[k + 1] - swept[i], east);
        }
    }
    }, q.size() < 8192 ? lines : 16);
}
}

SphericalTracerTransportDiagnostics advectSphericalTracerFfsl(
    int columns, int rows, const std::vector<float>& source,
    const std::vector<float>& zonalWindMps, const std::vector<float>& meridionalWindMps,
    float timeStepSeconds, float planetRadiusMetres, const FfslOptions& options,
    std::vector<float>& destination)
{
    SphericalTracerTransportDiagnostics d;
    const std::size_t cells = static_cast<std::size_t>(std::max(0, columns)) * std::max(0, rows);
    if (!climategrid::validGlobalGridDimensions(columns, rows) || source.size() != cells ||
        zonalWindMps.size() != cells || meridionalWindMps.size() != cells ||
        (!options.endZonalWindMps.empty() && options.endZonalWindMps.size() != cells) ||
        (!options.endMeridionalWindMps.empty() && options.endMeridionalWindMps.size() != cells) ||
        timeStepSeconds <= 0.0f || planetRadiusMetres <= 0.0f)
    {
        destination = source;
        return d;
    }
    thread_local TransportWorkspace retained;
    TransportWorkspace fresh;
    auto& work = options.reuseWorkspace ? retained : fresh;
    auto& grid = work.grid;
    if (grid.columns != columns || grid.rows != rows || grid.radiusMetres != planetRadiusMetres)
        grid = climategrid::makeSphericalGrid(columns, rows, planetRadiusMetres);
    work.resize(cells);
    auto& volume = work.volume; auto& east0 = work.east0; auto& east1 = work.east1;
    auto& south0 = work.south0; auto& south1 = work.south1;
    const auto& endU = options.endZonalWindMps.empty() ? zonalWindMps : options.endZonalWindMps;
    const auto& endV = options.endMeridionalWindMps.empty() ? meridionalWindMps : options.endMeridionalWindMps;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const std::size_t i = static_cast<std::size_t>(y) * columns + x;
            const auto e = x + 1 < columns ? i + 1 : i + 1 - columns;
            if (!std::isfinite(source[i]) || source[i] < 0.0f ||
                !std::isfinite(zonalWindMps[i]) || !std::isfinite(meridionalWindMps[i]) ||
                !std::isfinite(endU[i]) || !std::isfinite(endV[i]))
                throw std::invalid_argument("FFSL needs finite winds and nonnegative column water");
            volume[i] = grid.cellAreasSquareMetres[y];
            const double eastScale = 0.5 * grid.zonalFaceLengthsMetres[y] * timeStepSeconds;
            east0[i] = eastScale * (static_cast<double>(zonalWindMps[i]) + zonalWindMps[e]);
            east1[i] = eastScale * (static_cast<double>(endU[i]) + endU[e]);
            if (y + 1 < rows)
            {
                const double southScale = 0.5 * grid.southFaceLengthsMetres[y] * timeStepSeconds;
                south0[i] = southScale * (static_cast<double>(meridionalWindMps[i]) + meridionalWindMps[i + columns]);
                south1[i] = southScale * (static_cast<double>(endV[i]) + endV[i + columns]);
            }
            else south0[i] = south1[i] = 0.0;
            d.initialAreaWeightedMass += source[i] * volume[i];
        }
    double deformation = 0.0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const std::size_t i = static_cast<std::size_t>(y) * columns + x;
            const auto w = x > 0 ? i - 1 : i + columns - 1;
            const double dx = std::max(std::abs(east0[i] - east0[w]), std::abs(east1[i] - east1[w]));
            const double dy = std::max(std::abs(south0[i] - (y ? south0[i - columns] : 0.0)),
                std::abs(south1[i] - (y ? south1[i - columns] : 0.0)));
            deformation = std::max(deformation, (dx + dy) / volume[i]);
            d.maximumZonalCourant = std::max(d.maximumZonalCourant,
                static_cast<float>(std::max(std::abs(east0[i]), std::abs(east1[i])) / volume[i]));
            d.maximumMeridionalCourant = std::max(d.maximumMeridionalCourant,
                static_cast<float>(std::max(std::abs(south0[i]), std::abs(south1[i])) /
                    std::min(volume[i], grid.cellAreasSquareMetres[std::min(y + 1, rows - 1)])));
        }
    d.substeps = std::max(1, static_cast<int>(std::ceil(deformation / std::clamp(options.maximumDeformation, 0.05, 0.8))));
    d.cellSubsteps = cells * d.substeps;
    d.maximumMultidimensionalCourant = d.maximumZonalCourant + d.maximumMeridionalCourant;
    d.eastIntegratedFlux.assign(cells, 0.0); d.southIntegratedFlux.assign(cells, 0.0);
    auto& q = work.q; q.assign(source.begin(), source.end());
    auto& east = work.east; auto& south = work.south;
    auto& fx = work.fx; auto& fy = work.fy; auto& gx = work.gx; auto& gy = work.gy;
    auto& qx = work.qx; auto& qy = work.qy; auto& vx = work.vx; auto& vy = work.vy;
    for (int step = 0; step < d.substeps; ++step)
    {
        const double t = (step + 0.5) / d.substeps;
        for (std::size_t i = 0; i < cells; ++i)
        {
            east[i] = ((1.0 - t) * east0[i] + t * east1[i]) / d.substeps;
            south[i] = ((1.0 - t) * south0[i] + t * south1[i]) / d.substeps;
        }
        // Bendall & Kent, SWIFT (2024), equations 35-40:
        // https://arxiv.org/html/2405.20006v3#S3.SS1
        // Transport unity as well as density; the outer sweeps use its updated
        // volumes. This removes the spurious compression of plain Strang sweeps.
        fluxes(columns, rows, true, q, volume, east, fx, options.reuseWorkspace);
        fluxes(columns, rows, false, q, volume, south, fy, options.reuseWorkspace);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const std::size_t i = static_cast<std::size_t>(y) * columns + x;
                const auto w = x > 0 ? i - 1 : i + columns - 1;
                vx[i] = volume[i] + east[w] - east[i];
                vy[i] = volume[i] + (y ? south[i - columns] : 0.0) - south[i];
                qx[i] = (q[i] * volume[i] + fx[w] - fx[i]) / vx[i];
                qy[i] = (q[i] * volume[i] + (y ? fy[i - columns] : 0.0) - fy[i]) / vy[i];
            }
        fluxes(columns, rows, true, qy, vy, east, gx, options.reuseWorkspace);
        fluxes(columns, rows, false, qx, vx, south, gy, options.reuseWorkspace);
        for (std::size_t i = 0; i < cells; ++i)
        {
            fx[i] = 0.5 * (fx[i] + gx[i]); fy[i] = 0.5 * (fy[i] + gy[i]);
            d.eastIntegratedFlux[i] += fx[i]; d.southIntegratedFlux[i] += fy[i];
        }
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const std::size_t i = static_cast<std::size_t>(y) * columns + x;
                const auto w = x > 0 ? i - 1 : i + columns - 1;
                q[i] += (fx[w] - fx[i] + (y ? fy[i - columns] : 0.0) - fy[i]) / volume[i];
                if (!std::isfinite(q[i]) || q[i] < -1.0e-8)
                    throw std::runtime_error("FFSL departure geometry lost positivity");
                q[i] = std::max(0.0, q[i]); // Roundoff only, never a mass fixer.
            }
    }
    destination.resize(cells);
    std::transform(q.begin(), q.end(), destination.begin(), [](double value) { return static_cast<float>(value); });
    d.minimumMixingRatio = *std::min_element(destination.begin(), destination.end());
    for (std::size_t i = 0; i < cells; ++i) d.finalAreaWeightedMass += destination[i] * volume[i];
    return d;
}
}
