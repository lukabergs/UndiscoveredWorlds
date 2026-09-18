#include "climate_hydrology.hpp"
#include "climate_transport.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
#include <string>

namespace
{
constexpr double pi = 3.14159265358979323846;
int failures = 0;
void expect(bool value, const char* message)
{
    if (!value) { std::cerr << "FAIL: " << message << '\n'; ++failures; }
}

void checkFluxes(int columns, double radius, const std::vector<float>& initial,
    const std::vector<float>& result, const climatehydrology::SphericalTracerTransportDiagnostics& d)
{
    const auto grid = climategrid::makeSphericalGrid(columns, columns / 2, radius);
    expect(d.minimumMixingRatio >= 0.0f && std::isfinite(d.finalAreaWeightedMass), "transport remains finite and positive");
    expect(std::abs(d.finalAreaWeightedMass / d.initialAreaWeightedMass - 1.0) < 2.0e-6, "transport conserves spherical mass");
    double maximumResidual = 0.0;
    for (int y = 0; y < columns / 2; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto i = grid.index(x, y);
            const double net = d.eastIntegratedFlux[grid.index(x - 1, y)] - d.eastIntegratedFlux[i] +
                (y ? d.southIntegratedFlux[grid.index(x, y - 1)] : 0.0) - d.southIntegratedFlux[i];
            maximumResidual = std::max(maximumResidual, std::abs(result[i] - initial[i] - net / grid.cellAreasSquareMetres[y]));
        }
    expect(maximumResidual < 2.0e-5, "exported local face fluxes reproduce every cell change");
}

double rotationError(int columns, bool ffsl = false, double tilt = pi / 4.0, double phase = 0.0)
{
    const int rows = columns / 2;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 1000.0);
    std::vector<float> initial(columns * rows), u(initial.size()), v(initial.size()), next;
    // Tilted solid-body rotation crosses the seam and both polar regions.
    const double omega = 2.0 * pi / 100.0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const double lon = 2.0 * pi * (x + 0.5) / columns;
            const double lat = grid.latitudeCentresRadians[y];
            const auto i = grid.index(x, y);
            initial[i] = static_cast<float>(1.0 + 0.4 * std::cos(lat) * std::cos(lon + phase));
            u[i] = static_cast<float>(1000.0 * omega * (std::cos(tilt) * std::cos(lat) - std::sin(tilt) * std::sin(lat) * std::cos(lon)));
            v[i] = static_cast<float>(-1000.0 * omega * std::sin(tilt) * std::sin(lon));
        }
    auto q = initial;
    climatehydrology::MpdataOptions options;
    for (int step = 0; step < 4 * columns; ++step)
    {
        const auto d = ffsl
            ? climatehydrology::advectSphericalTracerFfsl(columns, rows, q, u, v,
                100.0f / (4 * columns), 1000.0f, {}, next)
            : climatehydrology::advectSphericalTracerMpdata(columns, rows, q, u, v,
                100.0f / (4 * columns), 1000.0f, options, next);
        checkFluxes(columns, 1000.0, q, next, d);
        q.swap(next);
    }
    double error = 0.0, area = 0.0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            error += grid.cellAreasSquareMetres[y] * std::pow(q[grid.index(x, y)] - initial[grid.index(x, y)], 2);
            area += grid.cellAreasSquareMetres[y];
        }
    return std::sqrt(error / area);
}

int benchmark(const char* path)
{
    constexpr int columns = 128, rows = 64, repetitions = 5;
    std::ifstream input(path);
    if (!input) { std::cerr << "Cannot read " << path << '\n'; return 2; }
    std::string line;
    std::getline(input, line);
    std::vector<float> u[8], v[8];
    while (std::getline(input, line))
    {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream values(line);
        int season, layer;
        double lat, lon, samples, duration, effective, east, south;
        if (!(values >> season >> layer >> lat >> lon >> samples >> duration >> effective >> east >> south) ||
            season < 0 || season > 3 || layer < 0 || layer > 1) return 2;
        u[2 * season + layer].push_back(static_cast<float>(east));
        v[2 * season + layer].push_back(static_cast<float>(south));
    }
    std::vector<float> source(columns * rows), result;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            source[y * columns + x] = static_cast<float>(10.0 + 5.0 * std::cos(2.0 * pi * x / columns) * std::cos(pi * (y + 0.5) / rows - pi / 2.0));
    std::cout << "method,season,layer,repetitions,ms_per_call,substeps,cell_substeps,mass_relative_error\n";
    for (bool ffsl : {false, true})
        for (int field = 0; field < 8; ++field)
        {
            if (u[field].size() != source.size()) { std::cerr << "Expected 128x64 wind fields\n"; return 2; }
            climatehydrology::MpdataOptions options;
            const auto start = std::chrono::steady_clock::now();
            climatehydrology::SphericalTracerTransportDiagnostics d;
            for (int repeat = 0; repeat < repetitions; ++repeat)
                d = ffsl ? climatehydrology::advectSphericalTracerFfsl(columns, rows, source, u[field], v[field], 28800.0f, 6371000.0f, {}, result)
                    : climatehydrology::advectSphericalTracerMpdata(columns, rows, source, u[field], v[field], 28800.0f, 6371000.0f, options, result);
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / repetitions;
            checkFluxes(columns, 6371000.0, source, result, d);
            std::cout << (ffsl ? "ffsl," : "mpdata,") << field / 2 << ',' << field % 2 << ',' << repetitions << ',' << ms << ','
                << d.substeps << ',' << d.cellSubsteps << ',' << d.finalAreaWeightedMass / d.initialAreaWeightedMass - 1.0 << '\n';
        }
    return failures ? 1 : 0;
}
}

int main(int argc, char** argv)
{
    // Explicit tier-4 microbenchmark; never run from default CTest.
    if (argc == 3 && std::string(argv[1]) == "--benchmark") return benchmark(argv[2]);
    const auto workspaceCheck = [](int phase)
    {
        // Grow/shrink and change radius/forcing while two callers share the row
        // executor. Cached geometry and scratch must never retain another call's data.
        for (int width : {32, 128, 16, 256, 128})
        {
            const auto count = width * (width / 2);
            std::vector<float> q(count), u(count), v(count), fresh, reused;
            for (int i = 0; i < count; ++i)
            {
                q[i] = static_cast<float>(2.0 + std::sin(0.03 * i + phase));
                u[i] = static_cast<float>(0.03 * std::cos(0.02 * i + phase));
                v[i] = static_cast<float>(0.01 * std::sin(0.01 * i));
            }
            climatehydrology::FfslOptions options;
            options.endZonalWindMps = v; options.endMeridionalWindMps = u;
            options.reuseWorkspace = false;
            const auto a = climatehydrology::advectSphericalTracerFfsl(width, width / 2, q, u, v,
                1.0f, 1000.0f + 10.0f * width + phase, options, fresh);
            options.reuseWorkspace = true;
            const auto b = climatehydrology::advectSphericalTracerFfsl(width, width / 2, q, u, v,
                1.0f, 1000.0f + 10.0f * width + phase, options, reused);
            if (fresh != reused || a.eastIntegratedFlux != b.eastIntegratedFlux ||
                a.southIntegratedFlux != b.southIntegratedFlux || a.substeps != b.substeps ||
                a.initialAreaWeightedMass != b.initialAreaWeightedMass || a.finalAreaWeightedMass != b.finalAreaWeightedMass)
                return false;
        }
        return true;
    };
    auto concurrent = std::async(std::launch::async, workspaceCheck, 1);
    const bool local = workspaceCheck(2);
    expect(local && concurrent.get(), "workspace reuse matches fresh storage bitwise across resizing and concurrent callers");
    const double coarse = rotationError(16), fine = rotationError(32);
    std::cout << "Tilted rotation RMS: " << coarse << " -> " << fine << '\n';
    expect(fine < 0.5 * coarse && fine < 0.008, "tilted rotation retains the original solver's accuracy under refinement");
    const double ffslCoarse = rotationError(16, true), ffslFine = rotationError(32, true);
    std::cout << "FFSL tilted rotation RMS: " << ffslCoarse << " -> " << ffslFine << '\n';
    expect(ffslFine < 0.5 * ffslCoarse && ffslFine < 0.008, "FFSL resolves tilted rotation under refinement");
    for (double tilt : {0.0, pi / 2.0})
    {
        // Keep the tracer perpendicular to the rotation axis so the polar
        // case actually transports it through both hemispheres.
        const double phase = tilt == 0.0 ? 0.0 : pi / 2.0;
        const double coarseAxis = rotationError(16, true, tilt, phase);
        const double fineAxis = rotationError(32, true, tilt, phase);
        std::cout << "FFSL rotation tilt=" << tilt << " RMS: " << coarseAxis << " -> " << fineAxis << '\n';
        expect(fineAxis < 0.6 * coarseAxis && fineAxis < 0.015,
            "both zonal and pole-crossing FFSL transport converge under refinement");
    }
    constexpr int columns = 32, rows = 16;
    std::vector<float> q(columns * rows), u(q.size()), v(q.size()), result;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto i = y * columns + x;
            q[i] = (x >= 8 && x < 16 && y > 3 && y < 12) ? 10.0f : 0.0f;
            u[i] = static_cast<float>(40.0 * std::sin(2.0 * pi * x / columns));
            v[i] = static_cast<float>(10.0 * std::sin(2.0 * pi * y / rows));
        }
    climatehydrology::MpdataOptions options;
    options.endZonalWindMps = u;
    options.endMeridionalWindMps = v;
    for (auto& value : options.endZonalWindMps) value = -value;
    for (auto& value : options.endMeridionalWindMps) value = -value;
    auto d = climatehydrology::advectSphericalTracerMpdata(columns, rows, q, u, v, 10.0f, 1000.0f, options, result);
    checkFluxes(columns, 1000.0, q, result, d);
    options.endZonalWindMps.clear(); options.endMeridionalWindMps.clear();
    std::fill(u.begin(), u.end(), 20.0f); std::fill(v.begin(), v.end(), 0.0f);
    d = climatehydrology::advectSphericalTracerMpdata(columns, rows, q, u, v, 10.0f, 1000.0f, options, result);
    checkFluxes(columns, 1000.0, q, result, d);
    expect(*std::max_element(result.begin(), result.end()) <= 10.00001f, "nondivergent transport creates no new extrema");
    expect(d.cellSubsteps == q.size() * d.substeps, "transport profiling records the actual CFL work");
    // Exact integer-cell translation at Courant 7: cost must not scale with
    // distance crossed, and the seam must preserve sharp fronts and dry cells.
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 1000.0);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            u[y * columns + x] = static_cast<float>(7.0 * grid.cellAreasSquareMetres[y] /
                grid.zonalFaceLengthsMetres[y]);
    d = climatehydrology::advectSphericalTracerFfsl(columns, rows, q, u, v, 1.0f, 1000.0f, {}, result);
    checkFluxes(columns, 1000.0, q, result, d);
    expect(d.substeps == 1, "FFSL crosses many cells without advective CFL subdivision");
    double shiftError = 0.0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            shiftError = std::max(shiftError, static_cast<double>(std::abs(result[grid.index(x, y)] - q[grid.index(x - 7, y)])));
    expect(shiftError < 1.0e-5, "large-Courant translation preserves the front at the correct longitude");
    expect(*std::max_element(result.begin(), result.end()) <= 10.00001f, "FFSL nondivergent transport creates no new extrema");
    // Compressing flow can exceed initial density extrema, but must remain
    // positive and locally conservative, including a reversal during the step.
    climatehydrology::FfslOptions ffsl;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            u[y * columns + x] = static_cast<float>(40.0 * std::sin(2.0 * pi * x / columns));
            v[y * columns + x] = static_cast<float>(10.0 * std::sin(2.0 * pi * y / rows));
        }
    d = climatehydrology::advectSphericalTracerFfsl(columns, rows, q, u, v, 10.0f, 1000.0f, ffsl, result);
    checkFluxes(columns, 1000.0, q, result, d);
    ffsl.endZonalWindMps = u; ffsl.endMeridionalWindMps = v;
    for (auto& value : ffsl.endZonalWindMps) value = -value;
    for (auto& value : ffsl.endMeridionalWindMps) value = -value;
    d = climatehydrology::advectSphericalTracerFfsl(columns, rows, q, u, v, 10.0f, 1000.0f, ffsl, result);
    checkFluxes(columns, 1000.0, q, result, d);
    return failures ? 1 : 0;
}
