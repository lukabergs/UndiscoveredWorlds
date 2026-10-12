// terrain_sandbox: the physics of the terrain sandbox (scripts/sandbox/terrain_sandbox.py). Routes and
// erodes a raster the sandbox writes, with the pipeline's own kernels (Priority-Flood, D8 flow graph,
// implicit stream power with deposition and hillslope diffusion), so procedural terrain can be judged
// after real drainage and erosion without running the pipeline. Not a pipeline stage: no world, no
// bundle, no tuning header. Files are raw little-endian float32, rows north first; elevations are metres
// relative to the sea level.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/glacial.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double secondsperyear = 365.25 * 86400.0;

std::vector<float> readraster(const std::string& path, std::size_t count)
{
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("cannot open " + path);
    std::vector<float> values(count);
    in.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(count * sizeof(float)));
    if (in.gcount() != static_cast<std::streamsize>(count * sizeof(float)))
        throw std::runtime_error(path + ": expected " + std::to_string(count) + " float32 values");
    return values;
}

void writeraster(const std::string& path, const std::vector<float>& values)
{
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(float)));
    if (!out)
        throw std::runtime_error("cannot write " + path);
}

// A window of square-in-degree pixels from latitude `north` down, `pixeldeg` per pixel. Longitude still
// wraps between the window's west and east columns: the caller makes the outer ring sea so nothing flows
// across that seam.
landscape::SphereMetric windowmetric(int columns, int rows, double north_deg, double pixel_deg, double radius_m)
{
    landscape::SphereMetric metric(columns, rows, radius_m);
    const double pixel = pixel_deg * pi / 180.0;
    metric.dy_m = radius_m * pixel;
    for (int j = 0; j < rows; j++)
    {
        const double north = north_deg * pi / 180.0 - j * pixel;
        const double south = north - pixel;
        const std::size_t row = static_cast<std::size_t>(j);
        metric.rowlatitude_rad[row] = 0.5 * (north + south);
        metric.rowdx_m[row] = radius_m * pixel * std::cos(metric.rowlatitude_rad[row]);
        metric.rowarea_m2[row] = radius_m * radius_m * pixel * (std::sin(north) - std::sin(south));
    }
    return metric;
}

// Threshold hillslopes: no land cell stands more than `slope` x distance above any neighbour, the sea
// counting at the sea level (0 m). A Dijkstra pass from the lowest cells up: each cell, once final, caps its
// land neighbours at its height plus slope x distance, so valley walls retreat as the rivers cut instead of
// leaving one-cell slots. Returns the volume removed (km^3); it leaves the system as landslide debris.
double limitslopes(const landscape::SphereMetric& metric, const std::vector<uint8_t>& sea, double slope,
                   std::vector<float>& elevation)
{
    using entry = std::pair<float, std::size_t>;
    std::priority_queue<entry, std::vector<entry>, std::greater<entry>> queue;
    std::vector<float> height(elevation.size());
    for (std::size_t i = 0; i < elevation.size(); i++)
    {
        height[i] = sea[i] ? 0.0f : elevation[i];
        queue.emplace(height[i], i);
    }
    std::vector<uint8_t> done(elevation.size(), 0);
    double removed = 0.0;
    while (!queue.empty())
    {
        const auto [value, cell] = queue.top();
        queue.pop();
        if (done[cell] || value != height[cell])
            continue;
        done[cell] = 1;
        for (int direction = 1; direction <= landscape::directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || done[next] || sea[next])
                continue;
            const float cap = static_cast<float>(value + slope * metric.distance_m(cell, direction));
            if (cap < height[next])
            {
                removed += (height[next] - cap) * metric.area_m2(next) * 1e-9;
                height[next] = cap;
                queue.emplace(cap, next);
            }
        }
    }
    for (std::size_t i = 0; i < elevation.size(); i++)
        if (!sea[i])
            elevation[i] = height[i];
    return removed;
}

// The sea floods every land cell at or below the sea level (0 m) connected to it: a bed cut below the sea
// by a glacier (a fjord) or a river mouth becomes sea, as the routing's mask requires (depressions.hpp).
void floodsea(const landscape::SphereMetric& metric, const std::vector<float>& elevation, std::vector<uint8_t>& sea)
{
    std::vector<std::size_t> stack;
    for (std::size_t i = 0; i < sea.size(); i++)
        if (sea[i])
            stack.push_back(i);
    while (!stack.empty())
    {
        const std::size_t cell = stack.back();
        stack.pop_back();
        for (int direction = 1; direction <= landscape::directioncount; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next != cell && !sea[next] && elevation[next] <= 0.0f)
            {
                sea[next] = 1;
                stack.push_back(next);
            }
        }
    }
}

void usage()
{
    std::cout << "terrain_sandbox --columns W --rows H --elevation FILE --sea FILE --out-elevation FILE --out-discharge FILE\n"
                 "  [--north DEG --pixel-deg DEG] (a window; default the whole sphere)\n"
                 "  [--erodibility FILE] (relative K per cell, default 1) [--uplift FILE] (m/Myr per cell, default 0)\n"
                 "  [--runoff M_PER_YR 0.3] [--duration MYR 10] [--step MYR 0.5] [--k 750] [--m 0.5]\n"
                 "  [--deposition G 1] [--diffusivity KM2_PER_MYR 1] [--delta-top M 0] [--threshold-slope S 0: off]\n"
                 "  [--ela FILE] (m per cell: glaciers) [--glacial-ratio 2] [--ice-full-range M 500] [--ice-accumulation 1]\n"
                 "  [--out-ice FILE]\n"
                 "  [--radius-km 6371]\n";
}
}

int main(int argc, char** argv)
{
    std::map<std::string, std::string> options;
    for (int i = 1; i + 1 < argc; i += 2)
        options[argv[i]] = argv[i + 1];
    const auto text = [&](const char* key, const char* fallback) {
        const auto found = options.find(key);
        return found == options.end() ? std::string(fallback) : found->second;
    };
    const auto number = [&](const char* key, double fallback) {
        const auto found = options.find(key);
        return found == options.end() ? fallback : std::stod(found->second);
    };
    if (argc < 2 || options.count("--columns") == 0 || options.count("--elevation") == 0)
    {
        usage();
        return 2;
    }
    try
    {
        const int columns = std::stoi(text("--columns", "0"));
        const int rows = std::stoi(text("--rows", "0"));
        const double radius_m = number("--radius-km", 6371.0) * 1000.0;
        const std::size_t count = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
        landscape::SphereMetric metric = options.count("--pixel-deg")
            ? windowmetric(columns, rows, number("--north", 90.0), number("--pixel-deg", 180.0 / rows), radius_m)
            : landscape::SphereMetric(columns, rows, radius_m);

        std::vector<float> elevation = readraster(text("--elevation", ""), count);
        const std::vector<float> seavalues = readraster(text("--sea", ""), count);
        std::vector<uint8_t> sea(count);
        for (std::size_t i = 0; i < count; i++)
            sea[i] = seavalues[i] > 0.5f ? 1 : 0;
        const double k = number("--k", 750.0);
        std::vector<float> erodibility = options.count("--erodibility") ? readraster(text("--erodibility", ""), count)
                                                                         : std::vector<float>(count, 1.0f);
        for (float& value : erodibility)
            value = static_cast<float>(value * k);
        const std::vector<float> uplift = options.count("--uplift") ? readraster(text("--uplift", ""), count)
                                                                     : std::vector<float>(count, 0.0f);
        const std::vector<float> runoff(count, static_cast<float>(number("--runoff", 0.3)));

        landscape::StreamPowerParameters parameters;
        parameters.areaexponent = number("--m", 0.5);
        parameters.depositioncoefficient = number("--deposition", 1.0);
        parameters.hillslopediffusivity_km2_myr = number("--diffusivity", 1.0);
        parameters.deltatop_m = number("--delta-top", 0.0);  // no sea cell filled above the sea level
        const double duration = number("--duration", 10.0);
        const double step = number("--step", 0.5);

        // Glaciers (erosion plan Phase G's kernel): an ELA per cell from the sandbox, the ice flux and thickness
        // on the ice surface's flow graph, the glacial term in the stream-power step (beds may be cut below the
        // sea: fjords).
        const bool glacial = options.count("--ela") != 0;
        const std::vector<float> ela = glacial ? readraster(text("--ela", ""), count) : std::vector<float>();
        landscape::GlacialParameters glacialparameters;
        glacialparameters.erodibilityratio = number("--glacial-ratio", 2.0);
        glacialparameters.fullicerange_m = number("--ice-full-range", 500.0);
        // The ice's own water supply: the runoff times --ice-accumulation (the wet coasts glaciers grow on; the
        // sandbox has no precipitation field). The ice flux may then exceed the rivers' discharge.
        const std::vector<float> snowfall(count, static_cast<float>(number("--runoff", 0.3) * number("--ice-accumulation", 1.0)));
        std::vector<double> iceflux(count, 0.0);
        std::vector<float> thickness(count, 0.0f), surface(count);

        std::vector<float> erosion(count), deposition(count);
        landscape::FlowGraph graph;
        const auto route = [&](const std::vector<float>& terrain) {
            floodsea(metric, elevation, sea);
            const landscape::Depressions depressions = landscape::filldepressions(metric, terrain, sea, 0.0f);
            landscape::RoutingOptions routing;
            routing.resolveflats = true;  // flats drain to their low edges, not in straight flood-order lines
            graph = landscape::buildflowgraph(metric, terrain, sea, 0.0f, depressions, routing);
            landscape::accumulateflow(metric, runoff, graph);
        };
        const double thresholdslope = number("--threshold-slope", 0.0);
        double eroded = 0.0, deposited = 0.0, marine = 0.0, landslides = 0.0;
        int steps = 0;
        for (double time = 0.0; time < duration - 1e-9; time += step)
        {
            parameters.step_myr = std::min(step, duration - time);
            for (std::size_t i = 0; i < count; i++)
                surface[i] = elevation[i] + thickness[i];
            route(surface);
            landscape::GlacialInput ice;
            if (glacial)
            {
                landscape::accumulateice(metric, graph, surface, ela, snowfall, glacialparameters.fullicerange_m, iceflux);
                landscape::icethickness(iceflux, glacialparameters, thickness);
                ice.iceflux_m3_yr = &iceflux;
                ice.thickness_m = &thickness;
                ice.erodibilityratio = glacialparameters.erodibilityratio;
                ice.flotationratio = glacialparameters.flotationratio;
            }
            const landscape::StreamPowerBudget budget = landscape::stepstreampower(
                metric, graph, erodibility, uplift, parameters, elevation, erosion, deposition, nullptr,
                glacial ? &ice : nullptr);
            eroded += budget.eroded_km3;
            deposited += budget.landdeposited_km3;
            marine += budget.marinedeposited_km3;
            if (thresholdslope > 0.0)
                landslides += limitslopes(metric, sea, thresholdslope, elevation);
            steps++;
        }
        route(elevation);
        std::vector<float> discharge(count, 0.0f);
        for (std::size_t i = 0; i < count; i++)
            discharge[i] = static_cast<float>(graph.discharge_m3_yr[i] / secondsperyear);
        writeraster(text("--out-elevation", ""), elevation);
        writeraster(text("--out-discharge", ""), discharge);
        if (options.count("--out-ice"))
            writeraster(text("--out-ice", ""), thickness);
        std::printf("steps %d eroded_km3 %.1f land_deposited_km3 %.1f marine_deposited_km3 %.1f landslides_km3 %.1f\n", steps,
                    eroded, deposited, marine, landslides);
    }
    catch (const std::exception& error)
    {
        std::cerr << "terrain_sandbox: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
