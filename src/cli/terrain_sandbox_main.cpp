// terrain_sandbox: the physics of the terrain sandbox (scripts/sandbox/terrain_sandbox.py). Routes and
// erodes a raster the sandbox writes, with the pipeline's own kernels (Priority-Flood, D8 flow graph,
// implicit stream power with deposition and hillslope diffusion), so procedural terrain can be judged
// after real drainage and erosion without running the pipeline. Not a pipeline stage: no world, no
// bundle, no tuning header. Files are raw little-endian float32, rows north first; elevations are metres
// relative to today's sea level (--sea-history moves the sea during the run).
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
    // The cones of successive channel cells leave a rib per cell down a steep wall (the bed drops in steps):
    // the cut depth is averaged over 3 x 3 neighbourhoods twice, on the cells the limiter cut only, so the
    // walls retreat without ribs and no channel bed is lowered.
    std::vector<float> cut(elevation.size()), smoothed(elevation.size());
    for (std::size_t i = 0; i < elevation.size(); i++)
        cut[i] = sea[i] ? 0.0f : elevation[i] - height[i];
    for (int pass = 0; pass < 2; pass++)
    {
        for (std::size_t i = 0; i < elevation.size(); i++)
        {
            if (cut[i] <= 0.0f)
            {
                smoothed[i] = 0.0f;
                continue;
            }
            double sum = cut[i];
            int count = 1;
            for (int direction = 1; direction <= landscape::directioncount; direction++)
            {
                const std::size_t next = metric.neighbour(i, direction);
                if (next != i)
                {
                    sum += cut[next];
                    count++;
                }
            }
            smoothed[i] = static_cast<float>(sum / count);
        }
        cut.swap(smoothed);
    }
    // The averaged cut may exceed a cell's own: never below the lowest capped height around it (itself
    // included), so the smoothing cannot dig a pit. Without this floor, cuts averaged in from glacier-steep
    // walls took cells below their neighbours step after step: closed pits down to -15 km in B-core at k 256.
    removed = 0.0;
    for (std::size_t i = 0; i < elevation.size(); i++)
        if (!sea[i] && cut[i] > 0.0f)
        {
            float floor = height[i];
            for (int direction = 1; direction <= landscape::directioncount; direction++)
            {
                const std::size_t next = metric.neighbour(i, direction);
                if (next != i)
                    floor = std::min(floor, sea[next] ? 0.0f : height[next]);
            }
            const float lowered = std::max(elevation[i] - cut[i], floor);
            removed += (elevation[i] - lowered) * metric.area_m2(i) * 1e-9;
            elevation[i] = lowered;
        }
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

// A sea-level history: points (Myr from the start, metres) in time order, piecewise linear between them; a
// repeated time is a jump (the later point holds from that time on); constant beyond the ends.
struct SeaPoint
{
    double time_myr;
    double level_m;
};

std::vector<SeaPoint> parseseahistory(const std::string& text)
{
    std::vector<SeaPoint> points;
    std::size_t start = 0;
    while (start <= text.size())
    {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos)
            end = text.size();
        const std::string item = text.substr(start, end - start);
        const std::size_t colon = item.find(':');
        if (colon == std::string::npos)
            throw std::runtime_error("--sea-history: expected t:level, got '" + item + "'");
        points.push_back({std::stod(item.substr(0, colon)), std::stod(item.substr(colon + 1))});
        if (points.size() > 1 && points.back().time_myr < points[points.size() - 2].time_myr)
            throw std::runtime_error("--sea-history: times must not decrease");
        start = end + 1;
    }
    return points;
}

double sealevel(const std::vector<SeaPoint>& points, double time_myr)
{
    if (time_myr < points.front().time_myr)
        return points.front().level_m;
    std::size_t last = 0;  // the last point at or before the time
    while (last + 1 < points.size() && points[last + 1].time_myr <= time_myr)
        last++;
    if (last + 1 == points.size())
        return points.back().level_m;
    const SeaPoint& a = points[last];
    const SeaPoint& b = points[last + 1];
    return a.level_m + (b.level_m - a.level_m) * (time_myr - a.time_myr) / (b.time_myr - a.time_myr);
}

// The sea under a moving level: the 8-connected flood from the seeds over the cells at or below the level
// (0 m: elevations here are relative to the level). Rebuilt each routing, so the sea retreats as well as
// advances.
void floodfromseeds(const landscape::SphereMetric& metric, const std::vector<float>& elevation,
                    const std::vector<uint8_t>& seeds, std::vector<uint8_t>& sea)
{
    sea = seeds;
    floodsea(metric, elevation, sea);
}

void usage()
{
    std::cout << "terrain_sandbox --columns W --rows H --elevation FILE --sea FILE --out-elevation FILE --out-discharge FILE\n"
                 "  [--north DEG --pixel-deg DEG] (a window; default the whole sphere)\n"
                 "  [--erodibility FILE] (relative K per cell, default 1) [--uplift FILE] (m/Myr per cell, default 0)\n"
                 "  [--runoff M_PER_YR 0.3 | --runoff-file FILE] [--duration MYR 10] [--step MYR 0.5] [--k 750] [--m 0.5]\n"
                 "  [--deposition G 1] [--diffusivity KM2_PER_MYR 1] [--delta-top M 0] [--threshold-slope S 0: off]\n"
                 "  [--ela FILE] (m per cell: glaciers) [--glacial-ratio 2] [--ice-full-range M 500] [--ice-accumulation 1]\n"
                 "  [--out-ice FILE] [--pinned FILE] (cells held at their starting height) [--sea-cap M 0]\n"
                 "  [--sea-history \"t:level,...\"] (Myr, m: the sea level in time, piecewise linear; default 0)\n"
                 "  [--random-routing 0|1] [--seed N 1] [--radius-km 6371]\n";
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
        // A sea-level history (--sea-history): the sea is flooded each routing from fixed seeds, the starting sea
        // deeper than the lowest level (and a window's outer ring), so a lowstand exposes the shelf above it.
        // Without the option the level is 0 and the sea only grows (floodsea), as before.
        const bool history = options.count("--sea-history") != 0;
        const std::vector<SeaPoint> seahistory = history ? parseseahistory(text("--sea-history", "")) : std::vector<SeaPoint>();
        double lowestlevel = 0.0;
        std::vector<uint8_t> seaseeds;
        if (history)
        {
            lowestlevel = seahistory.front().level_m;
            for (const SeaPoint& point : seahistory)
                lowestlevel = std::min(lowestlevel, point.level_m);
            const bool window = options.count("--pixel-deg") != 0;
            seaseeds.assign(count, 0);
            for (std::size_t i = 0; i < count; i++)
            {
                const int column = static_cast<int>(i % static_cast<std::size_t>(columns));
                const int row = static_cast<int>(i / static_cast<std::size_t>(columns));
                const bool ring = column == 0 || column == columns - 1 || row == 0 || row == rows - 1;
                if (sea[i] && (elevation[i] < lowestlevel || (window && ring)))
                    seaseeds[i] = 1;
            }
        }
        // Pinned cells (a window's inner boundary ring): held at their starting height after every step, so a
        // river leaving the window drains to the land's height there, not to the sea level of the outer ring.
        std::vector<std::size_t> pinned;
        std::vector<float> pinnedheight;
        if (options.count("--pinned"))
        {
            const std::vector<float> mask = readraster(text("--pinned", ""), count);
            for (std::size_t i = 0; i < count; i++)
                if (mask[i] > 0.5f && !sea[i])
                {
                    pinned.push_back(i);
                    pinnedheight.push_back(elevation[i]);
                }
        }
        const auto holdpinned = [&]() {
            for (std::size_t j = 0; j < pinned.size(); j++)
                elevation[pinned[j]] = pinnedheight[j];
        };
        const double k = number("--k", 750.0);
        std::vector<float> erodibility = options.count("--erodibility") ? readraster(text("--erodibility", ""), count)
                                                                         : std::vector<float>(count, 1.0f);
        for (float& value : erodibility)
            value = static_cast<float>(value * k);
        const std::vector<float> uplift = options.count("--uplift") ? readraster(text("--uplift", ""), count)
                                                                     : std::vector<float>(count, 0.0f);
        // Runoff per cell (m/yr): a file from the sandbox's precipitation proxy, else uniform --runoff.
        const std::vector<float> runoff = options.count("--runoff-file")
            ? readraster(text("--runoff-file", ""), count)
            : std::vector<float>(count, static_cast<float>(number("--runoff", 0.3)));

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
        std::vector<float> snowfall(runoff);
        for (float& value : snowfall)
            value = static_cast<float>(value * number("--ice-accumulation", 1.0));
        std::vector<double> iceflux(count, 0.0);
        std::vector<float> thickness(count, 0.0f), surface(count);

        std::vector<float> erosion(count), deposition(count);
        landscape::FlowGraph graph;
        const bool randomrouting = number("--random-routing", 0.0) > 0.5;
        const long long routingseed = static_cast<long long>(number("--seed", 1.0)) * 1000003LL;
        int steps = 0;
        const auto route = [&](const std::vector<float>& terrain) {
            if (history)
                floodfromseeds(metric, elevation, seaseeds, sea);
            else
                floodsea(metric, elevation, sea);
            const landscape::Depressions depressions = landscape::filldepressions(metric, terrain, sea, 0.0f);
            landscape::RoutingOptions routing;
            routing.resolveflats = true;  // flats drain to their low edges, not in straight flood-order lines
            // Randomised receivers (flow_graph.hpp, Rho8 on the sphere's metric): the mean flow follows the gradient
            // instead of the nearest of eight directions, so rivers wander instead of running along the grid's
            // axes and diagonals; a new draw every step (the seed plus the step), deterministic.
            routing.randomised = randomrouting;
            routing.seed = static_cast<std::uint64_t>(routingseed + steps);
            graph = landscape::buildflowgraph(metric, terrain, sea, 0.0f, depressions, routing);
            landscape::accumulateflow(metric, runoff, graph);
        };
        const double thresholdslope = number("--threshold-slope", 0.0);
        const float seacap = static_cast<float>(number("--sea-cap", 0.0));
        double eroded = 0.0, deposited = 0.0, marine = 0.0, landslides = 0.0;
        // With a history, each step works in the frame of the level at its midpoint: the elevation and the ELA
        // less the level before the step, the level added back after it (the pinned cells are held in the
        // absolute frame, after that).
        std::vector<float> elarelative(history ? ela : std::vector<float>());
        const auto shift = [&](float level) {
            for (float& value : elevation)
                value -= level;
        };
        double seafractionmin = 1.0;
        const auto seafraction = [&]() {
            std::size_t cells = 0;
            for (const uint8_t value : sea)
                cells += value;
            return count ? static_cast<double>(cells) / static_cast<double>(count) : 0.0;
        };
        for (double time = 0.0; time < duration - 1e-9; time += step)
        {
            parameters.step_myr = std::min(step, duration - time);
            const float level = history ? static_cast<float>(sealevel(seahistory, time + 0.5 * parameters.step_myr)) : 0.0f;
            if (history)
            {
                shift(level);
                if (glacial)
                    for (std::size_t i = 0; i < count; i++)
                        elarelative[i] = ela[i] - level;
            }
            for (std::size_t i = 0; i < count; i++)
                surface[i] = elevation[i] + thickness[i];
            route(surface);
            if (history)
                seafractionmin = std::min(seafractionmin, seafraction());
            landscape::GlacialInput ice;
            if (glacial)
            {
                landscape::accumulateice(metric, graph, surface, history ? elarelative : ela, snowfall,
                                         glacialparameters.fullicerange_m, iceflux);
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
            // No sea cell is filled above --sea-cap (m, default the sea level): with the delta top lifted (its
            // fill-and-spill search went quadratic), marine deposits otherwise stacked whole loads on one cell,
            // towers up to 12 km in glacially cut troughs the sea had flooded (F-troughs, G-coast at k 256).
            // With a history the cap is relative to the level (the frame here), seed cells above it included.
            for (std::size_t i = 0; i < count; i++)
                if (sea[i])
                    elevation[i] = std::min(elevation[i], seacap);
            if (history)
                shift(-level);
            holdpinned();
            steps++;
        }
        // The final routing and the outputs at the level at the end of the history.
        const float finallevel = history ? static_cast<float>(seahistory.back().level_m) : 0.0f;
        if (history)
            shift(finallevel);
        route(elevation);
        if (history)
            shift(-finallevel);
        std::vector<float> discharge(count, 0.0f);
        for (std::size_t i = 0; i < count; i++)
            discharge[i] = static_cast<float>(graph.discharge_m3_yr[i] / secondsperyear);
        writeraster(text("--out-elevation", ""), elevation);
        writeraster(text("--out-discharge", ""), discharge);
        if (options.count("--out-ice"))
            writeraster(text("--out-ice", ""), thickness);
        std::printf("steps %d eroded_km3 %.1f land_deposited_km3 %.1f marine_deposited_km3 %.1f landslides_km3 %.1f", steps,
                    eroded, deposited, marine, landslides);
        if (history)  // the schedule's lowest and final levels, the sea's least and final extent (cell fractions)
            std::printf(" sea_min_m %.1f sea_final_m %.1f sea_fraction_min %.4f sea_fraction_final %.4f", lowestlevel,
                        static_cast<double>(finallevel), seafractionmin, seafraction());
        std::printf("\n");
    }
    catch (const std::exception& error)
    {
        std::cerr << "terrain_sandbox: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
