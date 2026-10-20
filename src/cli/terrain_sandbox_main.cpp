// terrain_sandbox: the physics of the terrain sandbox (scripts/sandbox/terrain_sandbox.py). Routes and
// erodes a raster the sandbox writes, with the pipeline's own kernels (Priority-Flood, D8 flow graph,
// implicit stream power with deposition and hillslope diffusion), so procedural terrain can be judged
// after real drainage and erosion without running the pipeline. Not a pipeline stage: no world, no
// tuning header. Files are raw little-endian float32: on a lat-lon domain (--columns, --rows) rows north
// first; on a cube domain (--face-cells N, docs/architecture/sphere-mesh.md section 2) 6 N^2 values in id
// order, the kernels running on SphereMetric(makecubegrid(N)). Elevations are metres relative to today's
// sea level (--sea-history moves the sea during the run).
// Subcommands (the first argument): `grid` writes a cube's cell table (the bundle v2 layout, section 7),
// `export` remaps cube fields to the lat-lon grid (section 5) for the renders, `synthetic-base` writes a
// small cube base bundle from analytic fields (the checks' stand-in until the engine writes cube worlds).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_remap.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/glacial.hpp"
#include "simulations/landscape/refinement_grid.hpp"
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
    metric.grid.reset();  // the window's own rows, not the whole-sphere grid's tables (columns = 2 x rows)
    return metric;
}

// Threshold hillslopes: no land cell stands more than `slope` x distance above any neighbour, the sea
// counting at the sea level (0 m). A Dijkstra pass from the lowest cells up: each cell, once final, caps its
// land neighbours at its height plus slope x distance, so valley walls retreat as the rivers cut instead of
// leaving one-cell slots. Returns the volume removed (km^3); it leaves the system as landslide debris.
// `cut_m`, when given, is assigned each cell's lowering (m, 0 where none); the result does not depend on it.
double limitslopes(const landscape::SphereMetric& metric, const std::vector<uint8_t>& sea, double slope,
                   std::vector<float>& elevation, std::vector<float>* cut_m = nullptr)
{
    if (cut_m)
        cut_m->assign(elevation.size(), 0.0f);
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
            if (cut_m)
                (*cut_m)[i] = elevation[i] - lowered;
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

// An ice history (--ice-history, roadmap T1.2b): phases (Myr from the start, an offset added to the ELA file in
// m, the glacial erodibility ratio), each holding from its time on (the first one before its time too); a
// large offset (100000) means no ice. Times must not decrease.
struct IcePhase
{
    double time_myr;
    double elaoffset_m;
    double ratio;
};

std::vector<IcePhase> parseicehistory(const std::string& text)
{
    std::vector<IcePhase> phases;
    std::size_t start = 0;
    while (start <= text.size())
    {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos)
            end = text.size();
        const std::string item = text.substr(start, end - start);
        const std::size_t first = item.find(':');
        const std::size_t second = first == std::string::npos ? std::string::npos : item.find(':', first + 1);
        if (second == std::string::npos)
            throw std::runtime_error("--ice-history: expected t:ela_offset_m:ratio, got '" + item + "'");
        phases.push_back({std::stod(item.substr(0, first)), std::stod(item.substr(first + 1, second - first - 1)),
                          std::stod(item.substr(second + 1))});
        if (phases.size() > 1 && phases.back().time_myr < phases[phases.size() - 2].time_myr)
            throw std::runtime_error("--ice-history: times must not decrease");
        start = end + 1;
    }
    return phases;
}

const IcePhase& icephase(const std::vector<IcePhase>& phases, double time_myr)
{
    std::size_t last = 0;  // the last phase starting at or before the time (the first before every start)
    while (last + 1 < phases.size() && phases[last + 1].time_myr <= time_myr)
        last++;
    return phases[last];
}

// Step lengths in time (--step-after "T:STEP,..."): from each T on the step is STEP Myr. Times must not decrease.
std::vector<std::pair<double, double>> parsestepafter(const std::string& text)
{
    std::vector<std::pair<double, double>> changes;
    std::size_t start = 0;
    while (start <= text.size())
    {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos)
            end = text.size();
        const std::string item = text.substr(start, end - start);
        const std::size_t colon = item.find(':');
        if (colon == std::string::npos)
            throw std::runtime_error("--step-after: expected T:STEP, got '" + item + "'");
        changes.emplace_back(std::stod(item.substr(0, colon)), std::stod(item.substr(colon + 1)));
        if (changes.back().second <= 0.0)
            throw std::runtime_error("--step-after: a step must be positive");
        if (changes.size() > 1 && changes.back().first < changes[changes.size() - 2].first)
            throw std::runtime_error("--step-after: times must not decrease");
        start = end + 1;
    }
    return changes;
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

// Lakes trap sediment (--lake-trap; roadmap T1.3, the user's "too many lakes" of 2026-10-14). The kernel routes water
// and sediment over the routing's filled surface, so a load reaching a closed depression passed through it to the
// spill and only the Yuan share G a / (A + G a) of it settled on each lake cell: depressions survived 20 Myr of rivers.
// Here every closed depression of the routing (a lake: the cells the Priority-Flood raised) keeps the load that
// reaches it: its exits (lake cells whose receiver lies outside it) are made land roots for the fluvial step, so the
// kernel lays that load on them; the load is then taken off the exits and spread over the lake from its floor up (a
// level fill on the routing's surface, bed plus ice) to at most its spill level, and only the excess passes on down
// the receivers from the spill: into the next lake downstream, or out of the land (an outlet: it leaves the system and
// is reported, as the sea's own deposits are not judged). Upstream lakes are filled first (their exits lie later in
// the flow order). The hillslope diffusion runs after the fill (a second kernel call without erosion or uplift), so it
// never spreads a load that is then taken back.
struct LakeTrap
{
    std::vector<int32_t> spill;        // an exit's former receiver (the spill side), -1 elsewhere
    std::vector<std::size_t> exits;
    double dumped_m3 = 0.0;            // the loads the kernel laid on the exits
    double placed_m3 = 0.0;            // what the fills laid in the lakes
    double overflow_m3 = 0.0;          // what passed out of the land
};

// Makes the exits of every closed, non-terminal depression of the routing land roots of `graph` (in place; restored by
// releaselakes). The order stays valid: a root needs no receiver before it. Depressions flagged in `skip` (per
// depression, when given) are left open: their load passes on as without the trap.
void traplakes(landscape::FlowGraph& graph, const landscape::Depressions& depressions, LakeTrap& trap,
               const std::vector<uint8_t>* skip = nullptr)
{
    const std::vector<int32_t>& id = depressions.depressionid;
    trap.spill.assign(graph.receiver.size(), -1);
    trap.exits.clear();
    for (int32_t entry : graph.order)
    {
        const std::size_t i = static_cast<std::size_t>(entry);
        const int32_t lake = id[i];
        if (lake == 0 || depressions.depressions[static_cast<std::size_t>(lake - 1)].terminal
            || (skip && (*skip)[static_cast<std::size_t>(lake - 1)]))
            continue;
        const int32_t receiver = graph.receiver[i];
        if (static_cast<std::size_t>(receiver) == i || id[static_cast<std::size_t>(receiver)] == lake)
            continue;
        trap.spill[i] = receiver;
        trap.exits.push_back(i);
        graph.receiver[i] = static_cast<int32_t>(i);
    }
}

void releaselakes(landscape::FlowGraph& graph, const LakeTrap& trap)
{
    for (std::size_t exit : trap.exits)
        graph.receiver[exit] = trap.spill[exit];
}

// A sloped lake fill (--lake-fill-gradient G, --lake-fill-noise FILE; roadmap T1 v18): a flat fill left flat lake
// floors along whose routing lines the ice cut straight parallel grooves (comb-like lakes). The full fill of a lake is
// its spill level plus an offset that rises from its exits: a Dijkstra pass within each lake from its exits, each move
// adding max(0.1 G d, G d + n_next - n_cell) (d the centre distance, n the noise in m, as the sandbox's sloped_fill
// of the procedural pits), so every filled cell has a lower neighbour towards an exit. Cells the pass does not reach
// keep offset 0.
std::vector<double> filloffsets(const landscape::SphereMetric& metric, const std::vector<int32_t>& id,
                                const std::vector<std::size_t>& exits, double gradient, const std::vector<float>& noise)
{
    std::vector<double> offset(id.size(), std::numeric_limits<double>::infinity());
    using Entry = std::pair<double, std::size_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
    for (std::size_t exit : exits)
    {
        offset[exit] = 0.0;
        queue.emplace(0.0, exit);
    }
    while (!queue.empty())
    {
        const auto [value, cell] = queue.top();
        queue.pop();
        if (value > offset[cell])
            continue;
        for (int direction = 1; direction <= 8; direction++)
        {
            const std::size_t next = metric.neighbour(cell, direction);
            if (next == cell || id[next] != id[cell])
                continue;
            const double run = gradient * metric.distance_m(cell, direction);
            const double rise = noise.empty() ? run
                : run + static_cast<double>(noise[next]) - static_cast<double>(noise[cell]);
            const double candidate = value + std::max(0.1 * run, rise);
            if (candidate < offset[next])
            {
                offset[next] = candidate;
                queue.emplace(candidate, next);
            }
        }
    }
    for (double& value : offset)
        if (!std::isfinite(value))
            value = 0.0;
    return offset;
}

// After the fluvial step on the trapped graph (graph restored): takes the loads off the exits and fills the lakes
// (level, or sloped by `gradient` > 0 with `noise`: filloffsets).
void filllakes(const landscape::SphereMetric& metric, const landscape::FlowGraph& graph,
               const landscape::Depressions& depressions, const std::vector<float>& thickness, LakeTrap& trap,
               std::vector<float>& elevation, std::vector<float>& deposition, double gradient,
               const std::vector<float>& noise)
{
    const std::vector<int32_t>& id = depressions.depressionid;
    const std::size_t lakes = depressions.depressions.size();
    if (trap.exits.empty())
        return;
    const std::vector<double> offset = gradient > 0.0 ? filloffsets(metric, id, trap.exits, gradient, noise)
                                                       : std::vector<double>();
    // Per lake: its cells, its main spill (the spill of its exit with the most water) and its place in the order.
    std::vector<std::vector<std::size_t>> cells(lakes + 1);
    for (std::size_t i = 0; i < id.size(); i++)
        if (id[i] != 0)
            cells[static_cast<std::size_t>(id[i])].push_back(i);
    std::vector<std::size_t> position(id.size(), 0);
    for (std::size_t k = 0; k < graph.order.size(); k++)
        position[static_cast<std::size_t>(graph.order[k])] = k;
    std::vector<uint8_t> trapped(lakes + 1, 0), processed(lakes + 1, 0);
    std::vector<int32_t> mainspill(lakes + 1, -1);
    std::vector<double> mainwater(lakes + 1, -1.0), pending(lakes + 1, 0.0);
    std::vector<std::size_t> last(lakes + 1, 0);
    for (std::size_t exit : trap.exits)
    {
        const std::size_t lake = static_cast<std::size_t>(id[exit]);
        trapped[lake] = 1;
        const double load = static_cast<double>(deposition[exit]) * metric.area_m2(exit);
        pending[lake] += load;
        trap.dumped_m3 += load;
        elevation[exit] -= deposition[exit];
        deposition[exit] = 0.0f;
        if (graph.discharge_m3_yr[exit] > mainwater[lake])
        {
            mainwater[lake] = graph.discharge_m3_yr[exit];
            mainspill[lake] = trap.spill[exit];
        }
        last[lake] = std::max(last[lake], position[exit]);
    }
    // A level fill of `volume` (m^3) into lake `lake` up to its spill level; returns what did not fit. Sloped: the
    // same fill on the surface less the offset (the fill's level `height` lies at height + offset on each cell).
    std::vector<std::pair<double, std::size_t>> floor;
    const auto fill = [&](std::size_t lake, double volume) {
        if (volume <= 0.0)
            return 0.0;
        const double level = static_cast<double>(depressions.depressions[lake - 1].spillelevation_m);
        floor.clear();
        for (std::size_t i : cells[lake])
        {
            double surface = static_cast<double>(elevation[i]) + static_cast<double>(thickness[i]);
            if (!offset.empty())
                surface -= offset[i];
            if (surface < level)
                floor.emplace_back(surface, i);
        }
        std::sort(floor.begin(), floor.end());
        double area = 0.0, moment = 0.0, height = level;
        bool full = true;
        for (std::size_t k = 0; k < floor.size(); k++)
        {
            area += metric.area_m2(floor[k].second);
            moment += metric.area_m2(floor[k].second) * floor[k].first;
            const double next = k + 1 < floor.size() ? floor[k + 1].first : level;
            if (area * next - moment >= volume)
            {
                height = (volume + moment) / area;
                full = false;
                break;
            }
        }
        double placed = 0.0;
        for (const auto& [surface, i] : floor)
        {
            if (surface >= height)
                break;
            const float bed = offset.empty() ? static_cast<float>(height - static_cast<double>(thickness[i]))
                                             : static_cast<float>(height + offset[i] - static_cast<double>(thickness[i]));
            if (bed <= elevation[i])
                continue;
            const float raise = bed - elevation[i];
            elevation[i] = bed;
            deposition[i] += raise;
            placed += static_cast<double>(raise) * metric.area_m2(i);
        }
        trap.placed_m3 += placed;
        return full ? std::max(0.0, volume - placed) : 0.0;
    };
    // The excess of a lake down the receivers from `cell`: into the next trapped lake (pending, or filled at once if it
    // was already processed) or out of the land.
    const auto carry = [&](double volume, int32_t start) {
        std::size_t cell = static_cast<std::size_t>(start);
        for (std::size_t steps = 0; volume > 0.0 && steps <= graph.order.size(); steps++)
        {
            const std::size_t lake = static_cast<std::size_t>(id[cell]);
            if (lake != 0 && trapped[lake])
            {
                if (!processed[lake])
                {
                    pending[lake] += volume;
                    return;
                }
                volume = fill(lake, volume);
                if (mainspill[lake] < 0)
                    break;
                cell = static_cast<std::size_t>(mainspill[lake]);
                continue;
            }
            const std::size_t receiver = static_cast<std::size_t>(graph.receiver[cell]);
            if (receiver == cell)
                break;
            cell = receiver;
        }
        trap.overflow_m3 += volume;
    };
    std::vector<std::size_t> sequence;
    for (std::size_t lake = 1; lake <= lakes; lake++)
        if (trapped[lake])
            sequence.push_back(lake);
    std::sort(sequence.begin(), sequence.end(), [&](std::size_t a, std::size_t b) {
        return last[a] != last[b] ? last[a] > last[b] : a < b;
    });
    for (std::size_t lake : sequence)
    {
        processed[lake] = 1;
        const double excess = fill(lake, pending[lake]);
        pending[lake] = 0.0;
        if (excess > 0.0)
            carry(excess, mainspill[lake]);
    }
}

// ---------------------------------------------------------------- cube grids, exports, the synthetic base

using Options = std::map<std::string, std::string>;

std::string option(const Options& options, const char* key, const char* fallback)
{
    const auto found = options.find(key);
    return found == options.end() ? std::string(fallback) : found->second;
}

double numberoption(const Options& options, const char* key, double fallback)
{
    const auto found = options.find(key);
    return found == options.end() ? fallback : std::stod(found->second);
}

std::uint64_t fnv1a(const void* data, std::size_t size)  // the bundle writer's hash (src/io/run_bundle.cpp)
{
    std::uint64_t hash = 1469598103934665603ULL;
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; index++)
    {
        hash ^= bytes[index];
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string hexhash(std::uint64_t hash)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

std::string shapetext(const std::vector<std::size_t>& shape, bool json)
{
    std::ostringstream out;
    out << (json ? '[' : '(');
    for (std::size_t index = 0; index < shape.size(); index++)
        out << (index ? ", " : "") << shape[index];
    if (!json && shape.size() == 1)
        out << ',';
    out << (json ? ']' : ')');
    return out.str();
}

// A version 1.0 .npy file, as the bundle writer's (src/io/run_bundle.cpp, writenpy).
void writenpy(const std::filesystem::path& path, const std::string& dtype, const std::vector<std::size_t>& shape,
              const void* data, std::size_t bytes)
{
    std::filesystem::create_directories(path.parent_path());
    std::string header = "{'descr': '" + dtype + "', 'fortran_order': False, 'shape': " + shapetext(shape, false) + ", }";
    const std::size_t prefix = 10;
    const std::size_t total = ((prefix + header.size() + 1 + 63) / 64) * 64;
    header.append(total - prefix - header.size() - 1, ' ');
    header += '\n';
    std::ofstream out(path, std::ios::binary);
    out.write("\x93NUMPY", 6);
    const unsigned char version[2] = { 1, 0 };
    out.write(reinterpret_cast<const char*>(version), 2);
    const auto length = static_cast<std::uint16_t>(header.size());
    const unsigned char lengthbytes[2] = { static_cast<unsigned char>(length & 0xff), static_cast<unsigned char>(length >> 8) };
    out.write(reinterpret_cast<const char*>(lengthbytes), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    if (bytes)
        out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out)
        throw std::runtime_error("cannot write " + path.string());
}

// One array written into a bundle directory, for its manifest entry.
struct WrittenArray
{
    std::string id, dtype, file;
    std::vector<std::size_t> shape;
    std::uint64_t hash = 0;
};

template <typename T>
WrittenArray savearray(const std::filesystem::path& directory, const std::string& relative, const std::string& id,
                       const std::string& dtype, const std::vector<std::size_t>& shape, const std::vector<T>& values)
{
    writenpy(directory / relative, dtype, shape, values.data(), values.size() * sizeof(T));
    return { id, dtype, relative, shape, fnv1a(values.data(), values.size() * sizeof(T)) };
}

// The cube's cell table in id order (bundle v2, sphere-mesh.md section 7; as world/planet_bundle.cpp writes it):
// centre (unit vectors), latitude_deg, longitude_deg, area_m2, neighbours (slots 1..8, -1 none), reverse (0 none).
std::vector<WrittenArray> writecelltable(const mesh::CellGrid& grid, const std::filesystem::path& directory,
                                         const std::string& prefix)
{
    const std::size_t count = grid.cellcount();
    constexpr double degrees = 180.0 / pi;
    std::vector<double> centre(3 * count), latitude(count), longitude(count), area(count);
    std::vector<std::int32_t> neighbours(mesh::slotcount * count);
    std::vector<std::int8_t> reverse(mesh::slotcount * count);
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        centre[3 * c] = p.x;
        centre[3 * c + 1] = p.y;
        centre[3 * c + 2] = p.z;
        latitude[c] = grid.latitude_rad(c) * degrees;
        longitude[c] = grid.longitude_rad(c) * degrees;
        area[c] = grid.area_m2(c);
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            neighbours[mesh::slotcount * c + static_cast<std::size_t>(k - 1)] = grid.neighbour(c, k);
            reverse[mesh::slotcount * c + static_cast<std::size_t>(k - 1)] = static_cast<std::int8_t>(grid.reverse(c, k));
        }
    }
    const std::size_t slots = static_cast<std::size_t>(mesh::slotcount);
    return { savearray(directory, prefix + "centre.npy", "centre", "<f8", { count, 3 }, centre),
             savearray(directory, prefix + "latitude_deg.npy", "latitude_deg", "<f8", { count }, latitude),
             savearray(directory, prefix + "longitude_deg.npy", "longitude_deg", "<f8", { count }, longitude),
             savearray(directory, prefix + "area_m2.npy", "area_m2", "<f8", { count }, area),
             savearray(directory, prefix + "neighbours.npy", "neighbours", "<i4", { count, slots }, neighbours),
             savearray(directory, prefix + "reverse.npy", "reverse", "|i1", { count, slots }, reverse) };
}

// `terrain_sandbox grid --face-cells N [--radius-km R] --out DIR`: the cell table of a cube of N cells per face edge
// as DIR/<name>.npy (the sandbox's cube domains read their neighbours from it; Python never builds the topology).
int gridcommand(const Options& options)
{
    const int facecells = std::stoi(option(options, "--face-cells", "0"));
    const auto grid = mesh::makecubegrid(facecells, numberoption(options, "--radius-km", 6371.0) * 1000.0);
    writecelltable(*grid, option(options, "--out", "."), "");
    std::printf("cells %zu\n", grid->cellcount());
    return 0;
}

// `terrain_sandbox export --face-cells N --columns W --input FILE --rules RULE,... --output FILE [--radius-km R]
// [--depth D]`: FILE holds F fields of 6 N^2 float32 each (id order), one after the other, and the output F lat-lon
// rasters of W x W/2 (rows north first), field f by RULE f: `bilinear` (mesh::samplebilinear's stencil at the
// lat-lon centres: smooth renders, section 5), `mean` (remapmean over makeremap's overlaps, quadtree depth D,
// default 4) or `nearest` (the cell containing the centre: masks, codes, angles, the discharge of the network).
int exportcommand(const Options& options)
{
    const double radius_m = numberoption(options, "--radius-km", 6371.0) * 1000.0;
    const int facecells = std::stoi(option(options, "--face-cells", "0"));
    const int columns = std::stoi(option(options, "--columns", "0"));
    const auto cube = mesh::makecubegrid(facecells, radius_m);
    const auto latlon = mesh::makelatlongrid(columns, radius_m);
    std::vector<std::string> rules;
    {
        const std::string text = option(options, "--rules", "");
        std::size_t start = 0;
        while (start < text.size())
        {
            std::size_t end = text.find(',', start);
            if (end == std::string::npos)
                end = text.size();
            rules.push_back(text.substr(start, end - start));
            start = end + 1;
        }
    }
    const std::size_t count = cube->cellcount(), destinations = latlon->cellcount();
    const std::vector<float> input = readraster(option(options, "--input", ""), count * rules.size());
    std::vector<float> output(destinations * rules.size());
    std::vector<mesh::BilinearStencil> stencils;
    std::vector<mesh::CellId> nearest;
    std::unique_ptr<mesh::RemapWeights> weights;
    for (std::size_t f = 0; f < rules.size(); f++)
    {
        const std::vector<float> values(input.begin() + static_cast<std::ptrdiff_t>(f * count),
                                        input.begin() + static_cast<std::ptrdiff_t>((f + 1) * count));
        float* result = output.data() + f * destinations;
        if (rules[f] == "bilinear")
        {
            if (stencils.empty())
            {
                stencils.resize(destinations);
                for (std::size_t d = 0; d < destinations; d++)
                    stencils[d] = mesh::bilinearstencil(*cube, latlon->centre(static_cast<mesh::CellId>(d)));
            }
            for (std::size_t d = 0; d < destinations; d++)
            {
                const mesh::BilinearStencil& stencil = stencils[d];
                const double first = static_cast<double>(values[stencil.cell[0]]);
                double sum = 0.0;  // relative to the first value, as mesh::samplebilinear
                for (int n = 1; n < stencil.count; n++)
                    sum += stencil.weight[n] * (static_cast<double>(values[stencil.cell[n]]) - first);
                result[d] = static_cast<float>(first + sum);
            }
        }
        else if (rules[f] == "nearest")
        {
            if (nearest.empty())
            {
                nearest.resize(destinations);
                for (std::size_t d = 0; d < destinations; d++)
                    nearest[d] = cube->locate(latlon->centre(static_cast<mesh::CellId>(d)));
            }
            for (std::size_t d = 0; d < destinations; d++)
                result[d] = values[nearest[d]];
        }
        else if (rules[f] == "mean")
        {
            if (!weights)
            {
                mesh::RemapOptions remap;
                remap.depth = static_cast<int>(numberoption(options, "--depth", 4.0));
                weights = std::make_unique<mesh::RemapWeights>(mesh::makeremap(*cube, *latlon, remap));
            }
            const std::vector<float> mean = mesh::remapmean(*weights, *cube, *latlon, values);
            std::copy(mean.begin(), mean.end(), result);
        }
        else
            throw std::runtime_error("export: unknown rule '" + rules[f] + "' (bilinear, mean, nearest)");
    }
    writeraster(option(options, "--output", ""), output);
    return 0;
}

// The level keys of a set of level-0 cells' descendants at levelcells = basecells x factor per face edge
// (sphere-mesh.md section 8: the per-face quadtree), ascending.
std::vector<mesh::LevelKey> descendantkeys(const std::vector<mesh::LevelKey>& base, std::int64_t basecells,
                                           std::int64_t factor)
{
    const std::int64_t level = basecells * factor;
    std::vector<mesh::LevelKey> keys;
    keys.reserve(base.size() * static_cast<std::size_t>(factor * factor));
    for (const mesh::LevelKey key : base)
    {
        const std::int64_t k = static_cast<std::int64_t>(key);
        const std::int64_t face = k / (basecells * basecells), row = (k / basecells) % basecells, column = k % basecells;
        for (std::int64_t b = 0; b < factor; b++)
            for (std::int64_t a = 0; a < factor; a++)
                keys.push_back(static_cast<mesh::LevelKey>(face * level * level + (row * factor + b) * level
                                                           + column * factor + a));
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

std::vector<int> intlist(const std::string& text)
{
    std::vector<int> values;
    std::size_t start = 0;
    while (start < text.size())
    {
        std::size_t end = text.find(',', start);
        if (end == std::string::npos)
            end = text.size();
        values.push_back(std::stoi(text.substr(start, end - start)));
        start = end + 1;
    }
    return values;
}

// `terrain_sandbox window --base-face-cells N --factor K --pixels x0,y0,x1,y1 --margin M [--radius-km R] --out DIR`:
// a sandbox window on a cube base (sphere-mesh.md section 8): the level-0 cells whose centres lie in the box of the
// base's lat-lon export pixels (landscape::cellsinpixelbox), the domain those dilated by M rings, and its cells at
// N K per face edge (their descendants) as mesh::cubewindowgrid (compact): DIR/keys.npy (u8, ascending), kept.npy
// (u1: the box's descendants) and the window's cell table (neighbours as local ids, -1 outside the window).
int windowcommand(const Options& options)
{
    const int basecells = std::stoi(option(options, "--base-face-cells", "0"));
    const int factor = std::stoi(option(options, "--factor", "0"));
    const double radius_m = numberoption(options, "--radius-km", 6371.0) * 1000.0;
    const std::vector<int> box = intlist(option(options, "--pixels", ""));
    if (box.size() != 4 || basecells < 1 || factor < 1)
        throw std::runtime_error("window: --base-face-cells, --factor and --pixels x0,y0,x1,y1 are required");
    const auto base = mesh::makecubegrid(basecells, radius_m);
    const std::vector<mesh::LevelKey> kept0 = landscape::cellsinpixelbox(*base, 4 * basecells, box[0], box[1], box[2], box[3]);
    const std::vector<mesh::LevelKey> domain0 = landscape::dilatecells(*base, kept0,
                                                                       std::stoi(option(options, "--margin", "1")));
    const std::vector<mesh::LevelKey> keys = descendantkeys(domain0, basecells, factor);
    const std::vector<mesh::LevelKey> keptkeys = descendantkeys(kept0, basecells, factor);
    const auto grid = mesh::cubewindowgrid(basecells * factor, radius_m, keys, true);
    std::vector<std::uint8_t> kept(keys.size(), 0);
    for (const mesh::LevelKey key : keptkeys)
        kept[static_cast<std::size_t>(std::lower_bound(keys.begin(), keys.end(), key) - keys.begin())] = 1;
    const std::filesystem::path directory = option(options, "--out", ".");
    savearray(directory, "keys.npy", "keys", "<u8", { keys.size() }, keys);
    savearray(directory, "kept.npy", "kept", "|u1", { kept.size() }, kept);
    writecelltable(*grid, directory, "");
    std::printf("window cells %zu kept %zu level-0 %zu kept level-0 %zu\n", keys.size(), keptkeys.size(), domain0.size(),
                kept0.size());
    return 0;
}

// `terrain_sandbox stencils --base-face-cells N --level-cells NL (--keys FILE | --all) --out PREFIX`: for each level
// cell (FILE: raw uint64 level keys; --all: every cell of the level) the G stencil of section 8 at its centre over the
// base cube (mesh::cubicstencil, true-position ghosts): PREFIX-offsets.i8 (count + 1), -ids.i4, -weights.f8 (CSR rows in
// the keys' order, entries in ascending base id), so the sandbox samples world fields by G.
int stencilscommand(const Options& options)
{
    const std::int64_t basecells = std::stoll(option(options, "--base-face-cells", "0"));
    const std::int64_t level = std::stoll(option(options, "--level-cells", "0"));
    std::vector<mesh::LevelKey> keys;
    if (options.count("--all"))
    {
        keys.resize(static_cast<std::size_t>(6 * level * level));
        for (std::size_t i = 0; i < keys.size(); i++)
            keys[i] = static_cast<mesh::LevelKey>(i);
    }
    else
    {
        std::ifstream in(option(options, "--keys", ""), std::ios::binary | std::ios::ate);
        if (!in)
            throw std::runtime_error("stencils: cannot open --keys");
        keys.resize(static_cast<std::size_t>(in.tellg()) / sizeof(mesh::LevelKey));
        in.seekg(0);
        in.read(reinterpret_cast<char*>(keys.data()), static_cast<std::streamsize>(keys.size() * sizeof(mesh::LevelKey)));
    }
    std::vector<std::int64_t> offsets(keys.size() + 1, 0);
    std::vector<std::int32_t> ids;
    std::vector<double> weights;
    ids.reserve(keys.size() * 16);
    weights.reserve(keys.size() * 16);
    for (std::size_t c = 0; c < keys.size(); c++)
    {
        const std::int64_t k = static_cast<std::int64_t>(keys[c]);
        const int face = static_cast<int>(k / (level * level));
        const int row = static_cast<int>((k / level) % level), column = static_cast<int>(k % level);
        const mesh::Vec3 p = mesh::cubepoint(face, mesh::cubecentrea(static_cast<int>(level), column),
                                             mesh::cubecentreb(static_cast<int>(level), row));
        const mesh::KeyStencil stencil = mesh::cubicstencil(basecells, p);
        for (std::size_t n = 0; n < stencil.key.size(); n++)
        {
            ids.push_back(static_cast<std::int32_t>(stencil.key[n]));
            weights.push_back(stencil.weight[n]);
        }
        offsets[c + 1] = static_cast<std::int64_t>(ids.size());
    }
    const std::string prefix = option(options, "--out", "");
    const auto write = [&](const std::string& path, const void* data, std::size_t bytes) {
        std::ofstream out(path, std::ios::binary);
        out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
        if (!out)
            throw std::runtime_error("cannot write " + path);
    };
    write(prefix + "-offsets.i8", offsets.data(), offsets.size() * sizeof(std::int64_t));
    write(prefix + "-ids.i4", ids.data(), ids.size() * sizeof(std::int32_t));
    write(prefix + "-weights.f8", weights.data(), weights.size() * sizeof(double));
    std::printf("cells %zu entries %zu\n", keys.size(), ids.size());
    return 0;
}

// The synthetic cube base (`terrain_sandbox synthetic-base --face-cells N [--radius-km R] --out DIR`): a bundle v2 of
// grid kind cube with the fields the sandbox reads, analytic functions of the unit vector: three continents (one over
// the north pole, one across a cube corner, one on the equator astride a face edge) joined by a collision belt along a
// great circle (shortening across it, uplift, thick crust, a convergent boundary) and a spreading ridge along another
// (young oceanic crust, a divergent boundary). Lat-lon exports at 4 N columns by each field's rule; seed 1; no steps.
int syntheticbasecommand(const Options& options)
{
    const int facecells = std::stoi(option(options, "--face-cells", "32"));
    const double radius_km = numberoption(options, "--radius-km", 6371.0);
    const std::filesystem::path directory = option(options, "--out", "");
    if (directory.empty())
        throw std::runtime_error("synthetic-base: --out DIR is required");
    if (std::filesystem::exists(directory / "manifest.json"))
        throw std::runtime_error("synthetic-base: " + directory.string() + " already holds a bundle");
    const auto grid = mesh::makecubegrid(facecells, radius_km * 1000.0);
    const auto latlon = mesh::makelatlongrid(4 * facecells, radius_km * 1000.0);
    const std::size_t count = grid->cellcount();
    const auto unit = [](double latitude_deg, double longitude_deg) {
        return mesh::unitvector(latitude_deg * pi / 180.0, longitude_deg * pi / 180.0);
    };
    const auto angle = [](const mesh::Vec3& p, const mesh::Vec3& q) {
        return std::acos(std::clamp(mesh::dot(p, q), -1.0, 1.0));
    };
    struct Continent
    {
        mesh::Vec3 centre;
        double radius_rad, height_m;
    };
    const Continent continents[] = { { unit(80.0, 30.0), 0.38, 900.0 },     // over the north pole
                                     { unit(35.26, 45.0), 0.42, 700.0 },    // across the cube corner of faces 0, 1, 4
                                     { unit(-10.0, 135.0), 0.45, 600.0 } }; // astride the edge of faces 1 and 2
    const mesh::Vec3 beltpole = mesh::normalise(mesh::cross(continents[1].centre, continents[2].centre));
    const mesh::Vec3 ridgepole = unit(10.0, -100.0);  // the ridge is the great circle 90 degrees from it
    std::vector<float> elevation(count), shortening(count), direction(count), uplift(count), sediment(count),
        crust(count), age(count), strain(count), zero(count, 0.0f);
    std::vector<std::uint8_t> crustclass(count), boundary(count), deforming(count);
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::Vec3 p = grid->centre(c);
        double land = -0.22;  // the land function: > 0 on land
        for (const Continent& continent : continents)
        {
            const double x = angle(p, continent.centre) / continent.radius_rad;
            land += std::exp(-x * x) * (continent.height_m / 900.0);
        }
        // The belt: within 6 degrees of the great circle through the second and third continents' centres, between them.
        const double offbelt = std::asin(std::clamp(mesh::dot(p, beltpole), -1.0, 1.0));
        const double between = mesh::dot(p, mesh::normalise(continents[1].centre + continents[2].centre));
        const double belt = between > 0.75 ? std::exp(-std::pow(offbelt / (6.0 * pi / 180.0), 2.0)) : 0.0;
        land += 0.6 * belt;
        const double ridgedistance = std::abs(std::asin(std::clamp(mesh::dot(p, ridgepole), -1.0, 1.0))) * radius_km;
        const bool sea = land <= 0.0;
        elevation[c] = static_cast<float>(sea ? std::max(-5500.0, -1500.0 - 14000.0 * (-land)) : 50.0 + 2200.0 * land + 3500.0 * belt);
        shortening[c] = static_cast<float>(120000.0 * belt);
        // The shortening axis across the belt (toward its pole), radians from east toward south in the cell's frame.
        const mesh::Vec3 across = mesh::normalise(beltpole - mesh::dot(beltpole, p) * p);
        direction[c] = static_cast<float>(std::atan2(-mesh::dot(across, grid->north(c)), mesh::dot(across, grid->east(c))));
        uplift[c] = static_cast<float>(250.0 * belt - (sea ? 0.0 : 20.0 * std::max(0.0, 0.3 - land)));
        crust[c] = static_cast<float>(sea && elevation[c] < -1500.0f ? 0.0 : 33000.0 + 22000.0 * belt);
        crustclass[c] = crust[c] > 0.0f ? 2 : 1;
        age[c] = static_cast<float>(crustclass[c] == 1 ? std::min(180.0, ridgedistance / 25.0) : 900.0);
        sediment[c] = static_cast<float>(sea ? 400.0 : 1500.0 * std::max(0.0, 0.25 - land) * 4.0);
        strain[c] = static_cast<float>(0.5 * belt);
        boundary[c] = belt > 0.8 ? 1 : (crustclass[c] == 1 && ridgedistance < 120.0 ? 2 : 0);
        deforming[c] = belt > 0.5 ? 2 : 0;
    }
    struct Field
    {
        std::string id, layout, dtype, units, rule;
        const void* data;
        std::size_t elementsize;
        bool code;
    };
    const std::vector<Field> fields = {
        { "geosphere.base_elevation", "grid", "<f4", "m above datum", "mean", elevation.data(), 4, false },
        { "geosphere.mountain_height", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.extra_elevation", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.crater_rims", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.crater_centres", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.volcano", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.ocean_ridge_height", "grid", "<f4", "m", "mean", zero.data(), 4, false },
        { "geosphere.crust_age", "grid", "<f4", "Myr", "mean", age.data(), 4, false },
        { "geosphere.crust_class", "grid", "|u1", "code", "nearest", crustclass.data(), 1, true },
        { "geosphere.rock_uplift_rate", "grid", "<f4", "m/Myr", "mean", uplift.data(), 4, false },
        { "geosphere.shortening", "grid", "<f4", "m", "mean", shortening.data(), 4, false },
        { "geosphere.shortening_direction", "grid", "<f4", "rad", "nearest", direction.data(), 4, false },
        { "geosphere.accumulated_strain", "grid", "<f4", "engine units", "mean", strain.data(), 4, false },
        { "geosphere.boundary_type", "grid", "|u1", "code", "nearest", boundary.data(), 1, true },
        { "geosphere.deforming_region_type", "grid", "|u1", "code", "nearest", deforming.data(), 1, true },
        { "geosphere.continental_crust", "grid", "<f4", "m", "mean", crust.data(), 4, false },
        { "geosphere.sediment", "grid", "<f4", "m", "mean", sediment.data(), 4, false },
    };
    mesh::RemapOptions remapoptions;
    remapoptions.depth = 5;
    const mesh::RemapWeights weights = mesh::makeremap(*grid, *latlon, remapoptions);
    const std::size_t n = static_cast<std::size_t>(facecells);
    const std::size_t columns = static_cast<std::size_t>(4 * facecells), rows = columns / 2;
    std::vector<WrittenArray> written, exports;
    std::vector<std::string> rules;
    for (const Field& field : fields)
    {
        const std::string name = field.id.substr(field.id.find('.') + 1);
        const std::size_t bytes = count * field.elementsize;
        writenpy(directory / ("fields/geosphere/" + name + ".npy"), field.dtype, { 6, n, n }, field.data, bytes);
        written.push_back({ field.id, field.dtype, "fields/geosphere/" + name + ".npy", { 6, n, n }, fnv1a(field.data, bytes) });
        std::vector<float> values(count);
        for (std::size_t i = 0; i < count; i++)
            values[i] = field.code ? static_cast<float>(static_cast<const std::uint8_t*>(field.data)[i])
                                   : static_cast<const float*>(field.data)[i];
        const std::vector<float> exported = field.rule == "mean" ? mesh::remapmean(weights, *grid, *latlon, values)
                                                                 : mesh::remapnearest(*grid, *latlon, values);
        const std::string file = "latlon/geosphere/" + name + ".npy";
        if (field.code)
        {
            std::vector<std::uint8_t> codes(exported.size());
            for (std::size_t i = 0; i < codes.size(); i++)
                codes[i] = static_cast<std::uint8_t>(exported[i]);
            exports.push_back(savearray(directory, file, field.id, field.dtype, { rows, columns }, codes));
        }
        else
            exports.push_back(savearray(directory, file, field.id, field.dtype, { rows, columns }, exported));
        rules.push_back(field.rule);
    }
    const std::int32_t sealevel = 0;
    writenpy(directory / "fields/geosphere/sea_level.npy", "<i4", {}, &sealevel, sizeof(sealevel));
    const WrittenArray sealevelentry{ "geosphere.sea_level", "<i4", "fields/geosphere/sea_level.npy", {},
                                      fnv1a(&sealevel, sizeof(sealevel)) };
    const std::vector<WrittenArray> cells = writecelltable(*grid, directory, "grid/");
    std::filesystem::create_directories(directory / "tables");

    std::ostringstream out;
    out << std::setprecision(6) << std::fixed;
    out << "{\n  \"schema\": \"uw-run-bundle/v2\",\n  \"created\": \"synthetic\",\n  \"seed\": 1,\n";
    out << "  \"grid\": {\"kind\": \"cube\", \"face_cells\": " << facecells
        << ", \"orientation\": \"uw-cube-1\", \"layout\": \"face, row (top first), column\", \"radius_km\": " << radius_km
        << ", \"export_columns\": " << columns << "},\n";
    out << "  \"recipe\": \"synthetic-cube-base\",\n  \"command\": [\"terrain_sandbox\", \"synthetic-base\"],\n";
    out << "  \"properties\": {\"note\": \"analytic fields of the unit vector (terrain_sandbox synthetic-base), not an engine world\"},\n";
    out << "  \"parameters\": {\"world.body.radius_km\": " << radius_km << "},\n  \"steps\": [],\n  \"fields\": [";
    const auto entry = [&](const WrittenArray& array, const std::string& layout, const std::string& units) {
        out << "\n    {\"id\": \"" << array.id << "\", \"component\": \"geosphere\", \"layout\": \"" << layout
            << "\", \"dtype\": \"" << array.dtype << "\", \"shape\": " << shapetext(array.shape, true) << ", \"units\": \""
            << units << "\", \"file\": \"" << array.file << "\", \"hash\": \"" << hexhash(array.hash)
            << "\", \"written_by\": \"terrain_sandbox synthetic-base\", \"description\": \"\"}";
    };
    entry(sealevelentry, "scalar", "m above datum");
    for (std::size_t i = 0; i < fields.size(); i++)
    {
        out << ",";
        entry(written[i], fields[i].layout, fields[i].units);
    }
    out << "\n  ],\n  \"exports\": [";
    for (std::size_t i = 0; i < exports.size(); i++)
        out << (i ? "," : "") << "\n    {\"id\": \"" << exports[i].id << "\", \"dtype\": \"" << exports[i].dtype
            << "\", \"shape\": " << shapetext(exports[i].shape, true) << ", \"file\": \"" << exports[i].file
            << "\", \"hash\": \"" << hexhash(exports[i].hash) << "\", \"rule\": \"" << rules[i] << "\"}";
    out << "\n  ],\n  \"cells\": [";
    for (std::size_t i = 0; i < cells.size(); i++)
        out << (i ? "," : "") << "\n    {\"id\": \"" << cells[i].id << "\", \"dtype\": \"" << cells[i].dtype
            << "\", \"shape\": " << shapetext(cells[i].shape, true) << ", \"file\": \"" << cells[i].file
            << "\", \"hash\": \"" << hexhash(cells[i].hash) << "\"}";
    out << "\n  ],\n  \"tables\": []\n}\n";
    std::ofstream manifest(directory / "manifest.json", std::ios::binary);
    manifest << out.str();
    if (!manifest)
        throw std::runtime_error("cannot write the manifest in " + directory.string());
    std::printf("synthetic cube base N %d: %zu cells in %s\n", facecells, count, directory.string().c_str());
    return 0;
}

void usage()
{
    std::cout << "terrain_sandbox --columns W --rows H --elevation FILE --sea FILE --out-elevation FILE --out-discharge FILE\n"
                 "  [--north DEG --pixel-deg DEG] (a window; default the whole sphere)\n"
                 "  [--face-cells N] (instead of --columns and --rows: a cube domain, 6 N^2 cells in id order)\n"
                 "  [--cells FILE] (with --face-cells: raw int32 neighbours, n x 8, checked against the cube's table)\n"
                 "  [--erodibility FILE] (relative K per cell, default 1) [--uplift FILE] (m/Myr per cell, default 0)\n"
                 "  [--runoff M_PER_YR 0.3 | --runoff-file FILE] [--duration MYR 10] [--step MYR 0.5] [--k 750] [--m 0.5]\n"
                 "  [--deposition G 1] [--diffusivity KM2_PER_MYR 1] [--delta-top M 0] [--threshold-slope S 0: off]\n"
                 "  [--ela FILE] (m per cell: glaciers) [--glacial-ratio 2] [--ice-full-range M 500] [--ice-accumulation 1]\n"
                 "  [--out-ice FILE] [--pinned FILE] (cells held at their starting height) [--sea-cap M 0]\n"
                 "  [--sea-history \"t:level,...\"] (Myr, m: the sea level in time, piecewise linear; default 0)\n"
                 "  [--no-infill-after MYR] (from then on no deposition inside closed depressions of the bed)\n"
                 "  [--lake-trap 0|1] (lakes keep the sediment reaching them, filled up to their spill; not in no-infill steps)\n"
                 "  [--lake-fill-gradient G 0] [--lake-fill-noise FILE] (m per cell) (the trap's fill sloped up from the exits)\n"
                 "  [--lake-trap-ice-share S 2] (a lake with that share of its cells under ice does not trap; above 1: all trap)\n"
                 "  [--ice-history \"t:ela_offset_m:ratio,...\"] (Myr, m added to the ELA, glacial ratio: piecewise constant)\n"
                 "  [--step-after \"T:STEP,...\"] (from T Myr on, steps of STEP Myr)\n"
                 "  [--aggrade-subsidence 0|1] (land cells' negative uplift filled by sediment as they subside)\n"
                 "  [--out-fields PREFIX] (per-cell totals: PREFIX-eroded_m.f32, -deposited_m, -glacial_eroded_m,\n"
                 "   -ice_time_myr, -ice_max_m, -landslide_m, -last_eroded_m; -basin_fill_m with --aggrade-subsidence 1)\n"
                 "  [--ice-threshold M 10]\n"
                 "  [--random-routing 0|1] [--seed N 1] [--radius-km 6371]\n"
                 "terrain_sandbox grid --face-cells N [--radius-km R] --out DIR (the cube's cell table as .npy)\n"
                 "terrain_sandbox export --face-cells N --columns W --input FILE --rules bilinear|mean|nearest,...\n"
                 "  --output FILE [--radius-km R] [--depth D 4] (F stacked cube fields to F stacked lat-lon rasters)\n"
                 "terrain_sandbox synthetic-base --face-cells N [--radius-km R] --out DIR (an analytic cube base bundle)\n"
                 "terrain_sandbox window --base-face-cells N --factor K --pixels x0,y0,x1,y1 --margin M [--radius-km R]\n"
                 "  --out DIR (a cube window's keys, kept mask and cell table)\n"
                 "terrain_sandbox stencils --base-face-cells N --level-cells NL (--keys FILE | --all) --out PREFIX\n"
                 "  (the G stencils of the level cells over the base, CSR)\n"
                 "physics on a cube window: --level-cells NL --window-keys FILE (raw uint64 keys) instead of --face-cells\n";
}
}

int main(int argc, char** argv)
{
    const std::string command = argc > 1 && argv[1][0] != '-' ? argv[1] : "";
    std::map<std::string, std::string> options;
    for (int i = command.empty() ? 1 : 2; i + 1 < argc; i += 2)
        options[argv[i]] = argv[i + 1];
    const auto text = [&](const char* key, const char* fallback) {
        const auto found = options.find(key);
        return found == options.end() ? std::string(fallback) : found->second;
    };
    const auto number = [&](const char* key, double fallback) {
        const auto found = options.find(key);
        return found == options.end() ? fallback : std::stod(found->second);
    };
    if (!command.empty())
    {
        try
        {
            if (command == "grid")
                return gridcommand(options);
            if (command == "export")
                return exportcommand(options);
            if (command == "synthetic-base")
                return syntheticbasecommand(options);
            if (command == "window")
                return windowcommand(options);
            if (command == "stencils")
                return stencilscommand(options);
        }
        catch (const std::exception& error)
        {
            std::cerr << "terrain_sandbox " << command << ": " << error.what() << "\n";
            return 1;
        }
        usage();
        return 2;
    }
    // A cube window (--level-cells NL --window-keys FILE: raw uint64 level keys, ascending): mesh::cubewindowgrid.
    const bool cubewindow = options.count("--window-keys") != 0;
    const bool cube = options.count("--face-cells") != 0 || cubewindow;
    if (argc < 2 || (options.count("--columns") == 0 && !cube) || options.count("--elevation") == 0)
    {
        usage();
        return 2;
    }
    try
    {
        const int columns = cube ? 0 : std::stoi(text("--columns", "0"));
        const int rows = cube ? 0 : std::stoi(text("--rows", "0"));
        const double radius_m = number("--radius-km", 6371.0) * 1000.0;
        if (cube && options.count("--pixel-deg"))
            throw std::runtime_error("--pixel-deg is a lat-lon window's; a cube domain is the whole sphere");
        // A cube domain: the kernels on SphereMetric(makecubegrid(N)); --cells checks the caller's neighbour table.
        std::vector<mesh::LevelKey> windowkeys;
        if (cubewindow)
        {
            std::ifstream in(text("--window-keys", ""), std::ios::binary | std::ios::ate);
            if (!in)
                throw std::runtime_error("cannot open --window-keys");
            windowkeys.resize(static_cast<std::size_t>(in.tellg()) / sizeof(mesh::LevelKey));
            in.seekg(0);
            in.read(reinterpret_cast<char*>(windowkeys.data()),
                    static_cast<std::streamsize>(windowkeys.size() * sizeof(mesh::LevelKey)));
        }
        landscape::SphereMetric metric = cubewindow
            ? landscape::SphereMetric(mesh::cubewindowgrid(std::stoi(text("--level-cells", "0")), radius_m, windowkeys, true))
            : cube
            ? landscape::SphereMetric(mesh::makecubegrid(std::stoi(text("--face-cells", "0")), radius_m))
            : options.count("--pixel-deg")
            ? windowmetric(columns, rows, number("--north", 90.0), number("--pixel-deg", 180.0 / rows), radius_m)
            : landscape::SphereMetric(columns, rows, radius_m);
        const std::size_t count = metric.size();
        if (cube && options.count("--cells"))
        {
            const std::vector<float> bits = readraster(text("--cells", ""), count * mesh::slotcount);  // int32 values
            std::vector<std::int32_t> table(bits.size());
            std::memcpy(table.data(), bits.data(), bits.size() * sizeof(float));
            for (std::size_t c = 0; c < count; c++)
                for (int k = 1; k <= mesh::slotcount; k++)
                {
                    const std::int32_t stored = table[c * mesh::slotcount + static_cast<std::size_t>(k - 1)];
                    if (stored != metric.grid->neighbour(static_cast<mesh::CellId>(c), k))
                        throw std::runtime_error("--cells: neighbour table differs from the cube's at cell " +
                                                 std::to_string(c) + " slot " + std::to_string(k));
                }
        }

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
            const bool window = options.count("--pixel-deg") != 0 || cubewindow;
            seaseeds.assign(count, 0);
            for (std::size_t i = 0; i < count; i++)
            {
                bool ring = false;
                if (cubewindow)  // the window's edge: a cell missing a cardinal slot (the full level has all four)
                {
                    for (int k = 1; k <= mesh::slotcount; k += 2)
                        ring = ring || metric.grid->neighbour(static_cast<mesh::CellId>(i), k) < 0;
                }
                else if (window)
                {
                    const int column = static_cast<int>(i % static_cast<std::size_t>(columns));
                    const int row = static_cast<int>(i / static_cast<std::size_t>(columns));
                    ring = column == 0 || column == columns - 1 || row == 0 || row == rows - 1;
                }
                if (sea[i] && (elevation[i] < lowestlevel || ring))
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
        // Aggrading basins (--aggrade-subsidence 1, roadmap T1 v19): the engine's negative uplift lowered land cells with
        // no floor (20 % of the land subsided, p5 -146 m over the run) and drowned lowlands. With the option, a cell that
        // is land at the step's routing keeps its height where it subsides: the step's uplift there is 0 and the
        // subsidence U dt is counted as basin fill (sediment filling the basin as it sinks; per cell basin_fill_m, the
        // volume basin_fill_km3). Sea cells are unchanged (the kernel reads no uplift on them unless upliftsea is set).
        const bool aggrade = number("--aggrade-subsidence", 0.0) > 0.5;
        std::vector<float> stepuplift(aggrade ? count : 0, 0.0f);
        std::vector<double> basinfill(aggrade ? count : 0, 0.0);
        double basinfilled = 0.0;
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
        // An ice history (--ice-history): each step's ELA offset and glacial ratio are the phase's at its midpoint.
        const bool icehistory = options.count("--ice-history") != 0;
        if (icehistory && !glacial)
            throw std::runtime_error("--ice-history needs --ela");
        const std::vector<IcePhase> icephases = icehistory ? parseicehistory(text("--ice-history", "")) : std::vector<IcePhase>();
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
        landscape::Depressions routeddepressions;  // the last routing's, for --lake-trap
        const auto route = [&](const std::vector<float>& terrain) {
            if (history)
                floodfromseeds(metric, elevation, seaseeds, sea);
            else
                floodsea(metric, elevation, sea);
            landscape::Depressions depressions = landscape::filldepressions(metric, terrain, sea, 0.0f);
            landscape::RoutingOptions routing;
            routing.resolveflats = true;  // flats drain to their low edges, not in straight flood-order lines
            // Randomised receivers (flow_graph.hpp, Rho8 on the sphere's metric): the mean flow follows the gradient
            // instead of the nearest of eight directions, so rivers wander instead of running along the grid's
            // axes and diagonals; a new draw every step (the seed plus the step), deterministic.
            routing.randomised = randomrouting;
            routing.seed = static_cast<std::uint64_t>(routingseed + steps);
            graph = landscape::buildflowgraph(metric, terrain, sea, 0.0f, depressions, routing);
            landscape::accumulateflow(metric, runoff, graph);
            routeddepressions = std::move(depressions);
        };
        // Lakes trap sediment (--lake-trap 1; LakeTrap above), outside the final pulse's no-infill steps (those keep
        // their own rule: no deposition in closed depressions, so the young glacial lakes stay open).
        const bool laketrap = number("--lake-trap", 0.0) > 0.5;
        // Lakes under ice do not trap (--lake-trap-ice-share S, roadmap T1.7 fjords): a depression with at least the
        // share S of its cells under more than --ice-threshold m of ice is left open, the ice carrying its sediment on
        // (a subglacial trough stays a trough; it is drowned at the rise as a fjord). Above 1: every lake traps.
        const double laketrapiceshare = number("--lake-trap-ice-share", 2.0);
        // The trap's fill sloped up from the exits (--lake-fill-gradient, --lake-fill-noise; filloffsets above).
        const double lakefillgradient = number("--lake-fill-gradient", 0.0);
        const std::vector<float> lakefillnoise = options.count("--lake-fill-noise")
            ? readraster(text("--lake-fill-noise", ""), count) : std::vector<float>();
        double laketrapped = 0.0, lakeoverflow = 0.0;
        const std::vector<float> noerosion(laketrap ? count : 0, 0.0f), nouplift(laketrap ? count : 0, 0.0f);
        std::vector<float> erosionafter(laketrap ? count : 0), depositionafter(laketrap ? count : 0);
        const double thresholdslope = number("--threshold-slope", 0.0);
        const float seacap = static_cast<float>(number("--sea-cap", 0.0));
        double eroded = 0.0, deposited = 0.0, marine = 0.0, landslides = 0.0;
        // A final glacial pulse (--no-infill-after, roadmap T1.3): from that time on, the land cells inside a closed
        // depression of the bed (the Priority-Flood surface above the bed, in the step's frame and sea) get no
        // deposition: the step's deposition there is taken off after it, so the ice and rivers keep cutting while
        // the overdeepenings stay open (Earth's glacial lakes are about 15 kyr old, never filled). The removed
        // volume leaves the system (carried out by the ice) and is reported apart.
        const bool noinfill = options.count("--no-infill-after") != 0;
        const double noinfillafter = number("--no-infill-after", 0.0);
        double infillremoved = 0.0;
        // With a history, each step works in the frame of the level at its midpoint: the elevation and the ELA
        // less the level before the step, the level added back after it (the pinned cells are held in the
        // absolute frame, after that).
        std::vector<float> elarelative(history || icehistory ? ela : std::vector<float>());
        const auto shift = [&](float level) {
            for (float& value : elevation)
                value -= level;
        };
        // Per-cell totals for later mineral work (--out-fields PREFIX, raw float32 PREFIX-<name>.f32): the kernel's
        // erosion (fluvial, glacial and hillslope lowering) and deposition (land and marine, less what the final
        // pulse took off closed depressions), m, summed over the steps; the erosion of the steps in which the cell
        // carried more than --ice-threshold m of ice, the time it did (Myr) and its thickest ice (m); the threshold
        // slope limiter's lowering (m); and the lowering (erosion plus landslides) from the last sea-history
        // breakpoint before the end on (the whole run without a history): the rock the final phase exposed.
        const bool fields = options.count("--out-fields") != 0;
        const double icethreshold = number("--ice-threshold", 10.0);
        double laststart = 0.0;
        for (const SeaPoint& point : seahistory)
            if (point.time_myr < duration - 1e-9)
                laststart = std::max(laststart, point.time_myr);
        std::vector<double> erodedtotal, depositedtotal, glacialtotal, icetime, landslidetotal, lasteroded;
        std::vector<float> icemax, cut;
        if (fields)
        {
            for (std::vector<double>* total : {&erodedtotal, &depositedtotal, &glacialtotal, &icetime, &landslidetotal, &lasteroded})
                total->assign(count, 0.0);
            icemax.assign(count, 0.0f);
        }
        double seafractionmin = 1.0;
        const auto seafraction = [&]() {
            std::size_t cells = 0;
            for (const uint8_t value : sea)
                cells += value;
            return count ? static_cast<double>(cells) / static_cast<double>(count) : 0.0;
        };
        // Steps end at the history's breakpoints and at --no-infill-after, so a phase shorter than --step (a
        // 0.25 Myr lowstand) is a step of its own instead of being sampled at a step's midpoint.
        std::vector<double> breakpoints;
        for (const SeaPoint& point : seahistory)
            breakpoints.push_back(point.time_myr);
        if (noinfill)
            breakpoints.push_back(noinfillafter);
        for (const IcePhase& phase : icephases)
            breakpoints.push_back(phase.time_myr);
        const std::vector<std::pair<double, double>> stepafter =
            options.count("--step-after") ? parsestepafter(text("--step-after", "")) : std::vector<std::pair<double, double>>();
        for (const auto& change : stepafter)
            breakpoints.push_back(change.first);
        const auto nextstep = [&](double time) {
            double base = step;
            for (const auto& change : stepafter)
                if (time + 1e-9 >= change.first)
                    base = change.second;
            double length = std::min(base, duration - time);
            for (double point : breakpoints)
                if (point > time + 1e-9)
                    length = std::min(length, point - time);
            return length;
        };
        for (double time = 0.0; time < duration - 1e-9; time += parameters.step_myr)
        {
            parameters.step_myr = nextstep(time);
            const float level = history ? static_cast<float>(sealevel(seahistory, time + 0.5 * parameters.step_myr)) : 0.0f;
            if (history)
                shift(level);
            if (icehistory)
            {
                const IcePhase& phase = icephase(icephases, time + 0.5 * parameters.step_myr);
                glacialparameters.erodibilityratio = phase.ratio;
                const float offset = static_cast<float>(phase.elaoffset_m);
                for (std::size_t i = 0; i < count; i++)
                    elarelative[i] = ela[i] + offset - level;
            }
            else if (history && glacial)
                for (std::size_t i = 0; i < count; i++)
                    elarelative[i] = ela[i] - level;
            for (std::size_t i = 0; i < count; i++)
                surface[i] = elevation[i] + thickness[i];
            route(surface);
            if (history)
                seafractionmin = std::min(seafractionmin, seafraction());
            if (aggrade)
                for (std::size_t i = 0; i < count; i++)
                {
                    const bool filled = !sea[i] && uplift[i] < 0.0f;
                    stepuplift[i] = filled ? 0.0f : uplift[i];
                    if (filled)
                    {
                        const double depth = -static_cast<double>(uplift[i]) * parameters.step_myr;
                        basinfill[i] += depth;
                        basinfilled += depth * metric.area_m2(i) * 1e-9;
                    }
                }
            const std::vector<float>& upliftnow = aggrade ? stepuplift : uplift;
            landscape::GlacialInput ice;
            if (glacial)
            {
                landscape::accumulateice(metric, graph, surface, history || icehistory ? elarelative : ela, snowfall,
                                         glacialparameters.fullicerange_m, iceflux);
                landscape::icethickness(iceflux, glacialparameters, thickness);
                ice.iceflux_m3_yr = &iceflux;
                ice.thickness_m = &thickness;
                ice.erodibilityratio = glacialparameters.erodibilityratio;
                ice.flotationratio = glacialparameters.flotationratio;
            }
            std::vector<uint8_t> closed;  // land cells inside a closed depression of the bed, during the pulse
            if (noinfill && time + 1e-9 >= noinfillafter)
            {
                const landscape::Depressions basins = landscape::filldepressions(metric, elevation, sea, 0.0f);
                closed.assign(count, 0);
                for (std::size_t i = 0; i < count; i++)
                    closed[i] = !sea[i] && basins.filledelevation[i] > elevation[i] ? 1 : 0;
            }
            if (laketrap && closed.empty())
            {
                // The fluvial step with the lakes' exits as roots and no diffusion, the fill, then the diffusion alone.
                LakeTrap trap;
                std::vector<uint8_t> iced;  // depressions under ice (--lake-trap-ice-share): left open
                if (glacial && laketrapiceshare <= 1.0)
                {
                    const std::size_t lakes = routeddepressions.depressions.size();
                    std::vector<double> lakecells(lakes, 0.0), undericecells(lakes, 0.0);
                    for (std::size_t i = 0; i < count; i++)
                    {
                        const int32_t lake = routeddepressions.depressionid[i];
                        if (lake == 0)
                            continue;
                        lakecells[static_cast<std::size_t>(lake - 1)] += 1.0;
                        if (thickness[i] > icethreshold)
                            undericecells[static_cast<std::size_t>(lake - 1)] += 1.0;
                    }
                    iced.assign(lakes, 0);
                    for (std::size_t k = 0; k < lakes; k++)
                        iced[k] = lakecells[k] > 0.0 && undericecells[k] >= laketrapiceshare * lakecells[k] ? 1 : 0;
                }
                traplakes(graph, routeddepressions, trap, iced.empty() ? nullptr : &iced);
                landscape::StreamPowerParameters fluvial = parameters;
                fluvial.hillslopediffusivity_km2_myr = 0.0;
                const landscape::StreamPowerBudget budget = landscape::stepstreampower(
                    metric, graph, erodibility, upliftnow, fluvial, elevation, erosion, deposition, nullptr,
                    glacial ? &ice : nullptr);
                releaselakes(graph, trap);
                filllakes(metric, graph, routeddepressions, thickness, trap, elevation, deposition, lakefillgradient,
                          lakefillnoise);
                const landscape::StreamPowerBudget slopes = landscape::stepstreampower(
                    metric, graph, noerosion, nouplift, parameters, elevation, erosionafter, depositionafter);
                for (std::size_t i = 0; i < count; i++)
                {
                    erosion[i] += erosionafter[i];
                    deposition[i] += depositionafter[i];
                }
                eroded += budget.eroded_km3 + slopes.eroded_km3;
                deposited += budget.landdeposited_km3 + slopes.landdeposited_km3
                    + (trap.placed_m3 - trap.dumped_m3) * 1e-9;
                marine += budget.marinedeposited_km3 + slopes.marinedeposited_km3;
                laketrapped += trap.placed_m3 * 1e-9;
                lakeoverflow += trap.overflow_m3 * 1e-9;
            }
            else
            {
                const landscape::StreamPowerBudget budget = landscape::stepstreampower(
                    metric, graph, erodibility, upliftnow, parameters, elevation, erosion, deposition, nullptr,
                    glacial ? &ice : nullptr);
                eroded += budget.eroded_km3;
                deposited += budget.landdeposited_km3;
                marine += budget.marinedeposited_km3;
            }
            for (std::size_t i = 0; i < closed.size(); i++)
                if (closed[i] && deposition[i] > 0.0f)
                {
                    elevation[i] -= deposition[i];
                    infillremoved += static_cast<double>(deposition[i]) * metric.area_m2(i) * 1e-9;
                }
            if (thresholdslope > 0.0)
                landslides += limitslopes(metric, sea, thresholdslope, elevation, fields ? &cut : nullptr);
            // No sea cell is filled above --sea-cap (m, default the sea level): with the delta top lifted (its
            // fill-and-spill search went quadratic), marine deposits otherwise stacked whole loads on one cell,
            // towers up to 12 km in glacially cut troughs the sea had flooded (F-troughs, G-coast at k 256).
            // With a history the cap is relative to the level (the frame here), seed cells above it included.
            // The deposit the cap takes off is not counted in deposited_m.
            for (std::size_t i = 0; i < count; i++)
                if (sea[i])
                {
                    if (fields && elevation[i] > seacap)
                        depositedtotal[i] -= std::min(static_cast<double>(deposition[i]), static_cast<double>(elevation[i] - seacap));
                    elevation[i] = std::min(elevation[i], seacap);
                }
            if (fields)
            {
                const bool last = time + 1e-9 >= laststart;
                for (std::size_t i = 0; i < count; i++)
                {
                    const bool removed = !closed.empty() && closed[i] && deposition[i] > 0.0f;
                    const double slid = cut.empty() ? 0.0 : static_cast<double>(cut[i]);
                    erodedtotal[i] += erosion[i];
                    depositedtotal[i] += removed ? 0.0 : static_cast<double>(deposition[i]);
                    landslidetotal[i] += slid;
                    if (thickness[i] > icethreshold)
                    {
                        glacialtotal[i] += erosion[i];
                        icetime[i] += parameters.step_myr;
                    }
                    icemax[i] = std::max(icemax[i], thickness[i]);
                    if (last)
                        lasteroded[i] += erosion[i] + slid;
                }
            }
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
        if (fields)
        {
            const std::string prefix = text("--out-fields", "");
            const auto save = [&](const char* name, const std::vector<double>& values) {
                std::vector<float> stored(values.size());
                for (std::size_t i = 0; i < values.size(); i++)
                    stored[i] = static_cast<float>(values[i]);
                writeraster(prefix + "-" + name + ".f32", stored);
            };
            save("eroded_m", erodedtotal);
            save("deposited_m", depositedtotal);
            save("glacial_eroded_m", glacialtotal);
            save("ice_time_myr", icetime);
            writeraster(prefix + "-ice_max_m.f32", icemax);
            save("landslide_m", landslidetotal);
            save("last_eroded_m", lasteroded);
            if (aggrade)
                save("basin_fill_m", basinfill);
        }
        std::printf("steps %d eroded_km3 %.1f land_deposited_km3 %.1f marine_deposited_km3 %.1f landslides_km3 %.1f", steps,
                    eroded, deposited, marine, landslides);
        if (history)  // the schedule's lowest and final levels, the sea's least and final extent (cell fractions)
            std::printf(" sea_min_m %.1f sea_final_m %.1f sea_fraction_min %.4f sea_fraction_final %.4f", lowestlevel,
                        static_cast<double>(finallevel), seafractionmin, seafraction());
        if (noinfill)  // the deposition taken off closed depressions during the final pulse
            std::printf(" no_infill_after_myr %.2f infill_removed_km3 %.1f", noinfillafter, infillremoved);
        if (laketrap)  // the sediment the lakes kept, and what passed out of the land over full lakes
            std::printf(" lake_trapped_km3 %.1f lake_overflow_km3 %.1f", laketrapped, lakeoverflow);
        if (aggrade)  // the subsidence of land cells filled by sediment instead of lowering them
            std::printf(" basin_fill_km3 %.1f", basinfilled);
        std::printf("\n");
    }
    catch (const std::exception& error)
    {
        std::cerr << "terrain_sandbox: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
