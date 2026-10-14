// terrain_sandbox: the physics of the terrain sandbox (scripts/sandbox/terrain_sandbox.py). Routes and
// erodes a raster the sandbox writes, with the pipeline's own kernels (Priority-Flood, D8 flow graph,
// implicit stream power with deposition and hillslope diffusion), so procedural terrain can be judged
// after real drainage and erosion without running the pipeline. Not a pipeline stage: no world, no
// bundle, no tuning header. Files are raw little-endian float32, rows north first; elevations are metres
// relative to today's sea level (--sea-history moves the sea during the run).
#include <algorithm>
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
// releaselakes). The order stays valid: a root needs no receiver before it.
void traplakes(landscape::FlowGraph& graph, const landscape::Depressions& depressions, LakeTrap& trap)
{
    const std::vector<int32_t>& id = depressions.depressionid;
    trap.spill.assign(graph.receiver.size(), -1);
    trap.exits.clear();
    for (int32_t entry : graph.order)
    {
        const std::size_t i = static_cast<std::size_t>(entry);
        const int32_t lake = id[i];
        if (lake == 0 || depressions.depressions[static_cast<std::size_t>(lake - 1)].terminal)
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

// After the fluvial step on the trapped graph (graph restored): takes the loads off the exits and fills the lakes.
void filllakes(const landscape::SphereMetric& metric, const landscape::FlowGraph& graph,
               const landscape::Depressions& depressions, const std::vector<float>& thickness, LakeTrap& trap,
               std::vector<float>& elevation, std::vector<float>& deposition)
{
    const std::vector<int32_t>& id = depressions.depressionid;
    const std::size_t lakes = depressions.depressions.size();
    if (trap.exits.empty())
        return;
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
    // A level fill of `volume` (m^3) into lake `lake` up to its spill level; returns what did not fit.
    std::vector<std::pair<double, std::size_t>> floor;
    const auto fill = [&](std::size_t lake, double volume) {
        if (volume <= 0.0)
            return 0.0;
        const double level = static_cast<double>(depressions.depressions[lake - 1].spillelevation_m);
        floor.clear();
        for (std::size_t i : cells[lake])
        {
            const double surface = static_cast<double>(elevation[i]) + static_cast<double>(thickness[i]);
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
            const float bed = static_cast<float>(height - static_cast<double>(thickness[i]));
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
                 "  [--no-infill-after MYR] (from then on no deposition inside closed depressions of the bed)\n"
                 "  [--lake-trap 0|1] (lakes keep the sediment reaching them, filled up to their spill; not in no-infill steps)\n"
                 "  [--ice-history \"t:ela_offset_m:ratio,...\"] (Myr, m added to the ELA, glacial ratio: piecewise constant)\n"
                 "  [--step-after \"T:STEP,...\"] (from T Myr on, steps of STEP Myr)\n"
                 "  [--out-fields PREFIX] (per-cell totals: PREFIX-eroded_m.f32, -deposited_m, -glacial_eroded_m,\n"
                 "   -ice_time_myr, -ice_max_m, -landslide_m, -last_eroded_m) [--ice-threshold M 10]\n"
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
                traplakes(graph, routeddepressions, trap);
                landscape::StreamPowerParameters fluvial = parameters;
                fluvial.hillslopediffusivity_km2_myr = 0.0;
                const landscape::StreamPowerBudget budget = landscape::stepstreampower(
                    metric, graph, erodibility, uplift, fluvial, elevation, erosion, deposition, nullptr,
                    glacial ? &ice : nullptr);
                releaselakes(graph, trap);
                filllakes(metric, graph, routeddepressions, thickness, trap, elevation, deposition);
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
                    metric, graph, erodibility, uplift, parameters, elevation, erosion, deposition, nullptr,
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
        std::printf("\n");
    }
    catch (const std::exception& error)
    {
        std::cerr << "terrain_sandbox: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
