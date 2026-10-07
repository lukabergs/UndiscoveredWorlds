#include "simulations/hydrology/rivers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <vector>

#include "core/climate_seasons.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/landscape/distributaries.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "world/planet.hpp"

namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr double pi = 3.14159265358979323846;
// Legacy flow law (erosion plan section 6, "Legacy flow units"): a cell of mean size
// under 1 mm/month contributes 1 / river_factor units; 0.012 m/yr is 1 mm/month.
constexpr double metresperyearpermmpermonth = 0.012;
constexpr int januaryseason = seasonjanuary;
constexpr int julyseason = seasonjuly;

// Seasonal runoff of one cell in mm/month (section 6, "Seasonal flows"). The snow
// rule (B1): the precipitation of a season at or above the snow temperature, plus
// what the cold seasons before it held as snow (cycling over the four seasons); a
// cell with no warm season releases its annual depth evenly (annual / 4 a season),
// so the seasons sum to the annual depth everywhere. With landscape.runoff.budyko on, the
// snow-rule seasons are scaled so that they sum to the routing's annual runoff
// depth (hydrosphere.runoff), so annual and seasonal flows rest on the same depth;
// off, they stay the precipitation, as in Phase B. The winter baseflow (erosion
// plan C2): a season's runoff = (1 - b) x the snow-rule runoff + b x the annual
// depth spread evenly, b = landscape.rivers.baseflow_fraction (b = 0 keeps B1 on
// cells with a warm season).
void seasonalrunoff(const planet& world, const worldgrid<float>* annualrunoff, int x, int y,
                    double runoff[CLIMATESEASONCOUNT])
{
    const double baseflow = tuning::landscape::rivers::baseflowFraction;
    double rain[CLIMATESEASONCOUNT];
    bool warm[CLIMATESEASONCOUNT];
    int firstwarm = -1;
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        rain[season] = std::max(0.0, static_cast<double>(world.seasonalrainfloat(season, x, y)));
        warm[season] = world.seasonaltempfloat(season, x, y) >= tuning::landscape::runoff::snowTemperatureC;
        runoff[season] = 0.0;
        if (warm[season] && firstwarm < 0)
            firstwarm = season;
    }
    double precipitation = 0.0;  // mm/month summed over the seasons
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        precipitation += rain[season];
    // The annual depth in the units of the seasons (mm/month summed over the four
    // seasons): the routing's depth (mm/yr, three months a season) or the precipitation.
    const double annual = annualrunoff != nullptr
        ? std::max(0.0, static_cast<double>((*annualrunoff)[x][y])) / 3.0
        : precipitation;
    if (firstwarm >= 0)
    {
        double held = 0.0;
        for (int step = 1; step <= CLIMATESEASONCOUNT; step++)
        {
            const int season = (firstwarm + step) % CLIMATESEASONCOUNT;
            if (warm[season])
            {
                runoff[season] = rain[season] + held;
                held = 0.0;
            }
            else
                held += rain[season];
        }
        // Every season's precipitation reaches a warm season, so the snow rule sums to the precipitation.
        const double scale = annualrunoff != nullptr ? (precipitation > 0.0 ? annual / precipitation : 0.0) : 1.0;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            runoff[season] *= scale;
    }
    else
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            runoff[season] = annual / CLIMATESEASONCOUNT;
    if (baseflow > 0.0)
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            runoff[season] = (1.0 - baseflow) * runoff[season] + baseflow * annual / CLIMATESEASONCOUNT;
}

// Distributaries on prograded land (erosion plan C5b), written into the legacy delta
// fields by landscape::splitdistributaries (distributaries.hpp: the split, its order
// and the layout the delta fields' readers expect) on the routing's filled surface,
// with the rivers' January and July flows and landscape.rivers.minimum_flow. Every
// cell's delta fields are written: 0 except on the delta cells of prograded land.
void writedistributaries(planet& world, const landscape::SphereMetric& metric, const std::vector<std::uint8_t>& land)
{
    const planet& view = world;
    const std::size_t size = metric.size();
    const worldgrid<std::uint8_t>& prograded = view.prograded();
    const worldgrid<std::uint8_t>& direction = view.flowdirection();
    const worldgrid<float>& filled = view.filledelevation();

    std::vector<std::uint8_t> sea(size), progradedcells(size), codes(size);
    std::vector<float> surface(size);
    std::vector<std::int32_t> january(size), july(size);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int x = metric.column(cell);
        const int y = metric.row(cell);
        sea[cell] = land[cell] ? 0 : 1;
        progradedcells[cell] = prograded[x][y];
        codes[cell] = direction[x][y];
        surface[cell] = filled[x][y];
        january[cell] = view.riverjan(x, y);
        july[cell] = view.riverjul(x, y);
    }
    const landscape::Distributaries split = landscape::splitdistributaries(
        metric, sea, progradedcells, surface, static_cast<double>(view.sealevel()), codes, january, july,
        static_cast<double>(tuning::landscape::rivers::minimumFlow));
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const int x = metric.column(cell);
        const int y = metric.row(cell);
        world.setdeltadir(x, y, split.direction[cell]);
        world.setdeltajan(x, y, split.january[cell]);
        world.setdeltajul(x, y, split.july[cell]);
    }
}
}

void applyriverstage(planet& world)
{
    const planet& view = world;
    const int columns = view.width() + 1;
    const int rows = view.height() + 1;
    const double radius = static_cast<double>(view.radius());
    const landscape::SphereMetric metric(columns, rows, radius);
    const std::size_t size = metric.size();
    const worldgrid<std::uint8_t>& direction = view.flowdirection();
    const worldgrid<float>& discharge = view.discharge();

    // Land cells and their receivers along flow_direction (none into the sea or past a pole). The
    // sea is the routing's (flowrouting::oceancells, erosion plan H3): an enclosed below-sea basin is
    // land here, so rivers run into its lake or pan; the distributaries put the sea level only on
    // these sea cells and read a basin cell at its filled surface, the lake level.
    const std::vector<std::uint8_t> ocean = flowrouting::oceancells(view);
    std::vector<std::uint8_t> land(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
        land[cell] = ocean[cell] != 0 ? 0 : 1;

    constexpr std::size_t none = static_cast<std::size_t>(-1);
    std::vector<std::size_t> receiver(size, none);
    std::vector<int> indegree(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (!land[cell])
            continue;
        const int code = direction[metric.column(cell)][metric.row(cell)];
        if (code < 1 || code > landscape::directioncount)
            continue;
        const std::size_t next = metric.neighbour(cell, code);
        if (next == cell || !land[next])
            continue;
        receiver[cell] = next;
        indegree[next]++;
    }

    // Seasonal runoff x cell area, accumulated downstream (Kahn's order, increasing cell index).
    const worldgrid<float>* annualrunoff = tuning::landscape::runoff::budyko ? &view.runoff() : nullptr;
    std::vector<double> accumulated(size * CLIMATESEASONCOUNT, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (!land[cell])
            continue;
        double runoff[CLIMATESEASONCOUNT];
        seasonalrunoff(view, annualrunoff, metric.column(cell), metric.row(cell), runoff);
        const double area = metric.area_m2(cell);
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            accumulated[cell * CLIMATESEASONCOUNT + season] = runoff[season] * area;
    }
    std::deque<std::size_t> queue;
    for (std::size_t cell = 0; cell < size; cell++)
        if (land[cell] && indegree[cell] == 0)
            queue.push_back(cell);
    while (!queue.empty())
    {
        const std::size_t cell = queue.front();
        queue.pop_front();
        const std::size_t next = receiver[cell];
        if (next == none)
            continue;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            accumulated[next * CLIMATESEASONCOUNT + season] += accumulated[cell * CLIMATESEASONCOUNT + season];
        if (--indegree[next] == 0)
            queue.push_back(next);
    }

    // Annual flow in legacy units, split into January and July by the seasonal runoff. The
    // seasonal shares come from the lossless accumulation above, while the annual flow is the
    // routing's discharge after its transmission losses: the losses are spread evenly over the
    // seasons, in proportion to each season's share.
    // runoff.coefficient is the unit of the conversion, not the runoff law: a legacy
    // flow unit is the discharge of a mean-sized cell under 1 mm/month at the Phase A
    // coefficient (0.4) divided by river_factor, so the units, and rivers.minimum_flow
    // with them, stay those of Phase B whichever law (Budyko or the coefficient) made
    // the discharge.
    const double meanarea = 4.0 * pi * radius * radius / (static_cast<double>(columns) * static_cast<double>(rows));
    const double unitdischarge = tuning::landscape::runoff::coefficient * metresperyearpermmpermonth *
        static_cast<double>(view.riverfactor()) * meanarea;
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            if (!land[cell])
            {
                world.setriverdir(x, y, 0);
                world.setriverjan(x, y, 0);
                world.setriverjul(x, y, 0);
                continue;
            }
            world.setriverdir(x, y, direction[x][y]);
            const double flow = static_cast<double>(discharge[x][y]) * secondsperyear / unitdischarge;
            int january = 0;
            int july = 0;
            if (flow >= static_cast<double>(tuning::landscape::rivers::minimumFlow))
            {
                double total = 0.0;
                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                    total += accumulated[cell * CLIMATESEASONCOUNT + season];
                const double januaryshare = total > 0.0 ? 4.0 * accumulated[cell * CLIMATESEASONCOUNT + januaryseason] / total : 1.0;
                const double julyshare = total > 0.0 ? 4.0 * accumulated[cell * CLIMATESEASONCOUNT + julyseason] / total : 1.0;
                january = static_cast<int>(std::lround(flow * januaryshare));
                july = static_cast<int>(std::lround(flow * julyshare));
            }
            world.setriverjan(x, y, january);
            world.setriverjul(x, y, july);
        }
    }

    // River mouths: as the legacy map does, the sea cell an outlet drains into takes the
    // outlet's direction and flows (the delta stage reads them there); where several
    // outlets share a sea cell, their flows add and the direction is the largest one's.
    std::vector<int> mouthbest(size, -1);
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            const int code = direction[x][y];
            if (!land[cell] || code < 1 || code > landscape::directioncount)
                continue;
            const std::size_t mouth = metric.neighbour(cell, code);
            if (mouth == cell || land[mouth])
                continue;
            const int mx = metric.column(mouth);
            const int my = metric.row(mouth);
            const int january = view.riverjan(x, y);
            const int july = view.riverjul(x, y);
            if (january == 0 && july == 0)
                continue;
            const int flow = january + july;
            if (flow > mouthbest[mouth])
            {
                mouthbest[mouth] = flow;
                world.setriverdir(mx, my, code);
            }
            world.setriverjan(mx, my, view.riverjan(mx, my) + january);
            world.setriverjul(mx, my, view.riverjul(mx, my) + july);
        }
    }
    if (view.erosionduration() > 0.0f)
        writedistributaries(world, metric, land);
    world.setmaxriverflow();
}
