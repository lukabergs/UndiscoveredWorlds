#include "simulations/landscape/regional_pass.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include "simulations/landscape/regional_cirques.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_ice.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_troughs.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/regional_water.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan step E3 (docs/architecture/erosion-hydrology-plan.md, Phase E, and
// section 6, "Decided for Phase E"); Phase R step R0 split the start surface and the
// sea (regional_surface.cpp) and the lakes and ice (regional_water.cpp) out of it;
// Phase R2 step R6 wired the incision and the valleys (regional_valleys.cpp), which R9
// implemented with the solve's uplift term; Phase E5 step E5.0 wired the glacial forms
// (regional_glacial.cpp, regional_troughs.cpp, regional_cirques.cpp, regional_ice.cpp),
// inert at the defaults. Multiscale plan step S2.3 split the pass into the parts the
// refinement loop calls (regionallevelsteps, finishregionallevel), unchanged in order and
// arithmetic.

namespace landscape
{
namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

template<typename T>
void checksize(const std::vector<T>& field, std::size_t size, const char* name)
{
    if (field.size() != size)
        throw std::invalid_argument(std::string("runregionalpass: ") + name + " must hold width x height values");
}

// The Phase R2 fields may be absent (empty: zero everywhere).
template<typename T>
void checkoptional(const std::vector<T>& field, std::size_t size, const char* name)
{
    if (!field.empty() && field.size() != size)
        throw std::invalid_argument(std::string("runregionalpass: ") + name + " must be empty or hold width x height values");
}

// The domain's values at the kept window's cells.
template<typename T>
std::vector<T> cut(const std::vector<T>& domain, const RegionalWindow& window)
{
    const int width = window.finewidth(), height = window.fineheight();
    const int domainwidth = window.domainfinewidth();
    const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
    std::vector<T> kept(static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            kept[at(x, y, width)] = domain[at(x0 + x, y0 + y, domainwidth)];
    return kept;
}
}

void checkregionalinputs(const RegionalGlobalInputs& inputs, const RegionalWindow& window)
{
    const int width = inputs.width, height = inputs.height;
    if (width < 1 || height < 1 || window.worldwidth != width || window.worldheight != height || window.scale < 1)
        throw std::invalid_argument("runregionalpass: the window belongs to another grid");
    if (!(inputs.radius_m > 0.0) || !std::isfinite(inputs.sealevel))
        throw std::invalid_argument("runregionalpass: the radius must be positive and the sea level finite");
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    checksize(inputs.elevation_m, globalsize, "elevation_m");
    checksize(inputs.network.sea, globalsize, "network.sea");
    checksize(inputs.runoff_mm_yr, globalsize, "runoff_mm_yr");
    checksize(inputs.lakesurface_m, globalsize, "lakesurface_m");
    checksize(inputs.januarytemperature_c, globalsize, "januarytemperature_c");
    checksize(inputs.julytemperature_c, globalsize, "julytemperature_c");
    checksize(inputs.januaryprecipitation_mm, globalsize, "januaryprecipitation_mm");
    checksize(inputs.julyprecipitation_mm, globalsize, "julyprecipitation_mm");
    checksize(inputs.icethickness_m, globalsize, "icethickness_m");
    checksize(inputs.climate, globalsize, "climate");
    checkoptional(inputs.upliftrate_m_myr, globalsize, "upliftrate_m_myr");
    checkoptional(inputs.glaciation, globalsize, "glaciation");
    checkoptional(inputs.denudation_m, globalsize, "denudation_m");
    checkoptional(inputs.island, globalsize, "island");
    checkoptional(inputs.mountainisland, globalsize, "mountainisland");
    checkoptional(inputs.equilibriumline_m, globalsize, "equilibriumline_m");
}

RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters)
{
    checkregionalinputs(inputs, window);

    // The block-mean correction reaches 2 global cells per pass beyond any fine cell it
    // changes, so the start surface and the fine sea agree between windows from
    // 2 x meanpasses global cells inside the domain's edge. R1: the channels read them
    // along their paths, which reach one global cell nearer the edge
    // (channel_network.hpp), so the channels agree from 2 x meanpasses + 1 cells; the
    // two cells kept here leave the solve one global cell between a kept cell and the
    // pinned cells it depends on (Phase E allowed it two). As in Phase E, the overlap
    // tests at this smallest halo (k 4, 5 and 8, low relief) are what establish it.
    if (parameters.meanpasses < 0 || window.halo < 2 * parameters.meanpasses + 2)
        throw std::invalid_argument("runregionalpass: the halo must be at least 2 x meanpasses + 2 global cells");

    // 1. The warp and the start surface (regional_surface.hpp).
    const FineWarp warp = regionalwarp(window, inputs.seed, parameters);
    StartSurface start = regionalstartsurface(inputs, window, warp, parameters);
    const GlobalIncision incision = globalincision(inputs, parameters);
    RegionalLevel level = regionallevelsteps(inputs, window, parameters, warp, std::move(start), incision, RegionalDetails{});
    return finishregionallevel(inputs, window, parameters, level);
}

RegionalLevel regionallevelsteps(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const RegionalParameters& parameters, const FineWarp& warp, StartSurface start,
                                 const GlobalIncision& incision, const RegionalDetails& details)
{
    const int width = inputs.width, height = inputs.height;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const std::int64_t planetfinewidth = static_cast<std::int64_t>(width) * window.scale;
    if (start.start_m.size() != count || start.lapsebase_m.size() != count)
        throw std::invalid_argument("regionallevelsteps: the start surface must cover the window's domain");

    RegionalLevel level;
    level.blockmeanerror_m = start.blockmeanerror_m;
    level.sharepasses = start.sharepasses;
    level.sharemisses = start.sharemisses;
    level.lapsebase_m = std::move(start.lapsebase_m);
    std::vector<float> surface = std::move(start.start_m);

    // 2. The sea by local rules on the start surface (regional_surface.hpp), before the
    // channels (R1), so the channels end at the fine coast and never change the mask.
    level.sea = regionalsea(inputs, window, warp, parameters, surface);
    const std::vector<uint8_t>& sea = level.sea;

    // 3. The pinned channels: they end at their first fine sea cell, and their beds
    // follow the start surface's running minimum (channel_network.hpp). Beds are
    // pinned on land; a link's last cell in the sea keeps its own (seabed) elevation,
    // since the solve pins every sea cell at the sea level anyway. R9 (regional_valleys.hpp):
    // with the incision active the bed law reads the network's incised filled elevations
    // and the start surface lowered by the upsampled incision, so the channels cut below
    // the land by the world's uplift over the elapsed time; inactive, the inputs as they are.
    GlobalNetwork incisednetwork;
    std::vector<float> bedsurface;
    if (incision.active)
    {
        incisednetwork = inputs.network;
        incisednetwork.filledelevation_m = incision.filledelevation_m;
        const std::vector<float> fineincision = upsample(incision.incision_m, width, height, window);
        bedsurface = surface;
        for (std::size_t cell = 0; cell < count; cell++)
            bedsurface[cell] = static_cast<float>(static_cast<double>(surface[cell]) - std::max(0.0f, fineincision[cell]));
    }
    level.channels = rasterisechannels(window, incision.active ? incisednetwork : inputs.network,
                                       parameters.channels, inputs.sealevel,
                                       incision.active ? bedsurface : surface, sea);
    const FineChannels& channels = level.channels;
    for (std::size_t cell = 0; cell < count; cell++)
        if (channels.channel[cell] != 0 && sea[cell] == 0)
            surface[cell] = channels.bed_m[cell];

    // 4. The analytical relief: channels, sea and the domain's outer ring pinned. A
    // sea cell enters the solve at the sea level (the base level of the rivers), so
    // the solve never lowers land below it through the sea; it keeps its own elevation.
    const PatchMetric metric = finepatchmetric(window, inputs.radius_m);
    const std::vector<float> runoff = upsample(inputs.runoff_mm_yr, width, height, window);
    AnalyticalReliefInput solveinput;
    solveinput.start_m = surface;
    solveinput.pinned.assign(count, 0);
    solveinput.runoff_m_yr.resize(count);
    solveinput.key.resize(count);
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const bool ring = x == 0 || y == 0 || x == domainwidth - 1 || y == domainheight - 1;
            solveinput.pinned[cell] = channels.channel[cell] != 0 || sea[cell] != 0 || ring ? 1 : 0;
            if (sea[cell] != 0)
                solveinput.start_m[cell] = inputs.sealevel;
            solveinput.runoff_m_yr[cell] = std::max(0.0f, runoff[cell]) / 1000.0f;
            solveinput.key[cell] = static_cast<std::uint64_t>(window.absolutefinecolumn(x) +
                                                              window.absolutefinerow(y) * planetfinewidth);
        }
    // R9: the solve's uplift term from the incision's rate (analytical_relief.hpp), only
    // with the incision active; empty otherwise (the solve without uplift exactly).
    if (incision.active)
    {
        solveinput.uplift_m_myr = upsample(incision.rate_m_myr, width, height, window);
        for (float& rate : solveinput.uplift_m_myr)
            rate = std::max(0.0f, rate);
    }
    AnalyticalReliefResult solved = solveanalyticalrelief(metric, solveinput, parameters.relief);
    level.reliefconverged = solved.converged;
    level.changedreceivers = solved.changedreceivers;
    std::vector<float>& elevation = level.elevation_m;
    elevation = surface;
    for (std::size_t cell = 0; cell < count; cell++)
        if (solveinput.pinned[cell] == 0)
            elevation[cell] = solved.elevation_m[cell];
    {
        // The lowering by the solve over the kept window's unpinned (so land) cells.
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double total = 0.0, largest = 0.0;
        std::size_t cells = 0;
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (solveinput.pinned[cell] != 0)
                    continue;
                const double lowering = static_cast<double>(surface[cell]) - static_cast<double>(elevation[cell]);
                total += lowering;
                largest = std::max(largest, lowering);
                cells++;
            }
        level.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        level.maxlowering_m = largest;
    }
    level.solveddischarge_m3_yr = std::move(solved.discharge_m3_yr);

    // 4b. The valley floors around the channels (regional_valleys.hpp; R9, a no-op while
    // valley_floor_km is 0), after the solve's diagnostics and before the lakes, which
    // read the result.
    const bool band = details.factor >= 2;
    const double bandlower = band ? details.detailcells / 2.0 : 0.0;
    const double bandupper = band ? details.detailcells * details.factor / 2.0 : 0.0;
    if (details.valleys)
        carvevalleys(inputs, window, parameters, channels, sea, elevation, bandlower, bandupper);

    // 4c. The glacial forms (Phase E5): the fine ELA and glaciation fields (inactive, and
    // no field read, unless a glacial feature is on: regional_glacial.hpp), then the
    // troughs (E5.1) and the cirques (E5.2), each a no-op at its defaults.
    level.glacial = regionalglacialfields(inputs, window, parameters);
    if (details.troughs)
        carvetroughs(inputs, window, parameters, level.glacial, channels, sea, elevation, bandlower, bandupper);
    if (details.cirques)
        carvecirques(inputs, window, parameters, level.glacial, channels, sea, elevation, bandlower, bandupper);

    // 5. Lakes and ice (regional_water.hpp).
    level.water = regionalwater(inputs, window, warp, parameters, sea, channels, elevation, details.globallakedepth);
    return level;
}

RegionalFields finishregionallevel(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                   const RegionalParameters& parameters, RegionalLevel& level)
{
    const int width = inputs.width, height = inputs.height;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const double sealevel = static_cast<double>(inputs.sealevel);
    if (level.elevation_m.size() != count || level.lapsebase_m.size() != count)
        throw std::invalid_argument("finishregionallevel: the level's fields must cover the window's domain");
    const std::vector<float>& base = level.lapsebase_m;
    const std::vector<float>& elevation = level.elevation_m;
    const FineChannels& channels = level.channels;

    // 5b. Ice above the present ELA (E5.3, a no-op at its defaults), added to the ice of 5.
    iceabovetheela(inputs, window, parameters, level.glacial, level.sea, elevation, level.water.icethickness_m);

    RegionalFields fields;
    fields.window = window;
    fields.blockmeanerror_m = level.blockmeanerror_m;
    fields.sharepasses = level.sharepasses;
    fields.sharemisses = level.sharemisses;
    fields.reliefconverged = level.reliefconverged;
    fields.changedreceivers = level.changedreceivers;
    fields.meanlowering_m = level.meanlowering_m;
    fields.maxlowering_m = level.maxlowering_m;

    // 6.-7. Rivers and climate on the fine grid.
    const std::vector<float> januarytemperature = upsample(inputs.januarytemperature_c, width, height, window);
    const std::vector<float> julytemperature = upsample(inputs.julytemperature_c, width, height, window);
    const std::vector<float> januaryprecipitation = upsample(inputs.januaryprecipitation_mm, width, height, window);
    const std::vector<float> julyprecipitation = upsample(inputs.julyprecipitation_mm, width, height, window);
    const double lapse = static_cast<double>(inputs.lapserate_c_per_km) / 1000.0;

    std::vector<float> discharge(count, 0.0f);
    std::vector<uint8_t> riverdirection(count, 0);
    std::vector<int32_t> januaryflow(count, 0), julyflow(count, 0);
    std::vector<float> januarytemp(count), julytemp(count), januaryrain(count), julyrain(count);
    std::vector<int16_t> climate(count);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            if (channels.channel[cell] != 0)
            {
                riverdirection[cell] = channels.direction[cell];
                discharge[cell] = channels.discharge_m3s[cell];
                januaryflow[cell] = channels.januaryflow[cell];
                julyflow[cell] = channels.julyflow[cell];
            }
            else
                discharge[cell] = static_cast<float>(level.solveddischarge_m3_yr[cell] / secondsperyear);
            // Heights above the sea level only (the climate's rule, climate_elevation.hpp):
            // a sea cell whose interpolated elevation is below the sea level keeps the
            // interpolated temperature exactly.
            const double lift = std::max(0.0, static_cast<double>(base[cell]) - sealevel) -
                                std::max(0.0, static_cast<double>(elevation[cell]) - sealevel);
            januarytemp[cell] = static_cast<float>(static_cast<double>(januarytemperature[cell]) + lapse * lift);
            julytemp[cell] = static_cast<float>(static_cast<double>(julytemperature[cell]) + lapse * lift);
            januaryrain[cell] = std::max(0.0f, januaryprecipitation[cell]);
            julyrain[cell] = std::max(0.0f, julyprecipitation[cell]);
            // The global cell's class: classifyMonthly needs twelve months, the world
            // holds two seasons.
            climate[cell] = inputs.climate[global];
        }
    }

    // 8. The kept window.
    fields.elevation_m = cut(elevation, window);
    fields.sea = cut(level.sea, window);
    fields.riverdirection = cut(riverdirection, window);
    fields.discharge_m3s = cut(discharge, window);
    fields.januaryflow = cut(januaryflow, window);
    fields.julyflow = cut(julyflow, window);
    fields.lakesurface_m = cut(level.water.lakesurface_m, window);
    fields.januarytemperature_c = cut(januarytemp, window);
    fields.julytemperature_c = cut(julytemp, window);
    fields.januaryprecipitation_mm = cut(januaryrain, window);
    fields.julyprecipitation_mm = cut(julyrain, window);
    fields.icethickness_m = cut(level.water.icethickness_m, window);
    fields.climate = cut(climate, window);
    return fields;
}
}
