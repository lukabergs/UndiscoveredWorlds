#include "simulations/landscape/regional_pass.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_water.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan step E3 (docs/architecture/erosion-hydrology-plan.md, Phase E, and
// section 6, "Decided for Phase E"); Phase R step R0 split the start surface and the
// sea (regional_surface.cpp) and the lakes and ice (regional_water.cpp) out of it.

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

RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters)
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

    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const double sealevel = static_cast<double>(inputs.sealevel);
    const std::int64_t planetfinewidth = static_cast<std::int64_t>(width) * window.scale;

    // 1. The warp and the start surface (regional_surface.hpp).
    const FineWarp warp = regionalwarp(window, inputs.seed, parameters);
    const StartSurface start = regionalstartsurface(inputs, window, warp, parameters);
    const std::vector<float>& base = start.lapsebase_m;
    std::vector<float> surface = start.start_m;

    RegionalFields fields;
    fields.window = window;
    fields.blockmeanerror_m = start.blockmeanerror_m;
    fields.sharepasses = start.sharepasses;
    fields.sharemisses = start.sharemisses;

    // 2. The sea by local rules on the start surface (regional_surface.hpp), before the
    // channels (R1), so the channels end at the fine coast and never change the mask.
    const std::vector<uint8_t> sea = regionalsea(inputs, window, warp, parameters, surface);

    // 3. The pinned channels: they end at their first fine sea cell, and their beds
    // follow the start surface's running minimum (channel_network.hpp). Beds are
    // pinned on land; a link's last cell in the sea keeps its own (seabed) elevation,
    // since the solve pins every sea cell at the sea level anyway.
    const FineChannels channels =
        rasterisechannels(window, inputs.network, parameters.channels, inputs.sealevel, surface, sea);
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
    const AnalyticalReliefResult solved = solveanalyticalrelief(metric, solveinput, parameters.relief);
    fields.reliefconverged = solved.converged;
    fields.changedreceivers = solved.changedreceivers;
    std::vector<float> elevation = surface;
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
        fields.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        fields.maxlowering_m = largest;
    }

    // 5. Lakes and ice (regional_water.hpp).
    const FineWater water = regionalwater(inputs, window, warp, parameters, sea, channels, elevation);

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
                discharge[cell] = static_cast<float>(solved.discharge_m3_yr[cell] / secondsperyear);
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
    fields.sea = cut(sea, window);
    fields.riverdirection = cut(riverdirection, window);
    fields.discharge_m3s = cut(discharge, window);
    fields.januaryflow = cut(januaryflow, window);
    fields.julyflow = cut(julyflow, window);
    fields.lakesurface_m = cut(water.lakesurface_m, window);
    fields.januarytemperature_c = cut(januarytemp, window);
    fields.julytemperature_c = cut(julytemp, window);
    fields.januaryprecipitation_mm = cut(januaryrain, window);
    fields.julyprecipitation_mm = cut(julyrain, window);
    fields.icethickness_m = cut(water.icethickness_m, window);
    fields.climate = cut(climate, window);
    return fields;
}
}
