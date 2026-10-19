#include "simulations/landscape/regional_pass.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/refinement_projection.hpp"
#include "simulations/landscape/regional_folds.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_ice.hpp"
#include "simulations/landscape/regional_solve.hpp"
#include "simulations/landscape/regional_surface.hpp"
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
// arithmetic; S3.0 moved regionallevelsteps' parts into one function per LevelStep in its
// unit (the hook interface, regional_pass.hpp), unchanged likewise.

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
    checkoptional(inputs.reliefamplitude_m, globalsize, "reliefamplitude_m");
    checkoptional(inputs.valleyspacing_km, globalsize, "valleyspacing_km");
    checkoptional(inputs.foldstrain, globalsize, "foldstrain");
    checkoptional(inputs.folddirection_rad, globalsize, "folddirection_rad");
    checkoptional(inputs.sediment_m, globalsize, "sediment_m");
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
    // Multiscale plan S3.3: with tributaries on, tributaryreach more (regionalhalo).
    if (parameters.meanpasses < 0 || window.halo < regionalhalo(parameters.meanpasses, parameters.tributaries.enabled != 0))
        throw std::invalid_argument("runregionalpass: the halo must be at least 2 x meanpasses + 3 global cells (plus "
                                    "tributaryreach with tributaries on)");

    // 1. The warp (regional_surface.hpp) and the incision (regional_valleys.hpp); the start
    // surface is the noise step's (regionalstartsurface).
    const FineWarp warp = regionalwarp(window, inputs.seed, parameters);
    const GlobalIncision incision = globalincision(inputs, parameters);
    const RegionalDetails details{};
    const std::shared_ptr<const mesh::CellGrid> grid = regionalwindowgrid(window, inputs.radius_m);
    const LevelContext context{inputs, window, *grid, parameters, warp, incision, details};
    LevelState level;
    regionallevelsteps(context, level);
    return finishregionallevel(inputs, window, parameters, level);
}

void regionallevelsteps(const LevelContext& context, LevelState& state)
{
    // The steps in LevelStep order (refinement_levels.hpp); each unit's function is the code
    // the pass ran here up to S2.8, unchanged in its operations and their order (S3.0).
    noiselevel(context, state);         // the start surface (level 1, the single pass) or the octave band
    presmoothinglevel(context, state);  // the parent-mean correction; hands the start surface on
    coastlevel(context, state);         // 2. the sea by local rules, before the channels (R1)
    incisionlevel(context, state);      // 3. the incision's bed surface (R9)
    channelslevel(context, state);      // 3. the pinned channels and their tributaries
    solvelevel(context, state);         // 4. the analytical relief
    foldfabriclevel(context, state);    // 4a. the fold fabric (S5.2), cut by the forms after it
    valleyfloorslevel(context, state);  // 4b. the valley floors
    glacialformslevel(context, state);  // 4c. the glacial troughs and cirques
    unfoldforlakes(context, state);     // S5.3b: the lakes decide on the surface without the fabric,
    lakeslevel(context, state);         // 5. lakes and ice
    refoldafterlakes(context, state);   // and the fabric goes back on the dry land they left
    smoothinglevel(context, state);     // the refinement loop only
    projectionlevel(context, state);    // the refinement loop only
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
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = window.globalcell(x, y);
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
