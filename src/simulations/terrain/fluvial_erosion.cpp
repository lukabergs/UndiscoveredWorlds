#include "simulations/terrain/fluvial_erosion.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <vector>

#include "core/climate_seasons.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/landscape/glacial.hpp"
#include "simulations/landscape/isostasy.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/relief_roughness.hpp"
#include "simulations/landscape/stream_power.hpp"
#include "world/planet.hpp"

namespace
{
// Salt of the relief roughness's seed (landscape::ReliefRoughness::seed).
constexpr int roughnessseedsalt = 0x5266;

// K per cell: the crust-class multiplier times the age factor; the field
// geosphere.erodibility stores this relative value, the kernel gets it times
// landscape.stream_power.erodibility.
float relativeerodibility(const planet& world, mesh::CellId c)
{
    namespace sp = tuning::landscape::stream_power;
    float multiplier = sp::crustClassMultiplierContinental;
    switch (world.tectoniccrustclass(c))
    {
    case CrustClass::oceanic: multiplier = sp::crustClassMultiplierOceanic; break;
    case CrustClass::transitional: multiplier = sp::crustClassMultiplierTransitional; break;
    default: break;
    }
    const float age = std::max(0.0f, world.tectoniccrustagemyr(c));
    return multiplier * (1.0f + sp::crustAgeFactor * age / 1000.0f);
}

landscape::StreamPowerParameters kernelparameters(double baselevel)
{
    namespace t = tuning::landscape;
    landscape::StreamPowerParameters parameters;
    parameters.step_myr = t::erosion::stepMyr;
    parameters.baselevel_m = baselevel;
    parameters.areaexponent = t::stream_power::areaExponent;
    parameters.slopeexponent = t::stream_power::slopeExponent;
    parameters.depositioncoefficient = t::stream_power::depositionCoefficient;
    parameters.depositiontolerance_m = t::stream_power::depositionToleranceM;
    parameters.depositionmaxsweeps = t::stream_power::depositionMaxSweeps;
    parameters.hillslopediffusivity_km2_myr = t::hillslope::diffusivityKm2Myr;
    parameters.criticalslope = t::hillslope::criticalSlope;
    parameters.sedimenttravel_km = t::marine::sedimentTravelKm;
    parameters.deltatop_m = t::marine::deltaTopM;
    parameters.upliftsea = true;  // erosion plan D1: the engine's uplift moves the coast both ways
    return parameters;
}

// The glacial term of the run (erosion plan Phase G; section 6, "Decided for Phase G"):
// off with landscape.glacial.enabled false, when no vector is filled and the run is the
// fluvial one.
struct GlacialRun
{
    bool enabled = false;
    landscape::GlacialParameters parameters;
    std::vector<float> ela;          // the run's ELA, with the cooling (m)
    std::vector<float> presentela;   // the ELA without the cooling: the present-day ice at the end (m)
    std::vector<float> thickness;    // h of the last step (m); the next step's sea mask and mass balance use the bed plus h
    std::vector<double> iceflux;     // m^3/yr of the step
    std::vector<float> stepglacial;  // the glacial part of the step's erosion (m)
    std::vector<double> denudation;  // the glacial part of the run's denudation (m)
    std::vector<double> icetime;     // Myr during which the cell carried ice (h > 0)
    // Diagnostics of the summary line (no field reads them): cells whose h went from 0 to
    // > 0 (appearances) or back (disappearances) from one step to the next.
    long long appearances = 0;
    long long disappearances = 0;
    long long reroutedappearances = 0;     // of those, on cells whose receiver changed between the two steps
    long long rerouteddisappearances = 0;
};

landscape::GlacialParameters glacialparameters()
{
    namespace g = tuning::landscape::glacial;
    landscape::GlacialParameters parameters;
    parameters.erodibilityratio = g::erodibilityRatio;
    parameters.cooling_c = g::coolingC;
    parameters.fullicerange_m = g::fullIceRangeM;
    parameters.widthreference_m = g::widthReferenceM;
    parameters.widthexponent = g::widthExponent;
    parameters.thicknessratio = g::thicknessRatio;
    parameters.maximumthickness_m = g::maximumThicknessM;
    parameters.flotationratio = g::flotationRatio;
    return parameters;
}

// The ELAs from the coarse climate on the terrain it saw (`inputs.elevation`, before the
// roughness): the warmest season's temperature, the annual precipitation (mm/yr, the
// seasonal mm/month x 3 clamped at 0, as gatherinputs) and the climate's lapse rate.
void startglacialrun(const planet& world, const flowrouting::Inputs& inputs, GlacialRun& glacial)
{
    const landscape::SphereMetric& metric = inputs.metric;
    const std::size_t size = metric.size();
    glacial.enabled = true;
    glacial.parameters = glacialparameters();
    std::vector<float> warmest(size), precipitation(size);
    // The metric's cell index is the world's cell id (both in the grid's id order).
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        float warm = world.seasonaltempfloat(0, c);
        double annual_mm = 0.0;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            warm = std::max(warm, world.seasonaltempfloat(season, c));
            annual_mm += 3.0 * std::max(0.0f, world.seasonalrainfloat(season, c));
        }
        warmest[cell] = warm;
        precipitation[cell] = static_cast<float>(annual_mm);
    }
    // z of the ELA is the height the climate's temperature belongs to: the climate lapses
    // only above the sea (max(0, mapfloat - sealevel), climate_energy.cpp), so a cell below
    // the sea level when the coarse climate ran has its temperature at the sea level, not at
    // its bed; with the bed, such a cell's ELA would be too low by its depth and it would
    // glaciate as soon as uplift or progradation made it land. Procedural worlds carry no
    // template mountains, so the base elevation (nomfloat) above the sea is the mapfloat the
    // climate saw.
    std::vector<float> climateelevation(size);
    for (std::size_t cell = 0; cell < size; cell++)
        climateelevation[cell] = std::max(inputs.elevation[cell], inputs.sealevel);
    const double lapse = static_cast<double>(world.tempdecrease()) / 1000.0;  // C per 1000 m -> C per m
    landscape::equilibriumline(metric, climateelevation, warmest, precipitation, lapse, glacial.parameters.cooling_c,
                               glacial.ela);
    landscape::equilibriumline(metric, climateelevation, warmest, precipitation, lapse, 0.0, glacial.presentela);
    glacial.thickness.assign(size, 0.0f);
    glacial.denudation.assign(size, 0.0);
    glacial.icetime.assign(size, 0.0);
}

// Whether the step from `previous` to `elapsed` (Myr) reaches or passes a requested checkpoint time.
bool reachescheckpoint(const FluvialErosionCheckpoints& checkpoints, double previous, double elapsed)
{
    constexpr double tolerance = 1.0e-9;
    if (checkpoints.everyMyr > 0.0 &&
        std::floor(elapsed / checkpoints.everyMyr + tolerance) > std::floor(previous / checkpoints.everyMyr + tolerance))
        return true;
    for (const double time : checkpoints.atMyr)
        if (time > previous + tolerance && time <= elapsed + tolerance)
            return true;
    return false;
}

// The run's state as of `elapsed` (a checkpoint, or the end with the duration):
// base elevation, the erosion fields with rates over `elapsed`, the volumes so far
// and the routing of the current terrain; the glacial fields (zero when the term is
// off), and at the end of the run (`end`) also the present-day ice: the ELA without
// the cooling on the final terrain and its routing.
void writeerosionstate(planet& world, const flowrouting::Inputs& inputs, const std::vector<float>& elevation,
    const std::vector<float>& erodibility, const std::vector<double>& denudation, const std::vector<double>& aggradation,
    const std::vector<double>& sediment, const landscape::StreamPowerBudget& total, double elapsed,
    const GlacialRun& glacial, bool end, std::vector<std::uint8_t>* finalsea = nullptr)
{
    namespace t = tuning::landscape;
    const flowrouting::Routing routing = flowrouting::route(inputs, &elevation);
    flowrouting::writeroutingfields(world, inputs, routing);
    if (finalsea != nullptr)
        *finalsea = routing.sea;

    const landscape::SphereMetric& metric = inputs.metric;
    std::vector<float> presentice;
    if (end && glacial.enabled)
    {
        std::vector<double> presentflux;
        landscape::accumulateice(metric, routing.graph, elevation, glacial.presentela, inputs.runoff_m_yr,
                                 glacial.parameters.fullicerange_m, presentflux);
        landscape::icethickness(presentflux, glacial.parameters, presentice);
    }
    worldgrid<float>& equilibriumline = world.equilibriumline();
    worldgrid<float>& glaciation = world.glaciation();
    worldgrid<float>& glacialdenudation = world.glacialdenudation();
    worldgrid<float>& icethickness = world.icethickness();
    worldgrid<float>& erodibilityfield = world.erodibility();
    worldgrid<float>& erosionrate = world.erosionrate();
    worldgrid<float>& depositionrate = world.depositionrate();
    worldgrid<float>& denudationfield = world.denudation();
    worldgrid<float>& aggradationfield = world.aggradation();
    worldgrid<float>& sedimentflux = world.sedimentflux();
    worldgrid<std::uint8_t>& prograded = world.prograded();
    const double megayears = std::max(elapsed, 1.0e-9);
    constexpr double yearspermegayear = 1.0e6;
    for (std::size_t cell = 0; cell < metric.size(); cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        world.setnomfloat(c, elevation[cell]);
        erodibilityfield[c] = routing.sea[cell] ? 0.0f : erodibility[cell] / t::stream_power::erodibility;
        erosionrate[c] = static_cast<float>(denudation[cell] / megayears);
        depositionrate[c] = static_cast<float>(aggradation[cell] / megayears);
        denudationfield[c] = static_cast<float>(denudation[cell]);
        aggradationfield[c] = static_cast<float>(aggradation[cell]);
        // hydrosphere.sediment_flux: the mean over the elapsed time in m^3/yr; 0 on the current sea.
        sedimentflux[c] = routing.sea[cell] ? 0.0f
            : static_cast<float>(sediment[cell] / (megayears * yearspermegayear));
        // geosphere.prograded (erosion plan C5b): sea when the run started and now above the sea
        // level (raised by deposition); an arm of the sea merely cut off from the ocean (H3) is not.
        prograded[c] = inputs.sea[cell] != 0 && !flowrouting::belowsea(elevation[cell], inputs.sealevel) ? 1 : 0;
        // The glacial fields (erosion plan Phase G): the run's ELA, the share of the elapsed
        // time with ice, the glacial part of the denudation; the present-day ice at the end.
        equilibriumline[c] = glacial.enabled ? glacial.ela[cell] : 0.0f;
        glaciation[c] = glacial.enabled
            ? static_cast<float>(std::min(1.0, glacial.icetime[cell] / megayears)) : 0.0f;
        glacialdenudation[c] = glacial.enabled ? static_cast<float>(glacial.denudation[cell]) : 0.0f;
        if (end)
            icethickness[c] = glacial.enabled ? presentice[cell] : 0.0f;
    }
    world.seterosionduration(static_cast<float>(elapsed));
    world.seterodedvolume(static_cast<float>(total.eroded_km3));
    world.setlanddepositedvolume(static_cast<float>(total.landdeposited_km3));
    world.setmarinedepositedvolume(static_cast<float>(total.marinedeposited_km3));
}
}

void applyfluvialerosionstage(planet& world, const FluvialErosionCheckpoints& checkpoints)
{
    namespace t = tuning::landscape;
    flowrouting::Inputs inputs = flowrouting::gatherinputs(world);
    // The evolver's sea is the eroded rule from the start (seamask: the world ocean, erosion plan
    // H3), as at every routing of the run, whatever erosion_duration says before the run writes it.
    inputs.sea = flowrouting::seamask(inputs.metric, inputs.elevation, inputs.sealevel);
    const landscape::SphereMetric& metric = inputs.metric;
    const std::size_t size = metric.size();

    std::vector<float> erodibility(size), uplift(size), elevation = inputs.elevation;
    std::vector<float> steperosion(size), stepdeposition(size);
    std::vector<double> denudation(size, 0.0), aggradation(size, 0.0);
    std::vector<double> stepflux, sediment(size, 0.0);  // m^3/Myr of a step; m^3 summed over the run
    // The metric's cell index is the world's cell id (both in the grid's id order).
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        erodibility[cell] = relativeerodibility(world, c) * t::stream_power::erodibility;
        // U (erosion plan D1): the uniform uplift plus the engine's rock uplift rate, m/Myr,
        // signed; the kernel applies it to land and sea cells (upliftsea).
        uplift[cell] = static_cast<float>(t::erosion::upliftMMyr
            + t::erosion::engineUpliftScale * world.tectonicrockupliftrate(c));
    }

    // Relief roughness (erosion plan C3): the land's sub-cell relief, added once to
    // the working elevation; it stays in the eroded terrain. It can raise land, so
    // the result is held to the datum ceiling geosphere.max_elevation.
    if (t::roughness::amplitudeM != 0.0)
    {
        landscape::ReliefRoughness roughness;
        roughness.amplitude_m = t::roughness::amplitudeM;
        roughness.wavelength_km = t::roughness::wavelengthKm;
        roughness.octaves = t::roughness::octaves;
        roughness.persistence = t::roughness::persistence;
        roughness.seed = deterministiccontextseed(world.seed(), roughnessseedsalt);
        landscape::addreliefroughness(metric, roughness, inputs.sea, inputs.sealevel, elevation);
        const float ceiling = static_cast<float>(world.maxelevation());
        for (float& height : elevation)
            height = std::min(height, ceiling);
    }

    // The glacial term (erosion plan Phase G): the ELA per cell for the whole run; at every
    // step the routing runs on the bed with the sea mask of the previous step's ice surface
    // (bed + h, h = 0 at the first), the ice flux (mass balance on that surface) and the
    // thickness come from that routing and the kernel erodes with them. Off, nothing below
    // differs from the fluvial run.
    GlacialRun glacial;
    if (t::glacial::enabled)
        startglacialrun(world, inputs, glacial);
    std::vector<float> surface;
    std::vector<std::uint8_t> hadice;       // the previous step's h > 0, for the summary line only
    std::vector<std::int32_t> previousreceiver;  // the previous step's receivers, for the summary line only

    const landscape::StreamPowerParameters parameters = kernelparameters(inputs.sealevel);
    const landscape::IsostasyParameters isostasy{t::isostasy::reboundFraction, t::isostasy::loadingFraction};
    const float ceiling = static_cast<float>(world.maxelevation());
    const double duration = t::erosion::durationMyr;
    landscape::StreamPowerBudget total;
    double elapsed = 0.0;
    flowrouting::Routing routing;
    while (elapsed < duration)
    {
        landscape::StreamPowerParameters step = parameters;
        step.step_myr = std::min(parameters.step_myr, duration - elapsed);
        landscape::StreamPowerBudget budget;
        if (glacial.enabled)
        {
            // The ice surface of the previous step gives the sea mask (a bed below the sea level
            // under enough ice is land: grounded ice, through which the ice reaches the coast)
            // and the mass balance; receivers, depressions and accumulation are the bed's, as in
            // Hergarten 2021. Routing on the ice surface lifted an iced cell above its donors,
            // which then drained elsewhere, so the ice came and went from step to step (G3:
            // 643 appearances and 510 disappearances on 305 cells at the minimal size, against
            // 122 and 13 on 209 cells with the bed's routing).
            surface.resize(size);
            for (std::size_t cell = 0; cell < size; cell++)
                surface[cell] = elevation[cell] + glacial.thickness[cell];
            const std::vector<std::uint8_t> icesea = flowrouting::seamask(metric, surface, inputs.sealevel);
            routing = flowrouting::route(inputs, &elevation, &icesea);
            landscape::accumulateice(metric, routing.graph, surface, glacial.ela, inputs.runoff_m_yr,
                                     glacial.parameters.fullicerange_m, glacial.iceflux);
            hadice.resize(size);
            for (std::size_t cell = 0; cell < size; cell++)
                hadice[cell] = glacial.thickness[cell] > 0.0f ? 1 : 0;
            landscape::icethickness(glacial.iceflux, glacial.parameters, glacial.thickness);
            if (elapsed > 0.0)  // between steps: the first step starts from no ice by construction
                for (std::size_t cell = 0; cell < size; cell++)
                {
                    const bool hasice = glacial.thickness[cell] > 0.0f;
                    const bool rerouted = routing.graph.receiver[cell] != previousreceiver[cell];
                    if (hasice && hadice[cell] == 0)
                    {
                        glacial.appearances++;
                        if (rerouted)
                            glacial.reroutedappearances++;
                    }
                    else if (!hasice && hadice[cell] != 0)
                    {
                        glacial.disappearances++;
                        if (rerouted)
                            glacial.rerouteddisappearances++;
                    }
                }
            previousreceiver = routing.graph.receiver;
            landscape::GlacialInput input;
            input.iceflux_m3_yr = &glacial.iceflux;
            input.thickness_m = &glacial.thickness;
            input.erodibilityratio = glacial.parameters.erodibilityratio;
            input.flotationratio = glacial.parameters.flotationratio;
            input.glacialerosion_m = &glacial.stepglacial;
            budget = landscape::stepstreampower(metric, routing.graph, erodibility, uplift, step, elevation,
                                                steperosion, stepdeposition, &stepflux, &input);
            for (std::size_t cell = 0; cell < size; cell++)
            {
                glacial.denudation[cell] += glacial.stepglacial[cell];
                if (glacial.thickness[cell] > 0.0f)
                    glacial.icetime[cell] += step.step_myr;
            }
        }
        else
        {
            routing = flowrouting::route(inputs, &elevation);
            budget = landscape::stepstreampower(
                metric, routing.graph, erodibility, uplift, step, elevation, steperosion, stepdeposition, &stepflux);
        }
        // Local isostasy (erosion plan C8): the step's removed rock lifts each column by
        // a fraction, its laid sediment presses it down. The erosion fields keep counting
        // the rock removed and the sediment laid (layers.finalize's exhumation reads
        // them); the net elevation change of the run now includes the response. A cell
        // that both lost and gained in the step can end above where the kernel left it,
        // so the result is held to the datum ceiling, as after the roughness.
        landscape::applyisostasy(isostasy, steperosion, stepdeposition, elevation);
        for (float& height : elevation)
            height = std::min(height, ceiling);
        total.eroded_km3 += budget.eroded_km3;
        total.landdeposited_km3 += budget.landdeposited_km3;
        total.marinedeposited_km3 += budget.marinedeposited_km3;
        total.depositionsweeps += budget.depositionsweeps;
        total.strandedloads += budget.strandedloads;
        total.stranded_km3 += budget.stranded_km3;
        for (std::size_t cell = 0; cell < size; cell++)
        {
            denudation[cell] += steperosion[cell];
            aggradation[cell] += stepdeposition[cell];
            sediment[cell] += stepflux[cell] * step.step_myr;
        }
        const double previous = elapsed;
        elapsed += step.step_myr;
        // Checkpoint hook (time model): the end of the run is not a checkpoint.
        if (checkpoints.checkpoint && elapsed < duration && reachescheckpoint(checkpoints, previous, elapsed))
        {
            writeerosionstate(world, inputs, elevation, erodibility, denudation, aggradation, sediment, total, elapsed,
                              glacial, false);
            checkpoints.checkpoint(elapsed);
        }
    }

    // The end of the run: the final routing on the eroded terrain, so that no second
    // routing step follows. geosphere.max_elevation is the datum ceiling the terrain
    // stages clamp to, not the summit; the roughness and every step's isostatic
    // response were held to it above (with uplift 0 the stream-power steps lower the
    // terrain or deposit below the terrain upstream, but the rebound can lift a cell
    // that both lost and gained above where the kernel left it), so it stays valid.
    std::vector<std::uint8_t> finalsea;
    writeerosionstate(world, inputs, elevation, erodibility, denudation, aggradation, sediment, total, duration,
                      glacial, true, &finalsea);
    if (total.strandedloads > 0)
        std::cout << "Fluvial erosion: " << total.strandedloads << " marine load(s) found no sea cell below the delta top ("
                  << total.stranded_km3 << " km3 of " << total.marinedeposited_km3 << " km3 marine)\n";
    if (glacial.enabled)
    {
        // Diagnostics only (erosion plan G3): no field depends on them. "Fjord" cells are land in
        // the evolver's starting sea mask, sea in the final routing's and glacially eroded.
        long long icedcells = 0, fjordcells = 0;
        double deepest = 0.0;  // m below the sea level
        for (std::size_t cell = 0; cell < size; cell++)
        {
            if (glacial.icetime[cell] > 0.0)
                icedcells++;
            if (inputs.sea[cell] == 0 && finalsea[cell] != 0 && glacial.denudation[cell] > 0.0)
            {
                fjordcells++;
                deepest = std::max(deepest, static_cast<double>(inputs.sealevel) - static_cast<double>(elevation[cell]));
            }
        }
        std::cout << "Glacial erosion: " << icedcells << " cells carried ice at some step, " << glacial.appearances
                  << " ice appearances and " << glacial.disappearances << " disappearances between steps ("
                  << glacial.reroutedappearances << " and " << glacial.rerouteddisappearances
                  << " on cells whose receiver changed), " << fjordcells
                  << " cells land at the start and sea at the end with glacial denudation, deepest bed "
                  << deepest << " m below the sea level among them\n";
    }
}
