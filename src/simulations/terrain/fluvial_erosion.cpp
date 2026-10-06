#include "simulations/terrain/fluvial_erosion.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <vector>

#include "core/deterministic_random.hpp"
#include "simulations/hydrology/flow_routing.hpp"
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
float relativeerodibility(const planet& world, int x, int y)
{
    namespace sp = tuning::landscape::stream_power;
    float multiplier = sp::crustClassMultiplierContinental;
    switch (world.tectoniccrustclass(x, y))
    {
    case CrustClass::oceanic: multiplier = sp::crustClassMultiplierOceanic; break;
    case CrustClass::transitional: multiplier = sp::crustClassMultiplierTransitional; break;
    default: break;
    }
    const float age = std::max(0.0f, world.tectoniccrustagemyr(x, y));
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
    return parameters;
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
// and the routing of the current terrain.
void writeerosionstate(planet& world, const flowrouting::Inputs& inputs, const std::vector<float>& elevation,
    const std::vector<float>& erodibility, const std::vector<double>& denudation, const std::vector<double>& aggradation,
    const std::vector<double>& sediment, const landscape::StreamPowerBudget& total, double elapsed)
{
    namespace t = tuning::landscape;
    const flowrouting::Routing routing = flowrouting::route(inputs, &elevation);
    flowrouting::writeroutingfields(world, inputs, routing);

    const landscape::SphereMetric& metric = inputs.metric;
    worldgrid<float>& erodibilityfield = world.erodibility();
    worldgrid<float>& erosionrate = world.erosionrate();
    worldgrid<float>& depositionrate = world.depositionrate();
    worldgrid<float>& denudationfield = world.denudation();
    worldgrid<float>& aggradationfield = world.aggradation();
    worldgrid<float>& sedimentflux = world.sedimentflux();
    worldgrid<std::uint8_t>& prograded = world.prograded();
    const double megayears = std::max(elapsed, 1.0e-9);
    constexpr double yearspermegayear = 1.0e6;
    for (int y = 0; y < metric.height; y++)
    {
        for (int x = 0; x < metric.width; x++)
        {
            const std::size_t cell = metric.index(x, y);
            world.setnomfloat(x, y, elevation[cell]);
            erodibilityfield[x][y] = routing.sea[cell] ? 0.0f : erodibility[cell] / t::stream_power::erodibility;
            erosionrate[x][y] = static_cast<float>(denudation[cell] / megayears);
            depositionrate[x][y] = static_cast<float>(aggradation[cell] / megayears);
            denudationfield[x][y] = static_cast<float>(denudation[cell]);
            aggradationfield[x][y] = static_cast<float>(aggradation[cell]);
            // hydrosphere.sediment_flux: the mean over the elapsed time in m^3/yr; 0 on the current sea.
            sedimentflux[x][y] = routing.sea[cell] ? 0.0f
                : static_cast<float>(sediment[cell] / (megayears * yearspermegayear));
            // geosphere.prograded (erosion plan C5b): sea when the run started, land now.
            prograded[x][y] = inputs.sea[cell] != 0 && routing.sea[cell] == 0 ? 1 : 0;
        }
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
    const landscape::SphereMetric& metric = inputs.metric;
    const int columns = metric.width;
    const int rows = metric.height;
    const std::size_t size = metric.size();

    std::vector<float> erodibility(size), uplift(size), elevation = inputs.elevation;
    std::vector<float> steperosion(size), stepdeposition(size);
    std::vector<double> denudation(size, 0.0), aggradation(size, 0.0);
    std::vector<double> stepflux, sediment(size, 0.0);  // m^3/Myr of a step; m^3 summed over the run
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            erodibility[cell] = relativeerodibility(world, x, y) * t::stream_power::erodibility;
            uplift[cell] = t::erosion::upliftMMyr;  // the kernel applies it to land cells of the step's own mask
        }
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
        routing = flowrouting::route(inputs, &elevation);
        const landscape::StreamPowerBudget budget = landscape::stepstreampower(
            metric, routing.graph, erodibility, uplift, step, elevation, steperosion, stepdeposition, &stepflux);
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
            writeerosionstate(world, inputs, elevation, erodibility, denudation, aggradation, sediment, total, elapsed);
            checkpoints.checkpoint(elapsed);
        }
    }

    // The end of the run: the final routing on the eroded terrain, so that no second
    // routing step follows. geosphere.max_elevation is the datum ceiling the terrain
    // stages clamp to, not the summit; the roughness and every step's isostatic
    // response were held to it above (with uplift 0 the stream-power steps lower the
    // terrain or deposit below the terrain upstream, but the rebound can lift a cell
    // that both lost and gained above where the kernel left it), so it stays valid.
    writeerosionstate(world, inputs, elevation, erodibility, denudation, aggradation, sediment, total, duration);
    if (total.strandedloads > 0)
        std::cout << "Fluvial erosion: " << total.strandedloads << " marine load(s) found no sea cell below the delta top ("
                  << total.stranded_km3 << " km3 of " << total.marinedeposited_km3 << " km3 marine)\n";
}
