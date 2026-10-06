#include "simulations/terrain/fluvial_erosion.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "simulations/hydrology/flow_routing.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/stream_power.hpp"
#include "world/planet.hpp"

namespace
{
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
    return parameters;
}
}

void applyfluvialerosionstage(planet& world)
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
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            erodibility[cell] = relativeerodibility(world, x, y) * t::stream_power::erodibility;
            uplift[cell] = t::erosion::upliftMMyr;  // the kernel applies it to land cells of the step's own mask
        }
    }

    const landscape::StreamPowerParameters parameters = kernelparameters(inputs.sealevel);
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
            metric, routing.graph, erodibility, uplift, step, elevation, steperosion, stepdeposition);
        total.eroded_km3 += budget.eroded_km3;
        total.landdeposited_km3 += budget.landdeposited_km3;
        total.marinedeposited_km3 += budget.marinedeposited_km3;
        total.depositionsweeps += budget.depositionsweeps;
        for (std::size_t cell = 0; cell < size; cell++)
        {
            denudation[cell] += steperosion[cell];
            aggradation[cell] += stepdeposition[cell];
        }
        elapsed += step.step_myr;
    }

    // The final routing on the eroded terrain, so that no second routing step follows.
    routing = flowrouting::route(inputs, &elevation);
    flowrouting::writeroutingfields(world, inputs, routing);

    worldgrid<float>& erodibilityfield = world.erodibility();
    worldgrid<float>& erosionrate = world.erosionrate();
    worldgrid<float>& depositionrate = world.depositionrate();
    worldgrid<float>& denudationfield = world.denudation();
    worldgrid<float>& aggradationfield = world.aggradation();
    const double megayears = std::max(duration, 1.0e-9);
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = metric.index(x, y);
            world.setnomfloat(x, y, elevation[cell]);
            erodibilityfield[x][y] = routing.sea[cell] ? 0.0f : erodibility[cell] / t::stream_power::erodibility;
            erosionrate[x][y] = static_cast<float>(denudation[cell] / megayears);
            depositionrate[x][y] = static_cast<float>(aggradation[cell] / megayears);
            denudationfield[x][y] = static_cast<float>(denudation[cell]);
            aggradationfield[x][y] = static_cast<float>(aggradation[cell]);
        }
    }
    // geosphere.max_elevation is the datum ceiling the terrain stages clamp to, not
    // the summit; erosion only lowers the terrain, so the ceiling stays valid.
    world.seterosionduration(static_cast<float>(duration));
    world.seterodedvolume(static_cast<float>(total.eroded_km3));
    world.setlanddepositedvolume(static_cast<float>(total.landdeposited_km3));
    world.setmarinedepositedvolume(static_cast<float>(total.marinedeposited_km3));
}
