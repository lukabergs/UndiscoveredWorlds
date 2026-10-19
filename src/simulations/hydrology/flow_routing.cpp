#include "simulations/hydrology/flow_routing.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

#include "core/climate_seasons.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/ocean_mask.hpp"
#include "simulations/landscape/runoff.hpp"
#include "world/planet.hpp"

namespace flowrouting
{
namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;
// Salt of the routing's per-cell draws (landscape::RoutingOptions::seed).
constexpr int routingseedsalt = 0x7207;
}

landscape::SphereMetric worldmetric(const planet& world)
{
    // The world's grid, shared with the world and its copies.
    return landscape::SphereMetric(world.sharedgrid());
}

Inputs gatherinputs(const planet& world)
{
    namespace runoff = tuning::landscape::runoff;
    namespace transmission = tuning::landscape::transmission;
    namespace lakes = tuning::landscape::lakes;
    Inputs inputs;
    inputs.metric = worldmetric(world);
    const std::size_t size = inputs.metric.size();
    inputs.elevation.resize(size);
    inputs.sea.resize(size);
    inputs.runoff_m_yr.resize(size);
    inputs.sealevel = static_cast<float>(world.sealevel());
    inputs.routingoptions.randomised = tuning::landscape::routing::randomised;
    inputs.routingoptions.resolveflats = tuning::landscape::routing::resolveFlats;
    if (inputs.routingoptions.randomised)  // the seed is read only when a draw needs it
        inputs.routingoptions.seed = deterministiccontextseed(world.seed(), routingseedsalt);
    // The imported lakes are read only when the losses need them.
    const bool budyko = runoff::budyko;
    inputs.transmissionenabled = transmission::enabled;
    inputs.evaporation_m_yr.resize(size);
    const worldgrid<std::uint8_t>* importedlake = nullptr;
    if (inputs.transmissionenabled)
    {
        importedlake = &world.importedlake();
        landscape::TransmissionLoss& loss = inputs.transmission;
        loss.deficit_m_yr.resize(size);
        loss.openwater.resize(size);
        loss.widthcoefficient = transmission::widthCoefficient;
        loss.widthexponent = transmission::widthExponent;
        loss.floodplainfactor = transmission::floodplainFactor;
        loss.floodplainslope = transmission::floodplainSlope;
        loss.sealevel = inputs.sealevel;
    }
    // The metric's cell index is the world's cell id (both in the grid's id order).
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        inputs.elevation[cell] = world.nomfloat(c);
        // Runoff is computed for every cell, sea included (precipitation exists over the
        // sea): a seabed that marine deposition raises above the sea level during the
        // erosion run then contributes to the next routing; accumulateflow ignores sea cells.
        // Seasonal precipitation is mm/month; four seasons of three months make the year. The
        // potential evaporation is the open-water law of landscape.lakes applied to each
        // season's temperature and averaged over the four (m/yr): the one law of the runoff,
        // the transmission losses and the lake balance.
        double annual_mm = 0.0;
        double evaporation = 0.0;
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            annual_mm += 3.0 * std::max(0.0f, world.seasonalrainfloat(season, c));
            evaporation += landscape::openwaterevaporation(world.seasonaltempfloat(season, c),
                                                           lakes::evaporationPerC, lakes::evaporationZeroC);
        }
        evaporation /= CLIMATESEASONCOUNT;
        inputs.evaporation_m_yr[cell] = static_cast<float>(evaporation);
        double depth;
        if (budyko)  // Fu's curve: a runoff ratio near 1 where PET is small, so no snow or ice term
            depth = landscape::budykorunoff(annual_mm / 1000.0, evaporation, runoff::budykoW);
        else
        {
            int coldseasons = 0;
            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                if (world.seasonaltempfloat(season, c) < runoff::snowTemperatureC)
                    coldseasons++;
            const double snowfraction = static_cast<double>(coldseasons) / CLIMATESEASONCOUNT;
            const double icefraction = coldseasons == CLIMATESEASONCOUNT ? 1.0 : 0.0;
            depth = runoff::coefficient * annual_mm / 1000.0 *
                (1.0 - runoff::snowReduction * snowfraction) * (1.0 - runoff::iceReduction * icefraction);
        }
        inputs.runoff_m_yr[cell] = static_cast<float>(std::max(0.0, depth));
        if (inputs.transmissionenabled)
        {
            inputs.transmission.deficit_m_yr[cell] =
                static_cast<float>(std::max(0.0, evaporation - annual_mm / 1000.0));
            inputs.transmission.openwater[cell] = (*importedlake)[c] != 0 ? 1 : 0;
        }
    }
    // Only the below-sea cells of the ocean regions are sea (erosion plan H3).
    inputs.sea = oceancells(world);
    return inputs;
}

bool belowsea(float elevation, float sealevel)
{
    return static_cast<float>(std::lround(elevation)) <= sealevel;
}

std::vector<uint8_t> seamask(const landscape::SphereMetric& metric, const std::vector<float>& elevation, float sealevel)
{
    std::vector<uint8_t> below(elevation.size());
    for (std::size_t cell = 0; cell < elevation.size(); cell++)
        below[cell] = belowsea(elevation[cell], sealevel) ? 1 : 0;
    return landscape::oceanmask(metric, below, tuning::landscape::routing::oceanMinimumAreaFraction);
}

std::vector<uint8_t> oceancells(const planet& world)
{
    const landscape::SphereMetric metric = worldmetric(world);
    std::vector<uint8_t> below(metric.size());
    for (std::size_t cell = 0; cell < below.size(); cell++)
        below[cell] = world.basebelowsea(static_cast<mesh::CellId>(cell)) ? 1 : 0;
    if (world.erosionduration() > 0.0f)
        return landscape::oceanmask(metric, below, tuning::landscape::routing::oceanMinimumAreaFraction);
    // No erosion (imported terrain: imported-world, earth-hydrology): the imported sea map decides,
    // so every below-sea cell is sea unless the imported lake mask marks it.
    const worldgrid<std::uint8_t>& importedlake = world.importedlake();
    for (std::size_t cell = 0; cell < below.size(); cell++)
        if (importedlake[static_cast<mesh::CellId>(cell)] != 0)
            below[cell] = 0;
    return below;
}

Routing route(const Inputs& inputs, const std::vector<float>* elevation, const std::vector<uint8_t>* sea)
{
    const std::vector<float>& terrain = elevation != nullptr ? *elevation : inputs.elevation;
    Routing routing;
    if (elevation != nullptr && sea != nullptr)
        routing.sea = *sea;
    else
        routing.sea = elevation != nullptr ? seamask(inputs.metric, terrain, inputs.sealevel) : inputs.sea;
    routing.depressions = landscape::filldepressions(inputs.metric, terrain, routing.sea,
                                                     tuning::landscape::depressions::minimumDepthM);
    routing.graph = landscape::buildflowgraph(inputs.metric, terrain, routing.sea, inputs.sealevel, routing.depressions,
                                              inputs.routingoptions);
    if (inputs.transmissionenabled)
    {
        // The cells of this routing's depressions lose nothing to the channel: the lake balance
        // (hydrology.lake_levels) owns their evaporation.
        landscape::TransmissionLoss loss = inputs.transmission;
        loss.depressionid = routing.depressions.depressionid;
        routing.transmissionloss_m3_yr = landscape::accumulateflowwithlosses(
            inputs.metric, inputs.runoff_m_yr, terrain, loss, routing.graph);
    }
    else
        landscape::accumulateflow(inputs.metric, inputs.runoff_m_yr, routing.graph);
    return routing;
}

void writeroutingfields(planet& world, const Inputs& inputs, const Routing& routing)
{
    worldgrid<std::uint8_t>& direction = world.flowdirection();
    worldgrid<float>& area = world.drainagearea();
    worldgrid<float>& discharge = world.discharge();
    worldgrid<float>& filled = world.filledelevation();
    worldgrid<std::int32_t>& depression = world.depressionid();
    worldgrid<float>& runoff = world.runoff();
    for (std::size_t cell = 0; cell < inputs.metric.size(); cell++)
    {
        const mesh::CellId c = static_cast<mesh::CellId>(cell);
        direction[c] = routing.graph.direction[cell];
        area[c] = static_cast<float>(routing.graph.area_m2[cell] / 1.0e6);
        discharge[c] = static_cast<float>(routing.graph.discharge_m3_yr[cell] / secondsperyear);
        filled[c] = routing.depressions.filledelevation[cell];
        depression[c] = routing.depressions.depressionid[cell];
        runoff[c] = inputs.runoff_m_yr[cell] * 1000.0f;
    }
}
}

void applyflowroutingstage(planet& world)
{
    const flowrouting::Inputs inputs = flowrouting::gatherinputs(world);
    const flowrouting::Routing routing = flowrouting::route(inputs);
    flowrouting::writeroutingfields(world, inputs, routing);
}
