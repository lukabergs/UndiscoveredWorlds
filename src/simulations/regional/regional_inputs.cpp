#include "simulations/regional/regional_inputs.hpp"

#include "core/deterministic_random.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "world/planet.hpp"

namespace
{
// Salt of the regional channels' meanders (landscape::ChannelParameters::seed).
constexpr int channelseedsalt = 0x7e91;
}

landscape::RegionalGlobalInputs gatherregionalinputs(const planet& world)
{
    const int columns = world.width() + 1;
    const int rows = world.height() + 1;
    const std::size_t size = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
    landscape::RegionalGlobalInputs inputs;
    inputs.width = columns;
    inputs.height = rows;
    inputs.radius_m = static_cast<double>(world.radius());
    inputs.sealevel = static_cast<float>(world.sealevel());
    inputs.seed = static_cast<std::uint64_t>(static_cast<std::int64_t>(world.seed()));
    inputs.lapserate_c_per_km = world.tempdecrease();

    landscape::GlobalNetwork& network = inputs.network;
    network.width = columns;
    network.height = rows;
    // The sea of the routing's own state: a cell without a receiver whose base lies
    // below the sea level (every other below-sea cell is a lake; probes of the minimal
    // bundle and erosion run 9 found this equal to flowrouting::oceancells), so the
    // channels' sea and receivers agree by construction.
    network.sea.resize(size);
    network.direction.resize(size);
    network.drainagearea_km2.resize(size);
    network.discharge_m3s.resize(size);
    network.filledelevation_m.resize(size);
    network.januaryflow.resize(size);
    network.julyflow.resize(size);
    inputs.elevation_m.resize(size);
    inputs.runoff_mm_yr.resize(size);
    inputs.lakesurface_m.resize(size);
    inputs.januarytemperature_c.resize(size);
    inputs.julytemperature_c.resize(size);
    inputs.januaryprecipitation_mm.resize(size);
    inputs.julyprecipitation_mm.resize(size);
    inputs.icethickness_m.resize(size);
    inputs.climate.resize(size);

    const worldgrid<std::uint8_t>& direction = world.flowdirection();
    const worldgrid<float>& area = world.drainagearea();
    const worldgrid<float>& discharge = world.discharge();
    const worldgrid<float>& filled = world.filledelevation();
    const worldgrid<float>& runoff = world.runoff();
    const worldgrid<float>& ice = world.icethickness();
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(columns) + static_cast<std::size_t>(x);
            inputs.elevation_m[cell] = world.nomfloat(x, y);
            network.direction[cell] = direction[x][y];
            network.sea[cell] = direction[x][y] == 0 && world.basebelowsea(x, y) ? 1 : 0;
            network.drainagearea_km2[cell] = area[x][y];
            network.discharge_m3s[cell] = discharge[x][y];
            network.filledelevation_m[cell] = filled[x][y];
            // Flows stay in the legacy units the world stores (plan section 6, B1).
            network.januaryflow[cell] = world.riverjan(x, y);
            network.julyflow[cell] = world.riverjul(x, y);
            inputs.runoff_mm_yr[cell] = runoff[x][y];
            inputs.lakesurface_m[cell] = world.lakesurfacefloat(x, y);
            inputs.januarytemperature_c[cell] = world.jantempfloat(x, y);
            inputs.julytemperature_c[cell] = world.jultempfloat(x, y);
            inputs.januaryprecipitation_mm[cell] = world.janrainfloat(x, y);
            inputs.julyprecipitation_mm[cell] = world.julrainfloat(x, y);
            inputs.icethickness_m[cell] = ice[x][y];
            inputs.climate[cell] = static_cast<int16_t>(world.climate(x, y));
        }
    }
    return inputs;
}

landscape::RegionalParameters regionalparameters(std::uint64_t seed)
{
    namespace regional = tuning::landscape::regional;
    namespace streampower = tuning::landscape::stream_power;
    landscape::RegionalParameters parameters;
    parameters.noiseamplitude_m = regional::noiseAmplitudeM;
    parameters.noiserelieffraction = regional::noiseReliefFraction;
    parameters.noisepersistence = regional::noisePersistence;
    parameters.coastband_m = regional::coastBandM;
    parameters.meanpasses = regional::meanPasses;
    parameters.channels.pinarea_km2 = regional::pinAreaKm2;
    parameters.channels.meanderamplitude = regional::meanderAmplitude;
    parameters.channels.seed = deterministiccontextseed(static_cast<long>(seed), channelseedsalt);
    parameters.relief.erodibility = streampower::erodibility;
    parameters.relief.areaexponent = streampower::areaExponent;
    parameters.relief.elapsed_myr = regional::elapsedMyr;
    parameters.relief.thresholdslope = regional::thresholdSlope;
    parameters.relief.maxiterations = regional::maxIterations;
    parameters.relief.damping = regional::damping;
    return parameters;
}
