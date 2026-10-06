#include "pipeline/climate_stages.hpp"

#include <vector>

#include "simulations/climate/climate_classification.hpp"
#include "simulations/climate/climate_coupling.hpp"
#include "simulations/climate/climate_fields.hpp"
#include "simulations/climate/climate_world.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/climate/climate_parameters.hpp"
#include "core/procedural_noise.hpp"
#include "pipeline/world_processes.hpp"

using namespace std;

namespace
{
vector<vector<int>> makeintgrid(int value = 0)
{
    return vector<vector<int>>(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, value));
}

void reseedclimatepass(const climateworld& world, int salt)
{
    fast_srand(deterministicfastseed(deterministiccontextseed(world.seed(), salt)));
}

vector<vector<int>> buildclimatefractal(const climateworld& world, int salt, int maximum)
{
    vector<vector<int>> fractal = makeintgrid();
    reseedclimatepass(world, salt);
    createfractal(fractal, world.width(), world.height(), 8, 0.2f, static_cast<float>(random(3, 6)), 1, maximum, 0, 0);
    warp(fractal, world.width(), world.height(), world.maxelevation(), random(20, 80), 1);
    return fractal;
}
}

bool isclimateonlystage(GenerationStageId stage)
{
    return isclimatefeedbackstage(stage) && stage != GenerationStageId::terrain_roughness;
}

void countclimatelandandsea(climateworld& world)
{
    const int width = world.width();
    const int height = world.height();

    int seatotal = 0;
    int landtotal = 0;

    for (int i = 1; i < width; i++) // Avoid the margins, to avoid false positives.
    {
        for (int j = 1; j < height; j++)
        {
            if (!world.basebelowsea(i, j))
                landtotal++;
            else
                seatotal++;
        }
    }

    if (landtotal > 0)
    {
        if ((float)seatotal / (float)landtotal < 0.05) // If it's this close to 0, assume it is 0.
            seatotal = 0;
    }

    climateinputs inputs = world.inputs();
    inputs.landtotal = landtotal;
    inputs.seatotal = seatotal;
    world.bindinputs(inputs);
}

void executeclimatestage(GenerationStageId stage, climateworld& world, const ClimateStageContext& context)
{
    seedprocessstream(world.seed(), stage);
    const int maxelev = world.maxelevation();

    switch (stage)
    {
    case GenerationStageId::global_temperature:
    {
        vector<vector<int>> fractal = makeintgrid();
        reseedclimatepass(world, 0x5002);
        createfractal(fractal, world.width(), world.height(), 8, 0.2f, static_cast<float>(random(1, 4)), 1, maxelev, 0, 0);
        warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
        createtemperaturemap(world, fractal);
        return;
    }

    case GenerationStageId::ocean_currents:
        if (world.seatotal() > 0)
            createoceancurrentmap(world);
        return;

    case GenerationStageId::sea_surface_temperatures:
        if (world.seatotal() > 0)
            createsurfacetemperaturemap(world);
        return;

    case GenerationStageId::pressure:
        if (world.seatotal() > 0)
            createpressuremap(world);
        return;

    case GenerationStageId::winds:
        if (world.seatotal() > 0)
        {
            createvectorwindmap(world);
            updatehorsebeltsfrompressure(world);
            world.updatelevels();
        }
        return;

    case GenerationStageId::sea_ice_and_tides:
        if (world.seatotal() > 0)
        {
            vector<vector<int>> fractal = buildclimatefractal(world, 0x5003, maxelev);
            createseaicemap(world, fractal);
        }
        return;

    case GenerationStageId::rainfall:
    {
        vector<vector<int>> fractal = makeintgrid();
        reseedclimatepass(world, 0x5004);
        createfractal(fractal, world.width(), world.height(), 8, 0.2f, 0.6f, 1, maxelev, 0, 0);
        warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
        warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
        countclimatelandandsea(world);
        createrainmap(world, fractal, world.landtotal(), world.seatotal(), context.smalllake, context.landshape);
        if (world.seatotal() > 0 && tuning::climate::circulation::enableLaggedDiabaticCoupling)
        {
            vector<vector<int>> couplingfractal = buildclimatefractal(world, 0x5004, maxelev);
            convergeclimatecoupling(world, couplingfractal);
            world.updatelevels();
        }
        return;
    }

    case GenerationStageId::mountain_temperature_lapse:
        checkpoleclimates(world);
        world.syncseasonalclimatefromlegacy();
        return;

    case GenerationStageId::climates_and_biomes:
        createclimatemap(world);
        createbiomemap(world);
        return;

    default:
        return;
    }
}
