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
    // Phase R2 (R6): the erosion run's and the engine's fields for R7 and R9.
    inputs.upliftrate_m_myr.resize(size);
    inputs.glaciation.resize(size);
    inputs.denudation_m.resize(size);
    inputs.island.resize(size);
    inputs.mountainisland.resize(size);
    inputs.erosionduration_myr = static_cast<double>(world.erosionduration());
    // Phase E5 (E5.0): the run's ELA, and the present climate's ELA offset cooling_c over
    // the lapse rate (landscape.glacial.cooling_c is the run's when the region command
    // restored it from the world bundle's manifest, regional_command.cpp).
    inputs.equilibriumline_m.resize(size);
    inputs.presentelaoffset_m = inputs.lapserate_c_per_km > 0.0f
                                    ? static_cast<double>(tuning::landscape::glacial::coolingC) * 1000.0 /
                                          static_cast<double>(inputs.lapserate_c_per_km)
                                    : 0.0;

    const worldgrid<float>& equilibriumline = world.equilibriumline();
    const worldgrid<float>& glaciation = world.glaciation();
    const worldgrid<float>& denudation = world.denudation();
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
            inputs.upliftrate_m_myr[cell] = world.tectonicrockupliftrate(x, y);
            inputs.glaciation[cell] = glaciation[x][y];
            inputs.denudation_m[cell] = denudation[x][y];
            inputs.island[cell] = world.island(x, y) ? 1 : 0;
            inputs.mountainisland[cell] = world.mountainisland(x, y) ? 1 : 0;
            inputs.equilibriumline_m[cell] = equilibriumline[x][y];
        }
    }
    // A world without the glacial term (landscape.glacial.enabled false, or no erosion
    // stage) carries an ELA field of zeros: passed as present it would ice every land
    // cell above the datum's 881 m; left empty, every glacial hook is inert (the review of
    // E5.0).
    bool anyela = false;
    for (float value : inputs.equilibriumline_m)
        anyela = anyela || value != 0.0f;
    if (!anyela)
        inputs.equilibriumline_m.clear();

    // Multiscale plan S3.6: the sub-grid statistics of terrain.subgrid_statistics. A world
    // without the diagnoser (another recipe, a bundle written before S3.6) carries zeros:
    // each field is passed only when some cell is not zero, and is otherwise left empty
    // (absent: the pass's fallbacks apply).
    const worldgrid<float>& relief = world.reliefamplitude();
    const worldgrid<float>& spacing = world.valleyspacing();
    inputs.reliefamplitude_m.resize(size);
    inputs.valleyspacing_km.resize(size);
    bool anyrelief = false;
    bool anyspacing = false;
    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(columns) + static_cast<std::size_t>(x);
            inputs.reliefamplitude_m[cell] = relief[x][y];
            inputs.valleyspacing_km[cell] = spacing[x][y];
            anyrelief = anyrelief || relief[x][y] != 0.0f;
            anyspacing = anyspacing || spacing[x][y] != 0.0f;
        }
    }
    if (!anyrelief)
        inputs.reliefamplitude_m.clear();
    if (!anyspacing)
        inputs.valleyspacing_km.clear();
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
    parameters.noiserelieffield = regional::noiseReliefField;
    parameters.coastband_m = regional::coastBandM;
    parameters.meanpasses = regional::meanPasses;
    parameters.warp.amplitude_cells = regional::warpAmplitudeCells;
    parameters.warp.wavelength_cells = regional::warpWavelengthCells;
    parameters.warp.fade_m = regional::warpFadeM;
    parameters.coast.noiseamplitude_m = regional::coastNoiseAmplitudeM;
    parameters.coast.noisepersistence = regional::coastNoisePersistence;
    parameters.coast.noiserelieffraction = regional::coastNoiseReliefFraction;
    parameters.coast.seadrop_m = regional::coastSeaDropM;
    parameters.coast.landshare = regional::coastLandShare;
    parameters.coast.shorefloor_m = regional::coastShoreFloorM;
    parameters.coast.sharepasses = regional::sharePasses;
    parameters.coast.signedpasses = regional::signedPasses;
    parameters.coast.onesided = regional::coastOneSided;
    parameters.coast.glacialamplitude_m = regional::coastGlacialAmplitudeM;
    parameters.coast.glacialdropfraction = regional::coastGlacialDropFraction;
    parameters.coast.denudationamplitude_m = regional::coastDenudationAmplitudeM;
    parameters.coast.denudationreference_m = regional::coastDenudationReferenceM;
    parameters.coast.islandmargin = regional::islandMargin;
    parameters.edges.smoothgate = regional::gateSmooth;
    parameters.edges.straightmax = regional::straightMaxCells;
    parameters.edges.chunkmin = regional::straightChunkMinCells;
    parameters.edges.chunkmax = regional::straightChunkMaxCells;
    parameters.shores.shorenoise = regional::shoreNoise;
    parameters.shores.terrain_m = regional::shoreTerrainM;
    parameters.valleys.upliftscale = regional::upliftScale;
    parameters.valleys.denudationratescale = regional::denudationRateScale;
    parameters.valleys.floor_km = regional::valleyFloorKm;
    parameters.valleys.reference_m3s = regional::valleyReferenceM3s;
    parameters.valleys.floormax_km = regional::valleyFloorMaxKm;
    parameters.valleys.sideslope = regional::valleySideSlope;
    parameters.valleys.reachcells = regional::valleyReachCells;
    parameters.valleys.reach_km = regional::valleyReachKm;
    parameters.tributaries.enabled = regional::tributaryEnabled;
    parameters.tributaries.spacing_km = regional::tributarySpacingKm;
    parameters.tributaries.hackexponent = regional::tributaryHackExponent;
    parameters.tributaries.hackcoefficient = regional::tributaryHackCoefficient;
    parameters.tributaries.junctionangle_deg = regional::tributaryJunctionAngleDeg;
    parameters.tributaries.anglespread = regional::tributaryAngleSpread;
    parameters.tributaries.nodemargin = regional::tributaryNodeMargin;
    parameters.tributaries.lengthminfraction = regional::tributaryLengthMinFraction;
    parameters.troughs.width_km = regional::troughWidthKm;
    parameters.troughs.reference_km2 = regional::troughReferenceKm2;
    parameters.troughs.widthexponent = regional::troughWidthExponent;
    parameters.troughs.widthmax_km = regional::troughWidthMaxKm;
    parameters.troughs.wall_m = regional::troughWallM;
    parameters.troughs.exponent = regional::troughExponent;
    parameters.troughs.glaciationmin = regional::troughGlaciationMin;
    parameters.cirques.radius_km = regional::cirqueRadiusKm;
    parameters.cirques.depth_m = regional::cirqueDepthM;
    parameters.cirques.band_m = regional::cirqueBandM;
    parameters.cirques.candidates = regional::cirqueCandidates;
    parameters.cirques.glaciationmin = regional::cirqueGlaciationMin;
    parameters.ice.thicknessfraction = regional::iceElaThicknessFraction;
    parameters.ice.thicknessmax_m = regional::iceElaThicknessMaxM;
    parameters.ice.edgenoise_m = regional::iceElaEdgeNoiseM;
    parameters.channels.pinarea_km2 = regional::pinAreaKm2;
    parameters.channels.meanderamplitude = regional::meanderAmplitude;
    parameters.channels.meanderdischargeexponent = regional::meanderDischargeExponent;
    parameters.channels.meanderreference_m3s = regional::meanderReferenceM3s;
    parameters.channels.meanderroughness = regional::meanderRoughness;
    parameters.channels.meanderroughnesspersistence = regional::meanderRoughnessPersistence;
    parameters.channels.junctionjitter = regional::jitterJunction;
    parameters.channels.edgemargin = regional::jitterEdgeMargin;
    parameters.channels.springspercell = regional::springsPerCell;
    parameters.channels.diagonaltransit = regional::diagonalTransit;
    parameters.channels.transittolerance_m = regional::diagonalTransitToleranceM;
    parameters.channels.seed = deterministiccontextseed(static_cast<long>(seed), channelseedsalt);
    parameters.relief.erodibility = streampower::erodibility;
    parameters.relief.areaexponent = streampower::areaExponent;
    parameters.relief.elapsed_myr = regional::elapsedMyr;
    parameters.relief.thresholdslope = regional::thresholdSlope;
    parameters.relief.maxiterations = regional::maxIterations;
    parameters.relief.damping = regional::damping;
    parameters.relief.catchmentreach = regional::catchmentReachGlobalCells;
    return parameters;
}

landscape::RefinementParameters refinementparameters()
{
    namespace refinement = tuning::landscape::refinement;
    landscape::RefinementParameters parameters;
    parameters.factor = refinement::factor;
    parameters.projectionweightchannel = refinement::projectionWeightChannel;
    parameters.storedscalemax = refinement::storedScaleMax;
    parameters.detailcells = refinement::detailCells;
    parameters.tilecells = refinement::tileCells;
    parameters.smoothpasses = refinement::smoothPasses;
    parameters.smoothrelaxation = refinement::smoothRelaxation;
    parameters.smoothbasis = refinement::smoothBasis;
    parameters.prolongationlimit = refinement::prolongationLimit;
    parameters.channelcellslope = refinement::channelCellSlope;
    parameters.parentlift = refinement::parentLift;
    return parameters;
}
