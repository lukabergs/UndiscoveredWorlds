#include "simulations/landscape/regional_cube.hpp"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/window_grid.hpp"
#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/level_blocks.hpp"
#include "simulations/landscape/refinement_projection.hpp"
#include "simulations/landscape/regional_cirques.hpp"
#include "simulations/landscape/regional_coast.hpp"
#include "simulations/landscape/regional_folds.hpp"
#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_ice.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_troughs.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/regional_water.hpp"
#include "simulations/landscape/upsampling.hpp"

// T8 W2-W2f (regional_cube.hpp): the cube windows' drivers, the step functions of the lat-lon hook
// interface (regional_pass.hpp) rewritten on window grids, each calling its unit's cube form.

namespace landscape
{
namespace
{
constexpr double secondsperyear = 365.25 * 86400.0;

// The parent level of a cube level (level >= 2): its kept cells and their fields after its projection.
struct CubeParent
{
    std::shared_ptr<const mesh::CellGrid> kept;
    std::vector<float> elevation_m, lift_m, solved_m, lakesurface_m;
    std::vector<uint8_t> sea;
};

// What a cube level's steps read (LevelContext's analogue).
struct CubeLevel
{
    const RegionalGlobalInputs& inputs;
    const mesh::CellGrid& grid;           // the level's domain
    const std::vector<uint8_t>& kept;     // one flag per domain cell
    const RegionalParameters& parameters;
    const FineWarp& warp;
    const GlobalIncision& incision;
    const RegionalDetails& details;
    int level = 0;                        // 1 .. L in the loop, 0 in the single pass
    int factor = 0;                       // 0 in the single pass
    int k = 1;                            // level cells per level-0 cell edge
    const RefinementParameters* refinement = nullptr;
    const mesh::CellGrid* projectionparent = nullptr;     // level 1: level 0's kept cells; deeper: the parent's
    const std::vector<float>* projectionfield = nullptr;
    const CubeParent* parent = nullptr;                   // level >= 2
};

// The level-0 ancestor's id of every cell of a window at k.
std::vector<std::size_t> ancestors(const RegionalGlobalInputs& inputs, const mesh::CellGrid& grid, int k)
{
    std::vector<std::size_t> result(grid.cellcount());
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
        result[c] = inputs.grid->localofkey(parentkey(grid, grid.globalkey(c), k));
    return result;
}

// The descendants at k of level-0 keys (ascending).
std::vector<mesh::LevelKey> descendants(int facecells, int k, const std::vector<mesh::LevelKey>& keys)
{
    std::vector<mesh::LevelKey> result;
    result.reserve(keys.size() * static_cast<std::size_t>(k) * static_cast<std::size_t>(k));
    for (const mesh::LevelKey key : keys)
        for (const mesh::LevelKey child : mesh::cubechildkeys(facecells, k, key))
            result.push_back(child);
    std::sort(result.begin(), result.end());
    return result;
}

int keptchanged(const CubeLevel& c, const std::vector<float>& before, const std::vector<float>& after)
{
    int changed = 0;
    for (std::size_t cell = 0; cell < after.size(); cell++)
        changed += c.kept[cell] != 0 && before[cell] != after[cell] ? 1 : 0;
    return changed;
}

double channelcap(const CubeLevel& c)
{
    if (c.refinement == nullptr || !(c.refinement->channelcellslope > 0.0))
        return 0.0;
    return c.refinement->channelcellslope * northsouthspacing_m(c.grid);
}

void band(const CubeLevel& c, double& lower, double& upper)
{
    lower = 0.0;
    upper = 0.0;
    if (c.details.factor < 2)
        return;
    lower = c.details.detailcells / 2.0;
    upper = c.level == 1 ? std::numeric_limits<double>::infinity() : c.details.detailcells * c.details.factor / 2.0;
}

double nearestrank(std::vector<double>& values, double fraction)
{
    if (values.empty())
        return 0.0;
    const std::size_t rank = static_cast<std::size_t>(std::ceil(fraction * static_cast<double>(values.size())));
    const std::size_t index = std::min(values.size() - 1, rank > 0 ? rank - 1 : 0);
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
    return values[index];
}

// --- The steps (regional_pass.hpp's hook interface, in LevelStep order). ---

void noisestep(const CubeLevel& c, LevelState& state)
{
    int parentscale = 1;
    if (c.parent == nullptr)
    {
        state.start = regionalstartsurface(c.inputs, c.grid, c.warp, c.parameters, c.refinement != nullptr, &c.kept);
        if (state.diagnostics != nullptr)
        {
            state.diagnostics->pockets = state.start.pockets;
            state.diagnostics->pocketcells = state.start.pocketcells;
            state.diagnostics->pocketsfilled = state.start.pocketsfilled;
            state.diagnostics->pocketcellsfilled = state.start.pocketcellsfilled;
        }
    }
    else
    {
        parentscale = static_cast<int>(levelscale(c.factor, c.level - 1));
        addlevelnoise(c.inputs, c.grid, parentscale, c.parameters, state.start.start_m);
    }
    state.noiseband_m = levelnoiseband(c.inputs, c.grid, parentscale, c.parameters);
    recordlevelstep(state, LevelStep::noiseoctave);
}

void presmoothingstep(const CubeLevel& c, LevelState& state)
{
    StartSurface& start = state.start;
    if (c.parent != nullptr)
    {
        const CubeParent& parent = *c.parent;
        correctparentmeans(*parent.kept, parent.elevation_m, c.grid, c.factor, c.parameters.meanpasses, start.start_m);
        start.blockmeanerror_m = levelparentmeanerror(c.grid, start.start_m, *parent.kept, parent.elevation_m, c.factor, c.kept);
        start.lapsebase_m = upsample(c.inputs.elevation_m, *c.inputs.grid, c.grid);
    }
    if (state.diagnostics != nullptr)
        state.diagnostics->blockmeanerror_m = start.blockmeanerror_m;
    recordlevelstep(state, LevelStep::presmoothing);
    const std::size_t count = c.grid.cellcount();
    if (start.start_m.size() != count || start.lapsebase_m.size() != count)
        throw std::invalid_argument("cube level: the start surface must cover the window's domain");
    state.blockmeanerror_m = start.blockmeanerror_m;
    state.sharepasses = start.sharepasses;
    state.sharemisses = start.sharemisses;
    state.lapsebase_m = std::move(start.lapsebase_m);
    state.surface_m = std::move(start.start_m);
}

void coaststep(const CubeLevel& c, LevelState& state)
{
    const CoastParameters& coast = c.parameters.coast;
    if (c.parent != nullptr && coast.landshare > 0.0 && (coast.warpcells > 0.0 || coast.pocketreach > 0))
    {
        PocketCount filled, left;
        state.sea = refinedcoast(c.inputs, c.grid, c.factor, *c.parent->kept, c.parent->sea, c.parent->elevation_m, c.warp,
                                 c.parameters, c.details.parentseagate, state.surface_m, &c.kept, &filled, &left);
        if (state.diagnostics != nullptr && coast.pocketreach > 0)
        {
            state.diagnostics->pocketsfilled = filled.pockets;
            state.diagnostics->pocketcellsfilled = filled.pocketcells;
            state.diagnostics->pockets = left.pockets;
            state.diagnostics->pocketcells = left.pocketcells;
        }
    }
    else
        state.sea = regionalsea(c.inputs, c.grid, c.warp, c.parameters, state.surface_m, c.details.parentseagate,
                                c.refinement != nullptr);
    recordlevelstep(state, LevelStep::coast);
}

void incisionstep(const CubeLevel& c, LevelState& state)
{
    if (!c.incision.active)
        return;
    state.incisednetwork = c.inputs.network;
    state.incisednetwork.filledelevation_m = c.incision.filledelevation_m;
    const std::vector<float> fineincision = upsample(c.incision.incision_m, *c.inputs.grid, c.grid);
    state.bedsurface_m = state.surface_m;
    for (std::size_t cell = 0; cell < state.bedsurface_m.size(); cell++)
        state.bedsurface_m[cell] =
            static_cast<float>(static_cast<double>(state.surface_m[cell]) - std::max(0.0f, fineincision[cell]));
    recordlevelstep(state, LevelStep::incision);
}

void channelsstep(const CubeLevel& c, LevelState& state)
{
    const RegionalGlobalInputs& inputs = c.inputs;
    const bool incised = c.incision.active;
    ChannelParameters channelparameters = c.parameters.channels;
    channelparameters.pinsinks = c.refinement != nullptr && c.parameters.relief.catchmentreach > 0 ? 1 : 0;
    std::vector<float> relief;
    if (channelparameters.meanderroughness > 0.0 && channelparameters.meanderroughnessreliefgain > 0.0)
        relief = regionalnoiserelief(inputs, c.parameters);
    state.channels = rasterisechannels(*inputs.grid, c.grid, incised ? state.incisednetwork : inputs.network,
                                       channelparameters, inputs.sealevel, incised ? state.bedsurface_m : state.surface_m,
                                       state.sea, relief, &c.kept);
    const FineChannels& channels = state.channels;
    for (std::size_t cell = 0; cell < state.surface_m.size(); cell++)
        if (channels.channel[cell] != 0 && state.sea[cell] == 0)
            state.surface_m[cell] = channels.bed_m[cell];
    state.bed_m = channels.bed_m;
    state.sinkjunctions = channels.sinkjunctions;
    state.diagonallinks = channels.diagonallinks;
    state.cornerlinks = channels.cornerlinks;
    state.mergedpieces = channels.mergedpieces;
    state.drytails = channels.drytails;
    recordlevelstep(state, LevelStep::channels);
}

void solvestep(const CubeLevel& c, LevelState& state)
{
    const RegionalGlobalInputs& inputs = c.inputs;
    const mesh::CellGrid& grid = c.grid;
    const std::size_t count = grid.cellcount();
    const std::int64_t levelcells = grid.facecells();
    const std::vector<uint8_t>& sea = state.sea;
    const FineChannels& channels = state.channels;
    const std::vector<float>& surface = state.surface_m;

    // The domain's outer ring: a cell with a table slot missing that the full level has (not a cube
    // corner's missing diagonal).
    std::vector<uint8_t> ring(count, 0);
    for (mesh::CellId cell = 0; cell < count; cell++)
        for (int slot = 1; slot <= mesh::slotcount && ring[cell] == 0; slot++)
            if (grid.neighbour(cell, slot) < 0 && mesh::cubekeyneighbour(levelcells, grid.globalkey(cell), slot) != mesh::nokey)
                ring[cell] = 1;

    const std::vector<float> runoff = upsample(inputs.runoff_mm_yr, *inputs.grid, grid);
    AnalyticalReliefInput solveinput;
    solveinput.start_m = surface;
    solveinput.pinned.assign(count, 0);
    solveinput.runoff_m_yr.resize(count);
    solveinput.key.resize(count);
    for (mesh::CellId cell = 0; cell < count; cell++)
    {
        solveinput.pinned[cell] = channels.channel[cell] != 0 || sea[cell] != 0 || ring[cell] != 0 ? 1 : 0;
        if (sea[cell] != 0)
            solveinput.start_m[cell] = inputs.sealevel;
        solveinput.runoff_m_yr[cell] = std::max(0.0f, runoff[cell]) / 1000.0f;
        solveinput.key[cell] = grid.globalkey(cell);
    }
    AnalyticalReliefParameters relief = c.parameters.relief;
    if (c.refinement == nullptr)
        relief.catchmentreach = 0;
    if (relief.catchmentreach > 0)
    {
        const std::vector<std::size_t> global = ancestors(inputs, grid, c.k);
        solveinput.globalcell.resize(count);
        solveinput.globalgrid = inputs.grid.get();
        for (std::size_t cell = 0; cell < count; cell++)
            solveinput.globalcell[cell] = static_cast<mesh::CellId>(global[cell]);
    }
    if (c.incision.active)
    {
        solveinput.uplift_m_myr = upsample(c.incision.rate_m_myr, *inputs.grid, grid);
        for (float& rate : solveinput.uplift_m_myr)
            rate = std::max(0.0f, rate);
    }

    // S3.2: the solve against its parent (level >= 2).
    const RefinementParameters* refinement = c.refinement;
    std::vector<float> added;
    if (refinement != nullptr && refinement->parentlift != 0 && c.parent != nullptr)
    {
        const CubeParent& parent = *c.parent;
        std::vector<float> parentlift;
        if (refinement->parentlift == 1)
            parentlift = parent.lift_m;
        else
        {
            parentlift.resize(parent.elevation_m.size());
            for (std::size_t cell = 0; cell < parentlift.size(); cell++)
                parentlift[cell] = static_cast<float>(static_cast<double>(parent.elevation_m[cell]) -
                                                      static_cast<double>(parent.solved_m[cell]));
        }
        added = prolong(*parent.kept, parentlift, grid, c.factor, refinement->prolongationlimit == 1).values;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (solveinput.pinned[cell] == 0)
                solveinput.start_m[cell] =
                    static_cast<float>(static_cast<double>(solveinput.start_m[cell]) - static_cast<double>(added[cell]));
            else
                added[cell] = 0.0f;
        }
    }

    AnalyticalReliefResult solved = solveanalyticalrelief(grid, solveinput, relief);
    state.reliefconverged = solved.converged;
    state.changedreceivers = solved.changedreceivers;
    state.unreachedcells = 0;
    state.ringbasincells = 0;
    if (c.refinement != nullptr)
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (c.kept[cell] == 0)
                continue;
            if (!solved.unreached.empty() && solved.unreached[cell] != 0)
                state.unreachedcells++;
            const std::int32_t outlet = solved.outlet[cell];
            if (outlet >= 0 && ring[static_cast<std::size_t>(outlet)] != 0)
                state.ringbasincells++;
        }
    std::vector<float>& elevation = state.elevation_m;
    elevation = surface;
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    for (mesh::CellId cell = 0; cell < count; cell++)
    {
        if (solveinput.pinned[cell] != 0)
            continue;
        const float value = solved.elevation_m[cell];
        if (added.empty())
        {
            elevation[cell] = value;
            continue;
        }
        const CubeParent& parent = *c.parent;
        const float parentvalue = parent.elevation_m[parent.kept->localofkey(parentkey(grid, grid.globalkey(cell), c.factor))];
        elevation[cell] = std::max(static_cast<float>(static_cast<double>(value) + static_cast<double>(added[cell])),
                                   std::min({ value, abovesea, parentvalue }));
    }
    {
        double total = 0.0, largest = 0.0;
        std::size_t cells = 0;
        std::vector<double> lowerings;
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (c.kept[cell] == 0 || solveinput.pinned[cell] != 0)
                continue;
            const double lowering = static_cast<double>(surface[cell]) - static_cast<double>(elevation[cell]);
            total += lowering;
            largest = std::max(largest, lowering);
            cells++;
            if (state.diagnostics != nullptr)
                lowerings.push_back(lowering);
        }
        state.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        state.maxlowering_m = largest;
        if (state.diagnostics != nullptr && !lowerings.empty())
        {
            std::sort(lowerings.begin(), lowerings.end());
            const auto rank = [&](double share) {
                const std::size_t index = static_cast<std::size_t>(std::ceil(share * static_cast<double>(lowerings.size())));
                return lowerings[std::min(lowerings.size(), std::max<std::size_t>(index, 1)) - 1];
            };
            state.diagnostics->solvelowering50_m = rank(0.5);
            state.diagnostics->solvelowering99_m = rank(0.99);
        }
    }
    state.solveddischarge_m3_yr = std::move(solved.discharge_m3_yr);
    state.pinned = std::move(solveinput.pinned);

    if (refinement != nullptr && refinement->channelcellslope > 0.0)
    {
        // A land channel cell's cell mean b + min(h, theta c_l) over its table neighbours.
        const double cap = refinement->channelcellslope * northsouthspacing_m(grid);
        const std::vector<float> field = elevation;
        for (mesh::CellId cell = 0; cell < count; cell++)
        {
            if (channels.channel[cell] == 0 || sea[cell] != 0)
                continue;
            const double bed = static_cast<double>(state.bed_m[cell]);
            double sum = 0.0;
            int neighbours = 0;
            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const std::int32_t next = grid.neighbour(cell, slot);
                if (next < 0)
                    continue;
                const std::size_t other = static_cast<std::size_t>(next);
                if (channels.channel[other] != 0 || sea[other] != 0)
                    continue;
                sum += static_cast<double>(field[other]) - bed;
                neighbours++;
            }
            const double height_m = neighbours > 0 ? std::max(0.0, sum / neighbours) : 0.0;
            elevation[cell] = static_cast<float>(bed + std::min(height_m, cap));
        }
    }
    state.solved_m = elevation;
    if (!added.empty())
        for (std::size_t cell = 0; cell < count; cell++)
            if (state.pinned[cell] == 0)
                state.solved_m[cell] = solved.elevation_m[cell];
    recordlevelstep(state, LevelStep::solve);
}

void foldstep(const CubeLevel& c, LevelState& state)
{
    const RegionalGlobalInputs& inputs = c.inputs;
    if (!(c.parameters.folds.amplitudemax_m > 0.0) || inputs.foldstrain.size() != inputs.grid->cellcount())
        return;
    const int parentscale = c.parent != nullptr ? static_cast<int>(levelscale(c.factor, c.level - 1)) : 0;
    const bool record = state.diagnostics != nullptr;
    std::vector<float> before = state.elevation_m;
    std::vector<uint8_t> floorclamped;
    addfoldfabric(inputs, c.grid, c.parameters, state.channels, state.sea, state.elevation_m, parentscale,
                  record ? &floorclamped : nullptr);
    if (record)
    {
        std::vector<double> changes;
        for (std::size_t cell = 0; cell < before.size(); cell++)
        {
            if (c.kept[cell] == 0)
                continue;
            if (!floorclamped.empty() && floorclamped[cell] != 0)
                state.diagnostics->foldfloorclamped++;
            if (before[cell] != state.elevation_m[cell])
                changes.push_back(std::fabs(static_cast<double>(state.elevation_m[cell]) - static_cast<double>(before[cell])));
        }
        state.diagnostics->foldcells = static_cast<int>(changes.size());
        state.diagnostics->foldamplitude_p50_m = nearestrank(changes, 0.5);
        state.diagnostics->foldamplitude_p99_m = nearestrank(changes, 0.99);
    }
    channelcellmeans(c.grid, state.channels, state.sea, channelcap(c), before, state.elevation_m);
    state.folded_m = state.elevation_m;
    state.unfolded_m = std::move(before);
    recordlevelstep(state, LevelStep::foldfabric);
}

void valleystep(const CubeLevel& c, LevelState& state)
{
    if (!c.details.valleys)
        return;
    double lower = 0.0, upper = 0.0;
    band(c, lower, upper);
    const double cap = channelcap(c);
    std::vector<float> before;
    if (state.diagnostics != nullptr || cap > 0.0)
        before = state.elevation_m;
    carvevalleys(c.inputs, c.grid, c.parameters, state.channels, state.sea, state.elevation_m, lower, upper);
    channelcellmeans(c.grid, state.channels, state.sea, cap, before, state.elevation_m);
    if (state.diagnostics != nullptr)
        state.diagnostics->valleycells = keptchanged(c, before, state.elevation_m);
    recordlevelstep(state, LevelStep::valleyfloors);
}

void glacialstep(const CubeLevel& c, LevelState& state)
{
    double lower = 0.0, upper = 0.0;
    band(c, lower, upper);
    const double cap = channelcap(c);
    state.glacial = regionalglacialfields(c.inputs, c.grid, c.parameters);
    const bool record = state.diagnostics != nullptr;
    const bool keep = record || cap > 0.0;
    std::vector<float> before;
    if (c.details.troughs)
    {
        if (keep)
            before = state.elevation_m;
        carvetroughs(c.inputs, c.grid, c.parameters, state.glacial, state.channels, state.sea, state.elevation_m, lower, upper);
        channelcellmeans(c.grid, state.channels, state.sea, cap, before, state.elevation_m);
        if (record)
            state.diagnostics->troughcells = keptchanged(c, before, state.elevation_m);
    }
    if (c.details.cirques)
    {
        if (keep)
            before = state.elevation_m;
        carvecirques(c.inputs, c.grid, c.parameters, state.glacial, state.channels, state.sea, state.elevation_m, lower, upper);
        channelcellmeans(c.grid, state.channels, state.sea, cap, before, state.elevation_m);
        if (record)
            state.diagnostics->cirquecells = keptchanged(c, before, state.elevation_m);
    }
    if (c.details.troughs || c.details.cirques)
        recordlevelstep(state, LevelStep::glacialforms);
}

void lakesstep(const CubeLevel& c, LevelState& state)
{
    LakeInheritance lakes = c.details.lakes;
    if (lakes.flood)
    {
        lakes.startsurface_m = &state.surface_m;
        if (state.noiseband_m.size() == state.surface_m.size())
            lakes.rimcap_m = &state.noiseband_m;
        lakes.parentseagate = c.details.parentseagate;
    }
    std::vector<float>& elevation = state.elevation_m;
    const std::size_t count = elevation.size();
    const bool withbeds = state.bed_m.size() == count && state.channels.channel.size() == count && state.sea.size() == count;
    std::vector<float> means;
    if (withbeds)
    {
        means = elevation;
        for (std::size_t cell = 0; cell < count; cell++)
            if (state.channels.channel[cell] != 0 && state.sea[cell] == 0)
                elevation[cell] = state.bed_m[cell];
    }
    state.water = regionalwater(c.inputs, c.grid, c.warp, c.parameters, state.sea, state.channels, elevation,
                                c.details.globallakedepth, lakes);
    if (withbeds)
    {
        const std::vector<std::size_t> dried = channelmeansafterlakes(c.inputs, c.grid, state.channels, state.sea, state.bed_m,
                                                                      means, elevation, state.water, c.details.coarselevel);
        if (state.diagnostics != nullptr && c.refinement != nullptr)
            for (const std::size_t cell : dried)
                state.diagnostics->driedchannels += c.kept[cell] != 0 ? 1 : 0;
    }
    recordlevelstep(state, LevelStep::lakes);
}

void smoothingstep(const CubeLevel& c, LevelState& state, const std::vector<std::size_t>& global)
{
    if (c.refinement == nullptr)
        return;
    const RegionalGlobalInputs& inputs = c.inputs;
    const RefinementParameters& refinement = *c.refinement;
    const std::size_t count = c.grid.cellcount();
    std::vector<float>& elevation = state.elevation_m;
    const std::vector<float>& lakesurface = state.water.lakesurface_m;
    state.ownlake_m.assign(count, 0.0f);
    for (std::size_t cell = 0; cell < count; cell++)
        if (inputs.lakesurface_m[global[cell]] > 0.0f && inputs.network.sea[global[cell]] == 0)
            state.ownlake_m[cell] = inputs.lakesurface_m[global[cell]];
    static const std::vector<float> nobeds;
    const std::vector<float>& beds = refinement.channelcellslope > 0.0 ? state.bed_m : nobeds;
    state.bounds = projectionbounds(c.grid, state.channels.channel, state.sea, lakesurface, elevation, inputs.sealevel,
                                    refinement.projectionweightchannel, state.ownlake_m, beds);
    state.lift_m = elevation;
    const SmoothingReport smoothing =
        smoothparentmeans(c.grid, elevation, state.bounds, *c.projectionparent, *c.projectionfield, c.factor,
                          refinement.smoothpasses, refinement.smoothrelaxation, refinement.smoothbasis, c.kept);
    if (state.diagnostics != nullptr)
    {
        state.diagnostics->smoothing = smoothing;
        std::vector<double> moves;
        for (std::size_t cell = 0; cell < count; cell++)
            if (c.kept[cell] != 0)
                moves.push_back(std::fabs(static_cast<double>(elevation[cell]) - static_cast<double>(state.lift_m[cell])));
        if (!moves.empty())
        {
            const std::size_t rank = static_cast<std::size_t>(std::ceil(0.99 * static_cast<double>(moves.size()))) - 1;
            std::nth_element(moves.begin(), moves.begin() + static_cast<std::ptrdiff_t>(rank), moves.end());
            state.diagnostics->smoothmove99_m = moves[rank];
        }
    }
    if (refinement.smoothpasses > 0)
    {
        recordlevelstep(state, LevelStep::smoothing);
        state.bounds = projectionbounds(c.grid, state.channels.channel, state.sea, lakesurface, elevation, inputs.sealevel,
                                        refinement.projectionweightchannel, state.ownlake_m, beds);
    }
}

void projectionstep(const CubeLevel& c, LevelState& state)
{
    if (c.refinement == nullptr)
        return;
    const std::size_t count = c.grid.cellcount();
    std::vector<float>& elevation = state.elevation_m;
    std::vector<float>& lakesurface = state.water.lakesurface_m;
    std::vector<uint8_t> wet(count, 0);
    for (std::size_t cell = 0; cell < count; cell++)
        wet[cell] = lakesurface[cell] > 0.0f ? 1 : 0;
    const ProjectionReport projection =
        projectlevel(c.grid, elevation, state.bounds, *c.projectionparent, *c.projectionfield, c.factor, c.kept);
    if (state.diagnostics != nullptr)
        state.diagnostics->projection = projection;
    recordlevelstep(state, LevelStep::projection);
    const bool channelbounds = !state.bounds.bed.empty() && c.details.coarselevel;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (wet[cell] == 0 || elevation[cell] <= lakesurface[cell])
            continue;
        const bool channelcell = state.channels.channel[cell] != 0;
        if (channelcell && !(channelbounds && state.bounds.bed[cell] != 0))
            continue;
        lakesurface[cell] = 0.0f;
        if (c.kept[cell] != 0 && state.diagnostics != nullptr)
            (channelcell ? state.diagnostics->driedchannels : state.diagnostics->driedcells)++;
    }
    for (std::size_t cell = 0; cell < count; cell++)
        state.lift_m[cell] = static_cast<float>(static_cast<double>(elevation[cell]) - static_cast<double>(state.lift_m[cell]));
}

void levelsteps(const CubeLevel& c, LevelState& state, const std::vector<std::size_t>& global)
{
    noisestep(c, state);
    presmoothingstep(c, state);
    coaststep(c, state);
    incisionstep(c, state);
    channelsstep(c, state);
    solvestep(c, state);
    foldstep(c, state);
    valleystep(c, state);
    glacialstep(c, state);
    unfoldforlakes(state);
    lakesstep(c, state);
    refoldafterlakes(c.grid, state);
    smoothingstep(c, state, global);
    projectionstep(c, state);
}

// finishregionallevel on a cube window: the ice above the ELA, the rivers and the lapse-rate climate, the
// kept cells' fields in ascending key order.
RegionalFields finish(const CubeLevel& c, LevelState& level, const std::vector<std::size_t>& global)
{
    const RegionalGlobalInputs& inputs = c.inputs;
    const mesh::CellGrid& grid = c.grid;
    const std::size_t count = grid.cellcount();
    const double sealevel = static_cast<double>(inputs.sealevel);
    const std::vector<float>& base = level.lapsebase_m;
    const std::vector<float>& elevation = level.elevation_m;
    const FineChannels& channels = level.channels;
    iceabovetheela(inputs, grid, c.parameters, level.glacial, level.sea, elevation, level.water.icethickness_m);

    RegionalFields fields;
    fields.blockmeanerror_m = level.blockmeanerror_m;
    fields.sharepasses = level.sharepasses;
    fields.sharemisses = level.sharemisses;
    fields.reliefconverged = level.reliefconverged;
    fields.changedreceivers = level.changedreceivers;
    fields.meanlowering_m = level.meanlowering_m;
    fields.maxlowering_m = level.maxlowering_m;
    const mesh::CellGrid& levelzero = *inputs.grid;
    const std::vector<float> januarytemperature = upsample(inputs.januarytemperature_c, levelzero, grid);
    const std::vector<float> julytemperature = upsample(inputs.julytemperature_c, levelzero, grid);
    const std::vector<float> januaryprecipitation = upsample(inputs.januaryprecipitation_mm, levelzero, grid);
    const std::vector<float> julyprecipitation = upsample(inputs.julyprecipitation_mm, levelzero, grid);
    const double lapse = static_cast<double>(inputs.lapserate_c_per_km) / 1000.0;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (c.kept[cell] == 0)
            continue;
        if (channels.channel[cell] != 0)
        {
            fields.riverdirection.push_back(channels.direction[cell]);
            fields.discharge_m3s.push_back(channels.discharge_m3s[cell]);
            fields.januaryflow.push_back(channels.januaryflow[cell]);
            fields.julyflow.push_back(channels.julyflow[cell]);
        }
        else
        {
            fields.riverdirection.push_back(0);
            fields.discharge_m3s.push_back(static_cast<float>(level.solveddischarge_m3_yr[cell] / secondsperyear));
            fields.januaryflow.push_back(0);
            fields.julyflow.push_back(0);
        }
        const double lift = std::max(0.0, static_cast<double>(base[cell]) - sealevel) -
                            std::max(0.0, static_cast<double>(elevation[cell]) - sealevel);
        fields.januarytemperature_c.push_back(static_cast<float>(static_cast<double>(januarytemperature[cell]) + lapse * lift));
        fields.julytemperature_c.push_back(static_cast<float>(static_cast<double>(julytemperature[cell]) + lapse * lift));
        fields.januaryprecipitation_mm.push_back(std::max(0.0f, januaryprecipitation[cell]));
        fields.julyprecipitation_mm.push_back(std::max(0.0f, julyprecipitation[cell]));
        fields.climate.push_back(inputs.climate[global[cell]]);
        fields.elevation_m.push_back(elevation[cell]);
        fields.sea.push_back(level.sea[cell]);
        fields.lakesurface_m.push_back(level.water.lakesurface_m[cell]);
        fields.icethickness_m.push_back(level.water.icethickness_m[cell]);
    }
    return fields;
}

// The kept flags of a domain whose kept box is the descendants of sorted level-0 keys.
std::vector<uint8_t> keptflags(const std::vector<std::size_t>& global, const RegionalGlobalInputs& inputs,
                               const std::vector<mesh::LevelKey>& keptglobal)
{
    std::vector<uint8_t> kept(global.size(), 0);
    for (std::size_t cell = 0; cell < global.size(); cell++)
        kept[cell] = std::binary_search(keptglobal.begin(), keptglobal.end(), inputs.grid->globalkey(static_cast<mesh::CellId>(global[cell])))
                         ? 1
                         : 0;
    return kept;
}

void checkcube(const RegionalGlobalInputs& inputs, const std::vector<mesh::LevelKey>& kept0, const RegionalParameters& parameters,
               const char* name)
{
    if (inputs.grid == nullptr || inputs.grid->kind() != mesh::GridKind::cube || inputs.grid->issubgrid())
        throw std::invalid_argument(std::string(name) + ": the inputs need the world's full cube grid");
    if (!(inputs.radius_m > 0.0) || !std::isfinite(inputs.sealevel))
        throw std::invalid_argument(std::string(name) + ": the radius must be positive and the sea level finite");
    const std::size_t size = inputs.grid->cellcount();
    const std::vector<std::size_t> sizes = { inputs.elevation_m.size(), inputs.network.sea.size(), inputs.runoff_mm_yr.size(),
                                             inputs.lakesurface_m.size(), inputs.januarytemperature_c.size(),
                                             inputs.julytemperature_c.size(), inputs.januaryprecipitation_mm.size(),
                                             inputs.julyprecipitation_mm.size(), inputs.icethickness_m.size(),
                                             inputs.climate.size() };
    for (const std::size_t value : sizes)
        if (value != size)
            throw std::invalid_argument(std::string(name) + ": a global field does not hold one value per cell of the world");
    if (kept0.empty())
        throw std::invalid_argument(std::string(name) + ": the kept box holds no level-0 cell");
    if (parameters.tributaries.enabled != 0)
        throw std::invalid_argument(std::string(name) + ": the tributaries are not ported to cube windows "
                                                        "(landscape.regional.tributary_enabled must be 0)");
}
}

CubeRegionResult runcuberegionalpass(const RegionalGlobalInputs& inputs, const std::vector<mesh::LevelKey>& kept0, int scale,
                                     int halo, const RegionalParameters& parameters)
{
    checkcube(inputs, kept0, parameters, "runcuberegionalpass");
    if (scale < 1 || parameters.meanpasses < 0 || halo < regionalhalo(parameters.meanpasses, false))
        throw std::invalid_argument("runcuberegionalpass: the scale must be at least 1 and the halo at least 2 x meanpasses + 3 "
                                    "global cells");
    const mesh::CellGrid& base = *inputs.grid;
    const int n = base.facecells();
    std::vector<mesh::LevelKey> kept = kept0;
    std::sort(kept.begin(), kept.end());
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    const std::vector<mesh::LevelKey> domain0 = dilatecells(base, kept, halo);
    const auto started = std::chrono::steady_clock::now();
    const std::shared_ptr<const mesh::CellGrid> grid = mesh::cubewindowgrid(n * scale, inputs.radius_m, descendants(n, scale, domain0), false);
    const std::vector<std::size_t> global = ancestors(inputs, *grid, scale);
    const std::vector<uint8_t> keptcells = keptflags(global, inputs, kept);
    const FineWarp warp = regionalwarp(*grid, n, inputs.seed, parameters);
    const GlobalIncision incision = globalincision(inputs, parameters);
    const RegionalDetails details{};
    CubeLevel context{ inputs, *grid, keptcells, parameters, warp, incision, details };
    context.k = scale;
    LevelState state;
    levelsteps(context, state, global);
    CubeRegionResult result;
    result.fields = finish(context, state, global);
    result.kept = mesh::cubewindowgrid(n * scale, inputs.radius_m, descendants(n, scale, kept));
    result.keptcells.push_back(result.kept->cellcount());
    result.domaincells.push_back(grid->cellcount());
    (void)started;
    return result;
}

CubeRegionResult runcuberefinement(const RegionalGlobalInputs& inputs, const std::vector<mesh::LevelKey>& kept0, int halo,
                                   int levels, const RegionalParameters& parameters, const RefinementParameters& refinement)
{
    checkcube(inputs, kept0, parameters, "runcuberefinement");
    const int factor = refinement.factor;
    if (levels < 1 || factor < 2)
        throw std::invalid_argument("runcuberefinement: the level must be at least 1 and the factor at least 2");
    if (refinement.smoothpasses < 0 || !std::isfinite(refinement.smoothrelaxation) || !(refinement.smoothrelaxation > 0.0))
        throw std::invalid_argument("runcuberefinement: refinement.smooth_passes must not be negative and "
                                    "refinement.smooth_relaxation must be positive and finite");
    if (refinement.prolongationlimit != 0 && refinement.prolongationlimit != 1)
        throw std::invalid_argument("runcuberefinement: refinement.prolongation_limit must be 0 or 1");
    if (refinement.smoothbasis != 0)
        throw std::invalid_argument("runcuberefinement: the cube has the Catmull-Rom smoothing basis (0) only");
    if (refinement.parentlift < 0 || refinement.parentlift > 2)
        throw std::invalid_argument("runcuberefinement: refinement.parent_lift must be 0, 1 or 2");
    if (!std::isfinite(refinement.channelcellslope) || refinement.channelcellslope < 0.0)
        throw std::invalid_argument("runcuberefinement: refinement.channel_cell_slope must be finite and not negative");
    if (!std::isfinite(refinement.detailcells) || !(refinement.detailcells > 0.0))
        throw std::invalid_argument("runcuberefinement: refinement.detail_cells must be positive and finite");
    if (parameters.meanpasses < 0 || halo < refinementhalo(parameters.meanpasses, refinement.smoothpasses, false))
        throw std::invalid_argument("runcuberefinement: the halo must be at least refinementhalo, 2 x meanpasses + 7 + 2 x "
                                    "smooth_passes global cells");
    std::int64_t finest = 1;
    for (int level = 0; level < levels; level++)
    {
        finest *= factor;
        if (finest > INT_MAX)
            throw std::invalid_argument("runcuberefinement: factor^level must fit an int");
    }
    const mesh::CellGrid& base = *inputs.grid;
    const int n = base.facecells();
    const double radius = inputs.radius_m;

    CubeRegionResult result;
    result.chain = keyhalochain(base, factor, levels, kept0, halo, prolongationreach, solvemargin(parameters.relief.catchmentreach));
    const std::vector<KeyLevelWindow>& chain = result.chain;
    const GlobalIncision incision = globalincision(inputs, parameters);
    RegionalDetails details;
    details.valleys = parameters.valleys.floor_km > 0.0;
    details.troughs = parameters.troughs.width_km > 0.0;
    details.cirques = parameters.cirques.radius_km > 0.0;
    details.factor = factor;
    details.detailcells = refinement.detailcells;
    details.globallakedepth = true;
    details.lakes.flood = true;

    // Level 0: the global fields on chain[0]'s cells (level 1's projection parent and its masks).
    std::shared_ptr<const mesh::CellGrid> parentgrid = mesh::cubewindowgrid(n, radius, chain[0].keptglobal, false);
    std::vector<float> parentelevation(parentgrid->cellcount()), parentlake(parentgrid->cellcount());
    std::vector<uint8_t> parentsea(parentgrid->cellcount());
    for (mesh::CellId cell = 0; cell < parentgrid->cellcount(); cell++)
    {
        const std::size_t global = base.localofkey(parentgrid->globalkey(cell));
        const bool sea = inputs.network.sea[global] != 0;
        parentelevation[cell] = inputs.elevation_m[global];
        parentsea[cell] = sea ? 1 : 0;
        parentlake[cell] = !sea && inputs.lakesurface_m[global] > 0.0f ? inputs.lakesurface_m[global] : 0.0f;
    }
    CubeParent parentfields;
    int sharepasses = -1, sharemisses = 0;
    LevelState state;
    std::shared_ptr<const mesh::CellGrid> grid;
    std::vector<uint8_t> kept;
    std::vector<std::size_t> global;
    FineWarp warp;
    std::vector<uint8_t> seagate, ownseaparent, lakereach;
    std::vector<float> lakeparent;
    for (int l = 1; l <= levels; l++)
    {
        const auto started = std::chrono::steady_clock::now();
        const KeyLevelWindow& window = chain[static_cast<std::size_t>(l)];
        const int k = static_cast<int>(levelscale(factor, l));
        const std::int64_t parentcells = static_cast<std::int64_t>(n) * levelscale(factor, l - 1);
        LevelDiagnostics diagnostics;
        diagnostics.level = l;
        diagnostics.cellsize_km = cubecellarc_m(radius, n) / 1000.0 / static_cast<double>(k);
        diagnostics.window.factor = factor;
        diagnostics.window.level = l;
        diagnostics.window.halo = window.halo;

        grid = keywindowgrid(window, radius, false);
        global = ancestors(inputs, *grid, k);
        kept = keptflags(global, inputs, window.keptglobal);
        warp = regionalwarp(*grid, n, inputs.seed, parameters);
        state = LevelState{};
        state.diagnostics = &diagnostics;
        if (l >= 2)
        {
            Prolongation prolonged = prolong(*parentfields.kept, parentfields.elevation_m, *grid, factor,
                                             refinement.prolongationlimit == 1, &kept);
            state.start.start_m = std::move(prolonged.values);
            diagnostics.prolongationclamped = prolonged.clamped;
            diagnostics.maxprolongationclamp_m = prolonged.maxclamp_m;
        }

        // The parent's masks per cell (S2.8): its parent and the parent's 8 neighbours at the parent level.
        const std::size_t count = grid->cellcount();
        seagate.assign(count, 0);
        ownseaparent.assign(count, 0);
        lakereach.assign(count, 0);
        lakeparent.assign(count, 0.0f);
        for (mesh::CellId cell = 0; cell < count; cell++)
        {
            const mesh::LevelKey parent = parentkey(*grid, grid->globalkey(cell), factor);
            const mesh::CellId own = parentgrid->localofkey(parent);
            if (own == mesh::nocell)
                throw std::logic_error("runcuberefinement: a parent lies outside the parent level's kept cells");
            bool sea = parentsea[own] != 0, wet = parentlake[own] > 0.0f;
            for (int slot = 1; slot <= mesh::slotcount; slot++)
            {
                const mesh::LevelKey other = mesh::cubekeyneighbour(parentcells, parent, slot);
                const mesh::CellId local = other == mesh::nokey ? mesh::nocell : parentgrid->localofkey(other);
                if (local == mesh::nocell)
                    continue;
                sea = sea || parentsea[local] != 0;
                wet = wet || parentlake[local] > 0.0f;
            }
            seagate[cell] = sea ? 1 : 0;
            lakereach[cell] = wet ? 1 : 0;
            ownseaparent[cell] = parentsea[own];
            lakeparent[cell] = parentlake[own];
        }
        details.parentseagate = l >= 2 ? &seagate : nullptr;
        details.coarselevel = l < levels;
        details.lakes.parentsurface_m = l >= 2 ? &lakeparent : nullptr;
        details.lakes.parentreach = l >= 2 ? &lakereach : nullptr;

        CubeLevel context{ inputs, *grid, kept, parameters, warp, incision, details };
        context.level = l;
        context.factor = factor;
        context.k = k;
        context.refinement = &refinement;
        context.projectionparent = parentgrid.get();
        context.projectionfield = &parentelevation;
        context.parent = l >= 2 ? &parentfields : nullptr;
        levelsteps(context, state, global);
        state.diagnostics = nullptr;
        if (l == 1)
        {
            sharepasses = state.sharepasses;
            sharemisses = state.sharemisses;
        }
        details.parentseagate = nullptr;
        details.lakes.parentsurface_m = nullptr;
        details.lakes.parentreach = nullptr;
        diagnostics.ringbasincells = state.ringbasincells;
        diagnostics.unreachedcells = state.unreachedcells;
        diagnostics.sinkjunctions = state.sinkjunctions;
        diagnostics.diagonallinks = state.diagonallinks;
        diagnostics.cornerlinks = state.cornerlinks;
        diagnostics.mergedpieces = state.mergedpieces;
        diagnostics.drytails = state.drytails;
        diagnostics.reliefconverged = state.reliefconverged;
        diagnostics.changedreceivers = state.changedreceivers;

        // The loop's checks after the projection (refinement_levels.cpp), on the window's table.
        const std::vector<float>& elevation = state.elevation_m;
        const std::vector<float>& lakesurface = state.water.lakesurface_m;
        const std::vector<float> lakelevel = channellakelevels(inputs, *grid, state.channels);
        const bool withbeds = state.bed_m.size() == count;
        float largest = 0.0f;
        bool belowfound = false;
        int aligned = 0, roughcells = 0;
        double roughness = 0.0;
        const auto dryland = [&](std::size_t cell) { return state.sea[cell] == 0 && !(lakesurface[cell] > 0.0f); };
        for (mesh::CellId cell = 0; cell < count; cell++)
        {
            if (kept[cell] == 0)
                continue;
            const bool wet = lakesurface[cell] > 0.0f;
            const bool sea = state.sea[cell] != 0, channel = state.channels.channel[cell] != 0;
            largest = std::max(largest, std::fabs(elevation[cell]));
            if (!wet && !sea && !channel)
            {
                bool below = false, below4 = false;
                for (int slot = 1; slot <= mesh::slotcount; slot++)
                {
                    const std::int32_t other = grid->neighbour(cell, slot);
                    if (other < 0)
                        continue;
                    const float surface = lakesurface[static_cast<std::size_t>(other)];
                    below = below || (surface > 0.0f && elevation[cell] <= surface);
                    if (slot % 2 == 1)
                        below4 = below4 || (surface > 0.0f && elevation[cell] < surface);
                }
                diagnostics.shoresbelowlake += below ? 1 : 0;
                diagnostics.drybelowlake += below4 ? 1 : 0;
            }
            if (!sea && channel)
            {
                if (wet)
                {
                    if (elevation[cell] > lakesurface[cell])
                    {
                        diagnostics.wetchannelsabove++;
                        diagnostics.wetchannelexcess_m = std::max(diagnostics.wetchannelexcess_m,
                                                                  static_cast<double>(elevation[cell]) - lakesurface[cell]);
                    }
                }
                else
                {
                    const float bed = withbeds ? state.bed_m[cell] : elevation[cell];
                    if (lakelevel[cell] > 0.0f && bed <= lakelevel[cell])
                    {
                        bool guarded = false;
                        for (int slot = 1; slot <= mesh::slotcount; slot++)
                        {
                            const std::int32_t other = grid->neighbour(cell, slot);
                            if (other < 0)
                                continue;
                            const float surface = lakesurface[static_cast<std::size_t>(other)];
                            guarded = guarded || (surface > 0.0f && surface != lakelevel[cell]);
                        }
                        if (guarded)
                            diagnostics.channelbetweenlakes++;
                        else if (elevation[cell] <= lakelevel[cell])
                            diagnostics.drychannelinlake++;
                    }
                }
            }
            if (!wet && !sea && elevation[cell] <= inputs.sealevel)
            {
                diagnostics.landbelowsea++;
                const double depth = static_cast<double>(inputs.sealevel) - static_cast<double>(elevation[cell]);
                if (!belowfound || depth > diagnostics.maxlandbelowsea_m)
                    diagnostics.maxlandbelowsea_m = depth;
                belowfound = true;
            }
            if (!state.water.rimraise_m.empty() && state.water.rimraise_m[cell] > 0.0f)
            {
                diagnostics.rimraised++;
                diagnostics.maxrimraise_m = std::max(diagnostics.maxrimraise_m, static_cast<double>(state.water.rimraise_m[cell]));
            }
            if (sea && ownseaparent[cell] == 0)
            {
                diagnostics.seainlandparents++;
                diagnostics.seabeyondparent += seagate[cell] == 0 ? 1 : 0;
            }
            if (!sea)
            {
                diagnostics.landcells++;
                if (channel)
                {
                    diagnostics.channelcells++;
                    if (!wet)
                        diagnostics.maxchannelrise_m = std::max(diagnostics.maxchannelrise_m,
                                                                static_cast<double>(elevation[cell]) - state.bed_m[cell]);
                }
            }
            // The lake edges to the kept neighbours in slots 3 and 5 (east and south in the cell's frame), aligned
            // when the two cells lie in different level-0 cells; the roughness over the four cardinal neighbours.
            for (const int slot : { 3, 5 })
            {
                const std::int32_t other = grid->neighbour(cell, slot);
                if (other < 0 || kept[static_cast<std::size_t>(other)] == 0)
                    continue;
                if (wet != (lakesurface[static_cast<std::size_t>(other)] > 0.0f))
                {
                    diagnostics.lakeedges++;
                    aligned += global[static_cast<std::size_t>(other)] != global[cell] ? 1 : 0;
                }
            }
            if (!dryland(cell))
                continue;
            double sum = 0.0;
            bool all = true;
            for (const int slot : { 1, 3, 5, 7 })
            {
                const std::int32_t other = grid->neighbour(cell, slot);
                if (other < 0 || !dryland(static_cast<std::size_t>(other)))
                {
                    all = false;
                    break;
                }
                sum += static_cast<double>(elevation[static_cast<std::size_t>(other)]);
            }
            if (!all)
                continue;
            const double laplacian = sum - 4.0 * static_cast<double>(elevation[cell]);
            roughness += laplacian * laplacian;
            roughcells++;
        }
        diagnostics.meanbound_m = largest > 0.0f ? std::ldexp(1.0, std::ilogb(largest) - 24) : 0.0;
        diagnostics.drainagedensity_per_km =
            diagnostics.landcells > 0
                ? static_cast<double>(diagnostics.channelcells) / (static_cast<double>(diagnostics.landcells) * diagnostics.cellsize_km)
                : 0.0;
        diagnostics.lakeedgealigned =
            diagnostics.lakeedges > 0 ? static_cast<double>(aligned) / static_cast<double>(diagnostics.lakeedges) : 0.0;
        diagnostics.landroughness_m = roughcells > 0 ? std::sqrt(roughness / static_cast<double>(roughcells)) : 0.0;
        diagnostics.tributaries = state.tributaries;
        diagnostics.tributarycells = state.tributarycells;
        diagnostics.tributarybedsabove = state.tributarybedsabove;

        // The kept cells as the next level's parent (and the global-mean error over their level-0 blocks).
        const std::shared_ptr<const mesh::CellGrid> keptgrid =
            mesh::cubewindowgrid(static_cast<int>(window.levelcells()), radius, window.keptkeys(), l == levels);
        std::vector<float> keptelevation(keptgrid->cellcount());
        std::vector<std::size_t> keptlocal(keptgrid->cellcount());
        for (mesh::CellId cell = 0; cell < keptgrid->cellcount(); cell++)
        {
            keptlocal[cell] = grid->localofkey(keptgrid->globalkey(cell));
            keptelevation[cell] = elevation[keptlocal[cell]];
        }
        {
            const WindowBlocks blocks = windowblocks(*keptgrid, k);
            double error = 0.0;
            for (std::size_t b = 0; b < blocks.count(); b++)
                error = std::max(error, std::fabs(blockmean(blocks, b, keptelevation) -
                                                  static_cast<double>(inputs.elevation_m[base.localofkey(blocks.parents[b])])));
            diagnostics.globalmeanerror_m = error;
        }
        result.keptcells.push_back(keptgrid->cellcount());
        result.domaincells.push_back(count);
        if (l < levels)
        {
            parentfields.kept = keptgrid;
            parentfields.elevation_m = std::move(keptelevation);
            parentfields.lift_m.resize(keptlocal.size());
            parentfields.solved_m.resize(keptlocal.size());
            parentfields.sea.resize(keptlocal.size());
            parentfields.lakesurface_m.resize(keptlocal.size());
            for (std::size_t cell = 0; cell < keptlocal.size(); cell++)
            {
                parentfields.lift_m[cell] = state.lift_m[keptlocal[cell]];
                parentfields.solved_m[cell] = state.solved_m[keptlocal[cell]];
                parentfields.sea[cell] = state.sea[keptlocal[cell]];
                parentfields.lakesurface_m[cell] = lakesurface[keptlocal[cell]];
            }
            parentgrid = keptgrid;
            parentelevation = parentfields.elevation_m;
            parentsea = parentfields.sea;
            parentlake = parentfields.lakesurface_m;
        }
        else
        {
            diagnostics.steps.push_back(LevelStep::ice);
            result.kept = keptgrid;
        }
        diagnostics.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        result.levels.push_back(std::move(diagnostics));
    }

    // The ice above the ELA, the rivers and the climate at the final level.
    const auto started = std::chrono::steady_clock::now();
    CubeLevel context{ inputs, *grid, kept, parameters, warp, incision, details };
    context.level = levels;
    context.factor = factor;
    context.k = static_cast<int>(finest);
    context.refinement = &refinement;
    result.fields = finish(context, state, global);
    result.levels.back().seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    result.fields.sharepasses = sharepasses;
    result.fields.sharemisses = sharemisses;
    return result;
}
}
