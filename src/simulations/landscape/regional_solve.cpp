#include "simulations/landscape/regional_solve.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/upsampling.hpp"

// Multiscale plan S3.0: regional_pass.cpp's step 4 moved here unchanged; S3.2 adds the solve against
// its parent and the channel cells' means (regional_solve.hpp).

namespace landscape
{
namespace
{
std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}
}

void solvelevel(const LevelContext& context, LevelState& state)
{
    const RegionalGlobalInputs& inputs = context.inputs;
    const RegionalWindow& window = context.window;
    const GlobalIncision& incision = context.incision;
    const int width = inputs.width, height = inputs.height;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const std::int64_t planetfinewidth = static_cast<std::int64_t>(width) * window.scale;
    const std::vector<uint8_t>& sea = state.sea;
    const FineChannels& channels = state.channels;
    const std::vector<float>& surface = state.surface_m;

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
    // S3.9: bounded catchments in the refinement loop only (the single pass keeps the unbounded solve);
    // each cell's global cell from its absolute fine position.
    AnalyticalReliefParameters relief = context.parameters.relief;
    if (context.refinement == nullptr)
        relief.catchmentreach = 0;
    if (relief.catchmentreach > 0)
    {
        solveinput.globalcolumn.resize(count);
        solveinput.globalrow.resize(count);
        solveinput.globalwidth = width;
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                solveinput.globalcolumn[cell] = static_cast<std::int32_t>(window.absolutefinecolumn(x) / window.scale);
                solveinput.globalrow[cell] = static_cast<std::int32_t>(window.absolutefinerow(y) / window.scale);
            }
    }
    // R9: the solve's uplift term from the incision's rate (analytical_relief.hpp), only
    // with the incision active; empty otherwise (the solve without uplift exactly).
    if (incision.active)
    {
        solveinput.uplift_m_myr = upsample(incision.rate_m_myr, width, height, window);
        for (float& rate : solveinput.uplift_m_myr)
            rate = std::max(0.0f, rate);
    }

    // S3.2 (refinement loop, level >= 2): the solve against its parent. The parent's lift (1: its
    // lift_m, the final elevation less the elevation entering its smoothing; 2: its final elevation less
    // its solved elevation) prolonged as the start was (refinement.prolongation_limit): the solve starts
    // from the start surface less it on the unpinned cells (the parent's solved frame); the lift is added back to the
    // unpinned cells after it, F* = Z_l + lift (pinned cells keep their surface), so the level's solve
    // from its own beds keeps what the parent's later steps put on its solve.
    const RefinementParameters* refinement = context.refinement;
    std::vector<float> added;
    if (refinement != nullptr && refinement->parentlift != 0 && context.parent != nullptr && context.levelwindow != nullptr)
    {
        const ParentFields& parent = *context.parent;
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
        added = prolong(parent.window, parentlift, *context.levelwindow, refinement->prolongationlimit == 1).values;
        if (added.size() != count)
            throw std::invalid_argument("solvelevel: the parent's lift must prolong onto the window's domain");
        // Only the unpinned cells are de-lifted: de-lifting the pinned beds too (all of them, or those
        // below their link's filled elevation) raised level-2 lambda with the forms off from 1,172 to
        // 2,034 m and level-3 land below the sea level from 62 to 135 m on erosion run 9 (S3.2's review).
        for (std::size_t cell = 0; cell < count; cell++)
        {
            if (solveinput.pinned[cell] == 0)
                solveinput.start_m[cell] =
                    static_cast<float>(static_cast<double>(solveinput.start_m[cell]) - static_cast<double>(added[cell]));
            if (solveinput.pinned[cell] != 0)
                added[cell] = 0.0f;
        }
    }

    AnalyticalReliefResult solved = solveanalyticalrelief(metric, solveinput, relief);
    state.reliefconverged = solved.converged;
    state.changedreceivers = solved.changedreceivers;
    // S3.9's unreached cells stay unpinned here: the solve returned their start (with the parent's
    // lift, the de-lifted start), so their elevation and solved_m follow the unpinned cells' rules
    // below and the parent's lift of the next level is the one they really carry.
    // S3.9's diagnostics over the kept box (refinement loop only): the unreached cells, and the cells
    // whose outlet on the last pass's receivers is a cell of the domain's outer ring on a window edge
    // (the west and east columns; the north and south rows unless they lie on a pole, where every
    // window clips alike).
    state.unreachedcells = 0;
    state.ringbasincells = 0;
    if (context.levelwindow != nullptr)
    {
        const LevelWindow& level = *context.levelwindow;
        const bool north = level.domainrow0 > 0;
        const bool south = level.domainrow0 + domainheight < level.levelheight();
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (!solved.unreached.empty() && solved.unreached[cell] != 0)
                    state.unreachedcells++;
                const int32_t outlet = solved.outlet[cell];
                if (outlet < 0)
                    continue;
                const int ox = outlet % domainwidth, oy = outlet / domainwidth;
                if (ox == 0 || ox == domainwidth - 1 || (north && oy == 0) || (south && oy == domainheight - 1))
                    state.ringbasincells++;
            }
    }
    std::vector<float>& elevation = state.elevation_m;
    elevation = surface;
    // A negative lift takes no land cell below the sea level (the projection's floor for land) unless
    // the solve put it there or its parent cell lies there.
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (solveinput.pinned[cell] != 0)
                continue;
            const float value = solved.elevation_m[cell];
            if (added.empty())
            {
                elevation[cell] = value;
                continue;
            }
            const LevelWindow& parentwindow = context.parent->window;
            const LevelWindow& levelwindow = *context.levelwindow;
            const float parentvalue = context.parent->elevation_m[parentwindow.index(
                parentwindow.localcolumn(levelwindow.parentcolumn(x)),
                static_cast<int>(parentwindow.localrow(levelwindow.parentrow(y))))];
            elevation[cell] = std::max(static_cast<float>(static_cast<double>(value) + static_cast<double>(added[cell])),
                                       std::min({value, abovesea, parentvalue}));
        }
    {
        // The lowering by the solve over the kept window's unpinned (so land) cells (with the parent's
        // lift: the start surface less the result with the lift added back).
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double total = 0.0, largest = 0.0;
        std::size_t cells = 0;
        std::vector<double> lowerings;
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
                if (state.diagnostics != nullptr)
                    lowerings.push_back(lowering);
            }
        state.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        state.maxlowering_m = largest;
        // S3.9: the loop's record of the lowering's 50th and 99th percentiles (nearest rank).
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
        // A land channel cell's elevation becomes its cell mean m = b + min(h, theta c_l): b its bed
        // (kept in state.bed_m), h the mean height above b of its non-channel land 8-neighbours in the
        // domain (at least 0; 0 without such neighbours), c_l the level's north-south cell size. Read on
        // the field before any channel cell is written, so the order of the cells does not matter.
        const double cap = refinement->channelcellslope * metric.dy_m;
        const std::vector<float> field = elevation;
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (channels.channel[cell] == 0 || sea[cell] != 0)
                    continue;
                const double bed = static_cast<double>(state.bed_m[cell]);
                double sum = 0.0;
                int neighbours = 0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        const int nx = x + dx, ny = y + dy;
                        if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= domainwidth || ny >= domainheight)
                            continue;
                        const std::size_t other = at(nx, ny, domainwidth);
                        if (channels.channel[other] != 0 || sea[other] != 0)
                            continue;
                        sum += static_cast<double>(field[other]) - bed;
                        neighbours++;
                    }
                const double height_m = neighbours > 0 ? std::max(0.0, sum / neighbours) : 0.0;
                elevation[cell] = static_cast<float>(bed + std::min(height_m, cap));
            }
    }
    // The solved elevation: without the parent's lift (the next level's ParentFields::solved_m).
    state.solved_m = elevation;
    if (!added.empty())
        for (std::size_t cell = 0; cell < count; cell++)
            if (state.pinned[cell] == 0)
                state.solved_m[cell] = solved.elevation_m[cell];
    recordlevelstep(state, LevelStep::solve);
}
}
