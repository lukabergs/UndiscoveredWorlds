#include "simulations/landscape/regional_solve.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "simulations/landscape/analytical_relief.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/upsampling.hpp"

// Multiscale plan S3.0: regional_pass.cpp's step 4 moved here unchanged (S3.2 changes it).

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
    // R9: the solve's uplift term from the incision's rate (analytical_relief.hpp), only
    // with the incision active; empty otherwise (the solve without uplift exactly).
    if (incision.active)
    {
        solveinput.uplift_m_myr = upsample(incision.rate_m_myr, width, height, window);
        for (float& rate : solveinput.uplift_m_myr)
            rate = std::max(0.0f, rate);
    }
    AnalyticalReliefResult solved = solveanalyticalrelief(metric, solveinput, context.parameters.relief);
    state.reliefconverged = solved.converged;
    state.changedreceivers = solved.changedreceivers;
    std::vector<float>& elevation = state.elevation_m;
    elevation = surface;
    for (std::size_t cell = 0; cell < count; cell++)
        if (solveinput.pinned[cell] == 0)
            elevation[cell] = solved.elevation_m[cell];
    {
        // The lowering by the solve over the kept window's unpinned (so land) cells.
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double total = 0.0, largest = 0.0;
        std::size_t cells = 0;
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
            }
        state.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        state.maxlowering_m = largest;
    }
    state.solveddischarge_m3_yr = std::move(solved.discharge_m3_yr);
    state.pinned = std::move(solveinput.pinned);
    state.solved_m = elevation;
    recordlevelstep(state, LevelStep::solve);
}
}
