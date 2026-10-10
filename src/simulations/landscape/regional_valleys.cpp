#include "simulations/landscape/regional_valleys.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan Phase R2: the R6 contract, which step R9 implements (regional_valleys.hpp):
// the incision of the pinned channels by the world's uplift or denudation rate over the
// elapsed time, and the valley floors after the legacy expandrivers and pasterivercircle.

namespace landscape
{
namespace
{
constexpr std::size_t nocell = std::numeric_limits<std::size_t>::max();

// The receiver of a global land cell on the drainage tree (channel_network.cpp's
// globalreceiver): nocell for a sea cell, a sink and a direction beyond a pole.
std::size_t receiverof(const GlobalNetwork& network, std::size_t cell)
{
    const int direction = network.direction[cell];
    if (network.sea[cell] != 0 || direction < 1 || direction > directioncount)
        return nocell;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(network.width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(network.width)) + directionoffsety(direction);
    if (j < 0 || j >= network.height)
        return nocell;
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(network.width) +
           static_cast<std::size_t>((i + directionoffsetx(direction) + network.width) % network.width);
}

// The local global relief of one cell: the largest absolute difference between its
// elevation and its 8 neighbours' (columns wrap, rows beyond the poles do not exist), the
// expression of regional_surface.cpp's localrelief.
float localrelief(const std::vector<float>& elevation, int width, int height, std::size_t cell)
{
    const int i = static_cast<int>(cell % static_cast<std::size_t>(width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(width));
    const float own = elevation[cell];
    float largest = 0.0f;
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
        {
            const int column = ((i + di) % width + width) % width;
            const float other = elevation[static_cast<std::size_t>(row) * static_cast<std::size_t>(width) +
                                          static_cast<std::size_t>(column)];
            largest = std::max(largest, std::fabs(other - own));
        }
    }
    return largest;
}
}

GlobalIncision globalincision(const RegionalGlobalInputs& inputs, const RegionalParameters& parameters)
{
    GlobalIncision incision;
    const ValleyParameters& valleys = parameters.valleys;
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    const bool uplift = valleys.upliftscale > 0.0 && inputs.upliftrate_m_myr.size() == size;
    const bool denudation = valleys.denudationratescale > 0.0 && inputs.denudation_m.size() == size &&
                            inputs.erosionduration_myr > 0.0;
    const double elapsed = parameters.relief.elapsed_myr;
    const GlobalNetwork& network = inputs.network;
    if ((!uplift && !denudation) || !(elapsed > 0.0) || size == 0 || network.sea.size() != size ||
        network.direction.size() != size || network.filledelevation_m.size() != size)
        return incision;

    // The lowest incised level: the next float above the sea level, so a land channel at the
    // mouth never sits at the sea level itself (the bed law's lo, channel_network.hpp).
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    const bool withrelief = inputs.elevation_m.size() == size;

    // The rate and the capped depth I per land cell: at most F minus the lowest level and at
    // most the cell's local global relief, where the relief saturates once the plateau
    // itself lies within the knickpoint wave's reach (regional_valleys.hpp).
    std::vector<float> rate(size, 0.0f);
    std::vector<double> depth(size, 0.0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        if (network.sea[cell] != 0)
            continue;
        double value = 0.0;
        if (uplift)
            value = std::max(value, valleys.upliftscale * std::max(0.0, static_cast<double>(inputs.upliftrate_m_myr[cell])));
        if (denudation)
            value = std::max(value, valleys.denudationratescale * static_cast<double>(inputs.denudation_m[cell]) /
                                        inputs.erosionduration_myr);
        if (!std::isfinite(value))
            value = 0.0;
        rate[cell] = static_cast<float>(value);
        const double filled = static_cast<double>(network.filledelevation_m[cell]);
        double cap = filled - static_cast<double>(floorlevel);
        if (withrelief)
            cap = std::min(cap, static_cast<double>(localrelief(inputs.elevation_m, inputs.width, inputs.height, cell)));
        depth[cell] = std::max(0.0, std::min(value * elapsed, cap));
    }

    // Kahn's order over the forest of land links: a cell after all its donors.
    std::vector<std::size_t> receiver(size, nocell);
    std::vector<int> donors(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        receiver[cell] = receiverof(network, cell);
        if (receiver[cell] != nocell)
            donors[receiver[cell]]++;
    }
    // F'(r) = max(s+, min(F(r) - I(r), min over donors F'(u))) with s+ the next float above
    // the sea level, or F(r) where F(r) < s+; `lowest` holds the donors' minimum F' as they
    // are finished.
    std::vector<float> incised(network.filledelevation_m);
    std::vector<float> lowest(size, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> done(size, 0);
    std::vector<std::size_t> ready;
    for (std::size_t cell = 0; cell < size; cell++)
        if (donors[cell] == 0)
            ready.push_back(cell);
    const auto finish = [&](std::size_t cell, bool withdonors) {
        done[cell] = 1;
        if (network.sea[cell] != 0)
            return;
        const float filled = network.filledelevation_m[cell];
        if (!(filled >= floorlevel))
            return;  // at or below the sea level (or not a number): no incision
        float value = static_cast<float>(static_cast<double>(filled) - depth[cell]);
        if (withdonors)
            value = std::min(value, lowest[cell]);
        incised[cell] = std::min(filled, std::max(floorlevel, value));
    };
    while (!ready.empty())
    {
        const std::size_t cell = ready.back();
        ready.pop_back();
        finish(cell, true);
        const std::size_t next = receiver[cell];
        if (next == nocell)
            continue;
        lowest[next] = std::min(lowest[next], incised[cell]);
        if (--donors[next] == 0)
            ready.push_back(next);
    }
    // Cells on a receiver loop (bad input) and downstream of one: their own capped I.
    for (std::size_t cell = 0; cell < size; cell++)
        if (done[cell] == 0)
            finish(cell, false);

    bool any = false;
    std::vector<float> depthout(size, 0.0f);
    for (std::size_t cell = 0; cell < size; cell++)
    {
        const double value = static_cast<double>(network.filledelevation_m[cell]) - static_cast<double>(incised[cell]);
        depthout[cell] = static_cast<float>(std::max(0.0, value));
        any = any || incised[cell] < network.filledelevation_m[cell];
    }
    if (!any)
        return incision;
    incision.active = true;
    incision.incision_m = std::move(depthout);
    incision.filledelevation_m = std::move(incised);
    incision.rate_m_myr = std::move(rate);
    return incision;
}

void carvevalleys(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                  const FineChannels& channels, const std::vector<uint8_t>& sea, std::vector<float>& elevation_m,
                  double bandlower, double bandupper)
{
    const ValleyParameters& valleys = parameters.valleys;
    if (!(valleys.floor_km > 0.0))
        return;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || sea.size() != count || elevation_m.size() != count)
        return;
    const int k = window.scale;
    const double reference = valleys.reference_m3s > 0.0 ? valleys.reference_m3s : 1.0;
    const PatchMetric metric = finepatchmetric(window, inputs.radius_m);
    const double cell_m = metric.dy_m;
    // The refinement loop (bandupper > 0, multiscale plan S3.5): floors by Euclidean distance and,
    // with valleys.reach_km > 0, the sides' reach in km on the ground over the level's cell size,
    // its full width tested against the band; the single pass keeps R9's Chebyshev rule in cells.
    const bool loop = bandupper > 0.0;
    const bool reachkm = loop && valleys.reach_km > 0.0;
    const bool sides = valleys.sideslope > 0.0 && (reachkm || valleys.reachcells > 0);
    const double reach = !sides ? 0.0
                         : reachkm ? std::min(static_cast<double>(k), valleys.reach_km * 1000.0 / cell_m)
                                   : static_cast<double>(std::min(k, valleys.reachcells));
    // The floor in fine cells: kilometres on the ground over the fine north-south spacing.
    const double floorcells = valleys.floor_km * 1000.0 / cell_m;
    const double floormax = std::min(static_cast<double>(k), std::max(0.0, valleys.floormax_km * 1000.0 / cell_m));
    const double rise = valleys.sideslope * cell_m;  // per fine cell beyond the floor
    const float floorlevel = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    // S3.2: the beds are the channels' (FineChannels::bed_m, LevelState::bed_m), not the channel
    // cells' elevations, which the refinement loop sets to their cell means (regional_solve.hpp).
    const bool withbeds = channels.bed_m.size() == count;

    // The lowest target per cell over every channel that reaches it.
    std::vector<double> target(count, std::numeric_limits<double>::infinity());
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
            if (channels.channel[cell] == 0 || sea[cell] != 0)
                continue;
            const double discharge = std::max(0.0, static_cast<double>(channels.discharge_m3s[cell]));
            const double halfwidth = std::min(floormax, floorcells * std::sqrt(discharge / reference));
            // The band (multiscale plan S2.3): the full width (the floor or the sides' reach,
            // whichever is wider) decides the lower end, so a sub-cell floor with sides is
            // carved; the floor half-width the upper end (under one parent cell). S3.5: with the
            // reach in km the full width decides both ends (the whole form at one level).
            const double fullwidth = std::max(halfwidth, reach);
            if (loop && (fullwidth < bandlower || (reachkm ? fullwidth : halfwidth) >= bandupper))
                continue;  // carved at another level of the refinement
            const int extent = std::min(k, static_cast<int>(std::floor(fullwidth)));
            if (extent < 1)
                continue;
            const double bed = static_cast<double>(withbeds ? channels.bed_m[cell] : elevation_m[cell]);
            // Chebyshev in cells in the single pass (R9); in the loop Euclidean on the ground (S3.5), the
            // east-west offset scaled by the axis row's east-west over north-south cell size.
            const double ratio = loop ? eastwestratio(metric, y) : 1.0;
            const int columns = loop ? columnextent(fullwidth, ratio, k) : extent;
            for (int dy = -extent; dy <= extent; dy++)
            {
                const int ny = y + dy;
                if (ny < 0 || ny >= height)
                    continue;
                for (int dx = -columns; dx <= columns; dx++)
                {
                    const int nx = x + dx;
                    if (nx < 0 || nx >= width)
                        continue;
                    const double east = static_cast<double>(dx) * ratio;
                    const double distance = loop ? std::sqrt(east * east + static_cast<double>(dy * dy))
                                                 : static_cast<double>(std::max(std::abs(dx), std::abs(dy)));
                    double value;
                    if (distance <= halfwidth)
                        value = bed;
                    else if (distance <= reach)
                        value = bed + (distance - halfwidth) * rise;
                    else
                        continue;
                    const std::size_t next = static_cast<std::size_t>(ny) * static_cast<std::size_t>(width) + static_cast<std::size_t>(nx);
                    target[next] = std::min(target[next], value);
                }
            }
        }

    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (sea[cell] != 0 || !std::isfinite(target[cell]))
            continue;
        // S3.5: a channel cell keeps its value (its bed in the single pass; its cell mean in the
        // loop, which channelcellmeans lowers with the floor around it, never back to the bed).
        if (channels.channel[cell] != 0)
            continue;
        const float value = std::max(floorlevel, static_cast<float>(target[cell]));
        elevation_m[cell] = std::min(elevation_m[cell], value);
    }
}

double eastwestratio(const PatchMetric& metric, int y)
{
    if (y < 0 || static_cast<std::size_t>(y) >= metric.rowdx_m.size() || !(metric.dy_m > 0.0))
        return 1.0;
    return std::max(0.0, metric.rowdx_m[static_cast<std::size_t>(y)] / metric.dy_m);
}

int columnextent(double radius, double ratio, int k)
{
    if (!(radius > 0.0))
        return 0;
    if (!(ratio > 0.0))
        return k;
    return static_cast<int>(std::min(static_cast<double>(k), std::floor(radius / ratio)));
}

void channelcellmeans(const RegionalWindow& window, const FineChannels& channels, const std::vector<uint8_t>& sea,
                      double cap_m, const std::vector<float>& before, std::vector<float>& elevation_m)
{
    if (!(cap_m > 0.0))
        return;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (count == 0 || channels.channel.size() != count || channels.bed_m.size() != count || sea.size() != count ||
        before.size() != count || elevation_m.size() != count)
        return;
    const auto index = [width](int x, int y) {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
    };
    // Only non-channel cells are read and only channel cells written, so the order does not matter.
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = index(x, y);
            if (channels.channel[cell] == 0 || sea[cell] != 0)
                continue;
            const double bed = static_cast<double>(channels.bed_m[cell]);
            double sum = 0.0, carvedsum = 0.0;
            int neighbours = 0, carved = 0;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++)
                {
                    const int nx = x + dx, ny = y + dy;
                    if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= width || ny >= height)
                        continue;
                    const std::size_t other = index(nx, ny);
                    if (channels.channel[other] != 0 || sea[other] != 0)
                        continue;
                    const double above = static_cast<double>(elevation_m[other]) - bed;
                    sum += above;
                    neighbours++;
                    if (elevation_m[other] != before[other])
                    {
                        carvedsum += above;
                        carved++;
                    }
                }
            if (carved == 0)
                continue;
            // regional_solve.cpp's cell mean on the carved field, and never above the carved neighbours'
            // mean: under a form narrower than the diagonal (half-width below sqrt 2 cells) the uncarved
            // diagonal neighbours would otherwise hold the axis above its flanks.
            const double height_m = std::max(0.0, std::min(sum / neighbours, carvedsum / carved));
            const float mean = static_cast<float>(bed + std::min(height_m, cap_m));
            elevation_m[cell] = std::min(elevation_m[cell], mean);
        }
}

void incisionlevel(const LevelContext& context, LevelState& state)
{
    // R9 (regional_pass.cpp's step 3, moved here unchanged in S3.0): with the incision active
    // the bed law reads the network's incised filled elevations and the start surface lowered
    // by the upsampled incision; inactive, the channels read the inputs as they are.
    const GlobalIncision& incision = context.incision;
    if (!incision.active)
        return;
    const RegionalGlobalInputs& inputs = context.inputs;
    const std::vector<float>& surface = state.surface_m;
    const std::size_t count = surface.size();
    state.incisednetwork = inputs.network;
    state.incisednetwork.filledelevation_m = incision.filledelevation_m;
    const std::vector<float> fineincision = upsample(incision.incision_m, inputs.width, inputs.height, context.window);
    std::vector<float>& bedsurface = state.bedsurface_m;
    bedsurface = surface;
    for (std::size_t cell = 0; cell < count; cell++)
        bedsurface[cell] = static_cast<float>(static_cast<double>(surface[cell]) - std::max(0.0f, fineincision[cell]));
    recordlevelstep(state, LevelStep::incision);
}

void valleyfloorslevel(const LevelContext& context, LevelState& state)
{
    // 4b. The valley floors around the channels (R9, a no-op while valley_floor_km is 0),
    // after the solve's diagnostics and before the lakes, which read the result
    // (regional_pass.cpp's step 4b, moved here unchanged in S3.0).
    const RegionalDetails& details = context.details;
    if (!details.valleys)
        return;
    double bandlower = 0.0, bandupper = 0.0;
    formband(context, bandlower, bandupper);
    const double cap = channelcellcap_m(context);
    std::vector<float> before;
    if (state.diagnostics != nullptr || cap > 0.0)
        before = state.elevation_m;
    carvevalleys(context.inputs, context.window, context.parameters, state.channels, state.sea, state.elevation_m,
                 bandlower, bandupper);
    channelcellmeans(context.window, state.channels, state.sea, cap, before, state.elevation_m);
    if (state.diagnostics != nullptr)
        state.diagnostics->valleycells = keptcellschanged(context.window, before, state.elevation_m);
    recordlevelstep(state, LevelStep::valleyfloors);
}

void formband(const LevelContext& context, double& bandlower, double& bandupper)
{
    const RegionalDetails& details = context.details;
    bandlower = 0.0;
    bandupper = 0.0;
    if (details.factor < 2)
        return;
    bandlower = details.detailcells / 2.0;
    // S3.5: level 1 also carves the forms above its band (their extent clamped to one global cell
    // by each unit), since no coarser level of the loop carves them.
    bandupper = context.level == 1 ? std::numeric_limits<double>::infinity() : details.detailcells * details.factor / 2.0;
}

double channelcellcap_m(const LevelContext& context)
{
    if (context.refinement == nullptr || !(context.refinement->channelcellslope > 0.0))
        return 0.0;
    return context.refinement->channelcellslope * finepatchmetric(context.window, context.inputs.radius_m).dy_m;
}
}
