#include "simulations/landscape/regional_water.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/sphere_metric.hpp"

// Erosion plan steps E3, R0, R1 and R4 (docs/architecture/erosion-hydrology-plan.md,
// Phase E and Phase R): moved from regional_pass.cpp; R1 wets the channel cells in lakes;
// R4 draws the lake shores and ice edges as the zero contour of a signed indicator.

namespace landscape
{
namespace
{
// Salts of the indicators' hashed magnitudes and of their shore noises (R4).
constexpr std::uint64_t lakesalt = 0x1a4e;
constexpr std::uint64_t lakenoisesalt = 0x1a4f;
constexpr std::uint64_t icesalt = 0x1ce0;
constexpr std::uint64_t icenoisesalt = 0x1ce1;
// The indicator's magnitudes: 101 .. 900 thousandths of its scale (1), the legacy
// generator's range (plan section 6), so a zero contour lies between 0.1 and 0.9 of
// the way between two cell centres of different kinds. The form of the indicator, not
// a tuning value: the shore noise (shorenoise) is a fraction of the scale.
constexpr int magnitudeleast = 101;
constexpr int magnitudemost = 900;
constexpr double magnitudescale = 1000.0;
// The thickness a fine cell the ice indicator ices carries at least where the upsampled
// thickness is smaller (bicubic undershoot or the global ice edge): a floor that keeps
// the indicator's ice visible as ice (thickness above zero), not a tuning value.
constexpr float iceminimum_m = 1.0f;
constexpr std::size_t none = static_cast<std::size_t>(-1);

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

// The receiver of global cell `cell` on the global network (channel_network.cpp's
// globalreceiver): none for a sea cell, a sink and a direction beyond a pole.
std::size_t globalreceiver(const GlobalNetwork& network, std::size_t cell)
{
    if (cell >= network.direction.size() || cell >= network.sea.size())
        return none;
    const int direction = network.direction[cell];
    if (network.sea[cell] != 0 || direction < 1 || direction > directioncount)
        return none;
    const int i = static_cast<int>(cell % static_cast<std::size_t>(network.width));
    const int j = static_cast<int>(cell / static_cast<std::size_t>(network.width)) + directionoffsety(direction);
    if (j < 0 || j >= network.height)
        return none;
    return at(wrapcolumn(i + directionoffsetx(direction), network.width), j, network.width);
}

// True when global cell (i, j) or one of its 8 neighbours satisfies `test` (columns
// wrap, rows beyond the poles do not exist).
template<typename Test>
bool aroundcell(int width, int height, int i, int j, Test test)
{
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
            if (test(at(wrapcolumn(i + di, width), row, width)))
                return true;
    }
    return false;
}

// The signed indicator: -m on the cells `negative` marks, +m elsewhere, with m a hash of
// the cell's global indices, the seed and the salt in 0.101 .. 0.900.
std::vector<float> signedindicator(const std::vector<uint8_t>& negative, int width, int height, std::uint64_t seed,
                                   std::uint64_t salt)
{
    const std::uint64_t base = deterministiccombine(seed, salt);
    std::vector<float> indicator(negative.size());
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            const double magnitude = deterministicrandom(base, magnitudeleast, magnitudemost, i, j) / magnitudescale;
            indicator[cell] = static_cast<float>(negative[cell] != 0 ? -magnitude : magnitude);
        }
    return indicator;
}

// Phase E's rules with R1's channel rule (shores.shorenoise 0).
FineWater phaseewater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const std::vector<uint8_t>& sea,
                      const FineChannels& channels, std::vector<float>& elevation_m)
{
    const int width = inputs.width, height = inputs.height;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);

    std::vector<uint8_t> nearice(globalsize, 0);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
            nearice[at(i, j, width)] = aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.icethickness_m[other] > 0.0f;
            }) ? 1 : 0;
    const std::vector<float> ice = upsample(inputs.icethickness_m, width, height, window);

    FineWater water;
    water.lakesurface_m.assign(count, 0.0f);
    water.icethickness_m.assign(count, 0.0f);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            // R1: a channel cell carries the lake's water where its bed is at or below
            // the higher of the lake surface and the global cell's filled elevation,
            // and its elevation is lowered to at most the surface. Beds are pinned at
            // or below the filled elevation, which equals the surface of a lake that
            // spills and lies above that of a closed lake under its sill; the strict
            // rule left rivers dry through their lakes. A shore cell whose bed follows
            // a start surface above that level stays dry.
            const float surfacelevel = inputs.lakesurface_m[global];
            if (sea[cell] == 0 && surfacelevel > 0.0f)
            {
                if (channels.channel[cell] != 0)
                {
                    const float level = std::max(surfacelevel, inputs.network.filledelevation_m[global]);
                    if (elevation_m[cell] <= level)
                    {
                        water.lakesurface_m[cell] = surfacelevel;
                        elevation_m[cell] = std::min(elevation_m[cell], surfacelevel);
                    }
                }
                else if (elevation_m[cell] < surfacelevel)
                    water.lakesurface_m[cell] = surfacelevel;
            }
            water.icethickness_m[cell] = nearice[global] != 0 ? std::max(0.0f, ice[cell]) : 0.0f;
        }
    }
    return water;
}


// The Catmull-Rom weights of the four stencil points at the fraction t (upsampling.cpp's
// catmullrom written per point).
void catmullromweights(double t, double weights[4])
{
    const double t2 = t * t, t3 = t2 * t;
    weights[0] = 0.5 * (-t + 2.0 * t2 - t3);
    weights[1] = 0.5 * (2.0 - 5.0 * t2 + 3.0 * t3);
    weights[2] = 0.5 * (t + 4.0 * t2 - 3.0 * t3);
    weights[3] = 0.5 * (t3 - t2);
}

// R4's zone of a fine cell: the lake cell among its global cell and the 8 neighbours whose
// term in the bicubic indicator at the warped position (gx, gy), the stencil weight times
// the cell's magnitude, is the most negative; ties by the smaller global index (an absolute
// one). none: no lake cell around. A pure function of the fine cell's position.
std::size_t finezone(const RegionalGlobalInputs& inputs, const std::vector<uint8_t>& lakecell,
                     const std::vector<float>& magnitude, int i, int j, double gx, double gy)
{
    const int width = inputs.width, height = inputs.height;
    const double u = gx - 0.5, v = gy - 0.5;
    const double fu = std::floor(u), fv = std::floor(v);
    double wx[4], wy[4];
    catmullromweights(u - fu, wx);
    catmullromweights(v - fv, wy);
    const std::int64_t i0 = static_cast<std::int64_t>(fu), j0 = static_cast<std::int64_t>(fv);
    std::size_t best = none;
    double bestterm = 0.0;
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
        {
            const int column = wrapcolumn(i + di, width);
            const std::size_t cell = at(column, row, width);
            if (lakecell[cell] == 0)
                continue;
            // The cell's weight in the stencil (rows clamped at the poles may repeat it).
            double weight = 0.0;
            for (int b = 0; b < 4; b++)
            {
                const std::int64_t stencilrow = std::max<std::int64_t>(0, std::min<std::int64_t>(height - 1, j0 - 1 + b));
                if (stencilrow != row)
                    continue;
                for (int a = 0; a < 4; a++)
                    if (static_cast<int>(((i0 - 1 + a) % width + width) % width) == column)
                        weight += wx[a] * wy[b];
            }
            const double term = weight * std::fabs(static_cast<double>(magnitude[cell]));
            if (best == none || term > bestterm || (term == bestterm && cell < best))
            {
                best = cell;
                bestterm = term;
            }
        }
    }
    return best;
}

FineWater shorewater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                     const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                     const FineChannels& channels, std::vector<float>& elevation_m, bool globaldepth,
                     const LakeInheritance& inheritance)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const double shorenoise = parameters.shores.shorenoise;
    const double terrain = parameters.shores.terrain_m;
    const float above = std::numeric_limits<float>::infinity();
    const float below = -std::numeric_limits<float>::infinity();

    std::vector<uint8_t> lakecell(globalsize, 0), icecell(globalsize, 0), nearice(globalsize, 0);
    for (std::size_t global = 0; global < globalsize; global++)
    {
        lakecell[global] = inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0 ? 1 : 0;
        icecell[global] = inputs.icethickness_m[global] > 0.0f ? 1 : 0;
    }
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
            nearice[at(i, j, width)] =
                aroundcell(width, height, i, j, [&](std::size_t other) { return icecell[other] != 0; }) ? 1 : 0;

    // The indicators and the ice thickness at the fully warped positions; the noises are
    // evaluated only where a zone needs them.
    const std::vector<float> lakemagnitude = signedindicator(lakecell, width, height, inputs.seed, lakesalt);
    const std::vector<float> lakeindicator = upsample(lakemagnitude, width, height, window, warp);
    std::vector<float> lakemask(globalsize);
    for (std::size_t global = 0; global < globalsize; global++)
        lakemask[global] = lakecell[global] != 0 ? 1.0f : 0.0f;
    const std::vector<float> lakeweight = upsample(lakemask, width, height, window, warp);
    const std::vector<float> iceindicator =
        upsample(signedindicator(icecell, width, height, inputs.seed, icesalt), width, height, window, warp);
    const std::vector<float> ice = upsample(inputs.icethickness_m, width, height, window, warp);
    const double persistence = parameters.coast.noisepersistence;
    const RegionalNoise lakenoise(width, height, k, inputs.seed, lakenoisesalt, persistence);
    const RegionalNoise icenoise(width, height, k, inputs.seed, icenoisesalt, persistence);

    FineWater water;
    water.lakesurface_m.assign(count, 0.0f);
    water.icethickness_m.assign(count, 0.0f);

    // 1. Every fine cell that is not sea and has a lake cell around takes its zone's
    // surface; the sum is the indicator plus the shore noise plus the terrain term
    // clamp((z - surface) / terrain_m, -w, 1), w the lake cells' mask (1 on lake cells, 0
    // elsewhere) interpolated at the warped position and clamped to 0 .. 1: a cell far above
    // the surface is land whatever the indicator says, a cell far below it is water where the
    // lake cells dominate (the lake's deep part) and left to the indicator where they do not
    // (the valley beyond the lake's outlet, where forced water would run to the zone's edge
    // and stop at a wall); w is continuous, so the forcing makes no step at a block edge.
    // Water where the sum lies below zero and the elevation lies below the surface plus
    // terrain_m (so the lowering stays bounded).
    std::vector<std::size_t> zone(count, none);
    std::vector<float> surface(count, 0.0f);
    std::vector<double> lakesum(count, 0.0);
    std::vector<uint8_t> wet(count, 0);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        const std::int64_t finerow = window.absolutefinerow(y);
        const double gy = (static_cast<double>(finerow) + 0.5) / k;
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (sea[cell] != 0)
                continue;
            const std::int64_t finecolumn = window.absolutefinecolumn(x);
            const double gx = (static_cast<double>(finecolumn) + 0.5) / k;
            const std::size_t lake = finezone(inputs, lakecell, lakemagnitude, window.globalcolumn(x), row,
                                              gx + (warp.empty() ? 0.0 : warp.dx[cell]),
                                              gy + (warp.empty() ? 0.0 : warp.dy[cell]));
            if (lake == none)
                continue;
            zone[cell] = lake;
            surface[cell] = inputs.lakesurface_m[lake];
            const double height_m = static_cast<double>(elevation_m[cell]) - static_cast<double>(surface[cell]);
            const double forcing = std::clamp(static_cast<double>(lakeweight[cell]), 0.0, 1.0);
            const double term = terrain > 0.0 ? std::clamp(height_m / terrain, -forcing, 1.0) : 0.0;
            lakesum[cell] = static_cast<double>(lakeindicator[cell]) + shorenoise * lakenoise.at(finecolumn, finerow) + term;
            wet[cell] = lakesum[cell] < 0.0 && (terrain <= 0.0 || height_m < terrain) ? 1 : 0;
        }
    }

    // 1b. The refinement loop's flood (S2.8, LakeInheritance): the core (level 1: the water of
    // 1; level >= 2: the cells below S whose parent carries the zone's surface) and every cell
    // below S 4-connected to it through cells of the same zone below S, not sea, where the lake
    // cells dominate (the warped half-mask contour, lakeweight >= 0.5) or the global cell lies
    // above S, and at level >= 2 within one parent cell of the parent's water.
    if (inheritance.flood)
    {
        const bool inherited = inheritance.parentsurface_m != nullptr;
        // S3.4, the lake-sea contact at a sill: where the sea gate is open (regional_edges.hpp,
        // seagate, and at level >= 2 the parent's sea gate, so only where the level's own sea could
        // be), a cell between the sea level and S is passable whatever the lake mask says.
        const std::vector<uint8_t> gate = seagate(inputs, window, warp, parameters, &lakeweight);
        const auto passable = [&](std::size_t cell, int x, int y) {
            if (sea[cell] != 0 || zone[cell] == none || !(elevation_m[cell] < surface[cell]))
                return false;
            if (inheritance.startsurface_m != nullptr && !((*inheritance.startsurface_m)[cell] < surface[cell]))
                return false;
            if (inherited && (*inheritance.parentreach)[cell] == 0)
                return false;
            if (gate[cell] != 0 && (inheritance.parentseagate == nullptr || (*inheritance.parentseagate)[cell] != 0) &&
                elevation_m[cell] > inputs.sealevel)
                return true;
            const std::size_t global = at(window.globalcolumn(x), window.globalrow(y), width);
            return lakeweight[cell] >= 0.5f || inputs.elevation_m[global] > surface[cell];
        };
        std::vector<uint8_t> flooded(count, 0);
        std::vector<std::size_t> queue;
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                const bool core = inherited ? sea[cell] == 0 && zone[cell] != none &&
                                                  (*inheritance.parentsurface_m)[cell] == surface[cell] &&
                                                  elevation_m[cell] < surface[cell]
                                            : wet[cell] != 0;
                if (core)
                {
                    flooded[cell] = 1;
                    queue.push_back(cell);
                }
            }
        for (std::size_t next = 0; next < queue.size(); next++)
        {
            const std::size_t cell = queue[next];
            const int x = static_cast<int>(cell % static_cast<std::size_t>(domainwidth));
            const int y = static_cast<int>(cell / static_cast<std::size_t>(domainwidth));
            const int offsets[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
            for (const auto& offset : offsets)
            {
                const int nx = x + offset[0], ny = y + offset[1];
                if (nx < 0 || ny < 0 || nx >= domainwidth || ny >= domainheight)
                    continue;
                const std::size_t other = at(nx, ny, domainwidth);
                if (flooded[other] != 0 || zone[other] != zone[cell] || !passable(other, nx, ny))
                    continue;
                flooded[other] = 1;
                queue.push_back(other);
            }
        }
        wet = std::move(flooded);
    }

    // A wet 8-neighbour (inside the domain) of fine cell (x, y) passing `test`.
    auto neighbour = [&](int x, int y, auto test) {
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++)
            {
                const int nx = x + dx, ny = y + dy;
                if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= domainwidth || ny >= domainheight)
                    continue;
                if (test(at(nx, ny, domainwidth)))
                    return true;
            }
        return false;
    };

    // 2. No step between two lakes: a wet cell beside a wet cell of a lower surface is
    // dry (the dam between them, raised in 5), so adjacent water has one surface. Decided
    // on the cells of 1 only, so the order of the cells does not matter.
    std::vector<uint8_t> water2(count, 0);
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (wet[cell] != 0)
                water2[cell] = neighbour(x, y, [&](std::size_t other) {
                    return wet[other] != 0 && surface[other] < surface[cell];
                }) ? 0 : 1;
        }

    // 3. The channel cells of lake cells outside the lake of 1 and 2 (R4's replacement of
    // R1's rule, which wetted every channel cell at or below max(surface, filled) and drew
    // rivers through the dry part of a lake cell as one-cell lake lines): wet when the bed
    // lies below the surface plus terrain_m and an 8-neighbour that is not a channel is
    // wet, or when the bed lies strictly below the surface (no dry channel cell stands
    // below the water), unless a wet 8-neighbour has another surface (the outlet between
    // two lakes stays a river); a bed at or above the surface away from the water stays a
    // river. S3.13: in every block, a channel cell left dry whose lake level Λ (channellakelevels:
    // the lake its owner piece runs into or out of) is positive and whose bed lies at or below Λ
    // is wet at the surface Λ, under the same guard (no wet 8-neighbour of 2 with another surface);
    // step 4 lowers such a cell at Λ as every wet cell inside a lake cell's block whose zone has the
    // surface Λ, elsewhere to the next float below Λ only (outside the lake's block the zone's depth
    // would dig a trench; zone may be none there). Each candidate reads
    // the water of 2 and its own values only, so the order of the cells does not matter.
    const std::vector<float> lakelevel = channellakelevels(inputs, window, channels);
    std::vector<uint8_t> atlevel(count, 0);
    std::vector<uint8_t> water3 = water2;
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            if (water2[cell] != 0 || sea[cell] != 0 || channels.channel[cell] == 0)
                continue;
            if (zone[cell] != none && lakecell[global] != 0)
            {
                const double height_m = static_cast<double>(elevation_m[cell]) - static_cast<double>(surface[cell]);
                const bool meets = (terrain <= 0.0 || height_m < terrain) &&
                                   neighbour(x, y, [&](std::size_t other) { return water2[other] != 0 && channels.channel[other] == 0; });
                const bool otherlake = neighbour(x, y, [&](std::size_t other) {
                    return water2[other] != 0 && surface[other] != surface[cell];
                });
                water3[cell] = (meets || height_m < 0.0) && !otherlake ? 1 : 0;
                if (water3[cell] != 0)
                    continue;
            }
            const float level = lakelevel[cell];
            if (!(level > 0.0f) || !(elevation_m[cell] <= level))
                continue;
            if (neighbour(x, y, [&](std::size_t other) { return water2[other] != 0 && surface[other] != level; }))
                continue;
            // Inside a lake cell's block, a cell of the zone's own surface is lowered as every wet
            // cell of that block (step 4's depth); elsewhere only to the next float below Λ.
            atlevel[cell] = lakecell[global] != 0 && zone[cell] != none && surface[cell] == level ? 0 : 1;
            water3[cell] = 1;
            surface[cell] = level;
        }
    }

    // 3b. Step 3 checks a candidate against the water of 2 only, so two adjacent channel
    // cells it adds in lake cells of different zones (an outlet between two 8-adjacent lake
    // cells of different levels) could both be wet at different surfaces: every cell step 3
    // added that has a step-3-added 8-neighbour of a lower surface is dry again. Removal only,
    // decided on the cells of 3, so the order of the cells does not matter; the water of 2
    // is untouched, so no new step against it arises. Reach: one fine cell more than step 3
    // (whose candidates read the water of 2, itself one fine cell beyond the cells of 1), so
    // three fine cells in all, inside the halo.
    {
        std::vector<uint8_t> added(count, 0);
        for (std::size_t cell = 0; cell < count; cell++)
            added[cell] = water3[cell] != 0 && water2[cell] == 0 ? 1 : 0;
        for (int y = 0; y < domainheight; y++)
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (added[cell] != 0 && neighbour(x, y, [&](std::size_t other) {
                        return added[other] != 0 && surface[other] < surface[cell];
                    }))
                    water3[cell] = 0;
            }
    }

    // 3c. S3.4, the rim raise capped at the level's noise band (LakeInheritance::rimcap_m): a dry
    // cell (not sea, not a channel) beside the water of its own lake (8-neighbour, the zone's
    // surface S) that step 5 would raise by more than terrain_m and more than its cap joins the
    // water instead, within the flood's bounds (the lake's zone; below S on the start surface, so a
    // parent above S keeps a dry child; at level >= 2 within one parent cell of the parent's water)
    // and not beside water of a lower surface (steps 2 and 3b: adjacent water has one surface), in rounds until no cell joins (each round decided on the water
    // of the one before, so the order of the cells does not matter). Reach: inside the lake
    // cell's 3 x 3 zone, as the flood's. A deeper dry cell beyond those bounds keeps step 5's
    // raise (counted in rimraise_m).
    if (inheritance.flood && inheritance.rimcap_m != nullptr)
    {
        const std::vector<float>& cap = *inheritance.rimcap_m;
        for (;;)
        {
            std::vector<std::size_t> joins;
            for (int y = 0; y < domainheight; y++)
                for (int x = 0; x < domainwidth; x++)
                {
                    const std::size_t cell = at(x, y, domainwidth);
                    if (water3[cell] != 0 || sea[cell] != 0 || channels.channel[cell] != 0 || zone[cell] == none ||
                        (inheritance.parentreach != nullptr && (*inheritance.parentreach)[cell] == 0) ||
                        (inheritance.startsurface_m != nullptr && !((*inheritance.startsurface_m)[cell] < surface[cell])))
                        continue;
                    const double rise = static_cast<double>(surface[cell]) - static_cast<double>(elevation_m[cell]);
                    if (!(rise > terrain) || !(rise > static_cast<double>(cap[cell])))
                        continue;
                    if (neighbour(x, y, [&](std::size_t other) { return water3[other] != 0 && surface[other] < surface[cell]; }))
                        continue;
                    if (neighbour(x, y, [&](std::size_t other) { return water3[other] != 0 && surface[other] == surface[cell]; }))
                        joins.push_back(cell);
                }
            if (joins.empty())
                break;
            for (const std::size_t cell : joins)
                water3[cell] = 1;
        }
    }

    // 4. The beds: a wet cell at or above its surface is lowered to the surface minus a
    // depth that grows from 0 where the sum is 0 (the shore) to min(the global lake depth
    // (surface minus the zone's base elevation), terrain_m) where the sum reaches minus the
    // zone's magnitude, and to at least the next float below the surface; a wet cell
    // already below its surface keeps its elevation. So no cell is lowered by more than
    // 2 x terrain_m (it lay below the surface plus terrain_m). With `globaldepth` (the
    // refinement loop, multiscale plan S2.3) the depth is the global lake depth uncapped.
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (water3[cell] == 0)
            continue;
        const std::size_t lake = zone[cell];
        const float level = surface[cell];
        water.lakesurface_m[cell] = level;
        if (elevation_m[cell] < level)
            continue;
        if (atlevel[cell] != 0)
        {
            elevation_m[cell] = std::nextafter(level, below);  // S3.13: a channel cell at its lake level
            continue;
        }
        const double magnitude = std::fabs(static_cast<double>(lakemagnitude[lake]));
        double depth = std::max(0.0, static_cast<double>(level) - static_cast<double>(inputs.elevation_m[lake]));
        if (terrain > 0.0 && !globaldepth)
            depth = std::min(depth, terrain);
        const double share = std::clamp(-lakesum[cell] / magnitude, 0.0, 1.0);
        elevation_m[cell] = std::min(std::nextafter(level, below), static_cast<float>(static_cast<double>(level) - share * depth));
    }

    // 5. No dry land just below an adjacent lake's surface: a dry fine cell (not sea, not a
    // channel) at or below the level L, the highest surface of its wet 8-neighbours and of its
    // own global cell if that is a lake cell, and at most terrain_m below it, is raised to L
    // plus its depth below L (mirrored about L, at least the next float above it), so the
    // raised shore keeps its relief instead of a flat at the surface and rises by at most
    // 2 x terrain_m; a dry cell of a lake cell with no water beside it is raised to the mirror
    // of its own lake's surface the same way. A cell more than terrain_m below L is left: the
    // deep bed of a pit lake outside the indicator's outline (the pit's interpolated bed and
    // inner slope, its non-lake neighbours at or above the surface), left dry rather than
    // dammed by hundreds of metres (forcing water on all of a lake cell's own fine cells is a
    // later option).
    // The wet cells are final, so the order of the cells does not matter; a neighbour beyond
    // the domain is not seen (within the halo).
    if (inheritance.flood)
        water.rimraise_m.assign(count, 0.0f);
    std::vector<float> raised(count, 0.0f);
    std::vector<uint8_t> raise(count, 0);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (water.lakesurface_m[cell] > 0.0f || sea[cell] != 0 || channels.channel[cell] != 0)
                continue;
            const std::size_t global = at(window.globalcolumn(x), row, width);
            float level = below, beside = below;
            if (lakecell[global] != 0)
                level = inputs.lakesurface_m[global];
            neighbour(x, y, [&](std::size_t other) {
                if (water.lakesurface_m[other] > 0.0f)
                    beside = std::max(beside, water.lakesurface_m[other]);
                return false;
            });
            level = std::max(level, beside);
            if (!(elevation_m[cell] <= level))
                continue;
            const double rise = static_cast<double>(level) - static_cast<double>(elevation_m[cell]);
            if (terrain > 0.0 && rise > terrain)
            {
                // S2.8's flood: a deeper dry cell beside the water (outside the flood's reach) is
                // raised to just above its surface, so no dry cell stands below the water beside it.
                if (inheritance.flood && elevation_m[cell] <= beside)
                {
                    raised[cell] = std::nextafter(beside, above);
                    raise[cell] = 1;
                    water.rimraise_m[cell] =
                        static_cast<float>(static_cast<double>(raised[cell]) - static_cast<double>(elevation_m[cell]));
                }
                continue;
            }
            raised[cell] = std::max(std::nextafter(level, above), static_cast<float>(static_cast<double>(level) + rise));
            raise[cell] = 1;
        }
    }
    for (std::size_t cell = 0; cell < count; cell++)
        if (raise[cell] != 0)
            elevation_m[cell] = raised[cell];

    // The ice: where the ice indicator's sum lies below zero on the fine cells of ice cells
    // and their neighbours, the thickness at the warped position, at least iceminimum_m; on a
    // fine sea cell only where the global cell is iced sea (the fine coast is the calving
    // front).
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        const std::int64_t finerow = window.absolutefinerow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            if (nearice[global] == 0 || (sea[cell] != 0 && !(inputs.network.sea[global] != 0 && icecell[global] != 0)))
                continue;
            const double sum = static_cast<double>(iceindicator[cell]) +
                               shorenoise * icenoise.at(window.absolutefinecolumn(x), finerow);
            if (sum < 0.0)
                water.icethickness_m[cell] = std::max(iceminimum_m, ice[cell]);
        }
    }
    return water;
}
}

FineWater regionalwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const FineWarp& warp,
                        const RegionalParameters& parameters, const std::vector<uint8_t>& sea,
                        const FineChannels& channels, std::vector<float>& elevation_m, bool globaldepth,
                        const LakeInheritance& inheritance)
{
    if (!(parameters.shores.shorenoise > 0.0))
        return phaseewater(inputs, window, sea, channels, elevation_m);
    const std::size_t count =
        static_cast<std::size_t>(window.domainfinewidth()) * static_cast<std::size_t>(window.domainfineheight());
    if ((inheritance.parentsurface_m == nullptr) != (inheritance.parentreach == nullptr) ||
        (inheritance.parentsurface_m != nullptr &&
         (inheritance.parentsurface_m->size() != count || inheritance.parentreach->size() != count)))
        throw std::invalid_argument("regionalwater: the parent's lake masks must both hold one value per domain cell");
    if ((inheritance.parentseagate != nullptr && inheritance.parentseagate->size() != count) ||
        (inheritance.rimcap_m != nullptr && inheritance.rimcap_m->size() != count))
        throw std::invalid_argument("regionalwater: the parent sea gate and the rim cap must hold one value per domain cell");
    if (inheritance.startsurface_m != nullptr && inheritance.startsurface_m->size() != count)
        throw std::invalid_argument("regionalwater: the start surface must hold one value per domain cell");
    return shorewater(inputs, window, warp, parameters, sea, channels, elevation_m, globaldepth, inheritance);
}

std::vector<float> channellakelevels(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                     const FineChannels& channels)
{
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const std::size_t globalsize = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    std::vector<float> level(count, 0.0f);
    if (channels.channel.size() != count || channels.globalcell.size() != count || inputs.lakesurface_m.size() != globalsize ||
        inputs.network.sea.size() != globalsize || inputs.network.direction.size() != globalsize)
        return level;
    const auto lake = [&](std::size_t global) {
        return inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
    };
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (channels.channel[cell] == 0 || channels.globalcell[cell] < 0 ||
                static_cast<std::size_t>(channels.globalcell[cell]) >= globalsize)
                continue;
            const std::size_t upstream = static_cast<std::size_t>(channels.globalcell[cell]);
            const std::size_t next = globalreceiver(inputs.network, upstream);
            if (next != none && lake(next))
                level[cell] = inputs.lakesurface_m[next];
            else if (lake(upstream) && at(window.globalcolumn(x), row, inputs.width) == upstream)
                level[cell] = inputs.lakesurface_m[upstream];
        }
    }
    return level;
}

std::vector<std::size_t> channelmeansafterlakes(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                                const FineChannels& channels, const std::vector<uint8_t>& sea,
                                                const std::vector<float>& bed_m, const std::vector<float>& means_m,
                                                std::vector<float>& elevation_m, FineWater& water, bool drymeans)
{
    std::vector<std::size_t> dried;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    if (channels.channel.size() != count || sea.size() != count || bed_m.size() != count || means_m.size() != count ||
        elevation_m.size() != count || water.lakesurface_m.size() != count)
        return dried;
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            if (channels.channel[cell] == 0 || sea[cell] != 0 || means_m[cell] == bed_m[cell])
                continue;
            const float surface = water.lakesurface_m[cell];
            if (surface > 0.0f)
            {
                const std::size_t global = at(window.globalcolumn(x), row, inputs.width);
                if (inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                    continue;  // a lake bed of its own lake cell's block
            }
            elevation_m[cell] = means_m[cell];
            if (drymeans && surface > 0.0f && means_m[cell] > surface)
            {
                water.lakesurface_m[cell] = 0.0f;  // land at this level: its mean stands above the water
                dried.push_back(cell);
            }
        }
    }
    return dried;
}

void lakeslevel(const LevelContext& context, LevelState& state)
{
    // 5. Lakes and ice (regional_pass.cpp's step 5, moved here unchanged in S3.0).
    LakeInheritance lakes = context.details.lakes;
    if (lakes.flood)
    {
        lakes.startsurface_m = &state.surface_m;  // the start surface with the beds pinned (S2.8)
        if (state.noiseband_m.size() == state.surface_m.size())
            lakes.rimcap_m = &state.noiseband_m;  // S3.4: the rim raise capped at the level's noise band
        lakes.parentseagate = context.details.parentseagate;  // S3.4: the sill rule within the level's sea gate
    }
    // S3.2: the lakes read a land channel cell at its bed (state.bed_m), not at its cell mean (the
    // refinement loop's refinement.channel_cell_slope, regional_solve.hpp), so the water decision reads
    // the bed; then channelmeansafterlakes (S3.13). With every land channel cell at its bed (the single
    // pass, theta 0) this is regionalwater on state.elevation_m itself.
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
    state.water = regionalwater(context.inputs, context.window, context.warp, context.parameters, state.sea,
                                state.channels, elevation, context.details.globallakedepth, lakes);
    if (withbeds)
    {
        const std::vector<std::size_t> dried = channelmeansafterlakes(context.inputs, context.window, state.channels, state.sea,
                                                                      state.bed_m, means, elevation, state.water,
                                                                      context.details.coarselevel);
        if (state.diagnostics != nullptr && context.levelwindow != nullptr)
        {
            const LevelWindow& level = *context.levelwindow;
            const int domainwidth = context.window.domainfinewidth();
            const int x0 = level.keptcolumnoffset(), y0 = level.keptrowoffset();
            for (const std::size_t cell : dried)
            {
                const int x = static_cast<int>(cell % static_cast<std::size_t>(domainwidth));
                const int y = static_cast<int>(cell / static_cast<std::size_t>(domainwidth));
                if (x >= x0 && x < x0 + level.columns && y >= y0 && y < y0 + level.rows)
                    state.diagnostics->driedchannels++;
            }
        }
    }
    recordlevelstep(state, LevelStep::lakes);
}
}
