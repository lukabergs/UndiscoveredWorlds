// Tests of regional_water.hpp (erosion plan Phase R, step R4), suite Landscape:
// TEST_F(Landscape, RegionalWater<What>) and RegionalOverlapShores on the synthetic world of
// regional_test_world.hpp with the warp, R2's coast and the shore noise on (test values,
// shoreparameters()), at the smallest halo the pass allows (`testhalo`, two block-mean
// passes): the shores are pure functions of the fine cell's absolute position, and the
// neighbourhood rules (the step between two lakes, the channel cells, the raised shore)
// read up to three fine cells. Multiscale plan S3.4 adds the refinement loop's rim cap
// (RegionalWaterRimCap), the lake-sea contact at a sill (RegionalWaterSillContact) and both at a kept
// edge cutting a lake's zone (RegionalWaterOverlapSillAndRimCap).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_edges.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_water.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <limits>
#include <utility>

namespace
{
constexpr double testterrain_m = 50.0;

RegionalParameters shoreparameters(double shorenoise)
{
    RegionalParameters parameters;
    parameters.warp.amplitude_cells = 0.5;
    parameters.warp.wavelength_cells = 2.0;
    parameters.warp.fade_m = 400.0;
    parameters.coast.noiseamplitude_m = 150.0;
    parameters.coast.noisepersistence = 0.5;
    parameters.coast.noiserelieffraction = 0.25;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.shores.shorenoise = shorenoise;
    parameters.shores.terrain_m = testterrain_m;
    return parameters;
}

bool lakecell(const RegionalGlobalInputs& inputs, std::size_t global)
{
    return inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
}

// True when the global cell of kept fine cell (x, y) or one of its 8 neighbours passes `test`.
template<typename Test>
bool aroundglobal(const RegionalGlobalInputs& inputs, const RegionalWindow& window, int x, int y, Test test)
{
    const int i = (window.column0 + x / window.scale) % inputs.width, j = window.row0 + y / window.scale;
    for (int dj = -1; dj <= 1; dj++)
        for (int di = -1; di <= 1; di++)
            if (j + dj >= 0 && j + dj < inputs.height && test(at((i + di + inputs.width) % inputs.width, j + dj, inputs.width)))
                return true;
    return false;
}

// True when a kept 8-neighbour of kept fine cell (x, y) passes `test`.
template<typename Test>
bool aroundfine(const RegionalFields& fields, int x, int y, Test test)
{
    const int width = fields.window.finewidth(), height = fields.window.fineheight();
    for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++)
        {
            const int nx = x + dx, ny = y + dy;
            if ((dx != 0 || dy != 0) && nx >= 0 && ny >= 0 && nx < width && ny < height && test(at(nx, ny, width)))
                return true;
        }
    return false;
}

bool channel(const RegionalFields& fields, std::size_t cell) { return fields.riverdirection[cell] != 0; }
bool wet(const RegionalFields& fields, std::size_t cell) { return fields.lakesurface_m[cell] > 0.0f; }

// What the kept window's lakes and ice hold and which of R4's rules they break.
struct WaterCheck
{
    int lake = 0, lakeoutside = 0, lakeonsea = 0, bedabove = 0, drybelow = 0, drychannel = 0, lonechannel = 0;
    int drymeeting = 0;   // dry channel cells of lake cells beside a wet cell that is not a channel (and no lower water)
    int drybeside = 0;    // dry land (not sea, not a channel) at or below a wet 8-neighbour's surface and at most
                          // terrain_m below it, in any block (the raise's reach)
    int perched = 0;      // dry land more than terrain_m below a wet 8-neighbour's or its lake cell's surface (left)
    int steps = 0;        // 8-adjacent pairs of wet cells with different surfaces
    int ice = 0, iceoutside = 0, iceonsea = 0, seaice = 0;
    int lakeblocks = 0;                  // global lake cells in the kept window
    int lakeedges = 0, iceedges = 0;     // 4-neighbour pairs of land cells across a shore or an ice edge
    double lakeshare = 0.0, iceshare = 0.0;  // the share of them on global block edges
};

WaterCheck checkwater(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalFields& fields)
{
    WaterCheck check;
    const int k = window.scale, width = window.finewidth(), height = window.fineheight();
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = at(x, y, width);
            const std::size_t global = globalof(window, inputs.width, x, y);
            const bool sea = fields.sea[cell] != 0;
            const float surface = fields.lakesurface_m[cell];
            const float elevation = fields.elevation_m[cell];
            if (x % k == 0 && y % k == 0 && lakecell(inputs, global))
                check.lakeblocks++;
            if (surface > 0.0f)
            {
                check.lake++;
                if (!aroundglobal(inputs, window, x, y, [&](std::size_t other) { return lakecell(inputs, other); }))
                    check.lakeoutside++;
                check.lakeonsea += sea ? 1 : 0;
                check.bedabove += elevation >= surface ? 1 : 0;
                check.steps += aroundfine(fields, x, y, [&](std::size_t other) {
                    return wet(fields, other) && fields.lakesurface_m[other] < surface;
                }) ? 1 : 0;
                if (channel(fields, cell))
                    check.lonechannel += aroundfine(fields, x, y, [&](std::size_t other) {
                        return wet(fields, other) && !channel(fields, other);
                    }) ? 0 : 1;
            }
            else if (!sea)
            {
                float level = -std::numeric_limits<float>::infinity();
                aroundfine(fields, x, y, [&](std::size_t other) {
                    if (wet(fields, other))
                        level = std::max(level, fields.lakesurface_m[other]);
                    return false;
                });
                const double below = static_cast<double>(level) - static_cast<double>(elevation);
                if (!channel(fields, cell))
                {
                    check.drybeside += below >= 0.0 && below <= testterrain_m ? 1 : 0;
                    check.perched += below > testterrain_m ? 1 : 0;
                    if (lakecell(inputs, global))
                    {
                        const double own = static_cast<double>(inputs.lakesurface_m[global]) - static_cast<double>(elevation);
                        check.drybelow += own >= 0.0 && own <= testterrain_m ? 1 : 0;
                        check.perched += own > testterrain_m && !(below > testterrain_m) ? 1 : 0;
                    }
                }
                else if (lakecell(inputs, global))
                {
                    // R4's channel rule: a dry channel cell of a lake cell lies at or above the
                    // surface, and is not beside the water, unless it is the outlet between two
                    // lakes: beside water of another surface, or among lake cells of different
                    // surfaces (its zone, the surface it would take, may be a neighbour's lake).
                    const float level = inputs.lakesurface_m[global];
                    const bool otherlake = aroundfine(fields, x, y, [&](std::size_t other) {
                        return wet(fields, other) && fields.lakesurface_m[other] != level;
                    }) || aroundglobal(inputs, window, x, y, [&](std::size_t other) {
                        return lakecell(inputs, other) && inputs.lakesurface_m[other] != level;
                    });
                    if (!otherlake)
                    {
                        check.drychannel += elevation < level ? 1 : 0;
                        check.drymeeting += elevation < level + testterrain_m && aroundfine(fields, x, y, [&](std::size_t other) {
                            return wet(fields, other) && !channel(fields, other);
                        }) ? 1 : 0;
                    }
                }
            }
            if (fields.icethickness_m[cell] > 0.0f)
            {
                check.ice++;
                if (!aroundglobal(inputs, window, x, y, [&](std::size_t other) { return inputs.icethickness_m[other] > 0.0f; }))
                    check.iceoutside++;
                if (sea && !(inputs.network.sea[global] != 0 && inputs.icethickness_m[global] > 0.0f))
                    check.iceonsea++;
                check.seaice += sea ? 1 : 0;
            }
        }
    int lakeblock = 0, iceblock = 0;
    auto pair = [&](std::size_t a, std::size_t b, bool blockedge) {
        if (fields.sea[a] != 0 || fields.sea[b] != 0)
            return;
        if (wet(fields, a) != wet(fields, b))
        {
            check.lakeedges++;
            lakeblock += blockedge ? 1 : 0;
        }
        if ((fields.icethickness_m[a] > 0.0f) != (fields.icethickness_m[b] > 0.0f))
        {
            check.iceedges++;
            iceblock += blockedge ? 1 : 0;
        }
    };
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            if (x + 1 < width)
                pair(at(x, y, width), at(x + 1, y, width), (x + 1) % k == 0);
            if (y + 1 < height)
                pair(at(x, y, width), at(x, y + 1, width), (y + 1) % k == 0);
        }
    check.lakeshare = check.lakeedges > 0 ? static_cast<double>(lakeblock) / check.lakeedges : 0.0;
    check.iceshare = check.iceedges > 0 ? static_cast<double>(iceblock) / check.iceedges : 0.0;
    std::printf("[regional] water k %d: lake %d fine cells (%d global lake cells, area ratio %.2f), %d outside the zone, %d on sea, "
                "%d beds at or above, %d dry below, %d dry channels below, %d dry channels beside the lake, %d lone wet "
                "channels, %d dry beside water, %d perched, %d steps; ice %d, %d outside, %d on sea (%d off iced sea cells); block-edge "
                "share lake %.3f of %d, ice %.3f of %d\n",
                k, check.lake, check.lakeblocks, check.lakeblocks > 0 ? static_cast<double>(check.lake) / (check.lakeblocks * k * k) : 0.0,
                check.lakeoutside, check.lakeonsea, check.bedabove, check.drybelow, check.drychannel, check.drymeeting,
                check.lonechannel, check.drybeside, check.perched, check.steps, check.ice, check.iceoutside, check.seaice, check.iceonsea,
                check.lakeshare, check.lakeedges, check.iceshare, check.iceedges);
    return check;
}

void expectrules(const WaterCheck& check)
{
    EXPECT_EQ(check.lakeoutside, 0);
    EXPECT_EQ(check.lakeonsea, 0);
    EXPECT_EQ(check.bedabove, 0);
    EXPECT_EQ(check.drybelow, 0);
    EXPECT_EQ(check.drychannel, 0);
    EXPECT_EQ(check.drymeeting, 0);
    EXPECT_EQ(check.drybeside, 0);
    EXPECT_EQ(check.steps, 0);
    EXPECT_EQ(check.iceoutside, 0);
    EXPECT_EQ(check.iceonsea, 0);
}

// The elevation the shores change against the pass with shorenoise 0 and otherwise the same
// parameters (the same solve): cells lowered and raised, the largest of each, the largest
// change of a global cell's block mean, the lowered wet cells' mean depth below their surface
// at the shore (a dry land 4-neighbour) and inside, and raised cells outside lake blocks.
struct Reshape
{
    int lowered = 0, raised = 0, raisedoutside = 0, shorecells = 0, innercells = 0;
    double largestlowering = 0.0, largestraise = 0.0, largestblock = 0.0, shoredepth = 0.0, innerdepth = 0.0;
};

Reshape reshape(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalFields& before,
                const RegionalFields& after)
{
    Reshape result;
    const int k = window.scale, width = window.finewidth(), height = window.fineheight();
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = at(x, y, width);
            const double change = static_cast<double>(after.elevation_m[cell]) - static_cast<double>(before.elevation_m[cell]);
            if (change < 0.0)
            {
                result.lowered++;
                result.largestlowering = std::max(result.largestlowering, -change);
                if (wet(after, cell) && !channel(after, cell))
                {
                    const double depth = static_cast<double>(after.lakesurface_m[cell]) - static_cast<double>(after.elevation_m[cell]);
                    bool shore = false;
                    for (const auto& [dx, dy] : {std::pair{1, 0}, std::pair{-1, 0}, std::pair{0, 1}, std::pair{0, -1}})
                    {
                        const int nx = x + dx, ny = y + dy;
                        if (nx >= 0 && ny >= 0 && nx < width && ny < height && !wet(after, at(nx, ny, width)) &&
                            after.sea[at(nx, ny, width)] == 0)
                            shore = true;
                    }
                    (shore ? result.shoredepth : result.innerdepth) += depth;
                    (shore ? result.shorecells : result.innercells)++;
                }
            }
            else if (change > 0.0)
            {
                result.raised++;
                result.largestraise = std::max(result.largestraise, change);
                result.raisedoutside += lakecell(inputs, globalof(window, inputs.width, x, y)) ? 0 : 1;
            }
        }
    for (int r = 0; r < window.rows; r++)
        for (int c = 0; c < window.columns; c++)
        {
            double sum = 0.0;
            for (int y = r * k; y < (r + 1) * k; y++)
                for (int x = c * k; x < (c + 1) * k; x++)
                    sum += static_cast<double>(after.elevation_m[at(x, y, width)]) - static_cast<double>(before.elevation_m[at(x, y, width)]);
            result.largestblock = std::max(result.largestblock, std::fabs(sum) / (k * k));
        }
    result.shoredepth = result.shorecells > 0 ? result.shoredepth / result.shorecells : 0.0;
    result.innerdepth = result.innercells > 0 ? result.innerdepth / result.innercells : 0.0;
    std::printf("[regional] reshape k %d: %d lowered (largest %.1f m), %d raised (largest %.1f m, %d outside lake blocks), "
                "block means up to %.2f m; lowered wet depth: shore %.2f m over %d, inside %.2f m over %d\n",
                k, result.lowered, result.largestlowering, result.raised, result.largestraise, result.raisedoutside,
                result.largestblock, result.shoredepth, result.shorecells, result.innerdepth, result.innercells);
    return result;
}

// The lake of the synthetic world, (62, 17) on the 64 x 32 world, and a second lake cell two
// cells away (a non-lake cell between them; window C holds all three): the first land cell,
// rows 15 .. 19 then columns 60 .. 64, whose lowest neighbour lies above the first lake's
// surface, made a pit 30 m deep whose surface is that lowest neighbour's elevation (its spill
// level, as a routing would give it). `second` receives its index (or -1: none found).
RegionalGlobalInputs twolakes(std::size_t& second)
{
    RegionalGlobalInputs inputs = syntheticworld(64);
    const std::size_t first = at(62, 17, 64);
    second = static_cast<std::size_t>(-1);
    for (int j = 15; j <= 19 && second == static_cast<std::size_t>(-1); j++)
        for (int i = 60; i <= 64 && second == static_cast<std::size_t>(-1); i++)
        {
            if (std::max(std::abs(i - 62), std::abs(j - 17)) != 2)
                continue;
            const std::size_t cell = at(i % 64, j, 64);
            if (inputs.network.sea[cell] != 0)
                continue;
            float lowest = std::numeric_limits<float>::infinity();
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++)
                    if (di != 0 || dj != 0)
                        lowest = std::min(lowest, inputs.elevation_m[at((i + di + 64) % 64, j + dj, 64)]);
            if (!(lowest > inputs.lakesurface_m[first] + 1.0f))
                continue;
            second = cell;
            inputs.lakesurface_m[cell] = lowest;
            inputs.elevation_m[cell] = lowest - 30.0f;
            inputs.network.filledelevation_m[cell] = lowest;
        }
    return inputs;
}

RegionalGlobalInputs twolakes()
{
    std::size_t second = 0;
    return twolakes(second);
}

// The synthetic world with one iced sea cell: the first sea cell of window C's kept rows
// (west to east, north to south) with land among its 8 neighbours, 100 m of ice.
RegionalGlobalInputs icedsea(std::size_t& iced)
{
    RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = windowc(8);
    iced = static_cast<std::size_t>(-1);
    for (int j = window.row0; j < window.row0 + window.rows && iced == static_cast<std::size_t>(-1); j++)
        for (int c = 0; c < window.columns && iced == static_cast<std::size_t>(-1); c++)
        {
            const int i = (window.column0 + c) % 64;
            const std::size_t cell = at(i, j, 64);
            if (inputs.network.sea[cell] == 0)
                continue;
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++)
                    if (inputs.network.sea[at((i + di + 64) % 64, j + dj, 64)] == 0)
                        iced = cell;
        }
    if (iced != static_cast<std::size_t>(-1))
        inputs.icethickness_m[iced] = 100.0f;
    return inputs;
}
}

TEST_F(Landscape, RegionalWaterShores)
{
    // Window C at k = 16 holds the lake (62, 17), ice on the ridge and a coast. With the
    // shore noise on: lake water only around lake cells, never on sea, beds below the
    // surface, no dry land at or below the surface of its lake cell or of adjacent water in
    // any block, no step between adjacent water, R4's channel rule; ice only around ice
    // cells; both outlines off the global block edges (the rules before R4 on the same window
    // are the reference); against those rules the elevation is reshaped by a bounded amount.
    const RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = windowc(16);
    const RegionalFields off = timedpass(inputs, window, "shores off", shoreparameters(0.0));
    const RegionalFields on = timedpass(inputs, window, "shores on", shoreparameters(0.3));
    const WaterCheck before = checkwater(inputs, window, off);
    const WaterCheck after = checkwater(inputs, window, on);
    ASSERT_GT(after.lake, 0);
    ASSERT_GT(after.ice, 0);
    expectrules(after);
    // The lake's area of the order of its global cell's (one cell with a hashed magnitude
    // among neighbours with their own); the judging windows check the faithfulness of many.
    const double area = static_cast<double>(after.lakeblocks) * window.scale * window.scale;
    EXPECT_GT(after.lake, area / 3.0);
    EXPECT_LT(after.lake, area * 2.0);
    // The share of the outlines on block edges (1/k = 0.0625 for an outline that ignores them):
    // the synthetic lake is a pit 60 m deep among higher cells, so the terrain term draws most
    // of its shore (0.180 against 0.234 at R4); the ice has no terrain term.
    EXPECT_LT(after.lakeshare, before.lakeshare);
    EXPECT_LT(after.iceshare, 0.75 * before.iceshare);
    // The bed law: a wet cell lay below the surface plus terrain_m and is lowered to at most
    // terrain_m below it, so no cell is lowered by 2 x terrain_m or more; the raised shore
    // rises by its depth below the level, at most terrain_m above it.
    const Reshape change = reshape(inputs, window, off, on);
    EXPECT_LT(change.largestlowering, 2.0 * testterrain_m);
    EXPECT_LT(change.largestraise, 2.0 * testterrain_m);
}

TEST_F(Landscape, RegionalWaterBedLaw)
{
    // The depth grows from the shore: on a wide shallow lake (the 3 x 3 cells around (62, 17)
    // lake cells at the lake's surface, their base 10 m below it, so the regional noise leaves
    // parts of them at or above the surface) the lowered wet cells beside dry land lie less
    // deep below the surface than the lowered ones inside the lake, and none deeper than the
    // global lake depth (window C at k = 16).
    RegionalGlobalInputs inputs = syntheticworld(64);
    const float level = inputs.lakesurface_m[at(62, 17, 64)];
    for (int dj = -1; dj <= 1; dj++)
        for (int di = -1; di <= 1; di++)
        {
            const std::size_t cell = at(62 + di, 17 + dj, 64);
            inputs.lakesurface_m[cell] = level;
            inputs.elevation_m[cell] = level - 10.0f;
            inputs.network.filledelevation_m[cell] = level;
        }
    const RegionalWindow window = windowc(16);
    const RegionalFields off = timedpass(inputs, window, "bed law off", shoreparameters(0.0));
    const RegionalFields on = timedpass(inputs, window, "bed law on", shoreparameters(0.15));
    const Reshape change = reshape(inputs, window, off, on);
    ASSERT_GT(change.shorecells, 0);
    ASSERT_GT(change.innercells, 0);
    EXPECT_LT(change.shoredepth, change.innerdepth);
    EXPECT_LE(change.innerdepth, 10.0 + 0.01);
    EXPECT_LT(change.largestlowering, 2.0 * testterrain_m);
}

TEST_F(Landscape, RegionalWaterTwoLakes)
{
    // Two lake cells 40 m apart in level with a non-lake cell between them: each fine cell
    // takes the lake whose indicator term is the strongest at it, and water of the higher
    // lake beside water of the lower one is dry (a dam), so no two adjacent wet cells carry
    // different surfaces anywhere (inside the non-lake block between them included); both
    // lakes keep water, and the dams are raised above the higher surface.
    std::size_t secondlake = 0;
    const RegionalGlobalInputs inputs = twolakes(secondlake);
    ASSERT_NE(secondlake, static_cast<std::size_t>(-1));
    const float lower = inputs.lakesurface_m[at(62, 17, 64)], higher = inputs.lakesurface_m[secondlake];
    std::printf("[regional] two lakes: (62, 17) at %.1f m, (%d, %d) at %.1f m\n", lower, static_cast<int>(secondlake % 64),
                static_cast<int>(secondlake / 64), higher);
    ASSERT_LT(lower, higher);
    for (int scale : {8, 16})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowc(scale);
        const RegionalFields fields = timedpass(inputs, window, "two lakes", shoreparameters(0.3));
        const WaterCheck check = checkwater(inputs, window, fields);
        expectrules(check);
        int first = 0, second = 0, between = 0;
        for (int y = 0; y < window.fineheight(); y++)
            for (int x = 0; x < window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, window.finewidth());
                first += fields.lakesurface_m[cell] == lower ? 1 : 0;
                second += fields.lakesurface_m[cell] == higher ? 1 : 0;
                const std::size_t global = globalof(window, 64, x, y);
                between += !lakecell(inputs, global) && wet(fields, cell) ? 1 : 0;
            }
        std::printf("[regional] two lakes k %d: %d cells at the lower surface, %d at the higher, %d wet outside lake cells\n",
                    scale, first, second, between);
        EXPECT_GT(first, 0);
        EXPECT_GT(second, 0);
    }
}

TEST_F(Landscape, RegionalWaterAdjacentLakes)
{
    // The lake (62, 17) and its receiver, 8-adjacent, both lake cells, the receiver's surface
    // 5 m lower; their terrain 300 m above the surfaces (so the terrain term leaves no water of
    // the indicator) and their filled elevations 10 m below their surfaces (so the outlet's
    // channel bed, pinned at or below them, lies below both surfaces): the channel rule wets
    // the channel cells of both lake cells, each at its own surface, so adjacent channel
    // cells across the cells' edge would carry different surfaces; the higher stay rivers
    // (pass 3b), and no two adjacent wet cells carry different surfaces (window D at k = 8 and
    // 16 holds both cells).
    RegionalGlobalInputs inputs = syntheticworld(64);
    const std::size_t first = at(62, 17, 64);
    const int direction = inputs.network.direction[first];
    ASSERT_NE(direction, 0);
    const int ri = (62 + landscape::directionoffsetx(direction) + 64) % 64, rj = 17 + landscape::directionoffsety(direction);
    const std::size_t receiver = at(ri, rj, 64);
    ASSERT_EQ(inputs.network.sea[receiver], 0);
    const float upperlevel = inputs.lakesurface_m[first];
    const float level = upperlevel - 5.0f;
    inputs.lakesurface_m[receiver] = level;
    inputs.network.filledelevation_m[first] = upperlevel - 10.0f;
    inputs.network.filledelevation_m[receiver] = level - 10.0f;
    inputs.elevation_m[first] = upperlevel + 300.0f;
    inputs.elevation_m[receiver] = upperlevel + 300.0f;
    for (int scale : {8, 16})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowd(scale);
        const RegionalFields fields = timedpass(inputs, window, "adjacent lakes", shoreparameters(0.3));
        const WaterCheck check = checkwater(inputs, window, fields);
        EXPECT_EQ(check.steps, 0);
        EXPECT_EQ(check.bedabove, 0);
        EXPECT_EQ(check.lakeonsea, 0);
        EXPECT_EQ(check.drychannel, 0);
        int upper = 0, lower = 0, channels = 0;
        for (int y = 0; y < window.fineheight(); y++)
            for (int x = 0; x < window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, window.finewidth());
                const std::size_t global = globalof(window, 64, x, y);
                upper += fields.lakesurface_m[cell] == inputs.lakesurface_m[first] ? 1 : 0;
                lower += fields.lakesurface_m[cell] == level ? 1 : 0;
                channels += (global == first || global == receiver) && channel(fields, cell) ? 1 : 0;
            }
        std::printf("[regional] adjacent lakes k %d: receiver (%d, %d), %d cells at the upper surface, %d at the lower, "
                    "%d channel cells in the two lake cells\n",
                    scale, ri, rj, upper, lower, channels);
        EXPECT_GT(upper, 0);
        EXPECT_GT(lower, 0);
        EXPECT_GT(channels, 0);
    }
}

TEST_F(Landscape, RegionalWaterShoresBelowSill)
{
    // The lake closed 20 m below its sill (R1's review), where R1's rule wetted every channel
    // cell of the lake cell up to the sill: with the shore noise on, the rules hold (no dry
    // channel cell strictly below the surface or beside the water).
    const RegionalGlobalInputs inputs = lakebelowsill(syntheticworld(64), 20.0f);
    const RegionalWindow window = windowc(8);
    const RegionalFields fields = timedpass(inputs, window, "shores below sill", shoreparameters(0.3));
    const WaterCheck check = checkwater(inputs, window, fields);
    ASSERT_GT(check.lake, 0);
    expectrules(check);
}

TEST_F(Landscape, RegionalWaterSeaIce)
{
    // One iced sea cell beside the coast: some of its fine sea cells carry ice, and no fine
    // sea cell of another global cell does (the fine coast is the calving front elsewhere).
    std::size_t iced = 0;
    const RegionalGlobalInputs inputs = icedsea(iced);
    ASSERT_NE(iced, static_cast<std::size_t>(-1));
    const RegionalWindow window = windowc(8);
    const RegionalFields fields = timedpass(inputs, window, "sea ice", shoreparameters(0.3));
    const WaterCheck check = checkwater(inputs, window, fields);
    EXPECT_GT(check.seaice, 0);
    EXPECT_EQ(check.iceonsea, 0);
    EXPECT_EQ(check.iceoutside, 0);
}

TEST_F(Landscape, RegionalWaterPhaseERules)
{
    // shorenoise 0 keeps the rules before R4 exactly, with the warp and R2's coast on: ice is
    // the unwarped upsampled thickness, at least 0, around the ice cells; lake water on the
    // non-channel fine cells of lake cells below the surface and on the channel cells R1 wets
    // (lowered to at most the surface), nowhere else.
    const RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = windowc(8);
    const RegionalFields fields = timedpass(inputs, window, "shores off", shoreparameters(0.0));
    const std::vector<float> ice = keptpart(landscape::upsample(inputs.icethickness_m, inputs.width, inputs.height, window), window);
    int lake = 0, wrong = 0, icewrong = 0;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
        {
            const std::size_t cell = at(x, y, window.finewidth());
            const std::size_t global = globalof(window, inputs.width, x, y);
            const float surface = inputs.lakesurface_m[global];
            const float found = fields.lakesurface_m[cell];
            const float elevation = fields.elevation_m[cell];
            if (fields.sea[cell] != 0 || !(surface > 0.0f))
                wrong += found != 0.0f ? 1 : 0;
            else if (fields.riverdirection[cell] == 0)
            {
                const float expected = elevation < surface ? surface : 0.0f;
                wrong += std::memcmp(&expected, &found, sizeof(float)) != 0 ? 1 : 0;
            }
            else if (found > 0.0f)
                wrong += found != surface || elevation > surface ? 1 : 0;  // R1: lowered to at most the surface
            else
                wrong += elevation <= std::max(surface, inputs.network.filledelevation_m[global]) ? 1 : 0;
            lake += found > 0.0f ? 1 : 0;
            const bool nearice = aroundglobal(inputs, window, x, y, [&](std::size_t other) { return inputs.icethickness_m[other] > 0.0f; });
            const float expectedice = nearice ? std::max(0.0f, ice[cell]) : 0.0f;
            icewrong += std::memcmp(&expectedice, &fields.icethickness_m[cell], sizeof(float)) != 0 ? 1 : 0;
        }
    EXPECT_GT(lake, 0);
    EXPECT_EQ(wrong, 0);
    EXPECT_EQ(icewrong, 0);
}

TEST_F(Landscape, RegionalOverlapShores)
{
    // The pass with the shore noise, the warp and R2's coast on is bit-identical on the
    // overlap of two windows, in columns across the date line and in rows clear of the
    // poles (the overlap of C and D holds the lake and ice), at the smallest halo; also with
    // two lakes of different levels and with the lake below its sill.
    const RegionalParameters parameters = shoreparameters(0.3);
    const RegionalGlobalInputs columns = syntheticworld(32);
    const RegionalGlobalInputs rows = syntheticworld(64);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts first = compareoverlap(columns, windowa(scale), windowb(scale), "shores A/B", parameters);
        EXPECT_EQ(first.shared, 4 * 5 * scale * scale);
        const OverlapCounts second = compareoverlap(rows, windowc(scale), windowd(scale), "shores C/D", parameters);
        EXPECT_EQ(second.shared, 12 * 4 * scale * scale);
        EXPECT_GT(second.lake, 0);
        EXPECT_GT(second.ice, 0);
        EXPECT_GT(second.channel, 0);
    }
    const OverlapCounts sill = compareoverlap(lakebelowsill(rows, 20.0f), windowc(8), windowd(8), "shores below sill C/D", parameters);
    EXPECT_GT(sill.lake, 0);
    const OverlapCounts two = compareoverlap(twolakes(), windowc(8), windowd(8), "shores two lakes C/D", parameters);
    EXPECT_GT(two.lake, 0);
}

namespace
{
// The 64-column synthetic world with its lake cell (62, 17) between a sea cell on its west, (61, 17)
// at -100 m, and a low land cell on its east, (63, 17) at 5 m; the lake's surface 30 m above the sea
// level, its bed at 10 m. The east cell lies below the surface, so beyond the lake's warped half-mask
// its fine cells are the flood's barrier; the west cell gates the sea (the sill rule).
constexpr float barrierlakesurface = 30.0f;

RegionalGlobalInputs barrierlake()
{
    RegionalGlobalInputs inputs = syntheticworld(64);
    const int width = inputs.width;
    inputs.network.sea[at(61, 17, width)] = 1;
    inputs.elevation_m[at(61, 17, width)] = -100.0f;
    inputs.elevation_m[at(63, 17, width)] = 5.0f;
    inputs.elevation_m[at(62, 17, width)] = 10.0f;
    inputs.lakesurface_m[at(62, 17, width)] = barrierlakesurface;
    return inputs;
}

// A window of the 64-column world keeping global columns c0 .. c1 (c1 may pass the date line) and rows
// r0 .. r1 at the refinement loop's smallest halo.
RegionalWindow barrierwindow(int scale, int c0, int c1, int r0, int r1)
{
    const auto longitude = [](int i) { return -180.0 + 360.0 * ((i % 64) + 0.5) / 64.0; };
    const auto latitude = [](int j) { return 90.0 - 180.0 * (j + 0.5) / 32.0; };
    return landscape::regionalwindow(64, scale, landscape::refinementhalo(2, 0), longitude(c0), latitude(r1),
                                     longitude(c1), latitude(r0));
}

// The loop's level-1 sea and lakes on one window (lakegate on, the flood with the start surface): the
// fine cells of (63, 17) with a lake weight below 0.5 and a start surface below S, not sea, are lowered
// 300 m below S (a pure function of the absolute position), so only 3c can wet them; a cap below 0:
// none.
struct BarrierPass
{
    std::vector<float> start, weight, elevation;
    std::vector<uint8_t> sea, gate, pit;
    landscape::FineWater water;
};

BarrierPass barrierpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window, double cap)
{
    RegionalParameters parameters = shoreparameters(0.3);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    BarrierPass pass;
    pass.start = landscape::regionalstartsurface(inputs, window, warp, parameters, true).start_m;
    pass.sea = landscape::regionalsea(inputs, window, warp, parameters, pass.start, nullptr, true);
    pass.weight = landscape::lakecellweight(inputs, window, warp);
    pass.gate = landscape::seagate(inputs, window, warp, parameters, &pass.weight);
    const std::size_t count = pass.start.size();
    pass.elevation = pass.start;
    pass.pit.assign(count, 0);
    const std::size_t east = at(63, 17, inputs.width);
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x < window.domainfinewidth(); x++)
        {
            const std::size_t cell = at(x, y, window.domainfinewidth());
            if (at(window.globalcolumn(x), window.globalrow(y), inputs.width) == east && pass.weight[cell] < 0.5f &&
                pass.sea[cell] == 0 && pass.start[cell] < barrierlakesurface)
            {
                pass.elevation[cell] = barrierlakesurface - 300.0f;
                pass.pit[cell] = 1;
            }
        }
    landscape::FineChannels channels;
    channels.channel.assign(count, 0);
    const std::vector<float> caps(count, static_cast<float>(cap));
    landscape::LakeInheritance flood;
    flood.flood = true;
    flood.startsurface_m = &pass.start;
    flood.rimcap_m = cap >= 0.0 ? &caps : nullptr;
    pass.water = landscape::regionalwater(inputs, window, warp, parameters, pass.sea, channels, pass.elevation, true, flood);
    return pass;
}

struct BarrierCounts
{
    int pits = 0, wetpits = 0, sill = 0, wet = 0;
    float largestraise = 0.0f;
};

// The lowered cells and how many are wet; the wet cells only the sill rule lets the flood reach (gated,
// lake weight below 0.5, between the sea level and S, in the sea cell's block); the wet cells; the
// largest rim raise of the lake's zone (global columns 61 .. 63, rows 16 .. 18).
BarrierCounts barriercounts(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const BarrierPass& pass)
{
    BarrierCounts counts;
    const std::size_t west = at(61, 17, inputs.width);
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x < window.domainfinewidth(); x++)
        {
            const std::size_t cell = at(x, y, window.domainfinewidth());
            const int i = window.globalcolumn(x), j = window.globalrow(y);
            const bool wet = pass.water.lakesurface_m[cell] > 0.0f;
            counts.wet += wet ? 1 : 0;
            counts.pits += pass.pit[cell];
            counts.wetpits += pass.pit[cell] != 0 && wet ? 1 : 0;
            counts.sill += wet && at(i, j, inputs.width) == west && pass.gate[cell] != 0 && pass.weight[cell] < 0.5f &&
                                   pass.elevation[cell] > sealevel
                               ? 1
                               : 0;
            if (i >= 61 && i <= 63 && j >= 16 && j <= 18 && !pass.water.rimraise_m.empty())
                counts.largestraise = std::max(counts.largestraise, pass.water.rimraise_m[cell]);
        }
    return counts;
}
}

// Multiscale plan S3.4, the rim cap (the refinement loop's flood, LakeInheritance::rimcap_m), on window
// C at k 8 (barrierlake): the lowered cells of (63, 17) lie below S on the start surface and behind the
// flood's barrier. Without a cap the rim raise lifts those beside the water about 300 m; a cap of 1e9 m
// changes nothing; with a cap of 100 m they join the water (step 3c) and no cell of the lake's zone is
// raised by more than the cap.
TEST_F(Landscape, RegionalWaterRimCap)
{
    const RegionalGlobalInputs inputs = barrierlake();
    const RegionalWindow window = windowc(8);
    const BarrierPass none = barrierpass(inputs, window, -1.0);
    const BarrierPass huge = barrierpass(inputs, window, 1.0e9);
    const BarrierPass capped = barrierpass(inputs, window, 100.0);
    const std::size_t count = none.elevation.size();
    EXPECT_EQ(std::memcmp(none.elevation.data(), huge.elevation.data(), count * sizeof(float)), 0);
    EXPECT_EQ(std::memcmp(none.water.lakesurface_m.data(), huge.water.lakesurface_m.data(), count * sizeof(float)), 0);
    EXPECT_EQ(std::memcmp(none.water.rimraise_m.data(), huge.water.rimraise_m.data(), count * sizeof(float)), 0);
    const BarrierCounts a = barriercounts(inputs, window, none), b = barriercounts(inputs, window, capped);
    std::printf("[regional] rim cap: %d lowered cells; uncapped: %d wet, rim raise %.1f m, %d wet cells; capped at 100 m: "
                "%d wet, rim raise %.1f m, %d wet cells; sill-only wet cells %d\n",
                a.pits, a.wetpits, a.largestraise, a.wet, b.wetpits, b.largestraise, b.wet, b.sill);
    ASSERT_GT(a.pits, 0);
    EXPECT_EQ(a.wetpits, 0);
    EXPECT_GT(a.largestraise, 250.0f);
    EXPECT_GT(b.wetpits, 0);
    EXPECT_LE(b.largestraise, 100.0f);
    EXPECT_GT(b.wet, a.wet);
}

// Multiscale plan S3.4 review: two windows whose kept edge cuts the lake's zone at the loop's smallest
// halo (refinementhalo), E keeping columns 58 .. 62 (the zone's column 63 in its halo) and F columns
// 62 .. 66 across the date line (column 61 in its halo), rows 14 .. 19 and 15 .. 20, with the sill
// rule and 3c both acting (barrierlake, cap 20 m): the elevation, the lake surface and the rim raise
// bit-identical on the shared kept cells.
TEST_F(Landscape, RegionalWaterOverlapSillAndRimCap)
{
    const RegionalGlobalInputs inputs = barrierlake();
    const RegionalWindow e = barrierwindow(8, 58, 62, 14, 19);
    const RegionalWindow f = barrierwindow(8, 62, 66, 15, 20);
    ASSERT_EQ(e.column0, 58);
    ASSERT_EQ(e.columns, 5);
    ASSERT_EQ(f.column0, 62);
    ASSERT_EQ(f.row0, 15);
    const BarrierPass a = barrierpass(inputs, e, 20.0), b = barrierpass(inputs, f, 20.0);
    const BarrierCounts ca = barriercounts(inputs, e, a), cb = barriercounts(inputs, f, b);
    const std::int64_t finewidth = static_cast<std::int64_t>(inputs.width) * 8;
    int shared = 0, differ = 0, wet = 0;
    for (int y = 0; y < e.fineheight(); y++)
        for (int x = 0; x < e.finewidth(); x++)
        {
            const int dx = x + e.keptfinecolumnoffset(), dy = y + e.keptfinerowoffset();
            const std::int64_t column = e.absolutefinecolumn(dx), row = e.absolutefinerow(dy);
            const std::int64_t fx = ((column - f.absolutefinecolumn(0)) % finewidth + finewidth) % finewidth;
            const std::int64_t fy = row - f.absolutefinerow(0);
            if (fx < f.keptfinecolumnoffset() || fx >= f.keptfinecolumnoffset() + f.finewidth() || fy < f.keptfinerowoffset() ||
                fy >= f.keptfinerowoffset() + f.fineheight())
                continue;
            const std::size_t cellf = at(static_cast<int>(fx), static_cast<int>(fy), f.domainfinewidth());
            const std::size_t celle = at(dx, dy, e.domainfinewidth());
            shared++;
            differ += std::memcmp(&a.elevation[celle], &b.elevation[cellf], sizeof(float)) != 0 ? 1 : 0;
            differ += std::memcmp(&a.water.lakesurface_m[celle], &b.water.lakesurface_m[cellf], sizeof(float)) != 0 ? 1 : 0;
            differ += std::memcmp(&a.water.rimraise_m[celle], &b.water.rimraise_m[cellf], sizeof(float)) != 0 ? 1 : 0;
            wet += a.water.lakesurface_m[celle] > 0.0f ? 1 : 0;
        }
    std::printf("[regional] sill and rim cap overlap E/F: %d shared cells, %d differences, %d wet; E: %d of %d lowered "
                "cells wet, %d sill-only wet; F: %d of %d lowered cells wet, %d sill-only wet\n",
                shared, differ, wet, ca.wetpits, ca.pits, ca.sill, cb.wetpits, cb.pits, cb.sill);
    EXPECT_EQ(shared, 8 * 5 * 8);
    EXPECT_EQ(differ, 0);
    EXPECT_GT(wet, 0);
    EXPECT_GT(ca.wetpits, 0);
    EXPECT_GT(cb.wetpits, 0);
    EXPECT_GT(ca.sill, 0);
    EXPECT_GT(cb.sill, 0);
}

namespace
{
// The 256-column synthetic world with a lake at its sill height beside the sea (S3.4): the first
// coastal land cell of rows 58 .. 66, columns 214 .. 230 (the continent's west coast, about 132 E)
// whose west neighbour is sea and east neighbour land becomes a lake cell 2 m above the sea level
// with its bed 10 m below it, as global cell (169, 60) of erosion run 9.
RegionalGlobalInputs coastallake(int& lakei, int& lakej)
{
    RegionalGlobalInputs inputs = syntheticworld(256);
    lakei = lakej = -1;
    for (int j = 58; j <= 66 && lakei < 0; j++)
        for (int i = 214; i <= 230 && lakei < 0; i++)
        {
            const std::size_t cell = at(i, j, inputs.width);
            if (inputs.network.sea[cell] == 0 && inputs.lakesurface_m[cell] <= 0.0f &&
                inputs.network.sea[at(i - 1, j, inputs.width)] != 0 && inputs.network.sea[at(i + 1, j, inputs.width)] == 0)
            {
                lakei = i;
                lakej = j;
                inputs.lakesurface_m[cell] = sealevel + 2.0f;
                inputs.elevation_m[cell] = sealevel - 10.0f;
            }
        }
    return inputs;
}

// A window of the 256-column world keeping global columns c0 .. c1 and rows r0 .. r1 (centres in degrees).
RegionalWindow keptwindow(int c0, int c1, int r0, int r1)
{
    const auto longitude = [](int i) { return -180.0 + 360.0 * (i + 0.5) / 256.0; };
    const auto latitude = [](int j) { return 90.0 - 180.0 * (j + 0.5) / 128.0; };
    return landscape::regionalwindow(256, 16, landscape::refinementhalo(2, 0), longitude(c0), latitude(r1), longitude(c1),
                                     latitude(r0));
}

landscape::RefinementResult refinecoast(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                        const RegionalParameters& parameters)
{
    landscape::LevelTarget target;
    target.window = window;
    target.level = 2;
    return landscape::runrefinement(inputs, target, parameters, landscape::RefinementParameters{});
}

// The shared kept cells of two level windows, as index pairs into their kept fields.
std::vector<std::pair<std::size_t, std::size_t>> sharedlevelcells(const landscape::LevelWindow& a, const landscape::LevelWindow& b)
{
    std::vector<std::pair<std::size_t, std::size_t>> shared;
    const std::int64_t width = a.levelwidth();
    for (int y = 0; y < a.rows; y++)
        for (int x = 0; x < a.columns; x++)
        {
            const std::int64_t column = landscape::wraplevel(a.column0 + x, width), row = a.row0 + y;
            const std::int64_t xb = landscape::wraplevel(column - b.column0, width), yb = row - b.row0;
            if (xb < b.columns && yb >= 0 && yb < b.rows)
                shared.push_back({at(x, y, a.columns), at(static_cast<int>(xb), static_cast<int>(yb), b.columns)});
        }
    return shared;
}
}

// Multiscale plan S3.4, the lake-sea contact at a sill, in the refinement loop (two levels at
// factor 4, k 16): a lake cell 2 m above the sea level beside a sea cell (coastallake) meets the sea
// off the block edge between them (the gate's lake exclusion on the warped half-mask, the flood
// across the sea-gated block): at most a quarter of the level-2 rows of the shared block edge hold
// a lake cell beside a sea cell across it (the lake edges on global cell boundaries are only printed:
// this lake fills its block, whose land neighbours rise above S at the block edge). Two windows with
// different column and row origins holding the lake and the sea are bit-identical on their overlap
// at both levels. (The overlap across the date line
// with the lakes and the sea on is RefinementLevelsOverlap's, P/Q and R/S, whose loop runs S3.4's
// rules too.)
TEST_F(Landscape, RegionalWaterSillContact)
{
    int lakei = -1, lakej = -1;
    const RegionalGlobalInputs inputs = coastallake(lakei, lakej);
    ASSERT_GE(lakei, 0);
    RegionalParameters parameters = shoreparameters(0.15);
    parameters.edges.smoothgate = 3;
    const RegionalWindow first = keptwindow(lakei - 4, lakei + 3, lakej - 3, lakej + 3);
    const RegionalWindow second = keptwindow(lakei - 3, lakei + 4, lakej - 2, lakej + 4);
    ASSERT_EQ(first.column0, lakei - 4);
    ASSERT_EQ(first.row0, lakej - 3);
    const landscape::RefinementResult a = refinecoast(inputs, first, parameters);
    const landscape::RefinementResult b = refinecoast(inputs, second, parameters);
    ASSERT_EQ(a.levels.size(), 2u);
    const landscape::LevelDiagnostics& two = a.levels[1];
    const RegionalFields& fields = a.fields;
    const int columns = two.window.columns, k = 16;
    const int x0 = (lakei - first.column0) * k, y0 = (lakej - first.row0) * k;
    int contact = 0;
    for (int y = y0; y < y0 + k; y++)
    {
        const std::size_t west = at(x0 - 1, y, columns), east = at(x0, y, columns);
        const bool wetwest = fields.lakesurface_m[west] > 0.0f, weteast = fields.lakesurface_m[east] > 0.0f;
        contact += (weteast && fields.sea[west] != 0) || (wetwest && fields.sea[east] != 0) ? 1 : 0;
    }
    int lake = 0;
    for (int y = y0; y < y0 + k; y++)
        for (int x = x0; x < x0 + k; x++)
            lake += fields.lakesurface_m[at(x, y, columns)] > 0.0f ? 1 : 0;
    std::printf("[regional] sill contact at global cell (%d, %d): level 2 lake edges %d, on boundaries %.4f (2 / k %.4f); "
                "%d wet cells in the lake cell; lake beside sea across the west block edge in %d of %d rows\n",
                lakei, lakej, two.lakeedges, two.lakeedgealigned, 2.0 / k, lake, contact, k);
    EXPECT_GT(lake, 0);
    EXPECT_GT(two.lakeedges, 0);
    EXPECT_LE(contact, k / 4);

    for (std::size_t level = 0; level < 2; level++)
    {
        const auto cells = sharedlevelcells(a.levels[level].window, b.levels[level].window);
        ASSERT_GT(cells.size(), 0u);
        int differ = 0, sea = 0, wet = 0;
        for (const auto& [ca, cb] : cells)
        {
            if (level == 0)
                differ += std::memcmp(&a.coarserelevation_m[0][ca], &b.coarserelevation_m[0][cb], sizeof(float)) != 0 ? 1 : 0;
            else
            {
                differ += std::memcmp(&a.fields.elevation_m[ca], &b.fields.elevation_m[cb], sizeof(float)) != 0 ? 1 : 0;
                differ += a.fields.sea[ca] != b.fields.sea[cb] ? 1 : 0;
                differ += std::memcmp(&a.fields.lakesurface_m[ca], &b.fields.lakesurface_m[cb], sizeof(float)) != 0 ? 1 : 0;
                sea += a.fields.sea[ca] != 0 ? 1 : 0;
                wet += a.fields.lakesurface_m[ca] > 0.0f ? 1 : 0;
            }
        }
        std::printf("[regional] sill contact overlap, level %zu: %zu shared cells, %d differences, %d sea, %d wet\n", level + 1,
                    cells.size(), differ, sea, wet);
        EXPECT_EQ(differ, 0) << "level " << level + 1;
        if (level == 1)
        {
            EXPECT_GT(sea, 0);
            EXPECT_GT(wet, 0);
        }
    }
}
