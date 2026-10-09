// Tests of regional_water.hpp (erosion plan Phase R, step R4), suite Landscape:
// TEST_F(Landscape, RegionalWater<What>) and RegionalOverlapShores on the synthetic world of
// regional_test_world.hpp with the warp, R2's coast and the shore noise on (test values,
// shoreparameters()), at the smallest halo the pass allows (`testhalo`, two block-mean
// passes): the shores are pure functions of the fine cell's absolute position, and the
// neighbourhood rules (the step between two lakes, the channel cells, the raised shore)
// read up to three fine cells.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/upsampling.hpp"

#include <limits>

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
