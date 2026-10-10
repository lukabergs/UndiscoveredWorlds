// Tests of landscape/regional_tributaries.hpp (multiscale plan S3.0: the hook; S3.3 grows the
// tributaries), suite Landscape: TEST_F(Landscape, RegionalTributaries<What>) on the synthetic
// world of regional_test_world.hpp: off, the hook changes nothing; on, at level 1 (through the
// hook interface) every tributary lies inside its junction's parent cell and its 8 neighbours,
// drains into a channel without touching another, with beds and discharge that never rise
// downstream and never stand above the start surface; at level 2 (a hand-made parent) the
// inherited tributaries drain without loops into links within two global cells; the loop at two
// and three levels and the single pass are bit-identical on two overlapping windows with
// tributaries on at their smallest halos, and the density rises with the level; the outlets keep
// their discharge except where a tributary crosses a global cell's edge; off, the loop ignores the
// other tributary keys.

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_tributaries.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/sphere_metric.hpp"

#include "refinement_test_windows.hpp"

namespace
{
using landscape::LevelDiagnostics;
using landscape::LevelTarget;
using landscape::LevelWindow;
using landscape::RefinementParameters;
using landscape::RefinementResult;

// The level's steps up to the rasterised channels, through the hook interface.
void runtochannels(const landscape::LevelContext& context, landscape::LevelState& state)
{
    landscape::noiselevel(context, state);
    landscape::presmoothinglevel(context, state);
    landscape::coastlevel(context, state);
    landscape::incisionlevel(context, state);
    // As channelslevel: the incised network and bed surface with the incision active.
    const bool incised = context.incision.active;
    state.channels = landscape::rasterisechannels(context.window, incised ? state.incisednetwork : context.inputs.network,
                                                  context.parameters.channels, context.inputs.sealevel,
                                                  incised ? state.bedsurface_m : state.surface_m, state.sea);
}

RegionalParameters tributarieson()
{
    RegionalParameters parameters;
    parameters.tributaries.enabled = 1;
    return parameters;
}

// R9's incision and valleys at level 1 (refinement_levels_tests.cpp's valleysatlevelone), with
// tributaries on.
RegionalParameters tributarieswithvalleys()
{
    RegionalParameters parameters = tributarieson();
    parameters.relief.elapsed_myr = 20.0;
    parameters.valleys.upliftscale = 1.0;
    parameters.valleys.denudationratescale = 1.0;
    parameters.valleys.floor_km = 100.0;
    parameters.valleys.reference_m3s = 1.0e3;
    parameters.valleys.floormax_km = 150.0;
    parameters.valleys.sideslope = 0.02;
    parameters.valleys.reachcells = 4;
    return parameters;
}

// The 256-column world with R9's fields (refinement_levels_tests.cpp's wideworld).
RegionalGlobalInputs wideworld()
{
    RegionalGlobalInputs inputs = syntheticworld(256);
    const std::size_t size = static_cast<std::size_t>(inputs.width) * static_cast<std::size_t>(inputs.height);
    inputs.upliftrate_m_myr.resize(size);
    inputs.denudation_m.resize(size);
    inputs.erosionduration_myr = 20.0;
    for (int j = 0; j < inputs.height; j++)
        for (int i = 0; i < inputs.width; i++)
        {
            const double c = std::cos(latitude_rad(j, inputs.height));
            inputs.upliftrate_m_myr[at(i, j, inputs.width)] = static_cast<float>(40.0 * c * c - 5.0);
            inputs.denudation_m[at(i, j, inputs.width)] = static_cast<float>(300.0 + 200.0 * std::sin(longitude_rad(i, inputs.width)));
        }
    return withglacialfields(inputs);
}

// P (160 E .. 180, 15 S .. 10 N) and Q (170 E .. 170 W across the date line, 10 S .. 15 N):
// different column and row origins (refinement_levels_tests.cpp's windows).
// The loop's smallest halo with tributaries on (refinementhalo).
constexpr int loophalo = landscape::refinementhalo(2, 0, true);
RegionalWindow windowp(int scale, int halo = loophalo) { return landscape::regionalwindow(256, scale, halo, 160.0, -15.0, 180.0, 10.0); }
RegionalWindow windowq(int scale, int halo = loophalo) { return landscape::regionalwindow(256, scale, halo, 170.0, -10.0, -170.0, 15.0); }
// Smaller boxes for k 64: P' (170 E .. 180, 5 S .. 5 N) and Q' (175 E .. 175 W across the date line,
// 0 .. 10 N), sharing 175 E .. 180 by 0 .. 5 N.
RegionalWindow windowpsmall(int scale) { return landscape::regionalwindow(256, scale, loophalo, 170.0, -5.0, 180.0, 5.0); }
RegionalWindow windowqsmall(int scale) { return landscape::regionalwindow(256, scale, loophalo, 175.0, 0.0, -175.0, 10.0); }

RefinementResult refine(const RegionalGlobalInputs& inputs, const RegionalWindow& window, int level,
                        const RegionalParameters& parameters, const std::string& name)
{
    LevelTarget target;
    target.window = window;
    target.level = level;
    const auto start = std::chrono::steady_clock::now();
    RefinementResult result = landscape::runrefinement(inputs, target, parameters, RefinementParameters{});
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("[tributaries] %s: L %d, %.3f s\n", name.c_str(), level, seconds);
    for (const LevelDiagnostics& diagnostics : result.levels)
        std::printf("[tributaries]   level %d (%.2f km): %d grown, %d tributary cells, %d channel of %d land cells, "
                    "drainage density %.5f per km; forced %d, border step max %.1f m\n",
                    diagnostics.level, diagnostics.cellsize_km, diagnostics.tributaries, diagnostics.tributarycells,
                    diagnostics.channelcells, diagnostics.landcells, diagnostics.drainagedensity_per_km,
                    diagnostics.projection.forcedparents, diagnostics.projection.maxborderstep_m);
    return result;
}

std::vector<std::pair<std::size_t, std::size_t>> sharedcells(const LevelWindow& a, const LevelWindow& b)
{
    std::vector<std::pair<std::size_t, std::size_t>> shared;
    const std::int64_t width = a.levelwidth();
    for (int y = 0; y < a.rows; y++)
        for (int x = 0; x < a.columns; x++)
        {
            const std::int64_t column = landscape::wraplevel(a.column0 + x, width), row = a.row0 + y;
            const std::int64_t xb = landscape::wraplevel(column - b.column0, width), yb = row - b.row0;
            if (xb >= b.columns || yb < 0 || yb >= b.rows)
                continue;
            shared.push_back({at(x, y, a.columns), at(static_cast<int>(xb), static_cast<int>(yb), b.columns)});
        }
    return shared;
}

template<typename T>
int differing(const std::vector<T>& a, const std::vector<T>& b, const std::vector<std::pair<std::size_t, std::size_t>>& cells)
{
    int count = 0;
    for (const auto& [ca, cb] : cells)
        count += std::memcmp(&a[ca], &b[cb], sizeof(T)) != 0 ? 1 : 0;
    return count;
}

// The final level's channel cells in the shared kept box (no comparison).
int sharedchannels(const RefinementResult& a, const RefinementResult& b)
{
    int channel = 0;
    for (const auto& [ca, cb] : sharedcells(a.levels.back().window, b.levels.back().window))
        channel += a.fields.riverdirection[ca] != 0 ? 1 : 0;
    static_cast<void>(b);
    return channel;
}

// Every level's projected elevation and every field of the final level bit-identical on the overlap;
// returns the final level's shared channel cells.
int compareoverlaplevels(const RefinementResult& a, const RefinementResult& b, const std::string& label)
{
    const int levels = static_cast<int>(a.levels.size());
    for (int level = 1; level < levels; level++)
    {
        const auto cells = sharedcells(a.levels[static_cast<std::size_t>(level - 1)].window,
                                       b.levels[static_cast<std::size_t>(level - 1)].window);
        EXPECT_GT(cells.size(), 0u) << label << ", level " << level;
        EXPECT_EQ(differing(a.coarserelevation_m[static_cast<std::size_t>(level - 1)],
                            b.coarserelevation_m[static_cast<std::size_t>(level - 1)], cells),
                  0)
            << label << ": level " << level << "'s elevation";
    }
    const auto cells = sharedcells(a.levels.back().window, b.levels.back().window);
    const RegionalFields& fa = a.fields;
    const RegionalFields& fb = b.fields;
    int channel = 0;
    for (const auto& [ca, cb] : cells)
        channel += fa.riverdirection[ca] != 0 ? 1 : 0;
    std::printf("[tributaries] %s: %zu shared cells, %d channel\n", label.c_str(), cells.size(), channel);
    EXPECT_GT(channel, 0) << label;
    EXPECT_EQ(differing(fa.elevation_m, fb.elevation_m, cells), 0) << label << ": elevation_m";
    EXPECT_EQ(differing(fa.sea, fb.sea, cells), 0) << label << ": sea";
    EXPECT_EQ(differing(fa.riverdirection, fb.riverdirection, cells), 0) << label << ": riverdirection";
    EXPECT_EQ(differing(fa.discharge_m3s, fb.discharge_m3s, cells), 0) << label << ": discharge_m3s";
    EXPECT_EQ(differing(fa.januaryflow, fb.januaryflow, cells), 0) << label << ": januaryflow";
    EXPECT_EQ(differing(fa.julyflow, fb.julyflow, cells), 0) << label << ": julyflow";
    EXPECT_EQ(differing(fa.lakesurface_m, fb.lakesurface_m, cells), 0) << label << ": lakesurface_m";
    EXPECT_EQ(differing(fa.icethickness_m, fb.icethickness_m, cells), 0) << label << ": icethickness_m";
    EXPECT_EQ(differing(fa.climate, fb.climate, cells), 0) << label << ": climate";
    return channel;
}
}

TEST_F(Landscape, RegionalTributariesOffLeaveTheChannelsUnchanged)
{
    // tributaries.enabled 0 (the default): the step changes no channel field and counts none.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    ASSERT_EQ(parameters.tributaries.enabled, 0);
    const RegionalWindow window = windowa(4);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const landscape::LevelContext context{inputs, window, parameters, warp, incision, details};
    landscape::LevelState state;
    runtochannels(context, state);
    ASSERT_GT(std::count(state.channels.channel.begin(), state.channels.channel.end(), uint8_t{1}), 0);
    const landscape::FineChannels before = state.channels;
    const std::vector<float> surface = state.surface_m;
    state.tributaries = 7;
    state.tributarycells = 9;

    landscape::tributarieslevel(context, state);
    EXPECT_EQ(state.channels.channel, before.channel);
    EXPECT_EQ(state.channels.direction, before.direction);
    EXPECT_EQ(state.channels.bed_m, before.bed_m);
    EXPECT_EQ(state.channels.discharge_m3s, before.discharge_m3s);
    EXPECT_EQ(state.channels.januaryflow, before.januaryflow);
    EXPECT_EQ(state.channels.julyflow, before.julyflow);
    EXPECT_EQ(state.channels.globalcell, before.globalcell);
    EXPECT_EQ(state.surface_m, surface);
    EXPECT_EQ(state.tributaries, 0);
    EXPECT_EQ(state.tributarycells, 0);
}

TEST_F(Landscape, RegionalTributariesGrowInsideTheirParentCell)
{
    // Level 1 of the loop at factor 4 (k 4: the parent cell is the global cell), through the hook:
    // every tributary cell drains in fewer than 4 steps into a land link cell J of the global
    // network whose global cell is its own or one of its 8 neighbours; it touches no channel cell
    // but its own path and, for the cell next to J, cells within one cell of J; no diagonal step of
    // a tributary crosses a channel's; beds and discharges never rise downstream; the orders are
    // Strahler's (1 on the grown tributaries, which no tributary joins at level 1).
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = tributarieswithvalleys();
    const RegionalWindow window = windowp(4);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const landscape::LevelContext context{inputs, window, parameters, warp, incision, details, 1, 4};
    landscape::LevelState state;
    runtochannels(context, state);
    const landscape::FineChannels links = state.channels;
    const landscape::TributaryReport report = landscape::growtributaries(context, state);
    std::printf("[tributaries] level 1 of P: %d candidates, %d slot winners, %d grown (%d cells), max order %d\n",
                report.candidates, report.slotwinners, report.grown, report.growncells, report.maxorder);
    ASSERT_GT(report.grown, 50);
    EXPECT_EQ(state.tributaries, report.grown);
    EXPECT_EQ(state.tributarycells, report.growncells);
    EXPECT_EQ(report.inheritedcells, 0);
    EXPECT_EQ(report.maxorder, 1);
    EXPECT_EQ(report.bedsabove, 0);
    EXPECT_EQ(state.tributarybedsabove, 0);

    const landscape::FineChannels& channels = state.channels;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const auto index = [&](int x, int y) { return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x); };
    const auto down = [&](int x, int y, int& nx, int& ny) {
        const int direction = channels.direction[index(x, y)];
        nx = x + landscape::directionoffsetx(direction);
        ny = y + landscape::directionoffsety(direction);
        return direction >= 1 && direction <= 8 && nx >= 0 && ny >= 0 && nx < width && ny < height;
    };
    const std::vector<uint8_t> order = landscape::tributaryorders(window, channels, state.sea);
    int tributarycells = 0, failures = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = index(x, y);
            // Links unchanged.
            if (links.channel[cell] != 0)
            {
                failures += channels.channel[cell] == 0 || channels.direction[cell] != links.direction[cell] ||
                                    channels.bed_m[cell] != links.bed_m[cell] || channels.globalcell[cell] != links.globalcell[cell]
                                ? 1
                                : 0;
                continue;
            }
            if (channels.channel[cell] == 0)
                continue;
            tributarycells++;
            EXPECT_EQ(channels.globalcell[cell], -1);
            EXPECT_EQ(state.sea[cell], 0);
            EXPECT_EQ(order[cell], 1);
            int nx = 0, ny = 0;
            ASSERT_TRUE(down(x, y, nx, ny)) << "a tributary cell without a downstream cell at " << x << ", " << y;
            const std::size_t next = index(nx, ny);
            ASSERT_NE(channels.channel[next], 0);
            EXPECT_GE(channels.bed_m[cell], channels.bed_m[next]);
            EXPECT_LE(channels.bed_m[cell], state.surface_m[cell]) << "a bed above the start surface";
            EXPECT_LE(channels.discharge_m3s[cell], channels.discharge_m3s[next]);
            EXPECT_GE(channels.discharge_m3s[cell], 0.0f);
            // Down to the junction.
            int jx = x, jy = y, steps = 0;
            while (links.channel[index(jx, jy)] == 0 && steps < 8)
            {
                int ax = 0, ay = 0;
                ASSERT_TRUE(down(jx, jy, ax, ay));
                jx = ax;
                jy = ay;
                steps++;
            }
            EXPECT_LT(steps, 4) << "a tributary longer than a parent cell at " << x << ", " << y;
            EXPECT_EQ(state.sea[index(jx, jy)], 0);
            EXPECT_LE(std::abs(x / 4 - jx / 4), 1) << "outside the junction's parent cell and its 8 neighbours";
            EXPECT_LE(std::abs(y / 4 - jy / 4), 1) << "outside the junction's parent cell and its 8 neighbours";
            // No touching: every channel 8-neighbour is the downstream cell, an upstream cell of the
            // same tributary, or (next to the junction) within one cell of the junction.
            for (int oy = -1; oy <= 1; oy++)
                for (int ox = -1; ox <= 1; ox++)
                {
                    const int ax = x + ox, ay = y + oy;
                    if ((ox == 0 && oy == 0) || ax < 0 || ay < 0 || ax >= width || ay >= height)
                        continue;
                    const std::size_t other = index(ax, ay);
                    if (channels.channel[other] == 0 || other == next)
                        continue;
                    int bx = 0, by = 0;
                    const bool upstream = links.channel[other] == 0 && down(ax, ay, bx, by) && index(bx, by) == cell;
                    const bool nearjunction = links.channel[next] != 0 && std::abs(ax - nx) <= 1 && std::abs(ay - ny) <= 1;
                    if (!upstream && !nearjunction)
                    {
                        failures++;
                        if (failures < 5)
                            ADD_FAILURE() << "a tributary cell at " << x << ", " << y << " touches a channel at " << ax << ", " << ay;
                    }
                }
            // No crossing: a diagonal step's other diagonal is not a channel step.
            if (nx != x && ny != y)
            {
                const std::size_t c1 = index(nx, y), c2 = index(x, ny);
                if (channels.channel[c1] != 0 && channels.channel[c2] != 0)
                {
                    int ax = 0, ay = 0, bx = 0, by = 0;
                    const bool step = (down(nx, y, ax, ay) && index(ax, ay) == c2) || (down(x, ny, bx, by) && index(bx, by) == c1);
                    EXPECT_FALSE(step) << "a tributary crosses a channel at " << x << ", " << y;
                }
            }
        }
    EXPECT_EQ(failures, 0);
    EXPECT_EQ(tributarycells, report.growncells);
}

TEST_F(Landscape, RegionalTributariesInheritedAtLevelTwo)
{
    // Level 2 through the hook, its parent level 1's channels and tributaries handed over by hand
    // (ParentFields on a level-1 window holding level 2's domain and 3 global cells more): level 2
    // redraws level 1's tributaries and grows its own; no tributary bed stands above the start
    // surface; every tributary cell drains, without a loop, into a link, the sea or a counted end;
    // the link it reaches lies within two global cells (32 level-2 cells: the level-1 tree's
    // global cell and the redrawing's parent cells); the Strahler order rises above 1.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = tributarieswithvalleys();
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const RegionalWindow parentwindow = landscape::regionalwindow(256, 4, loophalo + 3, 170.0, -5.0, 180.0, 5.0);
    const landscape::FineWarp parentwarp = landscape::regionalwarp(parentwindow, inputs.seed, parameters);
    const landscape::LevelContext parentcontext{inputs, parentwindow, parameters, parentwarp, incision, details, 1, 4};
    landscape::LevelState parentstate;
    runtochannels(parentcontext, parentstate);
    landscape::growtributaries(parentcontext, parentstate);
    landscape::ParentFields parent;
    parent.window = handlevelwindow(256, 4, 1, 0, static_cast<std::int64_t>(parentwindow.domaincolumn0) * 4,
                                    parentwindow.domainfinewidth(), static_cast<std::int64_t>(parentwindow.domainrow0) * 4,
                                    parentwindow.domainfineheight());
    parent.channels = parentstate.channels;
    parent.sea = parentstate.sea;

    const RegionalWindow window = windowpsmall(16);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::LevelContext start{inputs, window, parameters, warp, incision, details};
    landscape::LevelState state;
    runtochannels(start, state);
    const landscape::FineChannels links = state.channels;
    const landscape::LevelContext context{inputs, window, parameters, warp, incision, details, 2, 4,
                                          nullptr, nullptr, nullptr, nullptr, &parent};
    const landscape::TributaryReport report = landscape::growtributaries(context, state);
    std::printf("[tributaries] level 2 of P': %d parent tributary cells (%d not drawn), %d inherited cells, %d ends; "
                "%d grown (%d cells); max order %d; %d beds above the start surface\n",
                report.inheritedparents, report.deadparents, report.inheritedcells, report.inheritedends, report.grown,
                report.growncells, report.maxorder, report.bedsabove);
    ASSERT_GT(report.inheritedcells, 0);
    ASSERT_GT(report.grown, 0);
    EXPECT_GE(report.maxorder, 2);
    EXPECT_EQ(report.bedsabove, 0);

    const landscape::FineChannels& channels = state.channels;
    const int width = window.domainfinewidth(), height = window.domainfineheight();
    const auto index = [&](int x, int y) { return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x); };
    int farthest = 0, ends = 0, failures = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = index(x, y);
            if (channels.channel[cell] == 0 || links.channel[cell] != 0)
                continue;
            EXPECT_EQ(channels.globalcell[cell], -1);
            EXPECT_LE(channels.bed_m[cell], state.surface_m[cell]) << "a bed above the start surface at " << x << ", " << y;
            int cx = x, cy = y, steps = 0;
            bool reached = false;
            while (steps < 400)
            {
                const std::size_t current = index(cx, cy);
                if (links.channel[current] != 0 || state.sea[current] != 0)
                {
                    reached = links.channel[current] != 0;
                    break;
                }
                const int direction = channels.direction[current];
                if (direction == 0)
                {
                    ends += steps == 0 ? 1 : 0;
                    break;
                }
                const int nx = cx + landscape::directionoffsetx(direction), ny = cy + landscape::directionoffsety(direction);
                if (nx < 0 || ny < 0 || nx >= width || ny >= height)
                    break;
                if (channels.channel[index(nx, ny)] == 0 && state.sea[index(nx, ny)] == 0)
                {
                    failures++;
                    break;
                }
                cx = nx;
                cy = ny;
                steps++;
            }
            EXPECT_LT(steps, 400) << "a loop through " << x << ", " << y;
            if (reached)
                farthest = std::max(farthest, std::max(std::abs(cx - x), std::abs(cy - y)));
        }
    std::printf("[tributaries] level 2 of P': farthest link reached %d level-2 cells away, %d ends\n", farthest, ends);
    EXPECT_EQ(failures, 0) << "tributary cells draining onto dry land";
    EXPECT_EQ(ends, report.inheritedends);
    EXPECT_LE(farthest, 32);
}

TEST_F(Landscape, RegionalTributariesOverlap)
{
    // The loop at two levels (k 16) with tributaries on (and R9's incision and valleys) at its
    // smallest halo: P and Q, which differ in column and row origin (Q across the date line), are
    // bit-identical on their overlap at every level, which holds more channel cells than without
    // tributaries. Level 2 inherits level 1's tributaries and grows its own, so its drainage
    // density rises over the one without tributaries by more than level 1's does.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = tributarieswithvalleys();
    RegionalParameters off = parameters;
    off.tributaries.enabled = 0;
    const RefinementResult p = refine(inputs, windowp(16), 2, parameters, "tributaries P");
    const RefinementResult q = refine(inputs, windowq(16), 2, parameters, "tributaries Q");
    const int channels = compareoverlaplevels(p, q, "tributaries P/Q");
    const RefinementResult plain = refine(inputs, windowp(16), 2, off, "no tributaries P");
    const RefinementResult plainq = refine(inputs, windowq(16), 2, off, "no tributaries Q");
    EXPECT_GT(channels, sharedchannels(plain, plainq));
    ASSERT_EQ(p.levels.size(), 2u);
    double previousratio = 1.0;
    for (std::size_t l = 0; l < p.levels.size(); l++)
    {
        EXPECT_GT(p.levels[l].tributaries, 0) << "level " << l + 1;
        EXPECT_GT(p.levels[l].tributarycells, 0) << "level " << l + 1;
        EXPECT_EQ(p.levels[l].tributarybedsabove, 0) << "level " << l + 1;
        EXPECT_EQ(plain.levels[l].tributarycells, 0) << "level " << l + 1;
        const double ratio = p.levels[l].drainagedensity_per_km / plain.levels[l].drainagedensity_per_km;
        std::printf("[tributaries] level %zu: drainage density %.5f against %.5f per km without (x %.2f)\n", l + 1,
                    p.levels[l].drainagedensity_per_km, plain.levels[l].drainagedensity_per_km, ratio);
        EXPECT_GT(ratio, previousratio) << "level " << l + 1;
        previousratio = ratio;
        EXPECT_EQ(p.levels[l].projection.forcedparents, 0) << "level " << l + 1;
    }
}

TEST_F(Landscape, RegionalTributariesOverlapAtThreeLevels)
{
    // The loop at three levels (k 64) with tributaries on and off at its smallest halo with them on:
    // P' and Q' (Q' across the date line) bit-identical on their overlap at every level, with more
    // shared channel cells than without tributaries.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = tributarieswithvalleys();
    RegionalParameters off = parameters;
    off.tributaries.enabled = 0;
    const RefinementResult p = refine(inputs, windowpsmall(64), 3, parameters, "tributaries P'");
    const RefinementResult q = refine(inputs, windowqsmall(64), 3, parameters, "tributaries Q'");
    const int channels = compareoverlaplevels(p, q, "tributaries P'/Q' at k 64");
    const RefinementResult plainp = refine(inputs, windowpsmall(64), 3, off, "no tributaries P'");
    const RefinementResult plainq = refine(inputs, windowqsmall(64), 3, off, "no tributaries Q'");
    // S3.8's regression test: bit-identical with tributaries off too. The unbounded solve (catchment
    // reach 0) takes the margin of 2 (unboundedsolvemarginglobalcells); without it level 3 has 3 global
    // cells of halo here and 263 elevations and 82 discharges differ at P's north edge: the level-3
    // solve reaches farther (S3.9's bounded catchments: RegionalTributariesOverlapAtThreeLevelsBoundedCatchments).
    const int without = compareoverlaplevels(plainp, plainq, "no tributaries P'/Q' at k 64");
    std::printf("[tributaries] P'/Q' at k 64: %d shared channel cells, %d without tributaries\n", channels, without);
    EXPECT_GT(channels, without);
    for (const LevelDiagnostics& level : p.levels)
    {
        EXPECT_GT(level.tributarycells, 0) << "level " << level.level;
        EXPECT_EQ(level.tributarybedsabove, 0) << "level " << level.level;
    }
}

TEST_F(Landscape, RegionalTributariesOverlapAtThreeLevelsBoundedCatchments)
{
    // S3.9 (multiscale plan §3, "Seamless solve: a bounded catchment", verify (b)): the pair above
    // with bounded catchments (catchment reach 1) at their solve margin, the bound's 2 (level 3 keeps 5
    // global cells of halo here): P' and Q' bit-identical on their overlap at every level with
    // tributaries off and on, and no kept cell drains to a domain's ring (a guard: 0 whenever the halo
    // holds 2 global cells).
    const RegionalGlobalInputs inputs = wideworld();
    RegionalParameters parameters = tributarieswithvalleys();
    parameters.relief.catchmentreach = 1;
    RegionalParameters off = parameters;
    off.tributaries.enabled = 0;
    const int margin = landscape::solvemargin(parameters.relief.catchmentreach);
    const RefinementResult plainp = refine(inputs, windowpsmall(64), 3, off, "bounded, no tributaries P'");
    const RefinementResult plainq = refine(inputs, windowqsmall(64), 3, off, "bounded, no tributaries Q'");
    const int without = compareoverlaplevels(plainp, plainq, "bounded, no tributaries P'/Q' at k 64");
    const RefinementResult p = refine(inputs, windowpsmall(64), 3, parameters, "bounded, tributaries P'");
    const RefinementResult q = refine(inputs, windowqsmall(64), 3, parameters, "bounded, tributaries Q'");
    const int channels = compareoverlaplevels(p, q, "bounded, tributaries P'/Q' at k 64");
    EXPECT_GT(channels, without);
    for (const RefinementResult* result : {&plainp, &plainq, &p, &q})
        for (const LevelDiagnostics& level : result->levels)
        {
            std::printf("[tributaries]   bounded level %d: halo %d parent cells, ring-draining %d, unreached %d, "
                        "sink junctions %d\n",
                        level.level, level.window.halo, level.ringbasincells, level.unreachedcells, level.sinkjunctions);
            EXPECT_EQ(level.window.halo, landscape::levelhalo(loophalo, 4, level.level, margin)) << "level " << level.level;
            EXPECT_EQ(level.ringbasincells, 0) << "level " << level.level;
        }
    std::printf("[tributaries] bounded: margin %d, level-3 halo %d parent cells\n", margin, plainp.levels[2].window.halo);
}

TEST_F(Landscape, RegionalTributariesOverlapSinglePass)
{
    // The single pass (k 16, the parent cell the global cell) with tributaries on at its smallest
    // halo, regionalhalo(2, true): P and Q bit-identical on their overlap, with more shared
    // channel cells than without tributaries; one global cell less is refused.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = tributarieswithvalleys();
    RegionalParameters off = parameters;
    off.tributaries.enabled = 0;
    constexpr int halo = landscape::regionalhalo(2, true);
    const OverlapCounts on = compareoverlap(inputs, windowp(16, halo), windowq(16, halo), "tributaries, single pass", parameters);
    const OverlapCounts without =
        compareoverlap(inputs, windowp(16, halo), windowq(16, halo), "no tributaries, single pass", off);
    EXPECT_GT(on.channel, without.channel);
    EXPECT_THROW(landscape::runregionalpass(inputs, windowp(16, halo - 1), parameters), std::invalid_argument);
    EXPECT_NO_THROW(landscape::runregionalpass(inputs, windowp(4, landscape::regionalhalo(2, false)), off));
}

TEST_F(Landscape, RegionalTributariesOutletDischarge)
{
    // RefinementLevelsOutletDischarge with tributaries on, its exception stated: a channel cell
    // whose direction leaves its global cell carries that cell's discharge (a link's outlet), or,
    // where a tributary crosses the global cell's edge, a tributary's discharge, which is at most
    // the largest discharge among the global cell and its 8 neighbours (a tributary's is capped at
    // its junction's, which lies there).
    const RegionalGlobalInputs inputs = wideworld();
    const RefinementResult result = refine(inputs, windowp(16), 2, tributarieswithvalleys(), "outlets P, tributaries");
    const RegionalFields& fields = result.fields;
    const RegionalWindow& window = fields.window;
    const int width = window.finewidth(), height = window.fineheight();
    const auto global = [&](int x, int y) {
        return at((window.column0 + x / window.scale) % inputs.width, window.row0 + y / window.scale, inputs.width);
    };
    int linkoutlets = 0, tributaryoutlets = 0, exceeding = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            const std::size_t cell = at(x, y, width);
            const int direction = fields.riverdirection[cell];
            if (direction == 0)
                continue;
            const int nx = x + landscape::directionoffsetx(direction), ny = y + landscape::directionoffsety(direction);
            if (nx < 0 || ny < 0 || nx >= width || ny >= height)
                continue;
            const std::size_t own = global(x, y);
            if (global(nx, ny) == own)
                continue;
            if (fields.discharge_m3s[cell] == inputs.network.discharge_m3s[own])
            {
                linkoutlets++;
                continue;
            }
            tributaryoutlets++;
            const int i = static_cast<int>(own % static_cast<std::size_t>(inputs.width));
            const int j = static_cast<int>(own / static_cast<std::size_t>(inputs.width));
            float largest = 0.0f;
            for (int oj = -1; oj <= 1; oj++)
                for (int oi = -1; oi <= 1; oi++)
                    if (j + oj >= 0 && j + oj < inputs.height)
                        largest = std::max(largest, inputs.network.discharge_m3s[at((i + oi + inputs.width) % inputs.width, j + oj, inputs.width)]);
            exceeding += fields.discharge_m3s[cell] > largest ? 1 : 0;
        }
    std::printf("[tributaries] outlets P: %d of links, %d of tributaries, %d above the neighbourhood's largest discharge\n",
                linkoutlets, tributaryoutlets, exceeding);
    EXPECT_GT(linkoutlets, 0);
    EXPECT_GT(tributaryoutlets, 0);
    EXPECT_EQ(exceeding, 0);
}

TEST_F(Landscape, RegionalTributariesOffIgnoresTheirKeys)
{
    // enabled 0: the other tributary keys change nothing (the S2 loop exactly).
    const RegionalGlobalInputs inputs = wideworld();
    RegionalParameters parameters = tributarieswithvalleys();
    parameters.tributaries.enabled = 0;
    RegionalParameters changed = parameters;
    changed.tributaries.spacing_km = 50.0;
    changed.tributaries.hackexponent = 0.7;
    changed.tributaries.hackcoefficient = 2.0;
    changed.tributaries.junctionangle_deg = 45.0;
    changed.tributaries.anglespread = 0.5;
    changed.tributaries.nodemargin = 0.0;
    changed.tributaries.lengthminfraction = 0.25;
    const int halo = landscape::refinementhalo(2, 0);
    const RefinementResult a = refine(inputs, windowp(16, halo), 2, parameters, "off P");
    const RefinementResult b = refine(inputs, windowp(16, halo), 2, changed, "off, keys changed, P");
    EXPECT_EQ(a.fields.elevation_m, b.fields.elevation_m);
    EXPECT_EQ(a.fields.riverdirection, b.fields.riverdirection);
    EXPECT_EQ(a.fields.discharge_m3s, b.fields.discharge_m3s);
    EXPECT_EQ(a.coarserelevation_m, b.coarserelevation_m);
    for (const LevelDiagnostics& level : a.levels)
    {
        EXPECT_EQ(level.tributaries, 0);
        EXPECT_EQ(level.tributarycells, 0);
    }
}
