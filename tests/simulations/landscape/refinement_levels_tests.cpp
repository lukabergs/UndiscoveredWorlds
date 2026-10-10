// Tests of landscape/refinement_levels.hpp (docs/architecture/multiscale-plan.md, section
// 3, "Tests"), suite Landscape: Landscape.RefinementLevels*. S0.4 tests the cell sizes and
// the level assignment; S2.3 adds the loop's tests on the synthetic world of
// regional_test_world.hpp at factor 4: the overlap identity at every level with every
// feature on (the Phase R rule: on a world wider than the widened boxes), the parent-mean
// error against its bound, the outlet discharge, the level consistency (k = 16
// box-averaged equals k = 4) and the parity of one level with the single pass; S2.3's
// rework (the bounds, the detail band, the global lake depth) the valleys at level 2.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/refinement_levels.hpp"

#include <stdexcept>
#include <string>

namespace
{
using landscape::LevelDiagnostics;
using landscape::LevelTarget;
using landscape::LevelWindow;
using landscape::RefinementParameters;
using landscape::RefinementResult;

// R9's inputs (regional_pass_tests.cpp's withvalleyfields): a rock uplift rate smooth in
// latitude (35 m/Myr at the equator, -5 at the poles) and a denudation smooth in longitude
// (100 to 500 m over a 20 Myr erosion run).
RegionalGlobalInputs refinementvalleyfields(RegionalGlobalInputs inputs)
{
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
    return inputs;
}

// The synthetic world on 256 columns (global cells of 156 km; levels 1 and 2 of 39.1 and
// 9.77 km) with R9's and E5's fields.
RegionalGlobalInputs wideworld() { return withglacialfields(refinementvalleyfields(syntheticworld(256))); }

// R9's incision and valleys: floors of 100 km x sqrt(Q / 1e3 m3/s), 2.56 level-1 cells at the
// reference (in the band [1, 4) at level 1 for Q >= 153 m3/s, at level 2 from 9.6 to 153
// m3/s), at most 150 km (3.84 level-1 cells: never beyond the band at level 1, so the single
// pass's clamp to k never acts).
RegionalParameters valleysatlevelone()
{
    RegionalParameters parameters;
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

// The valleys at level 2: floors of 20 km x sqrt(Q / 1e3 m3/s), 2.05 level-2 cells at the
// reference (in the band [1, 4) for Q from 240 to 3,800 m3/s; level 1 carves the larger,
// whose half-width reaches one level-1 cell, at most 40 km = 1.02 cells).
RegionalParameters valleysatleveltwo()
{
    RegionalParameters parameters = valleysatlevelone();
    parameters.valleys.floor_km = 20.0;
    parameters.valleys.floormax_km = 40.0;
    parameters.valleys.reachcells = 3;
    return parameters;
}

// Every feature on: the coast terms of Phases R and R2 (warp, coast noise, land share,
// one-sided noise, the island sink, the glacial and denudation terms, the smooth gate, the
// straight-run breaker), the lake shores, R9's incision and valleys (both levels), the
// troughs (30 km) and cirques (50 km: 1.28 level-1 and 5.1 level-2 cells, so level 1 carves
// them) and the ice above the ELA. Six share and six sign-keeping passes keep the widened
// land-share box clear of the poles.
RegionalParameters everyfeature()
{
    RegionalParameters parameters = valleysatlevelone();
    parameters.warp.amplitude_cells = 0.5;
    parameters.warp.wavelength_cells = 2.0;
    parameters.warp.fade_m = 400.0;
    parameters.coast.noiseamplitude_m = 150.0;
    parameters.coast.noisepersistence = 0.5;
    parameters.coast.noiserelieffraction = 0.25;
    parameters.coast.seadrop_m = 300.0;
    parameters.coast.landshare = 0.5;
    parameters.coast.onesided = 1;
    parameters.coast.islandmargin = 0.5;
    parameters.coast.glacialamplitude_m = 100.0;
    parameters.coast.glacialdropfraction = 0.9;
    parameters.coast.denudationamplitude_m = 50.0;
    parameters.coast.denudationreference_m = 1000.0;
    parameters.coast.sharepasses = 6;
    parameters.coast.signedpasses = 6;
    parameters.edges.smoothgate = 3;
    parameters.edges.straightmax = 4;
    parameters.shores.shorenoise = 0.15;
    parameters.shores.terrain_m = 50.0;
    parameters.troughs.width_km = 30.0;
    parameters.cirques.radius_km = 50.0;
    parameters.ice.thicknessfraction = 0.3;
    return parameters;
}

// On the 256-column world: the continent's centre, P (160 E .. 180, 15 S .. 10 N) and Q
// (170 E .. 170 W across the date line, 10 S .. 15 N), sharing 170 E .. 180 by 10 S .. 10 N
// in whole global cells, with the lake (global cell 248, 71), ice and channels; the coast,
// coastwest (120 .. 140 E, 5 S .. 15 N) and coasteast (130 .. 150 E, 0 .. 20 N), the windows of
// RegionalOverlapIslandsWideWorld.
RegionalWindow windowp() { return landscape::regionalwindow(256, 16, testhalo, 160.0, -15.0, 180.0, 10.0); }
RegionalWindow windowq() { return landscape::regionalwindow(256, 16, testhalo, 170.0, -10.0, -170.0, 15.0); }
RegionalWindow coastwest() { return landscape::regionalwindow(256, 16, testhalo, 120.0, -5.0, 140.0, 15.0); }
RegionalWindow coasteast() { return landscape::regionalwindow(256, 16, testhalo, 130.0, 0.0, 150.0, 20.0); }

RefinementResult refine(const RegionalGlobalInputs& inputs, const RegionalWindow& window, int level,
                        const RegionalParameters& parameters, const std::string& name)
{
    LevelTarget target;
    target.window = window;
    target.level = level;
    const auto start = std::chrono::steady_clock::now();
    RefinementResult result = landscape::runrefinement(inputs, target, parameters, RefinementParameters{});
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("[refinement] %s: L %d, kept %d x %d global, %.3f s\n", name.c_str(), level, window.columns, window.rows,
                seconds);
    for (const LevelDiagnostics& diagnostics : result.levels)
    {
        const landscape::ProjectionReport& projection = diagnostics.projection;
        std::printf("[refinement]   level %d (%.2f km): domain %d x %d, %.3f s; block-mean %.4g m, parent-mean %.3g m "
                    "(bound %.3g m), global-mean %.3g m; lambda max %.1f m, border step max %.1f m p99 %.1f m; beds moved "
                    "%d, forced %d (%.3g m); dried %d, shores below lakes %d, land below the sea level %d\n",
                    diagnostics.level, diagnostics.cellsize_km, diagnostics.window.domaincolumns,
                    diagnostics.window.domainrows, diagnostics.seconds, diagnostics.blockmeanerror_m,
                    projection.maxparentmeanerror_m, diagnostics.meanbound_m, diagnostics.globalmeanerror_m,
                    projection.maxlambda_m, projection.maxborderstep_m, projection.borderstep99_m,
                    projection.bedsmovedparents, projection.forcedparents, projection.maxforcedmove_m,
                    diagnostics.driedcells, diagnostics.shoresbelowlake, diagnostics.landbelowsea);
    }
    return result;
}

// The cells two level windows' kept boxes share, as index pairs into their kept fields.
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

// Every level's parent-mean error within its bound (half a float ulp of the kept box's
// largest |F|, 2.4e-4 m on the synthetic world) and the accumulated global-mean error within
// the sum of the bounds; no parent forced past its headrooms.
void expectexact(const RefinementResult& result, const std::string& name)
{
    double bounds = 0.0;
    for (const LevelDiagnostics& diagnostics : result.levels)
    {
        bounds += diagnostics.meanbound_m;
        EXPECT_GT(diagnostics.meanbound_m, 0.0) << name << ", level " << diagnostics.level;
        EXPECT_LE(diagnostics.meanbound_m, 2.5e-4) << name << ", level " << diagnostics.level;
        EXPECT_LE(diagnostics.projection.maxparentmeanerror_m, diagnostics.meanbound_m) << name << ", level " << diagnostics.level;
        EXPECT_LE(diagnostics.globalmeanerror_m, bounds) << name << ", level " << diagnostics.level;
        EXPECT_GT(diagnostics.projection.parents, 0) << name << ", level " << diagnostics.level;
        EXPECT_EQ(diagnostics.projection.forcedparents, 0) << name << ", level " << diagnostics.level;
    }
}

// Two refinements compared on their overlap at every level: the coarser levels' projected
// elevation and every field of the final level bit-identical. Returns the final level's
// shared cells and how many are sea, lake, ice and channel.
OverlapCounts compareoverlaplevels(const RefinementResult& a, const RefinementResult& b, const std::string& label)
{
    const int levels = static_cast<int>(a.levels.size());
    for (int level = 1; level < levels; level++)
    {
        const auto cells = sharedcells(a.levels[static_cast<std::size_t>(level - 1)].window,
                                       b.levels[static_cast<std::size_t>(level - 1)].window);
        const int count = differing(a.coarserelevation_m[static_cast<std::size_t>(level - 1)],
                                    b.coarserelevation_m[static_cast<std::size_t>(level - 1)], cells);
        std::printf("[refinement] %s, level %d: %zu shared cells, %d elevations differ\n", label.c_str(), level,
                    cells.size(), count);
        EXPECT_GT(cells.size(), 0u) << label << ", level " << level;
        EXPECT_EQ(count, 0) << label << ": level " << level << "'s elevation";
    }
    const auto cells = sharedcells(a.levels.back().window, b.levels.back().window);
    const RegionalFields& fa = a.fields;
    const RegionalFields& fb = b.fields;
    const std::vector<std::pair<std::string, int>> counts = {
        {"elevation_m", differing(fa.elevation_m, fb.elevation_m, cells)},
        {"sea", differing(fa.sea, fb.sea, cells)},
        {"riverdirection", differing(fa.riverdirection, fb.riverdirection, cells)},
        {"discharge_m3s", differing(fa.discharge_m3s, fb.discharge_m3s, cells)},
        {"januaryflow", differing(fa.januaryflow, fb.januaryflow, cells)},
        {"julyflow", differing(fa.julyflow, fb.julyflow, cells)},
        {"lakesurface_m", differing(fa.lakesurface_m, fb.lakesurface_m, cells)},
        {"januarytemperature_c", differing(fa.januarytemperature_c, fb.januarytemperature_c, cells)},
        {"julytemperature_c", differing(fa.julytemperature_c, fb.julytemperature_c, cells)},
        {"januaryprecipitation_mm", differing(fa.januaryprecipitation_mm, fb.januaryprecipitation_mm, cells)},
        {"julyprecipitation_mm", differing(fa.julyprecipitation_mm, fb.julyprecipitation_mm, cells)},
        {"icethickness_m", differing(fa.icethickness_m, fb.icethickness_m, cells)},
        {"climate", differing(fa.climate, fb.climate, cells)}};
    OverlapCounts overlap;
    overlap.shared = static_cast<int>(cells.size());
    for (const auto& [ca, cb] : cells)
    {
        overlap.sea += fa.sea[ca] != 0 ? 1 : 0;
        overlap.lake += fa.lakesurface_m[ca] > 0.0f ? 1 : 0;
        overlap.ice += fa.icethickness_m[ca] > 0.0f ? 1 : 0;
        overlap.channel += fa.riverdirection[ca] != 0 ? 1 : 0;
    }
    std::printf("[refinement] %s, final level: %d shared cells: %d sea, %d lake, %d ice, %d channel\n", label.c_str(),
                overlap.shared, overlap.sea, overlap.lake, overlap.ice, overlap.channel);
    for (const auto& [name, count] : counts)
    {
        if (count > 0)
            std::printf("[refinement] %s, final level: %s differs on %d cells\n", label.c_str(), name.c_str(), count);
        EXPECT_EQ(count, 0) << label << ": " << name;
    }
    return overlap;
}

// The global cell of a final-level kept cell.
std::size_t globalcell(const RegionalWindow& window, int width, int x, int y)
{
    return at((window.column0 + x / window.scale) % width, window.row0 + y / window.scale, width);
}
}

// A 256-column world (128 rows): level 0 cells of pi R / 128 = 156.4 km; at factor 4 the
// levels 1 to 4 are 39.1, 9.77, 2.44 and 0.611 km.
TEST_F(Landscape, RefinementLevelsCellSizesAndAssignment)
{
    const double global_km = pi * radius_m / 128.0 / 1000.0;
    EXPECT_NEAR(landscape::levelcellsize_km(radius_m, 256, 4, 0), global_km, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius_m, 256, 4, 1), global_km / 4.0, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius_m, 256, 4, 4), global_km / 256.0, 1e-12);
    EXPECT_NEAR(landscape::levelcellsize_km(radius_m, 256, 2, 3), global_km / 8.0, 1e-9);
    EXPECT_NEAR(landscape::levelcellsize_km(radius_m, 256, 4, 2), 9.7727, 1e-3);

    // The coarsest level whose cells are at most size / detail_cells, within 1 .. levels;
    // detail_cells 2 (the default): a detail is resolved by two cells of its level.
    EXPECT_EQ(landscape::detaillevel(10.0, 2.0, radius_m, 256, 4, 4), 3);    // 2.44 <= 5 < 9.77
    EXPECT_EQ(landscape::detaillevel(80.0, 2.0, radius_m, 256, 4, 4), 1);    // 39.1 <= 40
    EXPECT_EQ(landscape::detaillevel(78.0, 2.0, radius_m, 256, 4, 4), 2);    // 39.1 > 39
    EXPECT_EQ(landscape::detaillevel(5.0, 1.0, radius_m, 256, 4, 4), 3);     // one cell per detail: 2.44 <= 5
    EXPECT_EQ(landscape::detaillevel(1000.0, 2.0, radius_m, 256, 4, 4), 1);  // coarser than a global cell: level 1
    EXPECT_EQ(landscape::detaillevel(0.1, 2.0, radius_m, 256, 4, 4), 4);     // finer than the finest level: clamped
    EXPECT_EQ(landscape::detaillevel(0.1, 2.0, radius_m, 256, 4, 2), 2);
    EXPECT_EQ(landscape::detaillevel(10.0, 2.0, radius_m, 256, 2, 8), 5);    // factor 2: 78.2, 39.1, 19.5, 9.77, 4.89
    EXPECT_THROW(landscape::detaillevel(0.0, 2.0, radius_m, 256, 4, 4), std::invalid_argument);
    EXPECT_THROW(landscape::detaillevel(5.0, 0.0, radius_m, 256, 4, 4), std::invalid_argument);
    EXPECT_THROW(landscape::detaillevel(5.0, 2.0, radius_m, 256, 4, 0), std::invalid_argument);
    EXPECT_THROW(landscape::levelcellsize_km(radius_m, 255, 4, 1), std::invalid_argument);
    EXPECT_THROW(landscape::levelcellsize_km(radius_m, 256, 1, 1), std::invalid_argument);

    // The steps' fixed order within a level (section 3, "S2.3's level placement", item 5).
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::noiseoctave), 0);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::presmoothing), 1);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::incision), 3);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::solve), 5);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::lakes), 9);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::projection), 10);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::ice), landscape::levelstepcount - 1);
    EXPECT_EQ(std::string(landscape::levelstepname(landscape::LevelStep::glacialforms)), "glacial_forms");
}

TEST_F(Landscape, RefinementLevelsOverlap)
{
    // Two levels at factor 4 (k = 16) with every feature on, the detail units at both
    // levels: P and Q (Q across the date line; they differ in column and row origin) hold
    // the lake, ice and channels on their overlap; R and S the coast. Every level's shared
    // kept cells are bit-identical, and every level is exact.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    {
        const RefinementResult p = refine(inputs, windowp(), 2, parameters, "every feature P");
        const RefinementResult q = refine(inputs, windowq(), 2, parameters, "every feature Q");
        ASSERT_EQ(p.levels.size(), 2u);
        ASSERT_EQ(p.coarserelevation_m.size(), 1u);
        expectexact(p, "P");
        expectexact(q, "Q");
        const OverlapCounts counts = compareoverlaplevels(p, q, "every feature P/Q");
        EXPECT_GT(counts.lake, 0);
        EXPECT_GT(counts.ice, 0);
        EXPECT_GT(counts.channel, 0);
        // The detail units run at every level (each form where its size is in the band).
        const auto ran = [](const LevelDiagnostics& level, landscape::LevelStep step) {
            return std::find(level.steps.begin(), level.steps.end(), step) != level.steps.end();
        };
        for (const LevelDiagnostics& level : p.levels)
        {
            EXPECT_TRUE(ran(level, landscape::LevelStep::valleyfloors));
            EXPECT_TRUE(ran(level, landscape::LevelStep::glacialforms));
        }
        EXPECT_FALSE(ran(p.levels[0], landscape::LevelStep::ice));
        EXPECT_TRUE(ran(p.levels[1], landscape::LevelStep::ice));
        // The features act: the result differs from the defaults'.
        const RefinementResult plain = refine(inputs, windowp(), 2, RegionalParameters{}, "defaults P");
        int changed = 0;
        for (std::size_t cell = 0; cell < p.fields.elevation_m.size(); cell++)
            changed += p.fields.elevation_m[cell] != plain.fields.elevation_m[cell] ? 1 : 0;
        EXPECT_GT(changed, 0);
    }
    {
        const RefinementResult r = refine(inputs, coastwest(), 2, parameters, "every feature R");
        const RefinementResult s = refine(inputs, coasteast(), 2, parameters, "every feature S");
        expectexact(r, "R");
        expectexact(s, "S");
        const OverlapCounts counts = compareoverlaplevels(r, s, "every feature R/S");
        EXPECT_GT(counts.sea, 0);
        EXPECT_LT(counts.sea, counts.shared);
        EXPECT_GT(counts.channel, 0);
    }
}

TEST_F(Landscape, RefinementLevelsOverlapValleysAtLevelTwo)
{
    // The valleys carved at level 2 (half-widths of 1 to 4 level-2 cells, the band rule) on P
    // and Q, which differ in column and row origin (Q across the date line): every level's
    // shared kept cells bit-identical; the valleys act at level 2 (level 1 is the result
    // without them where no form reaches one level-1 cell).
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = valleysatleveltwo();
    ASSERT_NE(windowp().row0, windowq().row0);
    ASSERT_NE(windowp().column0, windowq().column0);
    const RefinementResult p = refine(inputs, windowp(), 2, parameters, "valleys at level 2, P");
    const RefinementResult q = refine(inputs, windowq(), 2, parameters, "valleys at level 2, Q");
    expectexact(p, "P");
    expectexact(q, "Q");
    const OverlapCounts counts = compareoverlaplevels(p, q, "valleys at level 2, P/Q");
    EXPECT_GT(counts.channel, 0);
    RegionalParameters without = parameters;
    without.valleys.floor_km = 0.0;
    const RefinementResult plain = refine(inputs, windowp(), 2, without, "no valleys, P");
    int levelone = 0, leveltwo = 0;
    for (std::size_t cell = 0; cell < p.coarserelevation_m[0].size(); cell++)
        levelone += p.coarserelevation_m[0][cell] != plain.coarserelevation_m[0][cell] ? 1 : 0;
    for (std::size_t cell = 0; cell < p.fields.elevation_m.size(); cell++)
        leveltwo += p.fields.elevation_m[cell] != plain.fields.elevation_m[cell] ? 1 : 0;
    std::printf("[refinement] valleys at level 2: %d level-1 cells and %d level-2 cells differ from the run without "
                "valleys\n", levelone, leveltwo);
    EXPECT_GT(leveltwo, 0);
}

TEST_F(Landscape, RefinementLevelsOutletDischarge)
{
    // At the final level, a channel cell whose direction leaves its global cell is that
    // cell's pinned outlet: its discharge is the global cell's discharge exactly.
    const RegionalGlobalInputs inputs = wideworld();
    const RefinementResult result = refine(inputs, windowp(), 2, valleysatlevelone(), "outlets P");
    const RegionalFields& fields = result.fields;
    const RegionalWindow& window = fields.window;
    ASSERT_EQ(window.scale, 16);
    const int width = window.finewidth(), height = window.fineheight();
    int outlets = 0, mismatches = 0;
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
            const std::size_t own = globalcell(window, inputs.width, x, y);
            if (globalcell(window, inputs.width, nx, ny) == own)
                continue;
            outlets++;
            mismatches += fields.discharge_m3s[cell] != inputs.network.discharge_m3s[own] ? 1 : 0;
        }
    std::printf("[refinement] outlets P: %d outlet cells, %d discharges differ from the global cell's\n", outlets, mismatches);
    EXPECT_GT(outlets, 0);
    EXPECT_EQ(mismatches, 0);
}

TEST_F(Landscape, RefinementLevelsConsistency)
{
    // With every feature on, the k = 16 result of the two-level pass averaged over 4 x 4
    // blocks equals its own level 1 and the k = 4 (one-level) result of the same target within
    // the level-2 bound; the level 1 of the two-level pass is the one-level result bit for bit
    // (no level depends on L: the band places every form by its own size in level cells).
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    const RefinementResult one = refine(inputs, windowp(), 1, parameters, "consistency P, one level");
    const RefinementResult two = refine(inputs, windowp(), 2, parameters, "consistency P, two levels");
    ASSERT_EQ(two.coarserelevation_m.size(), 1u);
    expectexact(one, "one level");
    expectexact(two, "two levels");
    const LevelWindow& levelone = two.levels[0].window;
    const LevelWindow& target = one.levels[0].window;
    const auto cells = sharedcells(target, levelone);
    ASSERT_EQ(cells.size(), target.keptsize());
    EXPECT_EQ(differing(one.fields.elevation_m, two.coarserelevation_m[0], cells), 0);

    const int columns = target.columns, rows = target.rows;
    const int finewidth = two.fields.window.finewidth();
    ASSERT_EQ(finewidth, 4 * columns);
    double largestown = 0.0, largestone = 0.0;
    for (const auto& [index, inlevelone] : cells)
    {
        const int x = static_cast<int>(index % static_cast<std::size_t>(columns));
        const int y = static_cast<int>(index / static_cast<std::size_t>(columns));
        double sum = 0.0;
        for (int fy = 4 * y; fy < 4 * y + 4; fy++)
            for (int fx = 4 * x; fx < 4 * x + 4; fx++)
                sum += static_cast<double>(two.fields.elevation_m[at(fx, fy, finewidth)]);
        const double mean = sum / 16.0;
        largestown = std::max(largestown, std::fabs(mean - static_cast<double>(two.coarserelevation_m[0][inlevelone])));
        largestone = std::max(largestone, std::fabs(mean - static_cast<double>(one.fields.elevation_m[index])));
    }
    std::printf("[refinement] consistency P: k = 16 box-averaged against its level 1 %.3g m, against k = 4 %.3g m "
                "(%d x %d level-1 cells)\n", largestown, largestone, columns, rows);
    EXPECT_LE(largestown, two.levels[1].meanbound_m);
    EXPECT_LE(largestone, two.levels[1].meanbound_m);
}

TEST_F(Landscape, RefinementLevelsParityAtOneLevel)
{
    // One level (k = 4) against the single pass on the same window, every feature on, every
    // form in the band at level 1 (design item 4, with S2.3's rework): sea, river directions,
    // discharge, flows, lake surfaces, precipitation and climate bit-identical; the elevation
    // differs by the projection and the lakes' global depth, its block means exact; the
    // temperatures follow it by the lapse rate. The valley sides are on: a floor under one
    // cell with sides of 4 cells is in the band (its full width decides the lower end), as
    // the single pass carves it.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    for (const RegionalWindow& base : {windowp(), coastwest()})
    {
        RegionalWindow window = base;
        window.scale = 4;
        SCOPED_TRACE("window at column " + std::to_string(window.column0));
        const RegionalFields single = timedpass(inputs, window, "parity single pass", parameters);
        const RefinementResult refined = refine(inputs, window, 1, parameters, "parity one level");
        const RegionalFields& level = refined.fields;
        ASSERT_EQ(level.window.scale, 4);
        ASSERT_EQ(level.elevation_m.size(), single.elevation_m.size());
        expectexact(refined, "one level");
        std::vector<std::pair<std::size_t, std::size_t>> cells;
        for (std::size_t cell = 0; cell < single.elevation_m.size(); cell++)
            cells.push_back({cell, cell});
        EXPECT_EQ(differing(single.sea, level.sea, cells), 0);
        EXPECT_EQ(differing(single.riverdirection, level.riverdirection, cells), 0);
        EXPECT_EQ(differing(single.discharge_m3s, level.discharge_m3s, cells), 0);
        EXPECT_EQ(differing(single.januaryflow, level.januaryflow, cells), 0);
        EXPECT_EQ(differing(single.julyflow, level.julyflow, cells), 0);
        EXPECT_EQ(differing(single.lakesurface_m, level.lakesurface_m, cells), 0);
        EXPECT_EQ(differing(single.januaryprecipitation_mm, level.januaryprecipitation_mm, cells), 0);
        EXPECT_EQ(differing(single.julyprecipitation_mm, level.julyprecipitation_mm, cells), 0);
        EXPECT_EQ(differing(single.climate, level.climate, cells), 0);

        // The bounds: land channel beds stay (unless a parent moved its beds), sea cells stay
        // at or below the sea level, wet cells below their surface, dry land above the sea
        // level where it was; the rest is lifted or lowered; the temperatures follow the
        // elevation by the lapse rate.
        const landscape::ProjectionReport& projection = refined.levels[0].projection;
        const double lapse = static_cast<double>(inputs.lapserate_c_per_km) / 1000.0;
        int bedsmoved = 0, beds = 0, seaabove = 0, wetabove = 0, landsunk = 0, changed = 0;
        double temperature = 0.0;
        for (std::size_t cell = 0; cell < single.elevation_m.size(); cell++)
        {
            const float before = single.elevation_m[cell], after = level.elevation_m[cell];
            const bool sea = single.sea[cell] != 0, wet = single.lakesurface_m[cell] > 0.0f;
            if (single.riverdirection[cell] != 0 && !sea && !wet)
            {
                beds++;
                bedsmoved += before != after ? 1 : 0;
            }
            seaabove += sea && after > sealevel ? 1 : 0;
            wetabove += wet && !(after < level.lakesurface_m[cell]) && single.riverdirection[cell] == 0 ? 1 : 0;
            landsunk += !sea && !wet && before > sealevel && !(after > sealevel) ? 1 : 0;
            changed += before != after ? 1 : 0;
            const double sealevelheight = static_cast<double>(sealevel);
            const double lift = std::max(0.0, static_cast<double>(before) - sealevelheight) -
                                std::max(0.0, static_cast<double>(after) - sealevelheight);
            temperature = std::max(temperature, std::fabs(static_cast<double>(level.januarytemperature_c[cell]) -
                                                          static_cast<double>(single.januarytemperature_c[cell]) - lapse * lift));
        }
        std::printf("[refinement] parity: %d land channel beds (%d moved; beds moved in %d parents), %d cells changed; sea "
                    "above the sea level %d, wet at or above the surface %d, land sunk to the sea level %d; temperature "
                    "off the lapse rule by %.3g C\n",
                    beds, bedsmoved, projection.bedsmovedparents, changed, seaabove, wetabove, landsunk, temperature);
        EXPECT_GT(beds, 0);
        if (projection.bedsmovedparents == 0 && projection.forcedparents == 0)
        {
            EXPECT_EQ(bedsmoved, 0);
            EXPECT_EQ(landsunk, 0);
        }
        EXPECT_EQ(seaabove, 0);
        EXPECT_EQ(wetabove, 0);
        EXPECT_GT(changed, 0);
        EXPECT_LE(temperature, 1.0e-3);
    }
}

TEST_F(Landscape, RefinementLevelsRejects)
{
    const RegionalGlobalInputs inputs = syntheticworld(64);
    LevelTarget target;
    target.window = windowc(4);
    target.level = 0;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, RefinementParameters{}), std::invalid_argument);
    target.level = 1;
    RefinementParameters one;
    one.factor = 1;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, one), std::invalid_argument);
    target.window.halo = 5;  // below 2 x meanpasses + 2
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, RefinementParameters{}), std::invalid_argument);
}
