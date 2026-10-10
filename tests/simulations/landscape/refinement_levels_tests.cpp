// Tests of landscape/refinement_levels.hpp (docs/architecture/multiscale-plan.md, section
// 3, "Tests"), suite Landscape: Landscape.RefinementLevels*. S0.4 tests the cell sizes and
// the level assignment; S2.3 adds the loop's tests on the synthetic world of
// regional_test_world.hpp at factor 4: the overlap identity at every level with every
// feature on (the Phase R rule: on a world wider than the widened boxes), the parent-mean
// error against its bound, the outlet discharge, the level consistency (k = 16
// box-averaged equals k = 4) and the parity of one level with the single pass; S2.3's
// rework (the bounds, the detail band, the global lake depth) the valleys at level 2; S2.4
// the resample of the final level to a requested metres per pixel; S2.7 the overlap and the
// consistency with the smoothing before the projection on; S2.8 the inherited masks; S3.1 the
// limited prolongation; S3.2 the solve against its parent with the channel cells' means.

#include <gtest/gtest.h>

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
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
// The loop's smallest halo at meanpasses 2 without smoothing: testhalo plus the lake flood's reach
// (refinementhalo).
constexpr int loophalo = landscape::refinementhalo(2, 0);
RegionalWindow windowp() { return landscape::regionalwindow(256, 16, loophalo, 160.0, -15.0, 180.0, 10.0); }
RegionalWindow windowq() { return landscape::regionalwindow(256, 16, loophalo, 170.0, -10.0, -170.0, 15.0); }
RegionalWindow coastwest() { return landscape::regionalwindow(256, 16, loophalo, 120.0, -5.0, 140.0, 15.0); }
RegionalWindow coasteast() { return landscape::regionalwindow(256, 16, loophalo, 130.0, 0.0, 150.0, 20.0); }

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
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::smoothing), 10);
    EXPECT_EQ(static_cast<int>(landscape::LevelStep::projection), 11);
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
        // The lake surfaces differ from S2.8 on: the loop floods every cell below a lake's
        // surface connected to its core (RefinementLevelsInheritedMasks).
        EXPECT_EQ(differing(single.januaryprecipitation_mm, level.januaryprecipitation_mm, cells), 0);
        EXPECT_EQ(differing(single.julyprecipitation_mm, level.julyprecipitation_mm, cells), 0);
        EXPECT_EQ(differing(single.climate, level.climate, cells), 0);
        int wetsingle = 0, wetlevel = 0;
        for (std::size_t cell = 0; cell < single.lakesurface_m.size(); cell++)
        {
            wetsingle += single.lakesurface_m[cell] > 0.0f ? 1 : 0;
            wetlevel += level.lakesurface_m[cell] > 0.0f ? 1 : 0;
        }
        // The flood only adds water: every wet cell of the single pass is wet at level 1.
        int lost = 0;
        for (std::size_t cell = 0; cell < single.lakesurface_m.size(); cell++)
            lost += single.lakesurface_m[cell] > 0.0f && !(level.lakesurface_m[cell] > 0.0f) ? 1 : 0;
        std::printf("[refinement] parity: wet cells %d in the single pass, %d at level 1 (the flood), %d of the single "
                    "pass's dry at level 1\n", wetsingle, wetlevel, lost);
        EXPECT_EQ(lost, 0);

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
            const bool sea = single.sea[cell] != 0, wet = level.lakesurface_m[cell] > 0.0f;  // the level's lakes (S2.8)
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
    // The lake flood's reach (S2.8): 6 is below 2 x 2 + 2 + 2.
    target.window.halo = testhalo;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, RefinementParameters{}), std::invalid_argument);
    // The smoothing's passes widen the rule: 6 is below 2 x 2 + 2 + 2 x 1 + 2.
    target.window.halo = testhalo;
    RefinementParameters smoothing;
    smoothing.smoothpasses = 1;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, smoothing), std::invalid_argument);
    target.window.halo = testhalo + 2;
    smoothing.smoothpasses = -1;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, smoothing), std::invalid_argument);
    smoothing.smoothpasses = 1;
    smoothing.smoothrelaxation = 0.0;
    EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, smoothing), std::invalid_argument);
}

// Step S2.7: the smoothing before every projection (6 passes, relaxation 1.4) with every
// feature on, the halo widened by 2 x 6 global cells: the overlap identity at both levels
// (P and Q, Q across the date line), every level exact, level 1 of two levels equal to one
// level bit for bit, the smoothing's step recorded and lambda reduced.
TEST_F(Landscape, RefinementLevelsSmoothing)
{
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    RefinementParameters refinement;
    refinement.smoothpasses = 6;
    const int halo = landscape::refinementhalo(2, refinement.smoothpasses);
    const RegionalWindow p = landscape::regionalwindow(256, 16, halo, 160.0, -15.0, 180.0, 10.0);
    const RegionalWindow q = landscape::regionalwindow(256, 16, halo, 170.0, -10.0, -170.0, 15.0);
    const auto run = [&](const RegionalWindow& window, int level, const std::string& name) {
        LevelTarget target;
        target.window = window;
        target.level = level;
        const auto start = std::chrono::steady_clock::now();
        RefinementResult result = landscape::runrefinement(inputs, target, parameters, refinement);
        std::printf("[refinement] %s: L %d, %.3f s\n", name.c_str(), level,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        for (const LevelDiagnostics& diagnostics : result.levels)
            std::printf("[refinement]   level %d: lambda %.1f -> %.2f m (step %.1f -> %.2f m), largest move %.1f m; "
                        "border step max %.2f m p99 %.2f m; parent-mean %.3g m (bound %.3g m); beds moved %d, forced %d; "
                        "land below the sea level %d\n",
                        diagnostics.level, diagnostics.smoothing.maxlambdabefore_m, diagnostics.smoothing.maxlambdaafter_m,
                        diagnostics.smoothing.maxlambdastepbefore_m, diagnostics.smoothing.maxlambdastepafter_m,
                        diagnostics.smoothing.maxmove_m, diagnostics.projection.maxborderstep_m,
                        diagnostics.projection.borderstep99_m, diagnostics.projection.maxparentmeanerror_m,
                        diagnostics.meanbound_m, diagnostics.projection.bedsmovedparents,
                        diagnostics.projection.forcedparents, diagnostics.landbelowsea);
        return result;
    };
    const RefinementResult a = run(p, 2, "smoothing P");
    const RefinementResult b = run(q, 2, "smoothing Q");
    expectexact(a, "smoothing P");
    expectexact(b, "smoothing Q");
    const OverlapCounts counts = compareoverlaplevels(a, b, "smoothing P/Q");
    EXPECT_GT(counts.lake, 0);
    EXPECT_GT(counts.channel, 0);
    for (const LevelDiagnostics& level : a.levels)
    {
        EXPECT_NE(std::find(level.steps.begin(), level.steps.end(), landscape::LevelStep::smoothing), level.steps.end());
        EXPECT_EQ(level.smoothing.passes, 6);
        EXPECT_NEAR(level.smoothing.maxlambdaafter_m, level.projection.maxlambda_m, 1.0e-6);
        EXPECT_LT(level.smoothing.maxlambdaafter_m, level.smoothing.maxlambdabefore_m);
    }

    // The level consistency: level 1 of two levels is one level bit for bit.
    const RefinementResult one = run(p, 1, "smoothing P, one level");
    const auto cells = sharedcells(one.levels[0].window, a.levels[0].window);
    ASSERT_EQ(cells.size(), one.levels[0].window.keptsize());
    EXPECT_EQ(differing(one.fields.elevation_m, a.coarserelevation_m[0], cells), 0);
}

// Step S2.8: the masks inherited and refined within one parent cell, every feature on, on P (the
// lake) and R (the coast), two levels. Level 2's sea lies inside level 1's sea dilated by one
// level-1 cell; no dry cell (not sea, not channel) lies below the surface of a wet 4-neighbour at
// either level (the flood wets every cell below a lake's surface connected to its water, step 5
// raises the rest, the projection keeps dry land above the water beside it). The overlap identity
// with the inherited masks on is RefinementLevelsOverlap's (every feature on: P and Q across the
// date line, differing in column and row origin, with the lake; R and S with the coast).
TEST_F(Landscape, RefinementLevelsInheritedMasks)
{
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    // Dry cells of a kept field below the surface of a wet 4-neighbour.
    const auto drybelow = [](const std::vector<float>& elevation, const std::vector<uint8_t>& sea,
                             const std::vector<float>& lake, const std::vector<uint8_t>* direction, int columns, int rows) {
        int count = 0;
        for (int y = 0; y < rows; y++)
            for (int x = 0; x < columns; x++)
            {
                const std::size_t cell = at(x, y, columns);
                if (lake[cell] > 0.0f || sea[cell] != 0 || (direction != nullptr && (*direction)[cell] != 0))
                    continue;
                const int offsets[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
                bool below = false;
                for (const auto& offset : offsets)
                {
                    const int nx = x + offset[0], ny = y + offset[1];
                    if (nx < 0 || ny < 0 || nx >= columns || ny >= rows)
                        continue;
                    const float surface = lake[at(nx, ny, columns)];
                    below = below || (surface > 0.0f && elevation[cell] < surface);
                }
                count += below ? 1 : 0;
            }
        return count;
    };
    for (const RegionalWindow& window : {windowp(), coastwest()})
    {
        const std::string name = "inherited masks at column " + std::to_string(window.column0);
        const RefinementResult result = refine(inputs, window, 2, parameters, name);
        ASSERT_EQ(result.coarsersea.size(), 1u);
        ASSERT_EQ(result.coarserlakesurface_m.size(), 1u);
        expectexact(result, name);
        const LevelWindow& one = result.levels[0].window;
        const LevelWindow& two = result.levels[1].window;
        const std::vector<uint8_t>& parentsea = result.coarsersea[0];
        const RegionalFields& fields = result.fields;

        // Level 2's sea inside level 1's sea dilated by one level-1 cell.
        int sea = 0, landparent = 0, beyond = 0;
        const std::int64_t widthone = one.levelwidth();
        for (int y = 0; y < two.rows; y++)
            for (int x = 0; x < two.columns; x++)
            {
                if (fields.sea[at(x, y, two.columns)] == 0)
                    continue;
                sea++;
                const std::int64_t column = landscape::wraplevel(two.column0 + x, two.levelwidth()) / 4;
                const std::int64_t row = (two.row0 + y) / 4;
                bool own = false, nearby = false;
                for (int dr = -1; dr <= 1; dr++)
                    for (int dc = -1; dc <= 1; dc++)
                    {
                        const std::int64_t px = landscape::wraplevel(column + dc - one.column0, widthone);
                        const std::int64_t py = row + dr - one.row0;
                        ASSERT_TRUE(px < one.columns && py >= 0 && py < one.rows) << name;
                        const bool parent = parentsea[at(static_cast<int>(px), static_cast<int>(py), one.columns)] != 0;
                        nearby = nearby || parent;
                        own = own || (dr == 0 && dc == 0 && parent);
                    }
                landparent += own ? 0 : 1;
                beyond += nearby ? 0 : 1;
            }
        const int drytwo = drybelow(fields.elevation_m, fields.sea, fields.lakesurface_m, &fields.riverdirection,
                                    two.columns, two.rows);
        const int dryone = drybelow(result.coarserelevation_m[0], parentsea, result.coarserlakesurface_m[0], nullptr,
                                    one.columns, one.rows);
        int wet = 0;
        for (const float surface : fields.lakesurface_m)
            wet += surface > 0.0f ? 1 : 0;
        std::printf("[refinement] %s: level 2 sea %d (%d in land parents, %d beyond one parent cell), wet %d; dry below a "
                    "lake: level 1 %d (channels included), level 2 %d (diagnostics %d, %d)\n",
                    name.c_str(), sea, landparent, beyond, wet, dryone, drytwo, result.levels[0].drybelowlake,
                    result.levels[1].drybelowlake);
        EXPECT_EQ(beyond, 0) << name;
        EXPECT_EQ(result.levels[1].seabeyondparent, 0) << name;
        EXPECT_EQ(result.levels[1].seainlandparents, landparent) << name;
        EXPECT_EQ(drytwo, 0) << name;
        EXPECT_EQ(result.levels[0].drybelowlake, 0) << name;
        EXPECT_EQ(result.levels[1].drybelowlake, 0) << name;
        if (window.column0 == windowp().column0)
            EXPECT_GT(wet, 0) << name;
        else
            EXPECT_GT(sea, 0) << name;
    }
}

// Step S2.8 review: the lake flood at the loop's smallest halo with the lake cell (global cell 248,
// 71) on a kept edge of both windows: E keeps columns up to 248 (the zone's column 249 lies in its
// halo), F from 248 across the date line (the zone's column 247 in its halo). Every level's shared
// kept cells (the lake cell's column) bit-identical, with every feature on, and the lake wet there.
TEST_F(Landscape, RefinementLevelsOverlapLakeAtEdge)
{
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    const RegionalWindow e = landscape::regionalwindow(256, 16, loophalo, 155.0, -15.0, 169.5, -5.0);
    const RegionalWindow f = landscape::regionalwindow(256, 16, loophalo, 169.5, -15.0, -175.0, -5.0);
    ASSERT_EQ(e.column0 + e.columns - 1, 248);
    ASSERT_EQ(f.column0, 248);
    ASSERT_GT(inputs.lakesurface_m[at(248, 71, inputs.width)], 0.0f);
    const RefinementResult a = refine(inputs, e, 2, parameters, "lake at the edge E");
    const RefinementResult b = refine(inputs, f, 2, parameters, "lake at the edge F");
    expectexact(a, "E");
    expectexact(b, "F");
    const OverlapCounts counts = compareoverlaplevels(a, b, "lake at the edge E/F");
    EXPECT_GT(counts.lake, 0);
}

// Step S2.4: the source level of a request and the resample of a refinement's final level,
// planet-free: a result built by hand on window C's chain at factor 4, one level (kept 96 x
// 28 cells of pi R / 128).
TEST_F(Landscape, RefinementLevelsResampledFields)
{
    const double cell_m = pi * radius_m / 128.0;
    EXPECT_EQ(landscape::resamplesourcelevel(cell_m, radius_m, 64, 4, 4), 1);
    EXPECT_EQ(landscape::resamplesourcelevel(4.0 * cell_m, radius_m, 64, 4, 4), 0);
    EXPECT_EQ(landscape::resamplesourcelevel(cell_m * 0.999, radius_m, 64, 4, 4), 2);
    EXPECT_EQ(landscape::resamplesourcelevel(cell_m / 4.0, radius_m, 64, 4, 2), 2);
    EXPECT_EQ(landscape::resamplesourcelevel(cell_m / 4.5, radius_m, 64, 4, 2), -1);
    EXPECT_THROW(landscape::resamplesourcelevel(0.0, radius_m, 64, 4, 2), std::invalid_argument);

    RefinementResult result;
    result.chain = landscape::halochain(windowc(4), 4, 1, landscape::prolongationreach);
    const LevelWindow& level = result.chain.back();
    const int columns = level.columns, rows = level.rows;
    ASSERT_EQ(columns, 96);
    ASSERT_EQ(rows, 28);
    landscape::RegionalFields& f = result.fields;
    f.window = windowc(4);
    const std::size_t size = level.keptsize();
    for (std::vector<float>* field : {&f.elevation_m, &f.discharge_m3s, &f.lakesurface_m, &f.januarytemperature_c,
                                      &f.julytemperature_c, &f.januaryprecipitation_mm, &f.julyprecipitation_mm,
                                      &f.icethickness_m})
        field->resize(size);
    f.sea.resize(size);
    f.riverdirection.resize(size);
    f.januaryflow.resize(size);
    f.julyflow.resize(size);
    f.climate.resize(size);
    for (int y = 0; y < rows; y++)
        for (int x = 0; x < columns; x++)
        {
            const std::size_t cell = at(x, y, columns);
            f.elevation_m[cell] = static_cast<float>(100 * y + x);
            // Rows 0..9 a checkerboard of sea (a tie in every 2 x 2 block), below one cell in four.
            const bool sea = y < 10 ? (x + y) % 2 == 0 : x % 2 == 0 && y % 2 == 0;
            f.sea[cell] = sea ? 1 : 0;
            // Rows 20.. lakes on the odd columns (half of every 2 x 2 block, never sea).
            f.lakesurface_m[cell] = y >= 20 && x % 2 == 1 ? static_cast<float>(300 + x) : 0.0f;
            f.discharge_m3s[cell] = static_cast<float>((7 * x + 13 * y) % 17);
            f.riverdirection[cell] = static_cast<uint8_t>(1 + (x + 3 * y) % 8);
            f.januaryflow[cell] = x * y;
            f.julyflow[cell] = x + y;
            f.climate[cell] = static_cast<int16_t>(100 * y + x);
            f.januarytemperature_c[cell] = 10.0f + 0.5f * static_cast<float>(x);
            f.julytemperature_c[cell] = 20.0f - 0.25f * static_cast<float>(y);
            f.januaryprecipitation_mm[cell] = static_cast<float>(y);
            f.julyprecipitation_mm[cell] = static_cast<float>(x + y);
            f.icethickness_m[cell] = static_cast<float>(x % 3);
        }

    // n = 1: aligned, the level grid itself.
    const landscape::ResampledRegionalFields one = landscape::resampleregionalfields(result, radius_m, cell_m);
    EXPECT_FALSE(one.report.resampled());
    EXPECT_EQ(one.report.sourcelevel, 1);
    EXPECT_EQ(one.width, columns);
    EXPECT_EQ(one.height, rows);
    EXPECT_EQ(one.fields.elevation_m, f.elevation_m);
    EXPECT_EQ(one.fields.sea, f.sea);
    EXPECT_EQ(one.fields.lakesurface_m, f.lakesurface_m);
    EXPECT_EQ(one.fields.discharge_m3s, f.discharge_m3s);
    EXPECT_EQ(one.fields.riverdirection, f.riverdirection);
    EXPECT_EQ(one.fields.januaryflow, f.januaryflow);
    EXPECT_EQ(one.fields.julyflow, f.julyflow);
    EXPECT_EQ(one.fields.climate, f.climate);
    EXPECT_EQ(one.fields.januarytemperature_c, f.januarytemperature_c);
    EXPECT_EQ(one.fields.icethickness_m, f.icethickness_m);

    // n = 2 divides k = 4: aligned, 2 x 2 blocks.
    const landscape::ResampledRegionalFields two = landscape::resampleregionalfields(result, radius_m, 2.0 * cell_m);
    EXPECT_FALSE(two.report.resampled());
    ASSERT_EQ(two.width, 48);
    ASSERT_EQ(two.height, 14);
    for (int q = 0; q < two.height; q++)
        for (int p = 0; p < two.width; p++)
        {
            double elevation = 0.0, surface = 0.0;
            int seacells = 0, wetcells = 0, january = 0, july = 0;
            float discharge = -1.0f;
            uint8_t direction = 0;
            for (int j = 2 * q; j < 2 * q + 2; j++)
                for (int i = 2 * p; i < 2 * p + 2; i++)
                {
                    const std::size_t cell = at(i, j, columns);
                    elevation += f.elevation_m[cell];
                    seacells += f.sea[cell];
                    if (f.sea[cell] == 0 && f.lakesurface_m[cell] > 0.0f)
                    {
                        wetcells++;
                        surface += f.lakesurface_m[cell];
                    }
                    if (f.discharge_m3s[cell] > discharge)
                    {
                        discharge = f.discharge_m3s[cell];
                        direction = f.riverdirection[cell];
                    }
                    january = std::max(january, f.januaryflow[cell]);
                    july = std::max(july, f.julyflow[cell]);
                }
            const std::size_t pixel = at(p, q, two.width);
            const std::string where = "pixel " + std::to_string(p) + ", " + std::to_string(q);
            EXPECT_EQ(two.fields.elevation_m[pixel], static_cast<float>(elevation / 4.0)) << where;
            EXPECT_EQ(two.fields.sea[pixel], seacells >= 2 ? 1 : 0) << where;
            EXPECT_EQ(two.fields.sea[pixel], q < 5 ? 1 : 0) << where;
            const float lake = seacells < 2 && wetcells >= 2 ? static_cast<float>(surface / wetcells) : 0.0f;
            EXPECT_EQ(two.fields.lakesurface_m[pixel], lake) << where;
            EXPECT_EQ(two.fields.lakesurface_m[pixel], q >= 10 ? static_cast<float>(300 + 2 * p + 1) : 0.0f) << where;
            EXPECT_EQ(two.fields.discharge_m3s[pixel], discharge) << where;
            EXPECT_EQ(two.fields.riverdirection[pixel], direction) << where;
            EXPECT_EQ(two.fields.januaryflow[pixel], january) << where;
            EXPECT_EQ(two.fields.julyflow[pixel], july) << where;
            EXPECT_EQ(two.fields.climate[pixel], f.climate[at(2 * p + 1, 2 * q + 1, columns)]) << where;
        }

    // n = 1.5: resampled, partial pixels dropped; overlap weights 1 and 0.5 on the first pixel.
    const landscape::ResampledRegionalFields unaligned = landscape::resampleregionalfields(result, radius_m, 1.5 * cell_m);
    EXPECT_TRUE(unaligned.report.resampled());
    EXPECT_EQ(unaligned.width, 64);
    EXPECT_EQ(unaligned.height, 18);
    EXPECT_NEAR(unaligned.fields.elevation_m[0], (100.0 * 0.5 * 1.5 + 0.5 * 1.5) / 2.25, 1e-4);
    EXPECT_EQ(unaligned.fields.climate[at(1, 1, unaligned.width)], f.climate[at(2, 2, columns)]);  // centre 2.25, 2.25
    EXPECT_EQ(unaligned.fields.elevation_m.size(), unaligned.fields.climate.size());

    // Pixel (0, 0) at n = 1.5 covers cells (0, 0), (1, 0), (0, 1), (1, 1) with weights 1, 0.5,
    // 0.5, 0.25 (half the total: 1.125). Two cells of four are a tie by count, decided by weight.
    struct Corner
    {
        const char* name;
        int sea[4];         // cells (0, 0), (1, 0), (0, 1), (1, 1)
        float lake[4];
        uint8_t sea_expected;
        float lake_expected;
    };
    const Corner corners[] = {
        {"sea on the two half cells (1.0): not sea; lakes on the others (1.25): wet", {0, 1, 1, 0},
         {400.0f, 0.0f, 0.0f, 500.0f}, 0, (400.0f * 1.0f + 500.0f * 0.25f) / 1.25f},
        {"sea on the whole and the quarter cell (1.25): sea", {1, 0, 0, 1}, {0.0f, 0.0f, 0.0f, 0.0f}, 1, 0.0f},
        {"lakes on the two half cells (1.0): dry", {0, 0, 0, 0}, {0.0f, 450.0f, 460.0f, 0.0f}, 0, 0.0f}};
    for (const Corner& corner : corners)
    {
        RefinementResult crafted = result;
        landscape::RegionalFields& c = crafted.fields;
        const std::size_t cells[4] = {at(0, 0, columns), at(1, 0, columns), at(0, 1, columns), at(1, 1, columns)};
        for (int index = 0; index < 4; index++)
        {
            c.sea[cells[index]] = static_cast<uint8_t>(corner.sea[index]);
            c.lakesurface_m[cells[index]] = corner.lake[index];
            c.discharge_m3s[cells[index]] = 1.0f;
            c.riverdirection[cells[index]] = 1;
            c.januaryflow[cells[index]] = 0;
            c.julyflow[cells[index]] = 0;
        }
        c.discharge_m3s[cells[3]] = 50.0f;  // the largest discharge on the quarter-cell sliver
        c.riverdirection[cells[3]] = 7;
        c.januaryflow[cells[1]] = 900;      // the largest flows on half-cell slivers
        c.julyflow[cells[2]] = 800;
        const landscape::ResampledRegionalFields out = landscape::resampleregionalfields(crafted, radius_m, 1.5 * cell_m);
        EXPECT_EQ(out.fields.sea[0], corner.sea_expected) << corner.name;
        EXPECT_FLOAT_EQ(out.fields.lakesurface_m[0], corner.lake_expected) << corner.name;
        EXPECT_EQ(out.fields.discharge_m3s[0], 50.0f) << corner.name;
        EXPECT_EQ(out.fields.riverdirection[0], 7) << corner.name;
        EXPECT_EQ(out.fields.januaryflow[0], 900) << corner.name;
        EXPECT_EQ(out.fields.julyflow[0], 800) << corner.name;
    }

    // The source level must be the final level, and no output finer than the data.
    EXPECT_THROW(landscape::resampleregionalfields(result, radius_m, 5.0 * cell_m), std::invalid_argument);
    EXPECT_THROW(landscape::resampleregionalfields(result, radius_m, 0.5 * cell_m), std::invalid_argument);
    RefinementResult empty;
    EXPECT_THROW(landscape::resampleregionalfields(empty, radius_m, cell_m), std::invalid_argument);
}

namespace
{
// A parent field over the parent window's domain, a function of the local indices.
std::vector<float> parentfield(const LevelWindow& parent, const std::function<float(int, int)>& value)
{
    std::vector<float> field(parent.domainsize());
    for (int y = 0; y < parent.domainrows; y++)
        for (int x = 0; x < parent.domaincolumns; x++)
            field[parent.index(x, y)] = value(x, y);
    return field;
}
}

// Step S3.1: the limited prolongation (refinement.prolongation_limit 1). On window C's chain at
// factor 4 (level-1 parents, level-2 children): a 6 km notch (two parent columns of 1,000 m in a
// ridge of 7,000 m) undershoots below its floor and overshoots above the ridge without the limit,
// which is S2.3's Catmull-Rom bit for bit; with it every child lies in the range of its 2 x 2
// bracketing parents, children already inside it keep their value, and the report counts the
// kept cells changed and the largest move; linear fields are unchanged. Then the loop: a limit
// other than 0 or 1 is refused; with the limit on, every feature on, P and Q (Q across the date
// line) are bit-identical on their overlap at both levels and exact, and level 2 clamps.
TEST_F(Landscape, RefinementLevelsLimitedProlongation)
{
    const std::vector<LevelWindow> chain = landscape::halochain(windowc(16), 4, 2, landscape::prolongationreach);
    const LevelWindow& parent = chain[1];
    const LevelWindow& level = chain[2];
    const int factor = level.factor;
    // As prolong: the child at offset o lies t = (o + 0.5) / factor - 0.5 parent cells from its
    // parent's centre; the stencil's second point is the parent (t >= 0) or its predecessor.
    const auto lower = [factor](int o) { return (o + 0.5) / factor - 0.5 < 0.0 ? -1 : 0; };
    const auto fraction = [factor](int o) {
        const double t = (o + 0.5) / factor - 0.5;
        return t < 0.0 ? t + 1.0 : t;
    };
    const auto parentx = [&](int x) { return parent.localcolumn(level.parentcolumn(x)) + lower(x % factor); };
    const auto parenty = [&](int y) { return static_cast<int>(parent.localrow(level.parentrow(y))) + lower(y % factor); };
    const auto reference = [&](const std::vector<float>& field, int x, int y) {
        const int px = parentx(x), py = parenty(y);
        double rows[4];
        for (int b = 0; b < 4; b++)
        {
            const auto value = [&](int a) { return static_cast<double>(field[parent.index(px - 1 + a, py - 1 + b)]); };
            rows[b] = landscape::catmullrom(value(0), value(1), value(2), value(3), fraction(x % factor));
        }
        return static_cast<float>(landscape::catmullrom(rows[0], rows[1], rows[2], rows[3], fraction(y % factor)));
    };
    const auto bracket = [&](const std::vector<float>& field, int x, int y) {
        const int px = parentx(x), py = parenty(y);
        float lowest = field[parent.index(px, py)], highest = lowest;
        for (int b = 0; b < 2; b++)
            for (int a = 0; a < 2; a++)
            {
                lowest = std::min(lowest, field[parent.index(px + a, py + b)]);
                highest = std::max(highest, field[parent.index(px + a, py + b)]);
            }
        return std::make_pair(lowest, highest);
    };
    const int keptx0 = level.keptcolumnoffset(), kepty0 = level.keptrowoffset();
    const auto kept = [&](int x, int y) {
        return x >= keptx0 && x < keptx0 + level.columns && y >= kepty0 && y < kepty0 + level.rows;
    };

    constexpr float floor_m = 1000.0f, ridge_m = 7000.0f;
    const std::vector<float> notch =
        parentfield(parent, [&](int x, int) { return x % 16 == 7 || x % 16 == 8 ? floor_m : ridge_m; });
    const landscape::Prolongation plain = landscape::prolong(parent, notch, level, false);
    const landscape::Prolongation limited = landscape::prolong(parent, notch, level, true);
    ASSERT_EQ(plain.values.size(), level.domainsize());
    ASSERT_EQ(limited.values.size(), level.domainsize());
    int mismatched = 0, below = 0, above = 0, outside = 0, insidechanged = 0, keptchanged = 0;
    double keptmove = 0.0;
    float lowest = limited.values[0], highest = limited.values[0];
    for (int y = 0; y < level.domainrows; y++)
        for (int x = 0; x < level.domaincolumns; x++)
        {
            const std::size_t cell = level.index(x, y);
            const float before = plain.values[cell], after = limited.values[cell];
            mismatched += before != reference(notch, x, y) ? 1 : 0;
            below += before < floor_m ? 1 : 0;
            above += before > ridge_m ? 1 : 0;
            const auto [low, high] = bracket(notch, x, y);
            outside += after < low || after > high ? 1 : 0;
            insidechanged += before >= low && before <= high && after != before ? 1 : 0;
            lowest = std::min(lowest, after);
            highest = std::max(highest, after);
            if (kept(x, y) && after != before)
            {
                keptchanged++;
                keptmove = std::max(keptmove, std::fabs(static_cast<double>(after) - static_cast<double>(before)));
            }
        }
    std::printf("[refinement] limited prolongation, notch: without the limit %d children below the floor, %d above the "
                "ridge; with it %.1f .. %.1f m, %d kept cells clamped (largest %.1f m)\n",
                below, above, static_cast<double>(lowest), static_cast<double>(highest), limited.clamped, limited.maxclamp_m);
    EXPECT_EQ(mismatched, 0);
    EXPECT_GT(below, 0);
    EXPECT_GT(above, 0);
    EXPECT_EQ(plain.clamped, 0);
    EXPECT_EQ(plain.maxclamp_m, 0.0);
    EXPECT_EQ(outside, 0);
    EXPECT_EQ(insidechanged, 0);
    EXPECT_GE(lowest, floor_m);
    EXPECT_LE(highest, ridge_m);
    EXPECT_GT(limited.clamped, 0);
    EXPECT_EQ(limited.clamped, keptchanged);
    EXPECT_EQ(limited.maxclamp_m, keptmove);
    EXPECT_GT(limited.maxclamp_m, 0.0);

    // Linear fields (exact in float) and a constant: the limit changes nothing.
    const std::vector<std::function<float(int, int)>> linear = {
        [](int x, int y) { return static_cast<float>(31000.0 + 37.5 * x - 12.25 * y); },
        [](int x, int) { return static_cast<float>(-250.0 + 512.0 * x); },
        [](int, int y) { return static_cast<float>(4000.0 - 64.0 * y); },
        [](int, int) { return 31043.0f; }};
    for (std::size_t index = 0; index < linear.size(); index++)
    {
        const std::vector<float> field = parentfield(parent, linear[index]);
        const landscape::Prolongation a = landscape::prolong(parent, field, level, false);
        const landscape::Prolongation b = landscape::prolong(parent, field, level, true);
        EXPECT_EQ(a.values, b.values) << "linear field " << index;
        EXPECT_EQ(b.clamped, 0) << "linear field " << index;
    }

    // The loop refuses a limit other than 0 or 1.
    {
        const RegionalGlobalInputs inputs = syntheticworld(64);
        LevelTarget target;
        target.window = windowc(4);
        target.window.halo = landscape::refinementhalo(2, 0);
        target.level = 1;
        RefinementParameters refinement;
        refinement.prolongationlimit = 2;
        EXPECT_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, refinement), std::invalid_argument);
        refinement.prolongationlimit = 1;
        EXPECT_NO_THROW(landscape::runrefinement(inputs, target, RegionalParameters{}, refinement));
    }

    // The overlap identity and exactness with the limit on.
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    RefinementParameters refinement;
    refinement.prolongationlimit = 1;
    const auto run = [&](const RegionalWindow& window, const std::string& name) {
        LevelTarget target;
        target.window = window;
        target.level = 2;
        RefinementResult result = landscape::runrefinement(inputs, target, parameters, refinement);
        for (const LevelDiagnostics& diagnostics : result.levels)
            std::printf("[refinement] %s, level %d: %d kept cells clamped (largest %.1f m); land below the sea level %d "
                        "(deepest %.1f m); forced %d\n",
                        name.c_str(), diagnostics.level, diagnostics.prolongationclamped,
                        diagnostics.maxprolongationclamp_m, diagnostics.landbelowsea, diagnostics.maxlandbelowsea_m,
                        diagnostics.projection.forcedparents);
        return result;
    };
    const RefinementResult a = run(windowp(), "limited P");
    const RefinementResult b = run(windowq(), "limited Q");
    expectexact(a, "limited P");
    expectexact(b, "limited Q");
    const OverlapCounts counts = compareoverlaplevels(a, b, "limited prolongation P/Q");
    EXPECT_GT(counts.lake, 0);
    EXPECT_GT(counts.channel, 0);
    EXPECT_EQ(a.levels[0].prolongationclamped, 0);
    EXPECT_GT(a.levels[1].prolongationclamped, 0);
}

// Step S3.2: the solve against its parent (refinement.parent_lift 2) with the channel cells' means
// (refinement.channel_cell_slope 0.3) and the smoothing (6 passes), every feature on: P and Q (Q
// across the date line) bit-identical on their overlap at both levels and exact; level 1 of two
// levels is one level bit for bit; the land channel cells lie at or above their beds after the
// projection; level 2's lambda before the smoothing is below the run without the two keys.
TEST_F(Landscape, RefinementLevelsSolveAgainstParent)
{
    const RegionalGlobalInputs inputs = wideworld();
    const RegionalParameters parameters = everyfeature();
    RefinementParameters refinement;
    refinement.smoothpasses = 6;
    refinement.parentlift = 2;
    refinement.channelcellslope = 0.3;
    const int halo = landscape::refinementhalo(2, refinement.smoothpasses);
    const RegionalWindow p = landscape::regionalwindow(256, 16, halo, 160.0, -15.0, 180.0, 10.0);
    const RegionalWindow q = landscape::regionalwindow(256, 16, halo, 170.0, -10.0, -170.0, 15.0);
    const auto run = [&](const RegionalWindow& window, int level, const RefinementParameters& settings,
                         const std::string& name) {
        LevelTarget target;
        target.window = window;
        target.level = level;
        RefinementResult result = landscape::runrefinement(inputs, target, parameters, settings);
        for (const LevelDiagnostics& diagnostics : result.levels)
            std::printf("[refinement] %s, level %d: lambda %.1f -> %.2f m, smoothing move p99 %.1f m; border step max "
                        "%.2f m; land below the sea level %d (deepest %.1f m); channel rise max %.1f m; forced %d\n",
                        name.c_str(), diagnostics.level, diagnostics.smoothing.maxlambdabefore_m,
                        diagnostics.smoothing.maxlambdaafter_m, diagnostics.smoothmove99_m,
                        diagnostics.projection.maxborderstep_m, diagnostics.landbelowsea, diagnostics.maxlandbelowsea_m,
                        diagnostics.maxchannelrise_m, diagnostics.projection.forcedparents);
        return result;
    };
    const RefinementResult a = run(p, 2, refinement, "against the parent P");
    const RefinementResult b = run(q, 2, refinement, "against the parent Q");
    expectexact(a, "against the parent P");
    expectexact(b, "against the parent Q");
    const OverlapCounts counts = compareoverlaplevels(a, b, "against the parent P/Q");
    EXPECT_GT(counts.lake, 0);
    EXPECT_GT(counts.channel, 0);
    for (const LevelDiagnostics& level : a.levels)
        EXPECT_GE(level.maxchannelrise_m, 0.0);

    // The level consistency: level 1 of two levels is one level bit for bit.
    const RefinementResult one = run(p, 1, refinement, "against the parent P, one level");
    const auto cells = sharedcells(one.levels[0].window, a.levels[0].window);
    ASSERT_EQ(cells.size(), one.levels[0].window.keptsize());
    EXPECT_EQ(differing(one.fields.elevation_m, a.coarserelevation_m[0], cells), 0);

    // Against the S2 loop (both keys 0): level 2's lambda before the smoothing is smaller.
    RefinementParameters plain = refinement;
    plain.parentlift = 0;
    plain.channelcellslope = 0.0;
    const RefinementResult before = run(p, 2, plain, "S2 loop P");
    EXPECT_LT(a.levels[1].smoothing.maxlambdabefore_m, before.levels[1].smoothing.maxlambdabefore_m);

    // The keys' range.
    LevelTarget target;
    target.window = p;
    target.level = 1;
    RefinementParameters bad = refinement;
    bad.parentlift = 3;
    EXPECT_THROW(landscape::runrefinement(inputs, target, parameters, bad), std::invalid_argument);
    bad = refinement;
    bad.channelcellslope = -0.1;
    EXPECT_THROW(landscape::runrefinement(inputs, target, parameters, bad), std::invalid_argument);
    bad = refinement;
    bad.smoothbasis = 2;
    EXPECT_THROW(landscape::runrefinement(inputs, target, parameters, bad), std::invalid_argument);
}
