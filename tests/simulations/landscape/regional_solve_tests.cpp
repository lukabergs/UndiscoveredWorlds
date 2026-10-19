// Tests of landscape/regional_solve.hpp (multiscale plan S3.0: the solve step of the hook
// interface, regional_pass.cpp's solve moved unchanged; S3.2 adds the channel cells' means and the
// solve against its parent), suite Landscape:
// TEST_F(Landscape, RegionalSolve<What>) on the synthetic world of regional_test_world.hpp, the
// level's steps up to the channels run through the hook interface (regional_pass.hpp).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/channel_network.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_solve.hpp"
#include "simulations/landscape/regional_surface.hpp"
#include "simulations/landscape/regional_valleys.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

TEST_F(Landscape, RegionalSolveKeepsThePinnedCellsAndTheBeds)
{
    // At refinement.channel_cell_slope 0 (S3.2's term off) the solve leaves every pinned cell
    // (channels, sea, the outer ring) at the working surface, so a land channel cell keeps its
    // bed, writes solved_m as the result and leaves bed_m unchanged.
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    const RegionalWindow window = windowa(4);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const auto contextgrid = landscape::regionalwindowgrid(window, inputs.radius_m);
    const landscape::LevelContext context{inputs, window, *contextgrid, parameters, warp, incision, details};
    landscape::LevelState state;
    landscape::noiselevel(context, state);
    landscape::presmoothinglevel(context, state);
    landscape::coastlevel(context, state);
    landscape::incisionlevel(context, state);
    landscape::channelslevel(context, state);
    const std::vector<float> surface = state.surface_m;
    const std::vector<float> beds = state.bed_m;

    landscape::solvelevel(context, state);
    const std::size_t count = surface.size();
    ASSERT_EQ(state.elevation_m.size(), count);
    ASSERT_EQ(state.pinned.size(), count);
    EXPECT_EQ(state.solved_m, state.elevation_m);
    EXPECT_EQ(state.bed_m, beds);
    EXPECT_EQ(state.surface_m, surface);
    std::size_t channelcells = 0, unpinned = 0;
    for (std::size_t cell = 0; cell < count; cell++)
    {
        if (state.pinned[cell] != 0)
            EXPECT_EQ(state.elevation_m[cell], surface[cell]) << cell;
        else
            unpinned++;
        if (state.channels.channel[cell] != 0 && state.sea[cell] == 0)
        {
            EXPECT_EQ(state.elevation_m[cell], state.bed_m[cell]) << cell;
            channelcells++;
        }
    }
    EXPECT_GT(channelcells, 0u);
    EXPECT_GT(unpinned, 0u);
    EXPECT_EQ(state.solveddischarge_m3_yr.size(), count);
}

namespace
{
// The steps of a level up to the solve through the hook interface (regional_pass.hpp).
landscape::LevelState stepstothesolve(const landscape::LevelContext& context, landscape::LevelState state)
{
    landscape::noiselevel(context, state);
    landscape::presmoothinglevel(context, state);
    landscape::coastlevel(context, state);
    landscape::incisionlevel(context, state);
    landscape::channelslevel(context, state);
    landscape::solvelevel(context, state);
    return state;
}
}

TEST_F(Landscape, RegionalSolveChannelCellMeans)
{
    // S3.2: in the refinement loop with refinement.channel_cell_slope theta > 0 a land channel cell
    // takes its cell mean b + min(h, theta c_l) (h the mean height above b of its non-channel land
    // 8-neighbours, c_l the cell size): within [b, b + theta c_l], above b where its neighbours are;
    // every other cell, state.bed_m and the solve's diagnostics as with theta 0, which equals the
    // single pass bit for bit; solved_m is the result (no parent: no lift).
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalParameters parameters;
    const RegionalWindow window = windowa(4);
    const landscape::FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    const auto singlegrid = landscape::regionalwindowgrid(window, inputs.radius_m);
    const landscape::LevelContext single{inputs, window, *singlegrid, parameters, warp, incision, details};
    landscape::RefinementParameters off;
    const auto withoffgrid = landscape::regionalwindowgrid(window, inputs.radius_m);
    const landscape::LevelContext withoff{inputs, window, *withoffgrid, parameters, warp, incision, details, 1, 4, &off};
    const landscape::LevelState plain = stepstothesolve(single, {});
    const landscape::LevelState zero = stepstothesolve(withoff, {});
    EXPECT_EQ(zero.elevation_m, plain.elevation_m);
    EXPECT_EQ(zero.solved_m, plain.solved_m);
    // theta 0.3 (the cap far above any h at k = 4 on 32 columns) and 0.0001 (a cap of about 30 m).
    for (const double theta : {0.3, 0.0001})
    {
        landscape::RefinementParameters on;
        on.channelcellslope = theta;
        const auto withongrid = landscape::regionalwindowgrid(window, inputs.radius_m);
        const landscape::LevelContext withon{inputs, window, *withongrid, parameters, warp, incision, details, 1, 4, &on};
        const landscape::LevelState means = stepstothesolve(withon, {});
        EXPECT_EQ(means.bed_m, plain.bed_m);
        EXPECT_EQ(means.solved_m, means.elevation_m);
        EXPECT_EQ(means.meanlowering_m, plain.meanlowering_m);
        EXPECT_EQ(means.maxlowering_m, plain.maxlowering_m);

        const int width = window.domainfinewidth(), height = window.domainfineheight();
        const double cap = on.channelcellslope * landscape::northsouthspacing_m(*landscape::regionalwindowgrid(window, inputs.radius_m));
        int channelcells = 0, raised = 0, capped = 0, mismatched = 0;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                const std::size_t cell = at(x, y, width);
                const float value = means.elevation_m[cell];
                if (means.channels.channel[cell] == 0 || means.sea[cell] != 0)
                {
                    mismatched += value != plain.elevation_m[cell] ? 1 : 0;
                    continue;
                }
                channelcells++;
                const double bed = static_cast<double>(means.bed_m[cell]);
                double sum = 0.0;
                int neighbours = 0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        const int nx = x + dx, ny = y + dy;
                        if ((dx == 0 && dy == 0) || nx < 0 || ny < 0 || nx >= width || ny >= height)
                            continue;
                        const std::size_t other = at(nx, ny, width);
                        if (means.channels.channel[other] != 0 || means.sea[other] != 0)
                            continue;
                        sum += static_cast<double>(plain.elevation_m[other]) - bed;
                        neighbours++;
                    }
                const double h = neighbours > 0 ? std::max(0.0, sum / neighbours) : 0.0;
                EXPECT_EQ(value, static_cast<float>(bed + std::min(h, cap))) << cell;
                EXPECT_GE(value, means.bed_m[cell]) << cell;
                EXPECT_LE(static_cast<double>(value), bed + cap + 1.0e-3) << cell;
                raised += value > means.bed_m[cell] ? 1 : 0;
                capped += h > cap ? 1 : 0;
            }
        std::printf("[solve] channel means, theta %.4g: %d land channel cells, %d above their beds, %d at the cap (%.1f m)\n",
                    theta, channelcells, raised, capped, cap);
        EXPECT_GT(channelcells, 0);
        EXPECT_GT(raised, 0);
        EXPECT_EQ(mismatched, 0);
        if (theta < 0.01)
            EXPECT_GT(capped, 0);
    }
}

TEST_F(Landscape, RegionalSolveAgainstItsParent)
{
    // S3.2: the solve against its parent at level 2. Level 1 is a real level of the loop (every step
    // through the hook interface, the projection onto the global elevation, 6 smoothing passes) with a
    // long elapsed time, so its solve lowers the land and its projection lifts it back: its lift
    // (final less solved) is large. At level 2, whose solve rebuilds the land from its own beds, lambda
    // = v - mean(F*) after the solve over the kept parents that are not sea is smaller with
    // refinement.parent_lift 2 than without it (90th percentile and mean), and every unpinned cell is its solve
    // plus the prolonged lift where the floor does not hold it.
    const RegionalGlobalInputs inputs = syntheticworld(64);
    RegionalParameters parameters;
    parameters.relief.elapsed_myr = 200.0;  // the knickpoints reach most of the land: the solve rebuilds it from the beds
    landscape::RefinementParameters refinement;
    refinement.smoothpasses = 6;
    refinement.parentlift = 2;
    const std::vector<landscape::LevelWindow> chain = landscape::halochain(windowc(16), 4, 2, landscape::prolongationreach);
    const landscape::GlobalIncision incision = landscape::globalincision(inputs, parameters);
    const landscape::RegionalDetails details{};
    std::vector<float> levelzero(chain[0].domainsize());
    for (int y = 0; y < chain[0].domainrows; y++)
        for (int x = 0; x < chain[0].domaincolumns; x++)
            levelzero[chain[0].index(x, y)] = inputs.elevation_m[at(static_cast<int>(chain[0].absolutecolumn(x)),
                                                                    static_cast<int>(chain[0].absoluterow(y)), inputs.width)];

    // Level 1.
    const RegionalWindow windowone = landscape::regionalwindowof(chain[1]);
    const landscape::FineWarp warpone = landscape::regionalwarp(windowone, inputs.seed, parameters);
    const auto onegrid = landscape::regionalwindowgrid(windowone, inputs.radius_m);
    const landscape::LevelContext one{inputs, windowone, *onegrid, parameters, warpone, incision, details, 1, 4, &refinement,
                                      &chain[1], &chain[0], &levelzero, nullptr};
    landscape::LevelState first;
    landscape::regionallevelsteps(one, first);
    landscape::ParentFields parent;
    parent.window = chain[1];
    parent.elevation_m = first.elevation_m;
    parent.lift_m = first.lift_m;
    parent.solved_m = first.solved_m;
    parent.sea = first.sea;
    parent.lakesurface_m = first.water.lakesurface_m;
    parent.channels = first.channels;
    parent.bed_m = first.bed_m;
    std::vector<float> lift(parent.elevation_m.size());
    double largestlift = 0.0;
    for (std::size_t cell = 0; cell < lift.size(); cell++)
    {
        lift[cell] = static_cast<float>(static_cast<double>(parent.elevation_m[cell]) - static_cast<double>(parent.solved_m[cell]));
        largestlift = std::max(largestlift, std::fabs(static_cast<double>(lift[cell])));
    }

    // Level 2 with and without the parent's lift.
    const landscape::LevelWindow& level = chain[2];
    const RegionalWindow windowtwo = landscape::regionalwindowof(level);
    const landscape::FineWarp warptwo = landscape::regionalwarp(windowtwo, inputs.seed, parameters);
    const auto run = [&](int mode) {
        landscape::RefinementParameters settings = refinement;
        settings.parentlift = mode;
        const auto contextgrid = landscape::regionalwindowgrid(windowtwo, inputs.radius_m);
        const landscape::LevelContext context{inputs, windowtwo, *contextgrid, parameters, warptwo, incision, details, 2, 4, &settings,
                                              &level, &chain[1], &parent.elevation_m, &parent};
        landscape::LevelState state;
        state.start.start_m = landscape::prolong(chain[1], parent.elevation_m, level, false).values;
        return stepstothesolve(context, std::move(state));
    };
    const landscape::LevelState with = run(2);
    const landscape::LevelState without = run(0);

    const int f = level.factor;
    const int keptx0 = level.keptcolumnoffset(), kepty0 = level.keptrowoffset();
    double largestwith = 0.0, largestwithout = 0.0, sumwith = 0.0, sumwithout = 0.0;
    std::vector<double> lambdaswith, lambdaswithout;
    int parents = 0;
    for (int q = kepty0 / f; q < (kepty0 + level.rows) / f; q++)
        for (int p = keptx0 / f; p < (keptx0 + level.columns) / f; p++)
        {
            const std::size_t cell = chain[1].index(chain[1].localcolumn(level.parentcolumn(p * f)),
                                                    static_cast<int>(chain[1].localrow(level.parentrow(q * f))));
            if (parent.sea[cell] != 0)
                continue;
            double a = 0.0, b = 0.0;
            for (int y = q * f; y < (q + 1) * f; y++)
                for (int x = p * f; x < (p + 1) * f; x++)
                {
                    a += static_cast<double>(with.elevation_m[level.index(x, y)]);
                    b += static_cast<double>(without.elevation_m[level.index(x, y)]);
                }
            const double v = static_cast<double>(parent.elevation_m[cell]);
            const double lambdawith = std::fabs(v - a / (f * f)), lambdawithout = std::fabs(v - b / (f * f));
            largestwith = std::max(largestwith, lambdawith);
            largestwithout = std::max(largestwithout, lambdawithout);
            sumwith += lambdawith;
            sumwithout += lambdawithout;
            lambdaswith.push_back(lambdawith);
            lambdaswithout.push_back(lambdawithout);
            parents++;
        }
    const auto percentile90 = [](std::vector<double> values) {
        const std::size_t rank = (90 * values.size() + 99) / 100 - 1;  // ceil(0.9 n) - 1
        std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(rank), values.end());
        return values[rank];
    };
    const double p90with = percentile90(lambdaswith), p90without = percentile90(lambdaswithout);
    const std::vector<float> prolonged = landscape::prolong(chain[1], lift, level, false).values;
    int unpinned = 0, mismatched = 0;
    for (std::size_t cell = 0; cell < with.elevation_m.size(); cell++)
    {
        if (with.pinned[cell] != 0)
            continue;
        unpinned++;
        const float expected =
            static_cast<float>(static_cast<double>(with.solved_m[cell]) + static_cast<double>(prolonged[cell]));
        mismatched += with.elevation_m[cell] < expected ? 1 : 0;  // the floor only raises
    }
    std::printf("[solve] against its parent: level-1 lift up to %.1f m; level 2, %d kept land parents: lambda after the "
                "solve largest %.1f m, p90 %.1f m and mean %.2f m with parent_lift 2, %.1f, %.1f and %.2f m without; %d of %d unpinned "
                "cells below their solve plus the lift\n",
                largestlift, parents, largestwith, p90with, sumwith / parents, largestwithout, p90without, sumwithout / parents,
                mismatched,
                unpinned);
    EXPECT_GT(parents, 0);
    EXPECT_GT(largestlift, 100.0);
    // The largest is not asserted: a parent whose children's beds follow the lifted start carries the
    // lift twice (the pinned beds are not de-lifted, regional_solve.cpp); measured 2,601 against 2,521 m.
    EXPECT_LT(p90with, p90without);
    EXPECT_LT(sumwith, 0.75 * sumwithout);  // measured: mean 346 against 715 m
    EXPECT_EQ(mismatched, 0);
}
