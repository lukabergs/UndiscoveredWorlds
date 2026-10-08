// Tests of landscape/regional_pass.hpp (erosion plan Phase E, step E3), suite Landscape:
// TEST_F(Landscape, Regional<What>) on the synthetic world of regional_test_world.hpp (a smooth
// continent across the date line with a ridge, a sea, a network routed by steepest
// descent on the Priority-Flood surface, a lake on one cell, climate varying with
// latitude), at scales 4, 5 and 8, default parameters except the halo (`testhalo`; the
// pass requires at least 2 x meanpasses + 2, RegionalPassRejectsASmallHalo).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

TEST_F(Landscape, RegionalOverlap)
{
    const RegionalGlobalInputs inputs = syntheticworld(32);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow wa = windowa(scale), wb = windowb(scale);
        ASSERT_EQ(wa.column0, 26);
        ASSERT_EQ(wa.columns, 10);
        ASSERT_EQ(wb.column0, 0);
        ASSERT_EQ(wb.row0, 5);
        const OverlapCounts counts = compareoverlap(inputs, wa, wb, "columns A/B");
        EXPECT_EQ(counts.shared, 4 * 5 * scale * scale);
    }
}

TEST_F(Landscape, RegionalOverlapInRows)
{
    const RegionalGlobalInputs inputs = syntheticworld(64);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow wc = windowc(scale), wd = windowd(scale);
        ASSERT_EQ(wc.column0, 40);
        ASSERT_EQ(wc.columns, 24);
        ASSERT_EQ(wc.row0, 11);
        ASSERT_EQ(wc.rows, 7);
        ASSERT_EQ(wd.column0, 52);
        ASSERT_EQ(wd.columns, 20);
        ASSERT_EQ(wd.row0, 14);
        ASSERT_EQ(wd.rows, 7);
        ASSERT_GT(wc.domainrow0, 0);
        ASSERT_LT(wd.domainrow0 + wd.domainrows, 32);
        const OverlapCounts counts = compareoverlap(inputs, wc, wd, "rows C/D");
        EXPECT_EQ(counts.shared, 12 * 4 * scale * scale);
        EXPECT_GT(counts.sea, 0);
        EXPECT_GT(counts.lake, 0);
        EXPECT_GT(counts.ice, 0);
        EXPECT_GT(counts.channel, 0);
    }
}

TEST_F(Landscape, RegionalOverlapLowRelief)
{
    // The synthetic elevations scaled by 0.01 about the sea level (plains of 0 to 38 m,
    // neighbour relief 3 to 6 m), at k = 8: the solve's reach on plains. The noise at
    // its defaults (about 21 m) is tapered over coast_band_m plus that amplitude, so it
    // reaches about half its amplitude inland and fades towards the coast.
    const RegionalGlobalInputs columns = syntheticworld(32, 0.01);
    const OverlapCounts first = compareoverlap(columns, windowa(8), windowb(8), "low relief A/B");
    EXPECT_GT(first.channel, 0);
    const RegionalGlobalInputs rows = syntheticworld(64, 0.01);
    const OverlapCounts second = compareoverlap(rows, windowc(8), windowd(8), "low relief C/D");
    EXPECT_GT(second.sea, 0);
    EXPECT_GT(second.lake, 0);
    EXPECT_GT(second.channel, 0);
}

TEST_F(Landscape, RegionalPassProperties)
{
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const int width = inputs.width;
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowa(scale);
        const RegionalFields fields = landscape::runregionalpass(inputs, window, RegionalParameters{});
        const int fw = window.finewidth(), fh = window.fineheight();
        const std::size_t size = static_cast<std::size_t>(fw) * static_cast<std::size_t>(fh);

        // The kept window has the window's size, in every field.
        ASSERT_EQ(fields.elevation_m.size(), size);
        ASSERT_EQ(fields.sea.size(), size);
        ASSERT_EQ(fields.riverdirection.size(), size);
        ASSERT_EQ(fields.discharge_m3s.size(), size);
        ASSERT_EQ(fields.januaryflow.size(), size);
        ASSERT_EQ(fields.julyflow.size(), size);
        ASSERT_EQ(fields.lakesurface_m.size(), size);
        ASSERT_EQ(fields.januarytemperature_c.size(), size);
        ASSERT_EQ(fields.julytemperature_c.size(), size);
        ASSERT_EQ(fields.januaryprecipitation_mm.size(), size);
        ASSERT_EQ(fields.julyprecipitation_mm.size(), size);
        ASSERT_EQ(fields.icethickness_m.size(), size);
        ASSERT_EQ(fields.climate.size(), size);

        // No NaN, no negative precipitation or ice.
        int nonfinite = 0, negative = 0;
        for (std::size_t cell = 0; cell < size; cell++)
        {
            for (float value : {fields.elevation_m[cell], fields.discharge_m3s[cell], fields.lakesurface_m[cell],
                                fields.januarytemperature_c[cell], fields.julytemperature_c[cell],
                                fields.januaryprecipitation_mm[cell], fields.julyprecipitation_mm[cell], fields.icethickness_m[cell]})
                nonfinite += std::isfinite(value) ? 0 : 1;
            negative += fields.januaryprecipitation_mm[cell] < 0.0f || fields.julyprecipitation_mm[cell] < 0.0f ||
                                fields.icethickness_m[cell] < 0.0f
                            ? 1
                            : 0;
        }
        EXPECT_EQ(nonfinite, 0);
        EXPECT_EQ(negative, 0);

        // Faithful: every global cell's centre fine cell carries the global mask.
        int centres = 0, wrong = 0, seacentres = 0;
        for (int y = 0; y < fh; y++)
            for (int x = 0; x < fw; x++)
            {
                if (!centrecell(window, x, y))
                    continue;
                const std::size_t global = globalof(window, width, x, y);
                centres++;
                seacentres += inputs.network.sea[global];
                wrong += fields.sea[at(x, y, fw)] != inputs.network.sea[global] ? 1 : 0;
            }
        EXPECT_EQ(centres, window.columns * window.rows);
        EXPECT_GT(seacentres, 0);
        EXPECT_LT(seacentres, centres);
        EXPECT_EQ(wrong, 0);

        // Lakes: the global surface where the fine bed lies below it (on a channel cell
        // at or below it, R1), never sea.
        int lakecells = 0;
        for (int y = 0; y < fh; y++)
            for (int x = 0; x < fw; x++)
            {
                const std::size_t cell = at(x, y, fw);
                if (fields.lakesurface_m[cell] <= 0.0f)
                    continue;
                lakecells++;
                const std::size_t global = globalof(window, width, x, y);
                EXPECT_EQ(fields.lakesurface_m[cell], inputs.lakesurface_m[global]);
                if (fields.riverdirection[cell] != 0)
                    EXPECT_LE(fields.elevation_m[cell], fields.lakesurface_m[cell]);
                else
                    EXPECT_LT(fields.elevation_m[cell], fields.lakesurface_m[cell]);
                EXPECT_EQ(fields.sea[cell], 0);
            }
        EXPECT_GT(lakecells, 0);

        // Temperature over the sea: the interpolated value exactly where the
        // interpolated elevation lies at or below the sea level (no lapse correction
        // there), the interpolated value plus the lapse rate x the interpolated height
        // above the sea level elsewhere.
        const std::vector<float> interpolated =
            keptpart(landscape::upsample(inputs.januarytemperature_c, width, inputs.height, window), window);
        const std::vector<float> base = keptpart(landscape::upsample(inputs.elevation_m, width, inputs.height, window), window);
        int seaexact = 0, seaabove = 0;
        for (std::size_t cell = 0; cell < size; cell++)
        {
            if (fields.sea[cell] == 0)
                continue;
            if (base[cell] <= sealevel)
            {
                seaexact++;
                EXPECT_EQ(fields.januarytemperature_c[cell], interpolated[cell]);
            }
            else
            {
                seaabove++;
                EXPECT_NEAR(fields.januarytemperature_c[cell], interpolated[cell] + 6.5 / 1000.0 * (base[cell] - sealevel), 1e-4);
            }
        }
        EXPECT_GT(seaexact, 0);

        // Channel cells' directions lead to the sea, to a land sink's centre or out of
        // the window, without a loop.
        int channels = 0, toseas = 0, outside = 0, failures = 0;
        for (int y = 0; y < fh; y++)
            for (int x = 0; x < fw; x++)
            {
                if (fields.riverdirection[at(x, y, fw)] == 0)
                    continue;
                channels++;
                int cx = x, cy = y;
                bool ended = false;
                for (std::size_t step = 0; step <= size && !ended; step++)
                {
                    const int direction = fields.riverdirection[at(cx, cy, fw)];
                    if (fields.sea[at(cx, cy, fw)] != 0)
                    {
                        toseas++;
                        ended = true;
                        break;
                    }
                    if (direction == 0)
                    {
                        const std::size_t global = globalof(window, width, cx, cy);
                        const bool sink = centrecell(window, cx, cy) && inputs.network.direction[global] == 0;
                        failures += sink ? 0 : 1;
                        ended = true;
                        break;
                    }
                    cx += landscape::directionoffsetx(direction);
                    cy += landscape::directionoffsety(direction);
                    if (cx < 0 || cy < 0 || cx >= fw || cy >= fh)
                    {
                        outside++;
                        ended = true;
                    }
                }
                failures += ended ? 0 : 1;
            }
        EXPECT_GT(channels, 0);
        EXPECT_GT(toseas, 0);
        EXPECT_EQ(failures, 0);
        std::printf("[regional] k %d: %d channel cells, %d reach the sea, %d leave the window; %d lake cells; "
                    "sea cells at the interpolated temperature %d, lapse-corrected %d\n",
                    scale, channels, toseas, outside, lakecells, seaexact, seaabove);
    }
}

TEST_F(Landscape, RegionalPassChannelsEndAtTheFineCoast)
{
    // Phase R, R1: the sea is decided before the channels, so a channel cell that leads
    // on is land, its bed above the sea level and at or above the next cell's (a link's
    // last cell in the sea keeps its seabed); a channel cell of a global lake cell whose
    // bed lies at or below the lake surface carries the lake's water. Wet cells are left
    // out of the bed checks: on real data a lowered lake cell may lie below the sea level
    // (a lake below it) and below its outlet's sill (regional_water.hpp).
    const RegionalGlobalInputs small = syntheticworld(32);
    const RegionalGlobalInputs large = syntheticworld(64);
    int lakechannels = 0;
    for (int scale : {4, 5, 8})
        for (const bool wide : {false, true})
        {
            SCOPED_TRACE(testing::Message() << "k " << scale << (wide ? ", window C" : ", window A"));
            const RegionalGlobalInputs& inputs = wide ? large : small;
            const RegionalWindow window = wide ? windowc(scale) : windowa(scale);
            const RegionalFields fields = landscape::runregionalpass(inputs, window, RegionalParameters{});
            const int fw = window.finewidth(), fh = window.fineheight();
            int channels = 0;
            for (int y = 0; y < fh; y++)
                for (int x = 0; x < fw; x++)
                {
                    const std::size_t cell = at(x, y, fw);
                    const int direction = fields.riverdirection[cell];
                    if (direction == 0)
                        continue;
                    channels++;
                    EXPECT_EQ(fields.sea[cell], 0) << x << ", " << y;
                    EXPECT_GT(fields.elevation_m[cell], sealevel) << x << ", " << y;
                    const int nx = x + landscape::directionoffsetx(direction), ny = y + landscape::directionoffsety(direction);
                    // (Under a lake the elevation may rise along the directions: wet
                    // channel cells are lowered to the surface, regional_water.hpp.)
                    if (nx >= 0 && ny >= 0 && nx < fw && ny < fh && fields.lakesurface_m[cell] <= 0.0f)
                        EXPECT_LE(fields.elevation_m[at(nx, ny, fw)], fields.elevation_m[cell]) << x << ", " << y;
                    const float surface = inputs.lakesurface_m[globalof(window, inputs.width, x, y)];
                    if (surface > 0.0f && fields.elevation_m[cell] <= surface)
                    {
                        lakechannels++;
                        EXPECT_EQ(fields.lakesurface_m[cell], surface) << x << ", " << y;
                    }
                }
            EXPECT_GT(channels, 0);
        }
    EXPECT_GT(lakechannels, 0);
}

TEST_F(Landscape, RegionalPassLakeBelowItsSill)
{
    // R1 (review): a closed lake whose surface lies 20 m below its cell's filled
    // elevation, where the channel beds are pinned (the global cell's centre at the
    // filled elevation). Channel cells at or below the filled elevation carry the lake
    // and are lowered to its surface, so no lake surface lies below its cell; the
    // lake's centre is wet. The overlap stays bit-identical.
    const RegionalGlobalInputs inputs = lakebelowsill(syntheticworld(32), 20.0f);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowa(scale);
        const RegionalFields fields = landscape::runregionalpass(inputs, window, RegionalParameters{});
        const int fw = window.finewidth(), fh = window.fineheight();
        int centres = 0, wetchannels = 0;
        for (int y = 0; y < fh; y++)
            for (int x = 0; x < fw; x++)
            {
                const std::size_t cell = at(x, y, fw);
                const std::size_t global = globalof(window, inputs.width, x, y);
                const float surface = inputs.lakesurface_m[global];
                if (fields.lakesurface_m[cell] > 0.0f)
                {
                    EXPECT_EQ(fields.lakesurface_m[cell], surface) << x << ", " << y;
                    EXPECT_LE(fields.elevation_m[cell], fields.lakesurface_m[cell]) << x << ", " << y;
                    wetchannels += fields.riverdirection[cell] != 0 ? 1 : 0;
                }
                if (surface <= 0.0f || !centrecell(window, x, y))
                    continue;
                centres++;
                ASSERT_GT(surface, 0.0f);
                EXPECT_LT(surface, inputs.network.filledelevation_m[global]);
                EXPECT_EQ(fields.lakesurface_m[cell], surface) << "the lake's centre is dry";
                EXPECT_LE(fields.elevation_m[cell], surface);
            }
        EXPECT_EQ(centres, 1);
        EXPECT_GT(wetchannels, 0);
    }
    const OverlapCounts counts =
        compareoverlap(lakebelowsill(syntheticworld(64), 20.0f), windowc(4), windowd(4), "lake below its sill C/D");
    EXPECT_GT(counts.lake, 0);
}

TEST_F(Landscape, RegionalPassCoastSpeckle)
{
    // Reported, not asserted (E3 review): land fine cells in global sea blocks and sea
    // fine cells in global land blocks away from the coast.
    const RegionalGlobalInputs small = syntheticworld(32);
    const RegionalGlobalInputs large = syntheticworld(64);
    for (int scale : {4, 8})
    {
        int landinsea = 0, seainland = 0;
        speckle(small, windowa(scale), landscape::runregionalpass(small, windowa(scale), RegionalParameters{}), landinsea,
                seainland);
        std::printf("[regional] speckle, window A, k %d: %d land cells in sea blocks, %d sea cells in land blocks\n", scale,
                    landinsea, seainland);
        speckle(large, windowc(scale), landscape::runregionalpass(large, windowc(scale), RegionalParameters{}), landinsea,
                seainland);
        std::printf("[regional] speckle, window C, k %d: %d land cells in sea blocks, %d sea cells in land blocks\n", scale,
                    landinsea, seainland);
    }
}

TEST_F(Landscape, RegionalPassRejectsASmallHalo)
{
    // The halo must be at least 2 x meanpasses + 2 global cells (6 with the default 2).
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow small = landscape::regionalwindow(32, 4, 5, 112.5, -22.5, -135.0, 56.25);
    EXPECT_THROW(landscape::runregionalpass(inputs, small, RegionalParameters{}), std::invalid_argument);
    RegionalParameters onepass;
    onepass.meanpasses = 1;
    EXPECT_NO_THROW(landscape::runregionalpass(inputs, small, onepass));
    onepass.meanpasses = 2;
    EXPECT_NO_THROW(landscape::runregionalpass(inputs, windowa(4), onepass));
}

// The run time of one window of 32 x 16 global cells at k = 16 (plan E3 asks for the
// number, not an assertion): a 64 x 32 world, the default halo of 8. Run with
// world_tests --gtest_also_run_disabled_tests --gtest_filter=Landscape.DISABLED_RegionalPassTiming.
TEST_F(Landscape, DISABLED_RegionalPassTiming)
{
    const RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = landscape::regionalwindow(64, 16, 8, 90.0, -45.0, -90.0, 45.0);
    ASSERT_EQ(window.columns, 32);
    ASSERT_EQ(window.rows, 16);
    const RegionalFields fields = timedpass(inputs, window, "timing window");
    EXPECT_EQ(fields.elevation_m.size(), static_cast<std::size_t>(32 * 16) * 16 * 16);
}
