// Tests of regional_surface.hpp (erosion plan Phase R, step R2), suite Landscape:
// TEST_F(Landscape, RegionalSurface<What>) and RegionalOverlapCoast* on the synthetic world of
// regional_test_world.hpp with the warp, the coast noise and the land-share rule on (test
// values, coastparameters()), at the smallest halo the pass allows (`testhalo`, two
// block-mean passes): the land-share decision reads nothing of the window, so the rule
// halo >= 2 x meanpasses + 2 stays sufficient. Multiscale plan S3.4 adds the level's noise share and
// the relief-amplitude field (RegionalSurfaceLevelNoiseShare, RegionalSurfaceNoiseReliefField).

#include "landscape_fixture.hpp"

#include "regional_test_world.hpp"

#include "simulations/landscape/regional_surface.hpp"

namespace
{
using landscape::FineWarp;
using landscape::StartSurface;

RegionalParameters coastparameters()
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
    return parameters;
}

// The kept window's land share per global cell, checked against the land-share rule; the
// fine sea checked against Phase E's local rules; one-cell holes counted.
struct ShareCheck
{
    int blocks = 0, coastal = 0, outside = 0, seawithoutsea = 0, seainlake = 0, holes = 0, pinches = 0, inconsistent = 0;
    double leastlandshare = 1.0, mostseashare = 0.0;  // over global land cells / land share over global sea cells
};

ShareCheck checkshares(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters)
{
    const FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const StartSurface start = landscape::regionalstartsurface(inputs, window, warp, parameters);
    EXPECT_EQ(start.sharemisses, 0) << "blocks of the domain that miss the land share";
    EXPECT_GE(start.sharepasses, 0);
    std::printf("[regional] land-share offset passes until the domain met the rule: %d\n", start.sharepasses);
    std::vector<float> surface = start.start_m;
    const std::vector<uint8_t> domainsea = landscape::regionalsea(inputs, window, warp, parameters, surface);
    EXPECT_EQ(std::memcmp(surface.data(), start.start_m.data(), surface.size() * sizeof(float)), 0)
        << "regionalsea must not change the surface under the land-share rule";
    const std::vector<uint8_t> sea = [&] {
        std::vector<uint8_t> kept;
        for (int y = 0; y < window.fineheight(); y++)
            for (int x = 0; x < window.finewidth(); x++)
                kept.push_back(domainsea[at(window.keptfinecolumnoffset() + x, window.keptfinerowoffset() + y,
                                            window.domainfinewidth())]);
        return kept;
    }();
    const std::vector<float> kept = keptpart(start.start_m, window);
    const int k = window.scale, width = inputs.width, finewidth = window.finewidth(), fineheight = window.fineheight();
    ShareCheck check;
    const double share = parameters.coast.landshare;
    for (int r = 0; r < window.rows; r++)
        for (int c = 0; c < window.columns; c++)
        {
            const std::size_t global = at((window.column0 + c) % width, window.row0 + r, width);
            int land = 0;
            for (int y = r * k; y < (r + 1) * k; y++)
                for (int x = c * k; x < (c + 1) * k; x++)
                    land += sea[at(x, y, finewidth)] == 0 ? 1 : 0;
            const double landshare = static_cast<double>(land) / (k * k);
            check.blocks++;
            const bool seacell = inputs.network.sea[global] != 0;
            if (landshare > 0.0 && landshare < 1.0)
                check.coastal++;
            if (seacell)
            {
                check.mostseashare = std::max(check.mostseashare, landshare);
                check.outside += landshare > 1.0 - share + 1.0e-12 ? 1 : 0;
            }
            else
            {
                check.leastlandshare = std::min(check.leastlandshare, landshare);
                check.outside += landshare < share - 1.0e-12 ? 1 : 0;
            }
        }
    for (int y = 0; y < fineheight; y++)
        for (int x = 0; x < finewidth; x++)
        {
            const std::size_t cell = at(x, y, finewidth);
            const int i = (window.column0 + x / k) % width, j = window.row0 + y / k;
            const std::size_t global = at(i, j, width);
            bool nearsea = false;
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++)
                    if (j + dj >= 0 && j + dj < inputs.height)
                        nearsea = nearsea || inputs.network.sea[at((i + di + width) % width, j + dj, width)] != 0;
            if (sea[cell] != 0 && !nearsea)
                check.seawithoutsea++;
            if (sea[cell] != 0 && inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0)
                check.seainlake++;
            if (sea[cell] != 0 && kept[cell] > inputs.sealevel)
                check.inconsistent++;
            if (sea[cell] == 0 && nearsea && inputs.lakesurface_m[global] <= 0.0f && kept[cell] <= inputs.sealevel)
                check.inconsistent++;
            if (sea[cell] != 0 && x > 0 && y > 0 && x < finewidth - 1 && y < fineheight - 1 &&
                sea[at(x - 1, y, finewidth)] == 0 && sea[at(x + 1, y, finewidth)] == 0 && sea[at(x, y - 1, finewidth)] == 0 &&
                sea[at(x, y + 1, finewidth)] == 0)
                check.holes++;
            if (x + 1 < finewidth && y + 1 < fineheight && sea[cell] == sea[at(x + 1, y + 1, finewidth)] &&
                sea[at(x + 1, y, finewidth)] == sea[at(x, y + 1, finewidth)] && sea[cell] != sea[at(x + 1, y, finewidth)])
                check.pinches++;
        }
    std::printf("[regional] land share, k %d: %d blocks (%d coastal), %d outside the rule; least land share of a land "
                "cell %.3f, most of a sea cell %.3f; block-mean error %.4f m\n",
                k, check.blocks, check.coastal, check.outside, check.leastlandshare, check.mostseashare,
                start.blockmeanerror_m);
    return check;
}
}

TEST_F(Landscape, RegionalSurfaceWarpedUpsample)
{
    const RegionalGlobalInputs inputs = syntheticworld(32);
    const RegionalWindow window = windowa(4);  // the domain crosses the date line
    const std::vector<float> plain = landscape::upsample(inputs.elevation_m, 32, 16, window);
    const std::size_t count = plain.size();

    // A zero warp gives upsample exactly.
    FineWarp zero;
    zero.dx.assign(count, 0.0);
    zero.dy.assign(count, 0.0);
    const std::vector<float> same = landscape::upsample(inputs.elevation_m, 32, 16, window, zero);
    ASSERT_EQ(same.size(), count);
    EXPECT_EQ(std::memcmp(same.data(), plain.data(), count * sizeof(float)), 0);

    // A displacement of one global cell eastwards samples the fine cell k columns east,
    // across the date line too (column W - 1 to column 0).
    FineWarp east;
    east.dx.assign(count, 1.0);
    east.dy.assign(count, 0.0);
    const std::vector<float> shifted = landscape::upsample(inputs.elevation_m, 32, 16, window, east);
    const int domainwidth = window.domainfinewidth();
    int crossings = 0;
    double largest = 0.0;
    for (int y = 0; y < window.domainfineheight(); y++)
        for (int x = 0; x + window.scale < domainwidth; x++)
        {
            largest = std::max(largest, std::fabs(static_cast<double>(shifted[at(x, y, domainwidth)]) -
                                                  plain[at(x + window.scale, y, domainwidth)]));
            if (window.globalcolumn(x) == 31)
                crossings++;
        }
    EXPECT_GT(crossings, 0);
    EXPECT_LT(largest, 1.0e-3);

    // A warp of the wrong size throws.
    FineWarp wrong;
    wrong.dx.assign(count - 1, 0.0);
    wrong.dy.assign(count - 1, 0.0);
    EXPECT_THROW(landscape::upsample(inputs.elevation_m, 32, 16, window, wrong), std::invalid_argument);
}

TEST_F(Landscape, RegionalSurfaceWarpBounds)
{
    const RegionalWindow window = windowc(8);
    RegionalParameters parameters = coastparameters();
    parameters.warp.amplitude_cells = 0.0;
    EXPECT_TRUE(landscape::regionalwarp(window, 1, parameters).empty());

    parameters.warp.amplitude_cells = 0.6;
    const FineWarp warp = landscape::regionalwarp(window, 1, parameters);
    ASSERT_EQ(warp.dx.size(), static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight());
    double longest = 0.0, total = 0.0;
    for (std::size_t cell = 0; cell < warp.dx.size(); cell++)
    {
        const double length = std::hypot(warp.dx[cell], warp.dy[cell]);
        longest = std::max(longest, length);
        total += length;
    }
    std::printf("[regional] warp at amplitude 0.6: longest %.4f, mean %.4f global cells\n", longest,
                total / static_cast<double>(warp.dx.size()));
    EXPECT_LT(longest, 0.6);
    EXPECT_GT(longest, 0.3);
    // Another seed gives another warp.
    const FineWarp other = landscape::regionalwarp(window, 2, parameters);
    EXPECT_NE(std::memcmp(other.dx.data(), warp.dx.data(), warp.dx.size() * sizeof(double)), 0);
}

TEST_F(Landscape, RegionalSurfaceNoise)
{
    // The exported noise: a pure function of the absolute fine cell, standard deviation
    // about 1, inactive (0) below k = 2.
    const landscape::RegionalNoise noise(64, 32, 16, 1, 0xc0a5, 0.5);
    ASSERT_TRUE(noise.active());
    double sum = 0.0, squares = 0.0, largest = 0.0;
    int count = 0;
    for (std::int64_t j = 0; j < 32 * 16; j += 3)
        for (std::int64_t i = 0; i < 64 * 16; i += 5)
        {
            const double value = noise.at(i, j);
            const double weight = std::cos(pi / 2.0 - pi * (static_cast<double>(j) + 0.5) / (32 * 16));
            sum += weight * value;
            squares += weight * value * value;
            largest = std::max(largest, std::fabs(value));
            count++;
        }
    double weights = 0.0;
    for (std::int64_t j = 0; j < 32 * 16; j += 3)
        weights += (64 * 16 + 4) / 5 * std::cos(pi / 2.0 - pi * (static_cast<double>(j) + 0.5) / (32 * 16));
    const double mean = sum / weights, deviation = std::sqrt(squares / weights - mean * mean);
    std::printf("[regional] coast noise: %d cells, mean %.4f, deviation %.4f, largest %.4f\n", count, mean, deviation,
                largest);
    EXPECT_NEAR(deviation, 1.0, 0.2);
    EXPECT_LT(largest, 7.0);
    EXPECT_EQ(noise.at(100, 200), noise.at(100 + 64 * 16, 200));
    EXPECT_FALSE(landscape::RegionalNoise(64, 32, 1, 1, 0xc0a5, 0.5).active());
}

TEST_F(Landscape, RegionalSurfaceWindowsAgree)
{
    // Two windows give the same warp on their shared domain cells and the same start
    // surface, base and sea on their kept overlap, bit for bit.
    const RegionalParameters parameters = coastparameters();
    const RegionalGlobalInputs inputs = syntheticworld(64);
    for (int scale : {4, 8})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow wc = windowc(scale), wd = windowd(scale);
        const FineWarp warpc = landscape::regionalwarp(wc, inputs.seed, parameters);
        const FineWarp warpd = landscape::regionalwarp(wd, inputs.seed, parameters);
        const StartSurface startc = landscape::regionalstartsurface(inputs, wc, warpc, parameters);
        const StartSurface startd = landscape::regionalstartsurface(inputs, wd, warpd, parameters);
        std::vector<float> surfacec = startc.start_m, surfaced = startd.start_m;
        const std::vector<uint8_t> seac = landscape::regionalsea(inputs, wc, warpc, parameters, surfacec);
        const std::vector<uint8_t> sead = landscape::regionalsea(inputs, wd, warpd, parameters, surfaced);
        const std::int64_t planetwidth = 64LL * scale;
        int shareddomain = 0, sharedkept = 0, warps = 0, surfaces = 0, seas = 0, sea = 0;
        for (int y = 0; y < wc.domainfineheight(); y++)
            for (int x = 0; x < wc.domainfinewidth(); x++)
            {
                const std::int64_t column = wc.absolutefinecolumn(x), row = wc.absolutefinerow(y);
                const int xd = static_cast<int>(((column - static_cast<std::int64_t>(wd.domaincolumn0) * scale) % planetwidth +
                                                 planetwidth) % planetwidth);
                const int yd = static_cast<int>(row - static_cast<std::int64_t>(wd.domainrow0) * scale);
                if (xd >= wd.domainfinewidth() || yd < 0 || yd >= wd.domainfineheight())
                    continue;
                shareddomain++;
                const std::size_t a = at(x, y, wc.domainfinewidth()), b = at(xd, yd, wd.domainfinewidth());
                warps += warpc.dx[a] != warpd.dx[b] || warpc.dy[a] != warpd.dy[b] ? 1 : 0;
                const int xk = x - wc.keptfinecolumnoffset(), yk = y - wc.keptfinerowoffset();
                const int xdk = xd - wd.keptfinecolumnoffset(), ydk = yd - wd.keptfinerowoffset();
                if (xk < 0 || yk < 0 || xk >= wc.finewidth() || yk >= wc.fineheight() || xdk < 0 || ydk < 0 ||
                    xdk >= wd.finewidth() || ydk >= wd.fineheight())
                    continue;
                sharedkept++;
                sea += seac[a];
                surfaces += std::memcmp(&surfacec[a], &surfaced[b], sizeof(float)) != 0 ||
                                    std::memcmp(&startc.base_m[a], &startd.base_m[b], sizeof(float)) != 0 ||
                                    std::memcmp(&startc.lapsebase_m[a], &startd.lapsebase_m[b], sizeof(float)) != 0
                                ? 1
                                : 0;
                seas += seac[a] != sead[b] ? 1 : 0;
            }
        std::printf("[regional] windows C/D, k %d: %d shared domain cells, %d kept (%d sea); differing: warp %d, "
                    "surface %d, sea %d\n",
                    scale, shareddomain, sharedkept, sea, warps, surfaces, seas);
        EXPECT_EQ(sharedkept, 12 * 4 * scale * scale);
        EXPECT_GT(sea, 0);
        EXPECT_EQ(warps, 0);
        EXPECT_EQ(surfaces, 0);
        EXPECT_EQ(seas, 0);
    }
}

TEST_F(Landscape, RegionalSurfaceLandShare)
{
    // Every global cell's block keeps the land share within the rule, sea only where
    // Phase E allowed it, no one-cell holes, the surface on the right side of the sea level.
    const RegionalParameters parameters = coastparameters();
    const RegionalGlobalInputs small = syntheticworld(32);
    const RegionalGlobalInputs large = syntheticworld(64);
    for (int scale : {4, 8, 16})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        for (const auto& [inputs, window] : {std::make_pair(&small, windowa(scale)), std::make_pair(&large, windowc(scale))})
        {
            const ShareCheck check = checkshares(*inputs, window, parameters);
            EXPECT_GT(check.coastal, 0);
            EXPECT_EQ(check.outside, 0);
            EXPECT_EQ(check.seawithoutsea, 0);
            EXPECT_EQ(check.seainlake, 0);
            EXPECT_EQ(check.holes, 0);
            // Pinches: one sweep before the hole fill; the sweep can leave a pinch of its own making
            // (erosion run 9, final sea masks: 1, 1 and 1 on windows A, B and C), so a small bound.
            std::printf("[regional] diagonal pinches left: %d\n", check.pinches);
            EXPECT_LE(check.pinches, 2);
            EXPECT_EQ(check.inconsistent, 0);
        }
    }
}

TEST_F(Landscape, RegionalOverlapCoast)
{
    // The pass with the R2 features on is bit-identical on the overlap of two windows, in
    // columns across the date line and in rows clear of the poles, at the smallest halo.
    const RegionalParameters parameters = coastparameters();
    const RegionalGlobalInputs columns = syntheticworld(32);
    const RegionalGlobalInputs rows = syntheticworld(64);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts first = compareoverlap(columns, windowa(scale), windowb(scale), "coast A/B", parameters);
        EXPECT_EQ(first.shared, 4 * 5 * scale * scale);
        const OverlapCounts second = compareoverlap(rows, windowc(scale), windowd(scale), "coast C/D", parameters);
        EXPECT_EQ(second.shared, 12 * 4 * scale * scale);
        EXPECT_GT(second.sea, 0);
        EXPECT_GT(second.lake, 0);
        EXPECT_GT(second.channel, 0);
    }
}

TEST_F(Landscape, RegionalOverlapCoastLowRelief)
{
    // The low-relief variant (elevations scaled by 0.01) at k = 8: the coast noise and the
    // land share act over plains of a few metres.
    const RegionalParameters parameters = coastparameters();
    const RegionalGlobalInputs rows = syntheticworld(64, 0.01);
    const OverlapCounts counts = compareoverlap(rows, windowc(8), windowd(8), "coast low relief C/D", parameters);
    EXPECT_GT(counts.sea, 0);
    EXPECT_GT(counts.channel, 0);
    const ShareCheck check = checkshares(rows, windowc(8), parameters);
    EXPECT_EQ(check.outside, 0);
    EXPECT_EQ(check.holes, 0);
}

namespace
{
// The quantile q of |cell - mean of its 8 neighbours| over the kept window's land fine cells
// of coastal (3 x 3 holds land and sea) or inland global cells, on the start surface.
double landroughness(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                     bool coastal, double q)
{
    const FineWarp warp = landscape::regionalwarp(window, inputs.seed, parameters);
    const StartSurface start = landscape::regionalstartsurface(inputs, window, warp, parameters);
    const std::vector<float> kept = keptpart(start.start_m, window);
    const int k = window.scale, width = inputs.width, w = window.finewidth(), h = window.fineheight();
    std::vector<double> values;
    for (int y = 1; y < h - 1; y++)
        for (int x = 1; x < w - 1; x++)
        {
            if (kept[at(x, y, w)] <= inputs.sealevel)
                continue;
            const int i = (window.column0 + x / k) % width, j = window.row0 + y / k;
            bool sea = false, land = false;
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++)
                    if (j + dj >= 0 && j + dj < inputs.height)
                        (inputs.network.sea[at((i + di + width) % width, j + dj, width)] != 0 ? sea : land) = true;
            if ((sea && land) != coastal)
                continue;
            // Land cells with land all round: the step down to the sea (the steep coast) is
            // not roughness.
            double mean = 0.0;
            bool shore = false;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++)
                    if (dx != 0 || dy != 0)
                    {
                        mean += kept[at(x + dx, y + dy, w)] / 8.0;
                        shore = shore || kept[at(x + dx, y + dy, w)] <= inputs.sealevel;
                    }
            if (!shore)
                values.push_back(std::fabs(kept[at(x, y, w)] - mean));
        }
    if (values.empty())
        return 0.0;
    std::sort(values.begin(), values.end());
    return values[static_cast<std::size_t>(q * static_cast<double>(values.size() - 1))];
}
}

TEST_F(Landscape, RegionalSurfaceCoastNoiseRoughness)
{
    // The coast noise sets the coastline on the decision surface; with the land-share rule on
    // the elevation takes none of it (only the shore, capped at coast.shorefloor_m), so
    // coastal land away from the shore stays about as rough as without it, and as rough as
    // inland land, even at a strong amplitude.
    const RegionalGlobalInputs inputs = syntheticworld(64);
    RegionalParameters strong = coastparameters();
    strong.coast.noiseamplitude_m = 300.0;
    RegionalParameters none = strong;
    none.coast.noiseamplitude_m = 0.0;
    for (int scale : {8, 16})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow window = windowc(scale);
        const double coastal = landroughness(inputs, window, strong, true, 0.99);
        const double plain = landroughness(inputs, window, none, true, 0.99);
        const double inland = landroughness(inputs, window, strong, false, 0.99);
        std::printf("[regional] land roughness q99, k %d: coastal %.2f m (%.2f m without the coast noise), inland %.2f m\n",
                    scale, coastal, plain, inland);
        EXPECT_LE(coastal, 2.0 * plain + 1.0);
        EXPECT_LE(coastal, 3.0 * inland + 1.0);
    }
}

TEST_F(Landscape, RegionalOverlapCoastNoShare)
{
    // The branch with the warp and the coast noise on and the land-share rule off (the
    // centre rule): bit-identical on the overlaps at the smallest halo.
    RegionalParameters parameters = coastparameters();
    parameters.coast.landshare = 0.0;
    parameters.coast.seadrop_m = 0.0;
    const RegionalGlobalInputs columns = syntheticworld(32);
    const RegionalGlobalInputs rows = syntheticworld(64);
    for (int scale : {4, 5})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const OverlapCounts first = compareoverlap(columns, windowa(scale), windowb(scale), "coast no share A/B", parameters);
        EXPECT_EQ(first.shared, 4 * 5 * scale * scale);
        const OverlapCounts second = compareoverlap(rows, windowc(scale), windowd(scale), "coast no share C/D", parameters);
        EXPECT_GT(second.sea, 0);
    }
}

TEST_F(Landscape, RegionalOverlapCoastWideWorld)
{
    // On a 256 x 128 world the land-share box (the domain widened by its margin) is neither
    // a whole turn nor at a pole, so its edge is exercised: two windows across the western
    // coast of the continent (about 132 E), overlapping in 130 .. 140 E, bit-identical at
    // the smallest halo, at k = 4 (share reach 1 global cell per pass) and k = 2 (2 cells).
    const RegionalParameters parameters = coastparameters();
    const RegionalGlobalInputs inputs = syntheticworld(256);
    for (int scale : {4, 2})
    {
        SCOPED_TRACE("k = " + std::to_string(scale));
        const RegionalWindow wa = landscape::regionalwindow(256, scale, testhalo, 120.0, -5.0, 140.0, 15.0);
        const RegionalWindow wb = landscape::regionalwindow(256, scale, testhalo, 130.0, 0.0, 150.0, 20.0);
        const int margin = (scale <= 2 ? 2 : 1) * parameters.coast.sharepasses + 2 + 2 * parameters.coast.signedpasses;
        for (const RegionalWindow& window : {wa, wb})
        {
            ASSERT_LT(window.domaincolumns + 2 * margin, 256);
            ASSERT_GT(window.domainrow0 - margin, 0);
            ASSERT_LT(window.domainrow0 + window.domainrows + margin, 128);
        }
        const auto start = std::chrono::steady_clock::now();
        const OverlapCounts counts = compareoverlap(inputs, wa, wb, "coast wide A/B", parameters);
        std::printf("[regional] wide world, k %d: both windows in %.2f s\n", scale,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        EXPECT_GT(counts.sea, 0);
        EXPECT_LT(counts.sea, counts.shared);
    }
}

// Multiscale plan S3.4: the share of a level's octave band in the noise (persistence p: octaves
// 0 .. log2 k weigh p^o; level 1 adds octaves 1 and 2 at k 4, level 2 octaves 3 and 4 at k 16).
TEST_F(Landscape, RegionalSurfaceLevelNoiseShare)
{
    const double p = 0.6;
    const auto squares = [p](int first, int last) {
        double sum = 0.0;
        for (int o = first; o <= last; o++)
            sum += std::pow(p, 2 * o);
        return sum;
    };
    EXPECT_NEAR(landscape::levelnoiseshare(4, 1, p), std::sqrt(squares(1, 2) / squares(0, 2)), 1e-12);
    EXPECT_NEAR(landscape::levelnoiseshare(16, 4, p), std::sqrt(squares(3, 4) / squares(0, 4)), 1e-12);
    EXPECT_NEAR(landscape::levelnoiseshare(64, 16, p), std::sqrt(squares(5, 6) / squares(0, 6)), 1e-12);
    // A scale that is not a power of two: k 8 over 1 adds octaves 1 .. 3; the same scale adds none.
    EXPECT_NEAR(landscape::levelnoiseshare(8, 1, p), std::sqrt(squares(1, 3) / squares(0, 3)), 1e-12);
    EXPECT_EQ(landscape::levelnoiseshare(16, 16, p), 0.0);
    EXPECT_THROW(landscape::levelnoiseshare(4, 0, p), std::invalid_argument);
    EXPECT_THROW(landscape::levelnoiseshare(4, 8, p), std::invalid_argument);
}

// Multiscale plan S3.4: regional.noise_relief_field 1 takes Phase E's noise amplitude from
// inputs.reliefamplitude_m: without the field the noise is the local relief's bit for bit; with a
// zero field it is the noise of noise_relief_fraction 0 (the amplitude noise_amplitude_m alone), in
// the start surface and in a level's octave band.
TEST_F(Landscape, RegionalSurfaceNoiseReliefField)
{
    const RegionalGlobalInputs inputs = syntheticworld(64);
    const RegionalWindow window = windowc(16);
    RegionalParameters plain = coastparameters();
    RegionalParameters field = plain;
    field.noiserelieffield = 1;
    RegionalParameters flat = plain;
    flat.noiserelieffraction = 0.0;
    RegionalGlobalInputs zero = inputs;
    zero.reliefamplitude_m.assign(inputs.elevation_m.size(), 0.0f);
    const auto start = [&](const RegionalGlobalInputs& in, const RegionalParameters& parameters) {
        return landscape::regionalstartsurface(in, window, landscape::regionalwarp(window, in.seed, parameters), parameters).start_m;
    };
    const auto band = [&](const RegionalGlobalInputs& in, const RegionalParameters& parameters) {
        std::vector<float> surface(static_cast<std::size_t>(window.domainfinewidth()) * window.domainfineheight(), 100.0f);
        landscape::addlevelnoise(in, window, 4, parameters, surface);
        return surface;
    };
    const auto same = [](const std::vector<float>& a, const std::vector<float>& b) {
        return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
    };
    EXPECT_TRUE(same(start(inputs, field), start(inputs, plain)));
    EXPECT_TRUE(same(band(inputs, field), band(inputs, plain)));
    EXPECT_TRUE(same(start(zero, field), start(inputs, flat)));
    EXPECT_TRUE(same(band(zero, field), band(inputs, flat)));
    EXPECT_FALSE(same(band(zero, field), band(inputs, plain)));
}
