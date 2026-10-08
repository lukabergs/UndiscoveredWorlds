// Tests of landscape/regional_pass.hpp (erosion plan Phase E, step E3), suite Landscape:
// TEST_F(Landscape, Regional<What>) on a synthetic world built here (a smooth
// continent across the date line with a ridge, a sea, a network routed by steepest
// descent on the Priority-Flood surface, a lake on one cell, climate varying with
// latitude), at scales 4, 5 and 8, default parameters except the halo (`testhalo`; the
// pass requires at least 2 x meanpasses + 2, RegionalPassRejectsASmallHalo).

#include "landscape_fixture.hpp"

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/upsampling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace
{
using landscape::RegionalFields;
using landscape::RegionalGlobalInputs;
using landscape::RegionalParameters;
using landscape::RegionalWindow;

constexpr double pi = 3.14159265358979323846;
constexpr double radius_m = 6371.0e3;
constexpr double secondsperyear = 365.25 * 86400.0;
constexpr float sealevel = 0.0f;
// The halo of these tests: two block-mean passes reach 4 global cells, so 6 leaves
// two cells for the channels and the solve (the plan's default is 8).
constexpr int testhalo = 6;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

double latitude_rad(int j, int height) { return pi / 2.0 - pi * (j + 0.5) / height; }
double longitude_rad(int i, int width) { return -pi + 2.0 * pi * (i + 0.5) / width; }

// The synthetic world on a W x W/2 grid: a continent centred at longitude 180, latitude
// 5 N (it crosses the date line), elevation falling linearly with the angular distance
// from its centre, a ridge whose latitude wanders with longitude, the sea beyond; a
// one-cell pit that the routing fills into a lake; climate from latitude and height.
RegionalGlobalInputs syntheticworld(int width, double relief = 1.0)
{
    const int height = width / 2;
    const std::size_t size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    RegionalGlobalInputs inputs;
    inputs.width = width;
    inputs.height = height;
    inputs.radius_m = radius_m;
    inputs.sealevel = sealevel;
    inputs.seed = 1;
    inputs.lapserate_c_per_km = 6.5f;
    inputs.elevation_m.resize(size);

    const double centrelatitude = 5.0 * pi / 180.0, centrelongitude = pi;
    const double reach = 0.85;  // the continent's angular radius
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const double phi = latitude_rad(j, height), lambda = longitude_rad(i, width);
            const double cosd = std::sin(phi) * std::sin(centrelatitude) +
                                std::cos(phi) * std::cos(centrelatitude) * std::cos(lambda - centrelongitude);
            const double d = std::acos(std::clamp(cosd, -1.0, 1.0));
            double e = 2400.0 * (reach - d) / reach - 50.0;
            if (d < reach)
            {
                const double ridgelatitude = (15.0 + 8.0 * std::sin(2.0 * lambda)) * pi / 180.0;
                const double offset = (phi - ridgelatitude) / (6.0 * pi / 180.0);
                e += 1500.0 * std::exp(-offset * offset) * (1.0 - d / reach);
            }
            else
                e = std::max(e, -4000.0);
            inputs.elevation_m[at(i, j, width)] = static_cast<float>(e);
        }

    // The lake: a pit 60 m below the lowest neighbour of the land cell at longitude
    // 170 E, latitude 10 S.
    const int lakei = static_cast<int>((170.0 + 180.0) / 360.0 * width);
    const int lakej = static_cast<int>((90.0 + 10.0) / 180.0 * height);
    {
        float lowest = inputs.elevation_m[at(lakei, lakej, width)];
        for (int dj = -1; dj <= 1; dj++)
            for (int di = -1; di <= 1; di++)
                if (di != 0 || dj != 0)
                    lowest = std::min(lowest, inputs.elevation_m[at((lakei + di + width) % width, lakej + dj, width)]);
        inputs.elevation_m[at(lakei, lakej, width)] = std::max(1.0f, lowest - 60.0f);
    }

    // The low-relief variant: every elevation scaled by `relief` about the sea level (0).
    for (float& e : inputs.elevation_m)
        e = static_cast<float>(relief * e);

    landscape::GlobalNetwork& network = inputs.network;
    network.width = width;
    network.height = height;
    network.sea.resize(size);
    network.direction.assign(size, 0);
    network.filledelevation_m.resize(size);
    network.drainagearea_km2.assign(size, 0.0f);
    network.discharge_m3s.assign(size, 0.0f);
    network.januaryflow.assign(size, 0);
    network.julyflow.assign(size, 0);
    for (std::size_t cell = 0; cell < size; cell++)
        network.sea[cell] = inputs.elevation_m[cell] <= sealevel ? 1 : 0;

    // Priority-Flood from the sea (ties by index): the filled surface and, for every
    // land cell, the direction towards the cell the flood came from.
    auto neighbourof = [&](int i, int j, int direction, int& ni, int& nj) {
        ni = (i + landscape::directionoffsetx(direction) + width) % width;
        nj = j + landscape::directionoffsety(direction);
        return nj >= 0 && nj < height;
    };
    std::vector<uint8_t> flooddirection(size, 0), done(size, 0);
    using entry = std::pair<float, std::size_t>;
    std::priority_queue<entry, std::vector<entry>, std::greater<entry>> queue;
    for (std::size_t cell = 0; cell < size; cell++)
        if (network.sea[cell] != 0)
        {
            network.filledelevation_m[cell] = inputs.elevation_m[cell];
            done[cell] = 1;
            queue.push({network.filledelevation_m[cell], cell});
        }
    while (!queue.empty())
    {
        const std::size_t cell = queue.top().second;
        queue.pop();
        const int i = static_cast<int>(cell % width), j = static_cast<int>(cell / width);
        for (int direction = 1; direction <= 8; direction++)
        {
            int ni = 0, nj = 0;
            if (!neighbourof(i, j, direction, ni, nj))
                continue;
            const std::size_t next = at(ni, nj, width);
            if (done[next] != 0)
                continue;
            done[next] = 1;
            network.filledelevation_m[next] = std::max(inputs.elevation_m[next], network.filledelevation_m[cell]);
            flooddirection[next] = static_cast<uint8_t>(landscape::oppositedirection(direction));
            queue.push({network.filledelevation_m[next], next});
        }
    }

    // Steepest descent on the filled surface where a lower neighbour exists, the
    // flood's direction otherwise (the lake's flat).
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            if (network.sea[cell] != 0)
                continue;
            const double dx = std::cos(latitude_rad(j, height));
            double steepest = 0.0;
            int best = 0;
            for (int direction = 1; direction <= 8; direction++)
            {
                int ni = 0, nj = 0;
                if (!neighbourof(i, j, direction, ni, nj))
                    continue;
                const double drop = static_cast<double>(network.filledelevation_m[cell]) -
                                    static_cast<double>(network.filledelevation_m[at(ni, nj, width)]);
                const int ox = landscape::directionoffsetx(direction), oy = landscape::directionoffsety(direction);
                const double distance = ox == 0 ? 1.0 : (oy == 0 ? dx : std::hypot(dx, 1.0));
                if (drop > 0.0 && drop / distance > steepest)
                {
                    steepest = drop / distance;
                    best = direction;
                }
            }
            network.direction[cell] = static_cast<uint8_t>(best != 0 ? best : flooddirection[cell]);
        }

    // Runoff, drainage area and discharge accumulated down the tree (Kahn's order).
    inputs.runoff_mm_yr.assign(size, 0.0f);
    std::vector<double> area(size, 0.0), water(size, 0.0);
    std::vector<int> donors(size, 0);
    std::vector<std::size_t> receiver(size, size);
    const double dlambda = 2.0 * pi / width;
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            if (network.sea[cell] != 0)
                continue;
            const double north = pi / 2.0 - pi * j / height, south = pi / 2.0 - pi * (j + 1) / height;
            area[cell] = radius_m * radius_m * dlambda * (std::sin(north) - std::sin(south)) / 1.0e6;
            const double c = std::cos(latitude_rad(j, height));
            inputs.runoff_mm_yr[cell] = static_cast<float>(200.0 + 600.0 * c * c);
            water[cell] = inputs.runoff_mm_yr[cell] / 1000.0 * area[cell] * 1.0e6;
            int ni = 0, nj = 0;
            neighbourof(i, j, network.direction[cell], ni, nj);
            receiver[cell] = at(ni, nj, width);
            if (network.sea[receiver[cell]] == 0)
                donors[receiver[cell]]++;
        }
    std::vector<std::size_t> ready;
    for (std::size_t cell = 0; cell < size; cell++)
        if (network.sea[cell] == 0 && donors[cell] == 0)
            ready.push_back(cell);
    while (!ready.empty())
    {
        const std::size_t cell = ready.back();
        ready.pop_back();
        network.drainagearea_km2[cell] = static_cast<float>(area[cell]);
        network.discharge_m3s[cell] = static_cast<float>(water[cell] / secondsperyear);
        network.januaryflow[cell] = static_cast<int32_t>(10.0 * network.discharge_m3s[cell]);
        network.julyflow[cell] = static_cast<int32_t>(5.0 * network.discharge_m3s[cell]);
        const std::size_t next = receiver[cell];
        if (network.sea[next] != 0)
            continue;
        area[next] += area[cell];
        water[next] += water[cell];
        if (--donors[next] == 0)
            ready.push_back(next);
    }

    inputs.lakesurface_m.assign(size, 0.0f);
    inputs.lakesurface_m[at(lakei, lakej, width)] = network.filledelevation_m[at(lakei, lakej, width)];

    inputs.januarytemperature_c.resize(size);
    inputs.julytemperature_c.resize(size);
    inputs.januaryprecipitation_mm.resize(size);
    inputs.julyprecipitation_mm.resize(size);
    inputs.icethickness_m.resize(size);
    inputs.climate.resize(size);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t cell = at(i, j, width);
            const double phi = latitude_rad(j, height);
            const double land = std::max(0.0, static_cast<double>(inputs.elevation_m[cell]));
            const double mean = 27.0 * std::cos(phi) - 12.0 - 6.5 * land / 1000.0;
            inputs.januarytemperature_c[cell] = static_cast<float>(mean - 10.0 * std::sin(phi));
            inputs.julytemperature_c[cell] = static_cast<float>(mean + 10.0 * std::sin(phi));
            const double c = std::cos(phi);
            inputs.januaryprecipitation_mm[cell] = static_cast<float>(40.0 + 120.0 * c * c);
            inputs.julyprecipitation_mm[cell] = static_cast<float>(30.0 + 150.0 * c * c * c);
            inputs.icethickness_m[cell] = static_cast<float>(std::max(0.0, 0.5 * (land - 2000.0)));
            inputs.climate[cell] = static_cast<int16_t>(std::fabs(phi) * 180.0 / pi / 15.0 + (network.sea[cell] != 0 ? 100 : 0));
        }
    return inputs;
}

// On the 32 x 16 world: window A crosses the date line (global columns 26 .. 35 mod 32,
// rows 3 .. 9), window B does not (columns 0 .. 7, rows 5 .. 11); they share columns
// 0 .. 3 and rows 5 .. 9 (4 x 5 global cells). Both domains reach a pole.
RegionalWindow windowa(int scale) { return landscape::regionalwindow(32, scale, testhalo, 112.5, -22.5, -135.0, 56.25); }
RegionalWindow windowb(int scale) { return landscape::regionalwindow(32, scale, testhalo, -180.0, -45.0, -90.0, 33.75); }

// On the 64 x 32 world, domains short of both poles: window C does not cross the date
// line (columns 40 .. 63, rows 11 .. 17; domain rows 5 .. 23), window D does (columns
// 52 .. 71 mod 64, rows 14 .. 20; domain rows 8 .. 26); they share columns 52 .. 63 and
// rows 14 .. 17 (12 x 4 global cells), which hold the coast, the lake (62, 17), ice
// and channels.
RegionalWindow windowc(int scale) { return landscape::regionalwindow(64, scale, testhalo, 45.0, -11.25, 180.0, 28.125); }
RegionalWindow windowd(int scale) { return landscape::regionalwindow(64, scale, testhalo, 112.5, -28.125, -135.0, 11.25); }

RegionalFields timedpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const char* name,
                         const RegionalParameters& parameters = RegionalParameters{})
{
    const auto start = std::chrono::steady_clock::now();
    RegionalFields fields = landscape::runregionalpass(inputs, window, parameters);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("[regional] %s: k %d, kept %d x %d global, domain %d x %d fine, %.3f s, block-mean error %.4f m, "
                "converged %s, %d receivers changed, lowering mean %.3f m max %.3f m\n",
                name, window.scale, window.columns, window.rows, window.domainfinewidth(), window.domainfineheight(),
                seconds, fields.blockmeanerror_m, fields.reliefconverged ? "yes" : "no", fields.changedreceivers,
                fields.meanlowering_m, fields.maxlowering_m);
    return fields;
}

// The largest difference of one float field on the overlap and where it lies.
struct FloatDifference
{
    std::string name;
    double largest = 0.0;
    int differing = 0;
    int edgea = 0, edgeb = 0;  // distance of the largest difference from each window's edge, fine cells
};

void comparefloat(FloatDifference& difference, float a, float b, int edgea, int edgeb)
{
    if (std::memcmp(&a, &b, sizeof(float)) == 0)
        return;
    difference.differing++;
    const double value = std::isfinite(a) && std::isfinite(b) ? std::fabs(static_cast<double>(a) - b) : 1.0e30;
    if (value >= difference.largest)
    {
        difference.largest = value;
        difference.edgea = edgea;
        difference.edgeb = edgeb;
    }
}

int edgedistance(int x, int y, int width, int height) { return std::min(std::min(x, width - 1 - x), std::min(y, height - 1 - y)); }

// What the overlap of two windows holds.
struct OverlapCounts
{
    int shared = 0, sea = 0, lake = 0, ice = 0, channel = 0;
};

// Runs both windows and compares them on their overlap: every code and mask must be
// identical; a float difference is printed with its largest value and where, and fails
// above 0.01 (m, C, mm, m3/s).
OverlapCounts compareoverlap(const RegionalGlobalInputs& inputs, const RegionalWindow& wa, const RegionalWindow& wb,
                             const std::string& label)
{
    const RegionalFields a = timedpass(inputs, wa, (label + " first").c_str());
    const RegionalFields b = timedpass(inputs, wb, (label + " second").c_str());
    const int scale = wa.scale;
    const std::int64_t planetwidth = static_cast<std::int64_t>(inputs.width) * scale;
    const int wa_width = wa.finewidth(), wa_height = wa.fineheight();
    const int wb_width = wb.finewidth(), wb_height = wb.fineheight();
    std::vector<FloatDifference> differences = {{"elevation_m"},          {"discharge_m3s"},     {"lakesurface_m"},
                                                {"januarytemperature_c"}, {"julytemperature_c"}, {"januaryprecipitation_mm"},
                                                {"julyprecipitation_mm"}, {"icethickness_m"}};
    OverlapCounts counts;
    int codes = 0;
    for (int ya = 0; ya < wa_height; ya++)
        for (int xa = 0; xa < wa_width; xa++)
        {
            const std::int64_t column = (static_cast<std::int64_t>(wa.column0) * scale + xa) % planetwidth;
            const std::int64_t row = static_cast<std::int64_t>(wa.row0) * scale + ya;
            const int xb = static_cast<int>(((column - static_cast<std::int64_t>(wb.column0) * scale) % planetwidth + planetwidth) % planetwidth);
            const int yb = static_cast<int>(row - static_cast<std::int64_t>(wb.row0) * scale);
            if (xb >= wb_width || yb < 0 || yb >= wb_height)
                continue;
            counts.shared++;
            const std::size_t ca = at(xa, ya, wa_width), cb = at(xb, yb, wb_width);
            counts.sea += a.sea[ca] != 0 ? 1 : 0;
            counts.lake += a.lakesurface_m[ca] > 0.0f ? 1 : 0;
            counts.ice += a.icethickness_m[ca] > 0.0f ? 1 : 0;
            counts.channel += a.riverdirection[ca] != 0 ? 1 : 0;
            if (a.sea[ca] != b.sea[cb] || a.riverdirection[ca] != b.riverdirection[cb] ||
                a.januaryflow[ca] != b.januaryflow[cb] || a.julyflow[ca] != b.julyflow[cb] || a.climate[ca] != b.climate[cb])
                codes++;
            const int ea = edgedistance(xa, ya, wa_width, wa_height), eb = edgedistance(xb, yb, wb_width, wb_height);
            comparefloat(differences[0], a.elevation_m[ca], b.elevation_m[cb], ea, eb);
            comparefloat(differences[1], a.discharge_m3s[ca], b.discharge_m3s[cb], ea, eb);
            comparefloat(differences[2], a.lakesurface_m[ca], b.lakesurface_m[cb], ea, eb);
            comparefloat(differences[3], a.januarytemperature_c[ca], b.januarytemperature_c[cb], ea, eb);
            comparefloat(differences[4], a.julytemperature_c[ca], b.julytemperature_c[cb], ea, eb);
            comparefloat(differences[5], a.januaryprecipitation_mm[ca], b.januaryprecipitation_mm[cb], ea, eb);
            comparefloat(differences[6], a.julyprecipitation_mm[ca], b.julyprecipitation_mm[cb], ea, eb);
            comparefloat(differences[7], a.icethickness_m[ca], b.icethickness_m[cb], ea, eb);
        }
    EXPECT_EQ(codes, 0) << label << ": cells whose sea, river direction, flows or climate class differ";
    std::printf("[regional] %s, k %d: %d shared fine cells: %d sea, %d lake, %d ice, %d channel; %d codes differ\n",
                label.c_str(), scale, counts.shared, counts.sea, counts.lake, counts.ice, counts.channel, codes);
    for (const FloatDifference& difference : differences)
    {
        if (difference.differing > 0)
            std::printf("[regional] %s, k %d: %s differs on %d cells, largest %.6g at %d and %d fine cells from the edges\n",
                        label.c_str(), scale, difference.name.c_str(), difference.differing, difference.largest,
                        difference.edgea, difference.edgeb);
        EXPECT_LE(difference.largest, 0.01) << label << ": " << difference.name;
        // The plan expects bit identity (section 6, 'Seams of the solve'); the tolerance
        // above is the plan row's failure bound, this the claim (E3's re-review).
        EXPECT_EQ(difference.differing, 0) << label << ": " << difference.name;
    }
    return counts;
}

// A cell of the kept window is the global cell's centre fine cell.
bool centrecell(const RegionalWindow& window, int x, int y)
{
    return x % window.scale == window.scale / 2 && y % window.scale == window.scale / 2;
}

std::size_t globalof(const RegionalWindow& window, int width, int x, int y)
{
    return at((window.column0 + x / window.scale) % width, window.row0 + y / window.scale, width);
}

// The kept window's part of a domain field.
std::vector<float> keptpart(const std::vector<float>& domain, const RegionalWindow& window)
{
    std::vector<float> kept;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
            kept.push_back(domain[at(window.keptfinecolumnoffset() + x, window.keptfinerowoffset() + y, window.domainfinewidth())]);
    return kept;
}

// Coast speckle: land fine cells in global sea blocks and sea fine cells in global land
// blocks that lie more than one fine cell from every block edge bordering a global cell
// of the other kind (Chebyshev distance to such a block above 2; no such block: always).
void speckle(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalFields& fields, int& landinsea,
             int& seainland)
{
    landinsea = seainland = 0;
    const int k = window.scale, width = inputs.width;
    for (int y = 0; y < window.fineheight(); y++)
        for (int x = 0; x < window.finewidth(); x++)
        {
            const std::size_t global = globalof(window, width, x, y);
            const uint8_t kind = inputs.network.sea[global];
            if (fields.sea[at(x, y, window.finewidth())] == kind)
                continue;
            const int ox = x % k, oy = y % k;
            const int i = (window.column0 + x / k) % width, j = window.row0 + y / k;
            int nearest = 1 << 30;
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++)
                {
                    if ((di == 0 && dj == 0) || j + dj < 0 || j + dj >= inputs.height)
                        continue;
                    if (inputs.network.sea[at((i + di + width) % width, j + dj, width)] == kind)
                        continue;
                    const int dx = di < 0 ? ox + 1 : (di > 0 ? k - ox : 0);
                    const int dy = dj < 0 ? oy + 1 : (dj > 0 ? k - oy : 0);
                    nearest = std::min(nearest, std::max(dx, dy));
                }
            if (nearest > 2)
                (kind != 0 ? landinsea : seainland)++;
        }
}
}

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

        // Lakes: the global surface where the fine bed lies below it, never sea.
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
