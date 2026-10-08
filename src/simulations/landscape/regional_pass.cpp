#include "simulations/landscape/regional_pass.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/spherical_noise.hpp"
#include "simulations/landscape/upsampling.hpp"

// Erosion plan step E3 (docs/architecture/erosion-hydrology-plan.md, Phase E, and
// section 6, "Decided for Phase E").

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double secondsperyear = 365.25 * 86400.0;
// Salt of the regional noise (deterministiccombine(seed, salt)).
constexpr std::uint64_t noisesalt = 0x5e6a;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

template<typename T>
void checksize(const std::vector<T>& field, std::size_t size, const char* name)
{
    if (field.size() != size)
        throw std::invalid_argument(std::string("runregionalpass: ") + name + " must hold width x height values");
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

// The local global relief: the largest absolute difference between a cell and its 8
// neighbours (columns wrap, rows beyond the poles do not exist).
std::vector<float> localrelief(const std::vector<float>& elevation, int width, int height)
{
    std::vector<float> relief(elevation.size(), 0.0f);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const float own = elevation[at(i, j, width)];
            float largest = 0.0f;
            aroundcell(width, height, i, j, [&](std::size_t other) {
                largest = std::max(largest, std::fabs(elevation[other] - own));
                return false;
            });
            relief[at(i, j, width)] = largest;
        }
    return relief;
}

// The domain's values at the kept window's cells.
template<typename T>
std::vector<T> cut(const std::vector<T>& domain, const RegionalWindow& window)
{
    const int width = window.finewidth(), height = window.fineheight();
    const int domainwidth = window.domainfinewidth();
    const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
    std::vector<T> kept(static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            kept[at(x, y, width)] = domain[at(x0 + x, y0 + y, domainwidth)];
    return kept;
}
}

RegionalFields runregionalpass(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                               const RegionalParameters& parameters)
{
    const int width = inputs.width, height = inputs.height;
    if (width < 1 || height < 1 || window.worldwidth != width || window.worldheight != height || window.scale < 1)
        throw std::invalid_argument("runregionalpass: the window belongs to another grid");
    if (!(inputs.radius_m > 0.0) || !std::isfinite(inputs.sealevel))
        throw std::invalid_argument("runregionalpass: the radius must be positive and the sea level finite");
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    checksize(inputs.elevation_m, globalsize, "elevation_m");
    checksize(inputs.network.sea, globalsize, "network.sea");
    checksize(inputs.runoff_mm_yr, globalsize, "runoff_mm_yr");
    checksize(inputs.lakesurface_m, globalsize, "lakesurface_m");
    checksize(inputs.januarytemperature_c, globalsize, "januarytemperature_c");
    checksize(inputs.julytemperature_c, globalsize, "julytemperature_c");
    checksize(inputs.januaryprecipitation_mm, globalsize, "januaryprecipitation_mm");
    checksize(inputs.julyprecipitation_mm, globalsize, "julyprecipitation_mm");
    checksize(inputs.icethickness_m, globalsize, "icethickness_m");
    checksize(inputs.climate, globalsize, "climate");

    // The block-mean correction reaches 2 global cells per pass beyond any fine cell it
    // changes; two more keep the channels and the solve's edge out of the kept window.
    if (parameters.meanpasses < 0 || window.halo < 2 * parameters.meanpasses + 2)
        throw std::invalid_argument("runregionalpass: the halo must be at least 2 x meanpasses + 2 global cells");

    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    const double sealevel = static_cast<double>(inputs.sealevel);
    const std::int64_t planetfinewidth = static_cast<std::int64_t>(width) * k;

    // 1. The start surface: the base elevation upsampled at the fine centres, plus
    // noise of the band from twice the global north-south spacing (exclusive) down to
    // twice the fine one, with an amplitude that grows with the local global relief
    // and tapers to zero within coastband_m of the sea level; then the block means.
    const std::vector<float> base = upsample(inputs.elevation_m, width, height, window);
    const std::vector<float> relief = upsample(localrelief(inputs.elevation_m, width, height), width, height, window);
    FractalNoise noise;
    noise.seed = deterministiccombine(inputs.seed, noisesalt);
    noise.radius_km = inputs.radius_m / 1000.0;
    noise.longest_km = 2.0 * pi * inputs.radius_m / height / 1000.0;
    noise.wavelength_km = noise.longest_km;
    // Twice the fine north-south spacing, written as longest / k so that the octave
    // whose wavelength equals it (k a power of two) is not lost to rounding.
    noise.shortest_km = noise.longest_km / k;
    noise.octaves = 1;
    while ((1 << noise.octaves) <= k)
        noise.octaves++;
    noise.persistence = parameters.noisepersistence;
    const FractalNoiseField field(noise);

    std::vector<float> surface = base;
    if (field.resolved() > 0)
    {
        for (int y = 0; y < domainheight; y++)
        {
            const double latitude = window.finelatitude_rad(y);
            for (int x = 0; x < domainwidth; x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                double px = 0.0, py = 0.0, pz = 0.0;
                unitvector(latitude, window.finelongitude_rad(x), px, py, pz);
                const double amplitude =
                    parameters.noiseamplitude_m +
                    parameters.noiserelieffraction * std::max(0.0, static_cast<double>(relief[cell]));
                const double height_m = std::fabs(static_cast<double>(base[cell]) - sealevel);
                // The taper scales with the local amplitude, so a strong noise cannot
                // lift a seabed or sink a lowland across the sea level far from the coast.
                const double band = parameters.coastband_m + amplitude;
                const double taper = band > 0.0 ? std::clamp(height_m / band, 0.0, 1.0) : 1.0;
                surface[cell] = static_cast<float>(static_cast<double>(base[cell]) +
                                                   taper * amplitude * field.sum(px, py, pz) / field.normalisation());
            }
        }
    }
    correctblockmeans(inputs.elevation_m, width, height, window, parameters.meanpasses, surface);

    RegionalFields fields;
    fields.window = window;
    {
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double largest = 0.0;
        for (int r = 0; r < window.rows; r++)
            for (int c = 0; c < window.columns; c++)
            {
                double sum = 0.0;
                for (int y = y0 + r * k; y < y0 + (r + 1) * k; y++)
                    for (int x = x0 + c * k; x < x0 + (c + 1) * k; x++)
                        sum += static_cast<double>(surface[at(x, y, domainwidth)]);
                const float global = inputs.elevation_m[at(wrapcolumn(window.column0 + c, width), window.row0 + r, width)];
                largest = std::max(largest, std::fabs(sum / (static_cast<double>(k) * k) - static_cast<double>(global)));
            }
        fields.blockmeanerror_m = largest;
    }

    // 2. The pinned channels take their bed (above the start surface or not).
    const FineChannels channels = rasterisechannels(window, inputs.network, parameters.channels, inputs.sealevel);
    for (std::size_t cell = 0; cell < count; cell++)
        if (channels.channel[cell] != 0)
            surface[cell] = channels.bed_m[cell];

    // 3. The sea by local rules: at or below the sea level, with a global sea cell
    // among the global cell and its 8 neighbours. Then each global cell's centre fine
    // cell follows the global mask: a sea centre above the sea level is lowered to the
    // sea level, a land centre counted as sea is raised to the next float above it.
    // A fine cell of a global lake cell (lake surface above zero) is never sea unless
    // its global cell is sea: a lake below the sea level beside the ocean stays a lake.
    std::vector<uint8_t> nearsea(globalsize, 0), nearice(globalsize, 0);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t global = at(i, j, width);
            const bool lake = inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
            nearsea[global] = !lake && aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.network.sea[other] != 0;
            }) ? 1 : 0;
            nearice[global] = aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.icethickness_m[other] > 0.0f;
            }) ? 1 : 0;
        }
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> sea(count, 0);
    std::vector<std::size_t> globalof(count);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        const bool centrerow = window.absolutefinerow(y) % k == k / 2;
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            globalof[cell] = global;
            if (centrerow && window.absolutefinecolumn(x) % k == k / 2)
            {
                if (inputs.network.sea[global] != 0)
                    surface[cell] = std::min(surface[cell], inputs.sealevel);
                else if (surface[cell] <= inputs.sealevel && nearsea[global] != 0)
                    surface[cell] = abovesea;
            }
            sea[cell] = surface[cell] <= inputs.sealevel && nearsea[global] != 0 ? 1 : 0;
        }
    }

    // 4. The analytical relief: channels, sea and the domain's outer ring pinned. A
    // sea cell enters the solve at the sea level (the base level of the rivers), so
    // the solve never lowers land below it through the sea; it keeps its own elevation.
    const PatchMetric metric = finepatchmetric(window, inputs.radius_m);
    const std::vector<float> runoff = upsample(inputs.runoff_mm_yr, width, height, window);
    AnalyticalReliefInput solveinput;
    solveinput.start_m = surface;
    solveinput.pinned.assign(count, 0);
    solveinput.runoff_m_yr.resize(count);
    solveinput.key.resize(count);
    for (int y = 0; y < domainheight; y++)
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const bool ring = x == 0 || y == 0 || x == domainwidth - 1 || y == domainheight - 1;
            solveinput.pinned[cell] = channels.channel[cell] != 0 || sea[cell] != 0 || ring ? 1 : 0;
            if (sea[cell] != 0)
                solveinput.start_m[cell] = inputs.sealevel;
            solveinput.runoff_m_yr[cell] = std::max(0.0f, runoff[cell]) / 1000.0f;
            solveinput.key[cell] = static_cast<std::uint64_t>(window.absolutefinecolumn(x) +
                                                              window.absolutefinerow(y) * planetfinewidth);
        }
    const AnalyticalReliefResult solved = solveanalyticalrelief(metric, solveinput, parameters.relief);
    fields.reliefconverged = solved.converged;
    fields.changedreceivers = solved.changedreceivers;
    std::vector<float> elevation = surface;
    for (std::size_t cell = 0; cell < count; cell++)
        if (solveinput.pinned[cell] == 0)
            elevation[cell] = solved.elevation_m[cell];
    {
        // The lowering by the solve over the kept window's unpinned (so land) cells.
        const int x0 = window.keptfinecolumnoffset(), y0 = window.keptfinerowoffset();
        double total = 0.0, largest = 0.0;
        std::size_t cells = 0;
        for (int y = y0; y < y0 + window.fineheight(); y++)
            for (int x = x0; x < x0 + window.finewidth(); x++)
            {
                const std::size_t cell = at(x, y, domainwidth);
                if (solveinput.pinned[cell] != 0)
                    continue;
                const double lowering = static_cast<double>(surface[cell]) - static_cast<double>(elevation[cell]);
                total += lowering;
                largest = std::max(largest, lowering);
                cells++;
            }
        fields.meanlowering_m = cells > 0 ? total / static_cast<double>(cells) : 0.0;
        fields.maxlowering_m = largest;
    }

    // 5.-7. Lakes, rivers and climate on the fine grid.
    const std::vector<float> januarytemperature = upsample(inputs.januarytemperature_c, width, height, window);
    const std::vector<float> julytemperature = upsample(inputs.julytemperature_c, width, height, window);
    const std::vector<float> januaryprecipitation = upsample(inputs.januaryprecipitation_mm, width, height, window);
    const std::vector<float> julyprecipitation = upsample(inputs.julyprecipitation_mm, width, height, window);
    const std::vector<float> ice = upsample(inputs.icethickness_m, width, height, window);
    const double lapse = static_cast<double>(inputs.lapserate_c_per_km) / 1000.0;

    std::vector<float> lakesurface(count, 0.0f), discharge(count, 0.0f);
    std::vector<uint8_t> riverdirection(count, 0);
    std::vector<int32_t> januaryflow(count, 0), julyflow(count, 0);
    std::vector<float> januarytemp(count), julytemp(count), januaryrain(count), julyrain(count), icethickness(count);
    std::vector<int16_t> climate(count);
    for (std::size_t cell = 0; cell < count; cell++)
    {
        const std::size_t global = globalof[cell];
        const float surfacelevel = inputs.lakesurface_m[global];
        if (sea[cell] == 0 && surfacelevel > 0.0f && elevation[cell] < surfacelevel)
            lakesurface[cell] = surfacelevel;
        if (channels.channel[cell] != 0)
        {
            riverdirection[cell] = channels.direction[cell];
            discharge[cell] = channels.discharge_m3s[cell];
            januaryflow[cell] = channels.januaryflow[cell];
            julyflow[cell] = channels.julyflow[cell];
        }
        else
            discharge[cell] = static_cast<float>(solved.discharge_m3_yr[cell] / secondsperyear);
        // Heights above the sea level only (the climate's rule, climate_elevation.hpp):
        // a sea cell whose interpolated elevation is below the sea level keeps the
        // interpolated temperature exactly.
        const double lift = std::max(0.0, static_cast<double>(base[cell]) - sealevel) -
                            std::max(0.0, static_cast<double>(elevation[cell]) - sealevel);
        januarytemp[cell] = static_cast<float>(static_cast<double>(januarytemperature[cell]) + lapse * lift);
        julytemp[cell] = static_cast<float>(static_cast<double>(julytemperature[cell]) + lapse * lift);
        januaryrain[cell] = std::max(0.0f, januaryprecipitation[cell]);
        julyrain[cell] = std::max(0.0f, julyprecipitation[cell]);
        icethickness[cell] = nearice[global] != 0 ? std::max(0.0f, ice[cell]) : 0.0f;
        // The global cell's class: classifyMonthly needs twelve months, the world
        // holds two seasons.
        climate[cell] = inputs.climate[global];
    }

    // 8. The kept window.
    fields.elevation_m = cut(elevation, window);
    fields.sea = cut(sea, window);
    fields.riverdirection = cut(riverdirection, window);
    fields.discharge_m3s = cut(discharge, window);
    fields.januaryflow = cut(januaryflow, window);
    fields.julyflow = cut(julyflow, window);
    fields.lakesurface_m = cut(lakesurface, window);
    fields.januarytemperature_c = cut(januarytemp, window);
    fields.julytemperature_c = cut(julytemp, window);
    fields.januaryprecipitation_mm = cut(januaryrain, window);
    fields.julyprecipitation_mm = cut(julyrain, window);
    fields.icethickness_m = cut(icethickness, window);
    fields.climate = cut(climate, window);
    return fields;
}
}
