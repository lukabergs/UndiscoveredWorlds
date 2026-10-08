#include "simulations/landscape/regional_surface.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "core/deterministic_random.hpp"
#include "simulations/landscape/spherical_noise.hpp"

// Erosion plan steps E3 and R0 (docs/architecture/erosion-hydrology-plan.md, Phase E
// and Phase R): moved from regional_pass.cpp unchanged.

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
// Salt of the regional noise (deterministiccombine(seed, salt)).
constexpr std::uint64_t noisesalt = 0x5e6a;

std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

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
}

FineWarp regionalwarp(const RegionalWindow& window, std::uint64_t seed, const RegionalParameters& parameters)
{
    // R0: the contract only; R2 draws the warp.
    (void)window;
    (void)seed;
    (void)parameters;
    return FineWarp{};
}

StartSurface regionalstartsurface(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const FineWarp& warp, const RegionalParameters& parameters)
{
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const double sealevel = static_cast<double>(inputs.sealevel);

    // The base elevation upsampled at the fine centres, plus noise of the band from
    // twice the global north-south spacing (exclusive) down to twice the fine one,
    // with an amplitude that grows with the local global relief and tapers to zero
    // within coastband_m of the sea level; then the block means.
    StartSurface result;
    result.base_m = upsample(inputs.elevation_m, width, height, window, warp);
    result.lapsebase_m = warp.empty() ? result.base_m : upsample(inputs.elevation_m, width, height, window);
    const std::vector<float>& base = result.base_m;
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

    std::vector<float>& surface = result.start_m;
    surface = base;
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
        result.blockmeanerror_m = largest;
    }
    return result;
}

std::vector<uint8_t> regionalsea(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                 const FineWarp& warp, const RegionalParameters& parameters,
                                 std::vector<float>& surface_m)
{
    (void)warp;
    (void)parameters;
    const int width = inputs.width, height = inputs.height;
    const int k = window.scale;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t globalsize = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);

    // At or below the sea level, with a global sea cell among the global cell and its
    // 8 neighbours. Then each global cell's centre fine cell follows the global mask:
    // a sea centre above the sea level is lowered to the sea level, a land centre
    // counted as sea is raised to the next float above it. A fine cell of a global
    // lake cell (lake surface above zero) is never sea unless its global cell is sea:
    // a lake below the sea level beside the ocean stays a lake.
    std::vector<uint8_t> nearsea(globalsize, 0);
    for (int j = 0; j < height; j++)
        for (int i = 0; i < width; i++)
        {
            const std::size_t global = at(i, j, width);
            const bool lake = inputs.lakesurface_m[global] > 0.0f && inputs.network.sea[global] == 0;
            nearsea[global] = !lake && aroundcell(width, height, i, j, [&](std::size_t other) {
                return inputs.network.sea[other] != 0;
            }) ? 1 : 0;
        }
    const float abovesea = std::nextafter(inputs.sealevel, std::numeric_limits<float>::infinity());
    std::vector<uint8_t> sea(count, 0);
    for (int y = 0; y < domainheight; y++)
    {
        const int row = window.globalrow(y);
        const bool centrerow = window.absolutefinerow(y) % k == k / 2;
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = at(x, y, domainwidth);
            const std::size_t global = at(window.globalcolumn(x), row, width);
            if (centrerow && window.absolutefinecolumn(x) % k == k / 2)
            {
                if (inputs.network.sea[global] != 0)
                    surface_m[cell] = std::min(surface_m[cell], inputs.sealevel);
                else if (surface_m[cell] <= inputs.sealevel && nearsea[global] != 0)
                    surface_m[cell] = abovesea;
            }
            sea[cell] = surface_m[cell] <= inputs.sealevel && nearsea[global] != 0 ? 1 : 0;
        }
    }
    return sea;
}
}
