#include "simulations/landscape/regional_ice.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <cstdint>

#include "simulations/landscape/regional_surface.hpp"

// Erosion plan Phase E5, step E5.3 (docs/architecture/erosion-hydrology-plan.md): the
// present-day ice above the present climate's ELA, the contract of regional_ice.hpp.

namespace landscape
{
namespace
{
// The salt of the ice edge's noise (its own: R4's ice and lake noises use 0x1ce1 and 0x1a4f).
constexpr std::uint64_t elanoisesalt = 0xe1a3;
// The thinnest ice the rule places: a floor that keeps a cell just above the threshold
// visible as ice (thickness above zero), as regional_water.cpp's iceminimum_m; not a
// tuning value.
constexpr double iceminimum_m = 1.0;
}

void iceabovetheela(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                    const FineGlacial& glacial, const std::vector<uint8_t>& sea, const std::vector<float>& elevation_m,
                    std::vector<float>& icethickness_m)
{
    if (!(parameters.ice.thicknessfraction > 0.0) || !glacial.active)
        return;
    const int domainwidth = window.domainfinewidth(), domainheight = window.domainfineheight();
    const std::size_t count = static_cast<std::size_t>(domainwidth) * static_cast<std::size_t>(domainheight);
    if (glacial.presentela_m.size() != count || sea.size() != count || elevation_m.size() != count ||
        icethickness_m.size() != count)
        return;
    const double fraction = parameters.ice.thicknessfraction;
    const double thickest = parameters.ice.thicknessmax_m;
    const double edgenoise = parameters.ice.edgenoise_m;
    const RegionalNoise noise(inputs.width, inputs.height, window.scale, inputs.seed, elanoisesalt,
                              parameters.coast.noisepersistence);
    const bool noisy = edgenoise > 0.0 && noise.active();
    for (int y = 0; y < domainheight; y++)
    {
        const std::int64_t finerow = window.absolutefinerow(y);
        for (int x = 0; x < domainwidth; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(domainwidth) +
                                     static_cast<std::size_t>(x);
            if (sea[cell] != 0)
                continue;
            const double z = static_cast<double>(elevation_m[cell]);
            const double ela = static_cast<double>(glacial.presentela_m[cell]);
            if (!std::isfinite(z) || !std::isfinite(ela))
                continue;
            const double threshold = noisy ? ela + edgenoise * noise.at(window.absolutefinecolumn(x), finerow) : ela;
            if (!(z >= threshold))
                continue;
            const double thickness = std::max(iceminimum_m, std::min(thickest, fraction * (z - ela)));
            icethickness_m[cell] = std::max(icethickness_m[cell], static_cast<float>(thickness));
        }
    }
}

void iceabovetheela(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const RegionalParameters& parameters,
                    const FineGlacial& glacial, const std::vector<uint8_t>& sea, const std::vector<float>& elevation_m,
                    std::vector<float>& icethickness_m)
{
    if (!(parameters.ice.thicknessfraction > 0.0) || !glacial.active)
        return;
    const std::size_t count = window.cellcount();
    if (glacial.presentela_m.size() != count || sea.size() != count || elevation_m.size() != count ||
        icethickness_m.size() != count)
        return;
    if (inputs.grid == nullptr || inputs.grid->kind() != mesh::GridKind::cube || window.kind() != mesh::GridKind::cube ||
        window.facecells() % inputs.grid->facecells() != 0)
        throw std::invalid_argument("iceabovetheela: a cube window refining the inputs' cube level 0");
    const int k = window.facecells() / inputs.grid->facecells();
    const int n0 = inputs.grid->facecells();
    const double fraction = parameters.ice.thicknessfraction;
    const double thickest = parameters.ice.thicknessmax_m;
    const double edgenoise = parameters.ice.edgenoise_m;
    const RegionalNoise noise(4 * n0, 2 * n0, k, inputs.seed, elanoisesalt, parameters.coast.noisepersistence);
    const bool noisy = edgenoise > 0.0 && noise.active();
    for (mesh::CellId cell = 0; cell < count; cell++)
    {
        if (sea[cell] != 0)
            continue;
        const double z = static_cast<double>(elevation_m[cell]);
        const double ela = static_cast<double>(glacial.presentela_m[cell]);
        if (!std::isfinite(z) || !std::isfinite(ela))
            continue;
        const mesh::Vec3 p = window.centre(cell);
        const double threshold = noisy ? ela + edgenoise * noise.atpoint(p.x, p.y, p.z) : ela;
        if (!(z >= threshold))
            continue;
        const double thickness = std::max(iceminimum_m, std::min(thickest, fraction * (z - ela)));
        icethickness_m[cell] = std::max(icethickness_m[cell], static_cast<float>(thickness));
    }
}
}
