#include "mineral_resources.hpp"
#include "physical_layers_internal.hpp"
#include "placer_transport.hpp"
#include <functional>

namespace physical_layers
{
using namespace detail;

namespace
{
int maxscoreinradius(const planet& world, int x, int y, int radius, const std::function<int(int, int)>& getter)
{
    const int width = world.width();
    const int height = world.height();
    int best = 0;

    for (int dy = -radius; dy <= radius; dy++)
    {
        const int ny = y + dy;

        if (ny < 0 || ny > height)
            continue;

        for (int dx = -radius; dx <= radius; dx++)
        {
            const int nx = wrap(x + dx, width);
            best = std::max(best, getter(nx, ny));
        }
    }

    return best;
}

float regimevolcanicbias(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::convergent_arc:
        return 1.0f;
    case GeologicRegime::divergent_rift:
        return 0.9f;
    case GeologicRegime::mid_ocean_ridge:
        return 0.85f;
    case GeologicRegime::continent_collision:
        return 0.3f;
    case GeologicRegime::trench_adjacent:
        return 0.25f;
    default:
        return 0.0f;
    }
}

float regimeorebias(GeologicRegime regime)
{
    switch (regime)
    {
    case GeologicRegime::continent_collision:
        return 1.0f;
    case GeologicRegime::convergent_arc:
        return 0.9f;
    case GeologicRegime::transform:
        return 0.55f;
    case GeologicRegime::divergent_rift:
        return 0.4f;
    case GeologicRegime::passive_margin:
        return 0.15f;
    default:
        return 0.0f;
    }
}

float boundaryvolcanicbias(BoundaryType boundarytype)
{
    switch (boundarytype)
    {
    case BoundaryType::convergent:
        return 1.0f;
    case BoundaryType::divergent:
        return 0.95f;
    case BoundaryType::transform:
        return 0.25f;
    case BoundaryType::passive_margin:
        return 0.12f;
    case BoundaryType::none:
    default:
        return 0.0f;
    }
}

float deformingvolcanicbias(DeformingRegionType regiontype)
{
    switch (regiontype)
    {
    case DeformingRegionType::continental_rift:
        return 0.85f;
    case DeformingRegionType::diffuse_collision:
        return 0.30f;
    case DeformingRegionType::none:
    default:
        return 0.0f;
    }
}

float crustvolcanicbias(CrustClass crustclass)
{
    switch (crustclass)
    {
    case CrustClass::transitional:
        return 1.0f;
    case CrustClass::oceanic:
        return 0.85f;
    case CrustClass::continental:
        return 0.75f;
    case CrustClass::none:
    default:
        return 0.0f;
    }
}

float boundaryorebias(BoundaryType boundarytype)
{
    switch (boundarytype)
    {
    case BoundaryType::convergent:
        return 1.0f;
    case BoundaryType::transform:
        return 0.80f;
    case BoundaryType::divergent:
        return 0.45f;
    case BoundaryType::passive_margin:
        return 0.25f;
    case BoundaryType::none:
    default:
        return 0.0f;
    }
}

float deformingorebias(DeformingRegionType regiontype)
{
    switch (regiontype)
    {
    case DeformingRegionType::diffuse_collision:
        return 1.0f;
    case DeformingRegionType::continental_rift:
        return 0.55f;
    case DeformingRegionType::none:
    default:
        return 0.0f;
    }
}

float crustorebias(CrustClass crustclass)
{
    switch (crustclass)
    {
    case CrustClass::continental:
        return 1.0f;
    case CrustClass::transitional:
        return 0.75f;
    case CrustClass::oceanic:
        return 0.35f;
    case CrustClass::none:
    default:
        return 0.0f;
    }
}

float boundaryproximityfactor(const planet& world, int x, int y, float maxdistance)
{
    if (maxdistance <= 0.0f)
        return 0.0f;

    return 1.0f - clampunit(static_cast<float>(world.tectonicboundarydistance(x, y)) / maxdistance);
}

void generate_volcanic_reserves(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const GeologicRegime regime = world.geologicregime(x, y);
                const BoundaryType boundarytype = world.tectonicboundarytype(x, y);
                const DeformingRegionType regiontype = world.tectonicdeformingregiontype(x, y);
                const CrustClass crustclass = world.tectoniccrustclass(x, y);
                const float regimemult = regimevolcanicbias(regime);
                const float divergencenorm = static_cast<float>(world.tectonicdivergence(x, y)) / 100.0f;
                const float convergencenorm = static_cast<float>(world.tectonicconvergence(x, y)) / 100.0f;
                const float uplift = clampunit(world.tectonicuplifttendency(x, y));
                const float subsidence = clampunit(world.tectonicsubsidencetendency(x, y));
                const float strain = clampunit(world.tectonicaccumulatedstrain(x, y));
                const float deformation = clampunit(world.tectonicdeformationrate(x, y));
                const float boundaryhistory = clampunit(world.tectonicboundaryhistory(x, y));
                const float boundaryproximity = boundaryproximityfactor(world, x, y, 10.0f);
                const float boundaryfocus = boundaryvolcanicbias(boundarytype) * boundaryproximity * (0.45f + 0.55f * boundaryhistory);
                const float regionbias = deformingvolcanicbias(regiontype);
                const float crustbias = crustvolcanicbias(crustclass);
                const float ridgebias = (world.oceanridges(x, y) != 0 || world.oceanrifts(x, y) != 0) ? 1.0f : 0.0f;
                const float volcanobias = (world.volcano(x, y) != 0 || world.strato(x, y)) ? 1.0f : 0.0f;
                const float nearbyvolcano = maxscoreinradius(world, x, y, 2, [&](int nx, int ny)
                {
                    return world.volcano(nx, ny) != 0 || world.strato(nx, ny) ? 100 : 0;
                }) / 100.0f;

                float score = 100.0f * (0.20f * regimemult + 0.18f * boundaryfocus + 0.12f * divergencenorm
                    + 0.08f * convergencenorm + 0.10f * uplift + 0.08f * strain + 0.06f * deformation
                    + 0.06f * regionbias + 0.06f * crustbias + 0.16f * volcanobias + 0.10f * nearbyvolcano);
                score += 12.0f * ridgebias;

                if (world.sea(x, y) != 0 && ridgebias == 0.0f && boundaryfocus < 0.35f && regimemult < 0.6f)
                    score *= 0.45f;

                if (subsidence > uplift && boundarytype == BoundaryType::passive_margin)
                    score *= 0.7f;

                world.setvolcanicreserve(x, y, clampscore(score));
            }
        }
    }, 16);
}

void generate_metal_ore_reserves(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const GeologicRegime regime = world.geologicregime(x, y);
                const BoundaryType boundarytype = world.tectonicboundarytype(x, y);
                const DeformingRegionType regiontype = world.tectonicdeformingregiontype(x, y);
                const CrustClass crustclass = world.tectoniccrustclass(x, y);
                const float regimebias = regimeorebias(regime);
                const float convergencenorm = static_cast<float>(world.tectonicconvergence(x, y)) / 100.0f;
                const float shearnorm = static_cast<float>(world.tectonicshear(x, y)) / 100.0f;
                const float uplift = clampunit(world.tectonicuplifttendency(x, y));
                const float strain = clampunit(world.tectonicaccumulatedstrain(x, y));
                const float deformation = clampunit(world.tectonicdeformationrate(x, y));
                const float boundaryhistory = clampunit(world.tectonicboundaryhistory(x, y));
                const float boundaryproximity = boundaryproximityfactor(world, x, y, 12.0f);
                const float boundaryfocus = boundaryorebias(boundarytype) * boundaryproximity * (0.40f + 0.60f * boundaryhistory);
                const float regionbias = deformingorebias(regiontype);
                const float crustbias = crustorebias(crustclass);
                const float volcanismnear = static_cast<float>(maxscoreinradius(world, x, y, 2, [&](int nx, int ny)
                {
                    return world.volcanicreserve(nx, ny);
                })) / 100.0f;
                const float mountainbias = clampunit(static_cast<float>(world.mountainheight(x, y)) / 5000.0f);

                float score = 100.0f * (0.18f * regimebias + 0.18f * boundaryfocus + 0.15f * convergencenorm
                    + 0.12f * shearnorm + 0.12f * strain + 0.08f * uplift + 0.05f * deformation
                    + 0.05f * regionbias + 0.05f * crustbias + 0.10f * mountainbias + 0.10f * volcanismnear);

                if (boundarytype == BoundaryType::transform)
                    score += 8.0f * boundaryhistory;

                if (regiontype == DeformingRegionType::diffuse_collision)
                    score += 12.0f * (0.5f + 0.5f * boundaryhistory);

                if (world.sea(x, y) != 0 && world.coast(x, y) == 0 && boundaryfocus < 0.4f)
                    score *= 0.30f;

                world.setmetalorereserve(x, y, clampscore(score));
            }
        }
    }, 16);
}

void generate_placer_reserves(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int maximumflow = std::max(1, world.maxriverflow());
    const size_t columns = size_t(width) + 1;
    const size_t count = columns * (size_t(height) + 1);
    std::vector<int32_t> receiver(count, -1);
    std::vector<uint8_t> source(count, 0);
    for (int y = 0; y <= height; ++y)
    {
        for (int x = 0; x <= width; ++x)
        {
            const size_t index = size_t(y) * columns + size_t(x);
            if (world.sea(x, y) != 0 && !isdeltacell(world, x, y)) continue;
            source[index] = static_cast<uint8_t>(clampscore(
                float(world.metalorereserve(x, y)) * float(world.erosionpotential(x, y)) / 100.0f));
            int direction = world.riverdir(x, y);
            if (world.deltadir(x, y) > 0)
                direction = world.deltadir(x, y);
            if (direction < 1 || direction > 8)
            {
                int lowest = world.nom(x, y);
                for (int candidate = 1; candidate <= 8; ++candidate)
                {
                    const auto next = getdestination(x, y, candidate);
                    if (next.y < 0 || next.y > height) continue;
                    const int elevation = world.nom(wrap(next.x, width), next.y);
                    if (elevation < lowest) { lowest = elevation; direction = candidate; }
                }
            }
            if (direction < 1 || direction > 8) continue;
            const auto next = getdestination(x, y, direction);
            if (next.y < 0 || next.y > height) continue;
            const size_t destination = size_t(next.y) * columns + size_t(wrap(next.x, width));
            if (destination != index) receiver[index] = static_cast<int32_t>(destination);
        }
    }
    const auto routedore = route_ore_potential(receiver, source);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const float depositionnorm = static_cast<float>(world.depositionpotential(x, y)) / 100.0f;
                const float flownorm = std::log1p(static_cast<float>(averageriverflow(world, x, y))) / std::log1p(static_cast<float>(maximumflow));
                const int upstreamore = routedore[size_t(y) * columns + size_t(x)];
                const float upstreamorenorm = static_cast<float>(upstreamore) / 100.0f;

                float score = 100.0f * upstreamorenorm * (0.35f + 0.32f * depositionnorm
                    + 0.20f * flownorm + 0.13f * (isdeltacell(world, x, y) ? 1.0f : 0.0f));

                if (world.sea(x, y) != 0 && isdeltacell(world, x, y) == false)
                    score *= 0.3f;

                world.setplacerreserve(x, y, clampscore(score));
            }
        }
    }, 16);
}

void generate_evaporite_reserves(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                if (world.sea(x, y) != 0)
                {
                    world.setevaporitereserve(x, y, 0);
                    continue;
                }

                const float aridity = 1.0f - clampunit(static_cast<float>(annualrainfall(world, x, y)) / 350.0f);
                const float endorheic = world.basinclass(x, y) == BasinClass::endorheic ? 1.0f : 0.0f;
                const float saltbias = issaltfeature(world, x, y) ? 1.0f : 0.0f;
                const float terminalwater = (islakecell(world, x, y) || issaltfeature(world, x, y)) ? 1.0f : 0.0f;

                float score = 100.0f * (0.45f * endorheic + 0.25f * aridity + 0.20f * saltbias + 0.10f * terminalwater);

                world.setevaporitereserve(x, y, clampscore(score));
            }
        }
    }, 16);
}
}

void generate_mineral_resources(planet& world)
{
    // Ore scores read volcanic reserves; placers read ore and deposition.
    generate_volcanic_reserves(world);
    generate_metal_ore_reserves(world);
    generate_placer_reserves(world);
    generate_evaporite_reserves(world);
}
}
