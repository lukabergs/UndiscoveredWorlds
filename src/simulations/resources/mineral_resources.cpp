#include "simulations/resources/mineral_resources.hpp"
#include "simulations/resources/physical_layers_internal.hpp"
#include "simulations/resources/placer_transport.hpp"

namespace physical_layers
{
using namespace detail;
using namespace cellqueries;

namespace
{
template <typename Getter>
int maxscoreinradius(const planet& world, CellId c, int radius, Getter&& getter)
{
    int best = 0;
    anyinbox(world, c, radius, [&](CellId n)
    {
        best = std::max(best, static_cast<int>(getter(n)));
        return false;
    });
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

float boundaryproximityfactor(const planet& world, CellId c, float maxdistance)
{
    if (maxdistance <= 0.0f)
        return 0.0f;

    return 1.0f - clampunit(static_cast<float>(world.tectonicboundarydistance(c)) / maxdistance);
}

void generate_volcanic_reserves(planet& world)
{
    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            const GeologicRegime regime = world.geologicregime(c);
            const BoundaryType boundarytype = world.tectonicboundarytype(c);
            const DeformingRegionType regiontype = world.tectonicdeformingregiontype(c);
            const CrustClass crustclass = world.tectoniccrustclass(c);
            const float regimemult = regimevolcanicbias(regime);
            const float divergencenorm = static_cast<float>(world.tectonicdivergence(c)) / 100.0f;
            const float convergencenorm = static_cast<float>(world.tectonicconvergence(c)) / 100.0f;
            const float uplift = clampunit(world.tectonicuplifttendency(c));
            const float subsidence = clampunit(world.tectonicsubsidencetendency(c));
            const float strain = clampunit(world.tectonicaccumulatedstrain(c));
            const float deformation = clampunit(world.tectonicdeformationrate(c));
            const float boundaryhistory = clampunit(world.tectonicboundaryhistory(c));
            const float boundaryproximity = boundaryproximityfactor(world, c, 10.0f);
            const float boundaryfocus = boundaryvolcanicbias(boundarytype) * boundaryproximity * (0.45f + 0.55f * boundaryhistory);
            const float regionbias = deformingvolcanicbias(regiontype);
            const float crustbias = crustvolcanicbias(crustclass);
            const float ridgebias = (world.oceanridges(c) != 0 || world.oceanrifts(c) != 0) ? 1.0f : 0.0f;
            const float volcanobias = (world.volcanofloat(c) != 0 || world.strato(c)) ? 1.0f : 0.0f;
            const float nearbyvolcano = maxscoreinradius(world, c, 2, [&](CellId n)
            {
                return world.volcanofloat(n) != 0 || world.strato(n) ? 100 : 0;
            }) / 100.0f;

            float score = 100.0f * (0.20f * regimemult + 0.18f * boundaryfocus + 0.12f * divergencenorm
                + 0.08f * convergencenorm + 0.10f * uplift + 0.08f * strain + 0.06f * deformation
                + 0.06f * regionbias + 0.06f * crustbias + 0.16f * volcanobias + 0.10f * nearbyvolcano);
            score += 12.0f * ridgebias;

            if (world.sea(c) != 0 && ridgebias == 0.0f && boundaryfocus < 0.35f && regimemult < 0.6f)
                score *= 0.45f;

            if (subsidence > uplift && boundarytype == BoundaryType::passive_margin)
                score *= 0.7f;

            world.setvolcanicreserve(c, clampscore(score));
        }
    });
}

void generate_metal_ore_reserves(planet& world)
{
    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            const GeologicRegime regime = world.geologicregime(c);
            const BoundaryType boundarytype = world.tectonicboundarytype(c);
            const DeformingRegionType regiontype = world.tectonicdeformingregiontype(c);
            const CrustClass crustclass = world.tectoniccrustclass(c);
            const float regimebias = regimeorebias(regime);
            const float convergencenorm = static_cast<float>(world.tectonicconvergence(c)) / 100.0f;
            const float shearnorm = static_cast<float>(world.tectonicshear(c)) / 100.0f;
            const float uplift = clampunit(world.tectonicuplifttendency(c));
            const float strain = clampunit(world.tectonicaccumulatedstrain(c));
            const float deformation = clampunit(world.tectonicdeformationrate(c));
            const float boundaryhistory = clampunit(world.tectonicboundaryhistory(c));
            const float boundaryproximity = boundaryproximityfactor(world, c, 12.0f);
            const float boundaryfocus = boundaryorebias(boundarytype) * boundaryproximity * (0.40f + 0.60f * boundaryhistory);
            const float regionbias = deformingorebias(regiontype);
            const float crustbias = crustorebias(crustclass);
            const float volcanismnear = static_cast<float>(maxscoreinradius(world, c, 2, [&](CellId n)
            {
                return world.volcanicreserve(n);
            })) / 100.0f;
            const float mountainbias = clampunit(static_cast<float>(world.mountainheightfloat(c)) / 5000.0f);

            float score = 100.0f * (0.18f * regimebias + 0.18f * boundaryfocus + 0.15f * convergencenorm
                + 0.12f * shearnorm + 0.12f * strain + 0.08f * uplift + 0.05f * deformation
                + 0.05f * regionbias + 0.05f * crustbias + 0.10f * mountainbias + 0.10f * volcanismnear);

            if (boundarytype == BoundaryType::transform)
                score += 8.0f * boundaryhistory;

            if (regiontype == DeformingRegionType::diffuse_collision)
                score += 12.0f * (0.5f + 0.5f * boundaryhistory);

            if (world.sea(c) != 0 && !iscoast(world, c) && boundaryfocus < 0.4f)
                score *= 0.30f;

            world.setmetalorereserve(c, clampscore(score));
        }
    });
}

void generate_placer_reserves(planet& world)
{
    const int maximumflow = std::max(1, world.maxriverflow());
    const std::size_t count = world.cellcount();
    std::vector<uint8_t> source(count, 0);
    for (CellId c = 0; c < count; ++c)
    {
        if (world.sea(c) != 0 && !isdeltacell(world, c)) continue;
        source[c] = static_cast<uint8_t>(clampscore(
            float(world.metalorereserve(c)) * float(world.erosionpotential(c)) / 100.0f));
    }
    const auto routedore = route_ore_potential(drainagereceivers(world), source);

    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            const float depositionnorm = static_cast<float>(world.depositionpotential(c)) / 100.0f;
            const float flownorm = std::log1p(static_cast<float>(averageriverflow(world, c))) / std::log1p(static_cast<float>(maximumflow));
            const int upstreamore = routedore[c];
            const float upstreamorenorm = static_cast<float>(upstreamore) / 100.0f;

            float score = 100.0f * upstreamorenorm * (0.35f + 0.32f * depositionnorm
                + 0.20f * flownorm + 0.13f * (isdeltacell(world, c) ? 1.0f : 0.0f));

            if (world.sea(c) != 0 && isdeltacell(world, c) == false)
                score *= 0.3f;

            world.setplacerreserve(c, clampscore(score));
        }
    });
}

void generate_evaporite_reserves(planet& world)
{
    parallelforcells(world, [&](CellId first, CellId last)
    {
        for (CellId c = first; c <= last; c++)
        {
            if (world.sea(c) != 0)
            {
                world.setevaporitereserve(c, 0);
                continue;
            }

            const float aridity = 1.0f - clampunit(static_cast<float>(annualrainfall(world, c)) / 350.0f);
            const float endorheic = world.basinclass(c) == BasinClass::endorheic ? 1.0f : 0.0f;
            const float saltbias = issaltfeature(world, c) ? 1.0f : 0.0f;
            const float terminalwater = (islakecell(world, c) || issaltfeature(world, c)) ? 1.0f : 0.0f;

            float score = 100.0f * (0.45f * endorheic + 0.25f * aridity + 0.20f * saltbias + 0.10f * terminalwater);

            world.setevaporitereserve(c, clampscore(score));
        }
    });
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
