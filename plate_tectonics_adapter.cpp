#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "platecapi.hpp"
#include "tectonic_contract.hpp"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"

namespace
{
template <int N>
struct prioritytag : prioritytag<N - 1>
{
};

template <>
struct prioritytag<0>
{
};

struct PlateTectonicsHandle
{
    void* pointer = nullptr;

    ~PlateTectonicsHandle()
    {
        if (pointer != nullptr)
            platec_api_destroy(pointer);
    }
};

float clamp01(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

float currentsearatio(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    int seacells = 0;
    const int totalcells = (width + 1) * (height + 1);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
            seacells += world.sea(x, y) ? 1 : 0;
    }

    return static_cast<float>(seacells) / static_cast<float>(std::max(1, totalcells));
}

void normalizeheightmap(std::vector<float>& heightmap)
{
    if (heightmap.empty())
        return;

    const auto bounds = std::minmax_element(heightmap.begin(), heightmap.end());
    const float lowest = *bounds.first;
    const float highest = *bounds.second;

    if (highest <= lowest)
    {
        std::fill(heightmap.begin(), heightmap.end(), 0.0f);
        return;
    }

    const float inverserange = 1.0f / (highest - lowest);

    for (float& value : heightmap)
        value = (value - lowest) * inverserange;
}

float findseathreshold(const std::vector<float>& heightmap, float searatio)
{
    float threshold = 0.5f;
    float step = 0.5f;
    const std::size_t cellcount = heightmap.size();

    while (step > 0.0005f)
    {
        std::size_t seacells = 0;

        for (float value : heightmap)
            seacells += value < threshold ? 1U : 0U;

        step *= 0.5f;

        if (static_cast<float>(seacells) / static_cast<float>(cellcount) < searatio)
            threshold += step;
        else
            threshold -= step;
    }

    return clamp01(threshold);
}

GeologicRegime translategeologicregime(std::uint8_t regime)
{
    using NativeRegime = platec::contract::GeologicRegime;

    switch (static_cast<NativeRegime>(regime))
    {
    case NativeRegime::Stable:
        return GeologicRegime::stable;
    case NativeRegime::ConvergentArc:
        return GeologicRegime::convergent_arc;
    case NativeRegime::ContinentCollision:
        return GeologicRegime::continent_collision;
    case NativeRegime::DivergentRift:
        return GeologicRegime::divergent_rift;
    case NativeRegime::Transform:
        return GeologicRegime::transform;
    case NativeRegime::PassiveMargin:
        return GeologicRegime::passive_margin;
    case NativeRegime::MidOceanRidge:
        return GeologicRegime::mid_ocean_ridge;
    case NativeRegime::TrenchAdjacent:
        return GeologicRegime::trench_adjacent;
    default:
        return GeologicRegime::stable;
    }
}

CrustClass translatecrustclass(std::uint8_t crustclass)
{
    using NativeCrustClass = platec::contract::CrustClass;

    switch (static_cast<NativeCrustClass>(crustclass))
    {
    case NativeCrustClass::Oceanic:
        return CrustClass::oceanic;
    case NativeCrustClass::Transitional:
        return CrustClass::transitional;
    case NativeCrustClass::Continental:
        return CrustClass::continental;
    case NativeCrustClass::None:
    default:
        return CrustClass::none;
    }
}

BoundaryType translateboundarytype(std::uint8_t boundarytype)
{
    using NativeBoundaryType = platec::contract::BoundaryType;

    switch (static_cast<NativeBoundaryType>(boundarytype))
    {
    case NativeBoundaryType::Convergent:
        return BoundaryType::convergent;
    case NativeBoundaryType::Divergent:
        return BoundaryType::divergent;
    case NativeBoundaryType::Transform:
        return BoundaryType::transform;
    case NativeBoundaryType::PassiveMargin:
        return BoundaryType::passive_margin;
    case NativeBoundaryType::None:
    default:
        return BoundaryType::none;
    }
}

DeformingRegionType translatedeformingregiontype(std::uint8_t regiontype)
{
    using NativeRegionType = platec::contract::DeformingRegionType;

    switch (static_cast<NativeRegionType>(regiontype))
    {
    case NativeRegionType::ContinentalRift:
        return DeformingRegionType::continental_rift;
    case NativeRegionType::DiffuseCollision:
        return DeformingRegionType::diffuse_collision;
    case NativeRegionType::None:
    default:
        return DeformingRegionType::none;
    }
}

float normalizeboundaryhistory(const platec::contract::BoundarySegment* segment)
{
    if (segment == nullptr)
        return 0.0f;

    const float agesignal = clamp01(static_cast<float>(segment->age_myr) / 40.0f);
    const float persistencesignal = clamp01(static_cast<float>(segment->persistence_steps) / 24.0f);
    return clamp01(0.55f * agesignal + 0.45f * persistencesignal);
}

float boundaryupliftbias(BoundaryType boundarytype)
{
    switch (boundarytype)
    {
    case BoundaryType::convergent:
        return 1.0f;
    case BoundaryType::transform:
        return 0.45f;
    case BoundaryType::divergent:
        return 0.28f;
    case BoundaryType::passive_margin:
        return 0.10f;
    case BoundaryType::none:
    default:
        return 0.0f;
    }
}

float deformingupliftbias(DeformingRegionType regiontype)
{
    switch (regiontype)
    {
    case DeformingRegionType::diffuse_collision:
        return 1.0f;
    case DeformingRegionType::continental_rift:
        return 0.45f;
    case DeformingRegionType::none:
    default:
        return 0.0f;
    }
}

template <typename World>
auto settectoniccyclecount(World& world, int value, prioritytag<1>) -> decltype(world.settectoniccyclecount(value), void())
{
    world.settectoniccyclecount(value);
}

template <typename World>
void settectoniccyclecount(World&, int, prioritytag<0>)
{
}

template <typename World>
auto settectonicplatecount(World& world, int value, prioritytag<1>) -> decltype(world.settectonicplatecount(value), void())
{
    world.settectonicplatecount(value);
}

template <typename World>
void settectonicplatecount(World&, int, prioritytag<0>)
{
}

template <typename World>
auto settectonicsealevelm(World& world, int value, prioritytag<1>) -> decltype(world.settectonicsealevelm(value), void())
{
    world.settectonicsealevelm(value);
}

template <typename World>
void settectonicsealevelm(World&, int, prioritytag<0>)
{
}

template <typename World, typename Segments>
auto settectonicboundarysegments(World& world, Segments&& segments, prioritytag<1>) -> decltype(world.settectonicboundarysegments(std::forward<Segments>(segments)), void())
{
    world.settectonicboundarysegments(std::forward<Segments>(segments));
}

template <typename World, typename Segments>
void settectonicboundarysegments(World&, Segments&&, prioritytag<0>)
{
}

template <typename World, typename Regions>
auto settectonicdeformingregions(World& world, Regions&& regions, prioritytag<1>) -> decltype(world.settectonicdeformingregions(std::forward<Regions>(regions)), void())
{
    world.settectonicdeformingregions(std::forward<Regions>(regions));
}

template <typename World, typename Regions>
void settectonicdeformingregions(World&, Regions&&, prioritytag<0>)
{
}
}

void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves)
{
    const int width = world.width();
    const int height = world.height();
    const int simwidth = width + 1;
    const int simheight = height + 1;
    const std::size_t cellcount = static_cast<std::size_t>(simwidth) * static_cast<std::size_t>(simheight);
    const int sealevel = world.sealevel();
    const int maxelev = world.maxelevation();
    const float searatio = std::clamp(currentsearatio(world) + tuning::terrain::platetectonics::seaLevelBias, 0.05f, 0.95f);

    world.cleartectonicprovenance();

    std::vector<float> inputheightmap(cellcount, 0.0f);
    std::vector<bool> originalsea(cellcount, false);
    int lowest = world.nom(0, 0);
    int highest = world.nom(0, 0);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const int value = world.nom(x, y);
            lowest = std::min(lowest, value);
            highest = std::max(highest, value);
        }
    }

    const float inverserange = highest > lowest ? 1.0f / static_cast<float>(highest - lowest) : 0.0f;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
            const int current = world.nom(x, y);
            inputheightmap[index] = inverserange > 0.0f ? static_cast<float>(current - lowest) * inverserange : 0.0f;
            originalsea[index] = world.sea(x, y) != 0;
        }
    }

    PlateTectonicsHandle simulation;
    simulation.pointer = platec_api_create(world.seed(), static_cast<uint32_t>(simwidth), static_cast<uint32_t>(simheight),
        searatio,
        tuning::terrain::platetectonics::erosionPeriod,
        tuning::terrain::platetectonics::foldingRatio,
        tuning::terrain::platetectonics::aggregationOverlapAbsolute,
        tuning::terrain::platetectonics::aggregationOverlapRelative,
        static_cast<uint32_t>(platetectonicscyclecount()),
        tuning::terrain::platetectonics::plateCount);

    if (simulation.pointer == nullptr)
        return;

    platec_api_load_heightmap(simulation.pointer, inputheightmap.data(), searatio);

    for (uint32_t step = 0; step < tuning::terrain::platetectonics::maximumSimulationSteps; step++)
    {
        if (platec_api_is_finished(simulation.pointer) != 0)
            break;

        platec_api_step(simulation.pointer);
    }

    float* outputheightmap = platec_api_get_heightmap(simulation.pointer);

    if (outputheightmap == nullptr)
        return;

    std::vector<float> transformed(outputheightmap, outputheightmap + cellcount);
    normalizeheightmap(transformed);

    if (tuning::terrain::platetectonics::outputBlend < 1.0f)
    {
        const float blend = clamp01(tuning::terrain::platetectonics::outputBlend);

        for (std::size_t index = 0; index < cellcount; index++)
            transformed[index] = inputheightmap[index] * (1.0f - blend) + transformed[index] * blend;
    }

    const float landbiasedsearatio = clamp01(searatio - tuning::terrain::platetectonics::landRetentionSeaBias);
    const float seathreshold = findseathreshold(transformed, landbiasedsearatio);
    const float searange = std::max(0.0001f, seathreshold);
    const float landrange = std::max(0.0001f, 1.0f - seathreshold);
    const int oceanceiling = std::max(tuning::terrain::platetectonics::minimumOceanDepth, sealevel - tuning::terrain::platetectonics::coastalOceanOffset);
    const int landfloor = std::min(maxelev - 1, sealevel + tuning::terrain::platetectonics::landStartOffset);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
                const float value = transformed[index];

                if (originalsea[index])
                {
                    const float oceanrelative = std::pow(clamp01(value / searange), tuning::terrain::platetectonics::oceanExponent);
                    const int newvalue = tuning::terrain::platetectonics::minimumOceanDepth +
                        static_cast<int>(std::round(static_cast<float>(oceanceiling - tuning::terrain::platetectonics::minimumOceanDepth) * oceanrelative));
                    world.setnom(x, y, std::max(tuning::terrain::platetectonics::minimumOceanDepth, std::min(oceanceiling, newvalue)));
                }
                else if (value <= seathreshold)
                {
                    const float oceanrelative = std::pow(clamp01(value / searange), tuning::terrain::platetectonics::oceanExponent);
                    const int newvalue = tuning::terrain::platetectonics::minimumOceanDepth +
                        static_cast<int>(std::round(static_cast<float>(oceanceiling - tuning::terrain::platetectonics::minimumOceanDepth) * oceanrelative));
                    world.setnom(x, y, std::max(tuning::terrain::platetectonics::minimumOceanDepth, std::min(oceanceiling, newvalue)));
                }
                else
                {
                    const float landrelative = std::pow(clamp01((value - seathreshold) / landrange), tuning::terrain::platetectonics::landExponent);
                    const int newvalue = landfloor +
                        static_cast<int>(std::round(static_cast<float>((maxelev - 1) - landfloor) * landrelative));
                    world.setnom(x, y, std::max(landfloor, std::min(maxelev - 1, newvalue)));
                }

                shelves[x][y] = false;
            }
        }
    });

    const std::uint8_t* convergencemap = platec_api_get_convergence_map(simulation.pointer);
    const std::uint8_t* divergencemap = platec_api_get_divergence_map(simulation.pointer);
    const std::uint8_t* shearmap = platec_api_get_shear_map(simulation.pointer);
    const std::uint8_t* regimemap = platec_api_get_geologic_regime_map(simulation.pointer);
    const float* crustagemyrmap = platec_api_get_crust_age_myr_map(simulation.pointer);
    const float* crustthicknessmap = platec_api_get_crust_thickness_map(simulation.pointer);
    const std::uint8_t* crustclassmap = platec_api_get_crust_class_map(simulation.pointer);
    const float* uplifttendencymap = platec_api_get_uplift_tendency_map(simulation.pointer);
    const float* subsidencetendencymap = platec_api_get_subsidence_tendency_map(simulation.pointer);
    const float* accumulatedstrainmap = platec_api_get_accumulated_strain_map(simulation.pointer);
    const std::uint8_t* boundarytypemap = platec_api_get_boundary_type_map(simulation.pointer);
    const std::uint16_t* boundarydistancemap = platec_api_get_boundary_distance_map(simulation.pointer);
    const std::uint32_t* boundarysegmentidmap = platec_api_get_boundary_segment_id_map(simulation.pointer);
    const std::uint32_t* nearestboundaryidmap = platec_api_get_nearest_boundary_id_map(simulation.pointer);
    const std::uint32_t* deformingregionidmap = platec_api_get_deforming_region_id_map(simulation.pointer);
    const std::uint8_t* deformingregiontypemap = platec_api_get_deforming_region_type_map(simulation.pointer);
    const float* deformationratemap = platec_api_get_deformation_rate_map(simulation.pointer);
    const float* deformationvelocityxmap = platec_api_get_deformation_velocity_x_map(simulation.pointer);
    const float* deformationvelocityymap = platec_api_get_deformation_velocity_y_map(simulation.pointer);

    if (convergencemap == nullptr || divergencemap == nullptr || shearmap == nullptr || regimemap == nullptr
        || crustagemyrmap == nullptr || crustthicknessmap == nullptr || crustclassmap == nullptr
        || uplifttendencymap == nullptr || subsidencetendencymap == nullptr || accumulatedstrainmap == nullptr
        || boundarytypemap == nullptr || boundarydistancemap == nullptr || boundarysegmentidmap == nullptr
        || nearestboundaryidmap == nullptr || deformingregionidmap == nullptr || deformingregiontypemap == nullptr
        || deformationratemap == nullptr || deformationvelocityxmap == nullptr || deformationvelocityymap == nullptr)
        return;

    world.settectonictimeoriginstep(static_cast<int>(platec_api_get_time_origin_step(simulation.pointer)));
    world.settectonictimemyr(static_cast<float>(platec_api_get_time_myr(simulation.pointer)));
    world.settectonicdeltatimemyr(static_cast<float>(platec_api_get_delta_time_myr(simulation.pointer)));
    settectoniccyclecount(world, static_cast<int>(platec_api_get_cycle_count(simulation.pointer)), prioritytag<1>{});
    settectonicplatecount(world, static_cast<int>(platec_api_get_plate_count(simulation.pointer)), prioritytag<1>{});
    settectonicsealevelm(world, static_cast<int>(platec_api_get_sea_level_m(simulation.pointer)), prioritytag<1>{});

    const std::uint32_t boundarysegmentcount = platec_api_get_boundary_segment_count(simulation.pointer);
    const platec::contract::BoundarySegment* boundarysegments = platec_api_get_boundary_segments(simulation.pointer);
    std::vector<platec::contract::BoundarySegment> boundarysegmentobjects;
    if (boundarysegments != nullptr && boundarysegmentcount > 0)
        boundarysegmentobjects.assign(boundarysegments, boundarysegments + boundarysegmentcount);
    settectonicboundarysegments(world, std::move(boundarysegmentobjects), prioritytag<1>{});

    const std::uint32_t deformingregioncount = platec_api_get_deforming_region_count(simulation.pointer);
    const platec::contract::DeformingRegion* deformingregions = platec_api_get_deforming_regions(simulation.pointer);
    std::vector<platec::contract::DeformingRegion> deformingregionobjects;
    if (deformingregions != nullptr && deformingregioncount > 0)
        deformingregionobjects.assign(deformingregions, deformingregions + deformingregioncount);
    settectonicdeformingregions(world, std::move(deformingregionobjects), prioritytag<1>{});

    std::vector<const platec::contract::BoundarySegment*> boundarysegmentlookup;
    if (boundarysegments != nullptr && boundarysegmentcount > 0)
    {
        std::size_t maxid = 0;
        for (std::uint32_t i = 0; i < boundarysegmentcount; i++)
            maxid = std::max(maxid, static_cast<std::size_t>(boundarysegments[i].id));

        boundarysegmentlookup.assign(maxid + 1, nullptr);
        for (std::uint32_t i = 0; i < boundarysegmentcount; i++)
            boundarysegmentlookup[boundarysegments[i].id] = &boundarysegments[i];
    }

    std::vector<std::vector<int>> rawmountains(ARRAYWIDTH, std::vector<int>(ARRAYHEIGHT, 0));
    bool havemountains = false;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
            const int convergencescore = static_cast<int>(convergencemap[index]);
            const int divergencescore = static_cast<int>(divergencemap[index]);
            const int shearscore = static_cast<int>(shearmap[index]);
            const CrustClass crustclass = translatecrustclass(crustclassmap[index]);
            const BoundaryType boundarytype = translateboundarytype(boundarytypemap[index]);
            const DeformingRegionType regiontype = translatedeformingregiontype(deformingregiontypemap[index]);
            const std::uint32_t nearestboundaryid = nearestboundaryidmap[index];
            const platec::contract::BoundarySegment* boundarysegment =
                nearestboundaryid < boundarysegmentlookup.size() ? boundarysegmentlookup[nearestboundaryid] : nullptr;
            // This history score is app-derived for presentation/stylization, not native contract data.
            const float boundaryhistory = normalizeboundaryhistory(boundarysegment);

            world.settectonicconvergence(x, y, convergencescore);
            world.settectonicdivergence(x, y, divergencescore);
            world.settectonicshear(x, y, shearscore);
            world.setgeologicregime(x, y, translategeologicregime(regimemap[index]));
            world.settectoniccrustagemyr(x, y, crustagemyrmap[index]);
            world.settectoniccrustthickness(x, y, crustthicknessmap[index]);
            world.settectoniccrustclass(x, y, crustclass);
            world.settectonicuplifttendency(x, y, uplifttendencymap[index]);
            world.settectonicsubsidencetendency(x, y, subsidencetendencymap[index]);
            world.settectonicaccumulatedstrain(x, y, accumulatedstrainmap[index]);
            world.settectonicboundarytype(x, y, boundarytype);
            world.settectonicboundarydistance(x, y, static_cast<int>(boundarydistancemap[index]));
            world.settectonicboundarysegmentid(x, y, static_cast<int>(boundarysegmentidmap[index]));
            world.settectonicnearestboundaryid(x, y, static_cast<int>(nearestboundaryid));
            world.settectonicboundaryhistory(x, y, boundaryhistory);
            world.settectonicdeformingregionid(x, y, static_cast<int>(deformingregionidmap[index]));
            world.settectonicdeformingregiontype(x, y, regiontype);
            world.settectonicdeformationrate(x, y, deformationratemap[index]);
            world.settectonicdeformationvelocityx(x, y, deformationvelocityxmap[index]);
            world.settectonicdeformationvelocityy(x, y, deformationvelocityymap[index]);

            if (originalsea[index] || world.nom(x, y) <= sealevel + 25)
                continue;

            const float convergencesignal = clamp01(static_cast<float>(convergencescore) / 100.0f);
            const float upliftsignal = clamp01(uplifttendencymap[index]);
            const float subsidencesignal = clamp01(subsidencetendencymap[index]);
            const float strainsignal = clamp01(accumulatedstrainmap[index]);
            const float deformationsignal = clamp01(deformationratemap[index]);
            const float boundaryproximity = 1.0f - clamp01(static_cast<float>(boundarydistancemap[index]) / 8.0f);
            const float crustbias = crustclass == CrustClass::continental ? 1.0f
                : crustclass == CrustClass::transitional ? 0.72f
                : crustclass == CrustClass::oceanic ? 0.28f
                : 0.0f;
            const float structuresignal = clamp01(boundaryproximity * boundaryupliftbias(boundarytype)
                + 0.55f * deformingupliftbias(regiontype));
            float mountainsource = clamp01(
                (0.44f * upliftsignal + 0.20f * strainsignal + 0.16f * convergencesignal
                    + 0.10f * boundaryhistory + 0.10f * deformationsignal) * (0.45f + 0.55f * crustbias)
                + 0.25f * structuresignal
                - 0.20f * subsidencesignal);

            if (regiontype == DeformingRegionType::continental_rift)
            {
                const float shouldersignal = clamp01(0.45f * upliftsignal + 0.30f * boundaryproximity
                    + 0.25f * deformationsignal - 0.10f * subsidencesignal);
                mountainsource = std::max(mountainsource, shouldersignal * 0.55f);
            }

            if (boundarytype == BoundaryType::transform)
                mountainsource *= 0.60f;
            else if (boundarytype == BoundaryType::passive_margin)
                mountainsource *= 0.25f;

            if (mountainsource <= 0.08f)
                continue;

            const bool riftshoulder = regiontype == DeformingRegionType::continental_rift;
            const float mountainsignal = std::pow(clamp01(mountainsource), riftshoulder ? 1.05f : 0.78f);
            const float upliftfactor = riftshoulder ? 0.55f : 1.0f;
            const int uplift = static_cast<int>(std::round(static_cast<float>(tuning::terrain::platetectonics::collisionUplift)
                * std::sqrt(clamp01(mountainsource)) * upliftfactor));
            const int minimumpeak = riftshoulder
                ? std::max(0, tuning::terrain::platetectonics::collisionMinimumPeak / 2)
                : tuning::terrain::platetectonics::collisionMinimumPeak;
            const int maximumpeak = riftshoulder
                ? std::max(minimumpeak + 1, (tuning::terrain::platetectonics::collisionMaximumPeak * 3) / 4)
                : tuning::terrain::platetectonics::collisionMaximumPeak;
            const int peakheight = minimumpeak +
                static_cast<int>(std::round(static_cast<float>(maximumpeak - minimumpeak) * mountainsignal));

            if (uplift > 0)
                world.setnom(x, y, std::min(maxelev - 1, world.nom(x, y) + uplift));

            rawmountains[x][y] = peakheight;
            havemountains = true;
        }
    }

    if (havemountains)
    {
        std::vector<std::vector<bool>> dummyok(ARRAYWIDTH, std::vector<bool>(ARRAYHEIGHT, false));
        createmountainsfromraw(world, rawmountains, dummyok);
    }
}
