#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "platecapi.hpp"
#include "tectonic_contract.hpp"
#include "topography_codec.hpp"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#include "functions.hpp"
#include "planet.hpp"

namespace
{
constexpr int kMinAggregationOverlapAbsolute = 64;
constexpr int kAggregationOverlapAreaDivisor = 1000;
constexpr auto kVisualizationInterval = std::chrono::milliseconds(50);
constexpr double kNativeHeightScale = 10000.0;

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

float maxnativeheight(const float* heightmap, std::size_t cellcount)
{
    if (heightmap == nullptr || cellcount == 0)
        return 0.0f;

    float maximum = heightmap[0];

    for (std::size_t index = 1; index < cellcount; index++)
        maximum = std::max(maximum, heightmap[index]);

    return maximum;
}

int encodenativeheight(float value)
{
    return static_cast<int>(std::llround(static_cast<double>(value) * kNativeHeightScale));
}

void applynativeheightmap(planet& world, std::vector<std::vector<bool>>& shelves, const float* outputheightmap, std::size_t cellcount)
{
    if (outputheightmap == nullptr)
        return;

    const int width = world.width();
    const int height = world.height();
    const int simwidth = width + 1;
    const float nativeheightmaximum = maxnativeheight(outputheightmap, cellcount);
    const int sealevel = encodenativeheight(TopographyCodec::kContinentalBase);
    const int maxelevation = std::max(encodenativeheight(TopographyCodec::kContinentalBase + 1.0f) + 1,
        encodenativeheight(nativeheightmaximum) + 1);

    world.setsealevel(sealevel);
    world.setmaxelevation(maxelevation);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const std::size_t index = static_cast<std::size_t>(y) * static_cast<std::size_t>(simwidth) + static_cast<std::size_t>(x);
                world.setnom(x, y, encodenativeheight(outputheightmap[index]));
                shelves[x][y] = false;
            }
        }
    });

    getlandandseatotals(world);
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

int defaultplatetectonicsaggregationoverlapabs(int width, int height)
{
    const std::uint64_t area = static_cast<std::uint64_t>(std::max(1, width)) * static_cast<std::uint64_t>(std::max(1, height));
    return std::max(kMinAggregationOverlapAbsolute, static_cast<int>(area / kAggregationOverlapAreaDivisor));
}

void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves, int cyclecount, int platecount)
{
    PlateTectonicsSimulationOptions options;
    options.cycleCount = std::max(1, cyclecount);
    options.plateCount = std::max(1, platecount);
    applyplatetectonicssimulation(world, shelves, options);
}

void applyplatetectonicssimulation(planet& world, std::vector<std::vector<bool>>& shelves, const PlateTectonicsSimulationOptions& options)
{
    const int width = world.width();
    const int height = world.height();
    const int simwidth = width + 1;
    const int simheight = height + 1;
    const std::size_t cellcount = static_cast<std::size_t>(simwidth) * static_cast<std::size_t>(simheight);
    const uint32_t normalizedcyclecount = static_cast<uint32_t>(std::max(1, options.cycleCount));
    const uint32_t normalizedcyclesteplimit = static_cast<uint32_t>(std::max(0, options.cycleStepLimit));
    const uint32_t normalizedplatecount = static_cast<uint32_t>(std::max(1, options.plateCount));
    const uint32_t normalizederosionperiod = static_cast<uint32_t>(std::max(1, options.erosionPeriod));
    const uint32_t aggregationoverlapabsolute = static_cast<uint32_t>(options.aggregationOverlapAbsolute >= 0
        ? options.aggregationOverlapAbsolute
        : defaultplatetectonicsaggregationoverlapabs(simwidth, simheight));
    const float aggregationoverlaprelative = clamp01(options.aggregationOverlapRelative);
    const float foldingratio = clamp01(options.foldingRatio);
    const float erosionstrength = std::max(0.0f, options.erosionStrength);
    const float landmassrotation = std::max(0.0f, options.landmassRotation);
    const float rotationstrength = std::max(0.0f, options.rotationStrength);
    const float subductionstrength = clamp01(options.subductionStrength);
    const float divergentcarvestrength = std::max(0.0f, options.divergentCarveStrength);
    const double deltatimemyr = std::max(0.001, static_cast<double>(options.deltaTimeMyr));
    const int32_t sealevelmeters = options.useSeaLevelMeters ? std::clamp(options.seaLevelMeters, 0, 65535) : TopographyCodec::kNoSeaLevelOverride;

    world.cleartectonicprovenance();

    platec::scenario::Scenario scenario;
    scenario.seed = world.seed();
    scenario.width = static_cast<uint32_t>(simwidth);
    scenario.height = static_cast<uint32_t>(simheight);
    scenario.erosion_period = normalizederosionperiod;
    scenario.folding_ratio = foldingratio;
    scenario.aggregation_overlap_abs = aggregationoverlapabsolute;
    scenario.aggregation_overlap_rel = aggregationoverlaprelative;
    scenario.cycle_count = normalizedcyclecount;
    scenario.cycle_step_limit = normalizedcyclesteplimit;
    scenario.plate_count = normalizedplatecount;
    scenario.erosion_strength = erosionstrength;
    scenario.crust_rotation_strength = landmassrotation;
    scenario.rotation_strength = rotationstrength;
    scenario.subduction_strength = subductionstrength;
    scenario.divergent_carve_strength = divergentcarvestrength;
    scenario.sea_level_m = sealevelmeters;
    scenario.delta_time_myr = deltatimemyr;

    PlateTectonicsHandle simulation;
    simulation.pointer = platec_api_create_from_scenario(&scenario);

    if (simulation.pointer == nullptr)
        return;

    const bool visualizeprogress = hasworldgenvisualizationcallback();
    auto nextvisualization = std::chrono::steady_clock::now();

    while (platec_api_is_finished(simulation.pointer) == 0)
    {
        platec_api_step(simulation.pointer);

        if (visualizeprogress == false || std::chrono::steady_clock::now() < nextvisualization)
            continue;

        float* previewheightmap = platec_api_get_heightmap(simulation.pointer);

        if (previewheightmap != nullptr)
        {
            applynativeheightmap(world, shelves, previewheightmap, cellcount);
            requestworldgenvisualization();
        }

        nextvisualization = std::chrono::steady_clock::now() + kVisualizationInterval;
    }

    float* outputheightmap = platec_api_get_heightmap(simulation.pointer);

    if (outputheightmap == nullptr)
        return;

    applynativeheightmap(world, shelves, outputheightmap, cellcount);

    if (visualizeprogress)
        requestworldgenvisualization();

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
        }
    }
}
