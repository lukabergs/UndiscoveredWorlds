#include "tectonic_frame_reader.hpp"

#include "tectonic_api.hpp"
#include "tectonic_contract.hpp"
#include "topography_codec.hpp"

#include <algorithm>

namespace mineralsystems
{
namespace
{
using tectonics::contract::BoundaryType;
using tectonics::contract::CrustClass;
using tectonics::contract::GeologicRegime;

constexpr bool matches(std::uint8_t code, CrustClass value) { return code == static_cast<std::uint8_t>(value); }
constexpr bool matches(std::uint8_t code, BoundaryType value) { return code == static_cast<std::uint8_t>(value); }
constexpr bool matches(std::uint8_t code, GeologicRegime value) { return code == static_cast<std::uint8_t>(value); }

static_assert(matches(crustcode::oceanic, CrustClass::Oceanic) && matches(crustcode::transitional, CrustClass::Transitional)
    && matches(crustcode::continental, CrustClass::Continental) && matches(crustcode::none, CrustClass::None));
static_assert(matches(boundarycode::convergent, BoundaryType::Convergent) && matches(boundarycode::divergent, BoundaryType::Divergent)
    && matches(boundarycode::transform, BoundaryType::Transform) && matches(boundarycode::passive_margin, BoundaryType::PassiveMargin));
static_assert(matches(regimecode::convergent_arc, GeologicRegime::ConvergentArc)
    && matches(regimecode::continent_collision, GeologicRegime::ContinentCollision)
    && matches(regimecode::divergent_rift, GeologicRegime::DivergentRift) && matches(regimecode::transform, GeologicRegime::Transform)
    && matches(regimecode::passive_margin, GeologicRegime::PassiveMargin)
    && matches(regimecode::mid_ocean_ridge, GeologicRegime::MidOceanRidge)
    && matches(regimecode::trench_adjacent, GeologicRegime::TrenchAdjacent));
}

const TectonicFrame& TectonicFrameReader::read(void* simulation)
{
    frame_.width = tectonics_get_map_width(simulation);
    frame_.height = tectonics_get_map_height(simulation);
    frame_.timeMyr = tectonics_get_time_myr(simulation);
    frame_.deltaTimeMyr = tectonics_get_delta_time_myr(simulation);
    frame_.plateIndex = tectonics_get_platesmap(simulation);
    frame_.crustClass = tectonics_get_crust_class_map(simulation);
    frame_.crustAgeMyr = tectonics_get_crust_age_myr_map(simulation);
    frame_.boundaryDistance = tectonics_get_boundary_distance_map(simulation);
    frame_.nearestBoundaryId = tectonics_get_nearest_boundary_id_map(simulation);
    frame_.geologicRegime = tectonics_get_geologic_regime_map(simulation);

    const float* heightmap = tectonics_get_heightmap(simulation);
    const std::uint16_t sealevel = tectonics_get_sea_level_m(simulation);
    elevationM_.resize(frame_.cellCount());
    for (std::size_t cell = 0; cell < elevationM_.size(); ++cell)
        elevationM_[cell] = static_cast<float>(TopographyCodec::internal_to_meters(heightmap[cell], sealevel)) - sealevel;
    frame_.elevationM = elevationM_.data();

    const tectonics::columns::Cell* columns = tectonics_get_column_cells(simulation);
    continentalCrustM_.clear();
    if (columns != nullptr)
    {
        continentalCrustM_.resize(frame_.cellCount());
        for (std::size_t cell = 0; cell < continentalCrustM_.size(); ++cell)
            continentalCrustM_[cell] = static_cast<float>(columns[cell].continent);
    }
    frame_.continentalCrustM = continentalCrustM_.empty() ? nullptr : continentalCrustM_.data();

    const tectonics::transport::Displacement* displacement = tectonics_get_material_displacement(simulation);
    materialShift_.clear();
    if (displacement != nullptr)
    {
        materialShift_.resize(frame_.cellCount());
        for (std::size_t cell = 0; cell < materialShift_.size(); ++cell)
            materialShift_[cell] = { displacement[cell].x, displacement[cell].y };
    }
    frame_.materialShift = materialShift_.empty() ? nullptr : materialShift_.data();

    const std::uint32_t platecount = tectonics_get_plate_count(simulation);
    frame_.plates.resize(platecount);
    for (std::uint32_t plate = 0; plate < platecount; ++plate)
    {
        frame_.plates[plate] = { tectonics_get_plate_uid(simulation, plate) };
    }

    const std::uint32_t segmentcount = tectonics_get_boundary_segment_count(simulation);
    const tectonics::contract::BoundarySegment* segments = tectonics_get_boundary_segments(simulation);
    std::uint32_t maximumid = 0;
    for (std::uint32_t i = 0; segments != nullptr && i < segmentcount; ++i)
        maximumid = std::max(maximumid, segments[i].id);
    frame_.boundaries.assign(static_cast<std::size_t>(maximumid) + 1, BoundarySummary{});
    for (std::uint32_t i = 0; segments != nullptr && i < segmentcount; ++i)
    {
        const tectonics::contract::BoundarySegment& segment = segments[i];
        frame_.boundaries[segment.id] = { static_cast<std::uint8_t>(segment.boundary_type),
            static_cast<std::uint8_t>(segment.geologic_regime), segment.average_convergence_score,
            segment.average_shear_score, static_cast<float>(segment.age_myr), segment.left_plate_id,
            segment.right_plate_id };
    }
    return frame_;
}
}
