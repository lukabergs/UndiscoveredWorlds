#ifndef TECTONIC_CONTRACT_HPP
#define TECTONIC_CONTRACT_HPP

#include "tectonic_scenario.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

class lithosphere;

namespace platec::contract {

inline constexpr uint32_t kPhase0SchemaVersion = 1;
inline constexpr uint32_t kPhase2SchemaVersion = 2;
inline constexpr uint32_t kPhase3SchemaVersion = 3;
inline constexpr uint32_t kPhase4SchemaVersion = 4;
inline constexpr uint32_t kPhase5SchemaVersion = 5;
inline constexpr float kAdapterProxyConvergentBoundaryThreshold = 0.18f;
inline constexpr uint32_t kNoBoundaryId = 0;
inline constexpr uint32_t kNoDeformingRegionId = 0;

enum class GeologicRegime : uint8_t {
    Stable = 0,
    ConvergentArc = 1,
    ContinentCollision = 2,
    DivergentRift = 3,
    Transform = 4,
    PassiveMargin = 5,
    MidOceanRidge = 6,
    TrenchAdjacent = 7,
};

enum class CrustClass : uint8_t {
    None = 0,
    Oceanic = 1,
    Transitional = 2,
    Continental = 3,
};

enum class BoundaryType : uint8_t {
    None = 0,
    Convergent = 1,
    Divergent = 2,
    Transform = 3,
    PassiveMargin = 4,
};

enum class DeformingRegionType : uint8_t {
    None = 0,
    ContinentalRift = 1,
    DiffuseCollision = 2,
};

struct PlateKinematics {
    float unit_velocity_x = 0.0f;
    float unit_velocity_y = 0.0f;
    float linear_velocity_x = 0.0f;
    float linear_velocity_y = 0.0f;
    float angular_velocity = 0.0f;
    float mass_center_x = 0.0f;
    float mass_center_y = 0.0f;
};

struct BoundarySegment {
    uint32_t id = kNoBoundaryId;
    uint32_t left_plate_id = 0;
    uint32_t right_plate_id = 0;
    uint32_t cell_count = 0;
    uint32_t persistence_steps = 0;
    float centroid_x = 0.0f;
    float centroid_y = 0.0f;
    float length_cells = 0.0f;
    float average_normal_motion = 0.0f;
    float average_shear_motion = 0.0f;
    uint8_t average_convergence_score = 0;
    uint8_t average_divergence_score = 0;
    uint8_t average_shear_score = 0;
    BoundaryType boundary_type = BoundaryType::None;
    GeologicRegime geologic_regime = GeologicRegime::Stable;
    double age_myr = 0.0;
};

struct DeformingRegion {
    uint32_t id = kNoDeformingRegionId;
    uint32_t boundary_segment_id = kNoBoundaryId;
    uint32_t primary_plate_id = 0;
    uint32_t secondary_plate_id = 0;
    uint32_t cell_count = 0;
    uint32_t persistence_steps = 0;
    float centroid_x = 0.0f;
    float centroid_y = 0.0f;
    float average_deformation_rate = 0.0f;
    float average_interpolated_velocity_x = 0.0f;
    float average_interpolated_velocity_y = 0.0f;
    float average_normal_motion = 0.0f;
    float average_shear_motion = 0.0f;
    DeformingRegionType type = DeformingRegionType::None;
    double age_myr = 0.0;
};

struct Snapshot {
    uint32_t schema_version = kPhase5SchemaVersion;
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t iteration_count = 0;
    uint32_t cycle_count = 0;
    uint32_t time_origin_step = 0;
    uint16_t sea_level_m = 0;
    double time_myr = 0.0;
    double delta_time_myr = platec::scenario::kDefaultDeltaTimeMyr;

    platec::scenario::Scenario run_scenario;

    std::vector<float> heightmap;
    std::vector<uint32_t> plate_id;
    std::vector<uint32_t> crust_age_steps;
    std::vector<float> crust_age_myr;
    std::vector<float> crust_thickness;
    std::vector<uint8_t> crust_class;
    std::vector<float> uplift_tendency;
    std::vector<float> subsidence_tendency;
    std::vector<float> accumulated_strain;
    std::vector<uint8_t> boundary_type;
    std::vector<uint16_t> boundary_distance;
    std::vector<uint32_t> boundary_segment_id;
    std::vector<uint32_t> nearest_boundary_id;
    std::vector<uint32_t> deforming_region_id;
    std::vector<uint8_t> deforming_region_type;
    std::vector<float> deformation_rate;
    std::vector<float> deformation_velocity_x;
    std::vector<float> deformation_velocity_y;

    std::vector<uint8_t> convergence_score;
    std::vector<uint8_t> divergence_score;
    std::vector<uint8_t> shear_score;
    std::vector<uint8_t> geologic_regime;

    std::vector<PlateKinematics> plates;
    std::vector<BoundarySegment> boundary_segments;
    std::vector<DeformingRegion> deforming_regions;

    size_t cell_count() const {
        return static_cast<size_t>(width) * static_cast<size_t>(height);
    }
};

Snapshot capture_snapshot(const lithosphere& lithosphere);
Snapshot capture_phase0_snapshot(const lithosphere& lithosphere);
const char* geologic_regime_name(GeologicRegime regime);
const char* crust_class_name(CrustClass crust_class);
const char* boundary_type_name(BoundaryType boundary_type);
const char* deforming_region_type_name(DeformingRegionType region_type);

} // namespace platec::contract

#endif
