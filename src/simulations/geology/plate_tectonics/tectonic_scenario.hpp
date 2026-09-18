#ifndef TECTONIC_SCENARIO_HPP
#define TECTONIC_SCENARIO_HPP

#include "topography_codec.hpp"

#include <cstdint>

namespace platec::scenario {

inline constexpr uint32_t kScenarioVersion = 1;
inline constexpr double kDefaultDeltaTimeMyr = 1.0;
inline constexpr float kDefaultGravityMps2 = 9.80665f;
inline constexpr float kDefaultGlacialErosionStrength = 0.35f;
inline constexpr uint32_t kDefaultGlacialErosionPeriod = 60;
inline constexpr float kDefaultCollisionUpliftRatio = 0.02f;
inline constexpr float kDefaultContinentalCompressionGain = 0.05f;
inline constexpr float kDefaultContinentalBoundaryFluidity = 0.95f;
inline constexpr int32_t kAutoCrustTypeBoundaryMeters = -1;
inline constexpr uint16_t kDefaultCrustTypeBoundaryOffsetMeters = 1200;
inline constexpr float kDefaultMovementEnergy = 1.0f;
inline constexpr uint32_t kDefaultPeriodicNoisePeriod = 0;
inline constexpr float kDefaultPeriodicNoiseStrength = 0.0f;

struct Scenario {
    uint32_t version = kScenarioVersion;
    long seed = 0;
    uint32_t width = 64;
    uint32_t height = 32;
    float sea_level = 0.65f;
    uint32_t erosion_period = 0;
    float folding_ratio = 0.02f;
    uint32_t aggregation_overlap_abs = 1000000;
    float aggregation_overlap_rel = 0.33f;
    uint32_t cycle_count = 2;
    uint32_t plate_count = 4;
    float erosion_strength = 1.0f;
    float crust_rotation_strength = 0.20f;
    float rotation_strength = 1.0f;
    float subduction_strength = 1.0f;
    int32_t sea_level_m = TopographyCodec::kNoSeaLevelOverride;
    int32_t crust_type_boundary_m = kAutoCrustTypeBoundaryMeters;
    uint16_t initial_min_height_m = TopographyCodec::kDefaultInitialMinHeightMeters;
    uint16_t initial_max_height_m = TopographyCodec::kDefaultInitialMaxHeightMeters;
    uint32_t cycle_step_limit = 600;
    double cycle_duration_myr = 0.0;
    float divergent_carve_strength = 0.015f;
    double delta_time_myr = kDefaultDeltaTimeMyr;
    float gravity_mps2 = kDefaultGravityMps2;
    float glacial_erosion_strength = kDefaultGlacialErosionStrength;
    uint32_t glacial_erosion_period = kDefaultGlacialErosionPeriod;
    float collision_uplift_ratio = kDefaultCollisionUpliftRatio;
    float continental_compression_gain = kDefaultContinentalCompressionGain;
    float continental_boundary_fluidity = kDefaultContinentalBoundaryFluidity;
    float movement_energy = kDefaultMovementEnergy;
    uint32_t hf_noise_period = kDefaultPeriodicNoisePeriod;
    float hf_noise_strength = kDefaultPeriodicNoiseStrength;
    uint32_t lf_noise_period = kDefaultPeriodicNoisePeriod;
    float lf_noise_strength = kDefaultPeriodicNoiseStrength;
    bool use_material_map_for_height_limit = false;
};

inline double steps_to_myr(uint32_t steps, double delta_time_myr) {
    return static_cast<double>(steps) * delta_time_myr;
}

} // namespace platec::scenario

#endif
