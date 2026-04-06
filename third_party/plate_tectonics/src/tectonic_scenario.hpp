#ifndef TECTONIC_SCENARIO_HPP
#define TECTONIC_SCENARIO_HPP

#include "topography_codec.hpp"

#include <cstdint>

namespace platec::scenario {

inline constexpr uint32_t kScenarioVersion = 1;
inline constexpr double kDefaultDeltaTimeMyr = 1.0;

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
    uint16_t initial_min_height_m = TopographyCodec::kDefaultInitialMinHeightMeters;
    uint16_t initial_max_height_m = TopographyCodec::kDefaultInitialMaxHeightMeters;
    uint32_t cycle_step_limit = 600;
    float divergent_carve_strength = 0.015f;
    double delta_time_myr = kDefaultDeltaTimeMyr;
};

inline double steps_to_myr(uint32_t steps, double delta_time_myr) {
    return static_cast<double>(steps) * delta_time_myr;
}

} // namespace platec::scenario

#endif
