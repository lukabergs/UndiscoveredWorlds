#ifndef TECTONIC_SCENARIO_HPP
#define TECTONIC_SCENARIO_HPP

#include "simulations/geology/plate_tectonics/topography_codec.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include "simulations/geology/plate_tectonics/tectonic_plumes.hpp"
#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"

#include <cstdint>

namespace tectonics::scenario {

inline constexpr uint32_t kScenarioVersion = 3;
inline constexpr double kDefaultDeltaTimeMyr = 1.0;
inline constexpr float kDefaultGravityMps2 = 9.80665f;
inline constexpr int32_t kAutoCrustTypeBoundaryMeters = -1;
inline constexpr uint16_t kDefaultCrustTypeBoundaryOffsetMeters = 1200;

struct Scenario {
    uint32_t version = kScenarioVersion;
    long seed = 0;
    uint32_t width = 64;
    uint32_t height = 32;
    float sea_level = 0.65f;
    uint32_t plate_count = 4;
    float erosion_strength = 1.0f;
    float thermal_subsidence_strength = 1.0f;
    int32_t sea_level_m = TopographyCodec::kNoSeaLevelOverride;
    int32_t crust_type_boundary_m = kAutoCrustTypeBoundaryMeters;
    // Positive magnitudes relative to sea level, independent of the uint16 datum.
    uint16_t initial_ocean_depth_m = TopographyCodec::kDefaultInitialOceanDepthMeters;
    uint16_t initial_land_relief_m = TopographyCodec::kDefaultInitialLandReliefMeters;
    // Finite positive duration of continuous material evolution.
    double evolution_duration_myr = 600.0;
    bool boundary_conditioned_terrain = false;
    double initial_boundary_influence = 1; // Local kinematic conditioning without changing the seeded column conversion.
    double initial_margin_width_m = 0; // Smooth continental/oceanic mixture around -2500 m; preserves seeded elevation.
    // Share of the surface seeded as continental crust (0: the legacy seeding, where
    // every column above -2500 m is continental). Above 0 the ocean side between the
    // land and this share is a shelf (0 to -200 m) and a slope to -2500 m; beyond it the
    // deep floor. At least 1 - sea_level; requires boundary-conditioned terrain.
    double continental_fraction = 0;
    bool sediment_load_feedback = false; // Include continental sediment in GPE and weak-crust redistribution.
    bool thermal_force_feedback = false; // Slab buoyancy follows advected thermal age, not chronological age.
    double world_width_km = 40000.0;
    double mantle_flow_speed_cm_yr = 4.0;
    uint32_t transport_subdivision = 2;
    double deformation_peak_km = 300.0;
    double deformation_variation = .2;
    double subduction_coupling = .2;
    double coastal_gpe_strength = 0; // Existing root-pressure resistance at ocean/continent contacts.
    double boundary_force_ramp_cm_yr = 0; // Smooth ridge/slab activation near zero normal motion; 0 selects instantaneous activation.
    bool consistent_motion = false; // Solve activity and velocity together; share subduction polarity.
    double boundary_memory_myr = 30; // Reduced attached-slab/ridge activity response, not mechanical inertia.
    double mountain_force_range_km = 0; // Inland GPE transmission strip; 0 samples contact cells only.
    double ocean_crust_variation_m = 0;
    double ocean_noise_wavelength_km = 600;
    double seamount_crust_m = 0;
    double seamount_radius_km = 350;
    double seamount_spacing_km = 3500;
    double seamount_roughness = 0;
    double submarine_erosion_ratio = 1;
    margins::Parameters margin_processes;
    topology::Parameters topology_processes;
    noise_pass::Parameters experimental_noise;
    plumes::Parameters plume_processes; // count 0 (default) leaves every run unchanged
    double delta_time_myr = kDefaultDeltaTimeMyr;
    // E-folding time of the exported rock uplift rate: an exponential mean of the
    // elevation change of the tectonic sub-steps, before the surface processes.
    double uplift_rate_window_myr = 10;
    float gravity_mps2 = kDefaultGravityMps2;
};

inline double steps_to_myr(uint32_t steps, double delta_time_myr) {
    return static_cast<double>(steps) * delta_time_myr;
}

} // namespace tectonics::scenario

#endif
