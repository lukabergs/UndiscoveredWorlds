/******************************************************************************
 *  plate-tectonics, a plate tectonics simulation library
 *  Copyright (C) 2012-2013 Lauri Viitanen
 *  Copyright (C) 2014-2015 Federico Tomassetti, Bret Curtis
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, see http://www.gnu.org/licenses/
 *****************************************************************************/

#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/plate_tectonics/noise.hpp"
#include "simulations/geology/plate_tectonics/plate.hpp"
#include "simulations/geology/plate_tectonics/simplexnoise.hpp"
#include "simulations/geology/plate_tectonics/sqrdmd.hpp"
#include "simulations/geology/plate_tectonics/tectonic_boundary_graph.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deforming_regions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junction_coupling.hpp"
#include "simulations/geology/plate_tectonics/tectonic_junctions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_deformation.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>


using namespace std;

static constexpr uint32_t MAX_BUOYANCY_AGE = 20;




uint32_t findBound(const uint32_t* map, uint32_t length, uint32_t x0, uint32_t y0, int dx, int dy);
uint32_t findPlate(plate** plates, float x, float y, uint32_t num_plates);

namespace {

struct FractalNoiseConfig {
    float octaves;
    float persistence;
    float scale;
    float noise_scale;
    float offset_a;
    float offset_b;
    float offset_c;
    float offset_d;
};


float wrap_coordinate(float value, float period) {
    float wrapped = std::fmod(value, period);
    if (wrapped < 0.0f) {
        wrapped += period;
    }
    return wrapped;
}

uint32_t wrap_index(int value, uint32_t period) {
    const int period_i = static_cast<int>(period);
    int wrapped = value % period_i;
    if (wrapped < 0) {
        wrapped += period_i;
    }
    return static_cast<uint32_t>(wrapped);
}

float wrapped_delta(float value, float reference, float period) {
    float delta = value - reference;
    if (delta > period * 0.5f) {
        delta -= period;
    } else if (delta < -period * 0.5f) {
        delta += period;
    }
    return delta;
}

FractalNoiseConfig make_noise_config(SimpleRandom& randsource, float octaves, float persistence,
                                     float scale, float noise_scale) {
    const uint32_t seed = randsource.next();
    return FractalNoiseConfig {
        octaves,
        persistence,
        scale,
        noise_scale,
        17.0f + static_cast<float>((seed * 13U) % 1024U) * 0.03125f,
        29.0f + static_cast<float>((seed * 29U) % 1024U) * 0.03125f,
        43.0f + static_cast<float>((seed * 43U) % 1024U) * 0.03125f,
        61.0f + static_cast<float>((seed * 61U) % 1024U) * 0.03125f,
    };
}

FractalNoiseConfig make_noise_config(uint32_t seed, float octaves, float persistence,
                                     float scale, float noise_scale) {
    SimpleRandom randsource(seed);
    return make_noise_config(randsource, octaves, persistence, scale, noise_scale);
}

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
}

float interpolate(float start, float end, float alpha) {
    return start + (end - start) * clamp_unit(alpha);
}

float smoothstep(float edge0, float edge1, float value) {
    if (std::fabs(edge1 - edge0) <= FLT_EPSILON) {
        return value >= edge1 ? 1.0f : 0.0f;
    }
    const float t = clamp_unit((value - edge0) / (edge1 - edge0));
    return t * t * (3.0f - 2.0f * t);
}

float sample_wrapped_octave_noise(float x, float y, const WorldDimension& dimension,
                                  float octaves, float persistence, float scale,
                                  float noise_scale, float lo_bound, float hi_bound,
                                  float offset_a, float offset_b, float offset_c,
                                  float offset_d) {
    const float width = static_cast<float>(dimension.getWidth());
    const float height = static_cast<float>(dimension.getHeight());
    const float wrapped_x = wrap_coordinate(x, width) / width;
    const float wrapped_y = wrap_coordinate(y, height) / height;
    const float theta = wrapped_x * 2.0f * PI;
    const float phi = wrapped_y * 2.0f * PI;
    const float sin_theta = std::sin(theta);
    const float cos_theta = std::cos(theta);
    const float sin_phi = std::sin(phi);
    const float cos_phi = std::cos(phi);

    return scaled_octave_noise_4d(octaves, persistence, scale, lo_bound, hi_bound,
                                  offset_a + sin_theta * noise_scale,
                                  offset_b + cos_theta * noise_scale,
                                  offset_c + sin_phi * noise_scale,
                                  offset_d + cos_phi * noise_scale);
}

float sample_toroidal_noise(float x, float y, const WorldDimension& dimension,
                            const FractalNoiseConfig& config) {
    return sample_wrapped_octave_noise(x, y, dimension, config.octaves, config.persistence,
                                       config.scale, config.noise_scale, -1.0f, 1.0f,
                                       config.offset_a, config.offset_b, config.offset_c,
                                       config.offset_d);
}

float logistic(float value, float sharpness) {
    return 1.0f / (1.0f + std::exp(-value * sharpness));
}

float remap_to_range(float value, float min_value, float max_value, float target_min,
                     float target_max) {
    if (max_value - min_value <= FLT_EPSILON) {
        return 0.5f * (target_min + target_max);
    }
    return target_min +
           ((value - min_value) / (max_value - min_value)) * (target_max - target_min);
}

struct GridStep {
    int dx;
    int dy;
};

GridStep make_grid_step(const Tectonics::FloatVector& velocity) {
    const float abs_x = std::fabs(velocity.x());
    const float abs_y = std::fabs(velocity.y());
    GridStep step {0, 0};
    if (abs_x >= 0.20f) {
        step.dx = velocity.x() >= 0.0f ? 1 : -1;
    }
    if (abs_y >= 0.20f) {
        step.dy = velocity.y() >= 0.0f ? 1 : -1;
    }
    if (step.dx == 0 && step.dy == 0) {
        if (abs_x >= abs_y) {
            step.dx = velocity.x() >= 0.0f ? 1 : -1;
        } else {
            step.dy = velocity.y() >= 0.0f ? 1 : -1;
        }
    }
    return step;
}

void advance_wrapped(uint32_t& x, uint32_t& y, int dx, int dy, const WorldDimension& dimension) {
    x = wrap_index(static_cast<int>(x) + dx, dimension.getWidth());
    y = wrap_index(static_cast<int>(y) + dy, dimension.getHeight());
}

float internal_to_physical_meters(float internal, uint16_t sea_level_m) {
    if (internal <= TopographyCodec::kOceanicBase) {
        return 0.0f;
    }

    if (TopographyCodec::is_oceanic_internal(internal)) {
        if (sea_level_m <= 1U) {
            return 0.0f;
        }

        const float oceanic_ceiling =
            std::nextafter(TopographyCodec::kContinentalBase, TopographyCodec::kOceanicBase);
        const float ratio =
            (std::min(oceanic_ceiling, internal) - TopographyCodec::kOceanicBase) /
            (oceanic_ceiling - TopographyCodec::kOceanicBase);
        return ratio * static_cast<float>(sea_level_m - 1U);
    }

    return static_cast<float>(sea_level_m) +
           (internal - TopographyCodec::kContinentalBase) *
               static_cast<float>(TopographyCodec::kMaxHeightMeters - sea_level_m);
}

// Inverse of the native physical conversion, without uint16 rounding or an
// upper image-format clamp. Relief processes must not use encoded columns as
// distances: their zero is the codec floor, not sea level.
float elevation_to_internal(float elevation_m, uint16_t sea_level_m) {
    if (elevation_m >= 0.0f) {
        return TopographyCodec::kContinentalBase + elevation_m /
            std::max(1.0f, static_cast<float>(TopographyCodec::kMaxHeightMeters - sea_level_m));
    }
    const float ceiling = std::nextafter(TopographyCodec::kContinentalBase,
                                        TopographyCodec::kOceanicBase);
    const float absolute_m = std::max(0.0f, elevation_m + sea_level_m);
    return TopographyCodec::kOceanicBase + std::min(1.0f,
        absolute_m / std::max(1.0f, static_cast<float>(sea_level_m) - 1.0f)) *
        (ceiling - TopographyCodec::kOceanicBase);
}

float trench_elevation(float elevation_m, float influence, float rate, double dt_myr,
                       double age_myr, float cooling_strength) {
    // Procedural relief scales, not an Earth-raster constraint. Relax toward
    // a 6 km basin / 10.5 km trench; repeated updates cannot dig toward the
    // codec floor. Existing deeper terrain is not clipped or lifted.
    // With cooling enabled the basin reference evolves from ridge depth to
    // mature bathymetry. Using the already-cooled 6 km reference for newborn
    // seafloor would count its thermal contraction a second time afterwards.
    const float basin_depth = cooling_strength > 0.0f
        ? 2500.0f + cooling_strength * static_cast<float>(tectonics::thermal::subsidence_m(age_myr + dt_myr))
        : 6000.0f;
    const float target_m = -basin_depth - 4500.0f * clamp_unit(influence);
    const float fraction = static_cast<float>(-std::expm1(-std::max(0.0f, rate) * dt_myr));
    return elevation_m - std::max(0.0f, elevation_m - target_m) * fraction;
}

float normalized_toroidal_latitude(uint32_t y, uint32_t height) {
    if (height == 0U) {
        return 0.0f;
    }

    const float equator = static_cast<float>(height) * 0.5f;
    const float sample = static_cast<float>(y) + 0.5f;
    const float delta = std::fabs(sample - equator);
    return clamp_unit(delta / std::max(1.0f, equator));
}

float snow_line_meters(uint32_t y, uint32_t height) {
    constexpr float kEquatorialSnowLineM = 4500.0f;
    constexpr float kTropicalSnowLineM = 6000.0f;
    constexpr float kTropicLatitude = 23.5f / 90.0f;

    const float latitude = normalized_toroidal_latitude(y, height);
    if (latitude <= kTropicLatitude) {
        const float tropical_alpha = smoothstep(0.0f, kTropicLatitude, latitude);
        return interpolate(kEquatorialSnowLineM, kTropicalSnowLineM, tropical_alpha);
    }

    const float polar_alpha = smoothstep(kTropicLatitude, 1.0f, latitude);
    return interpolate(kTropicalSnowLineM, 0.0f, polar_alpha);
}

float quantile_threshold(std::vector<float> values, float quantile) {
    if (values.empty()) {
        return 0.0f;
    }
    std::sort(values.begin(), values.end());
    const float clamped_quantile = std::max(0.0f, std::min(1.0f, quantile));
    const size_t index = static_cast<size_t>(
        std::floor(clamped_quantile * static_cast<float>(values.size() - 1)));
    return values[index];
}

tectonics::contract::CrustClass classify_crust_class(float total_thickness,
                                                  float continental_thickness) {
    using tectonics::contract::CrustClass;

    if (total_thickness <= FLT_EPSILON) {
        return CrustClass::None;
    }

    const float continental_ratio = continental_thickness / total_thickness;
    if (continental_ratio >= 0.70f) {
        return CrustClass::Continental;
    }
    if (continental_ratio <= 0.30f) {
        return CrustClass::Oceanic;
    }
    return CrustClass::Transitional;
}

tectonics::contract::CrustClass crust_class_from_material(uint8_t material_index) {
    using tectonics::contract::CrustClass;
    return tectonics::material::from_index(material_index) == tectonics::material::Type::Basalt
        ? CrustClass::Oceanic
        : CrustClass::Continental;
}

bool is_continental_material(uint8_t material_index) {
    return crust_class_from_material(material_index) == tectonics::contract::CrustClass::Continental;
}

uint16_t resolve_crust_type_boundary_m(const tectonics::scenario::Scenario& scenario,
                                       uint16_t sea_level_m) {
    if (scenario.crust_type_boundary_m != tectonics::scenario::kAutoCrustTypeBoundaryMeters) {
        return static_cast<uint16_t>(scenario.crust_type_boundary_m);
    }

    const uint32_t threshold =
        static_cast<uint32_t>(sea_level_m) +
        static_cast<uint32_t>(tectonics::scenario::kDefaultCrustTypeBoundaryOffsetMeters);
    return static_cast<uint16_t>(std::min<uint32_t>(TopographyCodec::kMaxHeightMeters, threshold));
}

} // namespace

static void push_unique_owner(std::vector<uint32_t>& candidates, uint32_t owner, uint32_t plate_count) {
    if (owner >= plate_count) {
        return;
    }
    for (uint32_t existing : candidates) {
        if (existing == owner) {
            return;
        }
    }
    candidates.push_back(owner);
}

WorldPoint lithosphere::randomPosition() {
    return WorldPoint(_randsource.next() % _worldDimension.getWidth(),
                      _randsource.next() % _worldDimension.getHeight(), _worldDimension);
}

void lithosphere::createNoise(float* tmp, const WorldDimension& tmpDim, bool useSimplex) {
    ::createNoise(tmp, tmpDim, _randsource, useSimplex);
}

void lithosphere::createSlowNoise(float* tmp, const WorldDimension& tmpDim) {
    ::createSlowNoise(tmp, tmpDim, _randsource);
}

lithosphere::lithosphere(const tectonics::scenario::Scenario& params) noexcept(false)
    : scenario(params), hmap(params.width, params.height), prev_hmap(params.width, params.height),
      display_hmap(params.width, params.height), prev_display_hmap(params.width, params.height),
      initial_hmap(params.width, params.height), imap(params.width, params.height),
      prev_imap(params.width, params.height), amap(params.width, params.height), plates(nullptr),
      plate_areas(params.plate_count), plate_indices_found(params.plate_count),
      iter_count(0), time_origin_step(0), max_plates(params.plate_count),
      plate_capacity(params.plate_count), num_plates(0),
      erosion_strength(std::max(0.0f, params.erosion_strength)),
      sea_level_m(TopographyCodec::default_import_sea_level_m()),
      initial_ocean_depth_m(params.initial_ocean_depth_m),
      initial_land_relief_m(params.initial_land_relief_m),
      convergence_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      divergence_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      shear_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      geologic_regime_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          static_cast<uint8_t>(tectonics::contract::GeologicRegime::Stable)),
      crust_age_myr_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      crust_thickness_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          0.0f),
      crust_class_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                      static_cast<uint8_t>(tectonics::contract::CrustClass::None)),
      uplift_tendency_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          0.0f),
      rock_uplift_rate_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                           0.0f),
      subsidence_tendency_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      accumulated_strain_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      boundary_type_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                        static_cast<uint8_t>(tectonics::contract::BoundaryType::None)),
      boundary_distance_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      boundary_segment_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          tectonics::contract::kNoBoundaryId),
      nearest_boundary_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          tectonics::contract::kNoBoundaryId),
      deforming_region_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          tectonics::contract::kNoDeformingRegionId),
      deforming_region_type_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None)),
      deformation_rate_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      deformation_velocity_x_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      deformation_velocity_y_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      junction_rift_influence_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      junction_transform_influence_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      junction_subduction_influence_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      surface_crust_type_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)),
      material_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                   tectonics::material::to_index(tectonics::material::kDefaultType)),
      continental_lock_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      rift_tension_x_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      rift_tension_y_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      rift_tension_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      _worldDimension(params.width, params.height), _randsource(params.seed), _steps(0) {
    if (params.width < 5 || params.height < 5) {
        throw runtime_error("Width and height should be >=5");
    }
    if (!std::isfinite(params.delta_time_myr) || params.delta_time_myr <= 0.0) {
        throw runtime_error("delta_time_myr must be finite and > 0");
    }
    if (!std::isfinite(params.uplift_rate_window_myr) || params.uplift_rate_window_myr <= 0.0) {
        throw runtime_error("uplift_rate_window_myr must be finite and > 0");
    }
    if (!std::isfinite(params.evolution_duration_myr) || params.evolution_duration_myr <= 0.0 ||
        params.evolution_duration_myr / params.delta_time_myr >
            static_cast<double>(std::numeric_limits<int>::max() - 10000)) {
        throw runtime_error("evolution_duration_myr must be finite, positive and fit the update counter");
    }
    if (!std::isfinite(params.world_width_km) || params.world_width_km <= 0.0 ||
        !std::isfinite(params.mantle_flow_speed_cm_yr) || params.mantle_flow_speed_cm_yr < 0.0) {
        throw runtime_error("world_width_km must be positive and mantle_flow_speed_cm_yr nonnegative and finite");
    }
    tectonics::ocean_relief::validate({params.ocean_crust_variation_m, params.ocean_noise_wavelength_km,
        params.seamount_crust_m,params.seamount_radius_km,params.seamount_spacing_km,params.seamount_roughness});
    if (!std::isfinite(params.submarine_erosion_ratio) || params.submarine_erosion_ratio<0 || params.submarine_erosion_ratio>1)
        throw runtime_error("Submarine erosion ratio must be in 0..1");
    tectonics::margins::validate(params.margin_processes);
    tectonics::noise_pass::validate(params.experimental_noise);
    tectonics::plumes::validate(params.plume_processes);
    if (params.plume_processes.weakening>0 && !params.topology_processes.progressive_rifting)
        throw runtime_error("Plume weakening acts through progressive rifting; enable it or set weakening to 0");
    if (!std::isfinite(params.mountain_force_range_km) || params.mountain_force_range_km<0)
        throw runtime_error("Inland mountain force requires isostatic columns and a nonnegative range");
    if (!std::isfinite(params.initial_margin_width_m) || params.initial_margin_width_m<0 || params.initial_margin_width_m>4000 ||
        (params.initial_margin_width_m>0 && !params.boundary_conditioned_terrain))
        throw runtime_error("Initial margin width requires boundary-conditioned terrain and 0..4000 m");
    if (params.margin_processes.process_variation>0 && !params.margin_processes.enabled)
        throw runtime_error("Process variation requires isostatic columns and active margins");
    if ((params.transport_subdivision != 1 && params.transport_subdivision != 2 && params.transport_subdivision != 4) ||
        !std::isfinite(params.deformation_peak_km) || params.deformation_peak_km < 0 ||
        !std::isfinite(params.deformation_variation) || params.deformation_variation < 0 || params.deformation_variation > .5 ||
        !std::isfinite(params.subduction_coupling) || params.subduction_coupling < 0 || params.subduction_coupling > 1 ||
        !std::isfinite(params.coastal_gpe_strength) || params.coastal_gpe_strength < 0 || params.coastal_gpe_strength > 1)
        throw runtime_error("Invalid transport subdivision or distributed deformation parameters");
    if (!std::isfinite(params.thermal_subsidence_strength) || params.thermal_subsidence_strength < 0.0f) {
        throw runtime_error("thermal_subsidence_strength must be finite and >= 0");
    }
    if (!std::isfinite(params.boundary_force_ramp_cm_yr) || params.boundary_force_ramp_cm_yr < 0)
        throw runtime_error("Boundary force activation speed must be finite and nonnegative");
    if (params.consistent_motion && params.boundary_force_ramp_cm_yr<=0)
        throw runtime_error("Consistent motion requires force balance, isostatic columns and a positive force ramp");
    if (!std::isfinite(params.boundary_memory_myr) || params.boundary_memory_myr<=0)
        throw runtime_error("Boundary memory must be finite and positive");
    if (params.crust_type_boundary_m < tectonics::scenario::kAutoCrustTypeBoundaryMeters ||
        params.crust_type_boundary_m > static_cast<int32_t>(TopographyCodec::kMaxHeightMeters)) {
        throw runtime_error("crust_type_boundary_m is out of range");
    }


    // Create default plates
    plates = new plate*[plate_capacity]();
    seedInitialTopography(params.sea_level, params.sea_level_m);
    seedInitialCrustTypes();
    seedMaterialMap();
    createPlates();
}

lithosphere::~lithosphere() throw() {
    clearPlates();
    delete[] plates;
    plates = 0;
}

void lithosphere::clearPlates() {
    for (uint32_t i = 0; i < num_plates; i++) {
        delete plates[i];
    }
    num_plates = 0;
}

void lithosphere::ensurePlateCapacity(uint32_t required_capacity) {
    if (required_capacity <= plate_capacity) {
        return;
    }

    uint32_t new_capacity = std::max(required_capacity, plate_capacity + std::max(4U, plate_capacity));
    plate** expanded = new plate*[new_capacity]();
    for (uint32_t i = 0; i < num_plates; ++i) {
        expanded[i] = plates[i];
    }

    delete[] plates;
    plates = expanded;
    plate_capacity = new_capacity;
    plate_areas.resize(plate_capacity);
    plate_indices_found.resize(plate_capacity);


}

void lithosphere::growPlates() {
    // "Grow" plates from their origins until surface is fully populated.
    uint32_t max_border = 1;
    uint32_t i;
    while (max_border) {
        for (max_border = i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];
            const uint32_t N = (uint32_t)area.border.size();
            max_border = max_border > N ? max_border : N;

            if (N == 0) {
                continue;
            }
            const uint32_t j = _randsource.next() % N;
            const uint32_t p = area.border[j];
            const uint32_t cy = _worldDimension.yFromIndex(p);
            const uint32_t cx = _worldDimension.xFromIndex(p);

            const uint32_t lft = cx > 0 ? cx - 1 : _worldDimension.getWidth() - 1;
            const uint32_t rgt = cx < _worldDimension.getWidth() - 1 ? cx + 1 : 0;
            const uint32_t top = cy > 0 ? cy - 1 : _worldDimension.getHeight() - 1;
            const uint32_t btm = cy < _worldDimension.getHeight() - 1 ? cy + 1 : 0;

            const uint32_t n = top * _worldDimension.getWidth() + cx; // North.
            const uint32_t s = btm * _worldDimension.getWidth() + cx; // South.
            const uint32_t w = cy * _worldDimension.getWidth() + lft; // West.
            const uint32_t e = cy * _worldDimension.getWidth() + rgt; // East.

            if (imap[n] >= num_plates) {
                imap[n] = i;
                area.border.push_back(n);

                if (area.top == _worldDimension.yMod(top + 1)) {
                    area.top = top;
                    area.hgt++;
                }
            }

            if (imap[s] >= num_plates) {
                imap[s] = i;
                area.border.push_back(s);

                if (btm == _worldDimension.yMod(area.btm + 1)) {
                    area.btm = btm;
                    area.hgt++;
                }
            }

            if (imap[w] >= num_plates) {
                imap[w] = i;
                area.border.push_back(w);

                if (area.lft == _worldDimension.xMod(lft + 1)) {
                    area.lft = lft;
                    area.wdt++;
                }
            }

            if (imap[e] >= num_plates) {
                imap[e] = i;
                area.border.push_back(e);

                if (rgt == _worldDimension.xMod(area.rgt + 1)) {
                    area.rgt = rgt;
                    area.wdt++;
                }
            }

            // Overwrite processed point with unprocessed one.
            area.border[j] = area.border.back();
            area.border.pop_back();
        }
    }
}

void lithosphere::rebuildPlateAreasFromOwnership() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        plateArea& area = plate_areas[i];
        area.border.clear();
        area.lft = _worldDimension.getWidth();
        area.rgt = 0;
        area.top = _worldDimension.getHeight();
        area.btm = 0;
        area.wdt = 0;
        area.hgt = 0;
    }

    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    std::vector<uint8_t> touches_boundary(map_area, 0);

    for (uint32_t index = 0; index < map_area; ++index) {
        const uint32_t owner = imap[index];
        ASSERT(owner < num_plates, "A point was not assigned to any plate");

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        plateArea& area = plate_areas[owner];

        area.lft = std::min(area.lft, x);
        area.rgt = std::max(area.rgt, x);
        area.top = std::min(area.top, y);
        area.btm = std::max(area.btm, y);

        const uint32_t left_x = x > 0 ? x - 1 : width - 1;
        const uint32_t right_x = x + 1 < width ? x + 1 : 0;
        const uint32_t top_y = y > 0 ? y - 1 : height - 1;
        const uint32_t bottom_y = y + 1 < height ? y + 1 : 0;
        const uint32_t neighbors[] = {
            _worldDimension.indexOf(left_x, y),
            _worldDimension.indexOf(right_x, y),
            _worldDimension.indexOf(x, top_y),
            _worldDimension.indexOf(x, bottom_y),
        };

        for (uint32_t neighbor : neighbors) {
            if (imap[neighbor] != owner) {
                touches_boundary[index] = 1;
                break;
            }
        }
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (touches_boundary[index]) {
            plate_areas[imap[index]].border.push_back(index);
        }
    }

    for (uint32_t i = 0; i < num_plates; ++i) {
        plateArea& area = plate_areas[i];
        area.wdt = area.rgt >= area.lft ? area.rgt - area.lft + 1U : 1U;
        area.hgt = area.btm >= area.top ? area.btm - area.top + 1U : 1U;
    }
}

void lithosphere::jaggedizePlateBoundaries() {
    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    std::vector<uint32_t> counts(num_plates, 0);
    for (uint32_t index = 0; index < map_area; ++index) {
        ASSERT(imap[index] < num_plates, "A point was not assigned to any plate");
        ++counts[imap[index]];
    }

    IndexMap current = imap;
    for (uint32_t pass = 0; pass < 2; ++pass) {
        IndexMap next = current;
        std::vector<uint32_t> next_counts = counts;

        for (uint32_t y = 0; y < height; ++y) {
            for (uint32_t x = 0; x < width; ++x) {
                const uint32_t index = _worldDimension.indexOf(x, y);
                const uint32_t owner = current[index];
                if (next_counts[owner] <= 1U) {
                    continue;
                }

                const uint32_t left_x = x > 0 ? x - 1 : width - 1;
                const uint32_t right_x = x + 1 < width ? x + 1 : 0;
                const uint32_t top_y = y > 0 ? y - 1 : height - 1;
                const uint32_t bottom_y = y + 1 < height ? y + 1 : 0;
                struct NeighborWeight {
                    uint32_t index;
                    float weight;
                };
                const NeighborWeight neighbors[] = {
                    {_worldDimension.indexOf(left_x, y), 1.0f},
                    {_worldDimension.indexOf(right_x, y), 1.0f},
                    {_worldDimension.indexOf(x, top_y), 1.0f},
                    {_worldDimension.indexOf(x, bottom_y), 1.0f},
                    {_worldDimension.indexOf(left_x, top_y), 0.7f},
                    {_worldDimension.indexOf(right_x, top_y), 0.7f},
                    {_worldDimension.indexOf(left_x, bottom_y), 0.7f},
                    {_worldDimension.indexOf(right_x, bottom_y), 0.7f},
                };

                bool touches_other_owner = false;
                std::vector<uint32_t> candidates;
                candidates.reserve(6);
                push_unique_owner(candidates, owner, num_plates);
                for (const NeighborWeight& neighbor : neighbors) {
                    const uint32_t neighbor_owner = current[neighbor.index];
                    touches_other_owner = touches_other_owner || neighbor_owner != owner;
                    push_unique_owner(candidates, neighbor_owner, num_plates);
                }
                if (!touches_other_owner) {
                    continue;
                }

                uint32_t best_owner = owner;
                float best_score = -FLT_MAX;
                for (uint32_t candidate : candidates) {
                    float score = candidate == owner ? 1.35f : 0.0f;
                    for (const NeighborWeight& neighbor : neighbors) {
                        if (current[neighbor.index] == candidate) {
                            score += neighbor.weight;
                        }
                    }

                    const float owner_bias = sample_wrapped_octave_noise(
                        static_cast<float>(x), static_cast<float>(y), _worldDimension, 4.0f,
                        0.57f, 0.10f, 1.9f, -1.0f, 1.0f,
                        static_cast<float>(candidate) * 13.0f + 17.0f,
                        static_cast<float>(candidate) * 7.0f + 29.0f,
                        static_cast<float>(pass) * 0.37f,
                        static_cast<float>(candidate) * 19.0f + 11.0f);
                    const float warp_bias = sample_wrapped_octave_noise(
                        static_cast<float>(x), static_cast<float>(y), _worldDimension, 2.0f,
                        0.50f, 0.23f, 2.5f, -1.0f, 1.0f, 5.0f, 13.0f,
                        static_cast<float>(candidate) * 0.41f,
                        static_cast<float>(pass) * 3.0f + 43.0f);
                    score += owner_bias * 1.15f + warp_bias * 0.75f;

                    if (score > best_score) {
                        best_score = score;
                        best_owner = candidate;
                    }
                }

                if (best_owner != owner && next_counts[owner] > 1U) {
                    next[index] = best_owner;
                    --next_counts[owner];
                    ++next_counts[best_owner];
                }
            }
        }

        current = next;
        counts.swap(next_counts);
    }

    imap = current;
    rebuildPlateAreasFromOwnership();
}

void lithosphere::createPlates() {
    try {
        const uint32_t map_area = _worldDimension.getArea();
        num_plates = max_plates;
        _steps = 0;
        if (num_plates > map_area) {
            throw runtime_error("plate_count cannot exceed map area");
        }
        for (uint32_t i = 0; i < num_plates; ++i) {
            plate_areas[i].border.clear();
            plate_areas[i].border.reserve(8);
        }

        // Initialize "Free plate center position" lookup table.
        // This way two plate centers will never be identical.
        std::vector<uint32_t> free_centers(map_area);
        for (uint32_t i = 0; i < map_area; ++i)
            free_centers[i] = i;

        // Select N plate centers from the global map.

        for (uint32_t i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];

            // Randomly select an unused plate origin.
            const uint32_t remaining = map_area - i;
            const uint32_t random_index = (uint32_t)_randsource.next() % remaining;
            const uint32_t p = free_centers[random_index];
            const uint32_t y = _worldDimension.yFromIndex(p);
            const uint32_t x = _worldDimension.xFromIndex(p);

            area.lft = area.rgt = x; // Save origin...
            area.top = area.btm = y;
            area.wdt = area.hgt = 1;

            area.border.clear();
            area.border.push_back(p); // ...and mark it as border.

            // Overwrite used entry with last unused entry in array.
            free_centers[random_index] = free_centers[remaining - 1U];
        }

        imap.set_all(0xFFFFFFFF);
        growPlates();

        // check all the points of the map are owned
        for (uint32_t i = 0; i < map_area; i++) {
            ASSERT(imap[i] < num_plates, "A point was not assigned to any plate");
        }

        // Extract and create plates from initial terrain.
        for (uint32_t i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];

            area.wdt = _worldDimension.xCap(area.wdt);
            area.hgt = _worldDimension.yCap(area.hgt);

            const uint32_t x0 = area.lft;
            const uint32_t x1 = 1 + x0 + area.wdt;
            const uint32_t y0 = area.top;
            const uint32_t y1 = 1 + y0 + area.hgt;
            const uint32_t width = x1 - x0;
            const uint32_t height = y1 - y0;
            float* pmap = new float[width * height];
            std::vector<uint8_t> pmat(
                static_cast<size_t>(width) * static_cast<size_t>(height),
                tectonics::material::to_index(tectonics::material::kDefaultType));

            // Copy plate's height data from global map into local map.
            for (uint32_t y = y0, j = 0; y < y1; ++y) {
                for (uint32_t x = x0; x < x1; ++x, ++j) {
                    uint32_t k = _worldDimension.normalizedIndexOf(x, y);
                    pmap[j] = hmap[k] * (imap[k] == i);
                    if (imap[k] == i) {
                        pmat[j] = material_map[k];
                    }
                }
            }
            // Create plate.
            // MK: The pmap array becomes owned by map, do not delete it
            plates[i] = new plate(_randsource.next(), pmap, width, height, x0, y0, i,
                                  _worldDimension, pmat.data(), 1.0f);
        }

        if (iter_count == 0U) {
            const double age_horizon_steps =
                std::ceil(260.0 / std::max(0.001, scenario.delta_time_myr));
            iter_count = std::max<uint32_t>(
                num_plates + MAX_BUOYANCY_AGE,
                static_cast<uint32_t>(std::min(
                    age_horizon_steps + static_cast<double>(MAX_BUOYANCY_AGE),
                    static_cast<double>(std::numeric_limits<uint32_t>::max() - 1U))));
        }
        const bool initial_creation = time_origin_step == 0U;
        if (initial_creation) {
            time_origin_step = iter_count;
        }


        updateDerivedMaps();
        // Repartitioning plates is not the creation of new crust. Restarts
        // restore these birth timestamps into the new owners below.
        if (initial_creation) {
            seedInitialCrustAges();
        }
        if (initial_creation) initializeColumns();
        updateDerivedMaps();

    } catch (const exception& e) {
        string msg = "Problem during createPlates: ";
        msg = msg + e.what();
        throw runtime_error(msg.c_str());
    }
}

uint32_t lithosphere::getPlateCount() const throw() {
    return num_plates;
}

uint32_t lithosphere::getProvenancePlateCount() const noexcept {
    return static_cast<uint32_t>(plate_kinematics.size());
}

const uint32_t* lithosphere::getAgeMap() const throw() {
    return amap.raw_data();
}

const float* lithosphere::getCrustAgeMyrMap() const noexcept {
    return crust_age_myr_map.data();
}

const float* lithosphere::getCrustThicknessMap() const noexcept {
    return crust_thickness_map.data();
}

const uint8_t* lithosphere::getCrustClassMap() const noexcept {
    return crust_class_map.data();
}

const uint8_t* lithosphere::getSurfaceCrustTypeMap() const noexcept {
    return surface_crust_type_map.data();
}

const float* lithosphere::getUpliftTendencyMap() const noexcept {
    return uplift_tendency_map.data();
}

const float* lithosphere::getRockUpliftRateMap() const noexcept {
    return rock_uplift_rate_map.data();
}

const float* lithosphere::getSubsidenceTendencyMap() const noexcept {
    return subsidence_tendency_map.data();
}

const float* lithosphere::getAccumulatedStrainMap() const noexcept {
    return accumulated_strain_map.data();
}

const float* lithosphere::getTensionMap() const noexcept {
    return rift_tension_map.data();
}

const uint8_t* lithosphere::getBoundaryTypeMap() const noexcept {
    return boundary_type_map.data();
}

const uint16_t* lithosphere::getBoundaryDistanceMap() const noexcept {
    return boundary_distance_map.data();
}

const uint32_t* lithosphere::getBoundarySegmentIdMap() const noexcept {
    return boundary_segment_id_map.data();
}

const uint32_t* lithosphere::getNearestBoundaryIdMap() const noexcept {
    return nearest_boundary_id_map.data();
}

const tectonics::contract::BoundarySegment* lithosphere::getBoundarySegments() const noexcept {
    return boundary_segments.data();
}

uint32_t lithosphere::getBoundarySegmentCount() const noexcept {
    return static_cast<uint32_t>(boundary_segments.size());
}

const uint32_t* lithosphere::getDeformingRegionIdMap() const noexcept {
    return deforming_region_id_map.data();
}

const uint8_t* lithosphere::getDeformingRegionTypeMap() const noexcept {
    return deforming_region_type_map.data();
}

const float* lithosphere::getDeformationRateMap() const noexcept {
    return deformation_rate_map.data();
}

const float* lithosphere::getDeformationVelocityXMap() const noexcept {
    return deformation_velocity_x_map.data();
}

const float* lithosphere::getDeformationVelocityYMap() const noexcept {
    return deformation_velocity_y_map.data();
}

const tectonics::contract::DeformingRegion* lithosphere::getDeformingRegions() const noexcept {
    return deforming_regions.data();
}

uint32_t lithosphere::getDeformingRegionCount() const noexcept {
    return static_cast<uint32_t>(deforming_regions.size());
}

const tectonics::contract::TectonicJunction* lithosphere::getTectonicJunctions() const noexcept {
    return tectonic_junctions.data();
}

uint32_t lithosphere::getTectonicJunctionCount() const noexcept {
    return static_cast<uint32_t>(tectonic_junctions.size());
}

const uint8_t* lithosphere::getMaterialMap() const noexcept {
    return material_map.data();
}

float* lithosphere::getTopography() const throw() {
    return display_hmap.raw_data();
}

const uint8_t* lithosphere::getConvergenceMap() const noexcept {
    return convergence_map.data();
}

const uint8_t* lithosphere::getDivergenceMap() const noexcept {
    return divergence_map.data();
}

const uint8_t* lithosphere::getShearMap() const noexcept {
    return shear_map.data();
}

const uint8_t* lithosphere::getGeologicRegimeMap() const noexcept {
    return geologic_regime_map.data();
}

const tectonics::contract::PlateKinematics* lithosphere::getPlateKinematics() const noexcept {
    return plate_kinematics.data();
}

double lithosphere::getTimeMyr() const noexcept {
    const uint32_t elapsed_steps =
        iter_count >= time_origin_step ? iter_count - time_origin_step : 0U;
    return tectonics::scenario::steps_to_myr(elapsed_steps, scenario.delta_time_myr);
}

bool lithosphere::isOceanic(float value) const noexcept {
    return TopographyCodec::is_oceanic_internal(value);
}

void lithosphere::resetSimulationState() {
    amap.set_all(0);
    plume_influence_map.clear();
    imap.set_all(0xFFFFFFFF);
    prev_imap.set_all(0xFFFFFFFF);
    prev_hmap = hmap;
    prev_display_hmap = display_hmap;
    std::fill(convergence_map.begin(), convergence_map.end(), static_cast<uint8_t>(0));
    std::fill(divergence_map.begin(), divergence_map.end(), static_cast<uint8_t>(0));
    std::fill(shear_map.begin(), shear_map.end(), static_cast<uint8_t>(0));
    std::fill(geologic_regime_map.begin(), geologic_regime_map.end(),
              static_cast<uint8_t>(tectonics::contract::GeologicRegime::Stable));
    std::fill(crust_age_myr_map.begin(), crust_age_myr_map.end(), 0.0f);
    std::fill(crust_thickness_map.begin(), crust_thickness_map.end(), 0.0f);
    std::fill(crust_class_map.begin(), crust_class_map.end(),
              static_cast<uint8_t>(tectonics::contract::CrustClass::None));
    std::fill(uplift_tendency_map.begin(), uplift_tendency_map.end(), 0.0f);
    std::fill(rock_uplift_rate_map.begin(), rock_uplift_rate_map.end(), 0.0f);
    std::fill(subsidence_tendency_map.begin(), subsidence_tendency_map.end(), 0.0f);
    std::fill(accumulated_strain_map.begin(), accumulated_strain_map.end(), 0.0f);
    std::fill(boundary_type_map.begin(), boundary_type_map.end(),
              static_cast<uint8_t>(tectonics::contract::BoundaryType::None));
    std::fill(boundary_distance_map.begin(), boundary_distance_map.end(), static_cast<uint16_t>(0));
    std::fill(boundary_segment_id_map.begin(), boundary_segment_id_map.end(),
              tectonics::contract::kNoBoundaryId);
    std::fill(nearest_boundary_id_map.begin(), nearest_boundary_id_map.end(),
              tectonics::contract::kNoBoundaryId);
    std::fill(deforming_region_id_map.begin(), deforming_region_id_map.end(),
              tectonics::contract::kNoDeformingRegionId);
    std::fill(deforming_region_type_map.begin(), deforming_region_type_map.end(),
              static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None));
    std::fill(deformation_rate_map.begin(), deformation_rate_map.end(), 0.0f);
    std::fill(deformation_velocity_x_map.begin(), deformation_velocity_x_map.end(), 0.0f);
    std::fill(deformation_velocity_y_map.begin(), deformation_velocity_y_map.end(), 0.0f);
    std::fill(junction_rift_influence_map.begin(), junction_rift_influence_map.end(), 0.0f);
    std::fill(junction_transform_influence_map.begin(), junction_transform_influence_map.end(), 0.0f);
    std::fill(junction_subduction_influence_map.begin(), junction_subduction_influence_map.end(), 0.0f);
    std::fill(continental_lock_map.begin(), continental_lock_map.end(), static_cast<uint8_t>(0));
    std::fill(rift_tension_x_map.begin(), rift_tension_x_map.end(), 0.0f);
    std::fill(rift_tension_y_map.begin(), rift_tension_y_map.end(), 0.0f);
    std::fill(rift_tension_map.begin(), rift_tension_map.end(), 0.0f);
    plate_kinematics.clear();
    boundary_segments.clear();
    deforming_regions.clear();
    tectonic_junctions.clear();
    next_boundary_segment_id = tectonics::contract::kNoBoundaryId + 1U;

    for (uint32_t i = 0; i < plate_capacity; ++i) {


    }

    clearPlates();

    iter_count = 0;
    time_origin_step = 0;


    _steps = 0;
    finished = false;
    motion_diagnostics = {};
    column_model = {};
}

void lithosphere::updateCrustAgeMyrMap() {
    const uint32_t map_area = _worldDimension.getArea();
    for (uint32_t i = 0; i < map_area; ++i) {
        const uint32_t age_steps = iter_count >= amap[i] ? iter_count - amap[i] : 0U;
        crust_age_myr_map[i] = static_cast<float>(
            tectonics::scenario::steps_to_myr(age_steps, scenario.delta_time_myr));
    }
}

void lithosphere::updateTectonicProvenance() {
    updateCrustAgeMyrMap();
    plate_kinematics.resize(num_plates);
    for (uint32_t plate_index = 0; plate_index < num_plates; ++plate_index) {
        const plate* current_plate = plates[plate_index];
        tectonics::contract::PlateKinematics& kinematics = plate_kinematics[plate_index];
        const Tectonics::FloatVector unit_velocity = current_plate->velocityUnitVector();
        const Tectonics::FloatVector linear_velocity = current_plate->linearVelocityVector();
        const FloatPoint mass_center = current_plate->worldMassCenter();

        kinematics.unit_velocity_x = unit_velocity.x();
        kinematics.unit_velocity_y = unit_velocity.y();
        kinematics.linear_velocity_x = linear_velocity.x();
        kinematics.linear_velocity_y = linear_velocity.y();
        kinematics.angular_velocity = current_plate->getAngularVelocity();
        kinematics.mass_center_x = mass_center.getX();
        kinematics.mass_center_y = mass_center.getY();
    }

    const tectonics::provenance::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        display_hmap.raw_data(),
        imap.raw_data(),
        plate_kinematics.data(),
        static_cast<uint32_t>(plate_kinematics.size()),
        crust_class_map.data(),
        motion_diagnostics.enabled
            ? scenario.world_width_km / (_worldDimension.getWidth() * scenario.delta_time_myr * 10.0)
            : 0.0,
    };
    tectonics::provenance::compute_maps(inputs, convergence_map.data(), divergence_map.data(),
                                     shear_map.data(), geologic_regime_map.data());
}

void lithosphere::updateBoundaryDistanceMap() {
    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    const uint16_t fallback_distance = static_cast<uint16_t>(
        std::min<uint32_t>(std::numeric_limits<uint16_t>::max(), width + height));

    std::fill(boundary_distance_map.begin(), boundary_distance_map.end(), fallback_distance);
    std::fill(nearest_boundary_id_map.begin(), nearest_boundary_id_map.end(),
              tectonics::contract::kNoBoundaryId);

    std::deque<uint32_t> frontier;
    frontier.clear();
    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != tectonics::contract::kNoBoundaryId) {
            boundary_distance_map[index] = 0U;
            nearest_boundary_id_map[index] = boundary_segment_id_map[index];
            frontier.push_back(index);
        }
    }

    if (frontier.empty()) {
        return;
    }

    while (!frontier.empty()) {
        const uint32_t index = frontier.front();
        frontier.pop_front();

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        const uint32_t left_x = x > 0 ? x - 1 : width - 1;
        const uint32_t right_x = x + 1 < width ? x + 1 : 0;
        const uint32_t top_y = y > 0 ? y - 1 : height - 1;
        const uint32_t bottom_y = y + 1 < height ? y + 1 : 0;
        const uint32_t neighbors[] = {
            _worldDimension.indexOf(left_x, y),
            _worldDimension.indexOf(right_x, y),
            _worldDimension.indexOf(x, top_y),
            _worldDimension.indexOf(x, bottom_y),
        };

        const uint16_t next_distance = static_cast<uint16_t>(boundary_distance_map[index] + 1U);
        const uint32_t boundary_id = nearest_boundary_id_map[index];
        for (uint32_t neighbor : neighbors) {
            if (next_distance < boundary_distance_map[neighbor] ||
                (next_distance == boundary_distance_map[neighbor] &&
                 boundary_id != tectonics::contract::kNoBoundaryId &&
                 (nearest_boundary_id_map[neighbor] == tectonics::contract::kNoBoundaryId ||
                  boundary_id < nearest_boundary_id_map[neighbor]))) {
                boundary_distance_map[neighbor] = next_distance;
                nearest_boundary_id_map[neighbor] = boundary_id;
                frontier.push_back(neighbor);
            }
        }
    }
}

void lithosphere::updateBoundaryGraphMaps() {
    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    const std::vector<uint8_t> raw_geologic_regime = geologic_regime_map;

    const tectonics::boundary_graph::Inputs inputs{
        width,
        height,
        imap.raw_data(),
        convergence_map.data(),
        divergence_map.data(),
        shear_map.data(),
        raw_geologic_regime.data(),
        crust_class_map.data(),
        boundary_segment_id_map.data(),
        nearest_boundary_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        next_boundary_segment_id,
        scenario.delta_time_myr,
    };
    tectonics::boundary_graph::Outputs outputs = tectonics::boundary_graph::build(inputs);
    boundary_segments = std::move(outputs.segments);
    boundary_segment_id_map = std::move(outputs.boundary_segment_id);
    boundary_type_map = std::move(outputs.boundary_type);
    geologic_regime_map = std::move(outputs.boundary_regime);
    next_boundary_segment_id = outputs.next_boundary_id;

    updateBoundaryDistanceMap();

    std::vector<const tectonics::contract::BoundarySegment*> segments_by_id(
        next_boundary_segment_id, nullptr);
    for (const tectonics::contract::BoundarySegment& segment : boundary_segments) {
        if (segment.id < segments_by_id.size()) {
            segments_by_id[segment.id] = &segment;
        }
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != tectonics::contract::kNoBoundaryId) {
            continue;
        }

        const uint32_t nearest_boundary_id = nearest_boundary_id_map[index];
        const auto raw_regime =
            static_cast<tectonics::contract::GeologicRegime>(raw_geologic_regime[index]);
        if (nearest_boundary_id == tectonics::contract::kNoBoundaryId ||
            nearest_boundary_id >= segments_by_id.size() ||
            segments_by_id[nearest_boundary_id] == nullptr) {
            geologic_regime_map[index] =
                static_cast<uint8_t>(raw_regime == tectonics::contract::GeologicRegime::PassiveMargin
                                         ? raw_regime
                                         : tectonics::contract::GeologicRegime::Stable);
            continue;
        }

        const tectonics::contract::BoundarySegment& segment =
            *segments_by_id[nearest_boundary_id];
        const uint16_t distance = boundary_distance_map[index];
        auto resolved_regime = tectonics::contract::GeologicRegime::Stable;
        if (distance <= 1U) {
            switch (segment.boundary_type) {
            case tectonics::contract::BoundaryType::Convergent:
                if (raw_regime == tectonics::contract::GeologicRegime::TrenchAdjacent ||
                    segment.geologic_regime == tectonics::contract::GeologicRegime::TrenchAdjacent) {
                    resolved_regime = tectonics::contract::GeologicRegime::TrenchAdjacent;
                } else {
                    resolved_regime = segment.geologic_regime;
                }
                break;
            case tectonics::contract::BoundaryType::Divergent:
                resolved_regime = segment.geologic_regime;
                break;
            case tectonics::contract::BoundaryType::Transform:
                resolved_regime = tectonics::contract::GeologicRegime::Transform;
                break;
            case tectonics::contract::BoundaryType::PassiveMargin:
                resolved_regime = tectonics::contract::GeologicRegime::PassiveMargin;
                break;
            case tectonics::contract::BoundaryType::None:
            default:
                resolved_regime = tectonics::contract::GeologicRegime::Stable;
                break;
            }
        } else if (distance <= 3U &&
                   segment.boundary_type == tectonics::contract::BoundaryType::PassiveMargin) {
            resolved_regime = tectonics::contract::GeologicRegime::PassiveMargin;
        } else if (distance <= 2U && raw_regime != tectonics::contract::GeologicRegime::Stable) {
            resolved_regime = raw_regime;
        }

        // Nearby collision belts cannot change an oceanic cell into continental crust.
        if (resolved_regime == tectonics::contract::GeologicRegime::ContinentCollision &&
            crust_class_map[index] == static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)) {
            resolved_regime = tectonics::contract::GeologicRegime::ConvergentArc;
        }
        geologic_regime_map[index] = static_cast<uint8_t>(resolved_regime);
    }
}

void lithosphere::updateTectonicJunctions() {
    const tectonics::tectonic_junctions::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        imap.raw_data(),
        boundary_segment_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        plate_kinematics.data(),
        static_cast<uint32_t>(plate_kinematics.size()),
    };
    tectonics::tectonic_junctions::Outputs outputs = tectonics::tectonic_junctions::build(inputs);
    tectonic_junctions = std::move(outputs.junctions);
}

void lithosphere::updateTectonicJunctionCouplingMaps() {
    const uint32_t map_area = _worldDimension.getArea();
    const tectonics::tectonic_junction_coupling::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        tectonic_junctions.data(),
        static_cast<uint32_t>(tectonic_junctions.size()),
    };
    tectonics::tectonic_junction_coupling::Outputs outputs =
        tectonics::tectonic_junction_coupling::build(inputs);
    junction_rift_influence_map = std::move(outputs.rift_influence);
    junction_transform_influence_map = std::move(outputs.transform_influence);
    junction_subduction_influence_map = std::move(outputs.subduction_influence);

    if (outputs.rift_tension_x.size() != map_area || outputs.rift_tension_y.size() != map_area ||
        rift_tension_x_map.size() != map_area || rift_tension_y_map.size() != map_area) {
        return;
    }

    constexpr float kJunctionRiftTensionGain = 18.0f * 0.35f;
    constexpr float kMaxTension = 100.0f * 2.5f;
    for (uint32_t index = 0; index < map_area; ++index) {
        rift_tension_x_map[index] += outputs.rift_tension_x[index] * kJunctionRiftTensionGain;
        rift_tension_y_map[index] += outputs.rift_tension_y[index] * kJunctionRiftTensionGain;
        const float magnitude =
            std::sqrt(rift_tension_x_map[index] * rift_tension_x_map[index] +
                      rift_tension_y_map[index] * rift_tension_y_map[index]);
        if (magnitude > kMaxTension) {
            const float scale = kMaxTension / magnitude;
            rift_tension_x_map[index] *= scale;
            rift_tension_y_map[index] *= scale;
        }
    }
    updateTensionMagnitudeMap();
}

void lithosphere::updateDeformingRegionMaps() {
    const tectonics::deforming_regions::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        imap.raw_data(),
        crust_class_map.data(),
        convergence_map.data(),
        divergence_map.data(),
        shear_map.data(),
        boundary_type_map.data(),
        geologic_regime_map.data(),
        boundary_distance_map.data(),
        nearest_boundary_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        plate_kinematics.data(),
        static_cast<uint32_t>(plate_kinematics.size()),
    };
    tectonics::deforming_regions::Outputs outputs =
        tectonics::deforming_regions::build(inputs);
    deforming_regions = std::move(outputs.regions);
    deforming_region_id_map = std::move(outputs.deforming_region_id);
    deforming_region_type_map = std::move(outputs.deforming_region_type);
    deformation_rate_map = std::move(outputs.deformation_rate);
    deformation_velocity_x_map = std::move(outputs.deformation_velocity_x);
    deformation_velocity_y_map = std::move(outputs.deformation_velocity_y);
}

void lithosphere::updateFieldStackMaps() {
    if (num_plates == 0) {
        return;
    }

    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    std::vector<float> continental_thickness(map_area, 0.0f);
    std::vector<float> oceanic_thickness(map_area, 0.0f);

    std::fill(crust_thickness_map.begin(), crust_thickness_map.end(), 0.0f);

    for (uint32_t plate_index = 0; plate_index < num_plates; ++plate_index) {
        const plate* current_plate = plates[plate_index];
        const uint32_t x0 = current_plate->getLeftAsUint();
        const uint32_t y0 = current_plate->getTopAsUint();
        const uint32_t x1 = x0 + current_plate->getWidth();
        const uint32_t y1 = y0 + current_plate->getHeight();

        const float* plate_map = nullptr;
        const uint32_t* plate_age = nullptr;
        const uint8_t* plate_material = nullptr;
        current_plate->getMap(&plate_map, &plate_age, &plate_material);
        (void)plate_age;

        const uint32_t x_mod_start = x0 % width;
        for (uint32_t y = y0, j = 0; y < y1; ++y) {
            const uint32_t y_mod = y % height;
            const uint32_t row_offset = y_mod * width;
            uint32_t x_mod = x_mod_start;
            for (uint32_t x = x0; x < x1;
                 ++x, ++j, x_mod = ++x_mod >= width ? x_mod - width : x_mod) {
                const float crust = plate_map[j];
                if (crust <= 2.0f * FLT_EPSILON) {
                    continue;
                }

                const uint32_t index = row_offset + x_mod;
                crust_thickness_map[index] += crust;
                if (is_continental_material(plate_material[j])) {
                    continental_thickness[index] += crust;
                } else {
                    oceanic_thickness[index] += crust;
                }
            }
        }
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        crust_thickness_map[index] =
            std::max(crust_thickness_map[index], std::max(0.0f, display_hmap[index]));
        crust_class_map[index] = static_cast<uint8_t>(classify_crust_class(
            crust_thickness_map[index], continental_thickness[index]));
        if (imap[index] >= num_plates || hmap[index] <= 2.0f * FLT_EPSILON) {
            surface_crust_type_map[index] =
                static_cast<uint8_t>(tectonics::contract::CrustClass::None);
        } else {
            surface_crust_type_map[index] =
                static_cast<uint8_t>(crust_class_from_material(material_map[index]));
        }
    }

    updateTectonicProvenance();
    updateBoundaryGraphMaps();
    updateTectonicJunctions();
    updateDeformingRegionMaps();
    updateTectonicJunctionCouplingMaps();

    constexpr float kSurfaceDeltaScale = 0.08f;
    constexpr float kYoungCrustMyr = 35.0f;
    constexpr float kOldCrustMyr = 80.0f;

    for (uint32_t index = 0; index < map_area; ++index) {
        const float total_thickness = crust_thickness_map[index];
        const float continental_ratio =
            total_thickness > FLT_EPSILON ? continental_thickness[index] / total_thickness : 0.0f;
        const float oceanic_ratio =
            total_thickness > FLT_EPSILON ? oceanic_thickness[index] / total_thickness : 0.0f;

        const float convergence = static_cast<float>(convergence_map[index]) / 100.0f;
        const float divergence = static_cast<float>(divergence_map[index]) / 100.0f;
        const float shear = static_cast<float>(shear_map[index]) / 100.0f;
        const auto regime =
            static_cast<tectonics::contract::GeologicRegime>(geologic_regime_map[index]);
        const auto boundary_type =
            static_cast<tectonics::contract::BoundaryType>(boundary_type_map[index]);
        const auto deforming_region_type = static_cast<tectonics::contract::DeformingRegionType>(
            deforming_region_type_map[index]);
        const float deformation_rate = deformation_rate_map[index];
        const float deformation_velocity = std::sqrt(
            deformation_velocity_x_map[index] * deformation_velocity_x_map[index] +
            deformation_velocity_y_map[index] * deformation_velocity_y_map[index]);
        const float junction_rift = junction_rift_influence_map[index];
        const float junction_transform = junction_transform_influence_map[index];
        const float junction_subduction = junction_subduction_influence_map[index];

        const float delta_surface = display_hmap[index] - prev_display_hmap[index];
        const float positive_surface = clamp_unit(delta_surface / kSurfaceDeltaScale);
        const float negative_surface = clamp_unit(-delta_surface / kSurfaceDeltaScale);
        const float age_old = clamp_unit(crust_age_myr_map[index] / kOldCrustMyr);
        const float age_young = 1.0f - clamp_unit(crust_age_myr_map[index] / kYoungCrustMyr);
        if (deforming_region_type ==
            tectonics::contract::DeformingRegionType::ContinentalRift) {
            crust_thickness_map[index] =
                std::max(display_hmap[index],
                         crust_thickness_map[index] * (1.0f - 0.12f * deformation_rate));
        } else if (deforming_region_type ==
                   tectonics::contract::DeformingRegionType::DiffuseCollision) {
            crust_thickness_map[index] +=
                deformation_rate * (0.12f + 0.18f * continental_ratio);
        }
        if (junction_rift > 0.0f &&
            (static_cast<tectonics::contract::CrustClass>(crust_class_map[index]) ==
                 tectonics::contract::CrustClass::Continental ||
             static_cast<tectonics::contract::CrustClass>(crust_class_map[index]) ==
                 tectonics::contract::CrustClass::Transitional)) {
            crust_thickness_map[index] =
                std::max(display_hmap[index], crust_thickness_map[index] * (1.0f - 0.045f * junction_rift));
        }

        float uplift =
            0.50f * convergence + 0.10f * shear + 0.15f * positive_surface +
            0.10f * continental_ratio + 0.08f * age_young;
        float subsidence =
            0.24f * divergence + 0.18f * negative_surface +
            0.18f * (oceanic_ratio * age_old) + 0.10f * (1.0f - continental_ratio) +
            0.06f * (boundary_type == tectonics::contract::BoundaryType::PassiveMargin ? 1.0f : 0.0f);

        switch (regime) {
        case tectonics::contract::GeologicRegime::ConvergentArc:
            uplift += 0.24f;
            break;
        case tectonics::contract::GeologicRegime::ContinentCollision:
            uplift += 0.30f;
            break;
        case tectonics::contract::GeologicRegime::TrenchAdjacent:
            uplift += 0.10f;
            subsidence += 0.16f;
            break;
        case tectonics::contract::GeologicRegime::DivergentRift:
            subsidence += 0.18f;
            break;
        case tectonics::contract::GeologicRegime::PassiveMargin:
            subsidence += 0.20f;
            break;
        case tectonics::contract::GeologicRegime::MidOceanRidge:
            uplift += 0.12f;
            subsidence -= 0.08f;
            break;
        case tectonics::contract::GeologicRegime::Transform:
            uplift += 0.05f;
            break;
        case tectonics::contract::GeologicRegime::Stable:
        default:
            break;
        }

        switch (deforming_region_type) {
        case tectonics::contract::DeformingRegionType::ContinentalRift:
            subsidence += 0.28f * deformation_rate;
            subsidence += 0.05f * deformation_velocity;
            break;
        case tectonics::contract::DeformingRegionType::DiffuseCollision:
            uplift += 0.30f * deformation_rate;
            uplift += 0.04f * deformation_velocity;
            subsidence -= 0.05f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::None:
        default:
            break;
        }
        subsidence += 0.20f * junction_rift;
        uplift += 0.05f * junction_rift;
        uplift += 0.04f * junction_transform;
        subsidence += 0.02f * junction_transform;
        uplift += 0.09f * junction_subduction;
        subsidence += 0.12f * junction_subduction;

        uplift_tendency_map[index] = clamp_unit(uplift);
        subsidence_tendency_map[index] = clamp_unit(subsidence);

        float strain_input =
            0.26f * convergence + 0.22f * divergence + 0.34f * shear +
            0.10f * std::max(positive_surface, negative_surface);
        switch (boundary_type) {
        case tectonics::contract::BoundaryType::Convergent:
            strain_input += 0.10f;
            break;
        case tectonics::contract::BoundaryType::Divergent:
            strain_input += 0.08f;
            break;
        case tectonics::contract::BoundaryType::Transform:
            strain_input += 0.12f;
            break;
        case tectonics::contract::BoundaryType::PassiveMargin:
            strain_input += 0.03f;
            break;
        case tectonics::contract::BoundaryType::None:
        default:
            break;
        }
        switch (deforming_region_type) {
        case tectonics::contract::DeformingRegionType::ContinentalRift:
            strain_input += 0.14f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::DiffuseCollision:
            strain_input += 0.18f * deformation_rate;
            break;
        case tectonics::contract::DeformingRegionType::None:
        default:
            break;
        }
        accumulated_strain_map[index] =
            clamp_unit(accumulated_strain_map[index] * 0.965f + clamp_unit(strain_input) * 0.22f);
    }
}

void lithosphere::updateDerivedMaps() {
    if (num_plates == 0) updateTectonicProvenance();
    updateFieldStackMaps();
}

void lithosphere::initializeHeightMapFromFloat(const float* heightmap) {
    if (heightmap == nullptr) {
        throw invalid_argument("Height map cannot be null");
    }

    const uint32_t map_area = _worldDimension.getArea();
    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = heightmap[i];
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
}

void lithosphere::initializeHeightMapFromMetric(const uint16_t* heightmap_m, uint16_t new_sea_level_m) {
    if (heightmap_m == nullptr) {
        throw invalid_argument("Metric height map cannot be null");
    }

    sea_level_m = new_sea_level_m;
    const uint32_t map_area = _worldDimension.getArea();
    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = TopographyCodec::meters_to_internal(heightmap_m[i], sea_level_m);
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
}

void lithosphere::seedInitialTopography(float ocean_coverage, int32_t sea_level_m_override) {
    const uint32_t map_area = _worldDimension.getArea();
    const bool has_override = sea_level_m_override != TopographyCodec::kNoSeaLevelOverride;
    if (sea_level_m_override < TopographyCodec::kNoSeaLevelOverride ||
        sea_level_m_override > TopographyCodec::kMaxHeightMeters) {
        throw invalid_argument("sea_level_m is out of range");
    }
    sea_level_m = has_override ? static_cast<uint16_t>(sea_level_m_override)
                               : TopographyCodec::default_import_sea_level_m();
    if (initial_ocean_depth_m == 0U || initial_ocean_depth_m > sea_level_m ||
        initial_land_relief_m > TopographyCodec::kMaxHeightMeters - sea_level_m) {
        throw invalid_argument("Initial depth/relief must fit around the sea-level datum");
    }
    SimpleRandom topo_rng(_randsource);

    const FractalNoiseConfig warp_x_noise = make_noise_config(topo_rng, 4.0f, 0.58f, 0.45f, 1.6f);
    const FractalNoiseConfig warp_y_noise = make_noise_config(topo_rng, 4.0f, 0.58f, 0.45f, 1.6f);
    const FractalNoiseConfig continent_macro_noise =
        make_noise_config(topo_rng, 5.0f, 0.57f, 0.16f, 0.95f);
    const FractalNoiseConfig continent_detail_noise =
        make_noise_config(topo_rng, 6.0f, 0.59f, 0.42f, 1.85f);
    const FractalNoiseConfig coast_breakup_noise =
        make_noise_config(topo_rng, 6.0f, 0.57f, 0.62f, 2.2f);
    const FractalNoiseConfig shelf_noise = make_noise_config(topo_rng, 6.0f, 0.61f, 0.85f, 2.9f);
    const FractalNoiseConfig sea_base_noise = make_noise_config(topo_rng, 7.0f, 0.63f, 1.15f, 3.9f);
    const FractalNoiseConfig sea_detail_noise =
        make_noise_config(topo_rng, 7.0f, 0.61f, 1.95f, 5.6f);
    const FractalNoiseConfig sea_micro_noise = make_noise_config(topo_rng, 5.0f, 0.48f, 3.95f, 10.2f);
    const FractalNoiseConfig land_noise = make_noise_config(topo_rng, 6.0f, 0.56f, 0.95f, 3.6f);
    const FractalNoiseConfig foothill_noise =
        make_noise_config(topo_rng, 6.0f, 0.53f, 1.55f, 4.7f);
    const FractalNoiseConfig mountain_noise =
        make_noise_config(topo_rng, 5.0f, 0.47f, 2.45f, 6.3f);
    const FractalNoiseConfig ridge_noise = make_noise_config(topo_rng, 5.0f, 0.50f, 2.6f, 7.4f);
    const FractalNoiseConfig trench_noise = make_noise_config(topo_rng, 5.0f, 0.48f, 3.8f, 9.8f);

    std::vector<float> warped_x(map_area);
    std::vector<float> warped_y(map_area);
    std::vector<float> continentality(map_area);
    std::vector<float> raw_height(map_area);
    std::vector<uint8_t> is_land(map_area, 0);

    for (uint32_t y = 0; y < _worldDimension.getHeight(); ++y) {
        for (uint32_t x = 0; x < _worldDimension.getWidth(); ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const float warp_x = sample_toroidal_noise(static_cast<float>(x), static_cast<float>(y),
                                                       _worldDimension, warp_x_noise) *
                                 18.0f;
            const float warp_y = sample_toroidal_noise(static_cast<float>(x), static_cast<float>(y),
                                                       _worldDimension, warp_y_noise) *
                                 18.0f;
            const float sample_x = static_cast<float>(x) + warp_x;
            const float sample_y = static_cast<float>(y) + warp_y;
            warped_x[index] = sample_x;
            warped_y[index] = sample_y;
            const float macro = sample_toroidal_noise(sample_x, sample_y, _worldDimension,
                                                      continent_macro_noise);
            const float detail = sample_toroidal_noise(sample_x * 1.15f, sample_y * 1.15f,
                                                       _worldDimension, continent_detail_noise);
            const float coast_breakup = sample_toroidal_noise(sample_x * 1.55f, sample_y * 1.55f,
                                                              _worldDimension,
                                                              coast_breakup_noise);
            continentality[index] = macro * 0.70f + detail * 0.22f + coast_breakup * 0.08f;
        }
    }

    const float continent_threshold = quantile_threshold(continentality, ocean_coverage);
    float raw_land_min = FLT_MAX;
    float raw_land_max = -FLT_MAX;
    float raw_ocean_min = FLT_MAX;
    float raw_ocean_max = -FLT_MAX;

    for (uint32_t y = 0; y < _worldDimension.getHeight(); ++y) {
        for (uint32_t x = 0; x < _worldDimension.getWidth(); ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const float sample_x = warped_x[index];
            const float sample_y = warped_y[index];
            const float detail = sample_toroidal_noise(sample_x * 1.15f, sample_y * 1.15f,
                                                       _worldDimension, continent_detail_noise);
            const float shelf = sample_toroidal_noise(sample_x * 1.4f, sample_y * 1.4f,
                                                      _worldDimension, shelf_noise);
            const float sea_base = sample_toroidal_noise(sample_x * 1.15f, sample_y * 1.15f,
                                                         _worldDimension, sea_base_noise);
            const float sea_detail = sample_toroidal_noise(sample_x * 2.3f, sample_y * 2.3f,
                                                           _worldDimension, sea_detail_noise);
            const float sea_micro = sample_toroidal_noise(sample_x * 4.8f, sample_y * 4.8f,
                                                          _worldDimension, sea_micro_noise);
            const float land_detail = sample_toroidal_noise(sample_x * 1.55f, sample_y * 1.55f,
                                                            _worldDimension, land_noise);
            const float foothills = sample_toroidal_noise(sample_x * 2.15f, sample_y * 2.15f,
                                                          _worldDimension, foothill_noise);
            const float mountains =
                1.0f - std::abs(sample_toroidal_noise(sample_x * 3.05f, sample_y * 3.05f,
                                                      _worldDimension, mountain_noise));
            const float ridge = 1.0f - std::abs(sample_toroidal_noise(sample_x * 2.1f, sample_y * 2.1f,
                                                                      _worldDimension, ridge_noise));
            const float trench = 1.0f - std::abs(sample_toroidal_noise(sample_x * 2.35f, sample_y * 2.35f,
                                                                       _worldDimension, trench_noise));
            const float distance_to_coast = continentality[index] - continent_threshold;
            const float coastal_transition =
                smoothstep(0.010f, 0.110f, std::abs(distance_to_coast));
            const float shoreline_profile =
                0.09f * std::tanh(distance_to_coast * 18.0f) + 0.07f * detail + 0.05f * shelf +
                0.04f * land_detail + 0.06f * ridge - 0.06f * trench + 0.03f * foothills;

            if (distance_to_coast >= 0.0f) {
                is_land[index] = 1U;
                const float inland_bias = logistic(distance_to_coast + 0.05f * land_detail, 8.0f);
                const float coastal_bias = 1.0f - inland_bias;
                const float inland_profile =
                    0.34f * inland_bias - 0.08f * coastal_bias + 0.18f * detail +
                    0.19f * land_detail + 0.15f * foothills + 0.17f * mountains +
                    0.10f * ridge - 0.05f * trench + 0.04f * shelf;
                raw_height[index] =
                    interpolate(shoreline_profile, inland_profile, coastal_transition);
                raw_land_min = std::min(raw_land_min, raw_height[index]);
                raw_land_max = std::max(raw_land_max, raw_height[index]);
            } else {
                const float shelf_weight = logistic(distance_to_coast + 0.10f * shelf + 0.05f, 11.0f);
                const float abyss_weight = 1.0f - shelf_weight;
                const float ocean_profile =
                    -0.62f * abyss_weight - 0.16f * shelf_weight + 0.17f * sea_base +
                    0.15f * sea_detail + 0.09f * sea_micro + 0.08f * shelf + 0.17f * ridge -
                    0.27f * trench + 0.04f * detail;
                raw_height[index] =
                    interpolate(shoreline_profile, ocean_profile, coastal_transition);
                raw_ocean_min = std::min(raw_ocean_min, raw_height[index]);
                raw_ocean_max = std::max(raw_ocean_max, raw_height[index]);
            }
        }
    }

    // Scale physical offsets, then encode once. Both profiles meet the sea datum
    // at the continentality threshold; independent rescaling must not create a
    // multi-kilometre jump between the nearest ocean and land samples.
    std::vector<float> profile(map_area);
    float land_peak = 0.0f;
    float ocean_peak = 0.0f;
    for (uint32_t i = 0; i < map_area; ++i) {
        const float coast_weight = smoothstep(
            0.0f, 0.110f, std::abs(continentality[i] - continent_threshold));
        const float amplitude = is_land[i]
            ? remap_to_range(raw_height[i], raw_land_min, raw_land_max, 0.0f, 1.0f)
            : remap_to_range(raw_height[i], raw_ocean_min, raw_ocean_max, 1.0f, 0.0f);
        profile[i] = coast_weight * amplitude;
        float& peak = is_land[i] ? land_peak : ocean_peak;
        peak = std::max(peak, profile[i]);
    }
    std::vector<uint16_t> seeded_heightmap(map_area);
    for (uint32_t i = 0; i < map_area; ++i) {
        const float peak = is_land[i] ? land_peak : ocean_peak;
        const float fraction = peak > FLT_EPSILON ? profile[i] / peak : 0.0f;
        const int offset_m = is_land[i]
            ? static_cast<int>(std::lround(initial_land_relief_m * fraction))
            : -std::max(1, static_cast<int>(std::lround(initial_ocean_depth_m * fraction)));
        seeded_heightmap[i] = static_cast<uint16_t>(sea_level_m + offset_m);
    }
    initializeHeightMapFromMetric(seeded_heightmap.data(), sea_level_m);
}

void lithosphere::fillEnclosedInitialOceans() {
    using tectonics::contract::CrustClass;

    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U) {
        return;
    }

    std::vector<uint8_t> visited(map_area, 0U);
    std::vector<uint32_t> largest_component;
    largest_component.reserve(map_area / 2U);

    for (uint32_t index = 0; index < map_area; ++index) {
        if (visited[index] != 0U ||
            surface_crust_type_map[index] != static_cast<uint8_t>(CrustClass::Oceanic)) {
            continue;
        }

        std::vector<uint32_t> component;
        component.reserve(256U);
        std::deque<uint32_t> frontier;
        frontier.push_back(index);
        visited[index] = 1U;

        while (!frontier.empty()) {
            const uint32_t current = frontier.front();
            frontier.pop_front();
            component.push_back(current);

            const uint32_t x = _worldDimension.xFromIndex(current);
            const uint32_t y = _worldDimension.yFromIndex(current);
            const uint32_t left_x = x > 0U ? x - 1U : width - 1U;
            const uint32_t right_x = x + 1U < width ? x + 1U : 0U;
            const uint32_t top_y = y > 0U ? y - 1U : height - 1U;
            const uint32_t bottom_y = y + 1U < height ? y + 1U : 0U;
            const uint32_t neighbors[] = {
                _worldDimension.indexOf(left_x, y),
                _worldDimension.indexOf(right_x, y),
                _worldDimension.indexOf(x, top_y),
                _worldDimension.indexOf(x, bottom_y),
            };

            for (uint32_t neighbor : neighbors) {
                if (visited[neighbor] != 0U ||
                    surface_crust_type_map[neighbor] !=
                        static_cast<uint8_t>(CrustClass::Oceanic)) {
                    continue;
                }
                visited[neighbor] = 1U;
                frontier.push_back(neighbor);
            }
        }

        if (component.size() > largest_component.size()) {
            largest_component.swap(component);
        }
    }

    std::vector<uint8_t> keep_ocean(map_area, 0U);
    for (uint32_t index : largest_component) {
        keep_ocean[index] = 1U;
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (surface_crust_type_map[index] == static_cast<uint8_t>(CrustClass::Oceanic) &&
            keep_ocean[index] == 0U) {
            surface_crust_type_map[index] = static_cast<uint8_t>(CrustClass::Continental);
        }
    }
}

void lithosphere::seedInitialCrustTypes() {
    using tectonics::contract::CrustClass;

    const uint32_t map_area = _worldDimension.getArea();
    const uint16_t threshold_m = resolve_crust_type_boundary_m(scenario, sea_level_m);
    for (uint32_t index = 0; index < map_area; ++index) {
        const uint16_t height_m = TopographyCodec::internal_to_meters(hmap[index], sea_level_m);
        surface_crust_type_map[index] =
            height_m >= threshold_m ? static_cast<uint8_t>(CrustClass::Continental)
                                    : static_cast<uint8_t>(CrustClass::Oceanic);
    }

    fillEnclosedInitialOceans();
}

void lithosphere::seedMaterialMap() {
    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U) {
        material_map.clear();
        return;
    }

    SimpleRandom material_rng(static_cast<uint32_t>(scenario.seed) ^ 0x6A09E667U);
    const FractalNoiseConfig warp_x_noise =
        make_noise_config(material_rng, 4.0f, 0.58f, 0.28f, 1.4f);
    const FractalNoiseConfig warp_y_noise =
        make_noise_config(material_rng, 4.0f, 0.58f, 0.28f, 1.4f);
    const FractalNoiseConfig material_macro_noise =
        make_noise_config(material_rng, 5.0f, 0.57f, 0.22f, 1.2f);
    const FractalNoiseConfig material_detail_noise =
        make_noise_config(material_rng, 6.0f, 0.61f, 0.75f, 2.4f);

    for (uint32_t y = 0; y < _worldDimension.getHeight(); ++y) {
        for (uint32_t x = 0; x < _worldDimension.getWidth(); ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const float warp_x =
                static_cast<float>(x) +
                sample_toroidal_noise(static_cast<float>(x), static_cast<float>(y),
                                      _worldDimension, warp_x_noise) *
                    2.4f;
            const float warp_y =
                static_cast<float>(y) +
                sample_toroidal_noise(static_cast<float>(x), static_cast<float>(y),
                                      _worldDimension, warp_y_noise) *
                    2.4f;
            const float macro =
                sample_toroidal_noise(warp_x, warp_y, _worldDimension, material_macro_noise);
            const float detail =
                sample_toroidal_noise(warp_x * 1.8f, warp_y * 1.8f, _worldDimension,
                                      material_detail_noise);
            const float composite = 0.68f * macro + 0.32f * detail;

            uint8_t material_index =
                tectonics::material::to_index(tectonics::material::Type::Granite);
            if (surface_crust_type_map[index] ==
                static_cast<uint8_t>(tectonics::contract::CrustClass::Oceanic)) {
                material_index = tectonics::material::to_index(tectonics::material::Type::Basalt);
            } else if (composite < -0.08f) {
                material_index =
                    tectonics::material::to_index(tectonics::material::Type::Sedimentary);
            } else if (composite > 0.48f) {
                material_index =
                    tectonics::material::to_index(tectonics::material::Type::Metamorphic);
            }

            material_map[index] = material_index;
        }
    }
}

void lithosphere::seedInitialCrustAges() {
    using tectonics::contract::BoundaryType;
    using tectonics::contract::CrustClass;
    using tectonics::contract::GeologicRegime;

    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U || iter_count == 0U) {
        return;
    }

    const double step_myr = std::max(0.001, scenario.delta_time_myr);
    std::vector<int32_t> ocean_distance(map_area, -1);
    std::deque<uint32_t> frontier;

    const auto is_oceanic_cell = [&](uint32_t index) {
        return index < map_area && imap[index] < num_plates && hmap[index] > FLT_EPSILON &&
               crust_class_map[index] == static_cast<uint8_t>(CrustClass::Oceanic);
    };
    const auto push_seed = [&](uint32_t index) {
        if (!is_oceanic_cell(index) || ocean_distance[index] >= 0) {
            return;
        }
        ocean_distance[index] = 0;
        frontier.push_back(index);
    };

    for (uint32_t index = 0; index < map_area; ++index) {
        if (!is_oceanic_cell(index)) {
            continue;
        }

        const auto boundary_type = static_cast<BoundaryType>(boundary_type_map[index]);
        const auto regime = static_cast<GeologicRegime>(geologic_regime_map[index]);
        if (boundary_type == BoundaryType::Divergent || regime == GeologicRegime::DivergentRift ||
            regime == GeologicRegime::MidOceanRidge) {
            push_seed(index);
        }
    }

    if (frontier.empty()) {
        for (uint32_t index = 0; index < map_area; ++index) {
            if (!is_oceanic_cell(index)) {
                continue;
            }
            if (static_cast<BoundaryType>(boundary_type_map[index]) != BoundaryType::None) {
                push_seed(index);
            }
        }
    }

    if (frontier.empty()) {
        for (uint32_t index = 0; index < map_area; ++index) {
            push_seed(index);
        }
    }

    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    while (!frontier.empty()) {
        const uint32_t current = frontier.front();
        frontier.pop_front();
        const int32_t next_distance = ocean_distance[current] + 1;
        const uint32_t x = _worldDimension.xFromIndex(current);
        const uint32_t y = _worldDimension.yFromIndex(current);
        const uint32_t neighbors[] = {
            _worldDimension.indexOf(x > 0U ? x - 1U : world_width - 1U, y),
            _worldDimension.indexOf(x + 1U < world_width ? x + 1U : 0U, y),
            _worldDimension.indexOf(x, y > 0U ? y - 1U : world_height - 1U),
            _worldDimension.indexOf(x, y + 1U < world_height ? y + 1U : 0U),
        };

        for (uint32_t neighbor : neighbors) {
            if (!is_oceanic_cell(neighbor) || ocean_distance[neighbor] >= 0) {
                continue;
            }
            ocean_distance[neighbor] = next_distance;
            frontier.push_back(neighbor);
        }
    }

    const FractalNoiseConfig continental_macro_noise =
        make_noise_config(static_cast<uint32_t>(scenario.seed) ^ 0x3C6EF372U, 5.0f, 0.57f, 0.26f, 1.3f);
    const FractalNoiseConfig continental_detail_noise =
        make_noise_config(static_cast<uint32_t>(scenario.seed) ^ 0xBB67AE85U, 6.0f, 0.60f, 0.72f, 2.5f);
    const FractalNoiseConfig ocean_age_noise =
        make_noise_config(static_cast<uint32_t>(scenario.seed) ^ 0xA54FF53AU, 4.0f, 0.58f, 0.34f, 1.6f);
    const FractalNoiseConfig ocean_detail_noise =
        make_noise_config(static_cast<uint32_t>(scenario.seed) ^ 0x510E527FU, 5.0f, 0.61f, 0.92f, 3.1f);

    uint32_t max_age_steps = 0U;
    std::vector<uint32_t> timestamps(map_area, 0U);

    for (uint32_t y = 0; y < world_height; ++y) {
        for (uint32_t x = 0; x < world_width; ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            if (imap[index] >= num_plates || hmap[index] <= FLT_EPSILON) {
                continue;
            }

            const auto crust_class = static_cast<CrustClass>(crust_class_map[index]);
            float age_myr = 0.0f;
            if (crust_class == CrustClass::Oceanic) {
                const int32_t distance = std::max(0, ocean_distance[index]);
                const float macro =
                    sample_toroidal_noise(static_cast<float>(x), static_cast<float>(y),
                                          _worldDimension, ocean_age_noise);
                const float detail =
                    sample_toroidal_noise(static_cast<float>(x) * 1.8f,
                                          static_cast<float>(y) * 1.8f, _worldDimension,
                                          ocean_detail_noise);
                age_myr = static_cast<float>(distance) * 0.55f;
                if (distance > 0) {
                    age_myr += 2.0f * macro + 1.5f * detail;
                }
                age_myr = std::max(0.0f, std::min(140.0f, age_myr));
            } else {
                const float macro =
                    0.5f + 0.5f * sample_toroidal_noise(static_cast<float>(x) * 0.75f,
                                                        static_cast<float>(y) * 0.75f,
                                                        _worldDimension, continental_macro_noise);
                const float detail =
                    0.5f + 0.5f * sample_toroidal_noise(static_cast<float>(x) * 1.7f,
                                                        static_cast<float>(y) * 1.7f,
                                                        _worldDimension, continental_detail_noise);
                const float normalized = clamp_unit(0.68f * macro + 0.32f * detail);
                const float min_age = crust_class == CrustClass::Transitional ? 55.0f : 110.0f;
                const float span = crust_class == CrustClass::Transitional ? 55.0f : 120.0f;
                age_myr = min_age + span * normalized;
            }

            const uint32_t age_steps =
                static_cast<uint32_t>(std::llround(age_myr / step_myr));
            max_age_steps = std::max(max_age_steps, age_steps);
            timestamps[index] = age_steps;
        }
    }

    if (max_age_steps + MAX_BUOYANCY_AGE > iter_count) {
        iter_count = max_age_steps + MAX_BUOYANCY_AGE;
        time_origin_step = iter_count;
    }

    for (uint32_t y = 0; y < world_height; ++y) {
        for (uint32_t x = 0; x < world_width; ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const uint32_t owner = imap[index];
            if (owner >= num_plates || hmap[index] <= FLT_EPSILON) {
                amap[index] = 0U;
                continue;
            }

            const uint32_t timestamp =
                timestamps[index] <= iter_count ? iter_count - timestamps[index] : 0U;
            amap[index] = timestamp;
            plates[owner]->setCrustTimestamp(x, y, timestamp);
        }
    }
}

float lithosphere::internalToPhysicalMeters(float internal_height) const noexcept {
    return internal_to_physical_meters(internal_height, sea_level_m);
}

float lithosphere::metersToInternalDelta(float delta_m) const noexcept {
    if (delta_m <= 0.0f) {
        return 0.0f;
    }

    const float land_range =
        static_cast<float>(std::max<uint16_t>(1U, TopographyCodec::kMaxHeightMeters - sea_level_m));
    return delta_m / land_range;
}

void lithosphere::importNormalizedHeightMap(const float* normalized_map, float sea_level) {
    if (normalized_map == nullptr) {
        throw invalid_argument("Normalized height map cannot be null");
    }

    const uint32_t map_area = _worldDimension.getArea();
    const float sea_threshold =
        TopographyCodec::infer_normalized_sea_threshold(normalized_map, map_area, sea_level);
    sea_level_m = static_cast<uint16_t>(std::lround(
        sea_threshold * static_cast<float>(TopographyCodec::kMaxHeightMeters)));

    for (uint32_t i = 0; i < map_area; ++i) {
        hmap[i] = TopographyCodec::normalized_to_internal(normalized_map[i], sea_threshold);
    }
    prev_hmap = hmap;
    display_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = hmap;
    seedInitialCrustTypes();
    seedMaterialMap();

    resetSimulationState();
    createPlates();
}

void lithosphere::importRawHeightMap(const float* normalized_map) {
    if (normalized_map == nullptr) {
        throw invalid_argument("Normalized height map cannot be null");
    }

    const uint32_t map_area = _worldDimension.getArea();
    std::vector<uint16_t> metric_heightmap(map_area);
    for (uint32_t i = 0; i < map_area; ++i) {
        metric_heightmap[i] = static_cast<uint16_t>(std::lround(
            TopographyCodec::clamp_normalized(normalized_map[i]) *
            static_cast<float>(TopographyCodec::kMaxHeightMeters)));
    }

    initializeHeightMapFromMetric(metric_heightmap.data(), TopographyCodec::default_import_sea_level_m());
    seedInitialCrustTypes();
    seedMaterialMap();
    resetSimulationState();
    createPlates();
}

void lithosphere::importMetricHeightMap(const uint16_t* heightmap_m, uint16_t new_sea_level_m) {
    initializeHeightMapFromMetric(heightmap_m, new_sea_level_m);
    seedInitialCrustTypes();
    seedMaterialMap();
    resetSimulationState();
    createPlates();
}

bool lithosphere::isFinished() const {
    return finished || getPlateCount() == 0;
}

// At least two plates are at same location.
// Move some crust from the SMALLER plate onto LARGER one.
// Update height and plate index maps.
// Doing it plate by plate is much faster than doing it index wise:
// Each plate's map's memory area is accessed sequentially and only
// once as opposed to calculating "num_plates" indices within plate
// maps in order to find out which plate(s) own current location.
void lithosphere::updateTensionMagnitudeMap() noexcept {
    const uint32_t map_area = _worldDimension.getArea();
    if (rift_tension_map.size() != map_area || rift_tension_x_map.size() != map_area ||
        rift_tension_y_map.size() != map_area) {
        return;
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        const float tx = rift_tension_x_map[index];
        const float ty = rift_tension_y_map[index];
        rift_tension_map[index] = std::sqrt(tx * tx + ty * ty);
    }
}

// Remove empty plates from the system.
void lithosphere::removeEmptyPlates() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        if (num_plates == 1)
            puts("ONLY ONE PLATE LEFT!");
        else if (plate_indices_found[i] == 0) {
            topology_model.removePlate(i);
            delete plates[i];
            plates[i] = plates[num_plates - 1];
            plate_indices_found[i] = plate_indices_found[num_plates - 1];

            // Life is seldom as simple as seems at first.
            // Replace the moved plate's index in the index map
            // to match its current position in the array!
            for (uint32_t j = 0; j < _worldDimension.getArea(); ++j)
                if (imap[j] == num_plates - 1)
                    imap[j] = i;

            --num_plates;
            --i;
        }
    }
}

void lithosphere::solvePlateMotion() {
    using namespace tectonics::dynamics;
    const uint32_t width = _worldDimension.getWidth(), height = _worldDimension.getHeight();
    const double cell_m = scenario.world_width_km * 1000.0 / width;
    const double dt_s = scenario.delta_time_myr * kSecondsPerMyr;
    std::vector<tectonics::dynamics::Plate> bodies(num_plates);
    std::vector<Contact> contacts;
    std::vector<size_t> contact_faces;
    std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::array<double,3>> activity;
    for (size_t face=0;face<motion_history.size();++face) {
        const auto& sample=motion_history[face];
        if (sample.a==sample.b || sample.a==UINT32_MAX) continue;
        auto& sum=activity[{sample.a,sample.b,static_cast<uint32_t>(face%2)}];
        sum[0]+=sample.ridge;sum[1]+=sample.slab;sum[2]+=1;
    }
    std::vector<uint32_t> pressure_owners;
    std::vector<double> pressure;
    if (scenario.mountain_force_range_km>0) {
        pressure_owners.resize(static_cast<size_t>(width)*height);
        for (uint32_t i=0;i<width*height;++i) pressure_owners[i]=imap[i];
        pressure.reserve(pressure_owners.size());
        for (const auto& c:column_model.cells()) pressure.push_back(tectonics::columns::mountain_force_n_m(c,scenario.sediment_load_feedback));
    }
    for (uint32_t y = 0; y < height; ++y) for (uint32_t x = 0; x < width; ++x) {
        const uint32_t index = y * width + x, owner = imap[index];
        if (owner >= num_plates) continue;
        auto& body = bodies[owner];
        const auto flow = mantle_flow((x + .5) / width, (y + .5) / height,
            getTimeMyr(), scenario.seed, scenario.mantle_flow_speed_cm_yr);
        body.area_m2 += cell_m * cell_m;
        body.mantle_velocity_m_s.x += flow.x * cell_m * cell_m;
        body.mantle_velocity_m_s.y += flow.y * cell_m * cell_m;
        // Each grid face is visited once, including the periodic seams.
        const uint32_t neighbors[] = {y * width + (x + 1) % width, ((y + 1) % height) * width + x};
        for (size_t axis = 0; axis < 2; ++axis) {
            const uint32_t neighbor = neighbors[axis], other = imap[neighbor];
            if (other == owner || other >= num_plates) continue;
            const auto a = plates[owner]->linearVelocityVector();
            const auto b = plates[other]->linearVelocityVector();
            contacts.push_back({owner, other, {axis == 0 ? 1.0 : 0.0, axis == 1 ? 1.0 : 0.0},
                cell_m, is_continental_material(material_map[index]),
                is_continental_material(material_map[neighbor]),
                (iter_count - std::min(iter_count, amap[index])) * scenario.delta_time_myr,
                (iter_count - std::min(iter_count, amap[neighbor])) * scenario.delta_time_myr,
                (axis == 0 ? a.x() - b.x() : a.y() - b.y()) * cell_m / dt_s});
            contacts.back().activation_speed_m_s = scenario.boundary_force_ramp_cm_yr / 3.15576e9;
            contact_faces.push_back(2*index+axis);
            if (scenario.consistent_motion) {
                const auto& stable=topology_model.diagnostics().stable_ids;
                const auto found=subduction_polarity.find(std::minmax(stable[owner],stable[other]));
                if (found!=subduction_polarity.end())
                    contacts.back().descending_side=found->second==stable[owner]?0:1;
                else {
                    const auto& ca=column_model.cells()[index];
                    const auto& cb=column_model.cells()[neighbor];
                    const double aa=ca.ocean>0?ca.oceanic_age_volume/ca.ocean:0;
                    const double ab=cb.ocean>0?cb.oceanic_age_volume/cb.ocean:0;
                    contacts.back().descending_side=aa>=ab?0:1;
                    if (!contacts.back().continental_a && !contacts.back().continental_b)
                        subduction_polarity.emplace(std::minmax(stable[owner],stable[other]),
                            contacts.back().descending_side==0?stable[owner]:stable[other]);
                }
            }
            {
                const auto& cells = column_model.cells();
                if (scenario.thermal_force_feedback) {
                    const auto thermal_age=[](const tectonics::columns::Cell& c) {
                        return c.ocean>0?std::max(0.0,
                            (c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean):0.0;
                    };
                    contacts.back().age_a_myr=thermal_age(cells[index]);
                    contacts.back().age_b_myr=thermal_age(cells[neighbor]);
                    contacts.back().slab_pull_scale=scenario.thermal_subsidence_strength;
                }
                contacts.back().mountain_force_n_m = .5 * (tectonics::columns::mountain_force_n_m(cells[index],scenario.sediment_load_feedback)
                    + tectonics::columns::mountain_force_n_m(cells[neighbor],scenario.sediment_load_feedback));
                if (!pressure.empty())
                    contacts.back().mountain_force_n_m=.5*(inland_pressure(width,height,cell_m,pressure_owners,pressure,
                        index,static_cast<uint32_t>(axis),-1,scenario.mountain_force_range_km*1000)
                        +inland_pressure(width,height,cell_m,pressure_owners,pressure,
                        neighbor,static_cast<uint32_t>(axis),1,scenario.mountain_force_range_km*1000));
                // The same excess crustal GPE can resist compression at an
                // active continental margin. Previously the force solver
                // discarded this pressure whenever the other side was oceanic.
                if (!contacts.back().continental_a || !contacts.back().continental_b)
                    contacts.back().mountain_force_n_m *= (contacts.back().continental_a || contacts.back().continental_b)
                        ? scenario.coastal_gpe_strength : 0;
            }
        }
    }
    for (auto& body : bodies) {
        if (body.area_m2 <= 0.0) throw runtime_error("Force balance found an empty plate owner");
        body.mantle_velocity_m_s.x /= body.area_m2;
        body.mantle_velocity_m_s.y /= body.area_m2;
    }
    if (scenario.consistent_motion) for (size_t k=0;k<contacts.size();++k) {
        auto& edge=contacts[k];
        const auto& stable=topology_model.diagnostics().stable_ids;
        const uint32_t a=stable[edge.a],b=stable[edge.b],axis=static_cast<uint32_t>(contact_faces[k]%2);
        const auto found=activity.find({a,b,axis});
        if (contact_faces[k]<motion_history.size() && motion_history[contact_faces[k]].a==a && motion_history[contact_faces[k]].b==b) {
            edge.ridge_activity=motion_history[contact_faces[k]].ridge;
            edge.slab_activity=motion_history[contact_faces[k]].slab;
        } else if (found!=activity.end()) {
            edge.ridge_activity=found->second[0]/found->second[2];
            edge.slab_activity=found->second[1]/found->second[2];
        } else {
            // Initial activity follows prescribed mantle forcing, not random
            // constructor velocities. Newly exposed faces inherit their pair.
            const auto va=bodies[edge.a].mantle_velocity_m_s,vb=bodies[edge.b].mantle_velocity_m_s;
            const double closing=(va.x-vb.x)*edge.normal.x+(va.y-vb.y)*edge.normal.y;
            const double strength=std::abs(closing)/std::hypot(closing,edge.activation_speed_m_s);
            edge.ridge_activity=closing<0?strength:0;edge.slab_activity=closing>0?strength:0;
        }
        edge.activity_response=scenario.delta_time_myr/(scenario.boundary_memory_myr+scenario.delta_time_myr);
    }
    std::vector<Vector> old_velocity(num_plates);
    for (uint32_t p=0;p<num_plates;++p) {
        const auto v=plates[p]->linearVelocityVector();
        old_velocity[p]={v.x()*cell_m/dt_s,v.y()*cell_m/dt_s};
    }
    const auto solution = scenario.consistent_motion?solve_consistent(bodies,contacts,old_velocity):solve(bodies, contacts);
    const auto previous_diagnostics=motion_diagnostics;
    motion_diagnostics = solution.diagnostics;
    motion_diagnostics.max_velocity_consistency_cm_yr=std::max(previous_diagnostics.max_velocity_consistency_cm_yr,
        motion_diagnostics.velocity_consistency_cm_yr);
    motion_diagnostics.persistent_face_steps=previous_diagnostics.persistent_face_steps;
    motion_diagnostics.reversing_face_steps=previous_diagnostics.reversing_face_steps;
    motion_diagnostics.strong_face_steps=previous_diagnostics.strong_face_steps;
    motion_diagnostics.strong_reversing_face_steps=previous_diagnostics.strong_reversing_face_steps;
    if (getTimeMyr()>motion_history_time) {
        motion_history.resize(static_cast<size_t>(width)*height*2);
        const auto& stable=topology_model.diagnostics().stable_ids;
        for (uint32_t i=0;i<width*height;++i) for (uint32_t axis=0;axis<2;++axis) {
            const uint32_t j=axis==0?(i/width)*width+(i%width+1)%width:(i+width)%(width*height);
            const auto a=imap[i],b=imap[j];
            auto& sample=motion_history[2*i+axis];
            const double closing=(axis==0?solution.velocity_m_s[a].x-solution.velocity_m_s[b].x:
                solution.velocity_m_s[a].y-solution.velocity_m_s[b].y)*3.15576e9;
            if (a!=b && sample.a==stable[a] && sample.b==stable[b]) {
                ++motion_diagnostics.persistent_face_steps;
                const bool reversal=sample.closing*closing<0;
                motion_diagnostics.reversing_face_steps+=reversal;
                if (std::abs(sample.closing)>.5 && std::abs(closing)>.5) {
                    ++motion_diagnostics.strong_face_steps;
                    motion_diagnostics.strong_reversing_face_steps+=reversal;
                }
            }
            sample={stable[a],stable[b],closing,0,0};
        }
        if (scenario.consistent_motion) for (size_t k=0;k<contacts.size();++k) {
            const auto& edge=contacts[k];auto& sample=motion_history[contact_faces[k]];
            const double strength=std::abs(sample.closing)/std::hypot(sample.closing,scenario.boundary_force_ramp_cm_yr);
            sample.ridge=edge.ridge_activity+edge.activity_response*((sample.closing<0?strength:0)-edge.ridge_activity);
            sample.slab=edge.slab_activity+edge.activity_response*((sample.closing>0?strength:0)-edge.slab_activity);
        }
        motion_history_time=getTimeMyr();
    }
    for (uint32_t p = 0; p < num_plates; ++p) {
        const auto velocity = solution.velocity_m_s[p];
        const double dx = velocity.x * dt_s / cell_m, dy = velocity.y * dt_s / cell_m;
        // Raster transport cannot resolve skipping cells. Fail explicitly rather
        // than hiding excessive velocities behind a minimum/maximum speed clamp.
        if (!std::isfinite(dx) || !std::isfinite(dy) || std::max(std::abs(dx), std::abs(dy)) > .95)
            throw runtime_error("Force-balance displacement exceeds 0.95 cell; reduce delta_time_myr");
        // Version 1 solves planar translations, not torques or spherical rotations.
        plates[p]->setKinematics(Tectonics::FloatVector(static_cast<float>(dx), static_cast<float>(dy)), 0.0f);
    }
}

void lithosphere::initializeColumns() {
    subduction_polarity.clear();
    motion_history.clear();motion_history_time=-1;
    updateCrustAgeMyrMap();
    std::vector<tectonics::columns::Cell> cells(_worldDimension.getArea());
    for (uint32_t i = 0; i < _worldDimension.getArea(); ++i)
        cells[i] = tectonics::columns::seed(internalToPhysicalMeters(hmap[i]) - sea_level_m,
            hmap[i] >= CONTINENTAL_BASE || is_continental_material(material_map[i]), crust_age_myr_map[i]);
    column_model.initialize(_worldDimension.getWidth(), _worldDimension.getHeight(),
        scenario.world_width_km * 1000 / _worldDimension.getWidth(), std::move(cells), scenario.transport_subdivision,
        scenario.seed, {scenario.ocean_crust_variation_m, scenario.ocean_noise_wavelength_km,
        scenario.seamount_crust_m,scenario.seamount_radius_km,scenario.seamount_spacing_km,scenario.seamount_roughness});
    std::vector<uint32_t> owners(_worldDimension.getArea());
    for (uint32_t i = 0; i < owners.size(); ++i) owners[i] = imap[i];
    topology_model.initialize(_worldDimension.getWidth(), _worldDimension.getHeight(), num_plates, owners);
    if (scenario.boundary_conditioned_terrain) {
        // Use the initial force-balanced motion when available, so terrain is
        // conditioned on opening contacts that will actually drive the first step.
        solvePlateMotion();
        std::vector<double> elevations(owners.size()),ages(owners.size());
        std::vector<tectonics::transport::Displacement> velocity(owners.size());
        for (uint32_t i=0;i<owners.size();++i) {
            elevations[i]=internalToPhysicalMeters(hmap[i])-sea_level_m;
            ages[i]=crust_age_myr_map[i];
            const auto v=plates[imap[i]]->surfaceVelocityAt(i%_worldDimension.getWidth(),i/_worldDimension.getWidth());
            velocity[i]={v.x(),v.y()};
        }
        auto conditioned=tectonics::initialization::terrain(_worldDimension.getWidth(),_worldDimension.getHeight(),
            scenario.world_width_km*1000/_worldDimension.getWidth(),scenario.delta_time_myr,
            elevations,ages,owners,velocity,scenario.thermal_subsidence_strength,scenario.initial_margin_width_m,
            scenario.initial_boundary_influence);
        column_model.initialize(_worldDimension.getWidth(),_worldDimension.getHeight(),
            scenario.world_width_km*1000/_worldDimension.getWidth(),std::move(conditioned),scenario.transport_subdivision,
            scenario.seed,{scenario.ocean_crust_variation_m,scenario.ocean_noise_wavelength_km,
            scenario.seamount_crust_m,scenario.seamount_radius_km,scenario.seamount_spacing_km,scenario.seamount_roughness},true);
    }
    publishColumns();
    // Conditioning changes buoyancy/GPE. Publish motion for that actual initial
    // surface, instead of retaining the solve from the preconditioned surface.
    if (scenario.consistent_motion) solvePlateMotion();
    prev_hmap = hmap; prev_display_hmap = display_hmap;
}

void lithosphere::publishColumns() {
    const auto& cells = column_model.cells();
    for (uint32_t i = 0; i < _worldDimension.getArea(); ++i) {
        const auto& c = cells[i];
        const double z = tectonics::columns::elevation_m(c, scenario.thermal_subsidence_strength);
        if (!std::isfinite(z) || z <= -sea_level_m)
            throw runtime_error("Isostatic elevation cannot be represented by the native height datum");
        hmap[i] = elevation_to_internal(static_cast<float>(z), sea_level_m);
        const auto material = c.continental_area >= .5 ? tectonics::material::Type::Granite
                                                     : tectonics::material::Type::Basalt;
        material_map[i] = tectonics::material::to_index(material);
        const double age_steps = tectonics::columns::age_myr(c) / scenario.delta_time_myr;
        amap[i] = iter_count - std::min(iter_count, static_cast<uint32_t>(std::round(age_steps)));
        const uint32_t x = i % _worldDimension.getWidth(), y = i / _worldDimension.getWidth();
        plates[imap[i]]->setCrust(x, y, hmap[i], amap[i], material_map[i]);
    }
    display_hmap = hmap;
    const auto& uplift_rate = column_model.upliftRate();
    for (size_t i = 0; i < rock_uplift_rate_map.size(); ++i)
        rock_uplift_rate_map[i] = i < uplift_rate.size() ? static_cast<float>(uplift_rate[i]) : 0.0f;
}

void lithosphere::updateColumns() {
    prev_hmap = hmap; prev_display_hmap = display_hmap; prev_imap.copy(imap);
    solvePlateMotion();
    std::vector<tectonics::columns::Displacement> displacement(_worldDimension.getArea());
    for (uint32_t i = 0; i < _worldDimension.getArea(); ++i) {
        const auto v = plates[imap[i]]->surfaceVelocityAt(i % _worldDimension.getWidth(), i / _worldDimension.getWidth());
        displacement[i] = {v.x(), v.y()};
    }
    std::vector<uint32_t> material_owners(displacement.size());
    std::vector<double> continental_fraction(displacement.size());
    for (uint32_t i=0; i<displacement.size(); ++i) {
        material_owners[i]=imap[i]; continental_fraction[i]=column_model.cells()[i].continental_area;
    }
    auto stable_owners=material_owners;
    for (auto& owner:stable_owners) owner=topology_model.diagnostics().stable_ids[owner];
    const auto margin_forcing=tectonics::margins::forcing(_worldDimension.getWidth(),_worldDimension.getHeight(),
        scenario.world_width_km*1000/_worldDimension.getWidth(),scenario.delta_time_myr,
        stable_owners,column_model.cells(),displacement,scenario.margin_processes,&subduction_polarity,
        scenario.seed,getTimeMyr());
    displacement=tectonics::deformation::distribute(_worldDimension.getWidth(),_worldDimension.getHeight(),
        scenario.world_width_km*1000/_worldDimension.getWidth(),material_owners,continental_fraction,
        displacement,scenario.seed,{scenario.deformation_peak_km,scenario.deformation_variation,scenario.subduction_coupling,
        scenario.margin_processes.process_variation,getTimeMyr(),scenario.margin_processes.process_wavelength_km,
        scenario.margin_processes.process_period_myr});
    if (scenario.topology_processes.progressive_rifting) {
        const auto rift=topology_model.riftDisplacement(material_owners,scenario.delta_time_myr,scenario.topology_processes);
        for (size_t i=0;i<displacement.size();++i) {displacement[i].x+=rift[i].x;displacement[i].y+=rift[i].y;}
    }
    material_displacement = displacement;
    tectonics::plumes::Source plume;
    if (scenario.plume_processes.count>0)
        plume=tectonics::plumes::source(_worldDimension.getWidth(),_worldDimension.getHeight(),
            scenario.world_width_km*1000/_worldDimension.getWidth(),scenario.seed,getTimeMyr(),
            getTimeMyr()+scenario.delta_time_myr,scenario.plume_processes);
    column_model.advance(displacement, scenario.delta_time_myr,
        scenario.erosion_strength, scenario.thermal_subsidence_strength, scenario.submarine_erosion_ratio,
        &margin_forcing,scenario.margin_processes,scenario.experimental_noise,scenario.sediment_load_feedback,
        plume.basalt_m.empty()?nullptr:&plume,scenario.uplift_rate_window_myr);
    plume_influence_map.resize(plume.influence.size());
    std::transform(plume.influence.begin(),plume.influence.end(),plume_influence_map.begin(),
        [](double value){return static_cast<float>(value);});
    std::vector<tectonics::transport::Displacement> plate_displacement(num_plates);
    for (uint32_t p = 0; p < num_plates; ++p) {
        const auto v = plates[p]->linearVelocityVector();
        plate_displacement[p] = {v.x(), v.y()};
    }
    auto owners = topology_model.advance(plate_displacement);
    const double cell_m = scenario.world_width_km * 1000 / _worldDimension.getWidth();
    const double velocity_scale = cell_m / (scenario.delta_time_myr * tectonics::dynamics::kSecondsPerMyr);
    std::vector<tectonics::dynamics::Vector> mantle(owners.size()), velocity(num_plates);
    std::vector<double> continent(owners.size());
    std::vector<double> strength(owners.size());
    for (uint32_t i = 0; i < owners.size(); ++i) {
        mantle[i] = tectonics::dynamics::mantle_flow(
            (i % _worldDimension.getWidth() + .5) / _worldDimension.getWidth(),
            (i / _worldDimension.getWidth() + .5) / _worldDimension.getHeight(),
            getTimeMyr() + scenario.delta_time_myr, scenario.seed, scenario.mantle_flow_speed_cm_yr);
        continent[i] = column_model.cells()[i].continental_area;
        const auto& c=column_model.cells()[i];
        const double thickness=c.continental_area>1e-6?c.continent/c.continental_area:35000;
        const double thermal_age=c.ocean>0?std::max(0.0,(c.oceanic_age_volume-c.thermal_rejuvenation_volume)/c.ocean):60;
        const double weakening=plume_influence_map.empty()?0:
            scenario.plume_processes.weakening*plume_influence_map[i];
        strength[i]=std::clamp((continent[i]*std::pow(thickness/35000,2)+
            (1-continent[i])*std::sqrt((thermal_age+5)/65))*(1-weakening),.15,4.0);
    }
    for (uint32_t p = 0; p < num_plates; ++p)
        velocity[p] = {plate_displacement[p].x * velocity_scale, plate_displacement[p].y * velocity_scale};
    const auto changes = topology_model.evolve(owners, mantle, velocity, continent,
        cell_m, scenario.delta_time_myr, getTimeMyr() + scenario.delta_time_myr,scenario.topology_processes,strength);
    for (const auto& change : changes) if (change.kind == tectonics::topology::EventKind::Split) {
        ensurePlateCapacity(num_plates + 1);
        auto* heights = new float[owners.size()]();
        for (uint32_t i = 0; i < owners.size(); ++i) if (owners[i] == change.b) heights[i] = hmap[i];
        // These are only per-plate lookup caches. Material columns remain untouched, and
        // the newborn inherits motion until the next force solve, with no kick.
        plates[change.b] = new plate(scenario.seed, heights, _worldDimension.getWidth(),
            _worldDimension.getHeight(), 0, 0, iter_count, _worldDimension, material_map.data(),
            1.0f);
        plates[change.b]->setKinematics(plates[change.a]->linearVelocityVector(), 0);
        ++num_plates;
    }
    std::fill(plate_indices_found.begin(), plate_indices_found.end(), 0);
    for (uint32_t i = 0; i < _worldDimension.getArea(); ++i) {
        // Plate-local raster copies are caches of the material columns. Moving
        // and re-rasterizing those caches a second time caused trailing islands.
        if (imap[i] != owners[i]) plates[imap[i]]->setCrust(
            i % _worldDimension.getWidth(), i / _worldDimension.getWidth(), 0, iter_count);
        imap[i] = owners[i];
        ++plate_indices_found[imap[i]];
    }
    ++iter_count;
    publishColumns();
    removeEmptyPlates();
    if (!changes.empty()) solvePlateMotion();
    updateDerivedMaps();
    // Preserve fractional material ages rather than rounding them to updates.
    for (size_t i = 0; i < displacement.size(); ++i)
        crust_age_myr_map[i] = static_cast<float>(tectonics::columns::age_myr(column_model.cells()[i]));
    finished = iter_count - time_origin_step >= static_cast<uint32_t>(
        std::ceil(scenario.evolution_duration_myr / scenario.delta_time_myr));
}

uint32_t lithosphere::getWidth() const {
    return _worldDimension.getWidth();
}

uint32_t lithosphere::getHeight() const {
    return _worldDimension.getHeight();
}

uint32_t* lithosphere::getPlatesMap() const throw() {
    return imap.raw_data();
}

const plate* lithosphere::getPlate(uint32_t index) const {
    ASSERT(index < num_plates, "invalid plate index");
    return plates[index];
}

void lithosphere::update() {
    if (isFinished()) return;
    ++_steps;
    updateColumns();
}
