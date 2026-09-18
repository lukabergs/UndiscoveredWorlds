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

#include "lithosphere.hpp"
#include "noise.hpp"
#include "plate.hpp"
#include "simplexnoise.hpp"
#include "sqrdmd.hpp"
#include "tectonic_boundary_graph.hpp"
#include "tectonic_deforming_regions.hpp"
#include "tectonic_junction_coupling.hpp"
#include "tectonic_junctions.hpp"
#include "tectonic_provenance.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#define BOOL_REGENERATE_CRUST 1

using namespace std;

static const float SUBDUCT_RATIO = 0.5f;

static const float BUOYANCY_BONUS_X = 3;
static const uint32_t MAX_BUOYANCY_AGE = 20;
static const float MULINV_MAX_BUOYANCY_AGE = 1.0f / (float)MAX_BUOYANCY_AGE;

static const float RESTART_ENERGY_RATIO = 0.15f;
static const float RESTART_SPEED_LIMIT = 2.0f;
static const uint32_t NO_COLLISION_TIME_LIMIT = 10;

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

constexpr float kPlateBirthTensionThreshold = 100.0f;
constexpr float kPlateBirthPropagationThreshold = 72.0f;
constexpr float kPlateBirthTensionDecay = 0.985f;
constexpr float kPlateBirthTensionGain = 18.0f;
constexpr float kPlateBirthKickStrength = 0.18f;

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

platec::scenario::Scenario make_scenario(
    long seed, uint32_t width, uint32_t height, float sea_level, uint32_t erosion_period,
    float folding_ratio, uint32_t aggr_overlap_abs, float aggr_overlap_rel,
    uint32_t cycle_count, uint32_t plate_count, float erosion_strength,
    float crust_rotation_strength, float rotation_strength, float subduction_strength,
    int32_t sea_level_m_override, uint16_t initial_min_height_m,
    uint16_t initial_max_height_m, uint32_t cycle_step_limit,
    float divergent_carve_strength) {
    platec::scenario::Scenario scenario;
    scenario.seed = seed;
    scenario.width = width;
    scenario.height = height;
    scenario.sea_level = sea_level;
    scenario.erosion_period = erosion_period;
    scenario.folding_ratio = folding_ratio;
    scenario.aggregation_overlap_abs = aggr_overlap_abs;
    scenario.aggregation_overlap_rel = aggr_overlap_rel;
    scenario.cycle_count = cycle_count;
    scenario.plate_count = plate_count;
    scenario.erosion_strength = erosion_strength;
    scenario.crust_rotation_strength = crust_rotation_strength;
    scenario.rotation_strength = rotation_strength;
    scenario.subduction_strength = subduction_strength;
    scenario.sea_level_m = sea_level_m_override;
    scenario.initial_min_height_m = initial_min_height_m;
    scenario.initial_max_height_m = initial_max_height_m;
    scenario.cycle_step_limit = cycle_step_limit;
    scenario.divergent_carve_strength = divergent_carve_strength;
    return scenario;
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

GridStep make_grid_step(const Platec::FloatVector& velocity) {
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

platec::contract::CrustClass classify_crust_class(float total_thickness,
                                                  float continental_thickness) {
    using platec::contract::CrustClass;

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

platec::contract::CrustClass crust_class_from_material(uint8_t material_index) {
    using platec::contract::CrustClass;
    return platec::material::from_index(material_index) == platec::material::Type::Basalt
        ? CrustClass::Oceanic
        : CrustClass::Continental;
}

bool is_continental_material(uint8_t material_index) {
    return crust_class_from_material(material_index) == platec::contract::CrustClass::Continental;
}

float subduction_wavelet(float x, float trench_scale, float arc_scale) {
    const float scale = x < 0.0f ? trench_scale : arc_scale;
    if (scale <= FLT_EPSILON) {
        return 0.0f;
    }
    return (x / scale) * std::exp(-(x * x) / (2.0f * scale * scale));
}

uint16_t resolve_crust_type_boundary_m(const platec::scenario::Scenario& scenario,
                                       uint16_t sea_level_m) {
    if (scenario.crust_type_boundary_m != platec::scenario::kAutoCrustTypeBoundaryMeters) {
        return static_cast<uint16_t>(scenario.crust_type_boundary_m);
    }

    const uint32_t threshold =
        static_cast<uint32_t>(sea_level_m) +
        static_cast<uint32_t>(platec::scenario::kDefaultCrustTypeBoundaryOffsetMeters);
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

lithosphere::lithosphere(long seed, uint32_t width, uint32_t height, float sea_level,
                         uint32_t _erosion_period, float _folding_ratio, uint32_t aggr_ratio_abs,
                         float aggr_ratio_rel, uint32_t num_cycles, uint32_t _max_plates,
                         float _erosion_strength, float _crust_rotation_strength,
                         float _rotation_strength, float _subduction_strength,
                         int32_t sea_level_m_override, uint16_t _initial_min_height_m,
                         uint16_t _initial_max_height_m, uint32_t _cycle_step_limit,
                         float _divergent_carve_strength) noexcept(false)
    : lithosphere(make_scenario(seed, width, height, sea_level, _erosion_period,
                                _folding_ratio, aggr_ratio_abs, aggr_ratio_rel, num_cycles,
                                _max_plates, _erosion_strength, _crust_rotation_strength,
                                _rotation_strength, _subduction_strength,
                                sea_level_m_override, _initial_min_height_m,
                                _initial_max_height_m, _cycle_step_limit,
                                _divergent_carve_strength)) {}

lithosphere::lithosphere(const platec::scenario::Scenario& params) noexcept(false)
    : scenario(params), hmap(params.width, params.height), prev_hmap(params.width, params.height),
      display_hmap(params.width, params.height), prev_display_hmap(params.width, params.height),
      initial_hmap(params.width, params.height), imap(params.width, params.height),
      prev_imap(params.width, params.height), amap(params.width, params.height), plates(nullptr),
      plate_areas(params.plate_count), plate_indices_found(params.plate_count),
      aggr_overlap_abs(params.aggregation_overlap_abs),
      aggr_overlap_rel(clamp_unit(params.aggregation_overlap_rel)), cycle_count(0),
      erosion_period(params.erosion_period), folding_ratio(clamp_unit(params.folding_ratio)),
      iter_count(0), time_origin_step(0), max_cycles(params.cycle_count),
      cycle_step_limit(params.cycle_step_limit), max_plates(params.plate_count),
      plate_capacity(params.plate_count), num_plates(0),
      erosion_strength(params.erosion_strength < 0.0f ? 0.0f : params.erosion_strength),
      crust_rotation_strength(params.crust_rotation_strength < 0.0f ? 0.0f
                                                                    : params.crust_rotation_strength),
      rotation_strength(params.rotation_strength < 0.0f ? 0.0f : params.rotation_strength),
      subduction_strength(clamp_unit(params.subduction_strength)),
      divergent_carve_strength(params.divergent_carve_strength < 0.0f ? 0.0f
                                                                      : params.divergent_carve_strength),
      sea_level_m(TopographyCodec::legacy_raw_sea_level_m()),
      initial_min_height_m(params.initial_min_height_m),
      initial_max_height_m(params.initial_max_height_m),
      convergence_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      divergence_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      shear_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      geologic_regime_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          static_cast<uint8_t>(platec::contract::GeologicRegime::Stable)),
      crust_age_myr_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      crust_thickness_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          0.0f),
      crust_class_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                      static_cast<uint8_t>(platec::contract::CrustClass::None)),
      uplift_tendency_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                          0.0f),
      subsidence_tendency_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      accumulated_strain_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0.0f),
      boundary_type_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                        static_cast<uint8_t>(platec::contract::BoundaryType::None)),
      boundary_distance_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height), 0U),
      boundary_segment_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          platec::contract::kNoBoundaryId),
      nearest_boundary_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          platec::contract::kNoBoundaryId),
      deforming_region_id_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          platec::contract::kNoDeformingRegionId),
      deforming_region_type_map(
          static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
          static_cast<uint8_t>(platec::contract::DeformingRegionType::None)),
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
          static_cast<uint8_t>(platec::contract::CrustClass::Oceanic)),
      material_map(static_cast<size_t>(params.width) * static_cast<size_t>(params.height),
                   platec::material::to_index(platec::material::kDefaultType)),
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
    if (!std::isfinite(params.cycle_duration_myr) || params.cycle_duration_myr < 0.0) {
        throw runtime_error("cycle_duration_myr must be finite and >= 0");
    }
    if (!std::isfinite(params.continental_compression_gain) ||
        params.continental_compression_gain < 0.0f) {
        throw runtime_error("continental_compression_gain must be finite and >= 0");
    }
    if (!std::isfinite(params.continental_boundary_fluidity) ||
        params.continental_boundary_fluidity < 0.0f ||
        params.continental_boundary_fluidity > 1.0f) {
        throw runtime_error("continental_boundary_fluidity must be finite and in [0, 1]");
    }
    if (!std::isfinite(params.movement_energy) || params.movement_energy < 0.0f) {
        throw runtime_error("movement_energy must be finite and >= 0");
    }
    if (params.crust_type_boundary_m < platec::scenario::kAutoCrustTypeBoundaryMeters ||
        params.crust_type_boundary_m > static_cast<int32_t>(TopographyCodec::kMaxHeightMeters)) {
        throw runtime_error("crust_type_boundary_m is out of range");
    }
    collisions.resize(plate_capacity);
    subductions.resize(plate_capacity);
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
    collisions.resize(plate_capacity);
    subductions.resize(plate_capacity);
}

float lithosphere::glacialCadenceScale() const noexcept {
    const uint32_t glacial_period = scenario.glacial_erosion_period;
    if (glacial_period == 0U) {
        return 0.0f;
    }
    return std::min(1.0f, 6.0f / static_cast<float>(glacial_period));
}

float lithosphere::effectiveGlacialErosionStrength() const noexcept {
    return std::max(0.0f, scenario.glacial_erosion_strength) * glacialCadenceScale();
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
                platec::material::to_index(platec::material::kDefaultType));

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
                                  _worldDimension, pmat.data(), erosion_strength,
                                  crust_rotation_strength, rotation_strength,
                                  scenario.movement_energy);
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
        if (time_origin_step == 0U) {
            time_origin_step = iter_count;
        }
        peak_Ek = 0;
        last_coll_count = 0;
        updateDerivedMaps();
        seedInitialCrustAges();
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

const platec::contract::BoundarySegment* lithosphere::getBoundarySegments() const noexcept {
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

const platec::contract::DeformingRegion* lithosphere::getDeformingRegions() const noexcept {
    return deforming_regions.data();
}

uint32_t lithosphere::getDeformingRegionCount() const noexcept {
    return static_cast<uint32_t>(deforming_regions.size());
}

const platec::contract::TectonicJunction* lithosphere::getTectonicJunctions() const noexcept {
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

const platec::contract::PlateKinematics* lithosphere::getPlateKinematics() const noexcept {
    return plate_kinematics.data();
}

double lithosphere::getTimeMyr() const noexcept {
    const uint32_t elapsed_steps =
        iter_count >= time_origin_step ? iter_count - time_origin_step : 0U;
    return platec::scenario::steps_to_myr(elapsed_steps, scenario.delta_time_myr);
}

bool lithosphere::isOceanic(float value) const noexcept {
    return TopographyCodec::is_oceanic_internal(value);
}

void lithosphere::resetSimulationState() {
    amap.set_all(0);
    imap.set_all(0xFFFFFFFF);
    prev_imap.set_all(0xFFFFFFFF);
    prev_hmap = hmap;
    prev_display_hmap = display_hmap;
    std::fill(convergence_map.begin(), convergence_map.end(), static_cast<uint8_t>(0));
    std::fill(divergence_map.begin(), divergence_map.end(), static_cast<uint8_t>(0));
    std::fill(shear_map.begin(), shear_map.end(), static_cast<uint8_t>(0));
    std::fill(geologic_regime_map.begin(), geologic_regime_map.end(),
              static_cast<uint8_t>(platec::contract::GeologicRegime::Stable));
    std::fill(crust_age_myr_map.begin(), crust_age_myr_map.end(), 0.0f);
    std::fill(crust_thickness_map.begin(), crust_thickness_map.end(), 0.0f);
    std::fill(crust_class_map.begin(), crust_class_map.end(),
              static_cast<uint8_t>(platec::contract::CrustClass::None));
    std::fill(uplift_tendency_map.begin(), uplift_tendency_map.end(), 0.0f);
    std::fill(subsidence_tendency_map.begin(), subsidence_tendency_map.end(), 0.0f);
    std::fill(accumulated_strain_map.begin(), accumulated_strain_map.end(), 0.0f);
    std::fill(boundary_type_map.begin(), boundary_type_map.end(),
              static_cast<uint8_t>(platec::contract::BoundaryType::None));
    std::fill(boundary_distance_map.begin(), boundary_distance_map.end(), static_cast<uint16_t>(0));
    std::fill(boundary_segment_id_map.begin(), boundary_segment_id_map.end(),
              platec::contract::kNoBoundaryId);
    std::fill(nearest_boundary_id_map.begin(), nearest_boundary_id_map.end(),
              platec::contract::kNoBoundaryId);
    std::fill(deforming_region_id_map.begin(), deforming_region_id_map.end(),
              platec::contract::kNoDeformingRegionId);
    std::fill(deforming_region_type_map.begin(), deforming_region_type_map.end(),
              static_cast<uint8_t>(platec::contract::DeformingRegionType::None));
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
    next_boundary_segment_id = platec::contract::kNoBoundaryId + 1U;

    for (uint32_t i = 0; i < plate_capacity; ++i) {
        collisions[i].clear();
        subductions[i].clear();
    }

    clearPlates();

    cycle_count = 0;
    iter_count = 0;
    time_origin_step = 0;
    peak_Ek = 0;
    last_coll_count = 0;
    _steps = 0;
    finished = false;
}

void lithosphere::updateCrustAgeMyrMap() {
    const uint32_t map_area = _worldDimension.getArea();
    for (uint32_t i = 0; i < map_area; ++i) {
        const uint32_t age_steps = iter_count >= amap[i] ? iter_count - amap[i] : 0U;
        crust_age_myr_map[i] = static_cast<float>(
            platec::scenario::steps_to_myr(age_steps, scenario.delta_time_myr));
    }
}

void lithosphere::updateTectonicProvenance() {
    updateCrustAgeMyrMap();
    plate_kinematics.resize(num_plates);
    for (uint32_t plate_index = 0; plate_index < num_plates; ++plate_index) {
        const plate* current_plate = plates[plate_index];
        platec::contract::PlateKinematics& kinematics = plate_kinematics[plate_index];
        const Platec::FloatVector unit_velocity = current_plate->velocityUnitVector();
        const Platec::FloatVector linear_velocity = current_plate->linearVelocityVector();
        const FloatPoint mass_center = current_plate->worldMassCenter();

        kinematics.unit_velocity_x = unit_velocity.x();
        kinematics.unit_velocity_y = unit_velocity.y();
        kinematics.linear_velocity_x = linear_velocity.x();
        kinematics.linear_velocity_y = linear_velocity.y();
        kinematics.angular_velocity = current_plate->getAngularVelocity();
        kinematics.mass_center_x = mass_center.getX();
        kinematics.mass_center_y = mass_center.getY();
    }

    const platec::provenance::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        display_hmap.raw_data(),
        imap.raw_data(),
        plate_kinematics.data(),
        static_cast<uint32_t>(plate_kinematics.size()),
        crust_class_map.data(),
    };
    platec::provenance::compute_maps(inputs, convergence_map.data(), divergence_map.data(),
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
              platec::contract::kNoBoundaryId);

    std::deque<uint32_t> frontier;
    frontier.clear();
    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != platec::contract::kNoBoundaryId) {
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
                 boundary_id != platec::contract::kNoBoundaryId &&
                 (nearest_boundary_id_map[neighbor] == platec::contract::kNoBoundaryId ||
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

    const platec::boundary_graph::Inputs inputs{
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
    platec::boundary_graph::Outputs outputs = platec::boundary_graph::build(inputs);
    boundary_segments = std::move(outputs.segments);
    boundary_segment_id_map = std::move(outputs.boundary_segment_id);
    boundary_type_map = std::move(outputs.boundary_type);
    geologic_regime_map = std::move(outputs.boundary_regime);
    next_boundary_segment_id = outputs.next_boundary_id;

    updateBoundaryDistanceMap();

    std::vector<const platec::contract::BoundarySegment*> segments_by_id(
        next_boundary_segment_id, nullptr);
    for (const platec::contract::BoundarySegment& segment : boundary_segments) {
        if (segment.id < segments_by_id.size()) {
            segments_by_id[segment.id] = &segment;
        }
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        if (boundary_segment_id_map[index] != platec::contract::kNoBoundaryId) {
            continue;
        }

        const uint32_t nearest_boundary_id = nearest_boundary_id_map[index];
        const auto raw_regime =
            static_cast<platec::contract::GeologicRegime>(raw_geologic_regime[index]);
        if (nearest_boundary_id == platec::contract::kNoBoundaryId ||
            nearest_boundary_id >= segments_by_id.size() ||
            segments_by_id[nearest_boundary_id] == nullptr) {
            geologic_regime_map[index] =
                static_cast<uint8_t>(raw_regime == platec::contract::GeologicRegime::PassiveMargin
                                         ? raw_regime
                                         : platec::contract::GeologicRegime::Stable);
            continue;
        }

        const platec::contract::BoundarySegment& segment =
            *segments_by_id[nearest_boundary_id];
        const uint16_t distance = boundary_distance_map[index];
        auto resolved_regime = platec::contract::GeologicRegime::Stable;
        if (distance <= 1U) {
            switch (segment.boundary_type) {
            case platec::contract::BoundaryType::Convergent:
                if (raw_regime == platec::contract::GeologicRegime::TrenchAdjacent ||
                    segment.geologic_regime == platec::contract::GeologicRegime::TrenchAdjacent) {
                    resolved_regime = platec::contract::GeologicRegime::TrenchAdjacent;
                } else {
                    resolved_regime = segment.geologic_regime;
                }
                break;
            case platec::contract::BoundaryType::Divergent:
                resolved_regime = segment.geologic_regime;
                break;
            case platec::contract::BoundaryType::Transform:
                resolved_regime = platec::contract::GeologicRegime::Transform;
                break;
            case platec::contract::BoundaryType::PassiveMargin:
                resolved_regime = platec::contract::GeologicRegime::PassiveMargin;
                break;
            case platec::contract::BoundaryType::None:
            default:
                resolved_regime = platec::contract::GeologicRegime::Stable;
                break;
            }
        } else if (distance <= 3U &&
                   segment.boundary_type == platec::contract::BoundaryType::PassiveMargin) {
            resolved_regime = platec::contract::GeologicRegime::PassiveMargin;
        } else if (distance <= 2U && raw_regime != platec::contract::GeologicRegime::Stable) {
            resolved_regime = raw_regime;
        }

        // Nearby collision belts cannot change an oceanic cell into continental crust.
        if (resolved_regime == platec::contract::GeologicRegime::ContinentCollision &&
            crust_class_map[index] == static_cast<uint8_t>(platec::contract::CrustClass::Oceanic)) {
            resolved_regime = platec::contract::GeologicRegime::ConvergentArc;
        }
        geologic_regime_map[index] = static_cast<uint8_t>(resolved_regime);
    }
}

void lithosphere::updateTectonicJunctions() {
    const platec::tectonic_junctions::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        imap.raw_data(),
        boundary_segment_id_map.data(),
        boundary_segments.data(),
        static_cast<uint32_t>(boundary_segments.size()),
        plate_kinematics.data(),
        static_cast<uint32_t>(plate_kinematics.size()),
    };
    platec::tectonic_junctions::Outputs outputs = platec::tectonic_junctions::build(inputs);
    tectonic_junctions = std::move(outputs.junctions);
}

void lithosphere::updateTectonicJunctionCouplingMaps() {
    const uint32_t map_area = _worldDimension.getArea();
    const platec::tectonic_junction_coupling::Inputs inputs{
        _worldDimension.getWidth(),
        _worldDimension.getHeight(),
        tectonic_junctions.data(),
        static_cast<uint32_t>(tectonic_junctions.size()),
    };
    platec::tectonic_junction_coupling::Outputs outputs =
        platec::tectonic_junction_coupling::build(inputs);
    junction_rift_influence_map = std::move(outputs.rift_influence);
    junction_transform_influence_map = std::move(outputs.transform_influence);
    junction_subduction_influence_map = std::move(outputs.subduction_influence);

    if (outputs.rift_tension_x.size() != map_area || outputs.rift_tension_y.size() != map_area ||
        rift_tension_x_map.size() != map_area || rift_tension_y_map.size() != map_area) {
        return;
    }

    constexpr float kJunctionRiftTensionGain = kPlateBirthTensionGain * 0.35f;
    constexpr float kMaxTension = kPlateBirthTensionThreshold * 2.5f;
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
    const platec::deforming_regions::Inputs inputs{
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
    platec::deforming_regions::Outputs outputs =
        platec::deforming_regions::build(inputs);
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
                static_cast<uint8_t>(platec::contract::CrustClass::None);
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
            static_cast<platec::contract::GeologicRegime>(geologic_regime_map[index]);
        const auto boundary_type =
            static_cast<platec::contract::BoundaryType>(boundary_type_map[index]);
        const auto deforming_region_type = static_cast<platec::contract::DeformingRegionType>(
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
            platec::contract::DeformingRegionType::ContinentalRift) {
            crust_thickness_map[index] =
                std::max(display_hmap[index],
                         crust_thickness_map[index] * (1.0f - 0.12f * deformation_rate));
        } else if (deforming_region_type ==
                   platec::contract::DeformingRegionType::DiffuseCollision) {
            crust_thickness_map[index] +=
                deformation_rate * (0.12f + 0.18f * continental_ratio);
        }
        if (junction_rift > 0.0f &&
            (static_cast<platec::contract::CrustClass>(crust_class_map[index]) ==
                 platec::contract::CrustClass::Continental ||
             static_cast<platec::contract::CrustClass>(crust_class_map[index]) ==
                 platec::contract::CrustClass::Transitional)) {
            crust_thickness_map[index] =
                std::max(display_hmap[index], crust_thickness_map[index] * (1.0f - 0.045f * junction_rift));
        }

        float uplift =
            0.50f * convergence + 0.10f * shear + 0.15f * positive_surface +
            0.10f * continental_ratio + 0.08f * age_young;
        float subsidence =
            0.24f * divergence + 0.18f * negative_surface +
            0.18f * (oceanic_ratio * age_old) + 0.10f * (1.0f - continental_ratio) +
            0.06f * (boundary_type == platec::contract::BoundaryType::PassiveMargin ? 1.0f : 0.0f);

        switch (regime) {
        case platec::contract::GeologicRegime::ConvergentArc:
            uplift += 0.24f;
            break;
        case platec::contract::GeologicRegime::ContinentCollision:
            uplift += 0.30f;
            break;
        case platec::contract::GeologicRegime::TrenchAdjacent:
            uplift += 0.10f;
            subsidence += 0.16f;
            break;
        case platec::contract::GeologicRegime::DivergentRift:
            subsidence += 0.18f;
            break;
        case platec::contract::GeologicRegime::PassiveMargin:
            subsidence += 0.20f;
            break;
        case platec::contract::GeologicRegime::MidOceanRidge:
            uplift += 0.12f;
            subsidence -= 0.08f;
            break;
        case platec::contract::GeologicRegime::Transform:
            uplift += 0.05f;
            break;
        case platec::contract::GeologicRegime::Stable:
        default:
            break;
        }

        switch (deforming_region_type) {
        case platec::contract::DeformingRegionType::ContinentalRift:
            subsidence += 0.28f * deformation_rate;
            subsidence += 0.05f * deformation_velocity;
            break;
        case platec::contract::DeformingRegionType::DiffuseCollision:
            uplift += 0.30f * deformation_rate;
            uplift += 0.04f * deformation_velocity;
            subsidence -= 0.05f * deformation_rate;
            break;
        case platec::contract::DeformingRegionType::None:
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
        case platec::contract::BoundaryType::Convergent:
            strain_input += 0.10f;
            break;
        case platec::contract::BoundaryType::Divergent:
            strain_input += 0.08f;
            break;
        case platec::contract::BoundaryType::Transform:
            strain_input += 0.12f;
            break;
        case platec::contract::BoundaryType::PassiveMargin:
            strain_input += 0.03f;
            break;
        case platec::contract::BoundaryType::None:
        default:
            break;
        }
        switch (deforming_region_type) {
        case platec::contract::DeformingRegionType::ContinentalRift:
            strain_input += 0.14f * deformation_rate;
            break;
        case platec::contract::DeformingRegionType::DiffuseCollision:
            strain_input += 0.18f * deformation_rate;
            break;
        case platec::contract::DeformingRegionType::None:
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
    sea_level_m = has_override ? static_cast<uint16_t>(sea_level_m_override)
                               : TopographyCodec::legacy_raw_sea_level_m();
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

    const float initial_domain_scale =
        static_cast<float>(initial_max_height_m) /
        static_cast<float>(TopographyCodec::kMaxHeightMeters);
    const float ocean_ceiling =
        TopographyCodec::kOceanicBase +
        initial_domain_scale *
            (std::nextafter(TopographyCodec::kContinentalBase, TopographyCodec::kOceanicBase) -
             TopographyCodec::kOceanicBase);
    const float land_ceiling = TopographyCodec::kContinentalBase + initial_domain_scale;

    std::vector<float> seeded_heightmap(map_area);
    for (uint32_t i = 0; i < map_area; ++i) {
        if (is_land[i]) {
            seeded_heightmap[i] =
                remap_to_range(raw_height[i], raw_land_min, raw_land_max,
                               TopographyCodec::kContinentalBase, land_ceiling);
        } else {
            seeded_heightmap[i] =
                remap_to_range(raw_height[i], raw_ocean_min, raw_ocean_max, 0.0f, ocean_ceiling);
        }
    }

    initializeHeightMapFromFloat(seeded_heightmap.data());
}

void lithosphere::fillEnclosedInitialOceans() {
    using platec::contract::CrustClass;

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
    using platec::contract::CrustClass;

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
                platec::material::to_index(platec::material::Type::Granite);
            if (surface_crust_type_map[index] ==
                static_cast<uint8_t>(platec::contract::CrustClass::Oceanic)) {
                material_index = platec::material::to_index(platec::material::Type::Basalt);
            } else if (composite < -0.08f) {
                material_index =
                    platec::material::to_index(platec::material::Type::Sedimentary);
            } else if (composite > 0.48f) {
                material_index =
                    platec::material::to_index(platec::material::Type::Metamorphic);
            }

            material_map[index] = material_index;
        }
    }
}

void lithosphere::seedInitialCrustAges() {
    using platec::contract::BoundaryType;
    using platec::contract::CrustClass;
    using platec::contract::GeologicRegime;

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

float lithosphere::boundaryDepthInPlate(const plate& target_plate, uint32_t world_x, uint32_t world_y,
                                        uint32_t max_radius) const noexcept {
    if (target_plate.getCrust(world_x, world_y) <= FLT_EPSILON) {
        return 0.0f;
    }

    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    for (uint32_t radius = 1U; radius <= max_radius; ++radius) {
        const int radius_i = static_cast<int>(radius);
        for (int dy = -radius_i; dy <= radius_i; ++dy) {
            for (int dx = -radius_i; dx <= radius_i; ++dx) {
                if (std::max(std::abs(dx), std::abs(dy)) != radius_i) {
                    continue;
                }

                const uint32_t sample_x =
                    wrap_index(static_cast<int>(world_x) + dx, world_width);
                const uint32_t sample_y =
                    wrap_index(static_cast<int>(world_y) + dy, world_height);
                if (target_plate.getCrust(sample_x, sample_y) <= FLT_EPSILON) {
                    return static_cast<float>(radius - 1U);
                }
            }
        }
    }

    return static_cast<float>(max_radius);
}

float lithosphere::supportHeightLimitMeters(uint32_t index) const noexcept {
    const uint8_t material_index =
        scenario.use_material_map_for_height_limit && index < material_map.size()
            ? material_map[index]
            : platec::material::to_index(platec::material::kDefaultType);
    const auto& material = platec::material::properties_from_index(material_index);
    const float gravity = scenario.gravity_mps2 > 0.0f ? scenario.gravity_mps2
                                                       : platec::scenario::kDefaultGravityMps2;
    return 2.0f * (material.tensile_strength_pa / (material.density_kg_m3 * gravity));
}

float lithosphere::effectiveHeightLimitMeters(uint32_t index) const noexcept {
    const float theoretical_limit_m = supportHeightLimitMeters(index);
    const float effective_glacial_erosion = effectiveGlacialErosionStrength();
    if (effective_glacial_erosion <= 0.0f) {
        return theoretical_limit_m;
    }

    const uint32_t y = _worldDimension.yFromIndex(index);
    const float snow_line_m = snow_line_meters(y, _worldDimension.getHeight());
    const float cap_noise =
        0.90f +
        0.20f *
            (0.5f +
             0.5f * sample_wrapped_octave_noise(
                         static_cast<float>(_worldDimension.xFromIndex(index)),
                         static_cast<float>(y), _worldDimension, 2.0f, 0.55f, 0.10f, 1.4f,
                         -1.0f, 1.0f, 131.0f, 137.0f, 139.0f, 149.0f));
    const float strong_glacier_buzzsaw =
        1.0f + 0.60f * std::max(0.0f, effective_glacial_erosion - 1.0f);
    const float buzzsaw_headroom_m =
        (4200.0f /
         ((1.0f + 0.35f * effective_glacial_erosion) * strong_glacier_buzzsaw)) *
        cap_noise;
    return std::min(theoretical_limit_m, snow_line_m + buzzsaw_headroom_m);
}

float lithosphere::softenHeight(uint32_t index, float height) const noexcept {
    if (height <= TopographyCodec::kContinentalBase ||
        TopographyCodec::is_oceanic_internal(height)) {
        return std::max(0.0f, height);
    }

    const float resistance_scale_m = std::max(750.0f, effectiveHeightLimitMeters(index));
    const float relief_m =
        std::max(0.0f, internalToPhysicalMeters(height) - static_cast<float>(sea_level_m));
    const float relief_ratio = relief_m / resistance_scale_m;
    if (relief_ratio <= 0.45f) {
        return height;
    }

    const float compressed_relief_m =
        resistance_scale_m * 0.85f * std::log1p(1.7f * relief_ratio);
    const float blend = smoothstep(0.45f, 1.10f, relief_ratio);
    const float relaxed_relief_m = interpolate(relief_m, compressed_relief_m, blend);
    return TopographyCodec::kContinentalBase + metersToInternalDelta(relaxed_relief_m);
}

float lithosphere::upliftCapacity(uint32_t index, float current_height,
                                  float candidate_height) const noexcept {
    if (candidate_height <= current_height) {
        return 1.0f;
    }

    const float resistance_scale_m = std::max(750.0f, effectiveHeightLimitMeters(index));
    const float current_elevation_m =
        std::max(0.0f, internalToPhysicalMeters(current_height) - static_cast<float>(sea_level_m));
    const float candidate_elevation_m =
        std::max(0.0f, internalToPhysicalMeters(candidate_height) - static_cast<float>(sea_level_m));
    const float current_ratio = current_elevation_m / resistance_scale_m;
    const float candidate_ratio = candidate_elevation_m / resistance_scale_m;
    const float current_drag = std::log1p(4.0f * current_ratio);
    const float candidate_drag = std::log1p(4.0f * candidate_ratio);
    float uplift_factor =
        1.0f / (1.0f + 2.00f * current_drag * current_drag +
                3.20f * candidate_drag * candidate_drag);

    const float effective_glacial_erosion = effectiveGlacialErosionStrength();
    if (effective_glacial_erosion > 0.0f) {
        const uint32_t y = _worldDimension.yFromIndex(index);
        const float local_snow_line_m = snow_line_meters(y, _worldDimension.getHeight());
        const float glacial_excess_m = std::max(0.0f, candidate_elevation_m - local_snow_line_m);
        uplift_factor /=
            1.0f + 2.0f * effective_glacial_erosion *
                        std::log1p(glacial_excess_m / 350.0f);
    }

    return std::max(0.001f, clamp_unit(uplift_factor));
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

    initializeHeightMapFromMetric(metric_heightmap.data(), TopographyCodec::legacy_raw_sea_level_m());
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
void lithosphere::resolveJuxtapositions(const uint32_t& i, const uint32_t& j, const uint32_t& k,
                                        const uint32_t& x_mod, const uint32_t& y_mod,
                                        const float*& this_map, const uint32_t*& this_age,
                                        const uint8_t*& this_material,
                                        uint32_t& continental_collisions) {
    ASSERT(i < num_plates, "Given invalid plate index");
    const uint32_t current_owner = imap[k];
    uint32_t surface_owner = current_owner;
    if (prev_imap[k] == i || prev_imap[k] == current_owner) {
        surface_owner = prev_imap[k];
    }
    const uint32_t penetrating_owner = surface_owner == i ? current_owner : i;

    plates[surface_owner]->addCollision(x_mod, y_mod);
    plates[penetrating_owner]->addCollision(x_mod, y_mod);

    const float surface_height = surface_owner == current_owner ? hmap[k] : this_map[j];
    const uint32_t surface_age = surface_owner == current_owner ? amap[k] : this_age[j];
    const uint8_t surface_material =
        surface_owner == current_owner ? material_map[k] : this_material[j];
    const float penetrating_height =
        penetrating_owner == current_owner ? hmap[k] : this_map[j];
    const uint32_t penetrating_age =
        penetrating_owner == current_owner ? amap[k] : this_age[j];

    const Platec::FloatVector surface_velocity =
        plates[surface_owner]->surfaceVelocityAt(x_mod, y_mod);
    const Platec::FloatVector penetrating_velocity =
        plates[penetrating_owner]->surfaceVelocityAt(x_mod, y_mod);
    const Platec::FloatVector relative_velocity =
        Platec::FloatVector(surface_velocity.x() - penetrating_velocity.x(),
                            surface_velocity.y() - penetrating_velocity.y());
    const float relative_speed = relative_velocity.length();
    const float surface_boundary_depth =
        boundaryDepthInPlate(*plates[surface_owner], x_mod, y_mod);
    const float penetrating_boundary_depth =
        boundaryDepthInPlate(*plates[penetrating_owner], x_mod, y_mod);
    const float uplift_ratio = clamp_unit(scenario.collision_uplift_ratio);
    const float compression_gain = std::max(0.0f, scenario.continental_compression_gain);
    const float boundary_fluidity = clamp_unit(scenario.continental_boundary_fluidity);
    const float boundary_rigidity = 1.0f - boundary_fluidity;
    const float lock_depth = static_cast<float>(continental_lock_map[k]);
    const uint8_t lock_increment =
        static_cast<uint8_t>(1U + static_cast<uint32_t>(std::lround(3.0f * boundary_rigidity)));
    continental_lock_map[k] = std::min<uint8_t>(
        12U, static_cast<uint8_t>(continental_lock_map[k] + lock_increment));

    const auto inward_scale = [&](float boundary_depth) {
        return std::exp(-(0.72f + 0.90f * boundary_rigidity) *
                        (boundary_depth + (0.25f + 0.65f * boundary_rigidity) * lock_depth));
    };
    const auto compression_noise = [&]() {
        return 0.88f + 0.24f * _randsource.next_float();
    };
    const auto find_anchor = [&](uint32_t owner, const Platec::FloatVector& velocity,
                                 bool include_boundary_cell) {
        uint32_t anchor_x = x_mod;
        uint32_t anchor_y = y_mod;
        uint32_t best_x = anchor_x;
        uint32_t best_y = anchor_y;
        const plate& target_plate = *plates[owner];
        bool found = include_boundary_cell &&
                     target_plate.getCrust(anchor_x, anchor_y) > FLT_EPSILON;
        const GridStep step = make_grid_step(velocity);
        const uint32_t walk_limit = 1U + continental_lock_map[k] / 3U +
                                    static_cast<uint32_t>(std::lround(3.0f * boundary_fluidity));

        for (uint32_t step_index = 0; step_index < walk_limit; ++step_index) {
            advance_wrapped(anchor_x, anchor_y, -step.dx, -step.dy, _worldDimension);
            if (target_plate.getCrust(anchor_x, anchor_y) > FLT_EPSILON) {
                best_x = anchor_x;
                best_y = anchor_y;
                found = true;
            }
        }

        if (!found) {
            best_x = x_mod;
            best_y = y_mod;
            if (!include_boundary_cell) {
                const GridStep fallback_step = make_grid_step(velocity);
                advance_wrapped(best_x, best_y, -fallback_step.dx, -fallback_step.dy,
                                _worldDimension);
            }
        }
        return std::pair<uint32_t, uint32_t>(best_x, best_y);
    };
    const auto reflect_world_cell = [&](uint32_t owner, uint32_t world_x, uint32_t world_y,
                                        float height, uint8_t material) {
        if (height <= FLT_EPSILON) {
            return;
        }
        const uint32_t world_index = _worldDimension.indexOf(world_x, world_y);
        if (imap[world_index] >= num_plates || imap[world_index] == owner ||
            height > hmap[world_index] + FLT_EPSILON) {
            hmap[world_index] = height;
            imap[world_index] = owner;
            amap[world_index] = plates[owner]->getCrustTimestamp(world_x, world_y);
            material_map[world_index] = material;
        }
    };

    const float speed_term = 0.75f + 0.35f * std::min(3.0f, relative_speed);
    const float surface_noise = compression_noise();
    const float surface_height_m = std::max(0.0f, internalToPhysicalMeters(surface_height));
    const float penetrating_height_m = std::max(0.0f, internalToPhysicalMeters(penetrating_height));
    const float boundary_depth =
        std::min(surface_boundary_depth, penetrating_boundary_depth);
    const float depth_factor = 0.45f + 0.55f * inward_scale(boundary_depth);
    const float consumed_crust_gain_m =
        uplift_ratio * (penetrating_height_m + 0.35f * surface_height_m) *
        speed_term * depth_factor * (0.55f + 0.45f * boundary_rigidity);
    const float uplifted_surface = softenHeight(
        k, surface_height + metersToInternalDelta(consumed_crust_gain_m * surface_noise) +
               metersToInternalDelta(_randsource.next_float() * 10.0f * uplift_ratio));

    const auto surface_anchor = find_anchor(surface_owner, surface_velocity, true);
    const uint32_t surface_anchor_index =
        _worldDimension.indexOf(surface_anchor.first, surface_anchor.second);
    const float surface_anchor_height =
        plates[surface_owner]->getCrust(surface_anchor.first, surface_anchor.second);

    const float surface_interior_gain =
        metersToInternalDelta(surface_height_m * compression_gain * uplift_ratio *
                              (0.35f + 0.65f * boundary_fluidity) * speed_term *
                              depth_factor);

    const float uplifted_surface_anchor = softenHeight(
        surface_anchor_index,
        surface_anchor_height + surface_interior_gain * surface_noise +
            metersToInternalDelta(_randsource.next_float() * 4.0f * uplift_ratio));

    plates[surface_owner]->setCrust(x_mod, y_mod, uplifted_surface,
                                    std::max(surface_age, penetrating_age), surface_material);
    plates[penetrating_owner]->setCrust(x_mod, y_mod, 0.0f, penetrating_age);

    if (surface_anchor.first != x_mod || surface_anchor.second != y_mod) {
        plates[surface_owner]->setCrust(surface_anchor.first, surface_anchor.second,
                                        uplifted_surface_anchor,
                                        std::max(surface_age, penetrating_age),
                                        surface_material);
        reflect_world_cell(surface_owner, surface_anchor.first, surface_anchor.second,
                           uplifted_surface_anchor, surface_material);
    }

    hmap[k] = uplifted_surface;
    imap[k] = surface_owner;
    amap[k] = plates[surface_owner]->getCrustTimestamp(x_mod, y_mod);
    material_map[k] = surface_material;

    const float surface_anchor_gain =
        (surface_anchor.first == x_mod && surface_anchor.second == y_mod)
            ? 0.0f
            : std::max(0.0f, uplifted_surface_anchor - surface_anchor_height);
    const float deformed_crust =
        std::max(0.0f, uplifted_surface - surface_height) + surface_anchor_gain +
        std::max(0.0f, penetrating_height);

    collisions[penetrating_owner].push_back(
        plateCollision(surface_owner, x_mod, y_mod, deformed_crust));
    ++continental_collisions;
}

// Update height and plate index maps.
// Doing it plate by plate is much faster than doing it index wise:
// Each plate's map's memory area is accessed sequentially and only
// once as opposed to calculating "num_plates" indices within plate
// maps in order to find out which plate(s) own current location.
void lithosphere::updateHeightAndPlateIndexMaps(uint32_t& oceanic_collisions,
                                                uint32_t& continental_collisions) {
    uint32_t world_width = _worldDimension.getWidth();
    uint32_t world_height = _worldDimension.getHeight();
    std::vector<uint32_t> plate_order(num_plates);
    for (uint32_t i = 0; i < num_plates; ++i) {
        plate_order[i] = i;
    }
    std::stable_sort(plate_order.begin(), plate_order.end(),
                     [&](uint32_t lhs, uint32_t rhs) {
                         const float lhs_buoyancy = plates[lhs]->buoyancy();
                         const float rhs_buoyancy = plates[rhs]->buoyancy();
                         if (std::fabs(lhs_buoyancy - rhs_buoyancy) > 1.0e-6f) {
                             return lhs_buoyancy < rhs_buoyancy;
                         }

                         const float lhs_score = plates[lhs]->continentalityScore();
                         const float rhs_score = plates[rhs]->continentalityScore();
                         if (std::fabs(lhs_score - rhs_score) > 1.0e-6f) {
                             return lhs_score < rhs_score;
                         }
                         return lhs < rhs;
                     });
    for (uint8_t& lock_value : continental_lock_map) {
        if (lock_value > 0U) {
            --lock_value;
        }
    }
    hmap.set_all(0);
    imap.set_all(0xFFFFFFFF);
    std::fill(material_map.begin(), material_map.end(),
              platec::material::to_index(platec::material::kDefaultType));
    for (uint32_t order_index = 0; order_index < num_plates; ++order_index) {
        const uint32_t i = plate_order[order_index];
        const uint32_t x0 = plates[i]->getLeftAsUint();
        const uint32_t y0 = plates[i]->getTopAsUint();
        const uint32_t x1 = x0 + plates[i]->getWidth();
        const uint32_t y1 = y0 + plates[i]->getHeight();

        const float* this_map;
        const uint32_t* this_age;
        const uint8_t* this_material;
        plates[i]->getMap(&this_map, &this_age, &this_material);
        const size_t plate_cell_count =
            static_cast<size_t>(x1 - x0) * static_cast<size_t>(y1 - y0);
        const std::vector<float> this_map_snapshot(this_map, this_map + plate_cell_count);
        const std::vector<uint32_t> this_age_snapshot(this_age, this_age + plate_cell_count);
        const std::vector<uint8_t> this_material_snapshot(
            this_material, this_material + plate_cell_count);
        this_map = this_map_snapshot.data();
        this_age = this_age_snapshot.data();
        this_material = this_material_snapshot.data();

        const uint32_t x_mod_start = x0 % world_width;

        // Copy first part of plate onto world map.
        // MK: These loops are ugly, but using modulus in here is a hog
        for (uint32_t y = y0, j = 0; y < y1; ++y) {
            const uint32_t y_mod = y % world_height;
            const uint32_t y_width = y_mod * world_width;
            uint32_t x_mod = x_mod_start;

            for (uint32_t x = x0; x < x1;
                 ++x, ++j, x_mod = ++x_mod >= world_width ? x_mod - world_width : x_mod) {
                const uint32_t k = x_mod + y_width;

                if (this_map[j] < 2 * FLT_EPSILON) // No crust here...
                    continue;

                if (imap[k] >= num_plates) // No one here yet?
                {
                    // This plate becomes the "owner" of current location
                    // if it is the first plate to have crust on it.
                    hmap[k] = this_map[j];
                    imap[k] = i;
                    amap[k] = this_age[j];
                    material_map[k] = this_material[j];

                    continue;
                }

                const uint32_t current_owner = imap[k];
                const float current_height = hmap[k];
                const uint32_t current_age = amap[k];
                const uint8_t current_material = material_map[k];
                const float incoming_height = this_map[j];
                const uint32_t incoming_age = this_age[j];
                const uint8_t incoming_material = this_material[j];
                const auto current_crust_class = crust_class_from_material(current_material);
                const auto incoming_crust_class = crust_class_from_material(incoming_material);
                const bool current_continental = is_continental_material(current_material);
                const bool incoming_continental = is_continental_material(incoming_material);

                if (current_owner == i) {
                    if (incoming_height > current_height + FLT_EPSILON ||
                        (std::fabs(incoming_height - current_height) <= FLT_EPSILON &&
                         incoming_age >= current_age)) {
                        hmap[k] = incoming_height;
                        amap[k] = incoming_age;
                        if (!(is_continental_material(current_material) &&
                              !is_continental_material(incoming_material))) {
                            material_map[k] = incoming_material;
                        }
                    }
                    continue;
                }

                plate& current_plate = *plates[current_owner];
                plate& incoming_plate = *plates[i];

                if (current_continental && incoming_continental) {
                    resolveJuxtapositions(i, j, k, x_mod, y_mod, this_map, this_age, this_material,
                                          continental_collisions);
                    continue;
                }

                uint32_t overriding_owner = current_owner;
                uint32_t subducting_owner = i;

                if (!current_continental && !incoming_continental) {
                    if (current_age > incoming_age ||
                        (current_age == incoming_age &&
                         (current_plate.buoyancy() + FLT_EPSILON < incoming_plate.buoyancy() ||
                          (std::fabs(current_plate.buoyancy() - incoming_plate.buoyancy()) <=
                               FLT_EPSILON &&
                           (current_height > incoming_height + FLT_EPSILON ||
                            (std::fabs(current_height - incoming_height) <= FLT_EPSILON &&
                             current_owner > i)))))) {
                        overriding_owner = i;
                        subducting_owner = current_owner;
                    }
                } else if (current_continental) {
                    overriding_owner = current_owner;
                    subducting_owner = i;
                } else {
                    overriding_owner = i;
                    subducting_owner = current_owner;
                }

                const float overriding_height =
                    overriding_owner == current_owner ? current_height : incoming_height;
                const uint32_t overriding_age =
                    overriding_owner == current_owner ? current_age : incoming_age;
                const uint8_t overriding_material =
                    overriding_owner == current_owner ? current_material : incoming_material;
                const auto overriding_crust_class =
                    overriding_owner == current_owner ? current_crust_class : incoming_crust_class;
                const float subducting_height =
                    subducting_owner == current_owner ? current_height : incoming_height;
                const uint32_t subducting_age =
                    subducting_owner == current_owner ? current_age : incoming_age;
                const auto subducting_crust_class =
                    subducting_owner == current_owner ? current_crust_class : incoming_crust_class;
                const uint8_t overriding_surface_material =
                    overriding_crust_class == platec::contract::CrustClass::Oceanic
                    ? overriding_material
                    : (is_continental_material(overriding_material)
                           ? overriding_material
                           : platec::material::to_index(platec::material::Type::Granite));
                current_plate.addCollision(x_mod, y_mod);
                incoming_plate.addCollision(x_mod, y_mod);

                const Platec::FloatVector overriding_velocity =
                    plates[overriding_owner]->surfaceVelocityAt(x_mod, y_mod);
                const Platec::FloatVector subducting_velocity =
                    plates[subducting_owner]->surfaceVelocityAt(x_mod, y_mod);
                const float uplift_ratio = clamp_unit(scenario.collision_uplift_ratio);
                const float subduction_ratio = clamp_unit(subduction_strength);
                const float overriding_depth =
                    boundaryDepthInPlate(*plates[overriding_owner], x_mod, y_mod);
                const float subducting_depth =
                    boundaryDepthInPlate(*plates[subducting_owner], x_mod, y_mod);
                const Platec::FloatVector relative_velocity =
                    Platec::FloatVector(subducting_velocity.x() - overriding_velocity.x(),
                                        subducting_velocity.y() - overriding_velocity.y());
                const float relative_speed = relative_velocity.length();
                const float uplift_scale = uplift_ratio * (0.35f + 2.15f * uplift_ratio);
                const float uplift_noise =
                    metersToInternalDelta(_randsource.next_float() * 6.0f * uplift_ratio);
                const bool subducting_oceanic =
                    subducting_crust_class == platec::contract::CrustClass::Oceanic;
                const bool overriding_oceanic =
                    overriding_crust_class == platec::contract::CrustClass::Oceanic;
                const float subduction_drive = clamp_unit(
                    subduction_ratio *
                    (0.30f + 0.18f * std::max(1.0f, subducting_depth) + 0.22f * relative_speed));
                float boundary_candidate = overriding_height;
                if (subducting_oceanic) {
                    const float trench_profile =
                        std::abs(std::min(0.0f, subduction_wavelet(-1.5f, 2.2f, 5.0f)));
                    const float trench_cut = subducting_height *
                        clamp_unit((0.20f + 0.62f * subduction_drive) * trench_profile * 1.35f);
                    boundary_candidate = overriding_oceanic
                        ? std::max(0.0f, overriding_height - trench_cut)
                        : overriding_height + uplift_noise * 0.05f;
                } else {
                    boundary_candidate =
                        overriding_height *
                            (1.0f + overriding_velocity.length() * overriding_depth * uplift_scale) +
                        uplift_noise;
                }
                const float uplifted_overriding = softenHeight(k, boundary_candidate);
                const float realized_uplift =
                    std::max(0.0f, uplifted_overriding - overriding_height);
                const float subducted_crust = subducting_oceanic
                    ? subducting_height * clamp_unit(0.55f + 0.45f * subduction_drive)
                    : subducting_height * subduction_drive;
                const float arc_fraction = subducting_oceanic
                    ? (overriding_oceanic ? 0.34f : 0.42f)
                    : 0.08f;
                const float arc_crust = subducted_crust * arc_fraction;

                plates[overriding_owner]->setCrust(
                    x_mod, y_mod, uplifted_overriding,
                    std::max(overriding_age, subducting_age), overriding_surface_material);
                plates[subducting_owner]->setCrust(x_mod, y_mod, 0.0f, subducting_age);

                hmap[k] = uplifted_overriding;
                imap[k] = overriding_owner;
                amap[k] = plates[overriding_owner]->getCrustTimestamp(x_mod, y_mod);
                material_map[k] = overriding_surface_material;

                if (arc_crust > FLT_EPSILON) {
                    subductions[overriding_owner].push_back(
                        plateCollision(subducting_owner, x_mod, y_mod, arc_crust));
                }
                collisions[subducting_owner].push_back(
                    plateCollision(overriding_owner, x_mod, y_mod,
                                   realized_uplift + subducted_crust));
                ++oceanic_collisions;
                ++continental_collisions;
            }
        }
    }
}

void lithosphere::updateCollisions() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        for (uint32_t j = 0; j < collisions[i].size(); ++j) {
            const plateCollision& coll = collisions[i][j];
            uint32_t coll_count, coll_count_i, coll_count_j;
            float coll_ratio, coll_ratio_i, coll_ratio_j;

            ASSERT(i != coll.index, "when colliding: SRC == DEST!");

            // Collision causes friction. Apply it to both plates.
            const bool continental_contact =
                is_continental_material(plates[i]->getMaterial(coll.wx, coll.wy)) &&
                is_continental_material(plates[coll.index]->getMaterial(coll.wx, coll.wy));
            const float contact_lock = continental_contact ? 1.0f : 1.0f;
            const float friction_mass = coll.crust * contact_lock;
            plates[i]->applyFriction(friction_mass);
            plates[coll.index]->applyFriction(friction_mass);

            plates[i]->getCollisionInfo(coll.wx, coll.wy, &coll_count_i, &coll_ratio_i);
            plates[coll.index]->getCollisionInfo(coll.wx, coll.wy, &coll_count_j, &coll_ratio_j);

            // Find the minimum count of collisions between two
            // continents on different plates.
            // It's minimum because large plate will get collisions
            // from all over whereas smaller plate will get just
            // a few. It's those few that matter between these two
            // plates, not what the big plate has with all the
            // other plates around it.
            coll_count = coll_count_i;
            coll_count -= (coll_count - coll_count_j) & -(coll_count > coll_count_j);

            // Find maximum amount of collided surface area between
            // two continents on different plates.
            // Like earlier, it's the "experience" of the smaller
            // plate that matters here.
            coll_ratio = coll_ratio_i;
            coll_ratio += (coll_ratio_j - coll_ratio) * (coll_ratio_j > coll_ratio);

            if ((coll_count > aggr_overlap_abs) | (coll_ratio > aggr_overlap_rel)) {
                plates[coll.index]->collide(*plates[i], coll.crust);
            }
        }

        collisions[i].clear();
    }
}

void lithosphere::resurfaceSubmergedCrust() {
    static constexpr uint32_t kResurfaceDelaySteps = 2U;
    static constexpr float kResurfaceFraction = 0.12f;

    const float min_resurface = std::max(0.005f, metersToInternalDelta(35.0f));
    for (uint32_t i = 0; i < num_plates; ++i) {
        if (!plates[i]->hasSubmergedCrust()) {
            continue;
        }

        const uint32_t x0 = plates[i]->getLeftAsUint();
        const uint32_t y0 = plates[i]->getTopAsUint();
        const uint32_t x1 = x0 + plates[i]->getWidth();
        const uint32_t y1 = y0 + plates[i]->getHeight();

        for (uint32_t y = y0; y < y1; ++y) {
            const uint32_t y_mod = _worldDimension.yMod(y);
            for (uint32_t x = x0; x < x1; ++x) {
                const uint32_t x_mod = _worldDimension.xMod(x);
                const float submerged = plates[i]->getSubmergedCrust(x_mod, y_mod);
                if (submerged <= FLT_EPSILON) {
                    continue;
                }

                const uint32_t hidden_since = plates[i]->getSubmergedSince(x_mod, y_mod);
                if (iter_count <= hidden_since + kResurfaceDelaySteps) {
                    continue;
                }

                const uint32_t index = _worldDimension.indexOf(x_mod, y_mod);
                if (imap[index] < num_plates && imap[index] != i) {
                    continue;
                }

                const float resurfaced = plates[i]->resurfaceSubmergedCrust(
                    x_mod, y_mod, kResurfaceFraction, min_resurface, iter_count);
                if (resurfaced <= FLT_EPSILON) {
                    continue;
                }

                const float surfaced_height = plates[i]->getCrust(x_mod, y_mod);
                if (surfaced_height <= FLT_EPSILON) {
                    continue;
                }

                if (imap[index] >= num_plates || imap[index] == i ||
                    surfaced_height > hmap[index] + FLT_EPSILON) {
                    hmap[index] = surfaced_height;
                    imap[index] = i;
                    amap[index] = plates[i]->getCrustTimestamp(x_mod, y_mod);
                    material_map[index] = plates[i]->getMaterial(x_mod, y_mod);
                }
            }
        }
    }
}

uint32_t lithosphere::chooseDivergentOwner(uint32_t x, uint32_t y, uint32_t index) const {
    std::vector<uint32_t> candidates;
    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    const uint32_t left_x = x > 0 ? x - 1 : world_width - 1;
    const uint32_t right_x = x + 1 < world_width ? x + 1 : 0;
    const uint32_t top_y = y > 0 ? y - 1 : world_height - 1;
    const uint32_t bottom_y = y + 1 < world_height ? y + 1 : 0;
    const uint32_t top_left = _worldDimension.indexOf(left_x, top_y);
    const uint32_t top_right = _worldDimension.indexOf(right_x, top_y);
    const uint32_t bottom_left = _worldDimension.indexOf(left_x, bottom_y);
    const uint32_t bottom_right = _worldDimension.indexOf(right_x, bottom_y);
    const uint32_t left = _worldDimension.indexOf(left_x, y);
    const uint32_t right = _worldDimension.indexOf(right_x, y);
    const uint32_t top = _worldDimension.indexOf(x, top_y);
    const uint32_t bottom = _worldDimension.indexOf(x, bottom_y);
    struct NeighborWeight {
        uint32_t index;
        float prev_weight;
        float current_weight;
    };
    const NeighborWeight neighbors[] = {
        {left, 0.75f, 2.6f},
        {right, 0.75f, 2.6f},
        {top, 0.75f, 2.6f},
        {bottom, 0.75f, 2.6f},
        {top_left, 0.35f, 1.5f},
        {top_right, 0.35f, 1.5f},
        {bottom_left, 0.35f, 1.5f},
        {bottom_right, 0.35f, 1.5f},
    };

    push_unique_owner(candidates, prev_imap[index], num_plates);
    for (const NeighborWeight& neighbor : neighbors) {
        push_unique_owner(candidates, prev_imap[neighbor.index], num_plates);
        push_unique_owner(candidates, imap[neighbor.index], num_plates);
    }

    if (candidates.empty()) {
        return prev_imap[index];
    }

    uint32_t best_owner = candidates[0];
    float best_score = -FLT_MAX;

    for (uint32_t owner : candidates) {
        float score = 0.0f;
        if (prev_imap[index] == owner) {
            score += 0.4f;
        }

        for (const NeighborWeight& neighbor : neighbors) {
            if (prev_imap[neighbor.index] == owner) {
                score += neighbor.prev_weight;
            }
            if (imap[neighbor.index] == owner) {
                score += neighbor.current_weight;
                if (amap[neighbor.index] == iter_count) {
                    score += neighbor.current_weight * 1.35f;
                }
            }
        }

        const float owner_noise = sample_wrapped_octave_noise(
            static_cast<float>(x), static_cast<float>(y), _worldDimension, 4.0f, 0.58f, 0.085f,
            1.5f, -1.0f, 1.0f, static_cast<float>(owner) * 13.0f + 17.0f,
            static_cast<float>(owner) * 7.0f + 29.0f, static_cast<float>(iter_count) * 0.09f,
            static_cast<float>(owner) * 19.0f + 11.0f);
        const float warp_noise = sample_wrapped_octave_noise(
            static_cast<float>(x), static_cast<float>(y), _worldDimension, 2.0f, 0.5f, 0.18f,
            2.2f, -0.75f, 0.75f, 5.0f, 13.0f,
            static_cast<float>(iter_count) * 0.05f + static_cast<float>(owner) * 0.37f,
            static_cast<float>(owner) * 31.0f + 43.0f);
        score += owner_noise * 2.1f + warp_noise * 1.5f;

        if (score > best_score) {
            best_score = score;
            best_owner = owner;
        }
    }

    return best_owner;
}

bool lithosphere::hasAssignedOwnerNeighbor(uint32_t x, uint32_t y) const {
    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    const uint32_t left_x = x > 0 ? x - 1 : world_width - 1;
    const uint32_t right_x = x + 1 < world_width ? x + 1 : 0;
    const uint32_t top_y = y > 0 ? y - 1 : world_height - 1;
    const uint32_t bottom_y = y + 1 < world_height ? y + 1 : 0;
    const uint32_t neighbors[] = {
        _worldDimension.indexOf(left_x, y),
        _worldDimension.indexOf(right_x, y),
        _worldDimension.indexOf(x, top_y),
        _worldDimension.indexOf(x, bottom_y),
        _worldDimension.indexOf(left_x, top_y),
        _worldDimension.indexOf(right_x, top_y),
        _worldDimension.indexOf(left_x, bottom_y),
        _worldDimension.indexOf(right_x, bottom_y),
    };

    for (uint32_t neighbor : neighbors) {
        if (imap[neighbor] < num_plates) {
            return true;
        }
    }

    return false;
}

void lithosphere::regenerateCrust() {
    const uint32_t map_area = _worldDimension.getArea();
    std::vector<uint32_t> frontier;
    std::vector<uint8_t> queued(map_area, 0);
    frontier.reserve(map_area / 32);

    auto enqueue = [&](uint32_t index, std::vector<uint32_t>& target) {
        if (imap[index] >= num_plates && !queued[index]) {
            target.push_back(index);
            queued[index] = 1;
        }
    };

    auto average_neighbor_height = [&](uint32_t x, uint32_t y) {
        const uint32_t world_width = _worldDimension.getWidth();
        const uint32_t world_height = _worldDimension.getHeight();
        const uint32_t left_x = x > 0 ? x - 1 : world_width - 1;
        const uint32_t right_x = x + 1 < world_width ? x + 1 : 0;
        const uint32_t top_y = y > 0 ? y - 1 : world_height - 1;
        const uint32_t bottom_y = y + 1 < world_height ? y + 1 : 0;
        const uint32_t neighbors[] = {
            _worldDimension.indexOf(left_x, y),
            _worldDimension.indexOf(right_x, y),
            _worldDimension.indexOf(x, top_y),
            _worldDimension.indexOf(x, bottom_y),
            _worldDimension.indexOf(left_x, top_y),
            _worldDimension.indexOf(right_x, top_y),
            _worldDimension.indexOf(left_x, bottom_y),
            _worldDimension.indexOf(right_x, bottom_y),
        };

        float sum = 0.0f;
        uint32_t count = 0U;
        for (uint32_t neighbor : neighbors) {
            if (hmap[neighbor] > FLT_EPSILON) {
                sum += hmap[neighbor];
                ++count;
            }
        }

        return count > 0U ? sum / static_cast<float>(count) : 0.0f;
    };

    auto choose_neighbor_owner = [&](uint32_t x, uint32_t y) {
        std::vector<uint32_t> owner_counts(num_plates, 0U);
        const uint32_t world_width = _worldDimension.getWidth();
        const uint32_t world_height = _worldDimension.getHeight();
        const uint32_t left_x = x > 0 ? x - 1 : world_width - 1;
        const uint32_t right_x = x + 1 < world_width ? x + 1 : 0;
        const uint32_t top_y = y > 0 ? y - 1 : world_height - 1;
        const uint32_t bottom_y = y + 1 < world_height ? y + 1 : 0;
        const uint32_t neighbors[] = {
            _worldDimension.indexOf(left_x, y),
            _worldDimension.indexOf(right_x, y),
            _worldDimension.indexOf(x, top_y),
            _worldDimension.indexOf(x, bottom_y),
            _worldDimension.indexOf(left_x, top_y),
            _worldDimension.indexOf(right_x, top_y),
            _worldDimension.indexOf(left_x, bottom_y),
            _worldDimension.indexOf(right_x, bottom_y),
        };

        for (uint32_t neighbor : neighbors) {
            const uint32_t owner = imap[neighbor];
            if (owner < num_plates) {
                ++owner_counts[owner];
            }
        }

        uint32_t best_owner = std::numeric_limits<uint32_t>::max();
        uint32_t best_count = 0U;
        for (uint32_t owner = 0; owner < num_plates; ++owner) {
            if (owner_counts[owner] > best_count) {
                best_count = owner_counts[owner];
                best_owner = owner;
            }
        }
        return best_owner;
    };

    auto assign_new_crust = [&](uint32_t index) {
        if (imap[index] < num_plates) {
            return;
        }

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        imap[index] = chooseDivergentOwner(x, y, index);
        amap[index] = iter_count;
        if (imap[index] >= num_plates) {
            imap[index] = choose_neighbor_owner(x, y);
            if (imap[index] >= num_plates) {
                const float neighbor_height = average_neighbor_height(x, y);
                if (neighbor_height > FLT_EPSILON) {
                    hmap[index] = neighbor_height;
                }
                return;
            }
        }

        struct BoundaryOwnerSample {
            uint32_t owner;
            float dir_x;
            float dir_y;
            float dir_weight;
            float retreat;
            float velocity_x;
            float velocity_y;
        };

        auto add_owner_sample = [&](std::vector<BoundaryOwnerSample>& samples, uint32_t owner,
                                    float dir_x, float dir_y, float weight) {
            if (owner >= num_plates || weight <= FLT_EPSILON) {
                return;
            }

            const float length = std::sqrt(dir_x * dir_x + dir_y * dir_y);
            if (length <= FLT_EPSILON) {
                return;
            }

            const float nx = dir_x / length;
            const float ny = dir_y / length;
            for (BoundaryOwnerSample& sample : samples) {
                if (sample.owner == owner) {
                    sample.dir_x += nx * weight;
                    sample.dir_y += ny * weight;
                    sample.dir_weight += weight;
                    return;
                }
            }

            BoundaryOwnerSample sample{};
            sample.owner = owner;
            sample.dir_x = nx * weight;
            sample.dir_y = ny * weight;
            sample.dir_weight = weight;
            sample.retreat = 0.0f;
            sample.velocity_x = 0.0f;
            sample.velocity_y = 0.0f;
            samples.push_back(sample);
        };

        auto neighbor_direction = [&](uint32_t neighbor_x, uint32_t neighbor_y, float& dx,
                                      float& dy) {
            dx = wrapped_delta(static_cast<float>(neighbor_x), static_cast<float>(x),
                               static_cast<float>(_worldDimension.getWidth()));
            dy = wrapped_delta(static_cast<float>(neighbor_y), static_cast<float>(y),
                               static_cast<float>(_worldDimension.getHeight()));
        };

        std::vector<BoundaryOwnerSample> boundary_samples;
        boundary_samples.reserve(8);

        const uint32_t world_width = _worldDimension.getWidth();
        const uint32_t world_height = _worldDimension.getHeight();
        const uint32_t left_x = x > 0 ? x - 1 : world_width - 1;
        const uint32_t right_x = x + 1 < world_width ? x + 1 : 0;
        const uint32_t top_y = y > 0 ? y - 1 : world_height - 1;
        const uint32_t bottom_y = y + 1 < world_height ? y + 1 : 0;
        const struct {
            uint32_t nx;
            uint32_t ny;
            float weight;
        } neighbors[] = {
            {left_x, y, 1.0f},       {right_x, y, 1.0f},       {x, top_y, 1.0f},
            {x, bottom_y, 1.0f},     {left_x, top_y, 0.65f},   {right_x, top_y, 0.65f},
            {left_x, bottom_y, 0.65f}, {right_x, bottom_y, 0.65f},
        };

        for (const auto& neighbor : neighbors) {
            const uint32_t neighbor_index = _worldDimension.indexOf(neighbor.nx, neighbor.ny);
            const uint32_t owner = imap[neighbor_index];
            float dir_x = 0.0f;
            float dir_y = 0.0f;
            neighbor_direction(neighbor.nx, neighbor.ny, dir_x, dir_y);
            add_owner_sample(boundary_samples, owner, dir_x, dir_y, neighbor.weight);
        }

        float total_opening = 0.0f;
        float owner_opening = 0.0f;
        float drift_x = 0.0f;
        float drift_y = 0.0f;
        float drift_weight = 0.0f;
        for (BoundaryOwnerSample& sample : boundary_samples) {
            const float normal_length =
                std::sqrt(sample.dir_x * sample.dir_x + sample.dir_y * sample.dir_y);
            if (normal_length <= FLT_EPSILON) {
                continue;
            }

            sample.dir_x /= normal_length;
            sample.dir_y /= normal_length;

            const Platec::FloatVector velocity = plates[sample.owner]->surfaceVelocityAt(x, y);
            sample.velocity_x = velocity.x();
            sample.velocity_y = velocity.y();
            sample.retreat =
                std::max(0.0f, -(sample.velocity_x * sample.dir_x + sample.velocity_y * sample.dir_y));
            total_opening += sample.retreat;
            if (sample.owner == imap[index]) {
                owner_opening = sample.retreat;
            }

            const float weight = 0.2f + sample.retreat;
            drift_x += sample.velocity_x * weight;
            drift_y += sample.velocity_y * weight;
            drift_weight += weight;
        }

        if (drift_weight > FLT_EPSILON) {
            drift_x /= drift_weight;
            drift_y /= drift_weight;
        }

        const float drift_magnitude = std::sqrt(drift_x * drift_x + drift_y * drift_y);
        const float owner_share =
            total_opening > FLT_EPSILON
                ? owner_opening / total_opening
                : (boundary_samples.empty() ? 1.0f : 1.0f / static_cast<float>(boundary_samples.size()));
        const float opening_strength = 1.0f - std::exp(-total_opening * 1.45f);
        const float random_signed = _randsource.next_float() * 2.0f - 1.0f;
        const float oceanic_ceiling = std::nextafter(CONTINENTAL_BASE, OCEANIC_BASE) - 0.01f;
        const float ridge_ratio =
            clamp_unit(0.18f + 0.18f * owner_share + 0.42f * opening_strength +
                       0.10f * std::min(1.5f, drift_magnitude) + 0.08f * random_signed);
        const float neighbor_height = average_neighbor_height(x, y);
        const bool using_neighbor_height = boundary_samples.empty() && neighbor_height > FLT_EPSILON;
        const float generated_crust = using_neighbor_height
            ? neighbor_height
            : OCEANIC_BASE + ridge_ratio * std::max(0.0f, oceanic_ceiling - OCEANIC_BASE);
        const uint8_t generated_material = using_neighbor_height &&
                generated_crust >= TopographyCodec::kContinentalBase
            ? platec::material::to_index(platec::material::Type::Granite)
            : platec::material::to_index(platec::material::Type::Basalt);

        hmap[index] = generated_crust;
        material_map[index] = generated_material;
        plates[imap[index]]->setCrust(x, y, generated_crust, iter_count, generated_material);
        ++plate_indices_found[imap[index]];
    };

    for (uint32_t index = 0; index < map_area; ++index) {
        if (imap[index] < num_plates) {
            if (hmap[index] > FLT_EPSILON) {
                ++plate_indices_found[imap[index]];
                continue;
            }

            imap[index] = 0xFFFFFFFFU;
            amap[index] = 0U;
            material_map[index] =
                platec::material::to_index(platec::material::kDefaultType);
        }

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        if (hasAssignedOwnerNeighbor(x, y)) {
            enqueue(index, frontier);
        }
    }

    if (frontier.empty()) {
        for (uint32_t index = 0; index < map_area; ++index) {
            enqueue(index, frontier);
        }
    }

    while (!frontier.empty()) {
        std::vector<uint32_t> next_frontier;
        next_frontier.reserve(frontier.size() * 2);

        while (!frontier.empty()) {
            const size_t pick = frontier.size() > 1 ? _randsource.next() % frontier.size() : 0;
            const uint32_t index = frontier[pick];
            frontier[pick] = frontier.back();
            frontier.pop_back();

            if (imap[index] < num_plates) {
                continue;
            }

            assign_new_crust(index);

            const uint32_t x = _worldDimension.xFromIndex(index);
            const uint32_t y = _worldDimension.yFromIndex(index);
            const uint32_t left_x = x > 0 ? x - 1 : _worldDimension.getWidth() - 1;
            const uint32_t right_x = x + 1 < _worldDimension.getWidth() ? x + 1 : 0;
            const uint32_t top_y = y > 0 ? y - 1 : _worldDimension.getHeight() - 1;
            const uint32_t bottom_y = y + 1 < _worldDimension.getHeight() ? y + 1 : 0;
            enqueue(_worldDimension.indexOf(left_x, y), next_frontier);
            enqueue(_worldDimension.indexOf(right_x, y), next_frontier);
            enqueue(_worldDimension.indexOf(x, top_y), next_frontier);
            enqueue(_worldDimension.indexOf(x, bottom_y), next_frontier);
            enqueue(_worldDimension.indexOf(left_x, top_y), next_frontier);
            enqueue(_worldDimension.indexOf(right_x, top_y), next_frontier);
            enqueue(_worldDimension.indexOf(left_x, bottom_y), next_frontier);
            enqueue(_worldDimension.indexOf(right_x, bottom_y), next_frontier);
        }

        frontier.swap(next_frontier);
    }

    for (uint32_t index = 0; index < map_area; ++index) {
        assign_new_crust(index);
    }
}

void lithosphere::applyPeriodicNoise() {
    if (num_plates == 0U || _steps <= 0) {
        return;
    }

    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    const uint32_t base_seed = static_cast<uint32_t>(scenario.seed);

    const auto apply_pass = [&](uint32_t seed_offset, uint32_t period, float strength,
                                float octaves, float persistence, float scale,
                                float noise_scale, float base_amplitude_m) {
        if (period == 0U || strength <= FLT_EPSILON ||
            static_cast<uint32_t>(_steps) % period != 0U) {
            return;
        }

        const uint32_t phase = static_cast<uint32_t>(_steps) / period;
        const FractalNoiseConfig config =
            make_noise_config(base_seed + seed_offset, octaves, persistence, scale, noise_scale);

        for (uint32_t y = 0; y < world_height; ++y) {
            for (uint32_t x = 0; x < world_width; ++x) {
                const uint32_t index = _worldDimension.indexOf(x, y);
                const uint32_t owner = imap[index];
                if (owner >= num_plates || hmap[index] <= 0.0f) {
                    continue;
                }

                const float sample = sample_wrapped_octave_noise(
                    static_cast<float>(x), static_cast<float>(y), _worldDimension,
                    config.octaves, config.persistence, config.scale, config.noise_scale,
                    -1.0f, 1.0f,
                    config.offset_a + static_cast<float>(phase) * 0.173f,
                    config.offset_b + static_cast<float>(phase) * 0.097f,
                    config.offset_c + static_cast<float>(phase) * 0.061f,
                    config.offset_d + static_cast<float>(phase) * 0.149f);
                const float crust_scale =
                    TopographyCodec::is_oceanic_internal(hmap[index]) ? 0.35f : 1.0f;
                const float delta_m = base_amplitude_m * strength * crust_scale * sample;
                const float delta_internal = metersToInternalDelta(std::fabs(delta_m));
                if (delta_internal <= FLT_EPSILON) {
                    continue;
                }

                hmap[index] = std::max(0.0f, hmap[index] +
                                                 (delta_m < 0.0f ? -delta_internal : delta_internal));
                plates[owner]->setCrust(x, y, hmap[index], amap[index]);
            }
        }
    };

    apply_pass(1U, scenario.hf_noise_period, scenario.hf_noise_strength, 5.0f, 0.56f, 2.4f,
               6.4f, 140.0f);
    apply_pass(2U, scenario.lf_noise_period, scenario.lf_noise_strength, 4.0f, 0.60f, 0.32f,
               1.2f, 360.0f);
}

void lithosphere::applySurfaceProcesses() {
    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U || num_plates == 0U) {
        return;
    }

    std::vector<float> base_height(map_area, 0.0f);
    std::vector<float> next_height(map_area, 0.0f);
    std::vector<float> deposited_height(map_area, 0.0f);
    for (uint32_t i = 0; i < map_area; ++i) {
        base_height[i] = hmap[i];
        next_height[i] = hmap[i];
    }

    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    const float bounded_erosion = std::max(0.0f, scenario.erosion_strength);
    const float cadence_scale = erosion_period > 0U
        ? std::min(1.0f, 6.0f / static_cast<float>(std::max<uint32_t>(1U, erosion_period)))
        : 1.0f;
    const float effective_erosion = bounded_erosion * cadence_scale;
    const float effective_glacial_erosion = effectiveGlacialErosionStrength();
    const float erosion_noise_scale = 0.55f + 0.45f * std::min(2.0f, effective_erosion);
    constexpr float kJunctionTerrainErosionBias = 0.15f;

    struct LowerNeighbor {
        uint32_t index = 0U;
        float drop_m = 0.0f;
    };
    const auto is_near_subduction_zone = [&](uint32_t index) {
        const auto boundary_type =
            static_cast<platec::contract::BoundaryType>(boundary_type_map[index]);
        const auto regime =
            static_cast<platec::contract::GeologicRegime>(geologic_regime_map[index]);
        return boundary_type == platec::contract::BoundaryType::Convergent ||
               regime == platec::contract::GeologicRegime::TrenchAdjacent ||
               (boundary_distance_map[index] <= 8U &&
                convergence_map[index] > divergence_map[index]);
    };

    for (uint32_t y = 0; y < world_height; ++y) {
        const uint32_t top_y = y > 0U ? y - 1U : world_height - 1U;
        const uint32_t bottom_y = y + 1U < world_height ? y + 1U : 0U;

        for (uint32_t x = 0; x < world_width; ++x) {
            const uint32_t left_x = x > 0U ? x - 1U : world_width - 1U;
            const uint32_t right_x = x + 1U < world_width ? x + 1U : 0U;
            const uint32_t index = _worldDimension.indexOf(x, y);
            const uint32_t owner = imap[index];

            if (owner >= num_plates || base_height[index] <= 0.0f) {
                continue;
            }

            const float center_m = internalToPhysicalMeters(base_height[index]);
            const float elevation_m = center_m - static_cast<float>(sea_level_m);
            if (elevation_m <= 0.0f) {
                if (crust_class_map[index] ==
                        static_cast<uint8_t>(platec::contract::CrustClass::Oceanic) &&
                    is_near_subduction_zone(index)) {
                    const float distance = static_cast<float>(boundary_distance_map[index]);
                    const float approach = std::exp(-distance / 2.4f);
                    const float loss_fraction =
                        clamp_unit(subduction_strength * (0.035f + 0.30f * approach));
                    next_height[index] =
                        std::max(0.0f, next_height[index] * (1.0f - loss_fraction));
                }
                continue;
            }

            const plate* owning_plate = plates[owner];
            const Platec::FloatVector velocity = owning_plate->surfaceVelocityAt(x, y);
            const float speed =
                std::sqrt(velocity.x() * velocity.x() + velocity.y() * velocity.y());
            const float slow_factor = clamp_unit(1.0f / (1.0f + 1.35f * speed));

            const LowerNeighbor neighbors[] = {
                {_worldDimension.indexOf(left_x, y), 0.0f},
                {_worldDimension.indexOf(right_x, y), 0.0f},
                {_worldDimension.indexOf(x, top_y), 0.0f},
                {_worldDimension.indexOf(x, bottom_y), 0.0f},
                {_worldDimension.indexOf(left_x, top_y), 0.0f},
                {_worldDimension.indexOf(right_x, top_y), 0.0f},
                {_worldDimension.indexOf(left_x, bottom_y), 0.0f},
                {_worldDimension.indexOf(right_x, bottom_y), 0.0f},
            };

            LowerNeighbor lower_neighbors[8];
            uint32_t lower_neighbor_count = 0U;
            float total_drop_m = 0.0f;
            float max_drop_m = 0.0f;

            for (const LowerNeighbor& neighbor : neighbors) {
                const float neighbor_m = internalToPhysicalMeters(base_height[neighbor.index]);
                const float drop_m = center_m - neighbor_m;
                if (drop_m <= 0.0f) {
                    continue;
                }

                lower_neighbors[lower_neighbor_count].index = neighbor.index;
                lower_neighbors[lower_neighbor_count].drop_m = drop_m;
                ++lower_neighbor_count;
                total_drop_m += drop_m;
                max_drop_m = std::max(max_drop_m, drop_m);
            }

            if (lower_neighbor_count == 0U) {
                continue;
            }

            const float mean_drop_m = total_drop_m / static_cast<float>(lower_neighbor_count);
            const float slope_norm =
                clamp_unit((0.65f * max_drop_m + 0.35f * mean_drop_m) / 1800.0f);
            const float support_noise =
                0.85f +
                0.30f *
                    (0.5f +
                     0.5f * sample_wrapped_octave_noise(
                                 static_cast<float>(x), static_cast<float>(y), _worldDimension,
                                 3.0f, 0.58f, 0.16f, 1.9f, -1.0f, 1.0f,
                                 static_cast<float>(iter_count) * 0.11f + 13.0f, 17.0f, 29.0f,
                                 43.0f));
            const float glacial_noise =
                0.80f +
                0.40f *
                    (0.5f +
                     0.5f * sample_wrapped_octave_noise(
                                 static_cast<float>(x), static_cast<float>(y), _worldDimension,
                                 3.0f, 0.62f, 0.18f, 2.2f, -1.0f, 1.0f, 23.0f,
                                 static_cast<float>(iter_count) * 0.09f + 31.0f, 37.0f, 47.0f));
            const float denudation_noise =
                0.75f +
                0.50f *
                    (0.5f +
                     0.5f * sample_wrapped_octave_noise(
                                 static_cast<float>(x), static_cast<float>(y), _worldDimension,
                                 4.0f, 0.56f, 0.14f, 1.7f, -1.0f, 1.0f, 41.0f, 53.0f,
                                 static_cast<float>(iter_count) * 0.07f + 59.0f, 67.0f));
            const float landslide_roll =
                0.5f +
                0.5f * sample_wrapped_octave_noise(
                            static_cast<float>(x), static_cast<float>(y), _worldDimension, 4.0f,
                            0.60f, 0.26f, 2.8f, -1.0f, 1.0f,
                            static_cast<float>(iter_count) * 0.13f + 71.0f, 79.0f, 83.0f, 97.0f);

            const float height_scale_m = std::max(750.0f, effectiveHeightLimitMeters(index));
            const float height_ratio = elevation_m / height_scale_m;
            const float overshoot_m = std::max(0.0f, elevation_m - height_scale_m);
            const float support_loss_m =
                support_noise *
                (height_scale_m * 0.0120f * smoothstep(0.40f, 0.90f, height_ratio) *
                     (0.45f + 0.55f * slope_norm) +
                 overshoot_m * (0.55f + 0.60f * slope_norm));

            const float snow_line_noise =
                0.88f +
                0.24f *
                    (0.5f +
                     0.5f * sample_wrapped_octave_noise(
                                 static_cast<float>(x), static_cast<float>(y), _worldDimension,
                                 2.0f, 0.54f, 0.10f, 1.2f, -1.0f, 1.0f, 101.0f, 107.0f,
                                 static_cast<float>(iter_count) * 0.05f + 109.0f, 113.0f));
            const float local_snow_line_m =
                snow_line_meters(y, world_height) * snow_line_noise;
            const float glacial_excess_m = std::max(0.0f, elevation_m - local_snow_line_m);
            const float glacial_ratio = clamp_unit(glacial_excess_m / 700.0f);
            const float strong_glacier_scale =
                1.0f + 0.45f * std::max(0.0f, effective_glacial_erosion - 1.0f);
            const float glacial_loss_m = std::min(
                elevation_m,
                effective_glacial_erosion * strong_glacier_scale * glacial_noise * glacial_ratio *
                    (1.5f + 3.5f * slope_norm +
                     0.0025f * std::min(2500.0f, glacial_excess_m)));

            const float denudation_loss_m =
                effective_erosion * erosion_noise_scale * denudation_noise * slow_factor *
                (55.0f + 280.0f * slope_norm + 140.0f * std::max(0.0f, height_ratio - 0.15f));
            float subduction_approach_loss_m = 0.0f;
            const float junction_subduction =
                junction_subduction_influence_map[index] * kJunctionTerrainErosionBias;
            if (crust_class_map[index] ==
                    static_cast<uint8_t>(platec::contract::CrustClass::Oceanic) &&
                (is_near_subduction_zone(index) || junction_subduction > 0.01f)) {
                const float distance = static_cast<float>(boundary_distance_map[index]);
                const float approach = std::max(std::exp(-distance / 2.4f), junction_subduction);
                subduction_approach_loss_m =
                    center_m * clamp_unit(subduction_strength * (0.025f + 0.22f * approach));
            }
            float transform_fault_loss_m = 0.0f;
            const float direct_junction_transform =
                junction_transform_influence_map[index] * kJunctionTerrainErosionBias;
            if (divergent_carve_strength > FLT_EPSILON &&
                (shear_map[index] > 0U || divergence_map[index] > 0U ||
                 direct_junction_transform > 0.01f)) {
                const int max_radius = 6;
                float transform_influence = direct_junction_transform;
                for (int dy = -max_radius; dy <= max_radius; ++dy) {
                    for (int dx = -max_radius; dx <= max_radius; ++dx) {
                        const float distance =
                            std::sqrt(static_cast<float>(dx * dx + dy * dy));
                        if (distance > static_cast<float>(max_radius)) {
                            continue;
                        }

                        const int wrapped_x =
                            (static_cast<int>(x) + dx + static_cast<int>(world_width)) %
                            static_cast<int>(world_width);
                        const int wrapped_y =
                            (static_cast<int>(y) + dy + static_cast<int>(world_height)) %
                            static_cast<int>(world_height);
                        const uint32_t nx = static_cast<uint32_t>(wrapped_x);
                        const uint32_t ny = static_cast<uint32_t>(wrapped_y);
                        const uint32_t neighbor_index = _worldDimension.indexOf(nx, ny);
                        if (static_cast<platec::contract::BoundaryType>(
                                boundary_type_map[neighbor_index]) !=
                            platec::contract::BoundaryType::Transform) {
                            continue;
                        }

                        const float strain = clamp_unit(accumulated_strain_map[neighbor_index]);
                        const float shear = clamp_unit(
                            static_cast<float>(shear_map[neighbor_index]) / 100.0f);
                        const float divergence = clamp_unit(
                            static_cast<float>(divergence_map[neighbor_index]) / 100.0f);
                        const float extent = 1.5f + 6.5f * strain;
                        if (distance > extent) {
                            continue;
                        }

                        const float falloff = 1.0f - distance / (extent + 0.5f);
                        const float local_influence =
                            falloff * (0.25f + 0.75f * shear) * (0.20f + 0.80f * divergence) *
                            (0.15f + 0.85f * strain);
                        transform_influence = std::max(transform_influence, local_influence);
                    }
                }

                transform_fault_loss_m =
                    divergent_carve_strength * (25.0f + 260.0f * transform_influence);
            }

            float direct_loss_internal =
                metersToInternalDelta(support_loss_m + glacial_loss_m + denudation_loss_m +
                                      subduction_approach_loss_m + transform_fault_loss_m);
            direct_loss_internal = std::min(direct_loss_internal, next_height[index]);
            next_height[index] -= direct_loss_internal;

            const float landslide_strength =
                clamp_unit((0.08f + 0.90f * slope_norm * slope_norm) *
                           (0.45f + 0.55f * slow_factor) *
                           (0.45f + 0.25f * std::min(2.0f, effective_erosion)));
            if (landslide_roll >= landslide_strength) {
                continue;
            }

            float landslide_loss_m =
                (60.0f + 420.0f * slope_norm + 0.050f * glacial_excess_m) *
                (0.65f + 0.35f * slow_factor) * (0.80f + 0.35f * erosion_noise_scale);
            if (height_ratio > 1.0f) {
                landslide_loss_m += (height_ratio - 1.0f) * 120.0f;
            }

            float landslide_internal =
                std::min(metersToInternalDelta(landslide_loss_m), next_height[index]);
            if (landslide_internal <= FLT_EPSILON) {
                continue;
            }

            next_height[index] -= landslide_internal;

            float lower_weight_sum = 0.0f;
            for (uint32_t n = 0; n < lower_neighbor_count; ++n) {
                lower_weight_sum += std::pow(lower_neighbors[n].drop_m, 1.35f);
            }
            if (lower_weight_sum <= FLT_EPSILON) {
                next_height[index] += landslide_internal;
                continue;
            }

            for (uint32_t n = 0; n < lower_neighbor_count; ++n) {
                const float weight = std::pow(lower_neighbors[n].drop_m, 1.35f) / lower_weight_sum;
                deposited_height[lower_neighbors[n].index] += landslide_internal * weight;
            }
        }
    }

    std::fill(plate_indices_found.begin(), plate_indices_found.end(), 0U);
    for (uint32_t y = 0; y < world_height; ++y) {
        for (uint32_t x = 0; x < world_width; ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const uint32_t owner = imap[index];

            hmap[index] = std::max(0.0f, next_height[index] + deposited_height[index]);
            hmap[index] = softenHeight(index, hmap[index]);
            if (owner >= num_plates) {
                continue;
            }

            plates[owner]->setCrust(x, y, hmap[index], amap[index]);
            if (hmap[index] > 2.0f * FLT_EPSILON) {
                material_map[index] = plates[owner]->getMaterial(x, y);
                ++plate_indices_found[owner];
            } else {
                imap[index] = 0xFFFFFFFFU;
                amap[index] = 0U;
                material_map[index] =
                    platec::material::to_index(platec::material::kDefaultType);
            }
        }
    }
}

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

bool lithosphere::triggerPlateBirth() {
    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U || num_plates == 0U || rift_tension_x_map.size() != map_area ||
        rift_tension_y_map.size() != map_area) {
        return false;
    }

    const uint32_t world_width = _worldDimension.getWidth();
    const uint32_t world_height = _worldDimension.getHeight();
    const auto tension_magnitude = [&](uint32_t index) {
        const float tx = rift_tension_x_map[index];
        const float ty = rift_tension_y_map[index];
        return std::sqrt(tx * tx + ty * ty);
    };
    const auto sample_opening_vector = [&](uint32_t x, uint32_t y, uint32_t owner) {
        Platec::FloatVector opening(0.0f, 0.0f);
        if (owner >= num_plates) {
            return opening;
        }

        const uint32_t left_x = x > 0U ? x - 1U : world_width - 1U;
        const uint32_t right_x = x + 1U < world_width ? x + 1U : 0U;
        const uint32_t top_y = y > 0U ? y - 1U : world_height - 1U;
        const uint32_t bottom_y = y + 1U < world_height ? y + 1U : 0U;
        const uint32_t left_index = _worldDimension.indexOf(left_x, y);
        const uint32_t right_index = _worldDimension.indexOf(right_x, y);
        const uint32_t top_index = _worldDimension.indexOf(x, top_y);
        const uint32_t bottom_index = _worldDimension.indexOf(x, bottom_y);

        if (imap[left_index] == owner && imap[right_index] == owner) {
            const Platec::FloatVector left_velocity = plates[owner]->surfaceVelocityAt(left_x, y);
            const Platec::FloatVector right_velocity = plates[owner]->surfaceVelocityAt(right_x, y);
            const float horizontal_opening = right_velocity.x() - left_velocity.x();
            if (horizontal_opening > FLT_EPSILON) {
                opening = Platec::FloatVector(opening.x() + horizontal_opening, opening.y());
            }
        }

        if (imap[top_index] == owner && imap[bottom_index] == owner) {
            const Platec::FloatVector top_velocity = plates[owner]->surfaceVelocityAt(x, top_y);
            const Platec::FloatVector bottom_velocity =
                plates[owner]->surfaceVelocityAt(x, bottom_y);
            const float vertical_opening = bottom_velocity.y() - top_velocity.y();
            if (vertical_opening > FLT_EPSILON) {
                opening = Platec::FloatVector(opening.x(), opening.y() + vertical_opening);
            }
        }

        return opening;
    };

    uint32_t seed_index = std::numeric_limits<uint32_t>::max();
    uint32_t seed_owner = std::numeric_limits<uint32_t>::max();
    float best_tension = kPlateBirthTensionThreshold;

    for (uint32_t y = 0; y < world_height; ++y) {
        for (uint32_t x = 0; x < world_width; ++x) {
            const uint32_t index = _worldDimension.indexOf(x, y);
            const uint32_t owner = imap[index];
            const float junction_rift = junction_rift_influence_map[index];

            Platec::FloatVector next_tension(rift_tension_x_map[index] * kPlateBirthTensionDecay,
                                             rift_tension_y_map[index] * kPlateBirthTensionDecay);
            bool eligible = false;
            if (owner < num_plates && hmap[index] > FLT_EPSILON && boundary_distance_map[index] > 1U) {
                const Platec::FloatVector opening = sample_opening_vector(x, y, owner);
                if (opening.length() > FLT_EPSILON) {
                    next_tension = Platec::FloatVector(
                        next_tension.x() + opening.x() * kPlateBirthTensionGain,
                        next_tension.y() + opening.y() * kPlateBirthTensionGain);
                    eligible = true;
                }
                if (junction_rift > 0.15f && next_tension.length() >= kPlateBirthPropagationThreshold) {
                    eligible = true;
                }
            }

            if (!eligible) {
                next_tension = Platec::FloatVector(next_tension.x() * 0.6f, next_tension.y() * 0.6f);
            }
            const float next_length = next_tension.length();
            if (next_length > kPlateBirthTensionThreshold * 2.5f) {
                const float scale = (kPlateBirthTensionThreshold * 2.5f) / next_length;
                next_tension = Platec::FloatVector(next_tension.x() * scale, next_tension.y() * scale);
            }
            rift_tension_x_map[index] = next_tension.x();
            rift_tension_y_map[index] = next_tension.y();

            if (eligible && next_length >= best_tension) {
                best_tension = next_length;
                seed_index = index;
                seed_owner = owner;
            }
        }
    }

    updateTensionMagnitudeMap();

    if (seed_index == std::numeric_limits<uint32_t>::max() || seed_owner >= num_plates) {
        return false;
    }

    std::vector<uint8_t> crack_mask(map_area, 0U);
    std::vector<uint32_t> crack_cells;
    std::vector<uint32_t> crack_frontier = {seed_index};
    std::vector<uint32_t> boundary_contact_cells;
    Platec::FloatVector accumulated_opening(0.0f, 0.0f);
    float opening_weight = 0.0f;

    while (!crack_frontier.empty()) {
        const uint32_t index = crack_frontier.back();
        crack_frontier.pop_back();
        if (crack_mask[index] != 0U || imap[index] != seed_owner || hmap[index] <= FLT_EPSILON) {
            continue;
        }
        if (index != seed_index && boundary_distance_map[index] > 1U &&
            tension_magnitude(index) < kPlateBirthPropagationThreshold) {
            continue;
        }

        crack_mask[index] = 1U;
        crack_cells.push_back(index);

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        const Platec::FloatVector opening_vector = sample_opening_vector(x, y, seed_owner);
        const float opening_strength = opening_vector.length();
        if (opening_strength > FLT_EPSILON) {
            accumulated_opening =
                Platec::FloatVector(accumulated_opening.x() + opening_vector.x(),
                                    accumulated_opening.y() + opening_vector.y());
            opening_weight += opening_strength;
        }
        const Platec::FloatVector tension_vector(rift_tension_x_map[index], rift_tension_y_map[index]);
        const float tension_strength = tension_vector.length();
        if (tension_strength > FLT_EPSILON) {
            accumulated_opening =
                Platec::FloatVector(accumulated_opening.x() + tension_vector.x(),
                                    accumulated_opening.y() + tension_vector.y());
            opening_weight += tension_strength;
        }
        if (boundary_distance_map[index] <= 1U) {
            boundary_contact_cells.push_back(index);
        }

        const int neighbor_dx[] = {-1, 1, 0, 0, -1, 1, -1, 1};
        const int neighbor_dy[] = {0, 0, -1, 1, -1, -1, 1, 1};
        for (size_t neighbor = 0; neighbor < 8U; ++neighbor) {
            const uint32_t nx = wrap_index(static_cast<int>(x) + neighbor_dx[neighbor], world_width);
            const uint32_t ny = wrap_index(static_cast<int>(y) + neighbor_dy[neighbor], world_height);
            const uint32_t neighbor_index = _worldDimension.indexOf(nx, ny);
            if (crack_mask[neighbor_index] != 0U || imap[neighbor_index] != seed_owner) {
                continue;
            }
            if (boundary_distance_map[neighbor_index] <= 1U ||
                tension_magnitude(neighbor_index) >= kPlateBirthPropagationThreshold) {
                crack_frontier.push_back(neighbor_index);
            }
        }
    }

    uint32_t max_contact_span = 0U;
    for (uint32_t a = 0; a < boundary_contact_cells.size(); ++a) {
        const uint32_t ax = _worldDimension.xFromIndex(boundary_contact_cells[a]);
        const uint32_t ay = _worldDimension.yFromIndex(boundary_contact_cells[a]);
        for (uint32_t b = a + 1U; b < boundary_contact_cells.size(); ++b) {
            const uint32_t bx = _worldDimension.xFromIndex(boundary_contact_cells[b]);
            const uint32_t by = _worldDimension.yFromIndex(boundary_contact_cells[b]);
            const uint32_t raw_dx = ax > bx ? ax - bx : bx - ax;
            const uint32_t raw_dy = ay > by ? ay - by : by - ay;
            const uint32_t dx = std::min(raw_dx, world_width - raw_dx);
            const uint32_t dy = std::min(raw_dy, world_height - raw_dy);
            max_contact_span = std::max(max_contact_span, dx + dy);
        }
    }
    const uint32_t min_contact_span = std::max<uint32_t>(
        3U, std::min(world_width, world_height) / 8U);

    if (crack_cells.size() < 4U || boundary_contact_cells.size() < 2U ||
        max_contact_span < min_contact_span) {
        return false;
    }

    const uint32_t min_split_area = std::max<uint32_t>(24U, map_area / 6000U);
    std::vector<uint8_t> component_visited(map_area, 0U);
    std::vector<uint8_t> chosen_component(map_area, 0U);
    uint32_t chosen_component_size = std::numeric_limits<uint32_t>::max();
    uint32_t owner_area = 0U;
    uint32_t component_count = 0U;

    for (uint32_t index = 0; index < map_area; ++index) {
        if (imap[index] != seed_owner || hmap[index] <= FLT_EPSILON) {
            continue;
        }
        if (crack_mask[index] != 0U) {
            ++owner_area;
            continue;
        }

        ++owner_area;
        if (component_visited[index] != 0U) {
            continue;
        }

        std::vector<uint32_t> frontier = {index};
        std::vector<uint32_t> component_cells;
        component_visited[index] = 1U;
        while (!frontier.empty()) {
            const uint32_t current = frontier.back();
            frontier.pop_back();
            component_cells.push_back(current);

            const uint32_t x = _worldDimension.xFromIndex(current);
            const uint32_t y = _worldDimension.yFromIndex(current);
            const uint32_t neighbors[] = {
                _worldDimension.indexOf(x > 0U ? x - 1U : world_width - 1U, y),
                _worldDimension.indexOf(x + 1U < world_width ? x + 1U : 0U, y),
                _worldDimension.indexOf(x, y > 0U ? y - 1U : world_height - 1U),
                _worldDimension.indexOf(x, y + 1U < world_height ? y + 1U : 0U),
            };

            for (uint32_t neighbor_index : neighbors) {
                if (component_visited[neighbor_index] != 0U || crack_mask[neighbor_index] != 0U ||
                    imap[neighbor_index] != seed_owner || hmap[neighbor_index] <= FLT_EPSILON) {
                    continue;
                }
                component_visited[neighbor_index] = 1U;
                frontier.push_back(neighbor_index);
            }
        }

        ++component_count;
        if (component_cells.size() < chosen_component_size) {
            std::fill(chosen_component.begin(), chosen_component.end(), static_cast<uint8_t>(0));
            for (uint32_t component_index : component_cells) {
                chosen_component[component_index] = 1U;
            }
            chosen_component_size = static_cast<uint32_t>(component_cells.size());
        }
    }

    if (component_count < 2U || chosen_component_size == std::numeric_limits<uint32_t>::max()) {
        return false;
    }

    std::vector<uint8_t> extract_mask = chosen_component;
    for (uint32_t crack_index : crack_cells) {
        extract_mask[crack_index] = 1U;
    }

    uint32_t extract_area = 0U;
    for (uint8_t selected : extract_mask) {
        extract_area += selected != 0U;
    }
    const uint32_t remaining_area = owner_area > extract_area ? owner_area - extract_area : 0U;
    if (extract_area < min_split_area || remaining_area < min_split_area) {
        return false;
    }

    uint32_t anchor_index = std::numeric_limits<uint32_t>::max();
    for (uint32_t index = 0; index < map_area; ++index) {
        if (extract_mask[index] != 0U) {
            anchor_index = index;
            break;
        }
    }
    if (anchor_index == std::numeric_limits<uint32_t>::max()) {
        return false;
    }

    std::vector<uint8_t> extract_visited(map_area, 0U);
    std::vector<int32_t> unwrapped_x(map_area, 0);
    std::vector<int32_t> unwrapped_y(map_area, 0);
    std::deque<uint32_t> extract_frontier;
    extract_frontier.push_back(anchor_index);
    extract_visited[anchor_index] = 1U;
    const int32_t anchor_x = static_cast<int32_t>(_worldDimension.xFromIndex(anchor_index));
    const int32_t anchor_y = static_cast<int32_t>(_worldDimension.yFromIndex(anchor_index));
    unwrapped_x[anchor_index] = anchor_x;
    unwrapped_y[anchor_index] = anchor_y;

    int32_t min_x = anchor_x;
    int32_t max_x = anchor_x;
    int32_t min_y = anchor_y;
    int32_t max_y = anchor_y;
    uint32_t visited_extract_cells = 0U;

    while (!extract_frontier.empty()) {
        const uint32_t current = extract_frontier.front();
        extract_frontier.pop_front();
        ++visited_extract_cells;

        min_x = std::min(min_x, unwrapped_x[current]);
        max_x = std::max(max_x, unwrapped_x[current]);
        min_y = std::min(min_y, unwrapped_y[current]);
        max_y = std::max(max_y, unwrapped_y[current]);

        const uint32_t x = _worldDimension.xFromIndex(current);
        const uint32_t y = _worldDimension.yFromIndex(current);
        const int neighbor_dx[] = {-1, 1, 0, 0};
        const int neighbor_dy[] = {0, 0, -1, 1};
        for (size_t neighbor = 0; neighbor < 4U; ++neighbor) {
            const uint32_t nx = wrap_index(static_cast<int>(x) + neighbor_dx[neighbor], world_width);
            const uint32_t ny = wrap_index(static_cast<int>(y) + neighbor_dy[neighbor], world_height);
            const uint32_t neighbor_index = _worldDimension.indexOf(nx, ny);
            if (extract_mask[neighbor_index] == 0U || extract_visited[neighbor_index] != 0U) {
                continue;
            }

            extract_visited[neighbor_index] = 1U;
            unwrapped_x[neighbor_index] = unwrapped_x[current] + neighbor_dx[neighbor];
            unwrapped_y[neighbor_index] = unwrapped_y[current] + neighbor_dy[neighbor];
            extract_frontier.push_back(neighbor_index);
        }
    }

    if (visited_extract_cells != extract_area) {
        return false;
    }

    const uint32_t plate_width = static_cast<uint32_t>(max_x - min_x + 1);
    const uint32_t plate_height = static_cast<uint32_t>(max_y - min_y + 1);
    if (plate_width == 0U || plate_height == 0U || plate_width > world_width ||
        plate_height > world_height) {
        return false;
    }

    std::vector<float> plate_map(static_cast<size_t>(plate_width) * static_cast<size_t>(plate_height),
                                 0.0f);
    std::vector<uint8_t> plate_material(
        static_cast<size_t>(plate_width) * static_cast<size_t>(plate_height),
        platec::material::to_index(platec::material::kDefaultType));

    for (uint32_t index = 0; index < map_area; ++index) {
        if (extract_mask[index] == 0U) {
            continue;
        }

        const uint32_t local_x = static_cast<uint32_t>(unwrapped_x[index] - min_x);
        const uint32_t local_y = static_cast<uint32_t>(unwrapped_y[index] - min_y);
        const size_t local_index =
            static_cast<size_t>(local_y) * static_cast<size_t>(plate_width) + local_x;
        plate_map[local_index] = hmap[index];
        plate_material[local_index] = material_map[index];
    }

    ensurePlateCapacity(num_plates + 1U);
    const uint32_t new_plate_index = num_plates;
    const uint32_t plate_left = wrap_index(min_x, world_width);
    const uint32_t plate_top = wrap_index(min_y, world_height);
    float* local_height = new float[plate_map.size()];
    std::copy(plate_map.begin(), plate_map.end(), local_height);

    plate* newborn = new plate(_randsource.next(), local_height, plate_width, plate_height, plate_left,
                               plate_top, iter_count, _worldDimension, plate_material.data(),
                               erosion_strength, crust_rotation_strength, rotation_strength,
                               scenario.movement_energy);

    Platec::FloatVector kick_direction = accumulated_opening;
    if (opening_weight > FLT_EPSILON) {
        kick_direction =
            Platec::FloatVector(kick_direction.x() / opening_weight, kick_direction.y() / opening_weight);
    }
    const float kick_length = kick_direction.length();
    if (kick_length > FLT_EPSILON) {
        kick_direction =
            Platec::FloatVector(kick_direction.x() / kick_length, kick_direction.y() / kick_length);
    } else {
        kick_direction = Platec::FloatVector(1.0f, 0.0f);
    }
    newborn->setKinematics(
        Platec::FloatVector(plates[seed_owner]->linearVelocityVector().x() +
                                kick_direction.x() * kPlateBirthKickStrength,
                            plates[seed_owner]->linearVelocityVector().y() +
                                kick_direction.y() * kPlateBirthKickStrength),
        plates[seed_owner]->getAngularVelocity());

    for (uint32_t index = 0; index < map_area; ++index) {
        if (extract_mask[index] == 0U) {
            continue;
        }

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        newborn->setCrustTimestamp(x, y, amap[index]);
        plates[seed_owner]->setCrust(x, y, 0.0f, amap[index]);
        imap[index] = new_plate_index;
        rift_tension_x_map[index] = 0.0f;
        rift_tension_y_map[index] = 0.0f;
    }

    plates[new_plate_index] = newborn;
    collisions[new_plate_index].clear();
    subductions[new_plate_index].clear();
    ++num_plates;
    updateTensionMagnitudeMap();
    return true;
}

void lithosphere::enforceSinglePlateRegions() {
    const uint32_t map_area = _worldDimension.getArea();
    if (map_area == 0U || num_plates == 0U) {
        return;
    }

    const uint32_t width = _worldDimension.getWidth();
    const uint32_t height = _worldDimension.getHeight();
    std::vector<uint8_t> visited(map_area, 0U);
    std::vector<std::vector<uint32_t>> components;
    std::vector<uint32_t> component_owner;
    components.reserve(num_plates);
    component_owner.reserve(num_plates);

    for (uint32_t index = 0; index < map_area; ++index) {
        const uint32_t owner = imap[index];
        if (visited[index] != 0U || owner >= num_plates || hmap[index] <= FLT_EPSILON) {
            continue;
        }

        std::vector<uint32_t> cells;
        std::deque<uint32_t> frontier;
        frontier.push_back(index);
        visited[index] = 1U;

        while (!frontier.empty()) {
            const uint32_t current = frontier.front();
            frontier.pop_front();
            cells.push_back(current);

            const uint32_t x = _worldDimension.xFromIndex(current);
            const uint32_t y = _worldDimension.yFromIndex(current);
            const uint32_t neighbors[] = {
                _worldDimension.indexOf(x > 0U ? x - 1U : width - 1U, y),
                _worldDimension.indexOf(x + 1U < width ? x + 1U : 0U, y),
                _worldDimension.indexOf(x, y > 0U ? y - 1U : height - 1U),
                _worldDimension.indexOf(x, y + 1U < height ? y + 1U : 0U),
            };

            for (uint32_t neighbor : neighbors) {
                if (visited[neighbor] != 0U || imap[neighbor] != owner ||
                    hmap[neighbor] <= FLT_EPSILON) {
                    continue;
                }
                visited[neighbor] = 1U;
                frontier.push_back(neighbor);
            }
        }

        component_owner.push_back(owner);
        components.push_back(std::move(cells));
    }

    std::vector<uint32_t> largest_component(num_plates, std::numeric_limits<uint32_t>::max());
    std::vector<uint32_t> largest_size(num_plates, 0U);
    for (uint32_t component_index = 0; component_index < components.size(); ++component_index) {
        const uint32_t owner = component_owner[component_index];
        const uint32_t size = static_cast<uint32_t>(components[component_index].size());
        if (size > largest_size[owner]) {
            largest_size[owner] = size;
            largest_component[owner] = component_index;
        }
    }

    bool changed = false;
    for (uint32_t component_index = 0; component_index < components.size(); ++component_index) {
        const uint32_t source_owner = component_owner[component_index];
        if (largest_component[source_owner] == component_index) {
            continue;
        }

        std::vector<uint32_t> contact_counts(num_plates, 0U);
        for (uint32_t cell : components[component_index]) {
            const uint32_t x = _worldDimension.xFromIndex(cell);
            const uint32_t y = _worldDimension.yFromIndex(cell);
            const uint32_t neighbors[] = {
                _worldDimension.indexOf(x > 0U ? x - 1U : width - 1U, y),
                _worldDimension.indexOf(x + 1U < width ? x + 1U : 0U, y),
                _worldDimension.indexOf(x, y > 0U ? y - 1U : height - 1U),
                _worldDimension.indexOf(x, y + 1U < height ? y + 1U : 0U),
            };

            for (uint32_t neighbor : neighbors) {
                const uint32_t neighbor_owner = imap[neighbor];
                if (neighbor_owner < num_plates && neighbor_owner != source_owner &&
                    hmap[neighbor] > FLT_EPSILON) {
                    ++contact_counts[neighbor_owner];
                }
            }
        }

        uint32_t target_owner = std::numeric_limits<uint32_t>::max();
        uint32_t best_contacts = 0U;
        for (uint32_t owner = 0; owner < num_plates; ++owner) {
            if (contact_counts[owner] > best_contacts) {
                best_contacts = contact_counts[owner];
                target_owner = owner;
            }
        }
        if (target_owner >= num_plates || best_contacts == 0U) {
            continue;
        }

        for (uint32_t cell : components[component_index]) {
            const uint32_t x = _worldDimension.xFromIndex(cell);
            const uint32_t y = _worldDimension.yFromIndex(cell);
            const float source_height = hmap[cell];
            const uint32_t source_age = amap[cell];
            const uint8_t source_material = material_map[cell];
            const float target_height = plates[target_owner]->getCrust(x, y);
            const uint8_t target_material = plates[target_owner]->getMaterial(x, y);
            const bool keep_target_material =
                target_height > FLT_EPSILON && is_continental_material(target_material) &&
                !is_continental_material(source_material);

            plates[source_owner]->setCrust(x, y, 0.0f, source_age);
            plates[target_owner]->setCrust(
                x, y, std::max(source_height, target_height), source_age,
                keep_target_material ? target_material : source_material);

            imap[cell] = target_owner;
            hmap[cell] = std::max(source_height, target_height);
            amap[cell] = source_age;
            material_map[cell] = keep_target_material ? target_material : source_material;
            rift_tension_x_map[cell] = 0.0f;
            rift_tension_y_map[cell] = 0.0f;
            changed = true;
        }
    }

    std::fill(plate_indices_found.begin(), plate_indices_found.end(), 0U);
    for (uint32_t index = 0; index < map_area; ++index) {
        const uint32_t owner = imap[index];
        if (owner < num_plates && hmap[index] > FLT_EPSILON) {
            ++plate_indices_found[owner];
        }
    }

    if (changed) {
        updateTensionMagnitudeMap();
    }
}

// Remove empty plates from the system.
void lithosphere::removeEmptyPlates() {
    for (uint32_t i = 0; i < num_plates; ++i) {
        if (num_plates == 1)
            puts("ONLY ONE PLATE LEFT!");
        else if (plate_indices_found[i] == 0) {
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

void lithosphere::update() {
    try {
        _steps++;
        float totalVelocity = 0;
        float systemKineticEnergy = 0;

        for (uint32_t i = 0; i < num_plates; ++i) {
            totalVelocity += plates[i]->getVelocity();
            systemKineticEnergy += plates[i]->getMomentum();
        }

        if (systemKineticEnergy > peak_Ek) {
            peak_Ek = systemKineticEnergy;
        }

        // If there's no continental collisions during past iterations,
        // then interesting activity has ceased and we should restart.
        // Also if the simulation has been going on for too long already,
        // restart, because interesting stuff has most likely ended.
        const bool fixed_cycle_duration = scenario.cycle_duration_myr > 0.0;
        const bool duration_elapsed =
            cycle_step_limit > 0 && _steps > static_cast<int>(cycle_step_limit);
        const bool early_restart =
            totalVelocity < RESTART_SPEED_LIMIT ||
            systemKineticEnergy / peak_Ek < RESTART_ENERGY_RATIO ||
            last_coll_count > NO_COLLISION_TIME_LIMIT;
        if ((fixed_cycle_duration && duration_elapsed) ||
            (!fixed_cycle_duration && (early_restart || duration_elapsed))) {
            const bool completing_last_cycle =
                max_cycles > 0U && cycle_count + 1U >= max_cycles;
            if (completing_last_cycle) {
                finishCurrentState();
                return;
            }
            restart();
            return;
        }

        const uint32_t map_area = _worldDimension.getArea();
        // Keep a copy of the previous index map
        prev_hmap.copy(hmap);
        prev_imap.copy(imap);
        prev_display_hmap.copy(display_hmap);

        // Realize accumulated external forces to each plate.
        for (uint32_t i = 0; i < num_plates; ++i) {
            plates[i]->resetSegments();

            if (erosion_period > 0 && _steps % erosion_period == 0)
                plates[i]->erode(CONTINENTAL_BASE, iter_count);

            plates[i]->move();
        }

        uint32_t oceanic_collisions = 0;
        uint32_t continental_collisions = 0;

        updateHeightAndPlateIndexMaps(oceanic_collisions, continental_collisions);
        display_hmap = hmap;

        // Update the counter of iterations since last continental collision.
        last_coll_count = (last_coll_count + 1) & -(continental_collisions == 0);

        for (uint32_t i = 0; i < num_plates; ++i) {
            for (uint32_t j = 0; j < subductions[i].size(); ++j) {
                const plateCollision& coll = subductions[i][j];

                ASSERT(i != coll.index, "when subducting: SRC == DEST!");

                // Do not apply friction to oceanic plates.
                // This is a very cheap way to emulate slab pull.
                // Just perform subduction and on our way we go!
                const Platec::FloatVector source_velocity =
                    plates[coll.index]->surfaceVelocityAt(coll.wx, coll.wy);
                const uint8_t source_material =
                    platec::material::to_index(platec::material::Type::Granite);
                plates[i]->addCrustBySubduction(coll.wx, coll.wy, coll.crust, iter_count,
                                                source_velocity.x(), source_velocity.y(),
                                                source_material);
            }

            subductions[i].clear();
        }

        updateCollisions();

        fill(plate_indices_found.begin(), plate_indices_found.end(), 0);

        // Fill divergent boundaries with new crustal material, molten magma.
        if (BOOL_REGENERATE_CRUST) {
            regenerateCrust();
        }

        applyPeriodicNoise();
        applySurfaceProcesses();
        enforceSinglePlateRegions();
        removeEmptyPlates();

        // delete[] indexFound;

        // Add some "virginity buoyancy" to all pixels for a visual boost! :)
        for (uint32_t i = 0; i < (BUOYANCY_BONUS_X > 0) * map_area; ++i) {
            // Calculate the inverted age of this piece of crust.
            // Force result to be minimum between inv. age and
            // max buoyancy bonus age.
            uint32_t crust_age = iter_count - amap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE);

            hmap[i] += isOceanic(hmap[i]) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age *
                       MULINV_MAX_BUOYANCY_AGE;
        }

        for (uint32_t i = 0; i < map_area; ++i) {
            float visual_height = hmap[i];
            if (imap[i] < num_plates && amap[i] > 0) {
                uint32_t crust_age = iter_count - amap[i];
                if (crust_age <= 12U) {
                    const float previous_visual =
                        prev_display_hmap[i] > 0.0f ? prev_display_hmap[i] : initial_hmap[i];
                    const float carry_ratio = 1.0f - static_cast<float>(crust_age) / 12.0f;
                    const uint32_t x = _worldDimension.xFromIndex(i);
                    const uint32_t y = _worldDimension.yFromIndex(i);
                    const float flux_noise = sample_wrapped_octave_noise(
                        static_cast<float>(x), static_cast<float>(y), _worldDimension, 3.0f,
                        0.52f, 0.11f, 1.6f, -1.0f, 1.0f,
                        static_cast<float>(imap[i]) * 17.0f + 17.0f,
                        static_cast<float>(imap[i]) * 7.0f + 29.0f,
                        static_cast<float>(iter_count) * 0.045f,
                        static_cast<float>(imap[i]) * 19.0f + 11.0f);
                    const float flux_delta = flux_noise * (0.006f + 0.016f * carry_ratio);
                    visual_height = std::max(visual_height, previous_visual + flux_delta);
                }
            }

            display_hmap[i] = visual_height;
        }

        ++iter_count;
        updateDerivedMaps();
        if (triggerPlateBirth()) {
            enforceSinglePlateRegions();
            removeEmptyPlates();
            updateDerivedMaps();
        }
    } catch (const exception& e) {
        string msg = "Problem during update: ";
        msg = msg + e.what();
        cerr << msg << endl;
        throw runtime_error(msg.c_str());
    }
}

void lithosphere::finishCurrentState() {
    cycle_count += max_cycles > 0U;
    prev_hmap = hmap;
    prev_display_hmap = display_hmap;
    initial_hmap = display_hmap;
    finished = true;
}

void lithosphere::restart() {
    try {
        const uint32_t map_area = _worldDimension.getArea();

        cycle_count += max_cycles > 0; // No increment if running forever.
        if (cycle_count > max_cycles)
            return;

        // Update height map to include all recent changes.
        hmap.set_all(0);
        for (uint32_t i = 0; i < num_plates; ++i) {
            const uint32_t x0 = plates[i]->getLeftAsUint();
            const uint32_t y0 = plates[i]->getTopAsUint();
            const uint32_t x1 = x0 + plates[i]->getWidth();
            const uint32_t y1 = y0 + plates[i]->getHeight();

            const float* this_map;
            const uint32_t* this_age;
            plates[i]->getMap(&this_map, &this_age);

            // Copy first part of plate onto world map.
            for (uint32_t y = y0, j = 0; y < y1; ++y) {
                for (uint32_t x = x0; x < x1; ++x, ++j) {
                    const uint32_t x_mod = _worldDimension.xMod(x);
                    const uint32_t y_mod = _worldDimension.yMod(y);
                    const float h0 = hmap[_worldDimension.indexOf(x_mod, y_mod)];
                    const float h1 = this_map[j];
                    const uint32_t a0 = amap[_worldDimension.indexOf(x_mod, y_mod)];
                    const uint32_t a1 = this_age[j];

                    const float h_sum = h0 + h1;
                    // Avoid division by zero: if both heights are zero, use the new age
                    amap[_worldDimension.indexOf(x_mod, y_mod)] =
                        (h_sum > 0.0f) ? static_cast<uint32_t>((h0 * a0 + h1 * a1) / h_sum) : a1;
                    hmap[_worldDimension.indexOf(x_mod, y_mod)] += this_map[j];
                }
            }
        }
        for (uint32_t i = 0; i < map_area; ++i) {
            hmap[i] = softenHeight(i, hmap[i]);
        }
        prev_hmap = hmap;
        display_hmap = hmap;
        prev_display_hmap = display_hmap;
        initial_hmap = display_hmap;

        // create new plates IFF there are cycles left to run!
        // However, if max cycle count is "ETERNITY", then 0 < 0 + 1 always.
        if (cycle_count < max_cycles + !max_cycles) {
            clearPlates();
            createPlates();

            // Restore the ages of plates' points of crust!
            for (uint32_t i = 0; i < num_plates; ++i) {
                const uint32_t x0 = plates[i]->getLeftAsUint();
                const uint32_t y0 = plates[i]->getTopAsUint();
                const uint32_t x1 = x0 + plates[i]->getWidth();
                const uint32_t y1 = y0 + plates[i]->getHeight();

                const float* this_map;
                const uint32_t* this_age_const;
                uint32_t* this_age;

                plates[i]->getMap(&this_map, &this_age_const);
                this_age = const_cast<uint32_t*>(this_age_const);

                for (uint32_t y = y0, j = 0; y < y1; ++y) {
                    for (uint32_t x = x0; x < x1; ++x, ++j) {
                        const uint32_t x_mod = _worldDimension.xMod(x);
                        const uint32_t y_mod = _worldDimension.yMod(y);

                        this_age[j] = amap[_worldDimension.indexOf(x_mod, y_mod)];
                    }
                }
            }

            updateDerivedMaps();
            return;
        }

        // Add some "virginity buoyancy" to all pixels for a visual boost.
        for (uint32_t i = 0; i < (BUOYANCY_BONUS_X > 0) * map_area; ++i) {
            uint32_t crust_age = iter_count - amap[i];
            crust_age = MAX_BUOYANCY_AGE - crust_age;
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE);

            hmap[i] += isOceanic(hmap[i]) * BUOYANCY_BONUS_X * OCEANIC_BASE * crust_age *
                       MULINV_MAX_BUOYANCY_AGE;
        }
        for (uint32_t i = 0; i < map_area; ++i) {
            hmap[i] = softenHeight(i, hmap[i]);
        }
        prev_hmap = hmap;
        display_hmap = hmap;
        prev_display_hmap = display_hmap;
        updateDerivedMaps();
        clearPlates();
    } catch (const exception& e) {
        std::string msg = "Problem during restart: ";
        msg = msg + e.what();
        throw runtime_error(msg.c_str());
    }
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
