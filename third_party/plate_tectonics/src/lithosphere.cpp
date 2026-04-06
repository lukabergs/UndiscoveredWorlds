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

float wrap_coordinate(float value, float period) {
    float wrapped = std::fmod(value, period);
    if (wrapped < 0.0f) {
        wrapped += period;
    }
    return wrapped;
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
      cycle_step_limit(params.cycle_step_limit), max_plates(params.plate_count), num_plates(0),
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
      _worldDimension(params.width, params.height), _randsource(params.seed), _steps(0) {
    if (params.width < 5 || params.height < 5) {
        throw runtime_error("Width and height should be >=5");
    }
    if (!std::isfinite(params.delta_time_myr) || params.delta_time_myr <= 0.0) {
        throw runtime_error("delta_time_myr must be finite and > 0");
    }
    if (initial_min_height_m >= initial_max_height_m) {
        throw runtime_error("Initial minimum height must be lower than the initial maximum height");
    }
    if (static_cast<uint32_t>(initial_max_height_m) - static_cast<uint32_t>(initial_min_height_m) <
        2U) {
        throw runtime_error(
            "Initial height range must leave room for both ocean and land values");
    }
    if (params.sea_level_m != TopographyCodec::kNoSeaLevelOverride &&
        (params.sea_level_m <= static_cast<int32_t>(initial_min_height_m) ||
         params.sea_level_m >= static_cast<int32_t>(initial_max_height_m))) {
        throw runtime_error("Sea level override must lie between the initial min and max heights");
    }

    collisions.resize(max_plates);
    subductions.resize(max_plates);

    // Create default plates
    plates = new plate*[max_plates];
    for (uint32_t i = 0; i < max_plates; i++) {
        plate_areas[i].border.reserve(8);
    }
    seedInitialTopography(params.sea_level, params.sea_level_m);
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

        // Initialize "Free plate center position" lookup table.
        // This way two plate centers will never be identical.
        for (uint32_t i = 0; i < map_area; ++i)
            imap[i] = i;

        // Select N plate centers from the global map.

        for (uint32_t i = 0; i < num_plates; ++i) {
            plateArea& area = plate_areas[i];

            // Randomly select an unused plate origin.
            const uint32_t p = imap[(uint32_t)_randsource.next() % (map_area - i)];
            const uint32_t y = _worldDimension.yFromIndex(p);
            const uint32_t x = _worldDimension.xFromIndex(p);

            area.lft = area.rgt = x; // Save origin...
            area.top = area.btm = y;
            area.wdt = area.hgt = 1;

            area.border.clear();
            area.border.push_back(p); // ...and mark it as border.

            // Overwrite used entry with last unused entry in array.
            imap[p] = imap[map_area - i - 1];
        }

        imap.set_all(0xFFFFFFFF);

        growPlates();
        jaggedizePlateBoundaries();

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

            // Copy plate's height data from global map into local map.
            for (uint32_t y = y0, j = 0; y < y1; ++y) {
                for (uint32_t x = x0; x < x1; ++x, ++j) {
                    uint32_t k = _worldDimension.normalizedIndexOf(x, y);
                    pmap[j] = hmap[k] * (imap[k] == i);
                }
            }
            // Create plate.
            // MK: The pmap array becomes owned by map, do not delete it
            plates[i] = new plate(_randsource.next(), pmap, width, height, x0, y0, i,
                                  _worldDimension, erosion_strength,
                                  crust_rotation_strength, rotation_strength);
        }

        iter_count = num_plates + MAX_BUOYANCY_AGE;
        if (time_origin_step == 0U) {
            time_origin_step = iter_count;
        }
        peak_Ek = 0;
        last_coll_count = 0;
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

const float* lithosphere::getUpliftTendencyMap() const noexcept {
    return uplift_tendency_map.data();
}

const float* lithosphere::getSubsidenceTendencyMap() const noexcept {
    return subsidence_tendency_map.data();
}

const float* lithosphere::getAccumulatedStrainMap() const noexcept {
    return accumulated_strain_map.data();
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
    plate_kinematics.clear();
    boundary_segments.clear();
    deforming_regions.clear();
    next_boundary_segment_id = platec::contract::kNoBoundaryId + 1U;

    for (uint32_t i = 0; i < max_plates; ++i) {
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

        geologic_regime_map[index] = static_cast<uint8_t>(resolved_regime);
    }
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
        const bool continental_plate = current_plate->isContinentalPlate();
        const uint32_t x0 = current_plate->getLeftAsUint();
        const uint32_t y0 = current_plate->getTopAsUint();
        const uint32_t x1 = x0 + current_plate->getWidth();
        const uint32_t y1 = y0 + current_plate->getHeight();

        const float* plate_map = nullptr;
        const uint32_t* plate_age = nullptr;
        current_plate->getMap(&plate_map, &plate_age);
        (void)plate_age;

        uint32_t x_mod_start = x0 % width;
        uint32_t y_mod = y0 % height;
        for (uint32_t y = y0, j = 0; y < y1;
             ++y, y_mod = ++y_mod >= height ? y_mod - height : y_mod) {
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
                if (continental_plate) {
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
    }

    updateBoundaryGraphMaps();
    updateDeformingRegionMaps();

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

        const float delta_surface = display_hmap[index] - prev_display_hmap[index];
        const float positive_surface = clamp_unit(delta_surface / kSurfaceDeltaScale);
        const float negative_surface = clamp_unit(-delta_surface / kSurfaceDeltaScale);
        const float age_old = clamp_unit(crust_age_myr_map[index] / kOldCrustMyr);
        const float age_young = 1.0f - clamp_unit(crust_age_myr_map[index] / kYoungCrustMyr);
        const float rift_shoulder = deforming_region_type ==
                                            platec::contract::DeformingRegionType::ContinentalRift &&
                                        boundary_distance_map[index] > 0U
                                    ? clamp_unit(static_cast<float>(boundary_distance_map[index]) /
                                                 3.0f)
                                    : 0.0f;

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
            uplift += 0.10f * deformation_rate * rift_shoulder;
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
    updateTectonicProvenance();
    updateFieldStackMaps();
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
    const uint16_t provisional_sea_level_m =
        has_override
            ? static_cast<uint16_t>(sea_level_m_override)
            : static_cast<uint16_t>(initial_min_height_m +
                                    (initial_max_height_m - initial_min_height_m) / 2U);

    const FractalNoiseConfig warp_x_noise = make_noise_config(_randsource, 4.0f, 0.58f, 0.45f, 1.6f);
    const FractalNoiseConfig warp_y_noise = make_noise_config(_randsource, 4.0f, 0.58f, 0.45f, 1.6f);
    const FractalNoiseConfig continent_macro_noise =
        make_noise_config(_randsource, 5.0f, 0.57f, 0.16f, 0.95f);
    const FractalNoiseConfig continent_detail_noise =
        make_noise_config(_randsource, 6.0f, 0.59f, 0.42f, 1.85f);
    const FractalNoiseConfig coast_breakup_noise =
        make_noise_config(_randsource, 6.0f, 0.57f, 0.62f, 2.2f);
    const FractalNoiseConfig shelf_noise = make_noise_config(_randsource, 6.0f, 0.61f, 0.85f, 2.9f);
    const FractalNoiseConfig sea_base_noise = make_noise_config(_randsource, 7.0f, 0.63f, 1.15f, 3.9f);
    const FractalNoiseConfig sea_detail_noise =
        make_noise_config(_randsource, 7.0f, 0.61f, 1.95f, 5.6f);
    const FractalNoiseConfig sea_micro_noise = make_noise_config(_randsource, 5.0f, 0.48f, 3.95f, 10.2f);
    const FractalNoiseConfig land_noise = make_noise_config(_randsource, 6.0f, 0.56f, 0.95f, 3.6f);
    const FractalNoiseConfig foothill_noise =
        make_noise_config(_randsource, 6.0f, 0.53f, 1.55f, 4.7f);
    const FractalNoiseConfig mountain_noise =
        make_noise_config(_randsource, 5.0f, 0.47f, 2.45f, 6.3f);
    const FractalNoiseConfig ridge_noise = make_noise_config(_randsource, 5.0f, 0.50f, 2.6f, 7.4f);
    const FractalNoiseConfig trench_noise = make_noise_config(_randsource, 5.0f, 0.48f, 3.8f, 9.8f);

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

    const float ocean_target_min = static_cast<float>(initial_min_height_m);
    const float ocean_target_max = static_cast<float>(provisional_sea_level_m - 1U);
    const float land_target_min = static_cast<float>(provisional_sea_level_m);
    const float land_target_max = static_cast<float>(initial_max_height_m);

    std::vector<uint16_t> metric_heightmap(map_area);
    uint16_t ocean_metric_max = 0;
    uint16_t land_metric_min = TopographyCodec::kMaxHeightMeters;
    uint32_t ocean_count = 0;
    uint32_t land_count = 0;
    for (uint32_t i = 0; i < map_area; ++i) {
        const float meters = is_land[i]
            ? remap_to_range(raw_height[i], raw_land_min, raw_land_max, land_target_min,
                             land_target_max)
            : remap_to_range(raw_height[i], raw_ocean_min, raw_ocean_max, ocean_target_min,
                             ocean_target_max);
        const uint16_t sample = static_cast<uint16_t>(std::lround(meters));
        metric_heightmap[i] = sample;
        if (is_land[i]) {
            land_metric_min = std::min(land_metric_min, sample);
            ++land_count;
        } else {
            ocean_metric_max = std::max(ocean_metric_max, sample);
            ++ocean_count;
        }
    }

    const uint16_t initial_sea_level_m =
        has_override
            ? provisional_sea_level_m
            : (ocean_count > 0U && land_count > 0U)
                ? static_cast<uint16_t>(std::min<uint32_t>(
                      land_metric_min,
                      std::max<uint32_t>(
                          static_cast<uint32_t>(ocean_metric_max) + 1U,
                          TopographyCodec::infer_metric_sea_level(metric_heightmap.data(),
                                                                  metric_heightmap.size(),
                                                                  ocean_coverage))))
                : provisional_sea_level_m;
    initializeHeightMapFromMetric(metric_heightmap.data(), initial_sea_level_m);
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
    resetSimulationState();
    createPlates();
}

void lithosphere::importMetricHeightMap(const uint16_t* heightmap_m, uint16_t new_sea_level_m) {
    initializeHeightMapFromMetric(heightmap_m, new_sea_level_m);
    resetSimulationState();
    createPlates();
}

bool lithosphere::isFinished() const {
    return getPlateCount() == 0;
}

// At least two plates are at same location.
// Move some crust from the SMALLER plate onto LARGER one.
void lithosphere::resolveJuxtapositions(const uint32_t& i, const uint32_t& j, const uint32_t& k,
                                        const uint32_t& x_mod, const uint32_t& y_mod,
                                        const float*& this_map, const uint32_t*& this_age,
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
    const float penetrating_height =
        penetrating_owner == current_owner ? hmap[k] : this_map[j];
    const uint32_t penetrating_age =
        penetrating_owner == current_owner ? amap[k] : this_age[j];

    const float folded_crust = std::min(surface_height, penetrating_height) * folding_ratio;
    const float uplifted_surface = surface_height + folded_crust;
    const float remaining_penetrating = std::max(0.0f, penetrating_height - folded_crust);

    plates[surface_owner]->setCrust(x_mod, y_mod, uplifted_surface,
                                    std::max(surface_age, penetrating_age));
    plates[penetrating_owner]->setCrust(x_mod, y_mod, remaining_penetrating, penetrating_age);

    hmap[k] = uplifted_surface;
    imap[k] = surface_owner;
    amap[k] = plates[surface_owner]->getCrustTimestamp(x_mod, y_mod);

    collisions[penetrating_owner].push_back(
        plateCollision(surface_owner, x_mod, y_mod, folded_crust));
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
    hmap.set_all(0);
    imap.set_all(0xFFFFFFFF);
    for (uint32_t i = 0; i < num_plates; ++i) {
        const uint32_t x0 = plates[i]->getLeftAsUint();
        const uint32_t y0 = plates[i]->getTopAsUint();
        const uint32_t x1 = x0 + plates[i]->getWidth();
        const uint32_t y1 = y0 + plates[i]->getHeight();

        const float* this_map;
        const uint32_t* this_age;
        plates[i]->getMap(&this_map, &this_age);

        uint32_t x_mod_start = (x0 + world_width) % world_width;
        uint32_t y_mod = (y0 + world_height) % world_height;

        // Copy first part of plate onto world map.
        // MK: These loops are ugly, but using modulus in here is a hog
        for (uint32_t y = y0, j = 0; y < y1;
             ++y, y_mod = ++y_mod >= world_height ? y_mod - world_height : y_mod) {
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

                    continue;
                }

                const uint32_t current_owner = imap[k];
                const float current_height = hmap[k];
                const uint32_t current_age = amap[k];
                const float incoming_height = this_map[j];
                const uint32_t incoming_age = this_age[j];
                plate& current_plate = *plates[current_owner];
                plate& incoming_plate = *plates[i];

                if (current_plate.isContinentalPlate() && incoming_plate.isContinentalPlate()) {
                    resolveJuxtapositions(i, j, k, x_mod, y_mod, this_map, this_age,
                                          continental_collisions);
                    continue;
                }

                uint32_t overriding_owner = current_owner;
                uint32_t subducting_owner = i;

                if (current_plate.isOceanicPlate() && incoming_plate.isOceanicPlate()) {
                    const float current_buoyancy = current_plate.buoyancy();
                    const float incoming_buoyancy = incoming_plate.buoyancy();
                    if (incoming_buoyancy > current_buoyancy + FLT_EPSILON ||
                        (std::fabs(incoming_buoyancy - current_buoyancy) <= FLT_EPSILON &&
                         (incoming_height > current_height + FLT_EPSILON ||
                          (std::fabs(incoming_height - current_height) <= FLT_EPSILON &&
                           incoming_age >= current_age)))) {
                        overriding_owner = i;
                        subducting_owner = current_owner;
                    }
                } else if (current_plate.isContinentalPlate()) {
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
                const float subducting_height =
                    subducting_owner == current_owner ? current_height : incoming_height;
                const uint32_t subducting_age =
                    subducting_owner == current_owner ? current_age : incoming_age;
                const float subducted_crust =
                    std::min(subducting_height, OCEANIC_BASE * subduction_strength);

                if (subducted_crust <= FLT_EPSILON) {
                    resolveJuxtapositions(i, j, k, x_mod, y_mod, this_map, this_age,
                                          continental_collisions);
                    continue;
                }

                current_plate.addCollision(x_mod, y_mod);
                incoming_plate.addCollision(x_mod, y_mod);

                const float convergence_fold_scale =
                    current_plate.isOceanicPlate() && incoming_plate.isOceanicPlate() ? 0.75f : 1.0f;
                const float folded_crust =
                    std::min(overriding_height, subducting_height) * folding_ratio *
                    convergence_fold_scale;
                const float sediment =
                    SUBDUCT_RATIO * subducted_crust *
                    (CONTINENTAL_BASE - std::min(subducting_height, CONTINENTAL_BASE)) /
                    CONTINENTAL_BASE;
                const float uplifted_overriding =
                    overriding_height + folded_crust + sediment * 0.5f;
                const float remaining_subducting =
                    std::max(0.0f, subducting_height - subducted_crust);

                plates[overriding_owner]->setCrust(
                    x_mod, y_mod, uplifted_overriding,
                    std::max(overriding_age, subducting_age));
                plates[subducting_owner]->setCrust(x_mod, y_mod, remaining_subducting,
                                                   subducting_age);

                hmap[k] = uplifted_overriding;
                imap[k] = overriding_owner;
                amap[k] = plates[overriding_owner]->getCrustTimestamp(x_mod, y_mod);

                subductions[overriding_owner].push_back(
                    plateCollision(subducting_owner, x_mod, y_mod, sediment));
                collisions[subducting_owner].push_back(
                    plateCollision(overriding_owner, x_mod, y_mod,
                                   folded_crust + subducted_crust * 0.35f));
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
            plates[i]->applyFriction(coll.crust);
            plates[coll.index]->applyFriction(coll.crust);

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

    auto assign_new_crust = [&](uint32_t index) {
        if (imap[index] < num_plates) {
            return;
        }

        const uint32_t x = _worldDimension.xFromIndex(index);
        const uint32_t y = _worldDimension.yFromIndex(index);
        imap[index] = chooseDivergentOwner(x, y, index);
        amap[index] = iter_count;
        if (imap[index] >= num_plates) {
            return;
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
        float oceanic_boundary_weight = 0.0f;
        float continental_boundary_weight = 0.0f;

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

            const float type_weight = sample.dir_weight * (0.35f + sample.retreat);
            if (plates[sample.owner]->isOceanicPlate()) {
                oceanic_boundary_weight += type_weight;
            } else {
                continental_boundary_weight += type_weight;
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
        const float previous_height =
            std::max(OCEANIC_BASE,
                     prev_hmap[index] > 0.0f ? prev_hmap[index] : initial_hmap[index]);
        const float random_signed = _randsource.next_float() * 2.0f - 1.0f;
        const float relief_step =
            divergent_carve_strength *
            (0.85f + 0.25f * owner_share + 0.20f * opening_strength +
             0.10f * std::min(1.5f, drift_magnitude) + 0.35f * random_signed);
        const float boundary_weight_sum = oceanic_boundary_weight + continental_boundary_weight;
        const float boundary_type_bias =
            boundary_weight_sum > FLT_EPSILON
                ? (oceanic_boundary_weight - continental_boundary_weight) / boundary_weight_sum
                : 0.0f;
        const float generated_candidate = previous_height + relief_step * boundary_type_bias;
        const float min_generated =
            boundary_type_bias >= 0.0f ? OCEANIC_BASE * 0.50f : OCEANIC_BASE * 0.10f;
        const float max_generated =
            boundary_type_bias > 0.0f ? CONTINENTAL_BASE - 0.05f : previous_height;
        const float generated_crust =
            std::max(min_generated, std::min(max_generated, generated_candidate));

        hmap[index] = generated_crust;
        plates[imap[index]]->setCrust(x, y, generated_crust, iter_count);
        ++plate_indices_found[imap[index]];
    };

    for (uint32_t index = 0; index < map_area; ++index) {
        if (imap[index] < num_plates) {
            if (++plate_indices_found[imap[index]] && hmap[index] <= 0) {
                puts("Occupied point has no land mass!");
                exit(1);
            }
            continue;
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
        if (totalVelocity < RESTART_SPEED_LIMIT ||
            systemKineticEnergy / peak_Ek < RESTART_ENERGY_RATIO ||
            last_coll_count > NO_COLLISION_TIME_LIMIT ||
            (cycle_step_limit > 0 && _steps > static_cast<int>(cycle_step_limit))) {
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
                plates[i]->addCrustBySubduction(coll.wx, coll.wy, coll.crust, iter_count,
                                                source_velocity.x(), source_velocity.y());
            }

            subductions[i].clear();
        }

        updateCollisions();

        fill(plate_indices_found.begin(), plate_indices_found.end(), 0);

        // Fill divergent boundaries with new crustal material, molten magma.
        if (BOOL_REGENERATE_CRUST) {
            regenerateCrust();
        }

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

            display_hmap[i] =
                std::max(0.0f, std::min(CONTINENTAL_BASE + 1.0f, visual_height));
        }

        ++iter_count;
        updateDerivedMaps();
    } catch (const exception& e) {
        string msg = "Problem during update: ";
        msg = msg + e.what();
        cerr << msg << endl;
        throw runtime_error(msg.c_str());
    }
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
