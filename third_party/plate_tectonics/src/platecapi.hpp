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

#ifndef PLATECAPI_H
#define PLATECAPI_H

#include "tectonic_contract.hpp"
#include "tectonic_scenario.hpp"
#include "topography_codec.hpp"
#include "utils.hpp"
#include <string.h> // For size_t.

void* platec_api_create(long seed, uint32_t width, uint32_t height, float sea_level,
                        uint32_t erosion_period, float folding_ratio, uint32_t aggr_overlap_abs,
                        float aggr_overlap_rel, uint32_t cycle_count, uint32_t num_plates,
                        float erosion_strength = 1.0f, float crust_rotation_strength = 0.20f,
                        float rotation_strength = 1.0f, float subduction_strength = 1.0f,
                        int32_t sea_level_m = -1,
                        uint16_t initial_min_height_m = TopographyCodec::kDefaultInitialMinHeightMeters,
                        uint16_t initial_max_height_m = TopographyCodec::kDefaultInitialMaxHeightMeters,
                        uint32_t cycle_step_limit = 600, float divergent_carve_strength = 0.015f);
void* platec_api_create_from_heightmap(long seed, uint32_t width, uint32_t height,
                                       const float* heightmap, float sea_level,
                                       uint32_t erosion_period, float folding_ratio,
                                       uint32_t aggr_overlap_abs, float aggr_overlap_rel,
                                       uint32_t cycle_count, uint32_t num_plates,
                                       float erosion_strength = 1.0f,
                                       float crust_rotation_strength = 0.20f,
                                       float rotation_strength = 1.0f,
                                       float subduction_strength = 1.0f,
                                       int32_t sea_level_m = -1,
                                       uint16_t initial_min_height_m = TopographyCodec::kDefaultInitialMinHeightMeters,
                                       uint16_t initial_max_height_m = TopographyCodec::kDefaultInitialMaxHeightMeters,
                                       uint32_t cycle_step_limit = 600,
                                       float divergent_carve_strength = 0.015f);
void* platec_api_create_from_scenario(const platec::scenario::Scenario* scenario);

void platec_api_destroy(void*);
const uint32_t* platec_api_get_agemap(uint32_t);
const float* platec_api_get_crust_age_myr_map(void*);
const float* platec_api_get_crust_thickness_map(void*);
const uint8_t* platec_api_get_crust_class_map(void*);
const float* platec_api_get_uplift_tendency_map(void*);
const float* platec_api_get_subsidence_tendency_map(void*);
const float* platec_api_get_accumulated_strain_map(void*);
const uint8_t* platec_api_get_boundary_type_map(void*);
const uint16_t* platec_api_get_boundary_distance_map(void*);
const uint32_t* platec_api_get_boundary_segment_id_map(void*);
const uint32_t* platec_api_get_nearest_boundary_id_map(void*);
uint32_t platec_api_get_boundary_segment_count(void*);
const platec::contract::BoundarySegment* platec_api_get_boundary_segments(void*);
const uint32_t* platec_api_get_deforming_region_id_map(void*);
const uint8_t* platec_api_get_deforming_region_type_map(void*);
const float* platec_api_get_deformation_rate_map(void*);
const float* platec_api_get_deformation_velocity_x_map(void*);
const float* platec_api_get_deformation_velocity_y_map(void*);
uint32_t platec_api_get_deforming_region_count(void*);
const platec::contract::DeformingRegion* platec_api_get_deforming_regions(void*);
float* platec_api_get_heightmap(void*);
uint32_t* platec_api_get_platesmap(void*);
const uint8_t* platec_api_get_convergence_map(void*);
const uint8_t* platec_api_get_divergence_map(void*);
const uint8_t* platec_api_get_shear_map(void*);
const uint8_t* platec_api_get_geologic_regime_map(void*);
uint32_t platec_api_get_cycle_count(void*);
uint32_t platec_api_get_iteration_count(void*);
uint32_t platec_api_get_time_origin_step(void*);
double platec_api_get_time_myr(void*);
double platec_api_get_delta_time_myr(void*);
uint32_t platec_api_is_finished(void*);
void platec_api_step(void*);
void platec_api_load_heightmap(void*, const float* normalized_heightmap, float sea_level);
void platec_api_load_heightmap_raw(void*, const float* normalized_heightmap);
void platec_api_load_heightmap_u16(void*, const uint16_t* heightmap_m, uint16_t sea_level_m);
uint16_t platec_api_get_sea_level_m(void*);
uint32_t platec_api_get_plate_count(void*);

float platec_api_velocity_unity_vector_x(void*, uint32_t plate_index);
float platec_api_velocity_unity_vector_y(void*, uint32_t plate_index);
float platec_api_velocity_vector_x(void*, uint32_t plate_index);
float platec_api_velocity_vector_y(void*, uint32_t plate_index);
float platec_api_angular_velocity(void*, uint32_t plate_index);
float platec_api_mass_center_x(void*, uint32_t plate_index);
float platec_api_mass_center_y(void*, uint32_t plate_index);

uint32_t lithosphere_getMapWidth(void* object);
uint32_t lithosphere_getMapHeight(void* object);

#endif
