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

#ifndef TECTONIC_API_HPP
#define TECTONIC_API_HPP

#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"
#include "simulations/geology/plate_tectonics/tectonic_scenario.hpp"
#include "simulations/geology/plate_tectonics/topography_codec.hpp"
#include "simulations/geology/plate_tectonics/utils.hpp"
#include "core/mesh/cell_grid.hpp"
#include <memory>
#include <string.h> // For size_t.

void* tectonics_create_from_scenario(const tectonics::scenario::Scenario* scenario);

void tectonics_destroy(void*);
const float* tectonics_get_crust_age_myr_map(void*);
const float* tectonics_get_crust_thickness_map(void*);
const uint8_t* tectonics_get_crust_class_map(void*);
const uint8_t* tectonics_get_surface_crust_type_map(void*);
const float* tectonics_get_uplift_tendency_map(void*);
const float* tectonics_get_rock_uplift_rate_map(void*); // m/Myr, signed; tectonic sub-steps only, in the plate's frame (semi-Lagrangian: carried with the rock by the engine's own remap, exact for translation where one column field varies, approximate at plate boundaries; columns::Model::upliftRate)
const float* tectonics_get_shortening_map(void*); // m across one native cell width, >= 0; cumulative over the run, in the plate's frame (columns::Model::strain, shortening_m)
const float* tectonics_get_shortening_direction_map(void*); // radians in [0, pi) from x toward y, the shortening axis (columns::shortening_direction)
const float* tectonics_get_subsidence_tendency_map(void*);
const float* tectonics_get_accumulated_strain_map(void*);
const float* tectonics_get_tension_map(void*);
const uint8_t* tectonics_get_boundary_type_map(void*);
const uint16_t* tectonics_get_boundary_distance_map(void*);
const uint32_t* tectonics_get_boundary_segment_id_map(void*);
const uint32_t* tectonics_get_nearest_boundary_id_map(void*);
uint32_t tectonics_get_boundary_segment_count(void*);
const tectonics::contract::BoundarySegment* tectonics_get_boundary_segments(void*);
const uint32_t* tectonics_get_deforming_region_id_map(void*);
const uint8_t* tectonics_get_deforming_region_type_map(void*);
const float* tectonics_get_deformation_rate_map(void*);
const float* tectonics_get_deformation_velocity_x_map(void*);
const float* tectonics_get_deformation_velocity_y_map(void*);
uint32_t tectonics_get_deforming_region_count(void*);
const tectonics::contract::DeformingRegion* tectonics_get_deforming_regions(void*);
uint32_t tectonics_get_tectonic_junction_count(void*);
const tectonics::contract::TectonicJunction* tectonics_get_tectonic_junctions(void*);
const uint8_t* tectonics_get_material_map(void*);
float* tectonics_get_heightmap(void*);
// Proximity (0..1) to an active mantle plume over the last update; nullptr when plumes are off.
const float* tectonics_get_plume_influence_map(void* pointer);
uint32_t* tectonics_get_platesmap(void*);
const uint8_t* tectonics_get_convergence_map(void*);
const uint8_t* tectonics_get_divergence_map(void*);
const uint8_t* tectonics_get_shear_map(void*);
const uint8_t* tectonics_get_geologic_regime_map(void*);
uint32_t tectonics_get_iteration_count(void*);
uint32_t tectonics_get_time_origin_step(void*);
double tectonics_get_time_myr(void*);
double tectonics_get_delta_time_myr(void*);
uint32_t tectonics_is_finished(void*);
void tectonics_step(void*);
void tectonics_load_heightmap(void*, const float* normalized_heightmap, float sea_level);
void tectonics_load_heightmap_raw(void*, const float* normalized_heightmap);
void tectonics_load_heightmap_u16(void*, const uint16_t* heightmap_m, uint16_t sea_level_m);
uint16_t tectonics_get_sea_level_m(void*);
uint32_t tectonics_get_plate_count(void*);
uint8_t tectonics_is_oceanic_plate(void*, uint32_t plate_index);

// The published plate kinematics (contract::PlateKinematics, schema 8: omega in rad/s, the
// centroid's longitude and latitude in degrees, its velocity (east, north) in cm/yr), one per
// plate index below tectonics_get_plate_kinematics_count; null when there are none.
uint32_t tectonics_get_plate_kinematics_count(void*);
const tectonics::contract::PlateKinematics* tectonics_get_plate_kinematics(void*);
// The engine's cube: the mesh::CellGrid of its cells (uw-cube-1, N = sphere_face_cells, the
// scenario's radius), shared so a hand-off can cache weights against it.
std::shared_ptr<const mesh::CellGrid> tectonics_get_cell_grid(void*);
// Deterministic topology identity survives index changes after plate removal.
uint64_t tectonics_get_plate_uid(void*, uint32_t plate_index);
// Per-cell material displacement of the latest isostatic-column update, (east, north) in cell
// spacings per update, in id order; null before the first physical update.
const tectonics::transport::Displacement* tectonics_get_material_displacement(void*);

// The id-order view of the cube's 6 N^2 cells: N columns x 6 N rows.
uint32_t tectonics_get_map_width(void*);
uint32_t tectonics_get_map_height(void*);
// Isostatic column state, one per cell in id order; null before columns exist.
const tectonics::columns::Cell* tectonics_get_column_cells(void*);
tectonics::contract::Snapshot tectonics_capture_snapshot(void*);

#endif
