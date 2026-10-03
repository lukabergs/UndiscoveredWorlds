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

#include "tectonic_api.hpp"
#include "lithosphere.hpp"
#include "plate.hpp"
#include <stdio.h>
#include <stdlib.h>

#include <vector>

void* tectonics_create_from_scenario(const tectonics::scenario::Scenario* scenario) {
    if (scenario == nullptr) {
        return nullptr;
    }

    lithosphere* litho = new lithosphere(*scenario);
    return litho;
}

void tectonics_destroy(void* pointer) {
    delete static_cast<lithosphere*>(pointer);
}

const float* tectonics_get_crust_age_myr_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getCrustAgeMyrMap();
}

const float* tectonics_get_crust_thickness_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getCrustThicknessMap();
}

const uint8_t* tectonics_get_crust_class_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getCrustClassMap();
}

const uint8_t* tectonics_get_surface_crust_type_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getSurfaceCrustTypeMap();
}

const float* tectonics_get_uplift_tendency_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getUpliftTendencyMap();
}

const float* tectonics_get_subsidence_tendency_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getSubsidenceTendencyMap();
}

const float* tectonics_get_accumulated_strain_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getAccumulatedStrainMap();
}

const float* tectonics_get_tension_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getTensionMap();
}

const uint8_t* tectonics_get_boundary_type_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getBoundaryTypeMap();
}

const uint16_t* tectonics_get_boundary_distance_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getBoundaryDistanceMap();
}

const uint32_t* tectonics_get_boundary_segment_id_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getBoundarySegmentIdMap();
}

const uint32_t* tectonics_get_nearest_boundary_id_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getNearestBoundaryIdMap();
}

uint32_t tectonics_get_boundary_segment_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getBoundarySegmentCount();
}

const tectonics::contract::BoundarySegment* tectonics_get_boundary_segments(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getBoundarySegments();
}

const uint32_t* tectonics_get_deforming_region_id_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformingRegionIdMap();
}

const uint8_t* tectonics_get_deforming_region_type_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformingRegionTypeMap();
}

const float* tectonics_get_deformation_rate_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformationRateMap();
}

const float* tectonics_get_deformation_velocity_x_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformationVelocityXMap();
}

const float* tectonics_get_deformation_velocity_y_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformationVelocityYMap();
}

uint32_t tectonics_get_deforming_region_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformingRegionCount();
}

const tectonics::contract::DeformingRegion* tectonics_get_deforming_regions(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeformingRegions();
}

uint32_t tectonics_get_tectonic_junction_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getTectonicJunctionCount();
}

const tectonics::contract::TectonicJunction* tectonics_get_tectonic_junctions(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getTectonicJunctions();
}

const uint8_t* tectonics_get_material_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getMaterialMap();
}

const float* tectonics_get_plume_influence_map(void* pointer) {
    const auto& map = static_cast<lithosphere*>(pointer)->getPlumeInfluenceMap();
    return map.empty() ? nullptr : map.data();
}

float* tectonics_get_heightmap(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    float* res = litho->getTopography();
    return res;
}

uint32_t* tectonics_get_platesmap(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    uint32_t* res = litho->getPlatesMap();
    return res;
}

const uint8_t* tectonics_get_convergence_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getConvergenceMap();
}

const uint8_t* tectonics_get_divergence_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDivergenceMap();
}

const uint8_t* tectonics_get_shear_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getShearMap();
}

const uint8_t* tectonics_get_geologic_regime_map(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getGeologicRegimeMap();
}

uint32_t tectonics_get_iteration_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getIterationCount();
}

uint32_t tectonics_get_time_origin_step(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getTimeOriginStep();
}

double tectonics_get_time_myr(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getTimeMyr();
}

double tectonics_get_delta_time_myr(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getDeltaTimeMyr();
}

uint32_t tectonics_is_finished(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    if (litho->isFinished()) {
        return 1;
    } else {
        return 0;
    }
}

void tectonics_step(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->update();
}

void tectonics_load_heightmap(void* pointer, const float* normalized_heightmap, float sea_level) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->importNormalizedHeightMap(normalized_heightmap, sea_level);
}

void tectonics_load_heightmap_raw(void* pointer, const float* normalized_heightmap) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->importRawHeightMap(normalized_heightmap);
}

void tectonics_load_heightmap_u16(void* pointer, const uint16_t* heightmap_m, uint16_t sea_level_m) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    litho->importMetricHeightMap(heightmap_m, sea_level_m);
}

uint16_t tectonics_get_sea_level_m(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getSeaLevelMeters();
}

uint32_t tectonics_get_plate_count(void* pointer) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlateCount();
}

uint8_t tectonics_is_oceanic_plate(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->isOceanicPlate() ? 1U : 0U;
}

uint32_t tectonics_get_map_width(void* pointer) {
    return static_cast<lithosphere*>(pointer)->getWidth();
}

uint32_t tectonics_get_map_height(void* pointer) {
    return static_cast<lithosphere*>(pointer)->getHeight();
}

const tectonics::columns::Cell* tectonics_get_column_cells(void* pointer) {
    const auto& cells = static_cast<const lithosphere*>(pointer)->getColumnModel().cells();
    return cells.empty() ? nullptr : cells.data();
}

tectonics::contract::Snapshot tectonics_capture_snapshot(void* pointer) {
    return tectonics::contract::capture_snapshot(*static_cast<const lithosphere*>(pointer));
}

float tectonics_velocity_unity_vector_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->velocityUnitVector().x();
}

float tectonics_velocity_unity_vector_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->velocityUnitVector().y();
}

float tectonics_velocity_vector_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->linearVelocityVector().x();
}

float tectonics_velocity_vector_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->linearVelocityVector().y();
}

float tectonics_angular_velocity(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->getAngularVelocity();
}

float tectonics_mass_center_x(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->worldMassCenter().getX();
}

float tectonics_mass_center_y(void* pointer, uint32_t plate_index) {
    lithosphere* litho = static_cast<lithosphere*>(pointer);
    return litho->getPlate(plate_index)->worldMassCenter().getY();
}

uint64_t tectonics_get_plate_uid(void* pointer, uint32_t index) {
    const auto& ids = static_cast<lithosphere*>(pointer)->getTopologyDiagnostics().stable_ids;
    return index < ids.size() ? ids[index] : 0;
}

const tectonics::transport::Displacement* tectonics_get_material_displacement(void* pointer) {
    const auto& displacement = static_cast<lithosphere*>(pointer)->getMaterialDisplacement();
    return displacement.empty() ? nullptr : displacement.data();
}
