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

#ifndef LITHOSPHERE_HPP
#define LITHOSPHERE_HPP

#include <array>
#include <cstdint>
#include <cstring> // For size_t.
#include <memory>
#include <stdexcept>
#include <vector>
#include "simulations/geology/plate_tectonics/material.hpp"
#include "simulations/geology/plate_tectonics/simplerandom.hpp"
#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"
#include "simulations/geology/plate_tectonics/tectonic_scenario.hpp"
#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
#include "simulations/geology/plate_tectonics/topography_codec.hpp"
#include <cmath>

using namespace std;

#define CONTINENTAL_BASE 1.0f
#define OCEANIC_BASE 0.1f

namespace tectonics {

// One plate (W2-E, E1): its identity, rotation, centroid, area and continental share. On the
// planar scaffold the motion is the torus engine's: a translation in cells per update kept as a
// unit direction and a speed (so the published velocities are unchanged); the rotation vector
// stays zero until the spherical kinematics (E2C).
struct PlateState {
    uint32_t stable_id = 0;
    mesh::Vec3 rotation_rad_s;   // omega
    mesh::Vec3 centroid;         // grid::TectonicGrid::centroids
    double area_m2 = 0;
    // Owned cells with crust whose owner material is continental, over owned cells with crust
    // (the torus engine's continentalAreaRatio).
    float continental_share = 0;
    float direction_x = 1, direction_y = 0;  // planar unit direction
    float speed = 0;                         // planar speed, cells per update

    // A new plate's motion drawn from its seed: speed 1 in a uniform direction (as the torus).
    static PlateState fromSeed(uint32_t seed);
    float velocityX() const { return direction_x * speed; }
    float velocityY() const { return direction_y * speed; }
    // Speed and unit direction of (x, y); below 1e-6 the plate stops facing +x.
    void setVelocity(float x, float y);
    bool isOceanic() const { return continental_share <= 0.10f; }
};

// Wall-clock seconds per engine phase, summed over the run.
struct PhaseTimings {
    static constexpr int count = 16;
    static const char* name(int phase);
    std::array<double, count> seconds{};
    uint64_t updates = 0;
};

} // namespace tectonics

/**
 * Lithosphere is the rigid outermost shell of a rocky planet.
 *
 * Lithosphere is divided into several rigid areas i.e. plates. In time passing
 * the topography of the planet evolves as the result of plate dynamics.
 * Lithosphere is the class responsible for creating and managing all the
 * plates. It updates heightmap to match the current setup of plates and thus
 * offers caller an easy access to system's topography.
 *
 */
class lithosphere {
  public:
    explicit lithosphere(const tectonics::scenario::Scenario& scenario) noexcept(false);

    ~lithosphere() noexcept; ///< Standard destructor.

    /**
     * Split the current topography into given number of (rigid) plates.
     *
     * Any previous set of plates is discarded.
     *
     * @param num_plates Number of areas the surface is divided into.
     */
    void createPlates();

    uint32_t getIterationCount() const noexcept {
        return iter_count;
    }
    uint32_t getTimeOriginStep() const noexcept {
        return time_origin_step;
    }
    double getTimeMyr() const noexcept;
    double getDeltaTimeMyr() const noexcept {
        return scenario.delta_time_myr;
    }
    const tectonics::scenario::Scenario& getScenario() const noexcept {
        return scenario;
    }
    const tectonics::dynamics::Diagnostics& getMotionDiagnostics() const noexcept {
        return motion_diagnostics;
    }
    const tectonics::grid::TectonicGrid& getGrid() const noexcept {
        return *grid_;
    }
    const tectonics::PhaseTimings& getPhaseTimings() const noexcept {
        return phase_timings;
    }
    uint32_t getPlateCount() const noexcept;    ///< Return number of plates.
    uint32_t getProvenancePlateCount() const noexcept;
    const uint32_t* getAgeMap() const noexcept; ///< Return surface age map.
    const float* getCrustAgeMyrMap() const noexcept;
    const float* getCrustThicknessMap() const noexcept;
    const uint8_t* getCrustClassMap() const noexcept;
    const uint8_t* getSurfaceCrustTypeMap() const noexcept;
    const float* getUpliftTendencyMap() const noexcept;
    /// Rock uplift rate, m/Myr, signed (columns::Model::upliftRate), republished each update.
    const float* getRockUpliftRateMap() const noexcept;
    /// Cumulative shortening, m across one native cell width, and its axis, radians in [0, pi)
    /// from x toward y (columns::Model::strain, shortening_m, shortening_direction), republished each update.
    const float* getShorteningMap() const noexcept;
    const float* getShorteningDirectionMap() const noexcept;
    const float* getSubsidenceTendencyMap() const noexcept;
    const float* getAccumulatedStrainMap() const noexcept;
    const float* getTensionMap() const noexcept;
    const uint8_t* getBoundaryTypeMap() const noexcept;
    const uint16_t* getBoundaryDistanceMap() const noexcept;
    const uint32_t* getBoundarySegmentIdMap() const noexcept;
    const uint32_t* getNearestBoundaryIdMap() const noexcept;
    const tectonics::contract::BoundarySegment* getBoundarySegments() const noexcept;
    uint32_t getBoundarySegmentCount() const noexcept;
    const uint32_t* getDeformingRegionIdMap() const noexcept;
    const uint8_t* getDeformingRegionTypeMap() const noexcept;
    const float* getDeformationRateMap() const noexcept;
    const float* getDeformationVelocityXMap() const noexcept;
    const float* getDeformationVelocityYMap() const noexcept;
    const tectonics::contract::DeformingRegion* getDeformingRegions() const noexcept;
    uint32_t getDeformingRegionCount() const noexcept;
    const tectonics::contract::TectonicJunction* getTectonicJunctions() const noexcept;
    uint32_t getTectonicJunctionCount() const noexcept;
    const uint8_t* getMaterialMap() const noexcept;
    float* getTopography() const noexcept;      ///< Return height map.
    /// Proximity (0..1) to an active mantle plume over the last update; empty when plumes are off.
    const std::vector<float>& getPlumeInfluenceMap() const noexcept { return plume_influence_map; }
    uint32_t* getPlatesMap() const noexcept;    ///< Return a map of the plates owning eaach point
    const uint8_t* getConvergenceMap() const noexcept;
    const uint8_t* getDivergenceMap() const noexcept;
    const uint8_t* getShearMap() const noexcept;
    const uint8_t* getGeologicRegimeMap() const noexcept;
    const tectonics::contract::PlateKinematics* getPlateKinematics() const noexcept;
    uint16_t getSeaLevelMeters() const noexcept {
        return sea_level_m;
    }
    void importNormalizedHeightMap(const float* normalized_map, float sea_level);
    void importRawHeightMap(const float* normalized_map);
    void importMetricHeightMap(const uint16_t* heightmap_m, uint16_t sea_level_m);
    void update();                              ///< Simulate one step of plate tectonics.
    const tectonics::columns::Model& getColumnModel() const { return column_model; }
    /// Per-cell material displacement (cells/update) applied by the latest column update;
    /// empty before the first physical update.
    const std::vector<tectonics::columns::Displacement>& getMaterialDisplacement() const { return material_displacement; }
    const tectonics::topology::Diagnostics& getTopologyDiagnostics() const { return topology_model.diagnostics(); }
    std::vector<double> getRiftDamage() const {
        return topology_model.damageField(imap);
    }
    uint32_t getWidth() const;
    uint32_t getHeight() const;
    bool isFinished() const;
    const tectonics::PlateState& getPlate(uint32_t index) const;

  protected:
  private:
    uint32_t cellCount() const noexcept { return static_cast<uint32_t>(grid_->cellcount()); }
    bool isOceanic(float value) const noexcept;
    void resetSimulationState();
    void initializeHeightMapFromFloat(const float* heightmap);
    void initializeHeightMapFromMetric(const uint16_t* heightmap_m, uint16_t sea_level_m);
    void updateCrustAgeMyrMap();
    void updateDerivedMaps();
    void updateFieldStackMaps();
    void updateBoundaryGraphMaps();
    void updateTectonicJunctions();
    void updateTectonicJunctionCouplingMaps();
    void updateDeformingRegionMaps();
    void updateBoundaryDistanceMap();
    void updateTectonicProvenance();
    void seedInitialTopography(float ocean_coverage, int32_t sea_level_m_override);
    void seedInitialCrustTypes();
    void fillEnclosedInitialOceans();
    void seedMaterialMap();
    void seedInitialCrustAges();
    void solvePlateMotion();
    void initializeColumns();
    void publishColumns();
    void updateColumns();
    void updateTensionMagnitudeMap() noexcept;
    void growPlates(std::vector<std::vector<uint32_t>>& borders);
    void removeEmptyPlates();
    float internalToPhysicalMeters(float internal_height) const noexcept;
    float metersToInternalDelta(float delta_m) const noexcept;

    tectonics::scenario::Scenario scenario;
    std::shared_ptr<const tectonics::grid::TectonicGrid> grid_;
    tectonics::PhaseTimings phase_timings;
    tectonics::dynamics::Diagnostics motion_diagnostics;
    struct MotionSample { uint32_t a=UINT32_MAX,b=UINT32_MAX; double closing=0,ridge=0,slab=0; };
    std::vector<MotionSample> motion_history;
    double motion_history_time = -1;
    tectonics::columns::Model column_model;
    std::vector<tectonics::columns::Displacement> material_displacement;
    tectonics::topology::Model topology_model;
    tectonics::margins::Polarity subduction_polarity;
    std::vector<float> hmap;              ///< Physical topography used by the simulation.
    std::vector<float> prev_hmap;         ///< Previous step's physical topography.
    std::vector<float> display_hmap;      ///< Visual/export topography exposed via the public API.
    std::vector<float> prev_display_hmap; ///< Previous frame's visual/export topography.
    std::vector<float> initial_hmap;      ///< Initial topography before plates begin moving.
    std::vector<uint32_t> imap;           ///< Plate index of the owner of each cell.
    std::vector<uint32_t> prev_imap;      ///< Plate index map from the last update
    std::vector<uint32_t> amap;           ///< Age map of the system's surface (topography).
    // The crust and material the owning plate last recorded at each cell (the torus engine's
    // per-plate rasters, which kept continental material until the cell changed owner): they
    // give the crust thickness split and each plate's continental share.
    std::vector<float> owner_crust_map;
    std::vector<uint8_t> owner_material_map;
    std::vector<tectonics::PlateState> plates; ///< The plates, by live index.
    vector<uint32_t> plate_indices_found; ///< Used in update loop to remove plates

    uint32_t iter_count;       ///< Whole-run iteration count. Used to timestamp new crust.
    uint32_t time_origin_step; ///< Whole-run timestamp corresponding to elapsed time 0.
    uint32_t max_plates;       ///< Number of plates in the initial setting.
    uint32_t num_plates;       ///< Number of plates in the current setting.
    float erosion_strength;    ///< Scales terrain erosion.
    uint16_t sea_level_m;      ///< Coastline threshold in metric export space.
    uint16_t initial_ocean_depth_m; ///< Maximum seeded depth below sea level, in metres.
    uint16_t initial_land_relief_m; ///< Maximum seeded relief above sea level, in metres.
    std::vector<uint8_t> convergence_map;
    std::vector<uint8_t> divergence_map;
    std::vector<uint8_t> shear_map;
    std::vector<uint8_t> geologic_regime_map;
    std::vector<float> crust_age_myr_map;
    std::vector<float> crust_thickness_map;
    std::vector<uint8_t> crust_class_map;
    std::vector<float> uplift_tendency_map;
    std::vector<float> rock_uplift_rate_map;
    std::vector<float> shortening_map;
    std::vector<float> shortening_direction_map;
    std::vector<float> subsidence_tendency_map;
    std::vector<float> accumulated_strain_map;
    std::vector<uint8_t> boundary_type_map;
    std::vector<uint16_t> boundary_distance_map;
    std::vector<uint32_t> boundary_segment_id_map;
    std::vector<uint32_t> nearest_boundary_id_map;
    std::vector<uint32_t> deforming_region_id_map;
    std::vector<uint8_t> deforming_region_type_map;
    std::vector<float> deformation_rate_map;
    std::vector<float> deformation_velocity_x_map;
    std::vector<float> deformation_velocity_y_map;
    std::vector<float> junction_rift_influence_map;
    std::vector<float> junction_transform_influence_map;
    std::vector<float> junction_subduction_influence_map;
    std::vector<uint8_t> surface_crust_type_map;
    std::vector<uint8_t> material_map;
    std::vector<uint8_t> continental_lock_map;
    std::vector<float> rift_tension_x_map;
    std::vector<float> rift_tension_y_map;
    std::vector<float> rift_tension_map;
    std::vector<float> plume_influence_map; // empty unless plumes are on
    std::vector<tectonics::contract::PlateKinematics> plate_kinematics;
    std::vector<tectonics::contract::BoundarySegment> boundary_segments;
    std::vector<tectonics::contract::DeformingRegion> deforming_regions;
    std::vector<tectonics::contract::TectonicJunction> tectonic_junctions;
    uint32_t next_boundary_segment_id = tectonics::contract::kNoBoundaryId + 1U;

    SimpleRandom _randsource;
    int _steps;
    bool finished = false;
};

#endif
