#include "tectonic_contract.hpp"

#include "lithosphere.hpp"
#include "tectonic_provenance.hpp"

#include <cstring>

namespace platec::contract {

Snapshot capture_snapshot(const lithosphere& litho) {
    Snapshot snapshot;
    snapshot.width = litho.getWidth();
    snapshot.height = litho.getHeight();
    snapshot.iteration_count = litho.getIterationCount();
    snapshot.cycle_count = litho.getCycleCount();
    snapshot.time_origin_step = litho.getTimeOriginStep();
    snapshot.sea_level_m = litho.getSeaLevelMeters();
    snapshot.time_myr = litho.getTimeMyr();
    snapshot.delta_time_myr = litho.getDeltaTimeMyr();
    snapshot.run_scenario = litho.getScenario();

    const size_t cell_count = snapshot.cell_count();
    snapshot.heightmap.resize(cell_count);
    snapshot.plate_id.resize(cell_count);
    snapshot.crust_age_steps.resize(cell_count);
    snapshot.crust_age_myr.resize(cell_count);
    snapshot.crust_thickness.resize(cell_count);
    snapshot.crust_class.resize(cell_count);
    snapshot.uplift_tendency.resize(cell_count);
    snapshot.subsidence_tendency.resize(cell_count);
    snapshot.accumulated_strain.resize(cell_count);
    snapshot.boundary_type.resize(cell_count);
    snapshot.boundary_distance.resize(cell_count);
    snapshot.boundary_segment_id.resize(cell_count);
    snapshot.nearest_boundary_id.resize(cell_count);
    snapshot.deforming_region_id.resize(cell_count);
    snapshot.deforming_region_type.resize(cell_count);
    snapshot.deformation_rate.resize(cell_count);
    snapshot.deformation_velocity_x.resize(cell_count);
    snapshot.deformation_velocity_y.resize(cell_count);
    snapshot.convergence_score.resize(cell_count);
    snapshot.divergence_score.resize(cell_count);
    snapshot.shear_score.resize(cell_count);
    snapshot.geologic_regime.resize(cell_count);

    std::memcpy(snapshot.heightmap.data(), litho.getTopography(), sizeof(float) * cell_count);
    std::memcpy(snapshot.plate_id.data(), litho.getPlatesMap(), sizeof(uint32_t) * cell_count);
    std::memcpy(snapshot.crust_age_steps.data(), litho.getAgeMap(),
                sizeof(uint32_t) * cell_count);
    std::memcpy(snapshot.crust_age_myr.data(), litho.getCrustAgeMyrMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.crust_thickness.data(), litho.getCrustThicknessMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.crust_class.data(), litho.getCrustClassMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.uplift_tendency.data(), litho.getUpliftTendencyMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.subsidence_tendency.data(), litho.getSubsidenceTendencyMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.accumulated_strain.data(), litho.getAccumulatedStrainMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.boundary_type.data(), litho.getBoundaryTypeMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.boundary_distance.data(), litho.getBoundaryDistanceMap(),
                sizeof(uint16_t) * cell_count);
    std::memcpy(snapshot.boundary_segment_id.data(), litho.getBoundarySegmentIdMap(),
                sizeof(uint32_t) * cell_count);
    std::memcpy(snapshot.nearest_boundary_id.data(), litho.getNearestBoundaryIdMap(),
                sizeof(uint32_t) * cell_count);
    std::memcpy(snapshot.deforming_region_id.data(), litho.getDeformingRegionIdMap(),
                sizeof(uint32_t) * cell_count);
    std::memcpy(snapshot.deforming_region_type.data(), litho.getDeformingRegionTypeMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.deformation_rate.data(), litho.getDeformationRateMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.deformation_velocity_x.data(), litho.getDeformationVelocityXMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.deformation_velocity_y.data(), litho.getDeformationVelocityYMap(),
                sizeof(float) * cell_count);
    std::memcpy(snapshot.convergence_score.data(), litho.getConvergenceMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.divergence_score.data(), litho.getDivergenceMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.shear_score.data(), litho.getShearMap(),
                sizeof(uint8_t) * cell_count);
    std::memcpy(snapshot.geologic_regime.data(), litho.getGeologicRegimeMap(),
                sizeof(uint8_t) * cell_count);

    snapshot.plates.resize(litho.getProvenancePlateCount());
    if (!snapshot.plates.empty()) {
        std::memcpy(snapshot.plates.data(), litho.getPlateKinematics(),
                    sizeof(PlateKinematics) * snapshot.plates.size());
    }

    snapshot.boundary_segments.resize(litho.getBoundarySegmentCount());
    if (!snapshot.boundary_segments.empty()) {
        std::memcpy(snapshot.boundary_segments.data(), litho.getBoundarySegments(),
                    sizeof(BoundarySegment) * snapshot.boundary_segments.size());
    }

    snapshot.deforming_regions.resize(litho.getDeformingRegionCount());
    if (!snapshot.deforming_regions.empty()) {
        std::memcpy(snapshot.deforming_regions.data(), litho.getDeformingRegions(),
                    sizeof(DeformingRegion) * snapshot.deforming_regions.size());
    }

    snapshot.tectonic_junctions.resize(litho.getTectonicJunctionCount());
    if (!snapshot.tectonic_junctions.empty()) {
        std::memcpy(snapshot.tectonic_junctions.data(), litho.getTectonicJunctions(),
                    sizeof(TectonicJunction) * snapshot.tectonic_junctions.size());
    }

    return snapshot;
}

Snapshot capture_phase0_snapshot(const lithosphere& lithosphere) {
    Snapshot snapshot = capture_snapshot(lithosphere);
    snapshot.schema_version = kPhase0SchemaVersion;
    if (!snapshot.heightmap.empty() && !snapshot.plate_id.empty() && !snapshot.plates.empty()) {
        const platec::provenance::Inputs inputs{
            snapshot.width,
            snapshot.height,
            snapshot.heightmap.data(),
            snapshot.plate_id.data(),
            snapshot.plates.data(),
            static_cast<uint32_t>(snapshot.plates.size()),
        };
        platec::provenance::compute_maps(inputs, snapshot.convergence_score.data(),
                                         snapshot.divergence_score.data(),
                                         snapshot.shear_score.data(),
                                         snapshot.geologic_regime.data());
    }
    return snapshot;
}

const char* geologic_regime_name(GeologicRegime regime) {
    switch (regime) {
    case GeologicRegime::Stable:
        return "stable";
    case GeologicRegime::ConvergentArc:
        return "convergent_arc";
    case GeologicRegime::ContinentCollision:
        return "continent_collision";
    case GeologicRegime::DivergentRift:
        return "divergent_rift";
    case GeologicRegime::Transform:
        return "transform";
    case GeologicRegime::PassiveMargin:
        return "passive_margin";
    case GeologicRegime::MidOceanRidge:
        return "mid_ocean_ridge";
    case GeologicRegime::TrenchAdjacent:
        return "trench_adjacent";
    default:
        return "unknown";
    }
}

const char* crust_class_name(CrustClass crust_class) {
    switch (crust_class) {
    case CrustClass::None:
        return "none";
    case CrustClass::Oceanic:
        return "oceanic";
    case CrustClass::Transitional:
        return "transitional";
    case CrustClass::Continental:
        return "continental";
    default:
        return "unknown";
    }
}

const char* boundary_type_name(BoundaryType boundary_type) {
    switch (boundary_type) {
    case BoundaryType::None:
        return "none";
    case BoundaryType::Convergent:
        return "convergent";
    case BoundaryType::Divergent:
        return "divergent";
    case BoundaryType::Transform:
        return "transform";
    case BoundaryType::PassiveMargin:
        return "passive_margin";
    default:
        return "unknown";
    }
}

const char* deforming_region_type_name(DeformingRegionType region_type) {
    switch (region_type) {
    case DeformingRegionType::None:
        return "none";
    case DeformingRegionType::ContinentalRift:
        return "continental_rift";
    case DeformingRegionType::DiffuseCollision:
        return "diffuse_collision";
    default:
        return "unknown";
    }
}

const char* junction_arm_type_name(JunctionArmType arm_type) {
    switch (arm_type) {
    case JunctionArmType::Unknown:
        return "unknown";
    case JunctionArmType::Ridge:
        return "ridge";
    case JunctionArmType::Trench:
        return "trench";
    case JunctionArmType::Fault:
        return "fault";
    case JunctionArmType::PassiveMargin:
        return "passive_margin";
    default:
        return "unknown";
    }
}

char junction_arm_type_code(JunctionArmType arm_type) {
    switch (arm_type) {
    case JunctionArmType::Ridge:
        return 'R';
    case JunctionArmType::Trench:
        return 'T';
    case JunctionArmType::Fault:
        return 'F';
    case JunctionArmType::PassiveMargin:
        return 'P';
    case JunctionArmType::Unknown:
    default:
        return '?';
    }
}

const char* junction_stability_name(JunctionStability stability) {
    switch (stability) {
    case JunctionStability::Unknown:
        return "unknown";
    case JunctionStability::Stable:
        return "stable";
    case JunctionStability::ConditionallyStable:
        return "conditionally_stable";
    case JunctionStability::Unstable:
        return "unstable";
    case JunctionStability::Indeterminate:
        return "indeterminate";
    default:
        return "unknown";
    }
}

} // namespace platec::contract
