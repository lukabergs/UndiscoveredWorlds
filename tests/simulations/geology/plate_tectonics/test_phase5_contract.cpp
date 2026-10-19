#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/lithosphere.hpp"
#include "simulations/geology/plate_tectonics/tectonic_api.hpp"
#include "simulations/geology/plate_tectonics/tectonic_contract.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

tectonics::scenario::Scenario make_phase5_scenario(long seed) {
    tectonics::scenario::Scenario scenario;
    scenario.seed = seed;
    scenario.sphere_face_cells = 24;  // 6 x 24^2 cells of about 417 km
    scenario.sea_level = 0.65f;
    scenario.plate_count = 4;
    scenario.delta_time_myr = 1.0;
    return scenario;
}

std::unordered_map<uint32_t, const tectonics::contract::DeformingRegion*> index_regions_by_id(
    const std::vector<tectonics::contract::DeformingRegion>& regions) {
    std::unordered_map<uint32_t, const tectonics::contract::DeformingRegion*> by_id;
    for (const tectonics::contract::DeformingRegion& region : regions) {
        by_id.emplace(region.id, &region);
    }
    return by_id;
}

std::unordered_map<uint32_t, const tectonics::contract::BoundarySegment*>
index_boundary_segments_by_id(const std::vector<tectonics::contract::BoundarySegment>& segments) {
    std::unordered_map<uint32_t, const tectonics::contract::BoundarySegment*> by_id;
    for (const tectonics::contract::BoundarySegment& segment : segments) {
        by_id.emplace(segment.id, &segment);
    }
    return by_id;
}

} // namespace

TEST(Phase5Contract, SnapshotExposesLimitedDeformingRegions) {
    const std::vector<long> seeds = {12345, 67890, 24680, 97531};

    size_t rift_regions = 0U;
    size_t collision_regions = 0U;
    size_t rift_cells = 0U;
    size_t collision_cells = 0U;
    double rift_uplift_total = 0.0;
    double rift_subsidence_total = 0.0;
    double collision_uplift_total = 0.0;
    double collision_subsidence_total = 0.0;

    for (long seed : seeds) {
        lithosphere simulation(make_phase5_scenario(seed));
        for (uint32_t i = 0; i < 6; ++i) {
            ASSERT_FALSE(simulation.isFinished());
            simulation.update();
        }

        const tectonics::contract::Snapshot snapshot = tectonics::contract::capture_snapshot(simulation);
        ASSERT_EQ(snapshot.schema_version, tectonics::contract::kSchemaVersion);
        ASSERT_EQ(snapshot.deforming_region_id.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deforming_region_type.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_rate.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_velocity_x.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_velocity_y.size(), snapshot.cell_count());

        const auto regions_by_id = index_regions_by_id(snapshot.deforming_regions);
        const auto boundary_segments_by_id =
            index_boundary_segments_by_id(snapshot.boundary_segments);

        for (const tectonics::contract::DeformingRegion& region : snapshot.deforming_regions) {
            EXPECT_NE(region.id, tectonics::contract::kNoDeformingRegionId);
            EXPECT_NE(region.boundary_segment_id, tectonics::contract::kNoBoundaryId);
            EXPECT_GT(region.cell_count, 0U);
            EXPECT_GT(region.persistence_steps, 0U);
            EXPECT_GT(region.average_deformation_rate, 0.0f);
            EXPECT_LE(region.average_deformation_rate, 1.0f);
            EXPECT_TRUE(std::isfinite(region.average_interpolated_velocity_x));
            EXPECT_TRUE(std::isfinite(region.average_interpolated_velocity_y));
            EXPECT_TRUE(boundary_segments_by_id.find(region.boundary_segment_id) !=
                        boundary_segments_by_id.end());

            if (region.type == tectonics::contract::DeformingRegionType::ContinentalRift) {
                ++rift_regions;
                EXPECT_EQ(boundary_segments_by_id.at(region.boundary_segment_id)->geologic_regime,
                          tectonics::contract::GeologicRegime::DivergentRift);
            } else if (region.type ==
                       tectonics::contract::DeformingRegionType::DiffuseCollision) {
                ++collision_regions;
                EXPECT_EQ(boundary_segments_by_id.at(region.boundary_segment_id)->geologic_regime,
                          tectonics::contract::GeologicRegime::ContinentCollision);
            }
        }

        for (size_t i = 0; i < snapshot.cell_count(); ++i) {
            const uint32_t region_id = snapshot.deforming_region_id[i];
            const auto region_type = static_cast<tectonics::contract::DeformingRegionType>(
                snapshot.deforming_region_type[i]);
            EXPECT_GE(snapshot.deformation_rate[i], 0.0f);
            EXPECT_LE(snapshot.deformation_rate[i], 1.0f);
            EXPECT_TRUE(std::isfinite(snapshot.deformation_velocity_x[i]));
            EXPECT_TRUE(std::isfinite(snapshot.deformation_velocity_y[i]));

            if (region_id == tectonics::contract::kNoDeformingRegionId) {
                EXPECT_EQ(region_type, tectonics::contract::DeformingRegionType::None);
                continue;
            }

            EXPECT_NE(region_type, tectonics::contract::DeformingRegionType::None);
            EXPECT_GT(snapshot.deformation_rate[i], 0.0f);
            EXPECT_TRUE(regions_by_id.find(region_id) != regions_by_id.end());

            if (region_type == tectonics::contract::DeformingRegionType::ContinentalRift) {
                ++rift_cells;
                rift_uplift_total += snapshot.uplift_tendency[i];
                rift_subsidence_total += snapshot.subsidence_tendency[i];
            } else if (region_type ==
                       tectonics::contract::DeformingRegionType::DiffuseCollision) {
                ++collision_cells;
                collision_uplift_total += snapshot.uplift_tendency[i];
                collision_subsidence_total += snapshot.subsidence_tendency[i];
            }
        }
    }

    EXPECT_GT(rift_regions, 0U);
    EXPECT_GT(collision_regions, 0U);
    EXPECT_GT(rift_cells, 0U);
    EXPECT_GT(collision_cells, 0U);
    EXPECT_GT(rift_subsidence_total / static_cast<double>(rift_cells),
              rift_uplift_total / static_cast<double>(rift_cells));
    EXPECT_GT(collision_uplift_total / static_cast<double>(collision_cells),
              collision_subsidence_total / static_cast<double>(collision_cells));
}

TEST(Phase5Contract, CApiDeformingRegionsPersistAfterFinishedState) {
    const tectonics::scenario::Scenario scenario = make_phase5_scenario(12345);
    void* simulation = tectonics_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    while (tectonics_is_finished(simulation) == 0U) {
        tectonics_step(simulation);
    }

    ASSERT_NE(tectonics_get_deforming_region_id_map(simulation), nullptr);
    ASSERT_NE(tectonics_get_deforming_region_type_map(simulation), nullptr);
    ASSERT_NE(tectonics_get_deformation_rate_map(simulation), nullptr);
    ASSERT_NE(tectonics_get_deformation_velocity_x_map(simulation), nullptr);
    ASSERT_NE(tectonics_get_deformation_velocity_y_map(simulation), nullptr);

    const uint32_t region_count = tectonics_get_deforming_region_count(simulation);
    ASSERT_GT(region_count, 0U);
    ASSERT_NE(tectonics_get_deforming_regions(simulation), nullptr);

    const auto* region_id = tectonics_get_deforming_region_id_map(simulation);
    const auto* region_type = tectonics_get_deforming_region_type_map(simulation);
    const auto* deformation_rate = tectonics_get_deformation_rate_map(simulation);
    const auto* velocity_x = tectonics_get_deformation_velocity_x_map(simulation);
    const auto* velocity_y = tectonics_get_deformation_velocity_y_map(simulation);
    const auto* regions = tectonics_get_deforming_regions(simulation);

    std::unordered_set<uint32_t> ids;
    for (uint32_t i = 0; i < region_count; ++i) {
        EXPECT_NE(regions[i].id, tectonics::contract::kNoDeformingRegionId);
        EXPECT_NE(regions[i].type, tectonics::contract::DeformingRegionType::None);
        EXPECT_GT(regions[i].cell_count, 0U);
        EXPECT_GT(regions[i].persistence_steps, 0U);
        EXPECT_GT(regions[i].age_myr, 0.0);
        ids.insert(regions[i].id);
    }

    const size_t cell_count =
        6 * static_cast<size_t>(scenario.sphere_face_cells) * scenario.sphere_face_cells;
    size_t active_cells = 0U;
    for (size_t i = 0; i < cell_count; ++i) {
        EXPECT_GE(deformation_rate[i], 0.0f);
        EXPECT_LE(deformation_rate[i], 1.0f);
        EXPECT_TRUE(std::isfinite(velocity_x[i]));
        EXPECT_TRUE(std::isfinite(velocity_y[i]));
        if (region_id[i] == tectonics::contract::kNoDeformingRegionId) {
            EXPECT_EQ(region_type[i],
                      static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None));
            continue;
        }

        ++active_cells;
        EXPECT_TRUE(ids.find(region_id[i]) != ids.end());
        EXPECT_NE(region_type[i],
                  static_cast<uint8_t>(tectonics::contract::DeformingRegionType::None));
        EXPECT_GT(deformation_rate[i], 0.0f);
    }

    EXPECT_GT(active_cells, 0U);

    tectonics_destroy(simulation);
}
