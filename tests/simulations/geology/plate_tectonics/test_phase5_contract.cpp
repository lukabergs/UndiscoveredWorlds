#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_contract.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

platec::scenario::Scenario make_phase5_scenario(long seed) {
    platec::scenario::Scenario scenario;
    scenario.seed = seed;
    scenario.width = 96;
    scenario.height = 48;
    scenario.sea_level = 0.65f;
    scenario.erosion_period = 0;
    scenario.folding_ratio = 0.02f;
    scenario.aggregation_overlap_abs = 1000000;
    scenario.aggregation_overlap_rel = 0.33f;
    scenario.cycle_count = 4;
    scenario.plate_count = 4;
    scenario.cycle_step_limit = 600;
    scenario.delta_time_myr = 1.0;
    return scenario;
}

std::unordered_map<uint32_t, const platec::contract::DeformingRegion*> index_regions_by_id(
    const std::vector<platec::contract::DeformingRegion>& regions) {
    std::unordered_map<uint32_t, const platec::contract::DeformingRegion*> by_id;
    for (const platec::contract::DeformingRegion& region : regions) {
        by_id.emplace(region.id, &region);
    }
    return by_id;
}

std::unordered_map<uint32_t, const platec::contract::BoundarySegment*>
index_boundary_segments_by_id(const std::vector<platec::contract::BoundarySegment>& segments) {
    std::unordered_map<uint32_t, const platec::contract::BoundarySegment*> by_id;
    for (const platec::contract::BoundarySegment& segment : segments) {
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

        const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
        ASSERT_EQ(snapshot.schema_version, platec::contract::kPhase6SchemaVersion);
        ASSERT_EQ(snapshot.deforming_region_id.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deforming_region_type.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_rate.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_velocity_x.size(), snapshot.cell_count());
        ASSERT_EQ(snapshot.deformation_velocity_y.size(), snapshot.cell_count());

        const auto regions_by_id = index_regions_by_id(snapshot.deforming_regions);
        const auto boundary_segments_by_id =
            index_boundary_segments_by_id(snapshot.boundary_segments);

        for (const platec::contract::DeformingRegion& region : snapshot.deforming_regions) {
            EXPECT_NE(region.id, platec::contract::kNoDeformingRegionId);
            EXPECT_NE(region.boundary_segment_id, platec::contract::kNoBoundaryId);
            EXPECT_GT(region.cell_count, 0U);
            EXPECT_GT(region.persistence_steps, 0U);
            EXPECT_GT(region.average_deformation_rate, 0.0f);
            EXPECT_LE(region.average_deformation_rate, 1.0f);
            EXPECT_TRUE(std::isfinite(region.average_interpolated_velocity_x));
            EXPECT_TRUE(std::isfinite(region.average_interpolated_velocity_y));
            EXPECT_TRUE(boundary_segments_by_id.find(region.boundary_segment_id) !=
                        boundary_segments_by_id.end());

            if (region.type == platec::contract::DeformingRegionType::ContinentalRift) {
                ++rift_regions;
                EXPECT_EQ(boundary_segments_by_id.at(region.boundary_segment_id)->geologic_regime,
                          platec::contract::GeologicRegime::DivergentRift);
            } else if (region.type ==
                       platec::contract::DeformingRegionType::DiffuseCollision) {
                ++collision_regions;
                EXPECT_EQ(boundary_segments_by_id.at(region.boundary_segment_id)->geologic_regime,
                          platec::contract::GeologicRegime::ContinentCollision);
            }
        }

        for (size_t i = 0; i < snapshot.cell_count(); ++i) {
            const uint32_t region_id = snapshot.deforming_region_id[i];
            const auto region_type = static_cast<platec::contract::DeformingRegionType>(
                snapshot.deforming_region_type[i]);
            EXPECT_GE(snapshot.deformation_rate[i], 0.0f);
            EXPECT_LE(snapshot.deformation_rate[i], 1.0f);
            EXPECT_TRUE(std::isfinite(snapshot.deformation_velocity_x[i]));
            EXPECT_TRUE(std::isfinite(snapshot.deformation_velocity_y[i]));

            if (region_id == platec::contract::kNoDeformingRegionId) {
                EXPECT_EQ(region_type, platec::contract::DeformingRegionType::None);
                continue;
            }

            EXPECT_NE(region_type, platec::contract::DeformingRegionType::None);
            EXPECT_GT(snapshot.deformation_rate[i], 0.0f);
            EXPECT_TRUE(regions_by_id.find(region_id) != regions_by_id.end());

            if (region_type == platec::contract::DeformingRegionType::ContinentalRift) {
                ++rift_cells;
                rift_uplift_total += snapshot.uplift_tendency[i];
                rift_subsidence_total += snapshot.subsidence_tendency[i];
            } else if (region_type ==
                       platec::contract::DeformingRegionType::DiffuseCollision) {
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
    const platec::scenario::Scenario scenario = make_phase5_scenario(12345);
    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    while (platec_api_is_finished(simulation) == 0U) {
        platec_api_step(simulation);
    }

    ASSERT_NE(platec_api_get_deforming_region_id_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_deforming_region_type_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_deformation_rate_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_deformation_velocity_x_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_deformation_velocity_y_map(simulation), nullptr);

    const uint32_t region_count = platec_api_get_deforming_region_count(simulation);
    ASSERT_GT(region_count, 0U);
    ASSERT_NE(platec_api_get_deforming_regions(simulation), nullptr);

    const auto* region_id = platec_api_get_deforming_region_id_map(simulation);
    const auto* region_type = platec_api_get_deforming_region_type_map(simulation);
    const auto* deformation_rate = platec_api_get_deformation_rate_map(simulation);
    const auto* velocity_x = platec_api_get_deformation_velocity_x_map(simulation);
    const auto* velocity_y = platec_api_get_deformation_velocity_y_map(simulation);
    const auto* regions = platec_api_get_deforming_regions(simulation);

    std::unordered_set<uint32_t> ids;
    for (uint32_t i = 0; i < region_count; ++i) {
        EXPECT_NE(regions[i].id, platec::contract::kNoDeformingRegionId);
        EXPECT_NE(regions[i].type, platec::contract::DeformingRegionType::None);
        EXPECT_GT(regions[i].cell_count, 0U);
        EXPECT_GT(regions[i].persistence_steps, 0U);
        EXPECT_GT(regions[i].age_myr, 0.0);
        ids.insert(regions[i].id);
    }

    const size_t cell_count =
        static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    size_t active_cells = 0U;
    for (size_t i = 0; i < cell_count; ++i) {
        EXPECT_GE(deformation_rate[i], 0.0f);
        EXPECT_LE(deformation_rate[i], 1.0f);
        EXPECT_TRUE(std::isfinite(velocity_x[i]));
        EXPECT_TRUE(std::isfinite(velocity_y[i]));
        if (region_id[i] == platec::contract::kNoDeformingRegionId) {
            EXPECT_EQ(region_type[i],
                      static_cast<uint8_t>(platec::contract::DeformingRegionType::None));
            continue;
        }

        ++active_cells;
        EXPECT_TRUE(ids.find(region_id[i]) != ids.end());
        EXPECT_NE(region_type[i],
                  static_cast<uint8_t>(platec::contract::DeformingRegionType::None));
        EXPECT_GT(deformation_rate[i], 0.0f);
    }

    EXPECT_GT(active_cells, 0U);

    platec_api_destroy(simulation);
}
