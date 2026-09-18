#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_contract.hpp"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

namespace {

platec::scenario::Scenario make_phase4_scenario() {
    platec::scenario::Scenario scenario;
    scenario.seed = 97531;
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

std::unordered_map<uint32_t, const platec::contract::BoundarySegment*>
index_segments_by_id(const std::vector<platec::contract::BoundarySegment>& segments) {
    std::unordered_map<uint32_t, const platec::contract::BoundarySegment*> by_id;
    for (const platec::contract::BoundarySegment& segment : segments) {
        by_id.emplace(segment.id, &segment);
    }
    return by_id;
}

} // namespace

TEST(Phase4Contract, SnapshotExposesBoundaryGraphAndPersistentIds) {
    lithosphere simulation(make_phase4_scenario());

    for (uint32_t i = 0; i < 5; ++i) {
        ASSERT_FALSE(simulation.isFinished());
        simulation.update();
    }

    const platec::contract::Snapshot before = platec::contract::capture_snapshot(simulation);
    ASSERT_EQ(before.schema_version, platec::contract::kPhase6SchemaVersion);
    ASSERT_EQ(before.boundary_segment_id.size(), before.cell_count());
    ASSERT_EQ(before.nearest_boundary_id.size(), before.cell_count());
    ASSERT_FALSE(before.boundary_segments.empty());

    std::unordered_set<uint32_t> ids_before;
    for (const platec::contract::BoundarySegment& segment : before.boundary_segments) {
        EXPECT_NE(segment.id, platec::contract::kNoBoundaryId);
        EXPECT_NE(segment.boundary_type, platec::contract::BoundaryType::None);
        EXPECT_GT(segment.cell_count, 0U);
        EXPECT_GT(segment.length_cells, 0.0f);
        EXPECT_GT(segment.persistence_steps, 0U);
        EXPECT_GT(segment.age_myr, 0.0);
        ids_before.insert(segment.id);
    }

    simulation.update();
    const platec::contract::Snapshot after = platec::contract::capture_snapshot(simulation);
    const auto segments_by_id = index_segments_by_id(after.boundary_segments);

    size_t boundary_cells = 0U;
    size_t interior_cells = 0U;
    size_t persisted_segments = 0U;

    for (const platec::contract::BoundarySegment& segment : after.boundary_segments) {
        if (ids_before.find(segment.id) != ids_before.end() && segment.persistence_steps > 1U) {
            ++persisted_segments;
        }
    }

    for (size_t i = 0; i < after.cell_count(); ++i) {
        const uint32_t boundary_id = after.boundary_segment_id[i];
        const uint32_t nearest_id = after.nearest_boundary_id[i];
        if (boundary_id != platec::contract::kNoBoundaryId) {
            ++boundary_cells;
            EXPECT_EQ(after.boundary_distance[i], 0U);
            EXPECT_EQ(nearest_id, boundary_id);
            EXPECT_NE(after.boundary_type[i],
                      static_cast<uint8_t>(platec::contract::BoundaryType::None));
            EXPECT_TRUE(segments_by_id.find(boundary_id) != segments_by_id.end());
        } else if (nearest_id != platec::contract::kNoBoundaryId) {
            ++interior_cells;
            EXPECT_GT(after.boundary_distance[i], 0U);
            EXPECT_TRUE(segments_by_id.find(nearest_id) != segments_by_id.end());
        }
    }

    EXPECT_GT(boundary_cells, 0U);
    EXPECT_GT(interior_cells, 0U);
    EXPECT_GT(persisted_segments, 0U);
}

TEST(Phase4Contract, CApiBoundaryGraphPersistsAfterFinishedState) {
    const platec::scenario::Scenario scenario = make_phase4_scenario();
    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    while (platec_api_is_finished(simulation) == 0U) {
        platec_api_step(simulation);
    }

    ASSERT_NE(platec_api_get_boundary_segment_id_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_nearest_boundary_id_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_boundary_type_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_boundary_distance_map(simulation), nullptr);

    const uint32_t segment_count = platec_api_get_boundary_segment_count(simulation);
    ASSERT_GT(segment_count, 0U);
    ASSERT_NE(platec_api_get_boundary_segments(simulation), nullptr);

    const auto* segment_id = platec_api_get_boundary_segment_id_map(simulation);
    const auto* nearest_id = platec_api_get_nearest_boundary_id_map(simulation);
    const auto* boundary_type = platec_api_get_boundary_type_map(simulation);
    const auto* boundary_distance = platec_api_get_boundary_distance_map(simulation);
    const auto* segments = platec_api_get_boundary_segments(simulation);

    std::unordered_set<uint32_t> ids;
    for (uint32_t i = 0; i < segment_count; ++i) {
        EXPECT_NE(segments[i].id, platec::contract::kNoBoundaryId);
        EXPECT_NE(segments[i].boundary_type, platec::contract::BoundaryType::None);
        EXPECT_GT(segments[i].cell_count, 0U);
        EXPECT_GT(segments[i].persistence_steps, 0U);
        EXPECT_GT(segments[i].age_myr, 0.0);
        ids.insert(segments[i].id);
    }

    const size_t cell_count =
        static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    size_t boundary_cells = 0U;
    size_t nearest_cells = 0U;
    for (size_t i = 0; i < cell_count; ++i) {
        boundary_cells += segment_id[i] != platec::contract::kNoBoundaryId ? 1U : 0U;
        nearest_cells += nearest_id[i] != platec::contract::kNoBoundaryId ? 1U : 0U;
        if (segment_id[i] != platec::contract::kNoBoundaryId) {
            EXPECT_EQ(boundary_distance[i], 0U);
            EXPECT_EQ(nearest_id[i], segment_id[i]);
            EXPECT_NE(boundary_type[i],
                      static_cast<uint8_t>(platec::contract::BoundaryType::None));
            EXPECT_TRUE(ids.find(segment_id[i]) != ids.end());
        } else if (nearest_id[i] != platec::contract::kNoBoundaryId) {
            EXPECT_GT(boundary_distance[i], 0U);
            EXPECT_TRUE(ids.find(nearest_id[i]) != ids.end());
        }
    }

    EXPECT_GT(boundary_cells, 0U);
    EXPECT_GT(nearest_cells, boundary_cells);

    platec_api_destroy(simulation);
}
