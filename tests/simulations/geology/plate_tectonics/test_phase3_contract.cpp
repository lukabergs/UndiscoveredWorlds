#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_contract.hpp"

#include <cstddef>
#include <cstdint>

namespace {

platec::scenario::Scenario make_phase3_scenario() {
    platec::scenario::Scenario scenario;
    scenario.seed = 13579;
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

} // namespace

TEST(Phase3Contract, SnapshotExposesNativeFieldStackMaps) {
    lithosphere simulation(make_phase3_scenario());

    for (uint32_t i = 0; i < 5; ++i) {
        ASSERT_FALSE(simulation.isFinished());
        simulation.update();
    }

    const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
    ASSERT_EQ(snapshot.schema_version, platec::contract::kPhase6SchemaVersion);
    ASSERT_EQ(snapshot.crust_thickness.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.crust_class.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.uplift_tendency.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.subsidence_tendency.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.accumulated_strain.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.boundary_type.size(), snapshot.cell_count());
    ASSERT_EQ(snapshot.boundary_distance.size(), snapshot.cell_count());

    size_t boundary_cells = 0;
    size_t interior_cells = 0;
    size_t oceanic_cells = 0;
    size_t continental_cells = 0;
    size_t strained_cells = 0;

    for (size_t i = 0; i < snapshot.cell_count(); ++i) {
        EXPECT_GE(snapshot.crust_thickness[i], snapshot.heightmap[i] - 1e-5f);
        EXPECT_GE(snapshot.uplift_tendency[i], 0.0f);
        EXPECT_LE(snapshot.uplift_tendency[i], 1.0f);
        EXPECT_GE(snapshot.subsidence_tendency[i], 0.0f);
        EXPECT_LE(snapshot.subsidence_tendency[i], 1.0f);
        EXPECT_GE(snapshot.accumulated_strain[i], 0.0f);
        EXPECT_LE(snapshot.accumulated_strain[i], 1.0f);

        if (snapshot.boundary_type[i] !=
            static_cast<uint8_t>(platec::contract::BoundaryType::None)) {
            ++boundary_cells;
            EXPECT_EQ(snapshot.boundary_distance[i], 0U);
        } else if (snapshot.boundary_distance[i] > 0U) {
            ++interior_cells;
        }

        if (snapshot.crust_class[i] ==
            static_cast<uint8_t>(platec::contract::CrustClass::Oceanic)) {
            ++oceanic_cells;
        } else if (snapshot.crust_class[i] ==
                   static_cast<uint8_t>(platec::contract::CrustClass::Continental)) {
            ++continental_cells;
        }

        strained_cells += snapshot.accumulated_strain[i] > 0.01f ? 1U : 0U;
    }

    EXPECT_GT(boundary_cells, 0U);
    EXPECT_GT(interior_cells, 0U);
    EXPECT_GT(oceanic_cells, 0U);
    EXPECT_GT(continental_cells, 0U);
    EXPECT_GT(strained_cells, 0U);
}

TEST(Phase3Contract, CApiFieldStackPersistsAfterFinishedState) {
    const platec::scenario::Scenario scenario = make_phase3_scenario();
    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    while (platec_api_is_finished(simulation) == 0U) {
        platec_api_step(simulation);
    }

    ASSERT_NE(platec_api_get_crust_thickness_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_crust_class_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_uplift_tendency_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_subsidence_tendency_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_accumulated_strain_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_boundary_type_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_boundary_distance_map(simulation), nullptr);

    const size_t cell_count =
        static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    const float* thickness = platec_api_get_crust_thickness_map(simulation);
    const float* uplift = platec_api_get_uplift_tendency_map(simulation);
    const float* subsidence = platec_api_get_subsidence_tendency_map(simulation);
    const float* strain = platec_api_get_accumulated_strain_map(simulation);
    const uint8_t* crust_class = platec_api_get_crust_class_map(simulation);
    const uint8_t* boundary_type = platec_api_get_boundary_type_map(simulation);
    const uint16_t* boundary_distance = platec_api_get_boundary_distance_map(simulation);

    size_t nonzero_thickness = 0;
    size_t zero_distance_boundaries = 0;
    size_t non_boundary_distance = 0;
    size_t classified_cells = 0;
    size_t active_tendency_cells = 0;

    for (size_t i = 0; i < cell_count; ++i) {
        nonzero_thickness += thickness[i] > 0.0f ? 1U : 0U;
        zero_distance_boundaries +=
            (boundary_type[i] != static_cast<uint8_t>(platec::contract::BoundaryType::None) &&
             boundary_distance[i] == 0U)
                ? 1U
                : 0U;
        non_boundary_distance +=
            (boundary_type[i] == static_cast<uint8_t>(platec::contract::BoundaryType::None) &&
             boundary_distance[i] > 0U)
                ? 1U
                : 0U;
        classified_cells +=
            crust_class[i] != static_cast<uint8_t>(platec::contract::CrustClass::None) ? 1U : 0U;
        active_tendency_cells +=
            (uplift[i] > 0.01f || subsidence[i] > 0.01f || strain[i] > 0.01f) ? 1U : 0U;
    }

    EXPECT_GT(nonzero_thickness, 0U);
    EXPECT_GT(zero_distance_boundaries, 0U);
    EXPECT_GT(non_boundary_distance, 0U);
    EXPECT_GT(classified_cells, 0U);
    EXPECT_GT(active_tendency_cells, 0U);

    platec_api_destroy(simulation);
}
