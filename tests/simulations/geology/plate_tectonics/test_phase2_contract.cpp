#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_contract.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>

namespace {

platec::scenario::Scenario make_phase2_scenario(double delta_time_myr) {
    platec::scenario::Scenario scenario;
    scenario.seed = 24680;
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
    scenario.delta_time_myr = delta_time_myr;
    return scenario;
}

} // namespace

TEST(Phase2Contract, SnapshotExposesScenarioAndExplicitTime) {
    const platec::scenario::Scenario scenario = make_phase2_scenario(2.5);
    lithosphere simulation(scenario);

    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 0.0);
    EXPECT_EQ(simulation.getTimeOriginStep(), simulation.getIterationCount());
    EXPECT_DOUBLE_EQ(simulation.getDeltaTimeMyr(), scenario.delta_time_myr);

    for (uint32_t i = 0; i < 3; ++i) {
        ASSERT_FALSE(simulation.isFinished());
        simulation.update();
    }

    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 7.5);

    const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
    EXPECT_EQ(snapshot.schema_version, platec::contract::kPhase6SchemaVersion);
    EXPECT_EQ(snapshot.time_origin_step, simulation.getTimeOriginStep());
    EXPECT_DOUBLE_EQ(snapshot.time_myr, simulation.getTimeMyr());
    EXPECT_DOUBLE_EQ(snapshot.delta_time_myr, scenario.delta_time_myr);
    EXPECT_EQ(snapshot.run_scenario.seed, scenario.seed);
    EXPECT_EQ(snapshot.run_scenario.width, scenario.width);
    EXPECT_EQ(snapshot.run_scenario.height, scenario.height);
    EXPECT_EQ(snapshot.run_scenario.plate_count, scenario.plate_count);
    ASSERT_EQ(snapshot.crust_age_myr.size(), snapshot.cell_count());

    for (size_t i = 0; i < snapshot.cell_count(); i += 257) {
        const uint32_t age_steps =
            snapshot.iteration_count >= snapshot.crust_age_steps[i]
                ? snapshot.iteration_count - snapshot.crust_age_steps[i]
                : 0U;
        const float expected =
            static_cast<float>(platec::scenario::steps_to_myr(age_steps, scenario.delta_time_myr));
        EXPECT_FLOAT_EQ(snapshot.crust_age_myr[i], expected);
    }
}

TEST(Phase2Contract, CApiExposesScenarioTimeAndAgeMaps) {
    const platec::scenario::Scenario scenario = make_phase2_scenario(1.25);
    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);

    EXPECT_DOUBLE_EQ(platec_api_get_time_myr(simulation), 0.0);
    EXPECT_DOUBLE_EQ(platec_api_get_delta_time_myr(simulation), scenario.delta_time_myr);
    const uint32_t origin_step = platec_api_get_time_origin_step(simulation);
    EXPECT_EQ(origin_step, platec_api_get_iteration_count(simulation));

    for (uint32_t i = 0; i < 2; ++i) {
        ASSERT_EQ(platec_api_is_finished(simulation), 0U);
        platec_api_step(simulation);
    }

    EXPECT_DOUBLE_EQ(platec_api_get_time_myr(simulation), 2.5);
    EXPECT_EQ(platec_api_get_iteration_count(simulation) - origin_step, 2U);
    EXPECT_EQ(platec_api_get_cycle_count(simulation), 0U);
    ASSERT_NE(platec_api_get_crust_age_myr_map(simulation), nullptr);

    const float* crust_age_myr = platec_api_get_crust_age_myr_map(simulation);
    const size_t cell_count = static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    for (size_t i = 0; i < cell_count; i += 509) {
        EXPECT_GE(crust_age_myr[i], 0.0f);
    }

    platec_api_destroy(simulation);
}

TEST(Phase2Contract, TimeMyrRemainsWholeRunMonotonicAcrossRestart) {
    platec::scenario::Scenario scenario = make_phase2_scenario(1.25);
    scenario.cycle_step_limit = 1;
    scenario.cycle_count = 3;

    lithosphere simulation(scenario);
    const uint32_t origin_step = simulation.getTimeOriginStep();

    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 0.0);
    EXPECT_EQ(simulation.getCycleCount(), 0U);

    ASSERT_FALSE(simulation.isFinished());
    simulation.update();
    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 1.25);
    EXPECT_EQ(simulation.getCycleCount(), 0U);
    EXPECT_EQ(simulation.getTimeOriginStep(), origin_step);

    ASSERT_FALSE(simulation.isFinished());
    simulation.update();
    EXPECT_EQ(simulation.getCycleCount(), 1U);
    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 1.25);
    EXPECT_EQ(simulation.getTimeOriginStep(), origin_step);

    ASSERT_FALSE(simulation.isFinished());
    simulation.update();
    EXPECT_EQ(simulation.getCycleCount(), 1U);
    EXPECT_DOUBLE_EQ(simulation.getTimeMyr(), 2.5);
    EXPECT_EQ(simulation.getTimeOriginStep(), origin_step);
    EXPECT_EQ(simulation.getIterationCount() - origin_step, 2U);

    const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
    EXPECT_EQ(snapshot.cycle_count, 1U);
    EXPECT_EQ(snapshot.time_origin_step, origin_step);
    EXPECT_DOUBLE_EQ(snapshot.time_myr, 2.5);
}

TEST(Phase2Contract, CApiTimeMyrRemainsWholeRunMonotonicAcrossRestart) {
    platec::scenario::Scenario scenario = make_phase2_scenario(0.5);
    scenario.cycle_step_limit = 1;
    scenario.cycle_count = 3;

    void* simulation = platec_api_create_from_scenario(&scenario);
    ASSERT_NE(simulation, nullptr);
    const uint32_t origin_step = platec_api_get_time_origin_step(simulation);

    ASSERT_EQ(platec_api_is_finished(simulation), 0U);
    platec_api_step(simulation);
    EXPECT_DOUBLE_EQ(platec_api_get_time_myr(simulation), 0.5);
    EXPECT_EQ(platec_api_get_cycle_count(simulation), 0U);
    EXPECT_EQ(platec_api_get_time_origin_step(simulation), origin_step);

    ASSERT_EQ(platec_api_is_finished(simulation), 0U);
    platec_api_step(simulation);
    EXPECT_DOUBLE_EQ(platec_api_get_time_myr(simulation), 0.5);
    EXPECT_EQ(platec_api_get_cycle_count(simulation), 1U);
    EXPECT_EQ(platec_api_get_time_origin_step(simulation), origin_step);

    ASSERT_EQ(platec_api_is_finished(simulation), 0U);
    platec_api_step(simulation);
    EXPECT_DOUBLE_EQ(platec_api_get_time_myr(simulation), 1.0);
    EXPECT_EQ(platec_api_get_cycle_count(simulation), 1U);
    EXPECT_EQ(platec_api_get_iteration_count(simulation) - origin_step, 2U);

    platec_api_destroy(simulation);
}
