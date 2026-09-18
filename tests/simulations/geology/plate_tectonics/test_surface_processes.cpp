#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "topography_codec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {

platec::scenario::Scenario make_surface_scenario(long seed) {
    platec::scenario::Scenario scenario;
    scenario.seed = seed;
    scenario.width = 96;
    scenario.height = 48;
    scenario.sea_level = 0.65f;
    scenario.erosion_period = 0;
    scenario.folding_ratio = 0.12f;
    scenario.aggregation_overlap_abs = 1000000;
    scenario.aggregation_overlap_rel = 0.33f;
    scenario.cycle_count = 3;
    scenario.plate_count = 8;
    scenario.cycle_step_limit = 140;
    scenario.rotation_strength = 1.0f;
    scenario.crust_rotation_strength = 0.20f;
    scenario.erosion_strength = 1.0f;
    scenario.delta_time_myr = 1.0;
    scenario.initial_max_height_m = TopographyCodec::kMaxHeightMeters;
    return scenario;
}

void run_updates(lithosphere& simulation, uint32_t updates) {
    for (uint32_t i = 0; i < updates && !simulation.isFinished(); ++i) {
        simulation.update();
    }
}

std::vector<float> copy_heightmap(const lithosphere& simulation) {
    const size_t map_size =
        static_cast<size_t>(simulation.getWidth()) * static_cast<size_t>(simulation.getHeight());
    return std::vector<float>(simulation.getTopography(), simulation.getTopography() + map_size);
}

float max_height(const lithosphere& simulation) {
    const std::vector<float> heightmap = copy_heightmap(simulation);
    return *std::max_element(heightmap.begin(), heightmap.end());
}

uint16_t max_height_meters(const lithosphere& simulation) {
    return TopographyCodec::internal_to_meters(max_height(simulation), simulation.getSeaLevelMeters());
}

float mean_absolute_difference(const std::vector<float>& lhs, const std::vector<float>& rhs) {
    if (lhs.size() != rhs.size() || lhs.empty()) {
        return 0.0f;
    }

    double total = 0.0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        total += std::fabs(static_cast<double>(lhs[i]) - static_cast<double>(rhs[i]));
    }
    return static_cast<float>(total / static_cast<double>(lhs.size()));
}

} // namespace

TEST(SurfaceProcesses, MaterialMapIsDeterministicAndNonUniform) {
    const platec::scenario::Scenario scenario = make_surface_scenario(12345);
    lithosphere primary(scenario);
    lithosphere reproducible(scenario);

    const size_t map_size =
        static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    const uint8_t* first_map = primary.getMaterialMap();
    const uint8_t* second_map = reproducible.getMaterialMap();

    ASSERT_NE(first_map, nullptr);
    ASSERT_NE(second_map, nullptr);
    EXPECT_EQ(0, std::memcmp(first_map, second_map, map_size));

    std::array<uint32_t, platec::material::kProperties.size()> counts = {};
    for (size_t i = 0; i < map_size; ++i) {
        ++counts[std::min<size_t>(first_map[i], counts.size() - 1U)];
    }

    size_t non_zero_buckets = 0;
    for (uint32_t count : counts) {
        non_zero_buckets += count > 0U ? 1U : 0U;
    }

    EXPECT_GE(non_zero_buckets, 2U);
}

TEST(SurfaceProcesses, InitialCrustAgeSeedingIsDeterministicAndRidgeBiased) {
    const platec::scenario::Scenario scenario = make_surface_scenario(67890);
    lithosphere first(scenario);
    lithosphere second(scenario);

    const size_t map_size =
        static_cast<size_t>(scenario.width) * static_cast<size_t>(scenario.height);
    const float* first_age = first.getCrustAgeMyrMap();
    const float* second_age = second.getCrustAgeMyrMap();
    const uint8_t* first_crust = first.getCrustClassMap();
    const uint8_t* first_boundary = first.getBoundaryTypeMap();

    ASSERT_NE(first_age, nullptr);
    ASSERT_NE(second_age, nullptr);
    ASSERT_NE(first_crust, nullptr);
    ASSERT_NE(first_boundary, nullptr);

    for (size_t i = 0; i < map_size; ++i) {
        EXPECT_FLOAT_EQ(first_age[i], second_age[i]);
    }

    double continental_total = 0.0;
    double oceanic_total = 0.0;
    double ridge_total = 0.0;
    size_t continental_count = 0U;
    size_t oceanic_count = 0U;
    size_t ridge_count = 0U;

    for (size_t i = 0; i < map_size; ++i) {
        const auto crust_class =
            static_cast<platec::contract::CrustClass>(first_crust[i]);
        if (crust_class == platec::contract::CrustClass::Oceanic) {
            oceanic_total += first_age[i];
            ++oceanic_count;
            if (static_cast<platec::contract::BoundaryType>(first_boundary[i]) ==
                platec::contract::BoundaryType::Divergent) {
                ridge_total += first_age[i];
                ++ridge_count;
            }
        } else if (crust_class != platec::contract::CrustClass::None) {
            continental_total += first_age[i];
            ++continental_count;
        }
    }

    ASSERT_GT(continental_count, 0U);
    ASSERT_GT(oceanic_count, 0U);
    ASSERT_GT(ridge_count, 0U);

    const double continental_mean = continental_total / static_cast<double>(continental_count);
    const double oceanic_mean = oceanic_total / static_cast<double>(oceanic_count);
    const double ridge_mean = ridge_total / static_cast<double>(ridge_count);

    EXPECT_GT(continental_mean, oceanic_mean + 20.0);
    EXPECT_LT(ridge_mean, 2.0);
    EXPECT_LT(ridge_mean, oceanic_mean);
}

TEST(SurfaceProcesses, MaxInitialHeightScalesSeededTopography) {
    platec::scenario::Scenario low = make_surface_scenario(4242);
    low.initial_max_height_m = 8000U;

    platec::scenario::Scenario high = low;
    high.initial_max_height_m = TopographyCodec::kMaxHeightMeters;

    lithosphere low_sim(low);
    lithosphere high_sim(high);

    EXPECT_LT(max_height(low_sim), max_height(high_sim) - 0.45f);
}

TEST(SurfaceProcesses, LowerGravityAllowsTallerRelief) {
    platec::scenario::Scenario low_gravity = make_surface_scenario(24680);
    low_gravity.gravity_mps2 = 3.71f;
    low_gravity.glacial_erosion_strength = 0.0f;
    low_gravity.erosion_strength = 0.0f;

    platec::scenario::Scenario high_gravity = low_gravity;
    high_gravity.gravity_mps2 = 19.62f;

    lithosphere low(low_gravity);
    lithosphere high(high_gravity);
    run_updates(low, 96);
    run_updates(high, 96);

    EXPECT_GT(max_height(low), max_height(high) + 0.05f);
}

TEST(SurfaceProcesses, GlacialBuzzsawCutsExtremePeaks) {
    platec::scenario::Scenario no_glaciers = make_surface_scenario(13579);
    no_glaciers.glacial_erosion_strength = 0.0f;

    platec::scenario::Scenario strong_glaciers = no_glaciers;
    strong_glaciers.glacial_erosion_strength = 2.5f;

    lithosphere baseline(no_glaciers);
    lithosphere glacial(strong_glaciers);
    run_updates(baseline, 96);
    run_updates(glacial, 96);

    EXPECT_GT(max_height(baseline), max_height(glacial) + 0.05f);
}

TEST(SurfaceProcesses, PeriodicNoiseIsDeterministicAndChangesTerrain) {
    platec::scenario::Scenario noisy = make_surface_scenario(31415);
    noisy.hf_noise_period = 2U;
    noisy.hf_noise_strength = 1.0f;
    noisy.lf_noise_period = 5U;
    noisy.lf_noise_strength = 1.0f;

    platec::scenario::Scenario baseline = noisy;
    baseline.hf_noise_period = 0U;
    baseline.hf_noise_strength = 0.0f;
    baseline.lf_noise_period = 0U;
    baseline.lf_noise_strength = 0.0f;

    lithosphere first(noisy);
    lithosphere second(noisy);
    lithosphere control(baseline);
    run_updates(first, 12);
    run_updates(second, 12);
    run_updates(control, 12);

    const std::vector<float> first_map = copy_heightmap(first);
    const std::vector<float> second_map = copy_heightmap(second);
    const std::vector<float> control_map = copy_heightmap(control);

    ASSERT_EQ(first_map.size(), second_map.size());
    for (size_t i = 0; i < first_map.size(); ++i) {
        EXPECT_FLOAT_EQ(first_map[i], second_map[i]);
    }

    EXPECT_GT(mean_absolute_difference(first_map, control_map), 0.0001f);
}

TEST(SurfaceProcesses, EarthGravityAvoidsCodecSaturation) {
    platec::scenario::Scenario scenario = make_surface_scenario(11223);
    scenario.gravity_mps2 = platec::scenario::kDefaultGravityMps2;
    scenario.glacial_erosion_strength = 1.0f;
    scenario.erosion_strength = 1.0f;

    lithosphere simulation(scenario);
    run_updates(simulation, 120);

    EXPECT_LT(max_height_meters(simulation), TopographyCodec::kMaxHeightMeters);
}
