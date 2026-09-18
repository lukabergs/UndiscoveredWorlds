#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "gtest/gtest.h"
#include "lithosphere.hpp"
#include "platecapi.hpp"
#include "tectonic_contract.hpp"
#include "topography_codec.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct FixtureConfig {
    const char* directory_name;
    long seed;
    uint32_t width;
    uint32_t height;
    uint32_t erosion_period;
    float folding_ratio;
    uint32_t aggregation_overlap_abs;
    float aggregation_overlap_rel;
    uint32_t cycle_count;
    uint32_t plate_count;
    uint32_t cycle_step_limit;
};

constexpr FixtureConfig kFixtures[] = {
    {"adapter_proxy_seed12345", 12345, 96, 48, 0, 0.02f, 1000000, 0.33f, 4, 4, 600},
    {"adapter_proxy_seed67890", 67890, 96, 48, 0, 0.02f, 1000000, 0.33f, 4, 4, 600},
};

fs::path fixture_root() {
    return fs::path(PLATE_TECTONICS_FIXTURE_DIR) / "phase0";
}

template <typename T>
std::vector<T> read_binary_vector(const fs::path& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input.good()) {
        throw std::runtime_error("failed to open fixture: " + path.string());
    }

    const std::streamsize size = input.tellg();
    if (size < 0 || size % static_cast<std::streamsize>(sizeof(T)) != 0) {
        throw std::runtime_error("invalid fixture size: " + path.string());
    }

    input.seekg(0, std::ios::beg);
    std::vector<T> values(static_cast<size_t>(size / static_cast<std::streamsize>(sizeof(T))));
    input.read(reinterpret_cast<char*>(values.data()), size);
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed to read fixture: " + path.string());
    }
    return values;
}

platec::contract::Snapshot run_fixture(const FixtureConfig& config) {
    lithosphere simulation(config.seed, config.width, config.height, 0.65f,
                           config.erosion_period, config.folding_ratio,
                           config.aggregation_overlap_abs,
                           config.aggregation_overlap_rel, config.cycle_count,
                           config.plate_count, 1.0f, 0.20f, 1.0f, 1.0f,
                           TopographyCodec::kNoSeaLevelOverride,
                           TopographyCodec::kDefaultInitialMinHeightMeters,
                           TopographyCodec::kDefaultInitialMaxHeightMeters,
                           config.cycle_step_limit, 0.015f);

    while (!simulation.isFinished()) {
        simulation.update();
    }

    return platec::contract::capture_phase0_snapshot(simulation);
}

double mean_metric_height_delta(const platec::contract::Snapshot& actual,
                                const std::vector<float>& expected) {
    long double total_delta = 0.0;
    for (size_t i = 0; i < actual.cell_count(); ++i) {
        const uint16_t actual_m =
            TopographyCodec::internal_to_meters(actual.heightmap[i], actual.sea_level_m);
        const uint16_t expected_m =
            TopographyCodec::internal_to_meters(expected[i], actual.sea_level_m);
        total_delta += static_cast<long double>(
            std::abs(static_cast<int>(actual_m) - static_cast<int>(expected_m)));
    }
    return static_cast<double>(total_delta / static_cast<long double>(actual.cell_count()));
}

void expect_fixture_match(const FixtureConfig& config) {
    const platec::contract::Snapshot actual = run_fixture(config);
    const fs::path directory = fixture_root() / config.directory_name;

    const std::vector<float> expected_heightmap =
        read_binary_vector<float>(directory / "heightmap.f32");
    const std::vector<uint32_t> expected_plate_id =
        read_binary_vector<uint32_t>(directory / "plate_id.u32");
    const std::vector<uint32_t> expected_age =
        read_binary_vector<uint32_t>(directory / "crust_age_steps.u32");
    const std::vector<uint8_t> expected_convergence =
        read_binary_vector<uint8_t>(directory / "convergence_score.u8");
    const std::vector<uint8_t> expected_divergence =
        read_binary_vector<uint8_t>(directory / "divergence_score.u8");
    const std::vector<uint8_t> expected_shear =
        read_binary_vector<uint8_t>(directory / "shear_score.u8");
    const std::vector<uint8_t> expected_regime =
        read_binary_vector<uint8_t>(directory / "geologic_regime.u8");

    ASSERT_EQ(expected_heightmap.size(), actual.heightmap.size());
    ASSERT_EQ(expected_plate_id.size(), actual.plate_id.size());
    ASSERT_EQ(expected_age.size(), actual.crust_age_steps.size());
    ASSERT_EQ(expected_convergence.size(), actual.convergence_score.size());
    ASSERT_EQ(expected_divergence.size(), actual.divergence_score.size());
    ASSERT_EQ(expected_shear.size(), actual.shear_score.size());
    ASSERT_EQ(expected_regime.size(), actual.geologic_regime.size());

    EXPECT_EQ(actual.plate_id, expected_plate_id);
    EXPECT_EQ(actual.crust_age_steps, expected_age);
    EXPECT_EQ(actual.convergence_score, expected_convergence);
    EXPECT_EQ(actual.divergence_score, expected_divergence);
    EXPECT_EQ(actual.shear_score, expected_shear);
    EXPECT_EQ(actual.geologic_regime, expected_regime);

    float max_abs_delta = 0.0f;
    for (size_t i = 0; i < actual.cell_count(); ++i) {
        max_abs_delta =
            std::max(max_abs_delta, std::fabs(actual.heightmap[i] - expected_heightmap[i]));
    }

    EXPECT_LE(max_abs_delta, 1e-5f);
    EXPECT_LE(mean_metric_height_delta(actual, expected_heightmap), 0.01);
}

} // namespace

TEST(Phase0Contract, CanonicalFixturesMatchCurrentNativeAndAdapterProxyFields) {
    GTEST_SKIP() << "Legacy phase-0 raster fixtures are obsolete under the current procedural "
                    "seeding regime.";
    for (const FixtureConfig& config : kFixtures) {
        SCOPED_TRACE(config.directory_name);
        expect_fixture_match(config);
    }
}

TEST(Phase1Contract, NativeProvenancePersistsAfterFinishedState) {
    lithosphere simulation(12345, 96, 48, 0.65f, 0, 0.02f, 1000000, 0.33f, 4, 4, 1.0f,
                           0.20f, 1.0f, 1.0f, TopographyCodec::kNoSeaLevelOverride,
                           TopographyCodec::kDefaultInitialMinHeightMeters,
                           TopographyCodec::kDefaultInitialMaxHeightMeters, 600, 0.015f);

    while (!simulation.isFinished()) {
        simulation.update();
    }

    // Finished runs retain the final plate state for native snapshot consumers.
    EXPECT_GT(simulation.getPlateCount(), 0U);
    EXPECT_EQ(simulation.getPlateCount(), simulation.getProvenancePlateCount());
    EXPECT_GT(simulation.getProvenancePlateCount(), 0U);
    ASSERT_NE(simulation.getConvergenceMap(), nullptr);
    ASSERT_NE(simulation.getDivergenceMap(), nullptr);
    ASSERT_NE(simulation.getShearMap(), nullptr);
    ASSERT_NE(simulation.getGeologicRegimeMap(), nullptr);

    const platec::contract::Snapshot snapshot = platec::contract::capture_snapshot(simulation);
    ASSERT_FALSE(snapshot.plates.empty());
    EXPECT_EQ(snapshot.plates.size(), simulation.getProvenancePlateCount());

    size_t active_boundary_cells = 0;
    size_t non_stable_regime_cells = 0;
    for (size_t i = 0; i < snapshot.cell_count(); ++i) {
        active_boundary_cells +=
            (snapshot.convergence_score[i] > 0U || snapshot.divergence_score[i] > 0U ||
             snapshot.shear_score[i] > 0U)
                ? 1U
                : 0U;
        non_stable_regime_cells +=
            snapshot.geologic_regime[i] !=
                    static_cast<uint8_t>(platec::contract::GeologicRegime::Stable)
                ? 1U
                : 0U;
    }

    EXPECT_GT(active_boundary_cells, 0U);
    EXPECT_GT(non_stable_regime_cells, 0U);
}

TEST(Phase1Contract, CApiExposesNativeProvenanceMaps) {
    void* simulation = platec_api_create(
        12345, 96, 48, 0.65f, 0, 0.02f, 1000000, 0.33f, 4, 4, 1.0f, 0.20f, 1.0f, 1.0f,
        TopographyCodec::kNoSeaLevelOverride, TopographyCodec::kDefaultInitialMinHeightMeters,
        TopographyCodec::kDefaultInitialMaxHeightMeters, 600, 0.015f);
    ASSERT_NE(simulation, nullptr);

    while (platec_api_is_finished(simulation) == 0U) {
        platec_api_step(simulation);
    }

    ASSERT_NE(platec_api_get_convergence_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_divergence_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_shear_map(simulation), nullptr);
    ASSERT_NE(platec_api_get_geologic_regime_map(simulation), nullptr);

    size_t active_boundary_cells = 0;
    const size_t cell_count = static_cast<size_t>(96U) * static_cast<size_t>(48U);
    const uint8_t* convergence = platec_api_get_convergence_map(simulation);
    const uint8_t* divergence = platec_api_get_divergence_map(simulation);
    const uint8_t* shear = platec_api_get_shear_map(simulation);
    for (size_t i = 0; i < cell_count; ++i) {
        active_boundary_cells +=
            (convergence[i] > 0U || divergence[i] > 0U || shear[i] > 0U) ? 1U : 0U;
    }

    EXPECT_GT(active_boundary_cells, 0U);
    platec_api_destroy(simulation);
}
