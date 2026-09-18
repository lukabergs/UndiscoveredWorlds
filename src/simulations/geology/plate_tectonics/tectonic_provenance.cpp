#include "tectonic_provenance.hpp"

#include "topography_codec.hpp"
#include "tectonic_kinematics.hpp"

#include <algorithm>
#include <cmath>

namespace platec::provenance {

namespace {

size_t map_index(int x, int y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
}

float clamp01(float value) {
    return std::clamp(value, 0.0f, 1.0f);
}

bool is_sea(const Inputs& inputs, int x, int y) {
    return TopographyCodec::is_oceanic_internal(
        inputs.heightmap[map_index(x, y, inputs.width)]);
}

bool is_continental(const Inputs& inputs, size_t index) {
    if (inputs.crust_class == nullptr) {
        return !TopographyCodec::is_oceanic_internal(inputs.heightmap[index]);
    }
    const auto crust = static_cast<contract::CrustClass>(inputs.crust_class[index]);
    return crust == contract::CrustClass::Continental ||
           crust == contract::CrustClass::Transitional;
}

bool has_adjacent_land(const Inputs& inputs, int x, int y) {
    const int width = static_cast<int>(inputs.width);
    const int height = static_cast<int>(inputs.height);

    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int ny = (y + dy + height) % height;
            const int nx = (x + dx + width) % width;
            if (!is_sea(inputs, nx, ny)) {
                return true;
            }
        }
    }

    return false;
}

bool has_adjacent_sea(const Inputs& inputs, int x, int y) {
    const int width = static_cast<int>(inputs.width);
    const int height = static_cast<int>(inputs.height);

    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int ny = (y + dy + height) % height;
            const int nx = (x + dx + width) % width;
            if (is_sea(inputs, nx, ny)) {
                return true;
            }
        }
    }

    return false;
}

bool is_coast_or_outline(const Inputs& inputs, int x, int y) {
    return is_sea(inputs, x, y) ? has_adjacent_land(inputs, x, y)
                               : has_adjacent_sea(inputs, x, y);
}

uint8_t normalize_tectonic_signal(float signal) {
    const float normalized = clamp01(
        (signal - contract::kAdapterProxyConvergentBoundaryThreshold) /
        (2.0f - contract::kAdapterProxyConvergentBoundaryThreshold));
    return static_cast<uint8_t>(std::lround(normalized * 100.0f));
}

contract::GeologicRegime classify_geologic_regime(const Inputs& inputs, int x, int y,
                                                  uint8_t convergence, uint8_t divergence,
                                                  uint8_t shear,
                                                  contract::GeologicRegime convergent_regime) {
    const uint8_t max_signal = std::max({convergence, divergence, shear});

    if (convergence >= 65U) {
        return convergent_regime;
    }

    if (divergence >= 55U) {
        return is_continental(inputs, map_index(x, y, inputs.width))
                   ? contract::GeologicRegime::DivergentRift
                   : contract::GeologicRegime::MidOceanRidge;
    }

    if (shear >= 50U && shear >= static_cast<uint8_t>(convergence + 10U) &&
        shear >= static_cast<uint8_t>(divergence + 10U)) {
        return contract::GeologicRegime::Transform;
    }

    if (is_coast_or_outline(inputs, x, y) && max_signal <= 20U) {
        return contract::GeologicRegime::PassiveMargin;
    }

    return contract::GeologicRegime::Stable;
}

void record_boundary_motion(const Inputs& inputs, int x, int y, int nx, int ny,
                            float normal_x, float normal_y,
                            std::vector<float>& convergence,
                            std::vector<float>& divergence,
                            std::vector<float>& shear,
                            std::vector<contract::GeologicRegime>& convergent_regime) {
    if (nx < 0 || nx >= static_cast<int>(inputs.width) ||
        ny < 0 || ny >= static_cast<int>(inputs.height)) {
        return;
    }

    const size_t index = map_index(x, y, inputs.width);
    const size_t neighbor_index = map_index(nx, ny, inputs.width);
    const uint32_t plate = inputs.plate_id[index];
    const uint32_t neighbor_plate = inputs.plate_id[neighbor_index];

    if (plate == neighbor_plate || plate >= inputs.plate_count ||
        neighbor_plate >= inputs.plate_count) {
        return;
    }

    const float boundary_x = static_cast<float>(x) + 0.5f * normal_x;
    const float boundary_y = static_cast<float>(y) + 0.5f * normal_y;
    const auto first = kinematics::local_velocity(inputs.plates[plate], boundary_x,
                                                 boundary_y, inputs.width, inputs.height);
    const auto second = kinematics::local_velocity(inputs.plates[neighbor_plate], boundary_x,
                                                  boundary_y, inputs.width, inputs.height);
    const float relative_x = first.x - second.x;
    const float relative_y = first.y - second.y;
    if (!std::isfinite(relative_x) || !std::isfinite(relative_y)) return;
    const float boundary_motion = relative_x * normal_x + relative_y * normal_y;
    const float tangential_motion =
        std::fabs(relative_x * (-normal_y) + relative_y * normal_x);

    if (boundary_motion > 0.0f) {
        const bool first_continent = is_continental(inputs, index);
        const bool second_continent = is_continental(inputs, neighbor_index);
        auto regime = [](bool own_continent, bool other_continent) {
            if (own_continent && other_continent)
                return contract::GeologicRegime::ContinentCollision;
            if (!own_continent && other_continent)
                return contract::GeologicRegime::TrenchAdjacent;
            // Ocean-ocean convergence is an arc system; polarity is not inferred here.
            return contract::GeologicRegime::ConvergentArc;
        };
        if (boundary_motion > convergence[index])
            convergent_regime[index] = regime(first_continent, second_continent);
        if (boundary_motion > convergence[neighbor_index])
            convergent_regime[neighbor_index] = regime(second_continent, first_continent);
        convergence[index] = std::max(convergence[index], boundary_motion);
        convergence[neighbor_index] = std::max(convergence[neighbor_index], boundary_motion);
    } else if (boundary_motion < 0.0f) {
        const float divergence_signal = -boundary_motion;
        divergence[index] = std::max(divergence[index], divergence_signal);
        divergence[neighbor_index] = std::max(divergence[neighbor_index], divergence_signal);
    }

    shear[index] = std::max(shear[index], tangential_motion);
    shear[neighbor_index] = std::max(shear[neighbor_index], tangential_motion);
}

} // namespace

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime) {
    if (inputs.width == 0 || inputs.height == 0 || inputs.heightmap == nullptr ||
        inputs.plate_id == nullptr || convergence_score == nullptr ||
        divergence_score == nullptr || shear_score == nullptr ||
        geologic_regime == nullptr) {
        return;
    }

    const size_t cell_count =
        static_cast<size_t>(inputs.width) * static_cast<size_t>(inputs.height);
    std::fill(convergence_score, convergence_score + cell_count, static_cast<uint8_t>(0));
    std::fill(divergence_score, divergence_score + cell_count, static_cast<uint8_t>(0));
    std::fill(shear_score, shear_score + cell_count, static_cast<uint8_t>(0));
    std::fill(geologic_regime, geologic_regime + cell_count,
              static_cast<uint8_t>(contract::GeologicRegime::Stable));

    if (inputs.plate_count == 0 || inputs.plates == nullptr) {
        return;
    }

    std::vector<float> convergence(cell_count, 0.0f);
    std::vector<float> divergence(cell_count, 0.0f);
    std::vector<float> shear(cell_count, 0.0f);
    std::vector<contract::GeologicRegime> convergent_regime(
        cell_count, contract::GeologicRegime::Stable);
    for (uint32_t y = 0; y < inputs.height; ++y) {
        for (uint32_t x = 0; x < inputs.width; ++x) {
            if (inputs.width > 1U) {
                record_boundary_motion(inputs, static_cast<int>(x), static_cast<int>(y),
                                       static_cast<int>((x + 1U) % inputs.width), static_cast<int>(y),
                                       1.0f, 0.0f, convergence, divergence, shear, convergent_regime);
            }
            if (inputs.height > 1U) {
                record_boundary_motion(inputs, static_cast<int>(x), static_cast<int>(y),
                                       static_cast<int>(x), static_cast<int>((y + 1U) % inputs.height),
                                       0.0f, 1.0f, convergence, divergence, shear, convergent_regime);
            }
        }
    }

    for (uint32_t y = 0; y < inputs.height; ++y) {
        for (uint32_t x = 0; x < inputs.width; ++x) {
            const size_t index = map_index(static_cast<int>(x), static_cast<int>(y),
                                           inputs.width);
            convergence_score[index] = normalize_tectonic_signal(convergence[index]);
            divergence_score[index] = normalize_tectonic_signal(divergence[index]);
            shear_score[index] = normalize_tectonic_signal(shear[index]);
            geologic_regime[index] = static_cast<uint8_t>(classify_geologic_regime(
                inputs, static_cast<int>(x), static_cast<int>(y), convergence_score[index],
                divergence_score[index], shear_score[index], convergent_regime[index]));
        }
    }
}

} // namespace platec::provenance
