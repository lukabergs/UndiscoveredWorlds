#include "tectonic_provenance.hpp"

#include "topography_codec.hpp"

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

bool has_adjacent_land(const Inputs& inputs, int x, int y) {
    const int width = static_cast<int>(inputs.width);
    const int height = static_cast<int>(inputs.height);

    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int ny = y + dy;
            if (ny < 0 || ny >= height) {
                continue;
            }

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

            const int ny = y + dy;
            if (ny < 0 || ny >= height) {
                continue;
            }

            const int nx = (x + dx + width) % width;
            if (is_sea(inputs, nx, ny)) {
                return true;
            }
        }
    }

    return false;
}

bool is_coast_or_outline(const Inputs& inputs, int x, int y) {
    const bool sea = is_sea(inputs, x, y);
    const int width = static_cast<int>(inputs.width);
    const int height = static_cast<int>(inputs.height);

    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            if (dx == 0 && dy == 0) {
                continue;
            }

            const int ny = y + dy;
            if (ny < 0 || ny >= height) {
                continue;
            }

            const int nx = (x + dx + width) % width;
            if (is_sea(inputs, nx, ny) != sea) {
                return true;
            }
        }
    }

    return false;
}

uint8_t normalize_tectonic_signal(float signal) {
    const float normalized = clamp01(
        (signal - contract::kAdapterProxyConvergentBoundaryThreshold) /
        (2.0f - contract::kAdapterProxyConvergentBoundaryThreshold));
    return static_cast<uint8_t>(std::lround(normalized * 100.0f));
}

contract::GeologicRegime classify_geologic_regime(const Inputs& inputs, int x, int y,
                                                  uint8_t convergence, uint8_t divergence,
                                                  uint8_t shear) {
    const bool sea = is_sea(inputs, x, y);
    const bool adjacent_land = has_adjacent_land(inputs, x, y);
    const bool adjacent_sea = has_adjacent_sea(inputs, x, y);
    const uint8_t max_signal = std::max({convergence, divergence, shear});

    if (convergence >= 65U) {
        if (sea && adjacent_land) {
            return contract::GeologicRegime::TrenchAdjacent;
        }

        if (!sea && adjacent_sea) {
            return contract::GeologicRegime::ConvergentArc;
        }

        return contract::GeologicRegime::ContinentCollision;
    }

    if (divergence >= 55U) {
        return sea ? contract::GeologicRegime::MidOceanRidge
                   : contract::GeologicRegime::DivergentRift;
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
                            std::vector<float>& convergence,
                            std::vector<float>& divergence,
                            std::vector<float>& shear) {
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

    float normal_x = static_cast<float>(nx - x);
    float normal_y = static_cast<float>(ny - y);
    const float length = std::sqrt(normal_x * normal_x + normal_y * normal_y);
    if (length <= 0.0f) {
        return;
    }

    normal_x /= length;
    normal_y /= length;

    const float relative_x =
        inputs.plates[plate].unit_velocity_x - inputs.plates[neighbor_plate].unit_velocity_x;
    const float relative_y =
        inputs.plates[plate].unit_velocity_y - inputs.plates[neighbor_plate].unit_velocity_y;
    const float boundary_motion = relative_x * normal_x + relative_y * normal_y;
    const float tangential_motion =
        std::fabs(relative_x * (-normal_y) + relative_y * normal_x);

    if (boundary_motion > 0.0f) {
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
    for (uint32_t y = 0; y < inputs.height; ++y) {
        for (uint32_t x = 0; x < inputs.width; ++x) {
            if (x + 1U < inputs.width) {
                record_boundary_motion(inputs, static_cast<int>(x), static_cast<int>(y),
                                       static_cast<int>(x + 1U), static_cast<int>(y),
                                       convergence, divergence, shear);
            }
            if (y + 1U < inputs.height) {
                record_boundary_motion(inputs, static_cast<int>(x), static_cast<int>(y),
                                       static_cast<int>(x), static_cast<int>(y + 1U),
                                       convergence, divergence, shear);
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
                divergence_score[index], shear_score[index]));
        }
    }
}

} // namespace platec::provenance
