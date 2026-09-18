#include "tectonic_junction_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace platec::tectonic_junction_coupling {
namespace {

using platec::contract::JunctionArmType;
using platec::contract::JunctionStability;
using platec::contract::TectonicJunction;

struct Vec2 {
    float x = 0.0f;
    float y = 0.0f;
};

size_t map_index(uint32_t x, uint32_t y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
}

uint32_t wrapped_add(uint32_t value, int offset, uint32_t limit) {
    int64_t shifted = static_cast<int64_t>(value) + static_cast<int64_t>(offset);
    const int64_t period = static_cast<int64_t>(limit);
    shifted %= period;
    if (shifted < 0) {
        shifted += period;
    }
    return static_cast<uint32_t>(shifted);
}

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
}

float wrapped_delta(float value, float center, uint32_t period) {
    float delta = value - center;
    const float period_f = static_cast<float>(period);
    const float half_period = period_f * 0.5f;
    if (delta > half_period) {
        delta -= period_f;
    } else if (delta < -half_period) {
        delta += period_f;
    }
    return delta;
}

float length(Vec2 value) {
    return std::sqrt(value.x * value.x + value.y * value.y);
}

Vec2 normalized(Vec2 value) {
    const float magnitude = length(value);
    if (magnitude <= 1e-6f || !std::isfinite(magnitude)) {
        return {};
    }
    return Vec2{value.x / magnitude, value.y / magnitude};
}

uint32_t count_arms(const TectonicJunction& junction, JunctionArmType type) {
    uint32_t count = 0U;
    for (JunctionArmType arm_type : junction.arm_types) {
        if (arm_type == type) {
            ++count;
        }
    }
    return count;
}

float stable_activity_score(const TectonicJunction& junction) {
    if (junction.stability == JunctionStability::Stable) {
        return clamp_unit(std::max(junction.stability_score,
                                   junction.geometry_quality_score * 0.65f));
    }
    if (junction.stability == JunctionStability::ConditionallyStable) {
        return clamp_unit(std::max(0.25f, junction.stability_score * 0.75f));
    }
    return 0.0f;
}

float unstable_activity_score(const TectonicJunction& junction) {
    if (junction.stability != JunctionStability::Unstable) {
        return 0.0f;
    }
    return clamp_unit(std::max(0.35f, 1.0f - junction.stability_score) *
                      std::max(0.35f, junction.quality_score));
}

Vec2 fallback_rift_direction(const TectonicJunction& junction) {
    for (size_t i = 0; i < 3U; ++i) {
        if (junction.arm_types[i] != JunctionArmType::Ridge ||
            junction.arm_orientation_quality[i] <= 0.0f) {
            continue;
        }
        Vec2 normal{-junction.arm_strike_y[i], junction.arm_strike_x[i]};
        normal = normalized(normal);
        if (normal.x < 0.0f || (std::fabs(normal.x) <= 1e-6f && normal.y < 0.0f)) {
            normal.x = -normal.x;
            normal.y = -normal.y;
        }
        if (length(normal) > 0.0f) {
            return normal;
        }
    }
    return Vec2{1.0f, 0.0f};
}

void clamp_tension_vectors(Outputs& outputs) {
    for (size_t i = 0; i < outputs.rift_tension_x.size(); ++i) {
        const Vec2 tension{outputs.rift_tension_x[i], outputs.rift_tension_y[i]};
        const float magnitude = length(tension);
        if (magnitude <= 1.0f) {
            continue;
        }
        outputs.rift_tension_x[i] = tension.x / magnitude;
        outputs.rift_tension_y[i] = tension.y / magnitude;
    }
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    const size_t cell_count =
        static_cast<size_t>(inputs.width) * static_cast<size_t>(inputs.height);
    outputs.rift_influence.assign(cell_count, 0.0f);
    outputs.transform_influence.assign(cell_count, 0.0f);
    outputs.subduction_influence.assign(cell_count, 0.0f);
    outputs.rift_tension_x.assign(cell_count, 0.0f);
    outputs.rift_tension_y.assign(cell_count, 0.0f);

    if (inputs.width == 0U || inputs.height == 0U || inputs.junctions == nullptr ||
        inputs.junction_count == 0U) {
        return outputs;
    }

    for (uint32_t i = 0; i < inputs.junction_count; ++i) {
        const TectonicJunction& junction = inputs.junctions[i];
        const float stable_score = stable_activity_score(junction);
        const float unstable_score = unstable_activity_score(junction);
        const uint32_t ridge_count = count_arms(junction, JunctionArmType::Ridge);
        const uint32_t fault_count = count_arms(junction, JunctionArmType::Fault);
        const uint32_t trench_count = count_arms(junction, JunctionArmType::Trench);

        const float rift_strength = stable_score * static_cast<float>(ridge_count) / 3.0f;
        const float subduction_strength =
            stable_score * static_cast<float>(trench_count) / 3.0f;
        const float transform_strength =
            stable_score * static_cast<float>(fault_count) / 4.0f +
            unstable_score * std::max(0.5f, static_cast<float>(fault_count) / 3.0f);
        if (rift_strength <= 0.0f && subduction_strength <= 0.0f &&
            transform_strength <= 0.0f) {
            continue;
        }

        const float radius =
            2.5f + 4.5f * clamp_unit(std::max(junction.quality_score,
                                              junction.geometry_quality_score));
        const int cell_radius = static_cast<int>(std::ceil(radius));
        const uint32_t center_x =
            wrapped_add(0U, static_cast<int>(std::floor(junction.centroid_x)), inputs.width);
        const uint32_t center_y =
            wrapped_add(0U, static_cast<int>(std::floor(junction.centroid_y)), inputs.height);
        const Vec2 fallback_direction = fallback_rift_direction(junction);

        for (int oy = -cell_radius; oy <= cell_radius; ++oy) {
            for (int ox = -cell_radius; ox <= cell_radius; ++ox) {
                const uint32_t x = wrapped_add(center_x, ox, inputs.width);
                const uint32_t y = wrapped_add(center_y, oy, inputs.height);
                const float dx =
                    wrapped_delta(static_cast<float>(x) + 0.5f, junction.centroid_x,
                                  inputs.width);
                const float dy =
                    wrapped_delta(static_cast<float>(y) + 0.5f, junction.centroid_y,
                                  inputs.height);
                const float distance = std::sqrt(dx * dx + dy * dy);
                if (distance > radius) {
                    continue;
                }

                const float falloff = clamp_unit(1.0f - distance / (radius + 0.5f));
                const float local_weight = falloff * falloff;
                const size_t index = map_index(x, y, inputs.width);
                outputs.rift_influence[index] =
                    std::max(outputs.rift_influence[index], local_weight * rift_strength);
                outputs.transform_influence[index] =
                    std::max(outputs.transform_influence[index],
                             local_weight * transform_strength);
                outputs.subduction_influence[index] =
                    std::max(outputs.subduction_influence[index],
                             local_weight * subduction_strength);

                if (rift_strength <= 0.0f) {
                    continue;
                }
                Vec2 direction = normalized(Vec2{dx, dy});
                if (length(direction) <= 0.0f) {
                    direction = fallback_direction;
                }
                const float tension = local_weight * rift_strength;
                outputs.rift_tension_x[index] += direction.x * tension;
                outputs.rift_tension_y[index] += direction.y * tension;
            }
        }
    }

    clamp_tension_vectors(outputs);
    return outputs;
}

} // namespace platec::tectonic_junction_coupling
