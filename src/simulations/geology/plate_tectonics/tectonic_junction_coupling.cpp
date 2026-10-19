#include "simulations/geology/plate_tectonics/tectonic_junction_coupling.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace tectonics::tectonic_junction_coupling {
namespace {

using tectonics::contract::JunctionArmType;
using tectonics::contract::JunctionStability;
using tectonics::contract::TectonicJunction;

struct Vec2 {
    float x = 0.0f;
    float y = 0.0f;
};

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
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

// One junction's influence at one cell, at distance (cells) and with the direction away from the
// junction (zero length: the fallback).
struct Strengths {
    float radius = 0.0f, rift = 0.0f, transform = 0.0f, subduction = 0.0f;
    Vec2 fallback_direction;
};

void apply(const Strengths& strengths, size_t index, float distance, Vec2 away, Outputs& outputs) {
    const float falloff = clamp_unit(1.0f - distance / (strengths.radius + 0.5f));
    const float local_weight = falloff * falloff;
    outputs.rift_influence[index] =
        std::max(outputs.rift_influence[index], local_weight * strengths.rift);
    outputs.transform_influence[index] =
        std::max(outputs.transform_influence[index],
                 local_weight * strengths.transform);
    outputs.subduction_influence[index] =
        std::max(outputs.subduction_influence[index],
                 local_weight * strengths.subduction);

    if (strengths.rift <= 0.0f) {
        return;
    }
    Vec2 direction = normalized(away);
    if (length(direction) <= 0.0f) {
        direction = strengths.fallback_direction;
    }
    const float tension = local_weight * strengths.rift;
    outputs.rift_tension_x[index] += direction.x * tension;
    outputs.rift_tension_y[index] += direction.y * tension;
}

// The cells within the radius (great circle over the cell spacing) of the junction (longitude,
// latitude), the tension (east, north) at each cell pointing away from it.
void apply_disc(const grid::TectonicGrid& grid, const TectonicJunction& junction,
                const Strengths& strengths, Outputs& outputs) {
    const double spacing = grid.cellspacing_m();
    const double cells_per_radian = grid.mesh()->radius_m() / spacing;
    const mesh::Vec3 p = boundary_cube::decodecentroid(junction.centroid_x, junction.centroid_y);
    for (const grid::CellId c :
         mesh::cellswithin(*grid.mesh(), p, static_cast<double>(strengths.radius) * spacing)) {
        const mesh::Vec3 q = grid.mesh()->centre(c);
        const float distance = static_cast<float>(mesh::greatcircle_rad(p, q) * cells_per_radian);
        if (distance > strengths.radius) {
            continue;
        }
        const auto away = boundary_cube::eastnorth(q, cells_per_radian * (mesh::dot(p, q) * q - p));
        apply(strengths, c, distance, Vec2{static_cast<float>(away[0]), static_cast<float>(away[1])},
              outputs);
    }
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    const size_t cell_count = inputs.grid != nullptr ? inputs.grid->cellcount() : 0U;
    outputs.rift_influence.assign(cell_count, 0.0f);
    outputs.transform_influence.assign(cell_count, 0.0f);
    outputs.subduction_influence.assign(cell_count, 0.0f);
    outputs.rift_tension_x.assign(cell_count, 0.0f);
    outputs.rift_tension_y.assign(cell_count, 0.0f);

    if (cell_count == 0U || inputs.junctions == nullptr || inputs.junction_count == 0U) {
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

        Strengths strengths;
        strengths.radius =
            2.5f + 4.5f * clamp_unit(std::max(junction.quality_score,
                                              junction.geometry_quality_score));
        strengths.rift = rift_strength;
        strengths.transform = transform_strength;
        strengths.subduction = subduction_strength;
        strengths.fallback_direction = fallback_rift_direction(junction);
        apply_disc(*inputs.grid, junction, strengths, outputs);
    }

    clamp_tension_vectors(outputs);
    return outputs;
}

} // namespace tectonics::tectonic_junction_coupling
