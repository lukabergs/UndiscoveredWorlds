#include "simulations/geology/plate_tectonics/tectonic_deforming_regions.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace tectonics::deforming_regions {

namespace {

using tectonics::contract::BoundarySegment;
using tectonics::contract::BoundaryType;
using tectonics::contract::CrustClass;
using tectonics::contract::DeformingRegion;
using tectonics::contract::DeformingRegionType;
using tectonics::contract::GeologicRegime;

float clamp_unit(float value) {
    return std::max(0.0f, std::min(1.0f, value));
}

CrustClass crust_class_from_byte(uint8_t value) {
    return static_cast<CrustClass>(value);
}

DeformingRegionType classify_region_type(const BoundarySegment& segment) {
    if (segment.boundary_type == BoundaryType::Divergent &&
        segment.geologic_regime == GeologicRegime::DivergentRift) {
        return DeformingRegionType::ContinentalRift;
    }
    if (segment.boundary_type == BoundaryType::Convergent &&
        segment.geologic_regime == GeologicRegime::ContinentCollision) {
        return DeformingRegionType::DiffuseCollision;
    }
    return DeformingRegionType::None;
}

uint16_t max_region_distance(DeformingRegionType region_type) {
    switch (region_type) {
    case DeformingRegionType::ContinentalRift:
        return 3U;
    case DeformingRegionType::DiffuseCollision:
        return 4U;
    case DeformingRegionType::None:
    default:
        return 0U;
    }
}

bool crust_supports_region(DeformingRegionType region_type, CrustClass crust_class) {
    switch (region_type) {
    case DeformingRegionType::ContinentalRift:
        return crust_class == CrustClass::Continental ||
               crust_class == CrustClass::Transitional;
    case DeformingRegionType::DiffuseCollision:
        return crust_class == CrustClass::Continental ||
               crust_class == CrustClass::Transitional;
    case DeformingRegionType::None:
    default:
        return false;
    }
}

float distance_falloff(uint16_t distance, uint16_t max_distance) {
    if (distance > max_distance) {
        return 0.0f;
    }
    return 1.0f - static_cast<float>(distance) /
                      static_cast<float>(max_distance + 1U);
}

float local_signal(const Inputs& inputs, size_t index, DeformingRegionType region_type,
                   const BoundarySegment& segment) {
    const float convergence =
        static_cast<float>(inputs.convergence_score[index]) / 100.0f;
    const float divergence =
        static_cast<float>(inputs.divergence_score[index]) / 100.0f;
    const float shear = static_cast<float>(inputs.shear_score[index]) / 100.0f;
    const float segment_convergence =
        static_cast<float>(segment.average_convergence_score) / 100.0f;
    const float segment_divergence =
        static_cast<float>(segment.average_divergence_score) / 100.0f;
    const float segment_shear =
        static_cast<float>(segment.average_shear_score) / 100.0f;

    switch (region_type) {
    case DeformingRegionType::ContinentalRift:
        return clamp_unit(0.72f * std::max(divergence, segment_divergence) +
                          0.28f * std::max(shear * 0.6f, segment_shear * 0.5f));
    case DeformingRegionType::DiffuseCollision:
        return clamp_unit(0.76f * std::max(convergence, segment_convergence) +
                          0.24f * std::max(shear * 0.5f, segment_shear * 0.45f));
    case DeformingRegionType::None:
    default:
        return 0.0f;
    }
}

float left_weight_of(const BoundarySegment& segment, uint32_t owner, DeformingRegionType region_type) {
    float left_weight = 0.5f;
    switch (region_type) {
    case DeformingRegionType::ContinentalRift:
        if (owner == segment.left_plate_id) {
            left_weight = 0.68f;
        } else if (owner == segment.right_plate_id) {
            left_weight = 0.32f;
        }
        break;
    case DeformingRegionType::DiffuseCollision:
        if (owner == segment.left_plate_id) {
            left_weight = 0.58f;
        } else if (owner == segment.right_plate_id) {
            left_weight = 0.42f;
        }
        break;
    case DeformingRegionType::None:
    default:
        break;
    }
    return left_weight;
}

// The two plates' omega x p at the cell, (east, north), weighted toward the owner.
void interpolate_velocity(const Inputs& inputs, const BoundarySegment& segment, size_t index,
                          DeformingRegionType region_type, float& velocity_x, float& velocity_y) {
    velocity_x = 0.0f;
    velocity_y = 0.0f;
    if (segment.left_plate_id >= inputs.plate_count || segment.right_plate_id >= inputs.plate_count) {
        return;
    }
    const float left_weight = left_weight_of(segment, inputs.plate_id[index], region_type);
    const float right_weight = 1.0f - left_weight;
    const mesh::Vec3 p = inputs.grid->mesh()->centre(static_cast<grid::CellId>(index));
    const auto left = boundary_cube::velocity(*inputs.grid, inputs.plate_rotation,
                                              segment.left_plate_id, inputs.plate_count, p);
    const auto right = boundary_cube::velocity(*inputs.grid, inputs.plate_rotation,
                                               segment.right_plate_id, inputs.plate_count, p);
    velocity_x = left[0] * left_weight + right[0] * right_weight;
    velocity_y = left[1] * left_weight + right[1] * right_weight;
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    const size_t cell_count = inputs.grid != nullptr ? inputs.grid->cellcount() : 0U;
    outputs.deforming_region_id.assign(cell_count, tectonics::contract::kNoDeformingRegionId);
    outputs.deforming_region_type.assign(
        cell_count, static_cast<uint8_t>(DeformingRegionType::None));
    outputs.deformation_rate.assign(cell_count, 0.0f);
    outputs.deformation_velocity_x.assign(cell_count, 0.0f);
    outputs.deformation_velocity_y.assign(cell_count, 0.0f);

    if (cell_count == 0U || inputs.plate_id == nullptr ||
        inputs.crust_class == nullptr || inputs.convergence_score == nullptr ||
        inputs.divergence_score == nullptr || inputs.shear_score == nullptr ||
        inputs.boundary_type == nullptr || inputs.geologic_regime == nullptr ||
        inputs.boundary_distance == nullptr || inputs.nearest_boundary_id == nullptr ||
        inputs.boundary_segments == nullptr || inputs.boundary_segment_count == 0U) {
        return outputs;
    }

    std::unordered_map<uint32_t, std::vector<uint32_t>> cells_by_boundary_id;
    cells_by_boundary_id.reserve(inputs.boundary_segment_count * 2U);
    for (uint32_t index = 0; index < cell_count; ++index) {
        const uint32_t boundary_id = inputs.nearest_boundary_id[index];
        if (boundary_id != tectonics::contract::kNoBoundaryId) {
            cells_by_boundary_id[boundary_id].push_back(index);
        }
    }

    for (uint32_t i = 0; i < inputs.boundary_segment_count; ++i) {
        const BoundarySegment& segment = inputs.boundary_segments[i];
        const DeformingRegionType region_type = classify_region_type(segment);
        if (region_type == DeformingRegionType::None) {
            continue;
        }

        const auto cells_it = cells_by_boundary_id.find(segment.id);
        if (cells_it == cells_by_boundary_id.end()) {
            continue;
        }

        const uint16_t max_distance = max_region_distance(region_type);
        uint32_t region_cell_count = 0U;
        mesh::Vec3 unit_sum;  // the area-weighted centre unit vectors
        double area_sum = 0.0;
        std::vector<mesh::Vec3> members;
        double deformation_total = 0.0;
        double velocity_x_total = 0.0;
        double velocity_y_total = 0.0;

        for (uint32_t index : cells_it->second) {
            const uint16_t distance = inputs.boundary_distance[index];
            if (distance > max_distance) {
                continue;
            }

            const CrustClass crust_class = crust_class_from_byte(inputs.crust_class[index]);
            if (!crust_supports_region(region_type, crust_class)) {
                continue;
            }

            const float signal = local_signal(inputs, index, region_type, segment);
            const float falloff = distance_falloff(distance, max_distance);
            const float deformation_rate = clamp_unit(signal * (0.35f + 0.65f * falloff));
            if (deformation_rate <= 0.01f) {
                continue;
            }

            float velocity_x = 0.0f;
            float velocity_y = 0.0f;
            interpolate_velocity(inputs, segment, index, region_type, velocity_x, velocity_y);

            outputs.deforming_region_id[index] = segment.id;
            outputs.deforming_region_type[index] =
                static_cast<uint8_t>(region_type);
            outputs.deformation_rate[index] = deformation_rate;
            outputs.deformation_velocity_x[index] = velocity_x;
            outputs.deformation_velocity_y[index] = velocity_y;

            unit_sum = unit_sum + inputs.grid->area_m2(index) * inputs.grid->position(index);
            area_sum += inputs.grid->area_m2(index);
            members.push_back(inputs.grid->position(index));
            deformation_total += deformation_rate;
            velocity_x_total += velocity_x;
            velocity_y_total += velocity_y;
            ++region_cell_count;
        }

        if (region_cell_count == 0U) {
            continue;
        }

        DeformingRegion region;
        region.id = segment.id;
        region.boundary_segment_id = segment.id;
        region.primary_plate_id = segment.left_plate_id;
        region.secondary_plate_id = segment.right_plate_id;
        region.cell_count = region_cell_count;
        region.persistence_steps = segment.persistence_steps;
        // A region around a great circle (a vanishing sum) takes its cell nearest the sum.
        const auto centroid = boundary_cube::encodecentroid(boundary_cube::meandirection(unit_sum, area_sum, members));
        region.centroid_x = centroid[0];
        region.centroid_y = centroid[1];
        region.average_deformation_rate =
            static_cast<float>(deformation_total / static_cast<double>(region_cell_count));
        region.average_interpolated_velocity_x =
            static_cast<float>(velocity_x_total / static_cast<double>(region_cell_count));
        region.average_interpolated_velocity_y =
            static_cast<float>(velocity_y_total / static_cast<double>(region_cell_count));
        region.average_normal_motion = segment.average_normal_motion;
        region.average_shear_motion = segment.average_shear_motion;
        region.type = region_type;
        region.age_myr = segment.age_myr;
        outputs.regions.push_back(region);
    }

    std::sort(outputs.regions.begin(), outputs.regions.end(),
              [](const DeformingRegion& lhs, const DeformingRegion& rhs) {
                  return lhs.id < rhs.id;
              });
    return outputs;
}

} // namespace tectonics::deforming_regions
