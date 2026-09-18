#include "tectonic_deforming_regions.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace platec::deforming_regions {

namespace {

using platec::contract::BoundarySegment;
using platec::contract::BoundaryType;
using platec::contract::CrustClass;
using platec::contract::DeformingRegion;
using platec::contract::DeformingRegionType;
using platec::contract::GeologicRegime;
using platec::contract::PlateKinematics;

constexpr double kTwoPi = 6.28318530717958647692;

size_t map_index(uint32_t x, uint32_t y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
}

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

void interpolate_velocity(const BoundarySegment& segment, const PlateKinematics* plates,
                          uint32_t plate_count, uint32_t owner, DeformingRegionType region_type,
                          float& velocity_x, float& velocity_y) {
    velocity_x = 0.0f;
    velocity_y = 0.0f;
    if (plates == nullptr || segment.left_plate_id >= plate_count ||
        segment.right_plate_id >= plate_count) {
        return;
    }

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

    const float right_weight = 1.0f - left_weight;
    velocity_x = plates[segment.left_plate_id].linear_velocity_x * left_weight +
                 plates[segment.right_plate_id].linear_velocity_x * right_weight;
    velocity_y = plates[segment.left_plate_id].linear_velocity_y * left_weight +
                 plates[segment.right_plate_id].linear_velocity_y * right_weight;
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    const size_t cell_count =
        static_cast<size_t>(inputs.width) * static_cast<size_t>(inputs.height);
    outputs.deforming_region_id.assign(cell_count, platec::contract::kNoDeformingRegionId);
    outputs.deforming_region_type.assign(
        cell_count, static_cast<uint8_t>(DeformingRegionType::None));
    outputs.deformation_rate.assign(cell_count, 0.0f);
    outputs.deformation_velocity_x.assign(cell_count, 0.0f);
    outputs.deformation_velocity_y.assign(cell_count, 0.0f);

    if (inputs.width == 0U || inputs.height == 0U || inputs.plate_id == nullptr ||
        inputs.crust_class == nullptr || inputs.convergence_score == nullptr ||
        inputs.divergence_score == nullptr || inputs.shear_score == nullptr ||
        inputs.boundary_type == nullptr || inputs.geologic_regime == nullptr ||
        inputs.boundary_distance == nullptr || inputs.nearest_boundary_id == nullptr ||
        inputs.boundary_segments == nullptr || inputs.boundary_segment_count == 0U) {
        return outputs;
    }

    std::unordered_map<uint32_t, std::vector<uint32_t>> cells_by_boundary_id;
    cells_by_boundary_id.reserve(inputs.boundary_segment_count * 2U);
    for (uint32_t y = 0; y < inputs.height; ++y) {
        for (uint32_t x = 0; x < inputs.width; ++x) {
            const uint32_t index =
                static_cast<uint32_t>(map_index(x, y, inputs.width));
            const uint32_t boundary_id = inputs.nearest_boundary_id[index];
            if (boundary_id != platec::contract::kNoBoundaryId) {
                cells_by_boundary_id[boundary_id].push_back(index);
            }
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
        double x_sin_sum = 0.0;
        double x_cos_sum = 0.0;
        double y_sin_sum = 0.0;
        double y_cos_sum = 0.0;
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
            interpolate_velocity(segment, inputs.plates, inputs.plate_count,
                                 inputs.plate_id[index], region_type, velocity_x, velocity_y);

            outputs.deforming_region_id[index] = segment.id;
            outputs.deforming_region_type[index] =
                static_cast<uint8_t>(region_type);
            outputs.deformation_rate[index] = deformation_rate;
            outputs.deformation_velocity_x[index] = velocity_x;
            outputs.deformation_velocity_y[index] = velocity_y;

            const uint32_t x = index % inputs.width;
            const uint32_t y = index / inputs.width;
            const double x_angle =
                kTwoPi * (static_cast<double>(x) / static_cast<double>(inputs.width));
            const double y_angle =
                kTwoPi * (static_cast<double>(y) / static_cast<double>(inputs.height));
            x_sin_sum += std::sin(x_angle);
            x_cos_sum += std::cos(x_angle);
            y_sin_sum += std::sin(y_angle);
            y_cos_sum += std::cos(y_angle);
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
        region.centroid_x = static_cast<float>(
            ((std::atan2(x_sin_sum, x_cos_sum) < 0.0
                  ? std::atan2(x_sin_sum, x_cos_sum) + kTwoPi
                  : std::atan2(x_sin_sum, x_cos_sum)) /
             kTwoPi) *
            static_cast<double>(inputs.width));
        region.centroid_y = static_cast<float>(
            ((std::atan2(y_sin_sum, y_cos_sum) < 0.0
                  ? std::atan2(y_sin_sum, y_cos_sum) + kTwoPi
                  : std::atan2(y_sin_sum, y_cos_sum)) /
             kTwoPi) *
            static_cast<double>(inputs.height));
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

} // namespace platec::deforming_regions
