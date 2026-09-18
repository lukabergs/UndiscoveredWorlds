#include "tectonic_junctions.hpp"
#include "tectonic_kinematics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace platec::tectonic_junctions {
namespace {

using platec::contract::BoundarySegment;
using platec::contract::BoundaryType;
using platec::contract::JunctionArmType;
using platec::contract::JunctionStability;
using platec::contract::PlateKinematics;
using platec::contract::TectonicJunction;

constexpr double kTwoPi = 6.28318530717958647692;
constexpr double kPi = 3.14159265358979323846;
constexpr int kProbePadding = 2;
constexpr int kOrientationRadius = 6;
constexpr float kMinStabilityOrientationQuality = 0.25f;
constexpr float kStableConstraintError = 0.25f;

size_t map_index(uint32_t x, uint32_t y, uint32_t width) {
    return static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
}

uint32_t wrap_next(uint32_t value, uint32_t limit) {
    return value + 1U < limit ? value + 1U : 0U;
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

bool add_unique_plate(std::array<uint32_t, 4>& plates, uint32_t& count, uint32_t plate_id) {
    for (uint32_t i = 0; i < count; ++i) {
        if (plates[i] == plate_id) {
            return true;
        }
    }
    if (count >= plates.size()) {
        return false;
    }
    plates[count++] = plate_id;
    return true;
}

int pair_slot(const BoundarySegment& segment, const std::array<uint32_t, 3>& plates) {
    const uint32_t left = std::min(segment.left_plate_id, segment.right_plate_id);
    const uint32_t right = std::max(segment.left_plate_id, segment.right_plate_id);
    if (left == plates[0] && right == plates[1]) {
        return 0;
    }
    if (left == plates[0] && right == plates[2]) {
        return 1;
    }
    if (left == plates[1] && right == plates[2]) {
        return 2;
    }
    return -1;
}

uint32_t choose_segment_id(const std::map<uint32_t, uint32_t>& votes) {
    uint32_t best_id = contract::kNoBoundaryId;
    uint32_t best_votes = 0U;
    for (const auto& vote : votes) {
        if (vote.second > best_votes ||
            (vote.second == best_votes && (best_id == contract::kNoBoundaryId ||
                                           vote.first < best_id))) {
            best_id = vote.first;
            best_votes = vote.second;
        }
    }
    return best_id;
}

JunctionArmType arm_type_from_boundary_type(BoundaryType boundary_type) {
    switch (boundary_type) {
    case BoundaryType::Divergent:
        return JunctionArmType::Ridge;
    case BoundaryType::Convergent:
        return JunctionArmType::Trench;
    case BoundaryType::Transform:
        return JunctionArmType::Fault;
    case BoundaryType::PassiveMargin:
        return JunctionArmType::PassiveMargin;
    case BoundaryType::None:
    default:
        return JunctionArmType::Unknown;
    }
}

struct JunctionKey {
    std::array<uint32_t, 3> plate_ids{};
    std::array<uint32_t, 3> boundary_segment_ids{};

    bool operator<(const JunctionKey& other) const {
        if (plate_ids != other.plate_ids) {
            return plate_ids < other.plate_ids;
        }
        return boundary_segment_ids < other.boundary_segment_ids;
    }
};

struct JunctionAccumulator {
    uint32_t candidate_count = 0;
    double x_sin_sum = 0.0;
    double x_cos_sum = 0.0;
    double y_sin_sum = 0.0;
    double y_cos_sum = 0.0;
};

struct OrientationEstimate {
    float strike_x = 0.0f;
    float strike_y = 0.0f;
    float quality = 0.0f;
};

struct Vec2 {
    float x = 0.0f;
    float y = 0.0f;
};

struct ConstraintLine {
    Vec2 point;
    Vec2 direction;
};

struct ArmKinematics {
    size_t plate_a = 0U;
    size_t plate_b = 0U;
    Vec2 velocity_a;
    Vec2 velocity_b;
    Vec2 relative_velocity;
    Vec2 strike;
    Vec2 normal;
    float normal_motion = 0.0f;
    float shear_motion = 0.0f;
    float normal_magnitude = 0.0f;
    float shear_magnitude = 0.0f;
    float speed = 0.0f;
    float speed_scale = 0.0f;
};

static constexpr std::array<std::array<size_t, 2>, 3> kArmPlatePairs{{
    std::array<size_t, 2>{0U, 1U},
    std::array<size_t, 2>{0U, 2U},
    std::array<size_t, 2>{1U, 2U},
}};

float dot(Vec2 a, Vec2 b) {
    return a.x * b.x + a.y * b.y;
}

float length(Vec2 value) {
    return std::sqrt(dot(value, value));
}

Vec2 normalized(Vec2 value) {
    const float magnitude = length(value);
    if (magnitude <= 1e-6f || !std::isfinite(magnitude)) {
        return {};
    }
    return Vec2{value.x / magnitude, value.y / magnitude};
}

Vec2 local_plate_velocity(const PlateKinematics& plate, float x, float y,
                          uint32_t width, uint32_t height) {
    const auto velocity = kinematics::local_velocity(plate, x, y, width, height);
    return Vec2{velocity.x, velocity.y};
}

Vec2 subtract(Vec2 a, Vec2 b) {
    return Vec2{a.x - b.x, a.y - b.y};
}

Vec2 midpoint(Vec2 a, Vec2 b) {
    return Vec2{(a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f};
}

void add_candidate_position(JunctionAccumulator& accumulator, uint32_t x, uint32_t y,
                            uint32_t width, uint32_t height) {
    const double x_angle =
        kTwoPi * ((static_cast<double>(x) + 0.5) / static_cast<double>(width));
    const double y_angle =
        kTwoPi * ((static_cast<double>(y) + 0.5) / static_cast<double>(height));
    accumulator.x_sin_sum += std::sin(x_angle);
    accumulator.x_cos_sum += std::cos(x_angle);
    accumulator.y_sin_sum += std::sin(y_angle);
    accumulator.y_cos_sum += std::cos(y_angle);
    ++accumulator.candidate_count;
}

float circular_centroid(double sin_sum, double cos_sum, uint32_t period) {
    double angle = std::atan2(sin_sum, cos_sum);
    if (angle < 0.0) {
        angle += kTwoPi;
    }
    return static_cast<float>((angle / kTwoPi) * static_cast<double>(period));
}

OrientationEstimate estimate_orientation(const Inputs& inputs, uint32_t segment_id,
                                         float centroid_x, float centroid_y) {
    OrientationEstimate estimate;
    if (segment_id == contract::kNoBoundaryId) {
        return estimate;
    }

    const auto center_x = static_cast<int>(std::floor(centroid_x));
    const auto center_y = static_cast<int>(std::floor(centroid_y));
    double x_sum = 0.0;
    double y_sum = 0.0;
    double xx_sum = 0.0;
    double xy_sum = 0.0;
    double yy_sum = 0.0;
    uint32_t sample_count = 0U;

    for (int oy = -kOrientationRadius; oy <= kOrientationRadius; ++oy) {
        for (int ox = -kOrientationRadius; ox <= kOrientationRadius; ++ox) {
            const uint32_t x = wrapped_add(static_cast<uint32_t>(center_x), ox, inputs.width);
            const uint32_t y = wrapped_add(static_cast<uint32_t>(center_y), oy, inputs.height);
            if (inputs.boundary_segment_id[map_index(x, y, inputs.width)] != segment_id) {
                continue;
            }

            const double dx =
                static_cast<double>(wrapped_delta(static_cast<float>(x) + 0.5f, centroid_x,
                                                  inputs.width));
            const double dy =
                static_cast<double>(wrapped_delta(static_cast<float>(y) + 0.5f, centroid_y,
                                                  inputs.height));
            x_sum += dx;
            y_sum += dy;
            xx_sum += dx * dx;
            xy_sum += dx * dy;
            yy_sum += dy * dy;
            ++sample_count;
        }
    }

    if (sample_count < 2U) {
        return estimate;
    }

    const double inv_count = 1.0 / static_cast<double>(sample_count);
    const double mean_x = x_sum * inv_count;
    const double mean_y = y_sum * inv_count;
    const double cov_xx = xx_sum * inv_count - mean_x * mean_x;
    const double cov_xy = xy_sum * inv_count - mean_x * mean_y;
    const double cov_yy = yy_sum * inv_count - mean_y * mean_y;
    const double trace = cov_xx + cov_yy;
    const double discriminant =
        std::sqrt(std::max(0.0, (cov_xx - cov_yy) * (cov_xx - cov_yy) + 4.0 * cov_xy * cov_xy));
    const double lambda_max = 0.5 * (trace + discriminant);
    const double lambda_min = 0.5 * (trace - discriminant);
    if (lambda_max <= 1e-6) {
        return estimate;
    }

    const double angle = 0.5 * std::atan2(2.0 * cov_xy, cov_xx - cov_yy);
    float strike_x = static_cast<float>(std::cos(angle));
    float strike_y = static_cast<float>(std::sin(angle));
    if (strike_x < 0.0f || (std::fabs(strike_x) <= 1e-6f && strike_y < 0.0f)) {
        strike_x = -strike_x;
        strike_y = -strike_y;
    }

    const double anisotropy = (lambda_max - std::max(0.0, lambda_min)) / lambda_max;
    estimate.strike_x = strike_x;
    estimate.strike_y = strike_y;
    estimate.quality = clamp_unit(static_cast<float>(anisotropy) *
                                  clamp_unit(static_cast<float>(sample_count) / 4.0f));
    return estimate;
}

bool gather_plate_kinematics(const Inputs& inputs, const TectonicJunction& junction,
                             std::array<const PlateKinematics*, 3>& plates) {
    if (inputs.plates == nullptr || inputs.plate_count == 0U) {
        return false;
    }

    for (size_t i = 0; i < plates.size(); ++i) {
        if (junction.plate_ids[i] >= inputs.plate_count) {
            return false;
        }
        plates[i] = &inputs.plates[junction.plate_ids[i]];
        if (!std::isfinite(plates[i]->linear_velocity_x) ||
            !std::isfinite(plates[i]->linear_velocity_y) ||
            !std::isfinite(plates[i]->angular_velocity) ||
            !std::isfinite(plates[i]->mass_center_x) ||
            !std::isfinite(plates[i]->mass_center_y)) {
            return false;
        }
    }
    return true;
}

std::array<Vec2, 3> local_plate_velocities(
    const Inputs& inputs, const TectonicJunction& junction,
    const std::array<const PlateKinematics*, 3>& plates) {
    std::array<Vec2, 3> velocities{};
    for (size_t i = 0; i < velocities.size(); ++i) {
        velocities[i] =
            local_plate_velocity(*plates[i], junction.centroid_x, junction.centroid_y,
                                 inputs.width, inputs.height);
    }
    return velocities;
}

bool build_arm_kinematics(const TectonicJunction& junction,
                          const std::array<Vec2, 3>& local_velocities,
                          size_t arm_index, ArmKinematics& arm) {
    const Vec2 strike =
        normalized(Vec2{junction.arm_strike_x[arm_index], junction.arm_strike_y[arm_index]});
    if (length(strike) <= 0.0f) {
        return false;
    }

    const auto pair = kArmPlatePairs[arm_index];
    arm.plate_a = pair[0];
    arm.plate_b = pair[1];
    arm.velocity_a = local_velocities[pair[0]];
    arm.velocity_b = local_velocities[pair[1]];
    arm.relative_velocity = subtract(arm.velocity_b, arm.velocity_a);
    arm.strike = strike;
    arm.normal = Vec2{-strike.y, strike.x};
    arm.normal_motion = dot(arm.relative_velocity, arm.normal);
    arm.shear_motion = dot(arm.relative_velocity, arm.strike);
    arm.normal_magnitude = std::fabs(arm.normal_motion);
    arm.shear_magnitude = std::fabs(arm.shear_motion);
    arm.speed = length(arm.relative_velocity);
    arm.speed_scale =
        std::max(arm.speed, std::max(length(arm.velocity_a), length(arm.velocity_b)));
    return true;
}

float arm_mode_residual(const TectonicJunction& junction, const ArmKinematics& arm,
                        size_t arm_index) {
    const float weak_normal =
        std::max(0.0f, arm.speed_scale - arm.normal_magnitude);
    const float weak_shear =
        std::max(0.0f, arm.speed_scale - arm.shear_magnitude);

    switch (junction.arm_types[arm_index]) {
    case JunctionArmType::Ridge:
        return std::sqrt(arm.shear_magnitude * arm.shear_magnitude +
                         weak_normal * weak_normal);
    case JunctionArmType::Fault:
        return std::sqrt(arm.normal_magnitude * arm.normal_magnitude +
                         weak_shear * weak_shear);
    case JunctionArmType::Trench:
        return arm.shear_magnitude;
    case JunctionArmType::PassiveMargin:
    case JunctionArmType::Unknown:
    default:
        return std::max(arm.speed_scale, 1.0f);
    }
}

void populate_velocity_metrics(const Inputs& inputs, TectonicJunction& junction) {
    std::array<const PlateKinematics*, 3> plates{};
    if (!gather_plate_kinematics(inputs, junction, plates)) {
        return;
    }

    const std::array<Vec2, 3> local_velocities =
        local_plate_velocities(inputs, junction, plates);
    junction.junction_velocity_x =
        (local_velocities[0].x + local_velocities[1].x + local_velocities[2].x) / 3.0f;
    junction.junction_velocity_y =
        (local_velocities[0].y + local_velocities[1].y + local_velocities[2].y) / 3.0f;

    float residual_sum = 0.0f;
    float speed_scale = 0.0f;
    for (size_t i = 0; i < kArmPlatePairs.size(); ++i) {
        ArmKinematics arm;
        float residual = 1.0f;
        if (build_arm_kinematics(junction, local_velocities, i, arm)) {
            residual = arm_mode_residual(junction, arm, i);
            speed_scale = std::max(speed_scale, arm.speed_scale);
        }
        residual_sum += residual * residual;
    }

    const float rms_residual =
        std::sqrt(residual_sum / static_cast<float>(kArmPlatePairs.size()));
    const float normalized_residual =
        speed_scale > 1e-6f ? rms_residual / speed_scale : rms_residual;
    junction.velocity_closure_error = normalized_residual;
    junction.velocity_quality_score = 1.0f / (1.0f + normalized_residual);
}

bool has_supported_stability_arms(const TectonicJunction& junction) {
    for (JunctionArmType arm_type : junction.arm_types) {
        if (arm_type == JunctionArmType::Unknown ||
            arm_type == JunctionArmType::PassiveMargin) {
            return false;
        }
    }
    return true;
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

float minimum_orientation_quality(const TectonicJunction& junction) {
    return std::min(junction.arm_orientation_quality[0],
                    std::min(junction.arm_orientation_quality[1],
                             junction.arm_orientation_quality[2]));
}

float undirected_axis_angle(Vec2 a, Vec2 b) {
    const float cosine = clamp_unit(std::fabs(dot(normalized(a), normalized(b))));
    return std::acos(cosine);
}

float geometry_quality_score(const TectonicJunction& junction) {
    const float min_quality = minimum_orientation_quality(junction);
    if (min_quality <= 0.0f) {
        return 0.0f;
    }
    if (count_arms(junction, JunctionArmType::Ridge) != 3U) {
        return min_quality;
    }

    const Vec2 strikes[] = {
        Vec2{junction.arm_strike_x[0], junction.arm_strike_y[0]},
        Vec2{junction.arm_strike_x[1], junction.arm_strike_y[1]},
        Vec2{junction.arm_strike_x[2], junction.arm_strike_y[2]},
    };
    const float target = static_cast<float>(kPi / 3.0);
    const float error =
        (std::fabs(undirected_axis_angle(strikes[0], strikes[1]) - target) +
         std::fabs(undirected_axis_angle(strikes[0], strikes[2]) - target) +
         std::fabs(undirected_axis_angle(strikes[1], strikes[2]) - target)) /
        (3.0f * target);
    return min_quality * clamp_unit(1.0f - error);
}

bool build_constraint_line(const TectonicJunction& junction,
                           const std::array<Vec2, 3>& local_velocities,
                           size_t arm_index, ConstraintLine& line,
                           float& velocity_scale) {
    if (junction.arm_orientation_quality[arm_index] < kMinStabilityOrientationQuality) {
        return false;
    }

    const Vec2 strike =
        normalized(Vec2{junction.arm_strike_x[arm_index], junction.arm_strike_y[arm_index]});
    if (length(strike) <= 0.0f) {
        return false;
    }

    const auto pair = kArmPlatePairs[arm_index];
    const Vec2 a = local_velocities[pair[0]];
    const Vec2 b = local_velocities[pair[1]];
    velocity_scale = std::max(velocity_scale, std::max(length(a), length(b)));
    velocity_scale = std::max(velocity_scale, length(subtract(a, b)));

    line.point = midpoint(a, b);
    line.direction = strike;
    return true;
}

bool solve_constraint_intersection(const std::array<ConstraintLine, 3>& lines,
                                   Vec2& intersection,
                                   float& error) {
    double a00 = 0.0;
    double a01 = 0.0;
    double a11 = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;

    for (const ConstraintLine& line : lines) {
        const double nx = -static_cast<double>(line.direction.y);
        const double ny = static_cast<double>(line.direction.x);
        const double c = nx * static_cast<double>(line.point.x) +
                         ny * static_cast<double>(line.point.y);
        a00 += nx * nx;
        a01 += nx * ny;
        a11 += ny * ny;
        b0 += nx * c;
        b1 += ny * c;
    }

    const double determinant = a00 * a11 - a01 * a01;
    if (std::fabs(determinant) <= 1e-8) {
        return false;
    }

    const double x = (b0 * a11 - a01 * b1) / determinant;
    const double y = (a00 * b1 - a01 * b0) / determinant;
    intersection = Vec2{static_cast<float>(x), static_cast<float>(y)};
    double residual_sum = 0.0;
    for (const ConstraintLine& line : lines) {
        const double nx = -static_cast<double>(line.direction.y);
        const double ny = static_cast<double>(line.direction.x);
        const double c = nx * static_cast<double>(line.point.x) +
                         ny * static_cast<double>(line.point.y);
        const double residual = nx * x + ny * y - c;
        residual_sum += residual * residual;
    }

    error = static_cast<float>(std::sqrt(residual_sum / static_cast<double>(lines.size())));
    return std::isfinite(intersection.x) && std::isfinite(intersection.y) &&
           std::isfinite(error);
}

void populate_stability_diagnostics(const Inputs& inputs, TectonicJunction& junction) {
    junction.geometry_quality_score = geometry_quality_score(junction);

    if (count_arms(junction, JunctionArmType::Fault) == 3U) {
        junction.stability = JunctionStability::Unstable;
        junction.stability_score = 0.0f;
        return;
    }

    if (!has_supported_stability_arms(junction)) {
        junction.stability = JunctionStability::Indeterminate;
        return;
    }

    std::array<const PlateKinematics*, 3> plates{};
    if (!gather_plate_kinematics(inputs, junction, plates) ||
        minimum_orientation_quality(junction) < kMinStabilityOrientationQuality) {
        junction.stability = JunctionStability::Indeterminate;
        return;
    }

    const std::array<Vec2, 3> local_velocities =
        local_plate_velocities(inputs, junction, plates);
    std::array<ConstraintLine, 3> lines{};
    float velocity_scale = 1.0f;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (!build_constraint_line(junction, local_velocities, i, lines[i], velocity_scale)) {
            junction.stability = JunctionStability::Indeterminate;
            return;
        }
    }

    Vec2 constraint_velocity;
    if (!solve_constraint_intersection(lines, constraint_velocity,
                                       junction.velocity_constraint_error)) {
        junction.velocity_constraint_error = velocity_scale;
        junction.stability = JunctionStability::Unstable;
        return;
    }

    junction.junction_velocity_x = constraint_velocity.x;
    junction.junction_velocity_y = constraint_velocity.y;

    velocity_scale = std::max(velocity_scale, 1.0f);
    float normalized_error = junction.velocity_constraint_error / velocity_scale;
    // Intersecting velocity lines alone can accept a "transform" arm whose
    // plates move across the fault. Check relative motion at the junction too.
    for (size_t i = 0; i < lines.size(); ++i) {
        ArmKinematics arm;
        if (junction.arm_types[i] == JunctionArmType::Fault &&
            build_arm_kinematics(junction, local_velocities, i, arm) && arm.speed > 1e-6f) {
            normalized_error = std::max(normalized_error, arm.normal_magnitude / arm.speed);
        }
    }
    const float velocity_score = 1.0f / (1.0f + normalized_error);
    junction.stability_score = clamp_unit(velocity_score * junction.geometry_quality_score);

    if (normalized_error > kStableConstraintError) {
        junction.stability = JunctionStability::Unstable;
    } else if (count_arms(junction, JunctionArmType::Ridge) == 2U &&
               count_arms(junction, JunctionArmType::Fault) == 1U) {
        junction.stability = JunctionStability::ConditionallyStable;
    } else {
        junction.stability = JunctionStability::Stable;
    }
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    if (inputs.width == 0U || inputs.height == 0U || inputs.plate_id == nullptr ||
        inputs.boundary_segment_id == nullptr || inputs.boundary_segments == nullptr ||
        inputs.boundary_segment_count == 0U) {
        return outputs;
    }

    std::unordered_map<uint32_t, const BoundarySegment*> segments_by_id;
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> unique_pair_segments;
    segments_by_id.reserve(inputs.boundary_segment_count);
    for (uint32_t i = 0; i < inputs.boundary_segment_count; ++i) {
        const BoundarySegment& segment = inputs.boundary_segments[i];
        if (segment.id != contract::kNoBoundaryId) {
            segments_by_id.emplace(segment.id, &segment);
            const auto pair = std::minmax(segment.left_plate_id, segment.right_plate_id);
            const auto inserted = unique_pair_segments.emplace(pair, segment.id);
            if (!inserted.second && inserted.first->second != segment.id) {
                inserted.first->second = contract::kNoBoundaryId;
            }
        }
    }

    std::unordered_set<uint32_t> represented_segments;
    const size_t cell_count = static_cast<size_t>(inputs.width) * inputs.height;
    for (size_t i = 0; i < cell_count; ++i) {
        if (inputs.boundary_segment_id[i] != contract::kNoBoundaryId) {
            represented_segments.insert(inputs.boundary_segment_id[i]);
        }
    }

    std::map<JunctionKey, JunctionAccumulator> candidates;
    for (uint32_t y = 0; y < inputs.height; ++y) {
        for (uint32_t x = 0; x < inputs.width; ++x) {
            std::array<uint32_t, 4> neighborhood_plates{};
            uint32_t neighborhood_plate_count = 0U;
            const uint32_t x1 = wrap_next(x, inputs.width);
            const uint32_t y1 = wrap_next(y, inputs.height);
            const uint32_t indices[] = {
                static_cast<uint32_t>(map_index(x, y, inputs.width)),
                static_cast<uint32_t>(map_index(x1, y, inputs.width)),
                static_cast<uint32_t>(map_index(x, y1, inputs.width)),
                static_cast<uint32_t>(map_index(x1, y1, inputs.width)),
            };
            for (uint32_t index : indices) {
                if (!add_unique_plate(neighborhood_plates, neighborhood_plate_count,
                                      inputs.plate_id[index])) {
                    break;
                }
            }
            if (neighborhood_plate_count != 3U) {
                continue;
            }

            std::array<uint32_t, 3> plate_ids{
                neighborhood_plates[0],
                neighborhood_plates[1],
                neighborhood_plates[2],
            };
            std::sort(plate_ids.begin(), plate_ids.end());

            std::array<std::map<uint32_t, uint32_t>, 3> segment_votes;
            for (int oy = -kProbePadding; oy <= kProbePadding + 1; ++oy) {
                for (int ox = -kProbePadding; ox <= kProbePadding + 1; ++ox) {
                    const uint32_t sx = wrapped_add(x, ox, inputs.width);
                    const uint32_t sy = wrapped_add(y, oy, inputs.height);
                    const uint32_t segment_id =
                        inputs.boundary_segment_id[map_index(sx, sy, inputs.width)];
                    if (segment_id == contract::kNoBoundaryId) {
                        continue;
                    }
                    const auto found = segments_by_id.find(segment_id);
                    if (found == segments_by_id.end()) {
                        continue;
                    }
                    const int slot = pair_slot(*found->second, plate_ids);
                    if (slot >= 0) {
                        segment_votes[static_cast<size_t>(slot)][segment_id] += 1U;
                    }
                }
            }

            JunctionKey key;
            key.plate_ids = plate_ids;
            for (size_t i = 0; i < key.boundary_segment_ids.size(); ++i) {
                key.boundary_segment_ids[i] = choose_segment_id(segment_votes[i]);
            }
            // A raster cell stores one segment, so a triple contact can lose
            // one arm. Recover only a unique graph edge absent from the entire
            // raster; never borrow a represented remote edge or guess between IDs.
            if (std::count(key.boundary_segment_ids.begin(), key.boundary_segment_ids.end(),
                           contract::kNoBoundaryId) == 1) {
                for (size_t i = 0; i < key.boundary_segment_ids.size(); ++i) {
                    if (key.boundary_segment_ids[i] != contract::kNoBoundaryId) {
                        continue;
                    }
                    const auto pair = kArmPlatePairs[i];
                    const auto found = unique_pair_segments.find(
                        {plate_ids[pair[0]], plate_ids[pair[1]]});
                    if (found != unique_pair_segments.end() &&
                        represented_segments.count(found->second) == 0U) {
                        key.boundary_segment_ids[i] = found->second;
                    }
                }
            }
            if (key.boundary_segment_ids[0] == contract::kNoBoundaryId ||
                key.boundary_segment_ids[1] == contract::kNoBoundaryId ||
                key.boundary_segment_ids[2] == contract::kNoBoundaryId) {
                continue;
            }

            add_candidate_position(candidates[key], x, y, inputs.width, inputs.height);
        }
    }

    outputs.junctions.reserve(candidates.size());
    uint32_t next_id = contract::kNoTectonicJunctionId + 1U;
    for (const auto& candidate : candidates) {
        TectonicJunction junction;
        junction.id = next_id++;
        junction.candidate_cell_count = candidate.second.candidate_count;
        junction.centroid_x =
            circular_centroid(candidate.second.x_sin_sum, candidate.second.x_cos_sum, inputs.width);
        junction.centroid_y =
            circular_centroid(candidate.second.y_sin_sum, candidate.second.y_cos_sum, inputs.height);
        junction.quality_score =
            std::min(1.0f, static_cast<float>(candidate.second.candidate_count) / 4.0f);

        for (size_t i = 0; i < candidate.first.plate_ids.size(); ++i) {
            junction.plate_ids[i] = candidate.first.plate_ids[i];
            junction.boundary_segment_ids[i] = candidate.first.boundary_segment_ids[i];
            const auto found = segments_by_id.find(junction.boundary_segment_ids[i]);
            if (found != segments_by_id.end()) {
                junction.arm_types[i] =
                    arm_type_from_boundary_type(found->second->boundary_type);
            }
            const OrientationEstimate orientation = estimate_orientation(
                inputs, junction.boundary_segment_ids[i], junction.centroid_x, junction.centroid_y);
            junction.arm_strike_x[i] = orientation.strike_x;
            junction.arm_strike_y[i] = orientation.strike_y;
            junction.arm_orientation_quality[i] = orientation.quality;
        }
        populate_velocity_metrics(inputs, junction);
        populate_stability_diagnostics(inputs, junction);

        outputs.junctions.push_back(junction);
    }

    return outputs;
}

} // namespace platec::tectonic_junctions
