#include "simulations/geology/plate_tectonics/tectonic_boundary_graph.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"

#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace tectonics::boundary_graph {

namespace {

using tectonics::contract::BoundarySegment;
using tectonics::contract::BoundaryType;
using tectonics::contract::CrustClass;
using tectonics::contract::GeologicRegime;
using tectonics::grid::TectonicGrid;

constexpr double kTwoPi = 6.28318530717958647692;

uint8_t average_score(uint8_t a, uint8_t b) {
    return static_cast<uint8_t>((static_cast<uint32_t>(a) + static_cast<uint32_t>(b) + 1U) / 2U);
}

float wrapped_distance(float a, float b, float period) {
    float delta = a - b;
    const float half_period = period * 0.5f;
    if (delta > half_period) {
        delta -= period;
    } else if (delta < -half_period) {
        delta += period;
    }
    return delta;
}

// Between two segment centroids in cells: the wrapped offset on the scaffold, the great circle
// over the cube's cell spacing on the cube (centroids as longitude and latitude).
float centroid_distance(const TectonicGrid& grid, float x0, float y0, float x1, float y1) {
    if (grid.kind() == TectonicGrid::Kind::cube) {
        const double arc = mesh::greatcircle_rad(boundary_cube::decodecentroid(x0, y0),
                                                 boundary_cube::decodecentroid(x1, y1));
        return static_cast<float>(arc * grid.mesh()->radius_m() / boundary_cube::cellspacing_m(grid));
    }
    const float dx = wrapped_distance(x0, x1, static_cast<float>(grid.width()));
    const float dy = wrapped_distance(y0, y1, static_cast<float>(grid.height()));
    return std::sqrt(dx * dx + dy * dy);
}

GeologicRegime regime_from_byte(const uint8_t value) {
    return static_cast<GeologicRegime>(value);
}

CrustClass crust_class_from_byte(const uint8_t value) {
    return static_cast<CrustClass>(value);
}

BoundaryType classify_boundary_type(GeologicRegime raw_regime, uint8_t convergence_score,
                                    uint8_t divergence_score, uint8_t shear_score,
                                    CrustClass crust_a, CrustClass crust_b) {
    switch (raw_regime) {
    case GeologicRegime::ConvergentArc:
    case GeologicRegime::ContinentCollision:
    case GeologicRegime::TrenchAdjacent:
        return BoundaryType::Convergent;
    case GeologicRegime::DivergentRift:
    case GeologicRegime::MidOceanRidge:
        return BoundaryType::Divergent;
    case GeologicRegime::Transform:
        return BoundaryType::Transform;
    case GeologicRegime::PassiveMargin:
        return BoundaryType::PassiveMargin;
    case GeologicRegime::Stable:
    default:
        break;
    }

    const uint8_t strongest =
        std::max({convergence_score, divergence_score, shear_score});
    const bool continental_mix =
        (crust_a == CrustClass::Continental && crust_b == CrustClass::Oceanic) ||
        (crust_a == CrustClass::Oceanic && crust_b == CrustClass::Continental) ||
        (crust_a == CrustClass::Continental && crust_b == CrustClass::Transitional) ||
        (crust_a == CrustClass::Transitional && crust_b == CrustClass::Continental);

    if (strongest <= 20U && continental_mix) {
        return BoundaryType::PassiveMargin;
    }
    if (strongest < 15U) {
        return BoundaryType::PassiveMargin;
    }
    if (convergence_score >= divergence_score && convergence_score >= shear_score) {
        return BoundaryType::Convergent;
    }
    if (divergence_score >= shear_score) {
        return BoundaryType::Divergent;
    }
    return BoundaryType::Transform;
}

GeologicRegime choose_regime_hint(BoundaryType boundary_type, CrustClass crust_a,
                                  CrustClass crust_b) {
    switch (boundary_type) {
    case BoundaryType::Convergent:
        if ((crust_a == CrustClass::Continental || crust_a == CrustClass::Transitional) &&
            (crust_b == CrustClass::Continental || crust_b == CrustClass::Transitional)) {
            return GeologicRegime::ContinentCollision;
        }
        if (crust_a == CrustClass::Oceanic &&
            (crust_b == CrustClass::Continental || crust_b == CrustClass::Transitional)) {
            return GeologicRegime::TrenchAdjacent;
        }
        return GeologicRegime::ConvergentArc;
    case BoundaryType::Divergent:
        if (crust_a == CrustClass::Oceanic && crust_b == CrustClass::Oceanic) {
            return GeologicRegime::MidOceanRidge;
        }
        return GeologicRegime::DivergentRift;
    case BoundaryType::Transform:
        return GeologicRegime::Transform;
    case BoundaryType::PassiveMargin:
        return GeologicRegime::PassiveMargin;
    case BoundaryType::None:
    default:
        return GeologicRegime::Stable;
    }
}

GeologicRegime finalize_segment_regime(BoundaryType boundary_type,
                                       const std::array<uint32_t, 8>& regime_counts,
                                       uint32_t continental_votes, uint32_t oceanic_votes) {
    auto count = [&](GeologicRegime regime) -> uint32_t {
        const size_t index = static_cast<size_t>(regime);
        return index < regime_counts.size() ? regime_counts[index] : 0U;
    };

    switch (boundary_type) {
    case BoundaryType::Convergent:
        if (count(GeologicRegime::ContinentCollision) > 0U &&
            count(GeologicRegime::ContinentCollision) >=
                std::max(count(GeologicRegime::ConvergentArc),
                         count(GeologicRegime::TrenchAdjacent))) {
            return GeologicRegime::ContinentCollision;
        }
        if (count(GeologicRegime::TrenchAdjacent) > count(GeologicRegime::ConvergentArc) &&
            oceanic_votes >= continental_votes) {
            return GeologicRegime::TrenchAdjacent;
        }
        return GeologicRegime::ConvergentArc;
    case BoundaryType::Divergent:
        if (count(GeologicRegime::MidOceanRidge) >= count(GeologicRegime::DivergentRift) &&
            oceanic_votes >= continental_votes) {
            return GeologicRegime::MidOceanRidge;
        }
        return GeologicRegime::DivergentRift;
    case BoundaryType::Transform:
        return GeologicRegime::Transform;
    case BoundaryType::PassiveMargin:
        return GeologicRegime::PassiveMargin;
    case BoundaryType::None:
    default:
        return GeologicRegime::Stable;
    }
}

struct CellSeed {
    bool active = false;
    uint32_t left_plate = 0;
    uint32_t right_plate = 0;
    BoundaryType boundary_type = BoundaryType::None;
    GeologicRegime regime_hint = GeologicRegime::Stable;
    uint8_t convergence_score = 0;
    uint8_t divergence_score = 0;
    uint8_t shear_score = 0;
    float normal_motion = 0.0f;
    float shear_motion = 0.0f;
};

CellSeed build_seed_for_cell(const Inputs& inputs, const TectonicGrid& grid, size_t index) {
    CellSeed best;
    const uint32_t plate = inputs.plate_id[index];
    // West, east, north, south: the torus stencil's x - 1, x + 1, y - 1, y + 1.
    const auto neighbours = grid.neighbours(static_cast<grid::CellId>(index));

    for (size_t i = 0; i < 4; ++i) {
        const size_t neighbor_index = neighbours[i];
        const uint32_t neighbor_plate = inputs.plate_id[neighbor_index];
        if (plate == neighbor_plate) {
            continue;
        }

        const uint8_t convergence_score =
            average_score(inputs.convergence_score[index], inputs.convergence_score[neighbor_index]);
        const uint8_t divergence_score =
            average_score(inputs.divergence_score[index], inputs.divergence_score[neighbor_index]);
        const uint8_t shear_score =
            average_score(inputs.shear_score[index], inputs.shear_score[neighbor_index]);
        const GeologicRegime regime_a = regime_from_byte(inputs.raw_geologic_regime[index]);
        const GeologicRegime regime_b =
            regime_from_byte(inputs.raw_geologic_regime[neighbor_index]);
        const CrustClass crust_a = crust_class_from_byte(inputs.crust_class[index]);
        const CrustClass crust_b = crust_class_from_byte(inputs.crust_class[neighbor_index]);
        const GeologicRegime raw_regime =
            regime_a != GeologicRegime::Stable ? regime_a : regime_b;
        const BoundaryType boundary_type = classify_boundary_type(
            raw_regime, convergence_score, divergence_score, shear_score, crust_a, crust_b);
        if (boundary_type == BoundaryType::None) {
            continue;
        }

        const float normal_motion =
            (static_cast<float>(convergence_score) - static_cast<float>(divergence_score)) /
            100.0f;
        const float shear_motion = static_cast<float>(shear_score) / 100.0f;
        const float strength = static_cast<float>(
            std::max({convergence_score, divergence_score, shear_score}));
        const float best_strength = static_cast<float>(std::max(
            {best.convergence_score, best.divergence_score, best.shear_score}));
        if (!best.active || strength > best_strength + 1.0f ||
            (std::fabs(strength - best_strength) <= 1.0f &&
             static_cast<uint8_t>(boundary_type) < static_cast<uint8_t>(best.boundary_type))) {
            best.active = true;
            best.left_plate = std::min(plate, neighbor_plate);
            best.right_plate = std::max(plate, neighbor_plate);
            best.boundary_type = boundary_type;
            best.regime_hint =
                choose_regime_hint(boundary_type, crust_a, crust_b);
            best.convergence_score = convergence_score;
            best.divergence_score = divergence_score;
            best.shear_score = shear_score;
            // Positive means convergence on both sides, independent of orientation.
            best.normal_motion = normal_motion;
            best.shear_motion = shear_motion;
        }
    }

    return best;
}

const BoundarySegment* find_previous_segment(
    const std::unordered_map<uint32_t, const BoundarySegment*>& previous_by_id, uint32_t id) {
    const auto it = previous_by_id.find(id);
    return it == previous_by_id.end() ? nullptr : it->second;
}

uint32_t choose_persistent_id(
    const Inputs& inputs, const TectonicGrid& grid,
    const std::unordered_map<uint32_t, const BoundarySegment*>& previous_by_id,
    const std::vector<uint32_t>& component_cells, BoundaryType boundary_type, float centroid_x,
    float centroid_y, const std::unordered_set<uint32_t>& claimed_ids) {
    std::unordered_map<uint32_t, uint32_t> overlap_votes;
    std::unordered_map<uint32_t, uint32_t> nearby_votes;

    for (uint32_t cell_index : component_cells) {
        const uint32_t overlap_id = inputs.previous_boundary_segment_id[cell_index];
        if (overlap_id != contract::kNoBoundaryId) {
            overlap_votes[overlap_id] += 1U;
        }
        const uint32_t nearby_id = inputs.previous_nearest_boundary_id[cell_index];
        if (nearby_id != contract::kNoBoundaryId) {
            nearby_votes[nearby_id] += 1U;
        }
    }

    auto compatible_previous_id =
        [&](const std::unordered_map<uint32_t, uint32_t>& votes, uint32_t min_votes,
            float max_distance) -> uint32_t {
        uint32_t best_id = contract::kNoBoundaryId;
        uint32_t best_votes = 0U;
        for (const auto& entry : votes) {
            if (entry.second < min_votes || claimed_ids.find(entry.first) != claimed_ids.end()) {
                continue;
            }
            const BoundarySegment* previous = find_previous_segment(previous_by_id, entry.first);
            if (previous == nullptr || previous->boundary_type != boundary_type) {
                continue;
            }
            if (centroid_distance(grid, centroid_x, centroid_y, previous->centroid_x,
                                  previous->centroid_y) > max_distance) {
                continue;
            }
            if (entry.second > best_votes) {
                best_votes = entry.second;
                best_id = entry.first;
            }
        }
        return best_id;
    };

    const uint32_t overlap_id = compatible_previous_id(overlap_votes, 1U, 12.0f);
    if (overlap_id != contract::kNoBoundaryId) {
        return overlap_id;
    }

    const uint32_t nearby_threshold =
        std::max<uint32_t>(2U, static_cast<uint32_t>(component_cells.size() / 8U));
    return compatible_previous_id(nearby_votes, nearby_threshold, 8.0f);
}

} // namespace

Outputs build(const Inputs& inputs) {
    Outputs outputs;
    outputs.next_boundary_id = inputs.next_boundary_id;

    std::shared_ptr<const TectonicGrid> owned;
    const TectonicGrid* grid = inputs.grid;
    if (grid == nullptr && inputs.width != 0U && inputs.height != 0U) {
        owned = TectonicGrid::planar(inputs.width, inputs.height, static_cast<double>(inputs.width));
        grid = owned.get();
    }
    if (grid == nullptr) {
        return outputs;
    }
    const bool cube = grid->kind() == TectonicGrid::Kind::cube;

    const size_t cell_count = grid->cellcount();
    outputs.boundary_segment_id.assign(cell_count, contract::kNoBoundaryId);
    outputs.boundary_type.assign(cell_count, static_cast<uint8_t>(BoundaryType::None));
    outputs.boundary_regime.assign(cell_count, static_cast<uint8_t>(GeologicRegime::Stable));

    if (inputs.plate_id == nullptr ||
        inputs.convergence_score == nullptr || inputs.divergence_score == nullptr ||
        inputs.shear_score == nullptr || inputs.raw_geologic_regime == nullptr ||
        inputs.crust_class == nullptr || inputs.previous_boundary_segment_id == nullptr ||
        inputs.previous_nearest_boundary_id == nullptr) {
        return outputs;
    }

    std::unordered_map<uint32_t, const BoundarySegment*> previous_by_id;
    previous_by_id.reserve(inputs.previous_segment_count);
    for (uint32_t i = 0; i < inputs.previous_segment_count; ++i) {
        previous_by_id.emplace(inputs.previous_segments[i].id, &inputs.previous_segments[i]);
    }

    std::vector<CellSeed> seeds(cell_count);
    for (size_t index = 0; index < cell_count; ++index) {
        seeds[index] = build_seed_for_cell(inputs, *grid, index);
    }

    std::vector<uint8_t> visited(cell_count, 0U);
    std::unordered_set<uint32_t> claimed_previous_ids;

    for (size_t start_index = 0; start_index < cell_count; ++start_index) {
        const CellSeed& start_seed = seeds[start_index];
        if (!start_seed.active || visited[start_index]) {
            continue;
        }

        std::vector<uint32_t> frontier{static_cast<uint32_t>(start_index)};
        std::vector<uint32_t> component_cells;
        component_cells.reserve(32);
        std::array<uint32_t, 8> regime_counts{};
        uint32_t continental_votes = 0U;
        uint32_t oceanic_votes = 0U;
        uint64_t convergence_total = 0U;
        uint64_t divergence_total = 0U;
        uint64_t shear_total = 0U;
        double normal_total = 0.0;
        double shear_total_f = 0.0;
        // The scaffold's circular means of columns and rows; the cube's area-weighted sum of
        // centre unit vectors.
        double x_sin_sum = 0.0;
        double x_cos_sum = 0.0;
        double y_sin_sum = 0.0;
        double y_cos_sum = 0.0;
        mesh::Vec3 unit_sum;

        while (!frontier.empty()) {
            const uint32_t current = frontier.back();
            frontier.pop_back();
            if (visited[current]) {
                continue;
            }

            const CellSeed& current_seed = seeds[current];
            if (!current_seed.active || current_seed.left_plate != start_seed.left_plate ||
                current_seed.right_plate != start_seed.right_plate ||
                current_seed.boundary_type != start_seed.boundary_type) {
                continue;
            }

            visited[current] = 1U;
            component_cells.push_back(current);
            convergence_total += current_seed.convergence_score;
            divergence_total += current_seed.divergence_score;
            shear_total += current_seed.shear_score;
            normal_total += static_cast<double>(current_seed.normal_motion);
            shear_total_f += static_cast<double>(current_seed.shear_motion);

            if (cube) {
                unit_sum = unit_sum + grid->area_m2(current) * grid->mesh()->centre(current);
            } else {
                const uint32_t cx = static_cast<uint32_t>(current % grid->width());
                const uint32_t cy = static_cast<uint32_t>(current / grid->width());
                const double x_angle =
                    kTwoPi * (static_cast<double>(cx) / static_cast<double>(grid->width()));
                const double y_angle =
                    kTwoPi * (static_cast<double>(cy) / static_cast<double>(grid->height()));
                x_sin_sum += std::sin(x_angle);
                x_cos_sum += std::cos(x_angle);
                y_sin_sum += std::sin(y_angle);
                y_cos_sum += std::cos(y_angle);
            }

            const GeologicRegime hint = current_seed.regime_hint;
            const size_t regime_index = static_cast<size_t>(hint);
            if (regime_index < regime_counts.size()) {
                regime_counts[regime_index] += 1U;
            }

            const CrustClass crust = crust_class_from_byte(inputs.crust_class[current]);
            if (crust == CrustClass::Oceanic) {
                oceanic_votes += 1U;
            } else if (crust == CrustClass::Continental ||
                       crust == CrustClass::Transitional) {
                continental_votes += 1U;
            }

            const auto neighbors = grid->neighbours(current);
            frontier.insert(frontier.end(), neighbors.begin(), neighbors.end());
        }

        if (component_cells.empty()) {
            continue;
        }

        float centroid_x = 0.0f;
        float centroid_y = 0.0f;
        if (cube) {
            const auto centroid = boundary_cube::encodecentroid(mesh::normalise(unit_sum));
            centroid_x = centroid[0];
            centroid_y = centroid[1];
        } else {
            const double x_angle = std::atan2(x_sin_sum, x_cos_sum);
            const double y_angle = std::atan2(y_sin_sum, y_cos_sum);
            centroid_x = static_cast<float>(
                (x_angle < 0.0 ? x_angle + kTwoPi : x_angle) / kTwoPi *
                static_cast<double>(grid->width()));
            centroid_y = static_cast<float>(
                (y_angle < 0.0 ? y_angle + kTwoPi : y_angle) / kTwoPi *
                static_cast<double>(grid->height()));
        }
        const uint32_t matched_id = choose_persistent_id(
            inputs, *grid, previous_by_id, component_cells, start_seed.boundary_type, centroid_x,
            centroid_y, claimed_previous_ids);

        BoundarySegment segment;
        segment.id = matched_id != contract::kNoBoundaryId ? matched_id : outputs.next_boundary_id++;
        segment.left_plate_id = start_seed.left_plate;
        segment.right_plate_id = start_seed.right_plate;
        segment.cell_count = static_cast<uint32_t>(component_cells.size());
        segment.centroid_x = centroid_x;
        segment.centroid_y = centroid_y;
        segment.length_cells = static_cast<float>(component_cells.size());
        segment.average_convergence_score = static_cast<uint8_t>(
            convergence_total / static_cast<uint64_t>(component_cells.size()));
        segment.average_divergence_score = static_cast<uint8_t>(
            divergence_total / static_cast<uint64_t>(component_cells.size()));
        segment.average_shear_score = static_cast<uint8_t>(
            shear_total / static_cast<uint64_t>(component_cells.size()));
        segment.average_normal_motion = static_cast<float>(
            normal_total / static_cast<double>(component_cells.size()));
        segment.average_shear_motion = static_cast<float>(
            shear_total_f / static_cast<double>(component_cells.size()));
        segment.boundary_type = start_seed.boundary_type;
        segment.geologic_regime = finalize_segment_regime(
            start_seed.boundary_type, regime_counts, continental_votes, oceanic_votes);

        const BoundarySegment* previous = matched_id != contract::kNoBoundaryId
                                              ? find_previous_segment(previous_by_id, matched_id)
                                              : nullptr;
        if (previous != nullptr) {
            segment.persistence_steps = previous->persistence_steps + 1U;
            segment.age_myr = previous->age_myr + inputs.delta_time_myr;
            claimed_previous_ids.insert(previous->id);
        } else {
            segment.persistence_steps = 1U;
            segment.age_myr = inputs.delta_time_myr;
        }

        outputs.segments.push_back(segment);
        for (uint32_t cell_index : component_cells) {
            outputs.boundary_segment_id[cell_index] = segment.id;
            outputs.boundary_type[cell_index] =
                static_cast<uint8_t>(segment.boundary_type);
            outputs.boundary_regime[cell_index] =
                static_cast<uint8_t>(seeds[cell_index].regime_hint);
        }
    }

    std::sort(outputs.segments.begin(), outputs.segments.end(),
              [](const BoundarySegment& lhs, const BoundarySegment& rhs) {
                  return lhs.id < rhs.id;
              });
    return outputs;
}

} // namespace tectonics::boundary_graph
