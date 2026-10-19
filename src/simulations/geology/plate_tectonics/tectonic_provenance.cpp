#include "simulations/geology/plate_tectonics/tectonic_provenance.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include "simulations/geology/plate_tectonics/topography_codec.hpp"

#include <algorithm>
#include <cmath>

namespace tectonics::provenance {

namespace {

using grid::TectonicGrid;

float clamp01(float value) {
    return std::clamp(value, 0.0f, 1.0f);
}

bool is_sea(const Inputs& inputs, size_t index) {
    return TopographyCodec::is_oceanic_internal(inputs.heightmap[index]);
}

bool is_continental(const Inputs& inputs, size_t index) {
    if (inputs.crust_class == nullptr) {
        return !TopographyCodec::is_oceanic_internal(inputs.heightmap[index]);
    }
    const auto crust = static_cast<contract::CrustClass>(inputs.crust_class[index]);
    return crust == contract::CrustClass::Continental ||
           crust == contract::CrustClass::Transitional;
}

// Any of the eight neighbours (the mesh's slots 1..8; a cube corner cell has seven) is sea
// (want_sea) or land.
bool has_adjacent(const Inputs& inputs, const TectonicGrid& grid, size_t index, bool want_sea) {
    for (int slot = 1; slot <= mesh::slotcount; ++slot) {
        const std::int32_t n = grid.mesh()->neighbour(static_cast<grid::CellId>(index), slot);
        if (n >= 0 && is_sea(inputs, static_cast<size_t>(n)) == want_sea) return true;
    }
    return false;
}

bool is_coast_or_outline(const Inputs& inputs, const TectonicGrid& grid, size_t index) {
    return is_sea(inputs, index) ? has_adjacent(inputs, grid, index, false)
                                 : has_adjacent(inputs, grid, index, true);
}

uint8_t normalize_tectonic_signal(float signal, double velocity_to_cm_yr) {
    if (velocity_to_cm_yr > 0.0) {
        // A 0..100 activity score over 0..5 cm/year, not an update displacement.
        return static_cast<uint8_t>(std::lround(std::clamp(signal * velocity_to_cm_yr / 5.0, 0.0, 1.0) * 100.0));
    }
    const float normalized = clamp01(
        (signal - contract::kConvergentBoundaryThreshold) /
        (2.0f - contract::kConvergentBoundaryThreshold));
    return static_cast<uint8_t>(std::lround(normalized * 100.0f));
}

contract::GeologicRegime classify_geologic_regime(const Inputs& inputs, const TectonicGrid& grid,
                                                  size_t index, uint8_t convergence,
                                                  uint8_t divergence, uint8_t shear,
                                                  contract::GeologicRegime convergent_regime) {
    const uint8_t max_signal = std::max({convergence, divergence, shear});

    if (inputs.velocity_to_cm_yr > 0.0 && max_signal > 0U) {
        if (convergence >= divergence && convergence >= shear) return convergent_regime;
        if (divergence >= shear)
            return is_continental(inputs, index)
                ? contract::GeologicRegime::DivergentRift : contract::GeologicRegime::MidOceanRidge;
        return contract::GeologicRegime::Transform;
    }

    if (convergence >= 65U) {
        return convergent_regime;
    }

    if (divergence >= 55U) {
        return is_continental(inputs, index)
                   ? contract::GeologicRegime::DivergentRift
                   : contract::GeologicRegime::MidOceanRidge;
    }

    if (shear >= 50U && shear >= static_cast<uint8_t>(convergence + 10U) &&
        shear >= static_cast<uint8_t>(divergence + 10U)) {
        return contract::GeologicRegime::Transform;
    }

    if (is_coast_or_outline(inputs, grid, index) && max_signal <= 20U) {
        return contract::GeologicRegime::PassiveMargin;
    }

    return contract::GeologicRegime::Stable;
}

struct Signals {
    std::vector<float> convergence, divergence, shear;
    std::vector<contract::GeologicRegime> convergent_regime;
};

// One contact's normal (positive: converging) and tangential motion, on both its cells.
void record_motion(const Inputs& inputs, size_t index, size_t neighbor_index,
                   float boundary_motion, float tangential_motion, Signals& signals) {
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
        if (boundary_motion > signals.convergence[index])
            signals.convergent_regime[index] = regime(first_continent, second_continent);
        if (boundary_motion > signals.convergence[neighbor_index])
            signals.convergent_regime[neighbor_index] = regime(second_continent, first_continent);
        signals.convergence[index] = std::max(signals.convergence[index], boundary_motion);
        signals.convergence[neighbor_index] = std::max(signals.convergence[neighbor_index], boundary_motion);
    } else if (boundary_motion < 0.0f) {
        const float divergence_signal = -boundary_motion;
        signals.divergence[index] = std::max(signals.divergence[index], divergence_signal);
        signals.divergence[neighbor_index] = std::max(signals.divergence[neighbor_index], divergence_signal);
    }

    signals.shear[index] = std::max(signals.shear[index], tangential_motion);
    signals.shear[neighbor_index] = std::max(signals.shear[neighbor_index], tangential_motion);
}

bool is_contact(const Inputs& inputs, size_t index, size_t neighbor_index) {
    const uint32_t plate = inputs.plate_id[index];
    const uint32_t neighbor_plate = inputs.plate_id[neighbor_index];
    return plate != neighbor_plate && plate < inputs.plate_count &&
           neighbor_plate < inputs.plate_count;
}

// Both plates' omega x p at the shared edge's midpoint, against the edge's normal.
void record_edge(const Inputs& inputs, const TectonicGrid& grid, const grid::Edge& edge, Signals& signals) {
    const size_t index = edge.a;
    const size_t neighbor_index = edge.b;
    if (!is_contact(inputs, index, neighbor_index)) return;
    const mesh::Vec3 p = edge.midpoint;
    const auto first = boundary_cube::velocity3(grid, inputs.plate_rotation, inputs.plate_id[index],
                                                inputs.plate_count, p);
    const auto second = boundary_cube::velocity3(grid, inputs.plate_rotation,
                                                 inputs.plate_id[neighbor_index], inputs.plate_count, p);
    const mesh::Vec3 relative{first[0] - second[0], first[1] - second[1], first[2] - second[2]};
    if (!std::isfinite(relative.x) || !std::isfinite(relative.y) || !std::isfinite(relative.z)) return;
    const mesh::Vec3 strike = mesh::cross(p, edge.normal);
    const float boundary_motion = static_cast<float>(mesh::dot(relative, edge.normal));
    const float tangential_motion = static_cast<float>(std::fabs(mesh::dot(relative, strike)));
    record_motion(inputs, index, neighbor_index, boundary_motion, tangential_motion, signals);
}

} // namespace

void compute_maps(const Inputs& inputs, uint8_t* convergence_score,
                  uint8_t* divergence_score, uint8_t* shear_score,
                  uint8_t* geologic_regime) {
    const TectonicGrid* grid = inputs.grid;
    if (grid == nullptr || inputs.heightmap == nullptr ||
        inputs.plate_id == nullptr || convergence_score == nullptr ||
        divergence_score == nullptr || shear_score == nullptr ||
        geologic_regime == nullptr) {
        return;
    }

    const size_t cell_count = grid->cellcount();
    std::fill(convergence_score, convergence_score + cell_count, static_cast<uint8_t>(0));
    std::fill(divergence_score, divergence_score + cell_count, static_cast<uint8_t>(0));
    std::fill(shear_score, shear_score + cell_count, static_cast<uint8_t>(0));
    std::fill(geologic_regime, geologic_regime + cell_count,
              static_cast<uint8_t>(contract::GeologicRegime::Stable));

    if (inputs.plate_count == 0) {
        return;
    }

    Signals signals;
    signals.convergence.assign(cell_count, 0.0f);
    signals.divergence.assign(cell_count, 0.0f);
    signals.shear.assign(cell_count, 0.0f);
    signals.convergent_regime.assign(cell_count, contract::GeologicRegime::Stable);
    // Every contact once.
    for (const grid::Edge& edge : grid->edges()) record_edge(inputs, *grid, edge, signals);

    for (size_t index = 0; index < cell_count; ++index) {
        convergence_score[index] = normalize_tectonic_signal(signals.convergence[index], inputs.velocity_to_cm_yr);
        divergence_score[index] = normalize_tectonic_signal(signals.divergence[index], inputs.velocity_to_cm_yr);
        shear_score[index] = normalize_tectonic_signal(signals.shear[index], inputs.velocity_to_cm_yr);
        geologic_regime[index] = static_cast<uint8_t>(classify_geologic_regime(
            inputs, *grid, index, convergence_score[index], divergence_score[index],
            shear_score[index], signals.convergent_regime[index]));
    }
}

} // namespace tectonics::provenance
