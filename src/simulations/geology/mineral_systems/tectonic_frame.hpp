#pragma once

#include "core/mesh/cell_grid.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

// Engine-neutral, read-only view of the native tectonic state after an update.
// Codes use the tectonic contract enums. Arrays are per cell of the engine's grid (the cube,
// mesh::CellGrid, in its id order); cells have their true areas. The frame does not own them.
namespace mineralsystems
{
// Codes from tectonics::contract; tectonic_frame_reader.cpp checks that they match.
namespace crustcode
{
inline constexpr std::uint8_t none = 0, oceanic = 1, transitional = 2, continental = 3;
}
namespace boundarycode
{
inline constexpr std::uint8_t none = 0, convergent = 1, divergent = 2, transform = 3, passive_margin = 4;
}
namespace regimecode
{
inline constexpr std::uint8_t stable = 0, convergent_arc = 1, continent_collision = 2, divergent_rift = 3,
    transform = 4, passive_margin = 5, mid_ocean_ridge = 6, trench_adjacent = 7;
}

inline constexpr std::uint32_t kNoPlate = 0xffffffffu;

struct PlateIdentity
{
    std::uint64_t uid = 0; // Stable topology identity; array indices can change.
};

struct BoundarySummary
{
    std::uint8_t boundaryType = boundarycode::none;
    std::uint8_t geologicRegime = regimecode::stable;
    std::uint8_t convergenceScore = 0; // Segment means, 0..100.
    std::uint8_t shearScore = 0;
    float ageMyr = 0.0f;
    std::uint32_t leftPlate = kNoPlate; // Plate indices on either side.
    std::uint32_t rightPlate = kNoPlate;
};

struct TectonicFrame
{
    // The engine's cells; every array below has grid->cellcount() values in its id order.
    const mesh::CellGrid* grid = nullptr;
    double timeMyr = 0.0;
    double deltaTimeMyr = 1.0;
    // The engine's cell spacing in kilometres (the mean cell side): the unit of
    // boundaryDistance.
    double spacingKm = 0.0;

    const std::uint32_t* plateIndex = nullptr; // >= plates.size(): no owner.
    const std::uint8_t* crustClass = nullptr;
    const float* crustAgeMyr = nullptr;
    const float* elevationM = nullptr;          // Relative to sea level.
    // Continental crust of the isostatic columns: volume per cell area, metres.
    const float* continentalCrustM = nullptr;
    const std::uint16_t* boundaryDistance = nullptr; // Engine cell spacings.
    const std::uint32_t* nearestBoundaryId = nullptr;
    const std::uint8_t* geologicRegime = nullptr;
    // Per-cell material displacement applied by the column solver: a 3-D tangent vector at
    // the cell's centre, metres per update.
    const mesh::Vec3* materialShift = nullptr;

    std::vector<PlateIdentity> plates;         // Indexed by plate index.
    std::vector<BoundarySummary> boundaries; // Indexed by boundary segment id.

    std::size_t cellCount() const { return grid != nullptr ? grid->cellcount() : 0; }
    double cellAreaKm2(std::size_t cell) const { return grid->area_m2(static_cast<mesh::CellId>(cell)) * 1.0e-6; }
};
}
