#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

// Engine-neutral, read-only view of the native tectonic state after an update.
// Codes use the tectonic contract enums. Arrays are row-major, width x height,
// on the native planar torus; the frame does not own them.
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

struct CellShift
{
    double x = 0.0; // Cells per update.
    double y = 0.0;
};

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
    std::uint32_t width = 0;
    std::uint32_t height = 0;
    double timeMyr = 0.0;
    double deltaTimeMyr = 1.0;
    // Physical width of the native planar domain in kilometres.
    double worldWidthKm = 40000.0;

    const std::uint32_t* plateIndex = nullptr; // >= plates.size(): no owner.
    const std::uint8_t* crustClass = nullptr;
    const float* crustAgeMyr = nullptr;
    const float* elevationM = nullptr;          // Relative to sea level.
    // Continental crust of the isostatic columns: volume per cell area, metres.
    const float* continentalCrustM = nullptr;
    const std::uint16_t* boundaryDistance = nullptr; // Native cells.
    const std::uint32_t* nearestBoundaryId = nullptr;
    const std::uint8_t* geologicRegime = nullptr;
    // Per-cell material displacement applied by the column solver.
    const CellShift* materialShift = nullptr;

    std::vector<PlateIdentity> plates;         // Indexed by plate index.
    std::vector<BoundarySummary> boundaries; // Indexed by boundary segment id.

    std::size_t cellCount() const { return static_cast<std::size_t>(width) * height; }
    std::size_t cellIndex(std::uint32_t x, std::uint32_t y) const
    {
        return static_cast<std::size_t>(y) * width + x;
    }
    // Native cell spacing along the torus axes.
    double cellKm() const { return width > 0 ? worldWidthKm / width : 0.0; }
    double cellAreaKm2() const { return cellKm() * cellKm(); }
};
}
