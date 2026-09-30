#pragma once

#include "tectonic_frame.hpp"

#include <cstdint>
#include <vector>

// Owns the arrays behind a synthetic TectonicFrame.
struct FrameFixture
{
    FrameFixture(std::uint32_t width, std::uint32_t height)
        : plate(static_cast<std::size_t>(width) * height, 0),
          crust(plate.size(), mineralsystems::crustcode::continental),
          age(plate.size(), 0.0f),
          elevation(plate.size(), 0.0f),
          continental(plate.size(), 35000.0f),
          distance(plate.size(), 0),
          nearest(plate.size(), 0),
          regime(plate.size(), mineralsystems::regimecode::stable)
    {
        frame.width = width;
        frame.height = height;
        frame.boundaries.resize(1);
    }

    std::size_t cell(std::uint32_t x, std::uint32_t y) const { return static_cast<std::size_t>(y) * frame.width + x; }

    // Refreshes the views; call after resizing any array.
    const mineralsystems::TectonicFrame& view()
    {
        frame.plateIndex = plate.data();
        frame.crustClass = crust.data();
        frame.crustAgeMyr = age.data();
        frame.elevationM = elevation.data();
        frame.continentalCrustM = continental.data();
        frame.boundaryDistance = distance.data();
        frame.nearestBoundaryId = nearest.data();
        frame.geologicRegime = regime.data();
        frame.materialShift = shift.empty() ? nullptr : shift.data();
        return frame;
    }

    std::vector<mineralsystems::CellShift> shift;
    void motion(double x, double y = 0) { shift.assign(plate.size(), {x, y}); }

    mineralsystems::TectonicFrame frame;
    std::vector<std::uint32_t> plate;
    std::vector<std::uint8_t> crust;
    std::vector<float> age;
    std::vector<float> elevation;
    std::vector<float> continental; // Column continental crust, metres.
    std::vector<std::uint16_t> distance;
    std::vector<std::uint32_t> nearest;
    std::vector<std::uint8_t> regime;
};
