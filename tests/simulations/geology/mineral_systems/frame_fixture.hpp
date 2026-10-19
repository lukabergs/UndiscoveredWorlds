#pragma once

#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

// Owns the arrays behind a synthetic TectonicFrame on a cell-centred width x width / 2
// latitude-longitude grid of Earth's radius (the frame takes any mesh::CellGrid; the engine's
// is the cube). Cells are addressed by column x and row y (row 0 north); positions by
// fractional (x, y) in the same units; motion and distances in equatorial cells.
struct FrameFixture
{
    static constexpr double kRadiusKm = 6371.0;
    static constexpr double kPi = 3.14159265358979323846;

    FrameFixture(std::uint32_t width, std::uint32_t height)
        : grid(mesh::makelatlongrid(static_cast<int>(width), kRadiusKm * 1000.0)),
          plate(grid->cellcount(), 0),
          crust(plate.size(), mineralsystems::crustcode::continental),
          age(plate.size(), 0.0f),
          elevation(plate.size(), 0.0f),
          continental(plate.size(), 35000.0f),
          distance(plate.size(), 0),
          nearest(plate.size(), 0),
          regime(plate.size(), mineralsystems::regimecode::stable),
          width_(width), height_(height)
    {
        frame.grid = grid.get();
        frame.spacingKm = cellkm();
        frame.boundaries.resize(1);
    }

    std::uint32_t width() const { return width_; }
    std::uint32_t height() const { return height_; }
    double cellkm() const { return 2.0 * kPi * kRadiusKm / width_; }  // equatorial
    std::size_t cell(std::uint32_t x, std::uint32_t y) const
    {
        return grid->cellof(static_cast<int>(x), static_cast<int>(y));
    }
    // The unit vector at fractional cell coordinates (x east from -180 degrees, y south from
    // the north pole), and back (x in [0, width)).
    mesh::Vec3 at(double x, double y) const
    {
        return mesh::unitvector(0.5 * kPi - kPi * y / height_, -kPi + 2.0 * kPi * x / width_);
    }
    std::array<double, 2> coordinates(const mesh::Vec3& p) const
    {
        double x = (mesh::longitudeof(p) + kPi) / (2.0 * kPi) * width_;
        if (x >= width_) x -= width_;
        return { x, (0.5 * kPi - mesh::latitudeof(p)) / kPi * height_ };
    }

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

    // A uniform displacement per update, (east, north) in equatorial cells at every cell.
    void motion(double east, double north = 0)
    {
        const double cell_m = cellkm() * 1000.0;
        shift.resize(plate.size());
        for (mesh::CellId c = 0; c < shift.size(); ++c)
            shift[c] = (east * cell_m) * grid->east(c) + (north * cell_m) * grid->north(c);
    }

    // A rigid rotation of every plate per update (a rotation vector, radians): each cell's
    // displacement is omega x R c and the frame carries the plates' rotations.
    void rotation(const mesh::Vec3& omega)
    {
        const double radius_m = kRadiusKm * 1000.0;
        shift.resize(plate.size());
        for (mesh::CellId c = 0; c < shift.size(); ++c)
            shift[c] = mesh::cross(omega, radius_m * grid->centre(c));
        frame.plateRotation.assign(frame.plates.size(), omega);
    }

    std::shared_ptr<const mesh::CellGrid> grid;
    std::vector<mesh::Vec3> shift;
    mineralsystems::TectonicFrame frame;
    std::vector<std::uint32_t> plate;
    std::vector<std::uint8_t> crust;
    std::vector<float> age;
    std::vector<float> elevation;
    std::vector<float> continental; // Column continental crust, metres.
    std::vector<std::uint16_t> distance;
    std::vector<std::uint32_t> nearest;
    std::vector<std::uint8_t> regime;

private:
    std::uint32_t width_, height_;
};
