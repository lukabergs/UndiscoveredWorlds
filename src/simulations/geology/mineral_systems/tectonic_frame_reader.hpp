#pragma once

#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <memory>
#include <vector>

// Reads the native tectonic engine state into a TectonicFrame on the engine's cube (true
// areas, displacements as tangent vectors in metres). The frame views engine arrays, so
// read() must be repeated after every update.
namespace mineralsystems
{
class TectonicFrameReader
{
public:
    const TectonicFrame& read(void* simulation);

private:
    TectonicFrame frame_;
    std::shared_ptr<const mesh::CellGrid> grid_;
    std::vector<float> elevationM_;
    std::vector<float> continentalCrustM_;
    std::vector<mesh::Vec3> materialShift_;
};
}
