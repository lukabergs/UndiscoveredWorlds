#pragma once

#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <vector>

// Reads the native tectonic engine state into a TectonicFrame. The frame views
// engine arrays, so read() must be repeated after every update.
namespace mineralsystems
{
class TectonicFrameReader
{
public:
    explicit TectonicFrameReader(double worldWidthKm) { frame_.worldWidthKm = worldWidthKm; }
    const TectonicFrame& read(void* simulation);

private:
    TectonicFrame frame_;
    std::vector<float> elevationM_;
    std::vector<float> continentalCrustM_;
    std::vector<CellShift> materialShift_;
};
}
