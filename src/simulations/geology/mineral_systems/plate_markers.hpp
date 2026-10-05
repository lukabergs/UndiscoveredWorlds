#pragma once

#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

// Approximate point tracers of the native column displacement on the torus.
// A marker occupies cell (floor(x), floor(y)). They do not reproduce the
// solver's conservative subcell mixing or track individual parcels of rock.
namespace mineralsystems
{
enum class MarkerState : std::uint8_t
{
    attached = 0,
    subducted = 1,
};

struct Marker
{
    double x = 0.0;
    double y = 0.0;
    std::uint64_t hostUid = 0;
    // Continental and transitional crust is not subducted by the engine; an
    // oceanic host is identified by the time its crust formed.
    bool continentalHost = true;
    double hostCrustBirthMyr = 0.0;
    MarkerState state = MarkerState::attached;
    double stateTimeMyr = 0.0;
};

struct MarkerStepSummary
{
    std::size_t moved = 0;
    std::size_t rehomed = 0;   // Ownership passed to another plate with the crust.
    std::size_t snapped = 0;   // Moved back onto host crust across a raster seam.
    std::size_t orphaned = 0;  // Host plate removed; no motion this update.
    std::size_t subducted = 0;
};

class MarkerTracker
{
public:
    // Birth times closer than the tolerance identify the same oceanic crust.
    explicit MarkerTracker(double birthToleranceMyr = 3.0) : birthToleranceMyr_(birthToleranceMyr) {}

    void begin(const TectonicFrame& frame);
    // Attaches the marker to the plate that owns its cell. Returns false when
    // no plate owns the cell; the marker is then left unchanged.
    bool attach(Marker& marker, const TectonicFrame& frame) const;
    MarkerStepSummary advance(std::vector<Marker>& markers, const TectonicFrame& frame);

private:
    double birthToleranceMyr_;
    std::vector<std::uint64_t> knownUids_; // Sorted plate identities before the update.
};

std::size_t markercell(const Marker& marker, const TectonicFrame& frame);
// Moves the marker through the frame's piecewise-constant material displacement.
void applymaterialshift(Marker& marker, const TectonicFrame& frame);
}
