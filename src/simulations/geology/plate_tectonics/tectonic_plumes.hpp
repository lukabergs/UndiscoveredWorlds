#ifndef TECTONIC_PLUMES_HPP
#define TECTONIC_PLUMES_HPP
#include "core/mesh/cell_grid.hpp"
#include <cstdint>
#include <vector>
namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::plumes {
// Mantle plumes fixed in the native (mantle) frame. A prescribed source, not a
// convection solution: each plume adds basalt at a steady tail flux, optionally
// a smooth head pulse at its birth, and reheats the oceanic lithosphere above it.
// Plates carrying that rock away build age-progressive chains. count 0 is off.
struct Parameters {
    uint32_t count = 0;
    double tail_flux_km3_yr = .05;     // basalt per plume and year after birth
    double tail_radius_km = 150;       // Gaussian sigma of the tail source
    double head_volume_km3 = 0;        // total head pulse per plume; 0 = no heads
    double head_radius_km = 800;
    double head_duration_myr = 10;     // smoothstep pulse from the birth time
    double birth_window_myr = 0;       // births spread over [0, window] of model time
    double swell_age_myr = 0;          // target thermal age of reheated ocean floor; 0 = no swell
    double swell_radius_km = 400;
    double weakening = 0;              // 0..1 strength reduction for progressive rifting
    double min_spacing_km = 3000;
};
// position: a unit vector (placesphere).
struct Plume { double birth_myr = 0; mesh::Vec3 position{}; };
struct Source {
    std::vector<double> basalt_m;      // added basalt, volume per cell area over the step
    std::vector<double> influence;     // 0..1 proximity to an active plume (swell, weakening, map)
    double swell_age_myr = 0;
};
void validate(const Parameters& parameters);
// Deterministic positions on the sphere of radius_km from the seed alone: stateless hashes of
// the seed read as a point uniform on the sphere (mesh::uniformpoint), spacing by great circle.
std::vector<Plume> placesphere(double radius_km, long seed, const Parameters& parameters);
// Source over [start_myr, end_myr]; volumes integrate exactly over any time step. The plumes
// are placed on the sphere, the Gaussians use great-circle distances and the basalt is the
// volume share per true cell area.
Source source(const grid::TectonicGrid& grid, long seed, double start_myr, double end_myr,
              const Parameters& parameters);
}
#endif
