#ifndef TECTONIC_NOISE_PASS_HPP
#define TECTONIC_NOISE_PASS_HPP
#include <cstdint>
#include <vector>
namespace tectonics::columns { struct Cell; class Metric; }
namespace tectonics::noise_pass {
// Deliberately synthetic detail experiment, not an erosion or tectonic law.
struct Parameters {
    double amplitude_m = 0; // Height-equivalent potential, not a guaranteed RMS.
    double wavelength_km = 600;
    double period_myr = 180;
    double coast_taper_m = 500; // Smooth suppression near sea level; zero reproduces the untapered experiment.
};
void validate(const Parameters& parameters);
// Conservative face transfers driven by the change in a bounded, continuous
// noise potential. No new material, per-step RNG draws, or initialization pass.
// The transfers run over the cell table (metric; null: the planar scaffold of width x
// height); the noise itself is planar until E4 (a cube metric with amplitude > 0 throws).
double apply(uint32_t width,uint32_t height,double cell_m,std::vector<columns::Cell>& cells,
             long seed,double start_myr,double end_myr,const Parameters& parameters,double thermal_strength = 1,
             double root_kg_m3 = 0, double root_max_thickness_m = 0, double root_delamination_width_m = 0,
             const columns::Metric* metric = nullptr);
}
#endif
