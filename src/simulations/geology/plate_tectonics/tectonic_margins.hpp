#ifndef TECTONIC_MARGINS_HPP
#define TECTONIC_MARGINS_HPP
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include <cstdint>
#include <vector>
#include <map>

namespace tectonics::columns { struct Cell; }
namespace tectonics::margins {
struct Parameters {
    bool enabled = false;
    double trench_load_n_m = 1.2e13; // Effective downward line load; not a depth cap.
    double flexural_length_km = 120;
    double trench_response_myr = 5;
    double arc_distance_km = 200;
    double arc_melt_ratio = .75; // Mantle-derived basalt mass / removed slab basalt mass.
    double river_erosion = 5; // Stream-power coefficient: m/Myr per km of sqrt(drainage area).
    double sediment_travel_km = 300;
    double alluvial_transport = 0; // 0: sediment bypass; >0: capacity / stream-power erosion potential.
    bool river_outlet_sea_level = false; // Clip river outlet slope at sea level; submarine routing still follows seabed.
    double rift_breakup_thickness_m = 12000; // Onset of gradual oceanization during extension.
    double process_variation = 0;
    double process_wavelength_km = 1800;
    double process_period_myr = 80;
    double thermal_age_rate = 1; // Effective thermal clock / chronological clock; keeps initial state and asymptotic contraction.
    double wave_erosion_m_myr = 0; // Unresolved coastal retreat, converted to regional rock removal.
    double wave_height_scale_m = 200; // Wave access decays inland/up coastal relief; not a height limit.
    double wave_retreat_km_myr = 0; // Subgrid cliff retreat converts coastal relief into removed rock volume.
    double root_mobility = 1; // Multiplier of lower-crust creep, inverse effective viscosity.
    double sediment_subduction_fraction = 0; // Remainder accretes; 0 preserves the old fully retained sediment budget.
};
struct Forcing {
    std::vector<double> thermal_age_myr; // -1: no active spreading heat source.
    std::vector<double> trench_m;
    std::vector<uint32_t> arc_receiver; // UINT32_MAX: no overriding plate destination.
    std::vector<uint32_t> owners;
    std::vector<double> arc_multiplier, erosion_multiplier, sediment_multiplier;
};
void validate(const Parameters& p);
using Polarity = std::map<std::pair<uint32_t,uint32_t>,uint32_t>; // Stable plate pair -> descending plate.
Forcing forcing(uint32_t w, uint32_t h, double cell_m, double dt_myr,
    const std::vector<uint32_t>& owners, const std::vector<columns::Cell>& cells,
    const std::vector<transport::Displacement>& rigid, const Parameters& p, Polarity* polarity = nullptr,
    long seed = 0, double time_myr = 0);
}
#endif
