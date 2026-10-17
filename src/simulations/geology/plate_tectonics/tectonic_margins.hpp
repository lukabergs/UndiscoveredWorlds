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
    // Rivers by the implicit landscape kernel on the torus; waves stay explicit.
    // Implicit rivers always end at sea level at an outlet (river_outlet_sea_level
    // has no effect), ignore the forcing's sediment_multiplier, require
    // alluvial_transport 0, and pass a depression's load on to its spill (G = 0:
    // closed basins do not trap sediment, unlike the explicit path's pits).
    bool river_implicit = false;
    double sediment_travel_km = 300;
    double alluvial_transport = 0; // 0: sediment bypass; >0: capacity / stream-power erosion potential.
    bool river_outlet_sea_level = false; // Clip river outlet slope at sea level; submarine routing still follows seabed.
    double rift_breakup_thickness_m = 12000; // Onset of gradual oceanization during extension.
    // Continental crust's least mean thickness, m (0: off, the legacy opening share).
    // Above 0 the opening's oceanic share rises from the onset to it, and a column's
    // continental area relaxes (5 Myr) to at most continent / this value: the freed
    // area becomes oceanic area, no rock moves.
    double continental_min_thickness_m = 0;
    // Continental crust's greatest mean thickness, m (0: off). Above 0 a column with
    // oceanic area whose continent exceeds continental area x this value gains
    // continental area (10 Myr) toward continent / this value from its oceanic area;
    // no rock moves. Requires margins and a value above the rift onset and the minimum.
    double continental_max_thickness_m = 0;
    // Orogenic collapse, km²/Myr (0: off): linear diffusion of the cell-mean continent
    // above continental_max_thickness_m into lower neighbours, interior and pure-ocean
    // cells included; continental age moves with the rock, sediment stays. Requires
    // that maximum.
    double orogenic_collapse_km2_myr = 0;
    // Onset of orogenic collapse, m of cell-mean continent (0, legacy: the collapse acts
    // above continental_max_thickness_m). Above 0 it replaces that threshold in the
    // collapse only, so margins keep the area cap while interiors thicken to the onset:
    // orogens spread gravitationally only above about 60-70 km of crust (Tibet ~70 km);
    // below it the crust is held by its strength. 0 or at least the maximum thickness.
    double orogenic_collapse_onset_thickness_m = 0;
    // With an onset: a face whose receiving cell has oceanic area (continental_area < 1:
    // a margin or ocean cell) collapses above continental_max_thickness_m as in legacy;
    // a face between two cells without oceanic area waits for the onset. Thick crust
    // spreads into its weak margins and forelands (the upper cap turns them into
    // continent); within a continent interior it is held until the onset. Requires an onset.
    bool orogenic_collapse_margin_cap = false;
    // With an onset, km²/Myr (0: off): weak creep between continental_max_thickness_m and
    // the onset, D_low (min(c_a, onset) - max(c_b, h_max))+, added to the full collapse
    // D (c_a - max(c_b, onset))+. With orogenic_collapse_margin_cap the creep applies to
    // interior faces only (a margin receiver takes the full collapse above h_max).
    // Requires an onset.
    double orogenic_collapse_below_onset_km2_myr = 0;
    // Share of arc melt emplaced as continental crust at age 0 (mass kept, volume
    // x 2900/2800); an outlet with mostly oceanic area also gains continental area
    // (added thickness / continental_min_thickness_m). Above 0 it requires that minimum.
    double arc_continent_fraction = 0;
    // Continental volume removed at an overlap cell per subducted oceanic volume,
    // capped by the continent there: subduction erosion, a sink of the reservoir.
    double subduction_erosion_ratio = 0;
    // A global feedback on the arc share in place of arc_continent_fraction: it holds
    // the felsic reservoir (continent and sediment mass) at its initial mass times
    // (1 + continental_growth_per_gyr x time in Gyr). Requires continental_min_thickness_m.
    bool continental_growth_control = false;
    double continental_growth_per_gyr = 0;
    // The controller's reservoir is the continental crust mass only (sediment out): its
    // target starts from the initial continent mass and only subduction erosion is booked
    // as a sink. False: the felsic reservoir.
    bool continental_growth_continent_only = false;
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
