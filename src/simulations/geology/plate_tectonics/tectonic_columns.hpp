#ifndef TECTONIC_COLUMNS_HPP
#define TECTONIC_COLUMNS_HPP

#include <cstdint>
#include <vector>
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include "simulations/geology/plate_tectonics/tectonic_plumes.hpp"

namespace tectonics::columns {
inline constexpr double mantle_density = 3300;
inline constexpr double continental_density = 2800;
inline constexpr double oceanic_density = 2900;
inline constexpr double sediment_density = 2700;
inline constexpr double water_density = 1030;
inline constexpr double reference_continent_m = 35000;
inline constexpr double ridge_crust_m = 7000;
inline constexpr double ridge_depth_m = 2500;
// Coarse regional sediment transport, m²/Myr. Not a resolved river model.
inline constexpr double transport_diffusivity = 2e8;

struct Cell {
    double continent = 0, ocean = 0, sediment = 0; // volume / cell area, metres
    double continental_area = 0, oceanic_area = 0; // normalized surface fractions
    double continental_age_volume = 0, oceanic_age_volume = 0; // m * Myr
    double thermal_rejuvenation_volume = 0; // Ocean volume * reduction in thermal age (Myr).
    double flexure_m = 0; // Transient downward response to subduction load, not lost rock.
};
using Displacement = transport::Displacement;
struct Budget {
    bool enabled = false;
    double initial_mass_kg = 0, mass_kg = 0;
    double generated_ocean_kg = 0, subducted_ocean_kg = 0;
    double subducted_sediment_kg = 0;
    double generated_arc_kg = 0; // Subset of generated_ocean_kg; mantle-derived basalt.
    double generated_plume_kg = 0; // Subset of generated_ocean_kg; mantle plume basalt.
    double river_eroded_kg = 0; // Subset of eroded_kg, routed downstream to sediment sinks.
    double eroded_kg = 0, deposited_kg = 0;
    double redistributed_continent_kg = 0;
    double redistributed_sediment_kg = 0;
    double noise_redistributed_kg = 0; // Synthetic experiment; internal transfer, not a mass source.
    double relative_mass_residual = 0;
};
double mass_per_area(const Cell& c);
double elevation_m(const Cell& c, double thermal_strength = 1);
double mountain_force_n_m(const Cell& c, bool sediment_feedback = false);
double age_myr(const Cell& c);
Cell seed(double elevation_m, bool continental, double age_myr);

class Model {
public:
    void initialize(uint32_t width, uint32_t height, double cell_m, std::vector<Cell> cells,
                    uint32_t transport_subdivision = 1, long seed = 0,
                    const ocean_relief::Parameters& relief = {}, bool preserve_initial = false);
    void advance(const std::vector<Displacement>& velocity, double dt_myr,
                 double erosion_strength, double thermal_strength, double submarine_erosion_ratio = 1,
                 const margins::Forcing* forcing = nullptr, const margins::Parameters& margins = {},
                 const noise_pass::Parameters& noise = {}, bool sediment_feedback = false,
                 const plumes::Source* plume = nullptr);
    const std::vector<Cell>& cells() const { return cells_; }
    const Budget& budget() const { return budget_; }
private:
    void spreadRoots(std::vector<Cell>& cells, double dt_myr, double thermal_strength, bool sediment_feedback);
    std::vector<Cell> transportColumns(const std::vector<Displacement>& velocity, double dt_myr);
    void reconcileTransport(const std::vector<Cell>& evolved);
    void routeSediment(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                       const margins::Parameters& parameters, const margins::Forcing* forcing);
    uint32_t width_ = 0, height_ = 0;
    double cell_m_ = 0;
    long seed_ = 0;
    double time_myr_ = 0;
    std::vector<Cell> cells_;
    uint32_t subdivision_ = 1;
    std::vector<Cell> transport_cells_; // Persistent rock subcells; never re-seeded between updates.
    std::vector<double> newborn_ocean_crust_; // Fixed mantle-space birth pattern; rock then advects.
    Budget budget_;
};
} // namespace tectonics::columns
#endif
