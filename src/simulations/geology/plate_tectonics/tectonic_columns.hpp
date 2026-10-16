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
    uint64_t river_unsupplied_cells = 0; // Implicit rivers: cell updates whose column could not supply the erosion.
    double relative_mass_residual = 0;
    // Continental growth (margins::Parameters). The felsic reservoir is the mass of
    // continent and sediment; it grows by arc continent and oceanic crust eroded into
    // sediment, and shrinks by subducted sediment and subduction erosion.
    double generated_continent_kg = 0; // Subset of generated_arc_kg; not in generated_ocean_kg.
    double subducted_continent_kg = 0; // Subduction erosion at overlap cells.
    double eroded_ocean_kg = 0; // Subset of eroded_kg taken from oceanic crust (becomes sediment).
    double oceanized_continental_area_m2 = 0; // Area freed by the minimum-thickness cap.
    double continentalized_oceanic_area_m2 = 0; // Area gained by the maximum-thickness cap.
    double collapsed_continent_kg = 0; // Orogenic collapse; internal transfer, not a mass source.
    double initial_continental_area_m2 = 0; // Sum of continental_area x cell area at initialize.
    double initial_felsic_kg = 0, felsic_kg = 0;
    double felsic_relative_residual = 0;
    double continental_growth_shortfall_kg = 0; // Arc continent the controller asked for beyond the arc melt.
    double arc_continent_fraction = 0; // Share of arc melt made continent in the last update.
};
// Cumulative horizontal strain of the rock: a symmetric 2D tensor on the grid's axes
// (x east, y the row direction), dimensionless, the sum of the per-update strains
// (negative is shortening). Measurement only.
struct Strain { double xx = 0, yy = 0, xy = 0; };
// Shortening, m: minus the most compressive principal component of the strain, times
// cell_m (the length lost across one native cell width, so it scales with the
// resolution; strain = shortening / cell_m); 0 where no axis shortens.
double shortening_m(const Strain& s, double cell_m);
// Axis of the most compressive principal component, radians in [0, pi), from the x
// axis toward the y axis; 0 where the two principal components are equal.
double shortening_direction(const Strain& s);
double mass_per_area(const Cell& c);
// Continent and sediment mass per unit area, kg/m²: the felsic reservoir.
double felsic_mass_per_area(const Cell& c);
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
                 const plumes::Source* plume = nullptr, double uplift_rate_window_myr = 10);
    const std::vector<Cell>& cells() const { return cells_; }
    // Rock uplift rate per cell, m/Myr, signed: an exponential mean (e-folding time
    // uplift_rate_window_myr) of the elevation change of the tectonic sub-steps
    // (transport, subduction overlap, rifting, trench flexure, arcs, plumes, roots,
    // thermal subsidence), before the diffusion, rivers and noise. In the plate's frame
    // (semi-Lagrangian: the elevation before the step and the previous mean are carried
    // with the rock by the engine's own remap, on the subcells with the material's
    // sub-steps and per-cell displacements), so plate motion over relief is neither
    // uplift nor subsidence and the mean follows the rock. Exact (to rounding) for a
    // translation when a single column field varies (land); at plate boundaries, from
    // the second sub-step on, the remap's limiters can act on the carried amounts and
    // the carried values are a bounded approximation. The remap smears the carried
    // mean: at most a(1 - a) subcells^2 of variance per axis and sub-step (a the
    // fractional subcell displacement; a flattened peak), a(1 - a) n / s^2 cells^2 per
    // step (n sub-steps, subdivision s), 9.5 times that in steady state at a 10 Myr
    // window and 1 Myr steps; less for resolved smooth relief. Converging rock mixes;
    // a divergent cell keeps the rock that stayed on its own side; new crust has no
    // history of its own. Zero at start.
    const std::vector<double>& upliftRate() const { return uplift_rate_m_myr_; }
    // Cumulative strain per cell (multiscale plan S5.1), in the plate's frame. Each
    // update first rotates the rock's tensor by the local spin and adds the strain of
    // the update's displacement field u (cells/update): E = sym(grad u), W = (duy/dx -
    // dux/dy) / 2. Each derivative along an axis uses the cell's two neighbours on that
    // axis: when only one is on the same plate (forcing owners; one plate without
    // forcing), the one-sided difference toward it; when both are, the smaller of the two
    // differences in magnitude (their mean when equal); when neither is, 0. So a velocity
    // jump to another plate (subduction, ridge, transform) is never strain; a jump inside
    // a plate (the rotation seam of plate::surfaceVelocityAt) is skipped only where the
    // other difference on that axis is the smaller one. Distributed deformation and
    // rift opening inside a plate count. Exact for a velocity linear on both sides of a
    // cell; a rigid rotation gives E = 0 and W = its angle. The tensor is then carried
    // with the rock as the uplift rate is (same remap, sub-steps and limits; converging
    // rock mixes); the ocean floor created in an opening has zero strain, so a divergent
    // cell keeps the tensor of the rock that stayed on the rest of its area only.
    // Small-strain sum: valid while the strain is small (it reaches about 1.4 at 256
    // columns). Cumulative over the run, no decay; zero at start.
    const std::vector<Strain>& strain() const { return strain_; }
    const Budget& budget() const { return budget_; }
private:
    void spreadRoots(std::vector<Cell>& cells, double dt_myr, double thermal_strength, bool sediment_feedback);
    void collapseOrogens(std::vector<Cell>& cells, double dt_myr, const margins::Parameters& margins);
    // carried: non-negative fields on the subcells, remapped in place with the material.
    std::vector<Cell> transportColumns(const std::vector<Displacement>& velocity, double dt_myr,
                                       std::vector<std::vector<double>>& carried);
    void reconcileTransport(const std::vector<Cell>& evolved);
    void routeSediment(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                       const margins::Parameters& parameters, const margins::Forcing* forcing);
    void routeRiversImplicit(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
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
    std::vector<double> uplift_rate_m_myr_; // Measurement only; never feeds back into the columns.
    std::vector<Strain> strain_; // Measurement only; never feeds back into the columns.
};
} // namespace tectonics::columns
#endif
