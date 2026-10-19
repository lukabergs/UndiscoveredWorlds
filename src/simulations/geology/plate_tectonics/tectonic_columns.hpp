#ifndef TECTONIC_COLUMNS_HPP
#define TECTONIC_COLUMNS_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include "simulations/geology/plate_tectonics/tectonic_plumes.hpp"

namespace tectonics::columns {
using Displacement = transport::Displacement;

// The engine's cell table as the column processes read it (W2-E, E2B): true cell areas,
// finite-volume face weights (shared length / centre distance), the eight-neighbour stencils
// and the displacement adapter. On the planar scaffold every weight and area ratio is
// exactly 1, every area cell_m², every distance cell_m, so the processes reproduce the
// torus bit for bit; E3 deletes the planar branch.
class Metric {
public:
    struct Neighbour { grid::CellId cell = grid::nocell; double distance_m = 0; };
    Metric() = default;
    // A null grid builds the planar scaffold of width x height cells of cell_m (unit tests);
    // a planar grid must have that size; on the cube the three are ignored.
    Metric(const grid::TectonicGrid* grid, uint32_t width, uint32_t height, double cell_m);

    const grid::TectonicGrid& grid() const { return *grid_; }
    bool planar() const { return planar_; }
    size_t cellcount() const { return grid_ ? grid_->cellcount() : 0; }
    double area_m2(grid::CellId c) const { return planar_ ? cell_area_m2_ : area_[c]; }
    // The finite-volume reference area: planar cell_m²; on the cube the least area / weight
    // over every cell and face, so an explicit coefficient per reference area bounds the
    // coefficient of every face in its source cell's thickness.
    double reference_area_m2() const { return reference_area_m2_; }
    double scale(grid::CellId c) const { return planar_ ? 1.0 : scale_[c]; } // reference area / area
    // The edge between c and its neighbour in direction k (an index of grid().edges()).
    size_t edge(grid::CellId c, uint8_t k) const;
    double weight(size_t e) const { return planar_ ? 1.0 : weight_[e]; }   // length / distance
    double distance_m(size_t e) const { return planar_ ? cell_m_ : grid_->edges()[e].distance_m; }
    double width_m(grid::CellId c) const { return planar_ ? cell_m_ : width_[c]; } // planar cell_m, cube sqrt(area)
    // A per-area quantity times the cell's area (planar: the torus' x cell_m cell_m).
    double mass(double per_area, grid::CellId c) const {
        return planar_ ? per_area * cell_m_ * cell_m_ : per_area * area_[c];
    }
    // A per-reference-area quantity times the reference area.
    double referencemass(double per_area) const {
        return planar_ ? per_area * cell_m_ * cell_m_ : per_area * reference_area_m2_;
    }
    // The eight neighbours in the planar stencils' (dy, dx) order: NW, N, NE, W, E, SW, S, SE
    // (mesh slots 8, 1, 2, 7, 3, 6, 5, 4); nocell for the diagonal at a cube corner.
    std::array<Neighbour, 8> eight(grid::CellId c) const;

    // The displacement adapter (identity in E3): a cell's displacement per update as a 3-D
    // tangent vector in displacement units. Planar: (x, y, 0) in cells along +x and +y
    // (unit cell_m). Cube: metres (unit 1) on the cell's frame, x along the unit tangent
    // toward the east neighbour (the face's u), y orthogonal to it toward the south one.
    double displacementunit_m() const { return planar_ ? cell_m_ : 1.0; }
    // An edge's centre distance in displacement units (planar 1).
    double displacementdistance(size_t e) const { return planar_ ? 1.0 : grid_->edges()[e].distance_m; }
    mesh::Vec3 vector(grid::CellId c, const Displacement& d) const;
    Displacement local(grid::CellId c, const mesh::Vec3& v) const;
    // A displacement of `from` on the frame of `to` (planar: unchanged).
    Displacement inframe(grid::CellId to, grid::CellId from, const Displacement& d) const;
    // The closing motion across edge e, displacement units: (d_a - d_b) along the normal from
    // a to b (planar: the axis component's difference).
    double closing(size_t e, const Displacement& a, const Displacement& b) const;

    // Position of a cell for position-placed variation, in [0, 1)²: planar the cell centre's
    // fraction of the torus; cube a lat-lon proxy ((longitude + pi) / 2 pi, (pi / 2 -
    // latitude) / pi) until E4. extent_km: the proxy's size (planar width and height in km,
    // cube 2 pi R and pi R).
    std::array<double, 2> unitposition(grid::CellId c) const;
    std::array<double, 2> extent_km() const;

private:
    std::shared_ptr<const grid::TectonicGrid> owned_;
    const grid::TectonicGrid* grid_ = nullptr;
    bool planar_ = true;
    uint32_t width_cells_ = 0, height_cells_ = 0;
    double cell_m_ = 0, cell_area_m2_ = 0, reference_area_m2_ = 0;
    std::vector<double> area_, scale_, width_, weight_;
    std::vector<uint32_t> edges_;          // 4 per cell (cube)
    std::vector<mesh::Vec3> east_, south_; // the displacement frame (cube)
};

inline constexpr double mantle_density = 3300;
inline constexpr double continental_density = 2800;
inline constexpr double oceanic_density = 2900;
inline constexpr double sediment_density = 2700;
inline constexpr double water_density = 1030;
inline constexpr double reference_continent_m = 35000;
// The root's delamination cap (scenario::Scenario::continental_root_max_thickness_m):
// 0 (no cap) or a crust thickness above reference_continent_m.
inline bool valid_root_max_thickness(double m) { return m == 0 || (std::isfinite(m) && m > reference_continent_m); }
// The width of the root's delamination ramp
// (scenario::Scenario::continental_root_delamination_width_m): 0 (the frozen cap), or
// finite and positive with a cap. The ramp never makes the elevation decrease with the
// crust (elevation_m), so no width needs a further bound.
inline bool valid_root_delamination_width(double max_thickness_m, double width_m) {
    return width_m == 0 || (std::isfinite(width_m) && width_m > 0 && max_thickness_m > 0);
}
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
    // Continental area terms, measurement only: area lost where coverage stays above 1
    // after the oceanic overlap sinks (collision), the continental share of an opening
    // (stretched crust) and the oceanic area converted at island-arc outlets. With the
    // caps: area change = stretched + arc + continentalized - collided - oceanized.
    double collided_continental_area_m2 = 0;
    double stretched_continental_area_m2 = 0;
    double arc_continental_area_m2 = 0;
    double collapsed_continent_kg = 0; // Orogenic collapse; internal transfer, not a mass source.
    double initial_continental_area_m2 = 0; // Sum of continental_area x cell area at initialize.
    double initial_felsic_kg = 0, felsic_kg = 0;
    // Continental crust mass: the reservoir of margins::Parameters::continental_growth_continent_only.
    double initial_continent_kg = 0, continent_kg = 0;
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
// root_kg_m3: the lithospheric root load per metre of crust away from the reference
// continent (scenario::Scenario::continental_root_kg_m3); 0 is the crust-only column.
// root_max_thickness_m: the crust thickness up to which the root load is full (0 no
// cap). root_delamination_width_m: 0, the root stops growing above the cap (frozen);
// positive, the root delaminates: its load (the whole root, root x (h - 35 km) per unit
// continental area) falls linearly to zero between the cap and cap + width, h the
// per-area crust thickness continent / continental_area (cells without continental
// area carry no root load, as with the frozen cap). Over the ramp d(buoyancy)/dh =
// (500 - root x (1 - (2h - cap - 35 km) / width)), at least 500 - root > 0: the
// elevation rises with the crust for every width. seed() inverts elevation_m() with
// both. mountain_force_n_m keeps the uncapped (500 - root) / 500 scale: neither the cap
// nor the ramp enters the GPE estimate.
double elevation_m(const Cell& c, double thermal_strength = 1, double root_kg_m3 = 0,
                   double root_max_thickness_m = 0, double root_delamination_width_m = 0);
double mountain_force_n_m(const Cell& c, bool sediment_feedback = false, double root_kg_m3 = 0);
double age_myr(const Cell& c);
Cell seed(double elevation_m, bool continental, double age_myr, double root_kg_m3 = 0,
          double root_max_thickness_m = 0, double root_delamination_width_m = 0);

// The column law a model evaluates: elevation_m with its arguments.
struct Isostasy {
    double thermal_strength = 1, root_kg_m3 = 0, root_max_thickness_m = 0, root_delamination_width_m = 0;
    double elevation(const Cell& c) const {
        return elevation_m(c, thermal_strength, root_kg_m3, root_max_thickness_m, root_delamination_width_m);
    }
};
// The column processes on the cell table, in the order Model::advance runs them (exposed
// so the cube is tested before E3 runs the whole update there). Each conserves the mass it
// moves in true areas and books it in the budget.
// Lower-crust root flow (n = 3 creep), sub-stepped below the explicit limit of the smallest
// reference face; dt_myr already scaled by margins::Parameters::root_mobility.
void spread_roots(const Metric& metric, std::vector<Cell>& cells, double dt_myr, const Isostasy& isostasy,
                  bool sediment_feedback, Budget& budget);
// Orogenic collapse (margins::Parameters::orogenic_collapse_*), sub-stepped below 0.1 per face.
void collapse_orogens(const Metric& metric, std::vector<Cell>& cells, double dt_myr,
                      const margins::Parameters& parameters, Budget& budget);
// The surface diffusion's explicit number per reference face, transport_diffusivity x
// erosion_strength x dt / reference area (Model::advance rejects above 0.2).
double surface_diffusion_number(const Metric& metric, double erosion_strength, double dt_myr);
// Conservative downhill transport of the elevation field (one value per cell) at that number.
void diffuse_surface(const Metric& metric, std::vector<Cell>& cells, const std::vector<double>& elevation,
                     double diffusion_number, double submarine_erosion_ratio, Budget& budget);

class Model {
public:
    void initialize(uint32_t width, uint32_t height, double cell_m, std::vector<Cell> cells,
                    uint32_t transport_subdivision = 1, long seed = 0,
                    const ocean_relief::Parameters& relief = {}, bool preserve_initial = false,
                    double root_kg_m3 = 0, double root_max_thickness_m = 0,
                    double root_delamination_width_m = 0, const grid::TectonicGrid* grid = nullptr);
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
    // a plate (a rotation seam, where the offset from the plate centre wraps) is skipped only where the
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
    // carried: non-negative fields on the subcells, remapped in place with the material.
    std::vector<Cell> transportColumns(const std::vector<Displacement>& velocity, double dt_myr,
                                       std::vector<std::vector<double>>& carried);
    void reconcileTransport(const std::vector<Cell>& evolved);
    void routeSediment(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                       const margins::Parameters& parameters, const margins::Forcing* forcing);
    void routeRiversImplicit(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                             const margins::Parameters& parameters, const margins::Forcing* forcing);
    uint32_t width_ = 0, height_ = 0;
    Metric metric_; // the engine's cell table (W2-E); a planar scaffold of its own in unit tests
    double cell_m_ = 0;
    long seed_ = 0;
    double root_kg_m3_ = 0; // Lithospheric root load in every elevation of the model.
    double root_max_thickness_m_ = 0; // The root's delamination cap (0 none), with root_kg_m3_.
    double root_delamination_width_m_ = 0; // The delamination ramp above the cap (0 frozen root).
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
