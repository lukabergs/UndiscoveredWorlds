#ifndef TECTONIC_COLUMNS_HPP
#define TECTONIC_COLUMNS_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_sphere_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include "simulations/geology/plate_tectonics/tectonic_plumes.hpp"

namespace tectonics::columns {
// A displacement per update: a 3-D tangent vector in metres at the cell's centre.
using Displacement = mesh::Vec3;

// The engine's cell table as the column processes read it (W2-E, E2B): true cell areas,
// finite-volume face weights (shared length / centre distance) with the non-orthogonal
// correction, the eight-neighbour stencils and the displacement adapter.
class Metric {
public:
    struct Neighbour { grid::CellId cell = grid::nocell; double distance_m = 0; };
    Metric() = default;
    // The grid must outlive the metric.
    explicit Metric(const grid::TectonicGrid& grid);

    const grid::TectonicGrid& grid() const { return *grid_; }
    size_t cellcount() const { return grid_ ? grid_->cellcount() : 0; }
    double area_m2(grid::CellId c) const { return area_[c]; }
    // The finite-volume reference area: the least area / weight over every cell and face, so
    // an explicit coefficient per reference area bounds the coefficient of every face in its
    // source cell's thickness.
    double reference_area_m2() const { return reference_area_m2_; }
    double scale(grid::CellId c) const { return scale_[c]; } // reference area / area
    // The edge between c and its neighbour in direction k (an index of grid().edges()).
    size_t edge(grid::CellId c, uint8_t k) const { return grid_->edge(c, k); }
    double weight(size_t e) const { return weight_[e]; }   // length / distance
    double distance_m(size_t e) const { return grid_->edges()[e].distance_m; }
    double width_m(grid::CellId c) const { return width_[c]; } // sqrt(area)
    // A per-area quantity times the cell's area.
    double mass(double per_area, grid::CellId c) const { return per_area * area_[c]; }
    // A per-reference-area quantity times the reference area.
    double referencemass(double per_area) const { return per_area * reference_area_m2_; }
    // The eight neighbours in the order NW, N, NE, W, E, SW, S, SE (mesh slots 8, 1, 2, 7, 3, 6,
    // 5, 4); nocell for the diagonal at a cube corner.
    std::array<Neighbour, 8> eight(grid::CellId c) const;

    // The displacement adapter: displacements per update are 3-D tangent vectors in metres
    // (E2D's convention), projected onto the stencil directions: the edge normals (closing)
    // and the cell's two stencil axes (axis 0 the unit tangent toward the east neighbour,
    // axis 1 orthogonal to it toward the south one).
    double component(grid::CellId c, const mesh::Vec3& v, uint8_t axis) const {
        return mesh::dot(v, axis == 0 ? east_[c] : south_[c]);
    }
    mesh::Vec3 axis(grid::CellId c, uint8_t axis) const { return axis == 0 ? east_[c] : south_[c]; }
    // A displacement of another cell on the tangent plane of `to`.
    mesh::Vec3 inframe(grid::CellId to, const mesh::Vec3& v) const;
    // The closing motion across edge e, metres: (d_a - d_b) along the face normal from a to b.
    double closing(size_t e, const mesh::Vec3& a, const mesh::Vec3& b) const;
    // The least-squares gradient weights of cell c (buildcorrection): grad phi = sum over the
    // four cardinal neighbours k of g_k (phi_k - phi_c), g_k as (axis 0, axis 1) components
    // in metres^-1.
    const std::array<std::array<double, 2>, 4>& gradientweights(grid::CellId c) const { return gradient_[c]; }

    // The difference phi_a - phi_b across edge e as the two-point flux reads it, with the
    // non-orthogonal correction (minimum correction from least-squares cell gradients), so the
    // flux (length / distance) x difference is consistent where the grid lines are skew:
    // (n . t)(phi_a - phi_b) - distance k . (grad_a + grad_b) / 2, n the face normal, t the
    // centre line, k = n - (n . t) t.
    double difference(size_t e, const std::vector<double>& phi) const;
    // The largest (n . t + sum |k . a| + sum |k . b|) / (n . t) over the edges: a corrected
    // flux reads its two cells and their eight neighbours with these weights, so an explicit
    // limit set for the two-point flux holds with the correction once divided by it.
    double correctionfactor() const { return correctionfactor_; }

private:
    void buildcorrection();
    const grid::TectonicGrid* grid_ = nullptr;
    double reference_area_m2_ = 0, correctionfactor_ = 1;
    std::vector<double> area_, scale_, width_, weight_;
    std::vector<mesh::Vec3> east_, south_; // the stencil axes
    std::vector<std::array<std::array<double, 2>, 4>> gradient_;
    struct Correction { double orthogonal = 1; std::array<double, 4> a{}, b{}; }; // per edge
    std::vector<Correction> correction_;
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
// The largest displacement per update, as a share of each cell's smallest edge.
inline constexpr double displacement_limit = 0.95;

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
// The transport of the latest update, measurement only (W2-E, E3 checklist item 9).
struct TransportDiagnostics {
    uint32_t substeps = 0;                // the stencil's substeps
    double max_lattice_displacement = 0;  // lattice cells per update (the stencil's maximum)
    double max_displacement_edges = 0;    // max over cells of |displacement| / smallest edge
    // |1 - coverage| (the plates' summed coverage after the transport) over rigid interiors:
    // cells whose cell, eight neighbours and two cardinal rings share one plate and have no
    // offsets; by class: face interiors, cells along a face edge, the 24 cube-corner cells.
    double coverage_interior_max = 0, coverage_edge_max = 0, coverage_corner_max = 0;
    double coverage_interior_mean = 0, coverage_edge_mean = 0, coverage_corner_mean = 0;
    uint32_t interior_cells = 0, edge_cells = 0, corner_cells = 0;
};
// Cumulative horizontal strain of the rock: a symmetric 2D tensor on the cell's stencil axes
// (Metric::axis: x toward the east neighbour, y orthogonal toward the south one),
// dimensionless, the sum of the per-update strains (negative is shortening). Measurement only.
struct Strain { double xx = 0, yy = 0, xy = 0; };
// Shortening, m: minus the most compressive principal component of the strain, times
// cell_m (the length lost across one cell spacing, so it scales with the resolution;
// strain = shortening / cell_m); 0 where no axis shortens.
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
// The column processes on the cell table, in the order Model::advance runs them. Each
// conserves the mass it moves in true areas and books it in the budget.
// Lower-crust root flow (n = 3 creep), sub-stepped below the explicit limit of the smallest
// reference face divided by Metric::correctionfactor; dt_myr already scaled by
// margins::Parameters::root_mobility.
void spread_roots(const Metric& metric, std::vector<Cell>& cells, double dt_myr, const Isostasy& isostasy,
                  bool sediment_feedback, Budget& budget);
// Orogenic collapse (margins::Parameters::orogenic_collapse_*), sub-stepped below 0.1 per face
// divided by Metric::correctionfactor.
void collapse_orogens(const Metric& metric, std::vector<Cell>& cells, double dt_myr,
                      const margins::Parameters& parameters, Budget& budget);
// The surface diffusion's explicit number per reference face, transport_diffusivity x
// erosion_strength x dt / reference area (Model::advance rejects above 0.2 divided by
// Metric::correctionfactor).
double surface_diffusion_number(const Metric& metric, double erosion_strength, double dt_myr);
// Conservative downhill transport of the elevation field (one value per cell) at that number.
void diffuse_surface(const Metric& metric, std::vector<Cell>& cells, const std::vector<double>& elevation,
                     double diffusion_number, double submarine_erosion_ratio, Budget& budget);

class Model {
public:
    // cells: one column per grid cell (the grid is kept); subdivision: the transport lattice's
    // subcells per cell edge (1, 2 or 4).
    void initialize(std::shared_ptr<const grid::TectonicGrid> grid, std::vector<Cell> cells,
                    uint32_t transport_subdivision = 1, long seed = 0,
                    const ocean_relief::Parameters& relief = {}, bool preserve_initial = false,
                    double root_kg_m3 = 0, double root_max_thickness_m = 0,
                    double root_delamination_width_m = 0);
    // motion: the update's motion (sphere_transport::Motion): per plate the rotation of the
    // update (radians), per cell its plate (noplate: at rest) and its offsets (the distributed
    // deformation and rift opening, metres; empty: none). A cell's displacement is its
    // plate's rotation of its centre plus its offset and must stay within displacement_limit
    // of its smallest edge.
    void advance(const sphere_transport::Motion& motion, double dt_myr,
                 double erosion_strength, double thermal_strength, double submarine_erosion_ratio = 1,
                 const margins::Forcing* forcing = nullptr, const margins::Parameters& margins = {},
                 const noise_pass::Parameters& noise = {}, bool sediment_feedback = false,
                 const plumes::Source* plume = nullptr, double uplift_rate_window_myr = 10);
    const std::vector<Cell>& cells() const { return cells_; }
    const Metric& metric() const { return metric_; }
    // Each cell's displacement of motion (its plate's rotation as w x r R, plus the offset's
    // tangent part), metres.
    static std::vector<Displacement> displacement(const grid::TectonicGrid& grid,
                                                  const sphere_transport::Motion& motion);
    // Rock uplift rate per cell, m/Myr, signed: an exponential mean (e-folding time
    // uplift_rate_window_myr) of the elevation change of the tectonic sub-steps
    // (transport, subduction overlap, rifting, trench flexure, arcs, plumes, roots,
    // thermal subsidence), before the diffusion, rivers and noise. In the plate's frame
    // (semi-Lagrangian: the elevation before the step and the previous mean are carried
    // with the rock by the engine's own transport on the subcells), so plate motion over
    // relief is neither uplift nor subsidence and the mean follows the rock. Exact (to
    // rounding) for a rigid motion when a single column field varies (land); at plate
    // boundaries the transport's margin rule moves every field by the rigid reference's
    // fractions and the carried values are a bounded approximation. Converging rock mixes;
    // a divergent cell keeps the rock that stayed on its own side; new crust has no
    // history of its own. Zero at start.
    const std::vector<double>& upliftRate() const { return uplift_rate_m_myr_; }
    // Cumulative strain per cell (multiscale plan S5.1), in the plate's frame. Each update
    // first rotates the rock's tensor by the local spin and adds the strain of the update's
    // displacement field u (metres per update, on the cell's stencil axes): E = sym(grad u),
    // W = (duy/dx - dux/dy) / 2, the gradient by least squares over the cardinal neighbours
    // on the same plate (Metric::gradientweights' offsets restricted to them; along one line
    // when they are collinear; zero without any), so a velocity jump to another plate
    // (subduction, ridge, transform) is never strain. Distributed deformation and rift
    // opening inside a plate count. Exact for a velocity linear around a cell; a rigid
    // rotation gives E = 0 to the rotation's curvature. The tensor is then carried with the
    // rock as the uplift rate is, as a 3-D tensor (so its components turn with the frames
    // across face edges); the ocean floor created in an opening has zero strain, so a
    // divergent cell keeps the tensor of the rock that stayed on the rest of its area only.
    // Small-strain sum; cumulative over the run, no decay; zero at start.
    const std::vector<Strain>& strain() const { return strain_; }
    const Budget& budget() const { return budget_; }
    const TransportDiagnostics& transportDiagnostics() const { return transport_diagnostics_; }
private:
    // carried: non-negative fields on the subcells, moved in place with the material.
    std::vector<Cell> transportColumns(const sphere_transport::Motion& motion, double dt_myr,
                                       std::vector<std::vector<double>>& carried);
    void reconcileTransport(const std::vector<Cell>& evolved);
    void measureCoverage(const std::vector<Cell>& moved, const sphere_transport::Motion& motion);
    void routeSediment(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                       const margins::Parameters& parameters, const margins::Forcing* forcing);
    void routeRiversImplicit(std::vector<Cell>& cells, double dt_myr, double thermal_strength,
                             const margins::Parameters& parameters, const margins::Forcing* forcing);
    std::shared_ptr<const grid::TectonicGrid> grid_;
    Metric metric_;
    std::shared_ptr<const sphere_transport::Lattice> lattice_;
    std::vector<std::vector<grid::CellId>> children_;  // per cell its subcells
    std::vector<uint8_t> cellclass_;                   // 0 face interior, 1 along a face edge, 2 cube corner
    std::vector<std::vector<grid::CellId>> neighbourhood_;  // eight neighbours and two cardinal rings
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
    TransportDiagnostics transport_diagnostics_;
    std::vector<double> uplift_rate_m_myr_; // Measurement only; never feeds back into the columns.
    std::vector<Strain> strain_; // Measurement only; never feeds back into the columns.
};
} // namespace tectonics::columns
#endif
