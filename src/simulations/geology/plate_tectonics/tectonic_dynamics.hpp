#ifndef TECTONIC_DYNAMICS_HPP
#define TECTONIC_DYNAMICS_HPP

#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace tectonics::dynamics {
using mesh::Vec3;

// The plate force balance as torques about the planet's centre (W2-E, E2C; the engine's only
// solver since E3). SI units throughout the solver.
// Basal Couette drag: eta / asthenosphere thickness = 1e20 / 200000 Pa s/m.
// Boundary dashpots are explicit effective rheology assumptions, not Earth fits.
inline constexpr double kBasalDrag = 5.0e14;
inline constexpr double kOceanBoundaryDrag = 1.0e21; // N s/m², times boundary length
inline constexpr double kContinentalBoundaryDrag = 8.0e21;
inline constexpr double kShearBoundaryDrag = 5.0e20;
inline constexpr double kSecondsPerMyr = 3.15576e13;

struct Diagnostics {
    bool enabled = false;
    double relative_force_residual = 0.0;
    double driving_power_w = 0.0;
    double dissipation_w = 0.0;
    double max_speed_cm_yr = 0.0;
    double mean_speed_cm_yr = 0.0;
    double ridge_length_km = 0.0;
    double subduction_length_km = 0.0;
    double continental_contact_length_km = 0.0;
    double mountain_boundary_force_n = 0.0; // Sum of one-flank force magnitudes.
    uint32_t nonlinear_iterations = 0;
    double velocity_consistency_cm_yr = 0;
    double max_velocity_consistency_cm_yr = 0;
    uint64_t persistent_face_steps = 0, reversing_face_steps = 0;
    uint64_t strong_face_steps = 0, strong_reversing_face_steps = 0;
    // Measures of the latest solve (sphere_measures; zero until the engine records them).
    double rms_speed_cm_yr = 0;       // area-RMS surface speed over all plates, sqrt(sum w^T D w / sum A)
    double net_rotation_cm_yr = 0;    // |omega_net| R of the area-weighted lithosphere rotation
    double condition_number = 0;      // largest over smallest eigenvalue of the torque balance
    double spin_share_mean = 0, spin_share_max = 0;  // |w . centroid| / |w| over the plates
};

// The slab pull per length of an attached slab of oceanic lithosphere of an age, N/m.
double slab_pull_per_length(double age_myr);

// ---- The sphere (W2-E, step E2C; docs/architecture/terrain-roadmap.md, "W2-E") ----
// Per plate a rotation vector w (rad/s, 3 degrees of freedom), a cell at unit vector c
// moving w x R c.

// The mantle flow at a point: its surface velocity and its local spin, half the radial
// vorticity (rad/s; a solid-body rotation W has spin W . c everywhere).
struct MantleSample {
    Vec3 velocity_m_s;
    double spin_rad_s = 0.0;
};

// The shape of the mantle flow (scenario keys mantle_flow_*): the noise's lattice spacing
// (wavelength) and octaves, the divergent pattern's wobble (radians of its cycle) and the
// rotational pattern's drift (patterns per 600 Myr), and the rotational part's amplitude
// relative to the divergent one (the mean square speed splits 0.5 : shear_ratio^2).
struct MantleFlowShape {
    double wavelength_km = 6371.0;
    int octaves = 2;
    double wobble_rad = 0.3;
    double drift_cycles = 1.0;
    double shear_ratio = 0.35;
};
void validate(const MantleFlowShape& shape);

// Smooth, seeded, periodic prescribed mantle flow on the sphere, v = grad chi + r x grad psi
// with chi and psi 3-D fractal noise at unit vectors. A forcing field, not a mantle-convection
// solution. Kept from the torus field: the 600 Myr clock, the seed's phase, the split of the
// mean square speed between the divergent and the rotational part (0.5 : 0.35^2 by default),
// the divergent pattern's small wobble and the rotational one's drift by one pattern per
// period, and the RMS speed, sqrt(0.5 + shear_ratio^2) x speed_cm_yr (measured per seed over a
// fixed point set).
class SphereMantleFlow {
public:
    SphereMantleFlow(long seed, double speed_cm_yr, const MantleFlowShape& shape = {}, double radius_km = 6371.0);
    // unit: a unit vector; the velocity is tangent there. radius_m scales the spin only.
    MantleSample at(const Vec3& unit, double time_myr, double radius_m) const;
    double poloidal_scale() const { return poloidal_; }   // m/s per unit surface gradient of chi
    double toroidal_scale() const { return toroidal_; }   // m/s per unit surface gradient of psi
private:
    struct Fields;
    std::shared_ptr<const Fields> fields_;
    double phase_ = 0.0, poloidal_ = 0.0, toroidal_ = 0.0;
    double lattice_ = 1.0;   // the noise's lattice unit on the unit sphere (wavelength / radius)
    double wobble_ = 0.0, drift_ = 0.0;  // offsets' radii on the unit sphere
};

// A plate's basal-drag moments, summed over its cells (kBasalDrag times them is the drag).
// Each cell is a disc of its area: besides the lever arm R c it resists spin about its own
// centre with the polar moment A^2 / (2 pi), so the drag of any plate, one cell included,
// is positive definite.
struct SpherePlate {
    double area_m2 = 0.0;
    std::array<double, 9> drag_m4{};   // sum A R^2 (I - c c^T) + A^2 / (2 pi) c c^T, row-major
    Vec3 mantle_m4_s;                  // sum A R c x v + A^2 / (2 pi) s c (v, s the mantle sample)
    Vec3 centre_sum_m2;                // sum A c; its direction is the plate's centroid
};
void add_sphere_cell(SpherePlate& plate, const Vec3& unit, double area_m2, double radius_m,
                     const MantleSample& mantle);
Vec3 sphere_centroid(const SpherePlate& plate);
// Per owner 0 .. count - 1 of the cube grid's cells (owners >= count are skipped).
std::vector<SpherePlate> sphere_plates(const grid::TectonicGrid& grid, const std::vector<uint32_t>& owners,
                                       uint32_t count, const std::vector<MantleSample>& mantle);

// One edge between plates a and b. The scalar fields mean what they mean in Contact.
struct SphereContact {
    uint32_t a = 0, b = 0;
    Vec3 position;  // unit vector of the edge's midpoint
    Vec3 normal;    // unit tangent there, from a toward b
    double length_m = 0.0;
    bool continental_a = false, continental_b = false;
    double age_a_myr = 0.0, age_b_myr = 0.0;
    double closing_m_s = 0.0;
    double mountain_force_n_m = 0.0;
    double activation_speed_m_s = 0.0;
    double slab_pull_scale = 1.0;
    int descending_side = -1;
    double ridge_activity = -1, slab_activity = -1;
    double activity_response = 1;
};
// The geometry of a grid edge (its midpoint, normal and length), between the owners of its
// two cells: a = owners[edge.a], b = owners[edge.b], so the normal points from a toward b.
SphereContact sphere_contact(const grid::TectonicGrid& grid, const grid::Edge& edge,
                             const std::vector<uint32_t>& owners);
// Normal convergence ((w_a - w_b) x R x) . n, m/s.
double sphere_closing_m_s(const SphereContact& contact, const Vec3& rotation_a, const Vec3& rotation_b,
                          double radius_m);

// The 3P x 3P torque balance K w = t: K row-major (N m s), t in N m.
struct SphereSystem {
    std::vector<double> matrix, rhs;
    Diagnostics diagnostics;
};
SphereSystem sphere_system(const std::vector<SpherePlate>& plates, const std::vector<SphereContact>& contacts,
                           double radius_m);
struct SphereSolution {
    std::vector<Vec3> rotation_rad_s;
    Diagnostics diagnostics;  // speeds at the plates' centroids
};
SphereSolution solve_sphere(const std::vector<SpherePlate>& plates, const std::vector<SphereContact>& contacts,
                            double radius_m);
// The consistent-motion Newton iteration: the activity-dependent forces re-evaluated at fixed
// geometry until the rotations reproduce themselves (3P columns).
SphereSolution solve_sphere_consistent(const std::vector<SpherePlate>& plates, std::vector<SphereContact> contacts,
                                       std::vector<Vec3> rotation, double radius_m);
// A plate's rotation under its basal drag alone (the mantle forcing's own motion of the plate).
Vec3 drag_rotation(const SpherePlate& plate, double radius_m);
// The measures of a solution (Diagnostics' last five fields): the area-RMS speed, the net
// rotation (lever moments of the plates, sum over the sphere 8 pi R^4 / 3), the torque balance's
// condition number (cyclic Jacobi) and the spin shares.
void sphere_measures(const std::vector<SpherePlate>& plates, const std::vector<SphereContact>& contacts,
                     const std::vector<Vec3>& rotation, double radius_m, Diagnostics& diagnostics);

// Effective GPE transmitted from a deforming strip on one side of a contact: from cell, walking
// straight along grid.ahead in the direction given in cell's frame, over the great-circle
// distances between centres, on the same plate only.
double inland_pressure(const grid::TectonicGrid& grid, const std::vector<uint32_t>& owners,
    const std::vector<double>& pressure, grid::CellId cell, std::uint8_t direction, double range_m);

} // namespace tectonics::dynamics
#endif
