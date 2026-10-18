#ifndef TECTONIC_DYNAMICS_HPP
#define TECTONIC_DYNAMICS_HPP

#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace tectonics::dynamics {
using mesh::Vec3;

struct Vector { double x = 0.0, y = 0.0; };

// Reduced Cartesian force balance, version 1. SI units throughout the solver.
// Basal Couette drag: eta / asthenosphere thickness = 1e20 / 200000 Pa s/m.
// Boundary dashpots are explicit effective rheology assumptions, not Earth fits.
inline constexpr double kBasalDrag = 5.0e14;
inline constexpr double kOceanBoundaryDrag = 1.0e21; // N s/m², times boundary length
inline constexpr double kContinentalBoundaryDrag = 8.0e21;
inline constexpr double kShearBoundaryDrag = 5.0e20;
inline constexpr double kSecondsPerMyr = 3.15576e13;

struct Plate {
    double area_m2 = 0.0;
    Vector mantle_velocity_m_s;
};

struct Contact {
    uint32_t a = 0, b = 0;
    Vector normal; // unit outward normal from a toward b
    double length_m = 0.0;
    bool continental_a = false, continental_b = false;
    double age_a_myr = 0.0, age_b_myr = 0.0;
    // Previous normal convergence in m/s. Local activity proxy, not slab history.
    double closing_m_s = 0.0;
    double mountain_force_n_m = 0.0; // Excess column GPE, outward on both flanks.
    double activation_speed_m_s = 0.0; // Continuous activation scale; zero preserves the historical switch.
    double slab_pull_scale = 1.0; // Thermal density anomaly multiplier, shared with surface compensation.
    int descending_side = -1; // -1: historical age choice; 0/1: shared persistent polarity.
    double ridge_activity = -1, slab_activity = -1; // Previous attached boundary state; negative selects instantaneous activity.
    double activity_response = 1; // Implicit history update dt/(tau+dt).
};

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
};

struct Solution {
    std::vector<Vector> velocity_m_s;
    Diagnostics diagnostics;
};

// Smooth, seeded, periodic prescribed mantle flow. This is a forcing field,
// not a mantle-convection solution. It changes on a 600 Myr period.
Vector mantle_flow(double x_fraction, double y_fraction, double time_myr,
                   long seed, double speed_cm_yr);
double slab_pull_per_length(double age_myr);
// Effective GPE transmitted from a deforming strip on one side of a contact.
// Samples only the same connected plate along its inward grid-face normal.
double inland_pressure(uint32_t width,uint32_t height,double cell_m,
    const std::vector<uint32_t>& owners,const std::vector<double>& pressure,
    uint32_t cell,uint32_t axis,int direction,double range_m,const grid::TectonicGrid* grid = nullptr);
Solution solve(const std::vector<Plate>& plates, const std::vector<Contact>& contacts);
Solution solve_consistent(const std::vector<Plate>& plates, std::vector<Contact> contacts,
                          std::vector<Vector> velocity);

// ---- The sphere (W2-E, step E2C; docs/architecture/terrain-roadmap.md, "W2-E") ----
// The same force balance as torques about the planet's centre: per plate a rotation
// vector w (rad/s, 3 degrees of freedom), a cell at unit vector c moving w x R c. The
// planar functions above stay until E3 wires the cube and deletes them.

// The mantle flow at a point: its surface velocity and its local spin, half the radial
// vorticity (rad/s; a solid-body rotation W has spin W . c everywhere).
struct MantleSample {
    Vec3 velocity_m_s;
    double spin_rad_s = 0.0;
};

// Smooth, seeded, periodic prescribed mantle flow on the sphere, v = grad chi + r x grad psi
// with chi and psi 3-D fractal noise at unit vectors. A forcing field, not a mantle-convection
// solution. Kept from the torus field: the 600 Myr clock, the seed's phase, the split of the
// mean square speed between the divergent and the rotational part (0.5 : 0.35^2), the
// divergent pattern's small wobble and the rotational one's drift by one pattern per period,
// and the RMS speed, sqrt(0.6225) x speed_cm_yr (measured per seed over a fixed point set).
class SphereMantleFlow {
public:
    SphereMantleFlow(long seed, double speed_cm_yr);
    // unit: a unit vector; the velocity is tangent there. radius_m scales the spin only.
    MantleSample at(const Vec3& unit, double time_myr, double radius_m) const;
    double poloidal_scale() const { return poloidal_; }   // m/s per unit surface gradient of chi
    double toroidal_scale() const { return toroidal_; }   // m/s per unit surface gradient of psi
private:
    struct Fields;
    std::shared_ptr<const Fields> fields_;
    double phase_ = 0.0, poloidal_ = 0.0, toroidal_ = 0.0;
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
// The geometry of a grid edge (the midpoint as the grid computes it, its normal and length).
SphereContact sphere_contact(const grid::TectonicGrid& grid, const grid::Edge& edge, uint32_t a, uint32_t b);
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
// The consistent-motion Newton iteration of solve_consistent with 3P columns.
SphereSolution solve_sphere_consistent(const std::vector<SpherePlate>& plates, std::vector<SphereContact> contacts,
                                       std::vector<Vec3> rotation, double radius_m);

// inland_pressure on the cube: from cell, walking straight along grid.ahead in the
// direction given in cell's frame, over the great-circle distances between centres.
double inland_pressure(const grid::TectonicGrid& grid, const std::vector<uint32_t>& owners,
    const std::vector<double>& pressure, grid::CellId cell, std::uint8_t direction, double range_m);

} // namespace tectonics::dynamics
#endif
