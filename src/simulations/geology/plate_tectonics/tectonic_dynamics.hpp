#ifndef TECTONIC_DYNAMICS_HPP
#define TECTONIC_DYNAMICS_HPP

#include <cstdint>
#include <vector>

namespace tectonics::grid { class TectonicGrid; }
namespace tectonics::dynamics {

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

} // namespace tectonics::dynamics
#endif
