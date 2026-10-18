#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"
#include "core/deterministic_random.hpp"
#include "simulations/landscape/spherical_noise.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace tectonics::dynamics {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double cm_yr_per_m_s = 3.15576e9;

double dot(Vector a, Vector b) { return a.x * b.x + a.y * b.y; }

// Deterministic Cholesky solve of the positive definite drag matrix.
std::vector<double> solve_spd(std::vector<double> matrix, std::vector<double> rhs) {
    const size_t n = rhs.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double value = matrix[i * n + j];
            for (size_t k = 0; k < j; ++k) value -= matrix[i * n + k] * matrix[j * n + k];
            if (i == j && (!(value > 0.0) || !std::isfinite(value)))
                throw std::runtime_error("Plate drag matrix is not positive definite");
            matrix[i * n + j] = i == j ? std::sqrt(value) : value / matrix[j * n + j];
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) rhs[i] -= matrix[i * n + j] * rhs[j];
        rhs[i] /= matrix[i * n + i];
    }
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i + 1; j < n; ++j) rhs[i] -= matrix[j * n + i] * rhs[j];
        rhs[i] /= matrix[i * n + i];
    }
    return rhs;
}
} // namespace

Vector mantle_flow(double x, double y, double time_myr, long seed, double speed_cm_yr) {
    // Seed changes phases only; no mutable RNG or dependence on update count.
    const auto bits = static_cast<uint32_t>(seed) * 1664525U + 1013904223U;
    const double phase = 2.0 * pi * static_cast<double>(bits) / 4294967296.0;
    const double t = 2.0 * pi * time_myr / 600.0;
    const double a = 2.0 * pi * x + phase + .3 * std::sin(t);
    const double b = 2.0 * pi * y - .7 * phase + .3 * std::cos(t);
    const double speed = speed_cm_yr / cm_yr_per_m_s;
    return {speed * (std::sin(a) * std::cos(b) + .35 * std::sin(2.0 * b + t)),
            speed * (std::cos(a) * std::sin(b) - .35 * std::sin(2.0 * a - t))};
}

double inland_pressure(uint32_t w,uint32_t h,double cell_m,const std::vector<uint32_t>& owners,
        const std::vector<double>& pressure,uint32_t cell,uint32_t axis,int direction,double range_m,
        const grid::TectonicGrid* grid) {
    grid::requireplanar(grid, "dynamics::inland_pressure");
    if (w<2 || h<2 || owners.size()!=static_cast<size_t>(w)*h || pressure.size()!=owners.size() ||
        cell>=owners.size() || axis>1 || (direction!=1 && direction!=-1) ||
        !std::isfinite(cell_m) || cell_m<=0 || !std::isfinite(range_m) || range_m<0)
        throw std::runtime_error("Invalid inland pressure strip");
    if (range_m==0) return pressure[cell];
    double sum=0,weight=0;
    const auto owner=owners[cell];
    // Integral of the parabolic shortening profile: peak halfway inland,
    // zero at both ends. This averages force, never terrain or crust thickness.
    const auto cdf=[](double u) { u=std::min(1.,u);return u*u*(3-2*u); };
    const uint32_t limit=axis==0?w:h;
    const double reach=std::min(range_m,limit*cell_m); // Never wrap around and count a strip twice.
    for (uint32_t step=0;step<limit && step*cell_m<reach;++step) {
        if (owners[cell]!=owner) break;
        const double k=cdf((step+1)*cell_m/reach)-cdf(step*cell_m/reach);
        sum+=k*pressure[cell];weight+=k;
        const uint32_t x=cell%w,y=cell/w;
        cell=axis==0?y*w+(x+w+direction)%w:((y+h+direction)%h)*w+x;
    }
    return sum/weight;
}

double slab_pull_per_length(double age_myr) {
    // Thermal buoyancy deficit per area (kg/m²) inferred consistently from
    // plate-cooling subsidence. Effective attached slab length 300 km; only
    // 20% transmitted horizontally after unresolved bending/mantle resistance.
    const double buoyancy_deficit = (3300.0 - 1030.0) * thermal::subsidence_m(age_myr);
    return .2 * buoyancy_deficit * 9.80665 * 300000.0;
}

Solution solve(const std::vector<Plate>& plates, const std::vector<Contact>& contacts) {
    const size_t n = plates.size() * 2;
    Solution result;
    result.velocity_m_s.resize(plates.size());
    if (n == 0) return result;
    auto& d = result.diagnostics;
    d.enabled = true;
    std::vector<double> matrix(n * n, 0.0), rhs(n, 0.0);
    for (size_t p = 0; p < plates.size(); ++p) {
        const auto& plate = plates[p];
        if (!(plate.area_m2 > 0.0) || !std::isfinite(plate.area_m2) ||
            !std::isfinite(plate.mantle_velocity_m_s.x) || !std::isfinite(plate.mantle_velocity_m_s.y))
            throw std::runtime_error("Invalid plate area or mantle velocity in force balance");
        const double drag = kBasalDrag * plate.area_m2;
        matrix[(2 * p) * n + 2 * p] = drag;
        matrix[(2 * p + 1) * n + 2 * p + 1] = drag;
        rhs[2 * p] = drag * plate.mantle_velocity_m_s.x;
        rhs[2 * p + 1] = drag * plate.mantle_velocity_m_s.y;
    }
    for (const auto& edge : contacts) {
        if (edge.a >= plates.size() || edge.b >= plates.size() || edge.a == edge.b ||
            !(edge.length_m > 0.0) || !std::isfinite(edge.length_m) ||
            !std::isfinite(edge.closing_m_s) || !std::isfinite(edge.age_a_myr) ||
            !std::isfinite(edge.mountain_force_n_m) || edge.mountain_force_n_m < 0 ||
            !std::isfinite(edge.activation_speed_m_s) || edge.activation_speed_m_s < 0 ||
            !std::isfinite(edge.slab_pull_scale) || edge.slab_pull_scale < 0 ||
            edge.descending_side < -1 || edge.descending_side > 1 ||
            !std::isfinite(edge.activity_response) || edge.activity_response<0 || edge.activity_response>1 ||
            !std::isfinite(edge.ridge_activity) || edge.ridge_activity < -1 || edge.ridge_activity>1 ||
            !std::isfinite(edge.slab_activity) || edge.slab_activity < -1 || edge.slab_activity>1 ||
            !std::isfinite(edge.age_b_myr) || std::abs(dot(edge.normal, edge.normal) - 1.0) > 1e-6 ||
            !std::isfinite(edge.normal.x) || !std::isfinite(edge.normal.y))
            throw std::runtime_error("Invalid contact in plate force balance");
        const bool continental = edge.continental_a && edge.continental_b;
        const double normal_drag = continental ? kContinentalBoundaryDrag : kOceanBoundaryDrag;
        const double normal[] = {edge.normal.x, edge.normal.y};
        // Symmetric pairwise viscous resistance; equal/opposite boundary forces.
        for (size_t i = 0; i < 2; ++i) for (size_t j = 0; j < 2; ++j) {
            const double drag = edge.length_m * (kShearBoundaryDrag * (i == j) +
                (normal_drag - kShearBoundaryDrag) * normal[i] * normal[j]);
            matrix[(2 * edge.a + i) * n + 2 * edge.a + j] += drag;
            matrix[(2 * edge.b + i) * n + 2 * edge.b + j] += drag;
            matrix[(2 * edge.a + i) * n + 2 * edge.b + j] -= drag;
            matrix[(2 * edge.b + i) * n + 2 * edge.a + j] -= drag;
        }
        if (continental) d.continental_contact_length_km += edge.length_m / 1000.0;
        if (edge.mountain_force_n_m > 0) {
            const double force = edge.mountain_force_n_m * edge.length_m;
            for (size_t axis = 0; axis < 2; ++axis) {
                rhs[2 * edge.a + axis] -= force * normal[axis];
                rhs[2 * edge.b + axis] += force * normal[axis];
            }
            d.mountain_boundary_force_n += force;
        }
        // A tiny change of normal velocity must not switch on a finite force.
        // This is a continuous activity proxy, not a resolved attached-slab history.
        const double threshold = edge.activation_speed_m_s > 0 ? 0 : 1e-12;
        const double activation = edge.activation_speed_m_s > 0
            ? std::abs(edge.closing_m_s) / std::hypot(edge.closing_m_s, edge.activation_speed_m_s) : 1;
        const double ridge_target=edge.closing_m_s < -threshold?activation:0;
        const double slab_target=edge.closing_m_s > threshold?activation:0;
        const double ridge=edge.ridge_activity<0?ridge_target:
            edge.ridge_activity+edge.activity_response*(ridge_target-edge.ridge_activity);
        const double slab=edge.slab_activity<0?slab_target:
            edge.slab_activity+edge.activity_response*(slab_target-edge.slab_activity);
        if (ridge>0 && !edge.continental_a && !edge.continental_b) {
            // Integrated ridge gravitational spreading, 2e12 N/m per flank.
            const double force = 2.0e12 * edge.length_m * ridge;
            for (size_t axis = 0; axis < 2; ++axis) {
                rhs[2 * edge.a + axis] -= force * normal[axis];
                rhs[2 * edge.b + axis] += force * normal[axis];
            }
            d.ridge_length_km += edge.length_m / 1000.0;
        }
        if (slab>0 && !continental) {
            // Continental lithosphere never receives slab pull. At ocean/ocean
            // contacts the older side sinks; exact-age ties use no slab force.
            const bool a_sinks = !edge.continental_a &&
                (edge.continental_b || edge.descending_side == 0 ||
                 (edge.descending_side < 0 && edge.age_a_myr > edge.age_b_myr));
            const bool b_sinks = !edge.continental_b &&
                (edge.continental_a || edge.descending_side == 1 ||
                 (edge.descending_side < 0 && edge.age_b_myr > edge.age_a_myr));
            if (a_sinks || b_sinks) {
                const uint32_t p = a_sinks ? edge.a : edge.b;
                const double force = slab_pull_per_length(a_sinks ? edge.age_a_myr : edge.age_b_myr)
                    * edge.length_m * (a_sinks ? 1.0 : -1.0) * slab * edge.slab_pull_scale;
                for (size_t axis = 0; axis < 2; ++axis) rhs[2 * p + axis] += force * normal[axis];
                d.subduction_length_km += edge.length_m / 1000.0;
            }
        }
    }
    const auto velocity = solve_spd(matrix, rhs);
    double residual_sq = 0.0, force_sq = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double drag_force = 0.0;
        for (size_t j = 0; j < n; ++j) drag_force += matrix[i * n + j] * velocity[j];
        residual_sq += (drag_force - rhs[i]) * (drag_force - rhs[i]);
        force_sq += rhs[i] * rhs[i];
        d.driving_power_w += rhs[i] * velocity[i];
        d.dissipation_w += drag_force * velocity[i];
    }
    d.relative_force_residual = std::sqrt(residual_sq / std::max(1.0, force_sq));
    for (size_t p = 0; p < plates.size(); ++p) {
        result.velocity_m_s[p] = {velocity[2 * p], velocity[2 * p + 1]};
        const double speed = std::hypot(velocity[2 * p], velocity[2 * p + 1]) * cm_yr_per_m_s;
        d.max_speed_cm_yr = std::max(d.max_speed_cm_yr, speed);
        d.mean_speed_cm_yr += speed / plates.size();
    }
    return result;
}

Solution solve_consistent(const std::vector<Plate>& plates, std::vector<Contact> contacts,
                          std::vector<Vector> velocity) {
    if (velocity.size()!=plates.size()) throw std::runtime_error("Invalid nonlinear motion seed");
    // Re-evaluate the activity-dependent forces at fixed geometry, before any
    // material moves. Relaxation is numerical iteration, not plate inertia.
    const size_t n=2*velocity.size();
    const auto evaluate=[&](const std::vector<Vector>& v) {
        for (auto& edge:contacts)
            edge.closing_m_s=dot({v[edge.a].x-v[edge.b].x,v[edge.a].y-v[edge.b].y},edge.normal);
        return solve(plates,contacts);
    };
    const auto residual=[&](const std::vector<Vector>& v,const Solution& s) {
        std::vector<double> r(n);
        for (size_t p=0;p<v.size();++p) {
            r[2*p]=(s.velocity_m_s[p].x-v[p].x)*cm_yr_per_m_s;
            r[2*p+1]=(s.velocity_m_s[p].y-v[p].y)*cm_yr_per_m_s;
        }
        return r;
    };
    const auto norm=[](const std::vector<double>& r) {
        double sum=0;for (double value:r) sum+=value*value;return std::sqrt(sum);
    };
    for (uint32_t iteration=1;iteration<=96;++iteration) {
        auto result=evaluate(velocity);
        const auto r=residual(velocity,result);
        double error=0;
        for (size_t p=0;p<plates.size();++p)
            error=std::max(error,std::hypot(result.velocity_m_s[p].x-velocity[p].x,
                                          result.velocity_m_s[p].y-velocity[p].y)*cm_yr_per_m_s);
        result.diagnostics.nonlinear_iterations=iteration;
        result.diagnostics.velocity_consistency_cm_yr=error;
        if (error<1e-4) return result;
        // Activity feedback is not always a contractive fixed point. A small
        // dense, damped Newton solve handles those contacts without timestep lag.
        std::vector<double> jac(n*n),step=r;
        constexpr double epsilon=1e-4; // cm/yr, finite-difference Jacobian.
        for (size_t j=0;j<n;++j) {
            auto perturbed=velocity;
            (j%2?perturbed[j/2].y:perturbed[j/2].x)+=epsilon/cm_yr_per_m_s;
            const auto q=residual(perturbed,evaluate(perturbed));
            for (size_t i=0;i<n;++i) jac[i*n+j]=(q[i]-r[i])/epsilon;
        }
        for (auto& value:step) value=-value;
        bool singular=false;
        for (size_t k=0;k<n;++k) {
            size_t pivot=k;
            for (size_t i=k+1;i<n;++i) if (std::abs(jac[i*n+k])>std::abs(jac[pivot*n+k])) pivot=i;
            if (std::abs(jac[pivot*n+k])<1e-12) {singular=true;break;}
            for (size_t j=k;j<n;++j) std::swap(jac[k*n+j],jac[pivot*n+j]);
            std::swap(step[k],step[pivot]);
            for (size_t i=k+1;i<n;++i) {
                const double factor=jac[i*n+k]/jac[k*n+k];
                for (size_t j=k;j<n;++j) jac[i*n+j]-=factor*jac[k*n+j];
                step[i]-=factor*step[k];
            }
        }
        if (!singular) for (size_t i=n;i-->0;) {
            for (size_t j=i+1;j<n;++j) step[i]-=jac[i*n+j]*step[j];
            step[i]/=jac[i*n+i];
        }
        double best=norm(r);
        auto next=velocity;
        for (int method=0;method<2;++method) {
            if (method==0 && singular) continue;
            for (double damping=1;damping>=1./256;damping*=.5) {
                auto candidate=velocity;
                for (size_t p=0;p<velocity.size();++p) {
                    candidate[p].x+=damping*(method?r[2*p]:step[2*p])/cm_yr_per_m_s;
                    candidate[p].y+=damping*(method?r[2*p+1]:step[2*p+1])/cm_yr_per_m_s;
                }
                const double candidate_error=norm(residual(candidate,evaluate(candidate)));
                if (candidate_error<best) {best=candidate_error;next=std::move(candidate);break;}
            }
        }
        if (best>=norm(r) || iteration==96)
            throw std::runtime_error("Activity-dependent plate motion did not converge; residual cm/yr="+std::to_string(error));
        velocity=std::move(next);
    }
    throw std::runtime_error("Unreachable nonlinear motion state");
}

// ---- The sphere ----
namespace {
// The shape of the spherical forcing, each kept from the torus field mantle_flow above:
// its mean square speed is 0.5 + 0.35^2 = 0.6225 of speed^2 (the divergent part
// sin a cos b, cos a sin b and the shear part 0.35 sin(..)), the divergent pattern
// wobbles by 0.3 rad of its 2 pi cycle and the shear pattern drifts by one cycle per
// 600 Myr. The noise has its lattice spacing at the planet's radius (patterns of about
// a radius, as the torus' one cycle per world) and two octaves.
constexpr double kDivergentMeanSquare = 0.5, kRotationalMeanSquare = 0.35 * 0.35;
constexpr double kDivergentWobble = 0.3 / (2.0 * pi), kRotationalDrift = 1.0 / (2.0 * pi); // lattice units
constexpr double kNoiseStep = 1e-4;           // lattice units, central differences
constexpr int kNormalisationPoints = 2048, kNormalisationTimes = 12;

landscape::FractalNoiseField mantle_noise(std::uint64_t seed) {
    landscape::FractalNoise noise;
    noise.seed = seed;
    noise.radius_km = 1.0;
    noise.wavelength_km = 1.0;
    noise.octaves = 2;
    noise.persistence = 0.5;
    noise.kernel = landscape::NoiseKernel::perlin;
    noise.rotateoctaves = 0;
    return landscape::FractalNoiseField(noise);
}

double noise_at(const landscape::FractalNoiseField& field, const Vec3& q) { return field.sum(q.x, q.y, q.z); }

// The surface gradient on the unit sphere at p of the field sampled at p + offset, and its
// surface Laplacian (Laplacian in space less the radial terms d2/dr2 + 2 d/dr at r = 1).
struct SurfaceDerivatives { Vec3 gradient; double laplacian = 0.0; };
SurfaceDerivatives surface_derivatives(const landscape::FractalNoiseField& field, const Vec3& p,
                                       const Vec3& offset, bool laplacian) {
    const double h = kNoiseStep;
    const Vec3 q = p + offset;
    const Vec3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double plus[3], minus[3];
    for (int i = 0; i < 3; ++i) {
        plus[i] = noise_at(field, q + h * axes[i]);
        minus[i] = noise_at(field, q - h * axes[i]);
    }
    const Vec3 gradient{(plus[0] - minus[0]) / (2 * h), (plus[1] - minus[1]) / (2 * h), (plus[2] - minus[2]) / (2 * h)};
    const double radial = mesh::dot(gradient, p);
    SurfaceDerivatives result{gradient - radial * p, 0.0};
    if (laplacian) {
        const double centre = noise_at(field, q);
        double space = 0;
        for (int i = 0; i < 3; ++i) space += (plus[i] - 2 * centre + minus[i]) / (h * h);
        const double outer = noise_at(field, q + h * p), inner = noise_at(field, q - h * p);
        result.laplacian = space - (outer - 2 * centre + inner) / (h * h) - 2 * radial;
    }
    return result;
}

double time_angle(double time_myr, double phase) { return 2.0 * pi * time_myr / 600.0 + phase; }
Vec3 divergent_offset(double angle) { return {kDivergentWobble * std::cos(angle), kDivergentWobble * std::sin(angle), 0}; }
Vec3 rotational_offset(double angle) { return {kRotationalDrift * std::cos(angle), 0, kRotationalDrift * std::sin(angle)}; }

void add3(std::vector<double>& matrix, size_t n, size_t row, size_t column, double sign, const double block[9]) {
    for (size_t i = 0; i < 3; ++i) for (size_t j = 0; j < 3; ++j)
        matrix[(3 * row + i) * n + 3 * column + j] += sign * block[3 * i + j];
}
void add_torque(std::vector<double>& rhs, uint32_t plate, double force, const Vec3& arm) {
    rhs[3 * plate] += force * arm.x;
    rhs[3 * plate + 1] += force * arm.y;
    rhs[3 * plate + 2] += force * arm.z;
}
bool finite3(const Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
} // namespace

struct SphereMantleFlow::Fields {
    landscape::FractalNoiseField divergent, rotational;
};

SphereMantleFlow::SphereMantleFlow(long seed, double speed_cm_yr) {
    if (!std::isfinite(speed_cm_yr) || speed_cm_yr < 0)
        throw std::runtime_error("Mantle flow speed must be finite and nonnegative");
    // The torus' seed phase; the seed also picks the two noise fields.
    const auto bits = static_cast<uint32_t>(seed) * 1664525U + 1013904223U;
    phase_ = 2.0 * pi * static_cast<double>(bits) / 4294967296.0;
    fields_ = std::make_shared<const Fields>(Fields{mantle_noise(deterministiccombine(bits, 1)),
                                                    mantle_noise(deterministiccombine(bits, 2))});
    // Mean square surface gradients over a Fibonacci point set and the 600 Myr cycle.
    double divergent = 0, rotational = 0;
    constexpr double golden = 2.39996322972865332; // pi (3 - sqrt 5)
    for (int t = 0; t < kNormalisationTimes; ++t) {
        const double angle = time_angle(600.0 * t / kNormalisationTimes, phase_);
        for (int i = 0; i < kNormalisationPoints; ++i) {
            const double z = 1.0 - (2.0 * i + 1.0) / kNormalisationPoints, r = std::sqrt(1.0 - z * z);
            const Vec3 p{r * std::cos(golden * i), r * std::sin(golden * i), z};
            const Vec3 g = surface_derivatives(fields_->divergent, p, divergent_offset(angle), false).gradient;
            const Vec3 s = surface_derivatives(fields_->rotational, p, rotational_offset(angle), false).gradient;
            divergent += mesh::dot(g, g);
            rotational += mesh::dot(s, s);
        }
    }
    const double samples = static_cast<double>(kNormalisationPoints) * kNormalisationTimes;
    const double speed = speed_cm_yr / cm_yr_per_m_s;
    poloidal_ = speed * std::sqrt(kDivergentMeanSquare * samples / divergent);
    toroidal_ = speed * std::sqrt(kRotationalMeanSquare * samples / rotational);
}

MantleSample SphereMantleFlow::at(const Vec3& unit, double time_myr, double radius_m) const {
    if (!finite3(unit) || std::abs(mesh::dot(unit, unit) - 1.0) > 1e-9 || !std::isfinite(time_myr) ||
        !(radius_m > 0) || !std::isfinite(radius_m))
        throw std::runtime_error("Invalid mantle flow sample point");
    const double angle = time_angle(time_myr, phase_);
    const auto chi = surface_derivatives(fields_->divergent, unit, divergent_offset(angle), false);
    const auto psi = surface_derivatives(fields_->rotational, unit, rotational_offset(angle), true);
    // v = grad chi + r x grad psi; only the rotational part has radial vorticity, Laplacian psi.
    return {poloidal_ * chi.gradient + toroidal_ * mesh::cross(unit, psi.gradient),
            toroidal_ * psi.laplacian / (2.0 * radius_m)};
}

void add_sphere_cell(SpherePlate& plate, const Vec3& c, double area, double radius, const MantleSample& mantle) {
    if (!finite3(c) || std::abs(mesh::dot(c, c) - 1.0) > 1e-9 || !(area > 0) || !std::isfinite(area) ||
        !(radius > 0) || !std::isfinite(radius) || !finite3(mantle.velocity_m_s) || !std::isfinite(mantle.spin_rad_s))
        throw std::runtime_error("Invalid cell in plate drag moments");
    const double lever = area * radius * radius, spin = area * area / (2.0 * pi);
    const double u[3] = {c.x, c.y, c.z};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
        plate.drag_m4[static_cast<size_t>(3 * i + j)] += lever * ((i == j ? 1.0 : 0.0) - u[i] * u[j]) + spin * u[i] * u[j];
    plate.mantle_m4_s = plate.mantle_m4_s + (area * radius) * mesh::cross(c, mantle.velocity_m_s)
        + (spin * mantle.spin_rad_s) * c;
    plate.centre_sum_m2 = plate.centre_sum_m2 + area * c;
    plate.area_m2 += area;
}

Vec3 sphere_centroid(const SpherePlate& plate) {
    const double norm = mesh::length(plate.centre_sum_m2);
    return norm > 0 ? (1.0 / norm) * plate.centre_sum_m2 : Vec3{};
}

std::vector<SpherePlate> sphere_plates(const grid::TectonicGrid& grid, const std::vector<uint32_t>& owners,
                                       uint32_t count, const std::vector<MantleSample>& mantle) {
    if (grid.kind() != grid::TectonicGrid::Kind::cube || owners.size() != grid.cellcount() ||
        mantle.size() != grid.cellcount())
        throw std::runtime_error("Plate drag moments need the cube grid, one owner and one mantle sample per cell");
    const double radius = grid.mesh()->radius_m();
    std::vector<SpherePlate> plates(count);
    for (grid::CellId c = 0; c < grid.cellcount(); ++c)
        if (owners[c] < count) add_sphere_cell(plates[owners[c]], grid.position(c), grid.area_m2(c), radius, mantle[c]);
    return plates;
}

SphereContact sphere_contact(const grid::TectonicGrid& grid, const grid::Edge& edge, uint32_t a, uint32_t b) {
    if (grid.kind() != grid::TectonicGrid::Kind::cube || edge.a >= grid.cellcount() || edge.b >= grid.cellcount())
        throw std::runtime_error("A sphere contact needs a cube grid edge");
    SphereContact contact;
    contact.a = a;
    contact.b = b;
    contact.position = mesh::normalise(grid.position(edge.a) + grid.position(edge.b));
    contact.normal = edge.normal;
    contact.length_m = edge.length_m;
    return contact;
}

double sphere_closing_m_s(const SphereContact& contact, const Vec3& wa, const Vec3& wb, double radius) {
    return mesh::dot(mesh::cross(wa - wb, radius * contact.position), contact.normal);
}

SphereSystem sphere_system(const std::vector<SpherePlate>& plates, const std::vector<SphereContact>& contacts,
                           double radius) {
    if (!(radius > 0) || !std::isfinite(radius)) throw std::runtime_error("Invalid planet radius in force balance");
    const size_t n = 3 * plates.size();
    SphereSystem system;
    system.matrix.assign(n * n, 0.0);
    system.rhs.assign(n, 0.0);
    auto& d = system.diagnostics;
    d.enabled = true;
    for (size_t p = 0; p < plates.size(); ++p) {
        const auto& plate = plates[p];
        bool finite = std::isfinite(plate.area_m2) && finite3(plate.mantle_m4_s);
        for (double value : plate.drag_m4) finite = finite && std::isfinite(value);
        if (!(plate.area_m2 > 0.0) || !finite)
            throw std::runtime_error("Invalid plate area or mantle drag in force balance");
        double block[9];
        for (size_t k = 0; k < 9; ++k) block[k] = kBasalDrag * plate.drag_m4[k];
        add3(system.matrix, n, p, p, 1.0, block);
        system.rhs[3 * p] = kBasalDrag * plate.mantle_m4_s.x;
        system.rhs[3 * p + 1] = kBasalDrag * plate.mantle_m4_s.y;
        system.rhs[3 * p + 2] = kBasalDrag * plate.mantle_m4_s.z;
    }
    for (const auto& edge : contacts) {
        if (edge.a >= plates.size() || edge.b >= plates.size() || edge.a == edge.b ||
            !(edge.length_m > 0.0) || !std::isfinite(edge.length_m) ||
            !std::isfinite(edge.closing_m_s) || !std::isfinite(edge.age_a_myr) || !std::isfinite(edge.age_b_myr) ||
            !std::isfinite(edge.mountain_force_n_m) || edge.mountain_force_n_m < 0 ||
            !std::isfinite(edge.activation_speed_m_s) || edge.activation_speed_m_s < 0 ||
            !std::isfinite(edge.slab_pull_scale) || edge.slab_pull_scale < 0 ||
            edge.descending_side < -1 || edge.descending_side > 1 ||
            !std::isfinite(edge.activity_response) || edge.activity_response<0 || edge.activity_response>1 ||
            !std::isfinite(edge.ridge_activity) || edge.ridge_activity < -1 || edge.ridge_activity>1 ||
            !std::isfinite(edge.slab_activity) || edge.slab_activity < -1 || edge.slab_activity>1 ||
            !finite3(edge.normal) || !finite3(edge.position) ||
            std::abs(mesh::dot(edge.normal, edge.normal) - 1.0) > 1e-6 ||
            std::abs(mesh::dot(edge.position, edge.position) - 1.0) > 1e-6 ||
            std::abs(mesh::dot(edge.normal, edge.position)) > 1e-6)
            throw std::runtime_error("Invalid contact in plate force balance");
        const bool continental = edge.continental_a && edge.continental_b;
        const double normal_drag = continental ? kContinentalBoundaryDrag : kOceanBoundaryDrag;
        // A force f along n at the edge exerts the torque f R (x x n) on a.
        const Vec3 x = edge.position, arm = radius * mesh::cross(x, edge.normal);
        const double xs[3] = {x.x, x.y, x.z}, ms[3] = {arm.x, arm.y, arm.z};
        // The planar dashpot of the relative velocity (w_a - w_b) x R x as a symmetric
        // torque block: L (shear R^2 (I - x x^T) + (normal - shear) m m^T), m = R x x n.
        double block[9];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
            block[3 * i + j] = edge.length_m * (kShearBoundaryDrag * radius * radius * ((i == j ? 1.0 : 0.0) - xs[i] * xs[j]) +
                (normal_drag - kShearBoundaryDrag) * ms[i] * ms[j]);
        add3(system.matrix, n, edge.a, edge.a, 1.0, block);
        add3(system.matrix, n, edge.b, edge.b, 1.0, block);
        add3(system.matrix, n, edge.a, edge.b, -1.0, block);
        add3(system.matrix, n, edge.b, edge.a, -1.0, block);
        if (continental) d.continental_contact_length_km += edge.length_m / 1000.0;
        if (edge.mountain_force_n_m > 0) {
            const double force = edge.mountain_force_n_m * edge.length_m;
            add_torque(system.rhs, edge.a, -force, arm);
            add_torque(system.rhs, edge.b, force, arm);
            d.mountain_boundary_force_n += force;
        }
        // Activity and slab choice exactly as in solve.
        const double threshold = edge.activation_speed_m_s > 0 ? 0 : 1e-12;
        const double activation = edge.activation_speed_m_s > 0
            ? std::abs(edge.closing_m_s) / std::hypot(edge.closing_m_s, edge.activation_speed_m_s) : 1;
        const double ridge_target=edge.closing_m_s < -threshold?activation:0;
        const double slab_target=edge.closing_m_s > threshold?activation:0;
        const double ridge=edge.ridge_activity<0?ridge_target:
            edge.ridge_activity+edge.activity_response*(ridge_target-edge.ridge_activity);
        const double slab=edge.slab_activity<0?slab_target:
            edge.slab_activity+edge.activity_response*(slab_target-edge.slab_activity);
        if (ridge>0 && !edge.continental_a && !edge.continental_b) {
            const double force = 2.0e12 * edge.length_m * ridge;
            add_torque(system.rhs, edge.a, -force, arm);
            add_torque(system.rhs, edge.b, force, arm);
            d.ridge_length_km += edge.length_m / 1000.0;
        }
        if (slab>0 && !continental) {
            const bool a_sinks = !edge.continental_a &&
                (edge.continental_b || edge.descending_side == 0 ||
                 (edge.descending_side < 0 && edge.age_a_myr > edge.age_b_myr));
            const bool b_sinks = !edge.continental_b &&
                (edge.continental_a || edge.descending_side == 1 ||
                 (edge.descending_side < 0 && edge.age_b_myr > edge.age_a_myr));
            if (a_sinks || b_sinks) {
                const uint32_t p = a_sinks ? edge.a : edge.b;
                const double force = slab_pull_per_length(a_sinks ? edge.age_a_myr : edge.age_b_myr)
                    * edge.length_m * (a_sinks ? 1.0 : -1.0) * slab * edge.slab_pull_scale;
                add_torque(system.rhs, p, force, arm);
                d.subduction_length_km += edge.length_m / 1000.0;
            }
        }
    }
    return system;
}

SphereSolution solve_sphere(const std::vector<SpherePlate>& plates, const std::vector<SphereContact>& contacts,
                            double radius) {
    SphereSolution result;
    result.rotation_rad_s.resize(plates.size());
    if (plates.empty()) return result;
    auto system = sphere_system(plates, contacts, radius);
    const size_t n = system.rhs.size();
    const auto rotation = solve_spd(system.matrix, system.rhs);
    auto& d = system.diagnostics;
    double residual_sq = 0.0, torque_sq = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double drag = 0.0;
        for (size_t j = 0; j < n; ++j) drag += system.matrix[i * n + j] * rotation[j];
        residual_sq += (drag - system.rhs[i]) * (drag - system.rhs[i]);
        torque_sq += system.rhs[i] * system.rhs[i];
        d.driving_power_w += system.rhs[i] * rotation[i];
        d.dissipation_w += drag * rotation[i];
    }
    d.relative_force_residual = std::sqrt(residual_sq / std::max(1.0, torque_sq));
    for (size_t p = 0; p < plates.size(); ++p) {
        const Vec3 w{rotation[3 * p], rotation[3 * p + 1], rotation[3 * p + 2]};
        result.rotation_rad_s[p] = w;
        const double speed = mesh::length(mesh::cross(w, radius * sphere_centroid(plates[p]))) * cm_yr_per_m_s;
        d.max_speed_cm_yr = std::max(d.max_speed_cm_yr, speed);
        d.mean_speed_cm_yr += speed / plates.size();
    }
    result.diagnostics = d;
    return result;
}

SphereSolution solve_sphere_consistent(const std::vector<SpherePlate>& plates, std::vector<SphereContact> contacts,
                                       std::vector<Vec3> rotation, double radius) {
    if (rotation.size()!=plates.size() || !(radius > 0) || !std::isfinite(radius))
        throw std::runtime_error("Invalid nonlinear motion seed");
    // As solve_consistent, in surface speed units: a rotation component times R in cm/yr.
    const size_t n=3*rotation.size();
    const double scale=radius*cm_yr_per_m_s;
    const auto component=[](Vec3& w,size_t k) -> double& { return k==0?w.x:k==1?w.y:w.z; };
    const auto evaluate=[&](const std::vector<Vec3>& w) {
        for (auto& edge:contacts) edge.closing_m_s=sphere_closing_m_s(edge,w[edge.a],w[edge.b],radius);
        return solve_sphere(plates,contacts,radius);
    };
    const auto residual=[&](const std::vector<Vec3>& w,const SphereSolution& s) {
        std::vector<double> r(n);
        for (size_t p=0;p<w.size();++p) {
            r[3*p]=(s.rotation_rad_s[p].x-w[p].x)*scale;
            r[3*p+1]=(s.rotation_rad_s[p].y-w[p].y)*scale;
            r[3*p+2]=(s.rotation_rad_s[p].z-w[p].z)*scale;
        }
        return r;
    };
    const auto norm=[](const std::vector<double>& r) {
        double sum=0;for (double value:r) sum+=value*value;return std::sqrt(sum);
    };
    for (uint32_t iteration=1;iteration<=96;++iteration) {
        auto result=evaluate(rotation);
        const auto r=residual(rotation,result);
        double error=0;
        for (size_t p=0;p<plates.size();++p)
            error=std::max(error,mesh::length(result.rotation_rad_s[p]-rotation[p])*scale);
        result.diagnostics.nonlinear_iterations=iteration;
        result.diagnostics.velocity_consistency_cm_yr=error;
        if (error<1e-4) return result;
        std::vector<double> jac(n*n),step=r;
        constexpr double epsilon=1e-4; // cm/yr, finite-difference Jacobian.
        for (size_t j=0;j<n;++j) {
            auto perturbed=rotation;
            component(perturbed[j/3],j%3)+=epsilon/scale;
            const auto q=residual(perturbed,evaluate(perturbed));
            for (size_t i=0;i<n;++i) jac[i*n+j]=(q[i]-r[i])/epsilon;
        }
        for (auto& value:step) value=-value;
        bool singular=false;
        for (size_t k=0;k<n;++k) {
            size_t pivot=k;
            for (size_t i=k+1;i<n;++i) if (std::abs(jac[i*n+k])>std::abs(jac[pivot*n+k])) pivot=i;
            if (std::abs(jac[pivot*n+k])<1e-12) {singular=true;break;}
            for (size_t j=k;j<n;++j) std::swap(jac[k*n+j],jac[pivot*n+j]);
            std::swap(step[k],step[pivot]);
            for (size_t i=k+1;i<n;++i) {
                const double factor=jac[i*n+k]/jac[k*n+k];
                for (size_t j=k;j<n;++j) jac[i*n+j]-=factor*jac[k*n+j];
                step[i]-=factor*step[k];
            }
        }
        if (!singular) for (size_t i=n;i-->0;) {
            for (size_t j=i+1;j<n;++j) step[i]-=jac[i*n+j]*step[j];
            step[i]/=jac[i*n+i];
        }
        double best=norm(r);
        auto next=rotation;
        for (int method=0;method<2;++method) {
            if (method==0 && singular) continue;
            for (double damping=1;damping>=1./256;damping*=.5) {
                auto candidate=rotation;
                for (size_t j=0;j<n;++j)
                    component(candidate[j/3],j%3)+=damping*(method?r[j]:step[j])/scale;
                const double candidate_error=norm(residual(candidate,evaluate(candidate)));
                if (candidate_error<best) {best=candidate_error;next=std::move(candidate);break;}
            }
        }
        if (best>=norm(r) || iteration==96)
            throw std::runtime_error("Activity-dependent plate motion did not converge; residual cm/yr="+std::to_string(error));
        rotation=std::move(next);
    }
    throw std::runtime_error("Unreachable nonlinear motion state");
}

double inland_pressure(const grid::TectonicGrid& grid, const std::vector<uint32_t>& owners,
        const std::vector<double>& pressure, grid::CellId cell, std::uint8_t direction, double range_m) {
    if (grid.kind() != grid::TectonicGrid::Kind::cube || owners.size() != grid.cellcount() ||
        pressure.size() != owners.size() || cell >= owners.size() || direction >= grid::directioncount ||
        !std::isfinite(range_m) || range_m < 0)
        throw std::runtime_error("Invalid inland pressure strip");
    if (range_m==0) return pressure[cell];
    const double radius = grid.mesh()->radius_m();
    const auto cdf=[](double u) { u=std::min(1.,u);return u*u*(3-2*u); };
    // A straight walk closes after 4 N cells; never count a strip twice.
    const size_t limit = 4 * static_cast<size_t>(grid.facecells());
    const double reach = std::min(range_m, 2.0 * pi * radius);
    const auto owner = owners[cell];
    double sum=0,weight=0,travelled=0;
    for (size_t step=0;step<limit && travelled<reach;++step) {
        if (owners[cell]!=owner) break;
        const grid::Step next = grid.ahead(cell, direction);
        const Vec3 p = grid.position(cell), q = grid.position(next.cell);
        const double width = radius * std::atan2(mesh::length(mesh::cross(p, q)), mesh::dot(p, q));
        const double k=cdf((travelled+width)/reach)-cdf(travelled/reach);
        sum+=k*pressure[cell];weight+=k;
        travelled += width;
        cell = next.cell;
        direction = next.direction;
    }
    return sum/weight;
}
} // namespace tectonics::dynamics
