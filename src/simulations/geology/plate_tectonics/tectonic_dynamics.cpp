#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"
#include "simulations/geology/plate_tectonics/tectonic_thermal.hpp"

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
        const std::vector<double>& pressure,uint32_t cell,uint32_t axis,int direction,double range_m) {
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
} // namespace tectonics::dynamics
