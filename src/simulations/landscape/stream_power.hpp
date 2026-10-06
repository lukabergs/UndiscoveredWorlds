#pragma once

// One time step of the stream power law on a FlowGraph, from the outlets upstream
// (docs/architecture/erosion-hydrology-plan.md, section 3, Kernels):
// - detachment-limited part (n = 1, Braun & Willett 2013): with
//   F_i = K_i Q_i^m dt / L_i, z_i <- (z_i + U_i dt + F_i z_r) / (1 + F_i), z_r the
//   receiver's already updated elevation;
// - deposition: the Davy-Lague term of Yuan et al. 2019 with coefficient G, solved by
//   Gauss-Seidel sweeps over the same order until the largest change is below the
//   tolerance (G = 0 reproduces the detachment-limited scheme exactly);
// - hillslope: linear diffusion with coefficient kappa, explicit, with the per-face
//   flux limiter of the former climate_erosion.cpp (a face moves at most a fixed
//   fraction of the drop times the smaller area);
// - sea cells: sediment that reaches an outlet is routed down the seabed's steepest
//   descent (the graph's sea receivers) and deposited with an e-folding length, the
//   remainder at a seabed pit, so eroded = land deposited + marine deposited.
// Units: elevation m, time Myr, discharge m^3/yr in the graph and km^3/yr in the
// law, erodibility m/Myr per (km^3/yr)^m, uplift m/Myr, kappa km^2/Myr. Storage is
// float, accumulation double. Serial over `order`; per-cell loops may be parallel
// with fixed-order reductions; the result is bit-identical for any thread count.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct StreamPowerParameters
{
    double step_myr = 1.0;                      // dt
    double areaexponent = 0.5;                  // m
    double slopeexponent = 1.0;                 // n; Phase A accepts 1 only (the implicit step is linear in the slope)
    double depositioncoefficient = 1.0;         // G of Yuan et al. 2019; 0 is detachment-limited
    double depositiontolerance_m = 1.0e-3;      // Gauss-Seidel stop: largest elevation change of a sweep
    int depositionmaxsweeps = 50;
    double hillslopediffusivity_km2_myr = 200.0;  // kappa
    double criticalslope = 0.0;                 // 0 off; else diffusion acts on max(0, gradient - criticalslope)
    double sedimenttravel_km = 300.0;           // marine e-folding length
};

struct StreamPowerBudget
{
    double eroded_km3 = 0.0;           // rock removed from land cells (fluvial incision and hillslope lowering)
    double landdeposited_km3 = 0.0;    // sediment laid down on land cells (fluvial deposition and hillslope raising)
    double marinedeposited_km3 = 0.0;  // sediment laid down on sea cells
    int depositionsweeps = 0;          // Gauss-Seidel sweeps the deposition solve used (0 when G = 0)
};

// Advances `elevation` by one step. `graph` was built and accumulated on the
// elevation passed in; `erodibility` (K) and `uplift_m_myr` (U) are per cell and
// read on land cells only. `erosion_m` and `deposition_m` are assigned (not
// accumulated) the thickness removed from and added to each cell by this step,
// both non-negative; the caller sums them over the run. Every vector has
// metric.size() elements. Throws std::invalid_argument on a size mismatch or
// a slope exponent other than 1.
StreamPowerBudget stepstreampower(const SphereMetric& metric, const FlowGraph& graph,
                                  const std::vector<float>& erodibility, const std::vector<float>& uplift_m_myr,
                                  const StreamPowerParameters& parameters, std::vector<float>& elevation,
                                  std::vector<float>& erosion_m, std::vector<float>& deposition_m);
}
