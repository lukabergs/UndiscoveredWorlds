#pragma once

// Analytical relief of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2c, and section 6, "The solve"): the transient solution of the
// linear stream power law (n = 1, uplift 0) on a fine patch whose pinned cells keep
// their elevation (Steer 2021, Earth Surface Dynamics 9, 1239; Tzathas, Gailleton,
// Steer and Cordonnier 2024: a non-flat starting terrain, fixed boundary heights, a
// prescribed network, elapsed time and a threshold-slope hillslope). With n = 1 and
// no uplift the solution travels upstream along the flow paths: the response time
// of a cell is tau(x) = sum over the links from x down to its first pinned cell of
// L / (K Q^m); after an elapsed time t a cell takes the start elevation of the point
// at response time tau(x) - t downstream of it (interpolated along the link), or
// the pinned cell's elevation when tau(x) <= t. Then every slope steeper than the
// threshold slope is lowered to it, from the pinned cells upwards. Receivers come
// from the result, so routing and solution are repeated until no receiver changes
// (or maxiterations). Units: elevation m, K in m/Myr per (km^3/yr)^m per unit slope
// (landscape.stream_power.erodibility), discharge m^3/yr inside. Serial in a
// canonical order (each step over the cells in flow order with ties by index), so
// the result is bit-identical for any thread count and, for a cell whose catchment
// and downstream path lie inside both of two domains, the same in both. Planet-free.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
struct AnalyticalReliefParameters
{
    double erodibility = 750.0;     // K, as the global kernel's
    double areaexponent = 0.5;      // m
    double elapsed_myr = 1.0;       // t
    double thresholdslope = 0.6;    // hillslope limit (rise over run)
    int maxiterations = 20;         // routing and solution repeated at most this often
};

struct AnalyticalReliefInput
{
    std::vector<float> start_m;          // starting surface (the upsampled terrain with noise)
    std::vector<uint8_t> pinned;         // 1: elevation fixed (channels, sea, the patch's outer ring)
    std::vector<float> runoff_m_yr;      // per cell; discharge = sum of runoff x area over the catchment
};

struct AnalyticalReliefResult
{
    std::vector<float> elevation_m;
    std::vector<uint8_t> direction;      // 1..8 towards the receiver, 0 on pinned cells
    std::vector<double> discharge_m3_yr; // accumulated over each unpinned catchment (0 on pinned cells)
    int iterations = 0;                  // routing passes used
    bool converged = false;              // no receiver changed in the last pass
};

// Throws std::invalid_argument when the input sizes differ from the metric's or a
// parameter is out of range (K <= 0, m < 0, t < 0, threshold slope <= 0, no
// iteration). Receivers: steepest descent where a lower neighbour exists, otherwise
// the neighbour from which a Priority-Flood seeded at the pinned cells reached the
// cell (Barnes' directions, as depressions.cpp gives them), so every unpinned cell
// has a path to a pinned cell and no receiver loop exists; the surface is not filled.
AnalyticalReliefResult solveanalyticalrelief(const PatchMetric& metric, const AnalyticalReliefInput& input,
                                             const AnalyticalReliefParameters& parameters);
}
