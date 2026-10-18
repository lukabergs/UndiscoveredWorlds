#pragma once

// Sub-grid statistics of the global grid (docs/architecture/multiscale-plan.md, step S3.6):
// the statistics of the unresolved terrain that every regional level reproduces, after the
// subgrid orography of weather models (Baines and Palmer 1990; Lott and Miller 1997), from
// the erosion run's own laws. Planet-free; the process (terrain/subgrid_statistics.cpp)
// passes landscape.stream_power.*, landscape.hillslope.* and terrain.subgrid.*.
//
// Valley spacing (Perron, Kirchner and Dietrich 2009): with dz/dt = D lap(z) - K_A A^m |grad z|
// + U (n = 1), the Peclet number K_A L^(2m+1) / D is 1 at Lc = (D / K_A)^(1/(2m+1)), and the
// first-order valley spacing is lambda = c Lc (c from 6.4 to 12.7; 7.0 to 9.5 at the field
// sites), independent of the erosion rate. The erosion's D and K are grid-scale values (Lc of
// order 100 km), so D/K is calibrated to a field site of the paper's Table 1 instead: D/K_ref
// (m^(2m_ref+1)) at a reference runoff and erodibility, and each cell's D/K is that times
// K_A(reference) / K_A(cell). With discharge Q = runoff x area in the site's law K Q^m_ref,
// K_A = K x relative erodibility x runoff^m_ref, so the ratio is (e_ref / e) x (r_ref / r)^m_ref
// and Lc = (D/K)^(1/(2 m_ref + 1)) varies as r^(-m_ref / (2 m_ref + 1)). A spacing wider than the
// cell means no sub-grid valleys: lambda is capped at the cell size (the square root of its
// area), so a stored value equal to the cell size reads "no sub-grid valleys resolved", and
// consumers compare it with the cell size.
//
// Relief amplitude: a grid-scale relief, the steady-state relief of the erosion run's own
// (grid-scale) diffusivity D over half the cell: a slope of half-length l eroding at the rate E
// everywhere is the parabola z = E (l^2 - x^2) / (2 D), of relief E l^2 / (2 D), with l half the
// cell size. It is held to the cell's local relief (the largest absolute elevation difference
// to its 8 neighbours, the expression of the regional pass's localrelief). The erosion's critical
// slope is a diffusion threshold, not a maximum slope, and bounds nothing here.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct SubgridParameters
{
    double diffusivity_km2_myr = 0.0;    // the relief's grid-scale D (landscape.hillslope.diffusivity_km2_myr)
    double spacingcoefficient = 8.0;     // c: lambda = c Lc (terrain.subgrid.spacing_coefficient)
    double referencedoverk = 0.0;        // D/K at the reference, m^(2 m_ref + 1) (terrain.subgrid.reference_d_over_k)
    double referencem = 0.35;            // m_ref of that D/K and of the runoff scaling (terrain.subgrid.reference_m)
    double referencerunoff_mm_yr = 0.0;  // runoff of the reference (terrain.subgrid.reference_runoff_mm_yr)
    double referenceerodibility = 1.0;   // relative erodibility of the reference (terrain.subgrid.reference_erodibility)
};

// Lc in km for a cell of relative erodibility and runoff (mm/yr):
// (D/K_ref x (e_ref / e) x (r_ref / r)^m_ref)^(1/(2 m_ref + 1)) m / 1000; 0 when the cell's or the
// reference's erodibility or runoff, or D/K_ref, is not positive (no channel: no length).
double characteristiclength_km(const SubgridParameters& parameters, double relativeerodibility, double runoff_mm_yr);

// The steady-state hillslope relief (m) of half-length l (km) at the erosion rate E (m/Myr):
// E l^2 / (2 D); 0 with D not positive.
double hillsloperelief_m(const SubgridParameters& parameters, double erosionrate_m_myr, double halflength_km);

// Both statistics on every cell of the metric's grid (in its cell order). `sea` 1
// marks the sea, where both are 0. `erodibility` is relative (geosphere.erodibility),
// `runoff_mm_yr` hydrosphere.runoff, `erosionrate_m_myr` geosphere.erosion_rate (negative
// counts as 0), `elevation_m` the base elevation. A spacing beyond the cell size is stored as
// the cell size (no sub-grid valleys resolved); returns the number of land cells so capped.
// Throws std::invalid_argument on a size mismatch.
std::size_t subgridstatistics(const SphereMetric& metric, const SubgridParameters& parameters,
                              const std::vector<float>& elevation_m, const std::vector<std::uint8_t>& sea,
                              const std::vector<float>& erodibility, const std::vector<float>& runoff_mm_yr,
                              const std::vector<float>& erosionrate_m_myr, std::vector<float>& relief_m,
                              std::vector<float>& spacing_km);
}
