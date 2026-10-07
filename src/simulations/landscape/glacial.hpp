#pragma once

// Glacial erosion in stream-power form after Hergarten 2021 (Earth Surface Dynamics
// 9, 937; docs/architecture/erosion-hydrology-plan.md, Phase G and section 6,
// "Decided for Phase G"): an equilibrium-line altitude per cell from the climate,
// an ice flux accumulated over the flow graph from an ice source that is positive
// above the ELA and negative below it, and a thickness from that flux. The erosion
// itself is the glacial term of stepstreampower (stream_power.hpp, GlacialInput).
// Units: elevation m, temperature C, precipitation mm/yr, runoff m/yr, ice flux
// m^3/yr of water equivalent. Serial over `order`, so bit-identical for any thread count.

#include <vector>

#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"

namespace landscape
{
struct GlacialParameters
{
    double erodibilityratio = 2.0;       // K_g / K (Hergarten's examples: 1 and 2)
    double cooling_c = 4.0;              // mean cooling of the run's climate below the climate passed in
    double fullicerange_m = 500.0;       // H_f - H_e: above ELA + this the whole runoff becomes ice
    double widthreference_m = 400.0;     // glacier width at an equivalent ice area of 1 km^2
    double widthexponent = 0.3;          // alpha: width ~ A_i^alpha
    double thicknessratio = 0.25;        // h / width
    double maximumthickness_m = 3000.0;  // cap of h (ice sheets)
    double flotationratio = 0.892;       // ice / sea water density: the deepest bed below the sea is this x h
};

// The summer (warmest season) air temperature at the ELA of a glacier receiving
// `precipitation_mm_yr`, from the relation of Ohmura, Kasser and Funk 1992 (Journal
// of Glaciology 38, 397) solved for the temperature on its increasing branch.
double equilibriumlinetemperature(double precipitation_mm_yr);

// ELA per cell: elevation + (warmest - cooling_c - equilibriumlinetemperature(P)) / lapse,
// `warmest` the warmest season's temperature at `elevation` (C), `lapse_c_per_m`
// the lapse rate that climate applied (C per m, positive). Every vector has
// metric.size() elements; `ela_m` is resized and assigned on every cell.
void equilibriumline(const SphereMetric& metric, const std::vector<float>& elevation,
                     const std::vector<float>& warmest_c, const std::vector<float>& precipitation_mm_yr,
                     double lapse_c_per_m, double cooling_c, std::vector<float>& ela_m);

// Ice flux on `graph` (built on `surface`, the bed plus the ice of the previous
// step): each land cell's source is runoff x area x min((surface - ela) /
// fullicerange_m, 1), negative below the ELA; over `order` reversed a cell's flux
// is max(0, source + the fluxes of its donors), so ablation ends a glacier where
// the flux reaches zero (its terminus). Sea cells carry none (ice reaching an
// outlet calves). With `runoff_m_yr` the routing's runoff, the flux never exceeds
// graph.discharge_m3_yr without transmission losses. `iceflux_m3_yr` is resized
// and assigned on every cell.
void accumulateice(const SphereMetric& metric, const FlowGraph& graph, const std::vector<float>& surface,
                   const std::vector<float>& ela_m, const std::vector<float>& runoff_m_yr, double fullicerange_m,
                   std::vector<double>& iceflux_m3_yr);

// Thickness per cell: min(maximumthickness_m, thicknessratio x widthreference_m x
// (A_i / 1 km^2)^widthexponent) with A_i = iceflux / (1 m/yr); 0 where the flux is not positive.
// `thickness_m` is resized and assigned on every cell.
void icethickness(const std::vector<double>& iceflux_m3_yr, const GlacialParameters& parameters,
                  std::vector<float>& thickness_m);
}
