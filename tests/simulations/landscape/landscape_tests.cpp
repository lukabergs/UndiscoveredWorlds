// Suite Landscape: one fixture per kernel of src/simulations/landscape, on
// synthetic grids of a few hundred cells with exact expected numbers
// (docs/architecture/erosion-hydrology-plan.md, section 3, Tests). Fixtures
// SphereMetric, FlowGraph and Depressions are step A1; StreamPower is step A2.

#include "simulations/landscape/depressions.hpp"
#include "simulations/landscape/flow_graph.hpp"
#include "simulations/landscape/sphere_metric.hpp"
#include "simulations/landscape/stream_power.hpp"

#include <gtest/gtest.h>
