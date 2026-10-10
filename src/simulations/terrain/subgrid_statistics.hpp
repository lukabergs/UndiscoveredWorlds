#pragma once

// terrain.subgrid_statistics (diagnoser, docs/architecture/multiscale-plan.md, step S3.6):
// the sub-grid statistics the regional levels reproduce, geosphere.relief_amplitude (m) and
// geosphere.valley_spacing (km), from the erosion run's stream-power and hillslope laws
// (landscape.stream_power.*, landscape.hillslope.*), its erodibility and erosion rate, the
// routing's runoff and the final base elevation; the formulas are the kernel's
// (simulations/landscape/subgrid_statistics.hpp). Runs after the last step that edits the
// base elevation.

class planet;

void applysubgridstatisticsstage(planet& world);
