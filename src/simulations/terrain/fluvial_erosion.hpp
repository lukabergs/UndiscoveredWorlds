#pragma once

// terrain.fluvial_erosion (evolver): post-tectonic denudation of the base
// elevation with the landscape stream power kernel over
// landscape.erosion.duration_myr in steps of step_myr; routes at every step
// and writes the final routing fields itself (docs/architecture/erosion-hydrology-plan.md,
// section 3, Processes).

class planet;

void applyfluvialerosionstage(planet& world);
