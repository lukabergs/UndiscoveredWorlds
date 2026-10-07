#pragma once

// terrain.fluvial_erosion (evolver): post-tectonic denudation of the base
// elevation with the landscape stream power kernel over
// landscape.erosion.duration_myr in steps of step_myr; routes at every step
// and writes the final routing fields itself (docs/architecture/erosion-hydrology-plan.md,
// section 3, Processes). With landscape.glacial.enabled the kernel also erodes with the
// ice of an ELA from the coarse climate, routing on the bed with the ice surface's sea, and
// the stage writes the glacial fields and the present-day ice (erosion plan Phase G).

#include <functional>
#include <vector>

class planet;

// Checkpoints inside the run, in Myr of erosion time. After the step that reaches
// or passes a requested time (one step may serve several), the world holds the
// state as of that step, as the end of the run writes it (base elevation, the
// erosion fields with rates over the elapsed time, erosion_duration = the elapsed
// time, the volumes so far and the routing of the terrain then), and the callback
// runs with the elapsed time. The end of the run is not a checkpoint; checkpoints
// leave the result unchanged.
struct FluvialErosionCheckpoints
{
    double everyMyr = 0;
    std::vector<double> atMyr;
    std::function<void(double elapsedMyr)> checkpoint;
};

void applyfluvialerosionstage(planet& world, const FluvialErosionCheckpoints& checkpoints = {});
