#pragma once
#include <memory>
#include <vector>

#include "world/world_grid.hpp"
class planet;

// Overrides of configs/tectonics-baseline.json; zero selects its value.
struct WorldTectonicsOptions {
    double durationMyr = 0;
    int plateCount = 0;
    double deltaTimeMyr = 0;
    bool mineralHistory = true;
};

// Tectonics engine and mineral history carried from the Geosphere initializer to
// the tectonics evolver. Solver state, not world state: copies never carry it.
struct GeologySolverState;

// Geosphere initializer: creates the engine (which seeds the starting topography)
// and the mineral history's inherited provinces, clears tectonic provenance and
// writes the starting elevation. Expects a cleared world.
std::shared_ptr<GeologySolverState> initializegeosphere(planet& world, worldgrid<bool>& shelves,
                                                        const WorldTectonicsOptions& options = {});

// Tectonics evolver: runs the engine to the end and writes the final Geosphere.
// State made for another seed or other options, or no state, is re-created first;
// creation is deterministic, so the result does not depend on where it happened.
void evolvegeosphere(planet& world, worldgrid<bool>& shelves, const WorldTectonicsOptions& options,
                     std::shared_ptr<GeologySolverState> state);

// Both steps.
void generateworldtectonics(planet& world, worldgrid<bool>& shelves,
                            const WorldTectonicsOptions& options = {});
