#pragma once

// Planet adapter of the regional pass (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E3): the global fields the planet-free pass reads, through the
// float getters only, and its parameters from tuning::landscape::regional and
// tuning::landscape::stream_power.

#include "simulations/landscape/regional_pass.hpp"

class planet;

landscape::RegionalGlobalInputs gatherregionalinputs(const planet& world);
landscape::RegionalParameters regionalparameters(std::uint64_t seed);
