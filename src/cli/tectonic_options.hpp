#pragma once

// Tectonic scenario options of world_pipeline: the engine scenario starts from
// configs/tectonics-baseline.json and each --key value overrides one member.

#include <ostream>
#include <string>

#include "tectonic_recipe.hpp"

// False when the argument is not a scenario option; throws std::runtime_error
// on a value the option cannot take.
bool parsetectonicoption(tectonics::scenario::Scenario& scenario, const std::string& argument, const std::string& value);
void printtectonicoptionsusage(std::ostream& out);
