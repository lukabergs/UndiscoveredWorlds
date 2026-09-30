#pragma once
#include "tectonic_scenario.hpp"

namespace tectonics::scenario {
// Compiled from configs/tectonics-baseline.json. Every production entry point
// starts here, then applies explicit run overrides. Rebuilt when the recipe changes.
Scenario baseline();
}
