#include "tectonic_recipe.hpp"

namespace tectonics::scenario {
Scenario baseline() {
    Scenario scenario;
#include "tectonic_baseline.inc"
    return scenario;
}
}
