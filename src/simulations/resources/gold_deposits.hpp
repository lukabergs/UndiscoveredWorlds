#pragma once

class planet;

namespace physical_layers
{
// Evaluates the deposits recorded during tectonics on the final world: present
// surface conditions at each deposit, and river transport of recently
// liberated gold into placers. Placers are rebuilt on every call.
void generate_gold_deposits(planet& world);
}
