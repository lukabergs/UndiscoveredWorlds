#pragma once

class climateworld;

// Writes the solver diagnostics of the last climate run (budgets, coupling,
// precipitation statistics, circulation precision) as CSV tables; does nothing
// outside a recorded generation run.
void writeclimatediagnostics(climateworld& world);
