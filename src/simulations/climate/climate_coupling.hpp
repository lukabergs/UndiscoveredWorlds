#pragma once

#include <vector>

class climateworld;

// World-grid adapters and retained iteration state for the numerical atmosphere,
// ocean and moisture solvers. Order: pressure -> winds -> ocean -> evaporation
// -> rainfall; convergence revisits these stages with lagged diagnosed heating.
void createpressuremap(climateworld& world);
void updatehorsebeltsfrompressure(climateworld& world);
void createvectorwindmap(climateworld& world);
void createoceancurrentmap(climateworld& world);
void createsurfacetemperaturemap(climateworld& world);
void createadvectedrainfall(climateworld& world, std::vector<std::vector<int>>& inland, std::vector<std::vector<int>>& fractal);
void convergeclimatecoupling(climateworld& world, std::vector<std::vector<int>>& fractal);
