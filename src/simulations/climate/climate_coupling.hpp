#pragma once

#include <vector>

class planet;

// World-grid adapters and retained iteration state for the numerical atmosphere,
// ocean and moisture solvers. Order: pressure -> winds -> ocean -> evaporation
// -> rainfall; convergence revisits these stages with lagged diagnosed heating.
void createpressuremap(planet& world);
void updatehorsebeltsfrompressure(planet& world);
void createvectorwindmap(planet& world);
void createoceancurrentmap(planet& world);
void createsurfacetemperaturemap(planet& world);
void createadvectedrainfall(planet& world, std::vector<std::vector<int>>& inland, std::vector<std::vector<int>>& fractal);
void convergeclimatecoupling(planet& world, std::vector<std::vector<int>>& fractal);
