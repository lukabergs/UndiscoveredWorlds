#pragma once
class climateworld;
class climatebenchmarkmapselection;

// Procedural worlds share field/map exporters without entering the Earth
// benchmark registry, workbook, ranking or observed-data calibration workflow.
bool exportproceduralclimateartifacts(climateworld& world, const climatebenchmarkmapselection& maps);
