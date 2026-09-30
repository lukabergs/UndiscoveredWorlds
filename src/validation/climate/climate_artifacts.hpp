#pragma once
class planet;
class climatebenchmarkmapselection;

// Procedural worlds share field/map exporters without entering the Earth
// benchmark registry, workbook, ranking or observed-data calibration workflow.
bool exportproceduralclimateartifacts(planet& world, const climatebenchmarkmapselection& maps);
