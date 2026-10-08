#pragma once

// The region command of world_pipeline (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2d): the regional pass on a window of a world bundle, written as a
// regional bundle (io/regional_bundle.hpp).

#include <ostream>
#include <string>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"

void printregionusage(std::ostream& out);

// The arguments after the command word; returns the exit code.
int runregioncommand(const std::vector<std::string>& arguments);

// The window of global cells x0..x1 (columns from the west edge at -180 degrees; x0 > x1
// crosses the date line) and rows y0..y1 (from the north), both inclusive: the cells'
// edge degrees through landscape::regionalwindow, checked to give exactly those cells.
// Throws std::invalid_argument on cells off the grid, y0 > y1 or a domain wider than
// the planet.
landscape::RegionalWindow regionalwindowfromcells(int worldwidth, int scale, int halo, int x0, int y0, int x1, int y1);
