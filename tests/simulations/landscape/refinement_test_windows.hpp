#pragma once

// Level windows built by hand for the refinement tests (docs/architecture/multiscale-plan.md,
// step S0.4), so that refinement_grid_tests.cpp and refinement_projection_tests.cpp need not
// wait for levelwindow and halochain (S2.1): the domain follows from the kept box by the
// rule of LevelWindow (refinement_grid.hpp). The caller passes an aligned kept box.

#include "simulations/landscape/refinement_grid.hpp"

#include <algorithm>
#include <cstdint>

inline landscape::LevelWindow handlevelwindow(int worldwidth, int factor, int level, int halo, std::int64_t column0,
                                              int columns, std::int64_t row0, int rows)
{
    landscape::LevelWindow window;
    window.worldwidth = worldwidth;
    window.worldheight = worldwidth / 2;
    window.factor = factor;
    window.level = level;
    window.halo = halo;
    window.column0 = column0;
    window.columns = columns;
    window.row0 = row0;
    window.rows = rows;
    const int cells = window.halocells();
    window.domaincolumn0 = landscape::wraplevel(column0 - cells, window.levelwidth());
    window.domaincolumns = columns + 2 * cells;
    window.domainrow0 = std::max<std::int64_t>(0, row0 - cells);
    window.domainrows = static_cast<int>(std::min(window.levelheight(), row0 + rows + cells) - window.domainrow0);
    return window;
}
