#pragma once

// Index helpers shared by the regional pass's translation units (erosion plan Phase R2,
// step R6): not part of the pass's interface. regional_surface.cpp keeps its own copies
// (its expressions are what the overlap tests established); the newer units use these.

#include <cstddef>

namespace landscape::regionaldetail
{
inline std::size_t at(int x, int y, int width)
{
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x);
}

// Columns wrap around the planet; rows beyond the poles do not exist.
inline int wrapcolumn(int column, int width) { return ((column % width) + width) % width; }

// Whether `test` holds for a cell of the 3 x 3 neighbourhood of global cell (i, j), the
// cell included (columns wrap, rows beyond the poles are skipped).
template<typename Test>
bool aroundcell(int width, int height, int i, int j, Test test)
{
    for (int dj = -1; dj <= 1; dj++)
    {
        const int row = j + dj;
        if (row < 0 || row >= height)
            continue;
        for (int di = -1; di <= 1; di++)
            if (test(at(wrapcolumn(i + di, width), row, width)))
                return true;
    }
    return false;
}
}
