#pragma once
#include <cmath>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

namespace climateio
{
// Explicit experimental boundary input, independent of saved-world formats.
// Cell-centred north-up rows from -180 to +180; dimensionless drag contrast.
inline std::vector<float> readSurfaceDragContrast(std::istream& stream, int columns, int rows)
{
    std::string magic;
    int inputColumns = 0, inputRows = 0;
    if (columns < 2 || columns % 2 || rows != columns / 2 ||
        !(stream >> magic >> inputColumns >> inputRows) || magic != "UWROUGH1" ||
        inputColumns != columns || inputRows != rows)
        throw std::runtime_error("Surface drag input has invalid grid or registration header");
    std::vector<float> values(static_cast<std::size_t>(columns) * rows);
    for (float& value : values)
        if (!(stream >> value) || !std::isfinite(value) || value < 0.25f || value > 4.0f)
            throw std::runtime_error("Surface drag input requires finite positive contrast in [0.25,4]");
    std::string extra;
    if (stream >> extra) throw std::runtime_error("Surface drag input contains extra cells");
    return values;
}
}
