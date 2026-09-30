#ifndef CLIMATE_TIFF_HPP
#define CLIMATE_TIFF_HPP

#include <cstdint>

namespace climateio
{
bool writefloat32geotiff(const char* filename, std::uint32_t width,
    std::uint32_t height, const float* samples, float nodatavalue = -9999.9f);
}

#endif
