#ifndef TOPOGRAPHY_TIFF_HPP
#define TOPOGRAPHY_TIFF_HPP

#include <cstdint>
#include <vector>

namespace platec::io {

bool writeGrayscaleTiff16(const char* filename, uint32_t width, uint32_t height,
                          const uint16_t* samples);

bool readGrayscaleTiff16(const char* filename, std::vector<uint16_t>& samples,
                         uint32_t& width, uint32_t& height);

} // namespace platec::io

#endif
