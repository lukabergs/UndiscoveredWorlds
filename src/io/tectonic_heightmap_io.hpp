#ifndef HEIGHTMAP_IO_HPP
#define HEIGHTMAP_IO_HPP

#include <stdio.h>
#include <png.h>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "topography_codec.hpp"

// This function actually writes out the PNG image file. The string 'title' is
// also written into the image file
int writeImageGray(const char* filename, int width, int height, float *heightmap, const char* title);

int writeImageColors(const char* filename, int width, int height, float *heightmap, const char* title);

int writeImageRgb(const char* filename, int width, int height, const png_byte* rgb, const char* title);
int writeImageRgba(const char* filename, int width, int height, const png_byte* rgba, const char* title);

void renderImageGrayRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb);

void renderImageColorsRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb);
void renderMaterialMapRgb(int width, int height, const uint8_t* material_map,
                          std::vector<png_byte>& rgb);
void renderCrustClassMapRgb(int width, int height, const uint8_t* crust_class_map,
                            std::vector<png_byte>& rgb);
void renderCrustAgeMapRgb(int width, int height, const float* crust_age_myr_map,
                          std::vector<png_byte>& rgb);
void renderPlateIdMapRgb(int width, int height, const uint32_t* plate_map,
                         std::vector<png_byte>& rgb);

int writeImageGrayTiff16(const char* filename, int width, int height, const uint16_t* heightmap);
int readImageGrayTiff16(const char* filename, std::vector<uint16_t>& heightmap, int& width,
                        int& height);
int writeTopographyMetadataJson(const char* filename, const TopographyCodec::Metadata& metadata);
int readTopographyMetadataJson(const char* filename, TopographyCodec::Metadata& metadata);

#endif
