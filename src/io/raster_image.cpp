#include "raster_image.hpp"

#include <cstdio>
#include <iostream>
#include <memory>

#include <png.h>

namespace rasterimage
{
namespace
{
struct FileCloser
{
    void operator()(std::FILE* file) const { if (file) std::fclose(file); }
};
using FileHandle = std::unique_ptr<std::FILE, FileCloser>;

FileHandle openfile(const std::string& path, const char* mode)
{
    std::FILE* file = nullptr;
    if (fopen_s(&file, path.c_str(), mode) != 0)
        return FileHandle(nullptr);
    return FileHandle(file);
}
}

void Image::create(unsigned int width, unsigned int height, const Color& color)
{
    width_ = width;
    height_ = height;
    pixels_.resize(static_cast<std::size_t>(width) * height * 4);
    for (std::size_t index = 0; index < pixels_.size(); index += 4)
    {
        pixels_[index] = color.r;
        pixels_[index + 1] = color.g;
        pixels_[index + 2] = color.b;
        pixels_[index + 3] = color.a;
    }
}

Color Image::getPixel(unsigned int x, unsigned int y) const
{
    const Uint8* pixel = &pixels_[(static_cast<std::size_t>(y) * width_ + x) * 4];
    return Color(pixel[0], pixel[1], pixel[2], pixel[3]);
}

void Image::setPixel(unsigned int x, unsigned int y, const Color& color)
{
    Uint8* pixel = &pixels_[(static_cast<std::size_t>(y) * width_ + x) * 4];
    pixel[0] = color.r;
    pixel[1] = color.g;
    pixel[2] = color.b;
    pixel[3] = color.a;
}

bool Image::loadFromFile(const std::string& path)
{
    FileHandle file = openfile(path, "rb");
    if (!file)
    {
        std::cerr << "Failed to load image \"" << path << "\". Reason: Unable to open file\n";
        return false;
    }

    unsigned char signature[8] = {};
    if (std::fread(signature, 1, 8, file.get()) != 8 || png_sig_cmp(signature, 0, 8) != 0)
    {
        std::cerr << "Failed to load image \"" << path << "\". Reason: Not a PNG file\n";
        return false;
    }

    // Declared before setjmp so a libpng error still releases them.
    std::vector<Uint8> pixels;
    std::vector<png_bytep> rows;
    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info = png ? png_create_info_struct(png) : nullptr;
    if (!png || !info)
    {
        png_destroy_read_struct(&png, &info, nullptr);
        return false;
    }
    if (setjmp(png_jmpbuf(png)))
    {
        png_destroy_read_struct(&png, &info, nullptr);
        std::cerr << "Failed to load image \"" << path << "\". Reason: Corrupt PNG\n";
        return false;
    }

    png_init_io(png, file.get());
    png_set_sig_bytes(png, 8);
    png_read_info(png, info);

    const png_byte colortype = png_get_color_type(png, info);
    const png_byte bitdepth = png_get_bit_depth(png, info);
    if (bitdepth == 16)
        png_set_strip_16(png);
    if (colortype == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);
    if (colortype == PNG_COLOR_TYPE_GRAY && bitdepth < 8)
        png_set_expand_gray_1_2_4_to_8(png);
    if (png_get_valid(png, info, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png);
    if (colortype == PNG_COLOR_TYPE_GRAY || colortype == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);
    if (colortype == PNG_COLOR_TYPE_RGB || colortype == PNG_COLOR_TYPE_GRAY || colortype == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png, 0xFF, PNG_FILLER_AFTER);
    png_set_interlace_handling(png);
    png_read_update_info(png, info);

    const png_uint_32 width = png_get_image_width(png, info);
    const png_uint_32 height = png_get_image_height(png, info);
    if (png_get_rowbytes(png, info) != static_cast<std::size_t>(width) * 4)
    {
        png_destroy_read_struct(&png, &info, nullptr);
        std::cerr << "Failed to load image \"" << path << "\". Reason: Unsupported PNG layout\n";
        return false;
    }

    pixels.resize(static_cast<std::size_t>(width) * height * 4);
    rows.resize(height);
    for (png_uint_32 y = 0; y < height; y++)
        rows[y] = pixels.data() + static_cast<std::size_t>(y) * width * 4;
    png_read_image(png, rows.data());
    png_read_end(png, nullptr);
    png_destroy_read_struct(&png, &info, nullptr);

    width_ = width;
    height_ = height;
    pixels_ = std::move(pixels);
    return true;
}

bool Image::saveToFile(const std::string& path) const
{
    if (width_ == 0 || height_ == 0)
        return false;
    FileHandle file = openfile(path, "wb");
    if (!file)
    {
        std::cerr << "Failed to save image \"" << path << "\"\n";
        return false;
    }

    std::vector<png_bytep> rows(height_);
    for (unsigned int y = 0; y < height_; y++)
        rows[y] = const_cast<png_bytep>(pixels_.data() + static_cast<std::size_t>(y) * width_ * 4);
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info = png ? png_create_info_struct(png) : nullptr;
    if (!png || !info)
    {
        png_destroy_write_struct(&png, &info);
        return false;
    }
    if (setjmp(png_jmpbuf(png)))
    {
        png_destroy_write_struct(&png, &info);
        std::cerr << "Failed to save image \"" << path << "\"\n";
        return false;
    }

    png_init_io(png, file.get());
    png_set_IHDR(png, info, width_, height_, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    png_write_image(png, rows.data());
    png_write_end(png, nullptr);
    png_destroy_write_struct(&png, &info);
    return true;
}
}
