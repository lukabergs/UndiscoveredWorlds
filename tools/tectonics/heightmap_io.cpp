#include "heightmap_io.hpp"
#include "material.hpp"
#include "tectonic_contract.hpp"
#include "topography_tiff.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

namespace {

bool extract_json_number(const std::string& content, const char* key, uint32_t& value)
{
    const std::string token = std::string("\"") + key + "\"";
    size_t pos = content.find(token);
    if (pos == std::string::npos) {
        return false;
    }

    pos = content.find(':', pos + token.size());
    if (pos == std::string::npos) {
        return false;
    }
    ++pos;
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])) != 0) {
        ++pos;
    }

    size_t end = pos;
    while (end < content.size() && std::isdigit(static_cast<unsigned char>(content[end])) != 0) {
        ++end;
    }

    if (end == pos) {
        return false;
    }

    char* parse_end = nullptr;
    const unsigned long parsed = std::strtoul(content.substr(pos, end - pos).c_str(), &parse_end, 10);
    if (parse_end == nullptr || *parse_end != '\0' ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        return false;
    }

    value = static_cast<uint32_t>(parsed);
    return true;
}

bool extract_json_string(const std::string& content, const char* key, std::string& value)
{
    const std::string token = std::string("\"") + key + "\"";
    size_t pos = content.find(token);
    if (pos == std::string::npos) {
        return false;
    }

    pos = content.find(':', pos + token.size());
    if (pos == std::string::npos) {
        return false;
    }
    ++pos;
    while (pos < content.size() && std::isspace(static_cast<unsigned char>(content[pos])) != 0) {
        ++pos;
    }
    if (pos >= content.size() || content[pos] != '"') {
        return false;
    }
    ++pos;

    const size_t end = content.find('"', pos);
    if (end == std::string::npos) {
        return false;
    }

    value = content.substr(pos, end - pos);
    return true;
}

} // namespace

inline void setGray(png_byte *ptr, int val)
{
    ptr[0] = static_cast<png_byte>(val);
    ptr[1] = static_cast<png_byte>(val);
    ptr[2] = static_cast<png_byte>(val);
}

inline void setColor(png_byte *ptr, png_byte r, png_byte g, png_byte b)
{
    ptr[0] = r;
    ptr[1] = g;
    ptr[2] = b;
}

void hsvColor(png_byte* ptr, float hue_cycle, float saturation, float value)
{
    while (hue_cycle < 0.0f) {
        hue_cycle += 1.0f;
    }
    while (hue_cycle >= 1.0f) {
        hue_cycle -= 1.0f;
    }

    const float h = hue_cycle * 6.0f;
    const int sector = static_cast<int>(std::floor(h));
    const float fraction = h - static_cast<float>(sector);
    const float p = value * (1.0f - saturation);
    const float q = value * (1.0f - saturation * fraction);
    const float t = value * (1.0f - saturation * (1.0f - fraction));

    float r = value;
    float g = value;
    float b = value;
    switch (sector % 6) {
    case 0:
        r = value;
        g = t;
        b = p;
        break;
    case 1:
        r = q;
        g = value;
        b = p;
        break;
    case 2:
        r = p;
        g = value;
        b = t;
        break;
    case 3:
        r = p;
        g = q;
        b = value;
        break;
    case 4:
        r = t;
        g = p;
        b = value;
        break;
    default:
        r = value;
        g = p;
        b = q;
        break;
    }

    setColor(ptr, static_cast<png_byte>(255.0f * r), static_cast<png_byte>(255.0f * g),
             static_cast<png_byte>(255.0f * b));
}

void gradient(png_byte *ptr, png_byte ra, png_byte ga, png_byte ba, png_byte rb, png_byte gb,
              png_byte bb, float h, float ha, float hb);

struct ColorStop {
    float position;
    png_byte r;
    png_byte g;
    png_byte b;
};

void spectralColor(png_byte* ptr, float h)
{
    static const ColorStop kStops[] = {
        {0.00f,  0,   8,  92},
        {0.12f,  0,  84, 196},
        {0.24f,  0, 232, 255},
        {0.38f,  0, 168,  74},
        {0.52f, 244, 235,  60},
        {0.66f, 255, 164,  32},
        {0.78f, 228,  42,  24},
        {0.88f, 118,  72,  42},
        {0.95f, 132,  86, 178},
        {1.00f, 255, 255, 255},
    };

    if (h <= kStops[0].position) {
        setColor(ptr, kStops[0].r, kStops[0].g, kStops[0].b);
        return;
    }

    for (size_t i = 1; i < sizeof(kStops) / sizeof(kStops[0]); ++i) {
        if (h <= kStops[i].position) {
            gradient(ptr, kStops[i - 1].r, kStops[i - 1].g, kStops[i - 1].b, kStops[i].r,
                     kStops[i].g, kStops[i].b, h, kStops[i - 1].position, kStops[i].position);
            return;
        }
    }

    setColor(ptr, 255, 255, 255);
}

int writeImageRgb(const char* filename, int width, int height, const png_byte* rgb, const char* title)
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = const_cast<char*>(title);
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    for (int y = 0; y < height; ++y) {
        png_write_row(png_ptr, const_cast<png_bytep>(rgb + static_cast<size_t>(y) * static_cast<size_t>(width) * 3U));
    }
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) fclose(fp);
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

int writeImageRgba(const char* filename, int width, int height, const png_byte* rgba, const char* title)
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;

#ifdef _WIN32
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        return 1;
    }

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fclose(fp);
        return 1;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        png_structp png_ptr_nv = png_ptr;
        png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        fclose(fp);
        return 1;
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGBA,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = const_cast<char*>(title);
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    for (int y = 0; y < height; ++y) {
        png_write_row(
            png_ptr,
            const_cast<png_bytep>(rgba + static_cast<size_t>(y) * static_cast<size_t>(width) * 4U));
    }
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) {
        fclose(fp);
    }
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

int writeImage(const char* filename, int width, int height, float *heightmap, const char* title,
               void (drawFunction)(png_structp&, png_bytep&, int, int, float*))
{
    volatile int code = 0;
    FILE * volatile fp = nullptr;
    png_structp volatile png_ptr = nullptr;
    png_infop volatile info_ptr = nullptr;
    png_bytep volatile row = nullptr;
    size_t row_bytes = 0;  // Declare early to avoid goto issues

    // Open file for writing (binary mode)
#ifdef _WIN32
    // fopen_s doesn't accept volatile pointer, so use a non-volatile temporary
    FILE* fp_temp = nullptr;
    errno_t err = fopen_s(&fp_temp, filename, "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(filename, "wb");
    if (fp == nullptr) {
#endif
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        code = 1;
        goto finalise;
    }

    // Initialize write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fprintf(stderr, "Could not allocate write struct\n");
        code = 1;
        goto finalise;
    }

    // Initialize info structure
    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        fprintf(stderr, "Could not allocate info struct\n");
        code = 1;
        goto finalise;
    }

    // Setup Exception handling
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)  // Interaction between setjmp and C++ object destruction is acceptable here
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during png creation\n");
        code = 1;
        goto finalise;
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);

    // Write header (8-bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    // Set title
    if (title != nullptr) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = const_cast<char*>("Title");
        title_text.text = (char*)title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }

    png_write_info(png_ptr, info_ptr);

    // Allocate memory for one row (3 bytes per pixel - RGB)
    row_bytes = 3 * width * sizeof(png_byte);
    row = (png_bytep) malloc(row_bytes);

    if (row == nullptr) {
        fprintf(stderr, "Could not allocate memory for one row\n");
        code = 1;
        goto finalise;
    }

    // Write image data
    // Need to create non-volatile references for function calls
    {
        png_structp png_ptr_nv = png_ptr;
        png_bytep row_nv = row;
        drawFunction(png_ptr_nv, row_nv, width, height, heightmap);
    }

    // End write
    png_write_end(png_ptr, nullptr);

finalise:
    if (fp != nullptr) fclose(fp);
    if (row != nullptr) free(row);
    if (png_ptr != nullptr) {
        if (info_ptr != nullptr) {
            png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
        }
        {
            png_structp png_ptr_nv = png_ptr;
            png_destroy_write_struct(&png_ptr_nv, (png_infopp)nullptr);
        }
    }

    return code;
}

float find_value_for_quantile(const float quantile, const float* array, const uint32_t size)
{
    float value = 0.5;
    float th_step = 0.5;

    while (th_step > 0.00001)
    {
        uint32_t count = 0;
        for (uint32_t i = 0; i < size; ++i)
            count += (array[i] < value);

        th_step *= 0.5;
        if (count / (float)size < quantile)
            value += th_step;
        else
            value -= th_step;
    }
    return value;
}

void gradient(png_byte *ptr, png_byte ra, png_byte ga, png_byte ba, png_byte rb, png_byte gb, png_byte bb, float h, float ha, float hb)
{
    if (ha>hb) {
        printf("BAD1\n");
        throw runtime_error("BAD1\n");
    }
    if (hb<h) {
        printf("BAD2\n");
        throw runtime_error("BAD2\n");
    }
    if (ha>h) {
        printf("BAD3\n");
        throw runtime_error("BAD3\n");
    }
    float h_delta = hb - ha;
    float simil_b = (h - ha)/h_delta;
    float simil_a = (1.0f - simil_b);
    setColor(ptr,
             static_cast<png_byte>((float)simil_a * ra + (float)simil_b * rb),
             static_cast<png_byte>((float)simil_a * ga + (float)simil_b * gb),
             static_cast<png_byte>((float)simil_a * ba + (float)simil_b * bb));
}

void drawGrayImage(png_structp& png_ptr, png_bytep& row, int width, int height, float *heightmap)
{
    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {

            float h = heightmap[(y*width + x)];
            float res = 0.0f;
            if (h <= 0.0f) {
                res = 0;
            } else if (h >= 1.0f) {
                res = 255;
            } else {
                res = (h * 255.0f);
            }

            setGray(&(row[x*3]), static_cast<int>(res));
        }
        png_write_row(png_ptr, row);
    }
}

void renderImageGrayRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float h = heightmap[y * width + x];
            float res = 0.0f;
            if (h <= 0.0f) {
                res = 0.0f;
            } else if (h >= 1.0f) {
                res = 255.0f;
            } else {
                res = h * 255.0f;
            }

            png_byte* ptr = &rgb[(static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * 3U];
            setGray(ptr, static_cast<int>(res));
        }
    }
}

void drawColorsImage(png_structp& png_ptr, png_bytep& row, int width, int height, float *heightmap)
{
    int x, y;
    for (y=0 ; y<height ; y++) {
        for (x=0 ; x<width ; x++) {
            float h = heightmap[(y*width + x)];
            spectralColor(&(row[x*3]), h);
        }
        png_write_row(png_ptr, row);
    }
}

void renderImageColorsRgb(int width, int height, float* heightmap, std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float h = heightmap[y * width + x];
            png_byte* ptr = &rgb[(static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * 3U];
            spectralColor(ptr, h);
        }
    }
}

void renderMaterialMapRgb(int width, int height, const uint8_t* material_map,
                          std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            const auto& material = platec::material::properties_from_index(material_map[pixel_index]);
            png_byte* ptr = &rgb[pixel_index * 3U];
            setColor(ptr, material.r, material.g, material.b);
        }
    }
}

void renderCrustClassMapRgb(int width, int height, const uint8_t* crust_class_map,
                            std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            png_byte* ptr = &rgb[pixel_index * 3U];
            switch (static_cast<platec::contract::CrustClass>(crust_class_map[pixel_index])) {
            case platec::contract::CrustClass::Oceanic:
                setColor(ptr, 0x01, 0xBE, 0xBE);
                break;
            case platec::contract::CrustClass::Transitional:
                setColor(ptr, 0x90, 0xC9, 0x9D);
                break;
            case platec::contract::CrustClass::Continental:
                setColor(ptr, 0xF4, 0x98, 0x80);
                break;
            case platec::contract::CrustClass::None:
            default:
                setColor(ptr, 0, 0, 0);
                break;
            }
        }
    }
}

void renderCrustAgeMapRgb(int width, int height, const float* crust_age_myr_map,
                          std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);
    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    float max_age = 0.0f;
    for (size_t i = 0; i < pixel_count; ++i) {
        max_age = std::max(max_age, crust_age_myr_map[i]);
    }
    const float denom = std::log1pf(std::max(1.0f, max_age));

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            png_byte* ptr = &rgb[pixel_index * 3U];
            const float age = crust_age_myr_map[pixel_index];
            if (age <= 0.0f) {
                setColor(ptr, 0, 0, 0);
                continue;
            }

            const float normalized =
                denom > 0.0f ? std::log1pf(std::max(0.0f, age)) / denom : 0.0f;
            spectralColor(ptr, normalized);
        }
    }
}

void renderPlateIdMapRgb(int width, int height, const uint32_t* plate_map,
                         std::vector<png_byte>& rgb)
{
    rgb.resize(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U);
    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<uint32_t> active_ids;
    active_ids.reserve(32);

    for (size_t i = 0; i < pixel_count; ++i) {
        const uint32_t owner = plate_map[i];
        if (owner == 0xFFFFFFFFU) {
            continue;
        }
        if (std::find(active_ids.begin(), active_ids.end(), owner) == active_ids.end()) {
            active_ids.push_back(owner);
        }
    }

    std::sort(active_ids.begin(), active_ids.end());
    const uint32_t active_plate_count = static_cast<uint32_t>(active_ids.size());

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            const uint32_t owner = plate_map[pixel_index];
            png_byte* ptr = &rgb[pixel_index * 3U];

            if (owner == 0xFFFFFFFFU || active_plate_count == 0U) {
                setColor(ptr, 0, 0, 0);
                continue;
            }

            const auto it = std::lower_bound(active_ids.begin(), active_ids.end(), owner);
            if (it == active_ids.end() || *it != owner) {
                setColor(ptr, 0, 0, 0);
                continue;
            }

            const uint32_t plate_order =
                static_cast<uint32_t>(std::distance(active_ids.begin(), it));
            const float hue = active_plate_count > 1U
                ? static_cast<float>(plate_order) / static_cast<float>(active_plate_count)
                : 0.0f;
            hsvColor(ptr, hue, 0.90f, 1.0f);
        }
    }
}

int writeImageGray(const char* filename, int width, int height, float *heightmap, const char* title)
{
    return writeImage(filename, width, height, heightmap, title, drawGrayImage);
}

int writeImageColors(const char* filename, int width, int height, float *heightmap, const char* title)
{
    return writeImage(filename, width, height, heightmap, title, drawColorsImage);
}

int writeImageGrayTiff16(const char* filename, int width, int height, const uint16_t* heightmap)
{
    if (width <= 0 || height <= 0) {
        return 1;
    }

    return platec::io::writeGrayscaleTiff16(filename, static_cast<uint32_t>(width),
                                            static_cast<uint32_t>(height), heightmap)
               ? 0
               : 1;
}

int readImageGrayTiff16(const char* filename, std::vector<uint16_t>& heightmap, int& width,
                        int& height)
{
    uint32_t loaded_width = 0;
    uint32_t loaded_height = 0;
    if (!platec::io::readGrayscaleTiff16(filename, heightmap, loaded_width, loaded_height)) {
        heightmap.clear();
        width = 0;
        height = 0;
        return 1;
    }

    if (loaded_width > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
        loaded_height > static_cast<uint32_t>(std::numeric_limits<int>::max())) {
        heightmap.clear();
        width = 0;
        height = 0;
        return 1;
    }

    width = static_cast<int>(loaded_width);
    height = static_cast<int>(loaded_height);
    return 0;
}

int writeTopographyMetadataJson(const char* filename, const TopographyCodec::Metadata& metadata)
{
    std::ofstream output(filename, std::ios::binary);
    if (!output) {
        return 1;
    }

    output << "{\n";
    output << "  \"width\": " << metadata.width << ",\n";
    output << "  \"height\": " << metadata.height << ",\n";
    output << "  \"sea_level_m\": " << metadata.sea_level_m << ",\n";
    output << "  \"max_height_m\": " << metadata.max_height_m << ",\n";
    output << "  \"format\": \"" << metadata.format << "\",\n";
    output << "  \"endianness\": \"" << metadata.endianness << "\",\n";
    output << "  \"layout\": \"" << metadata.layout << "\",\n";
    output << "  \"version\": " << metadata.version << "\n";
    output << "}\n";

    return output ? 0 : 1;
}

int readTopographyMetadataJson(const char* filename, TopographyCodec::Metadata& metadata)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        return 1;
    }

    const std::string content((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    uint32_t width = 0;
    uint32_t height = 0;
    uint32_t sea_level_m = 0;
    uint32_t max_height_m = 0;
    uint32_t version = 0;
    std::string format;
    std::string endianness;
    std::string layout;

    if (!extract_json_number(content, "width", width) ||
        !extract_json_number(content, "height", height) ||
        !extract_json_number(content, "sea_level_m", sea_level_m) ||
        !extract_json_number(content, "max_height_m", max_height_m) ||
        !extract_json_number(content, "version", version) ||
        !extract_json_string(content, "format", format) ||
        !extract_json_string(content, "endianness", endianness) ||
        !extract_json_string(content, "layout", layout) ||
        sea_level_m > TopographyCodec::kMaxHeightMeters ||
        max_height_m > TopographyCodec::kMaxHeightMeters) {
        return 1;
    }

    metadata.width = width;
    metadata.height = height;
    metadata.sea_level_m = static_cast<uint16_t>(sea_level_m);
    metadata.max_height_m = static_cast<uint16_t>(max_height_m);
    metadata.version = version;
    metadata.format = format;
    metadata.endianness = endianness;
    metadata.layout = layout;
    return 0;
}

