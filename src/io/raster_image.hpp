#pragma once

// RGBA8 image with PNG load/save through libpng. Headless replacement for the
// subset of sf::Image used by imports and validation; method names follow
// sf::Image so callers read the same. Decoding matches the desktop loader:
// palette, grey and tRNS expand to RGBA, 16-bit samples keep their high byte.

#include <cstdint>
#include <string>
#include <vector>

namespace rasterimage
{
using Uint8 = std::uint8_t;

struct Color
{
    Uint8 r = 0;
    Uint8 g = 0;
    Uint8 b = 0;
    Uint8 a = 255;

    constexpr Color() = default;
    constexpr Color(Uint8 red, Uint8 green, Uint8 blue, Uint8 alpha = 255) : r(red), g(green), b(blue), a(alpha) {}

    static const Color Black;
    static const Color White;

    friend constexpr bool operator==(const Color& left, const Color& right)
    {
        return left.r == right.r && left.g == right.g && left.b == right.b && left.a == right.a;
    }
    friend constexpr bool operator!=(const Color& left, const Color& right) { return !(left == right); }
};

inline constexpr Color Color::Black{ 0, 0, 0 };
inline constexpr Color Color::White{ 255, 255, 255 };

struct Size
{
    unsigned int x = 0;
    unsigned int y = 0;

    friend constexpr bool operator==(const Size& left, const Size& right) { return left.x == right.x && left.y == right.y; }
    friend constexpr bool operator!=(const Size& left, const Size& right) { return !(left == right); }
};

class Image
{
public:
    void create(unsigned int width, unsigned int height, const Color& color = Color(0, 0, 0));
    Size getSize() const { return { width_, height_ }; }
    Color getPixel(unsigned int x, unsigned int y) const;
    void setPixel(unsigned int x, unsigned int y, const Color& color);
    const Uint8* getPixelsPtr() const { return pixels_.empty() ? nullptr : pixels_.data(); }
    // PNG only; failures print the reason to standard error like the desktop loader.
    bool loadFromFile(const std::string& path);
    bool saveToFile(const std::string& path) const;

private:
    unsigned int width_ = 0;
    unsigned int height_ = 0;
    std::vector<Uint8> pixels_;
};
}
