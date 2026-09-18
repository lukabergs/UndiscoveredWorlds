#include "../../src/validation/climate/climate_palette.hpp"
#include <iostream>
#include <limits>

int main()
{
    using climatepalette::colour;
    int failures = 0;
    auto check = [&](bool condition) { if (!condition) ++failures; };
    check(colour(0, 0, 25, false) == std::array<unsigned char, 3>{40, 71, 199});
    check(colour(25, 0, 25, false) == std::array<unsigned char, 3>{207, 41, 41});
    check(colour(-100, -30, 30, true) == std::array<unsigned char, 3>{225, 225, 225});
    check(colour(0, -54, 40, true) == std::array<unsigned char, 3>{0, 220, 230});
    check(colour(0, -30, 30, true) == colour(0, -54, 40, true));
    check(colour(std::numeric_limits<double>::quiet_NaN(), 0, 1, false) == std::array<unsigned char, 3>{0, 0, 0});
    for (const auto& stop : climatepalette::positiveStops)
        check(colour(stop.position*25, 0, 25, false) == stop.rgb);
    std::cout << "Palette contract failures: " << failures << '\n';
    return failures ? 1 : 0;
}
