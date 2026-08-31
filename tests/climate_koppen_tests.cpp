#include "climate_koppen.hpp"

#include <cstdlib>
#include <iostream>

namespace
{
void expect(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}
}

int main()
{
    expect(
        climatekoppen::isWinterDry(9.9f, 100.0f),
        "winter precipitation below one tenth of the summer maximum must be dry-winter");
    expect(
        !climatekoppen::isWinterDry(10.0f, 100.0f),
        "the one-tenth boundary itself must not be dry-winter");
    expect(
        !climatekoppen::isWinterDry(25.0f, 100.0f),
        "the obsolete one-quarter temperate rule must not classify dry-winter");
    return 0;
}
