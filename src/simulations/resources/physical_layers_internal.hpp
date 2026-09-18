#ifndef physical_layers_internal_hpp
#define physical_layers_internal_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include "grid_coordinates.hpp"
#include "parallel_rows.hpp"
#include "planet.hpp"

// Shared planet-grid queries for derived physical layers. Coordinates use the
// planet's inclusive bounds; longitude wraps and latitude never wraps.
namespace physical_layers::detail
{
constexpr std::array<std::array<int, 2>, 8> drainageoffsets =
{ {
    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
} };

inline twointegers wrappeddestination(const planet& world, int x, int y, int dir)
{
    twointegers destination = getdestination(x, y, dir);
    destination.x = wrap(destination.x, world.width());
    destination.y = std::clamp(destination.y, 0, world.height());
    return destination;
}

inline float clampunit(float value)
{
    return std::clamp(value, 0.0f, 1.0f);
}

inline int clampscore(float value)
{
    return static_cast<int>(std::round(std::clamp(value, 0.0f, 100.0f)));
}

inline bool islakecell(const planet& world, int x, int y)
{
    return world.truelake(x, y) != 0 || world.riftlakesurface(x, y) != 0;
}

inline bool issaltfeature(const planet& world, int x, int y)
{
    const int special = world.special(x, y);
    return special == 100 || special == 110;
}

inline bool isdeltacell(const planet& world, int x, int y)
{
    return world.deltadir(x, y) != 0 || world.deltajan(x, y) != 0 || world.deltajul(x, y) != 0;
}

inline int annualrainfall(const planet& world, int x, int y)
{
    return (world.summerrain(x, y) + world.winterrain(x, y)) / 2;
}

inline int averageriverflow(const planet& world, int x, int y)
{
    const int riverflow = world.riveraveflow(x, y);
    const int deltaflow = (std::abs(world.deltajan(x, y)) + std::abs(world.deltajul(x, y))) / 2;
    return std::max(riverflow, deltaflow);
}

template<typename Predicate>
inline bool hasfeatureinradius(const planet& world, int x, int y, int radius, Predicate predicate)
{
    const int width = world.width();
    const int height = world.height();

    for (int dy = -radius; dy <= radius; dy++)
    {
        const int ny = y + dy;

        if (ny < 0 || ny > height)
            continue;

        for (int dx = -radius; dx <= radius; dx++)
        {
            const int nx = wrap(x + dx, width);

            if (predicate(nx, ny))
                return true;
        }
    }

    return false;
}
}

#endif
