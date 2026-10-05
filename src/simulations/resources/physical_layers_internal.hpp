#ifndef physical_layers_internal_hpp
#define physical_layers_internal_hpp

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>
#include "core/grid_coordinates.hpp"
#include "core/parallel_rows.hpp"
#include "world/planet.hpp"

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
    return world.truelake(x, y) != 0 || world.riftlakesurfacefloat(x, y) != 0;
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

// Downstream cell of each land or delta cell, row-major with the planet's
// inclusive bounds; -1 where water leaves the map or no lower neighbour exists.
// Rivers and deltas give the direction; elsewhere the steepest descent does.
inline std::vector<int32_t> drainagereceivers(const planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const size_t columns = size_t(width) + 1;
    std::vector<int32_t> receiver(columns * (size_t(height) + 1), -1);
    for (int y = 0; y <= height; ++y)
    {
        for (int x = 0; x <= width; ++x)
        {
            if (world.sea(x, y) != 0 && !isdeltacell(world, x, y)) continue;
            int direction = world.riverdir(x, y);
            if (world.deltadir(x, y) > 0)
                direction = world.deltadir(x, y);
            if (direction < 1 || direction > 8)
            {
                int lowest = world.nom(x, y);
                for (int candidate = 1; candidate <= 8; ++candidate)
                {
                    const auto next = getdestination(x, y, candidate);
                    if (next.y < 0 || next.y > height) continue;
                    const int elevation = world.nom(wrap(next.x, width), next.y);
                    if (elevation < lowest) { lowest = elevation; direction = candidate; }
                }
            }
            if (direction < 1 || direction > 8) continue;
            const auto next = getdestination(x, y, direction);
            if (next.y < 0 || next.y > height) continue;
            const size_t index = size_t(y) * columns + size_t(x);
            const size_t destination = size_t(next.y) * columns + size_t(wrap(next.x, width));
            if (destination != index) receiver[index] = static_cast<int32_t>(destination);
        }
    }
    return receiver;
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
