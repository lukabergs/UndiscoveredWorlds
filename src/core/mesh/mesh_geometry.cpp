#include "core/mesh/mesh_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace mesh
{
namespace
{
constexpr double pi = 3.14159265358979323846;

double localspacing(const CellGrid& grid, CellId c)
{
    double spacing = grid.distance_m(c, 1);
    for (int k = 3; k <= slotcount; k += 2)
        spacing = std::min(spacing, grid.distance_m(c, k));
    return spacing > 0.0 ? spacing : grid.minimumspacing_m();
}
}

double greatcircle_rad(const Vec3& p, const Vec3& q) { return std::atan2(length(cross(p, q)), dot(p, q)); }
double greatcircle_m(const Vec3& p, const Vec3& q, double radius_m) { return radius_m * greatcircle_rad(p, q); }

std::vector<CellId> cellswithin(const CellGrid& grid, const Vec3& centre, double radius_m)
{
    std::vector<CellId> result;
    const CellId start = grid.locate(centre);
    if (start == nocell || !(radius_m >= 0.0))
        return result;
    const double limit = radius_m / grid.radius_m();
    std::unordered_set<CellId> visited{ start };
    std::vector<CellId> frontier{ start };
    while (!frontier.empty())
    {
        const CellId c = frontier.back();
        frontier.pop_back();
        const bool inside = greatcircle_rad(centre, grid.centre(c)) <= limit;
        if (inside)
            result.push_back(c);
        else if (c != start)
            continue;
        for (int k = 1; k <= slotcount; k++)
        {
            const std::int32_t other = grid.neighbour(c, k);
            if (other < 0)
                continue;
            const CellId o = static_cast<CellId>(other);
            if (visited.insert(o).second)
                frontier.push_back(o);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<CellId> cellsalongarc(const CellGrid& grid, const Vec3& p, const Vec3& q)
{
    std::vector<CellId> result;
    const auto add = [&](CellId c) {
        if (c != nocell && (result.empty() || result.back() != c))
            result.push_back(c);
    };
    const double total = greatcircle_rad(p, q);
    const Vec3 w0 = q - dot(p, q) * p;
    const double wlength = length(w0);
    CellId current = grid.locate(p);
    add(current);
    if (total == 0.0 || wlength == 0.0)
    {
        add(grid.locate(q));
        return result;
    }
    const Vec3 w = (1.0 / wlength) * w0;
    double s = 0.0;
    const double radius = grid.radius_m();
    while (true)
    {
        const double spacing = current == nocell ? grid.minimumspacing_m() : localspacing(grid, current);
        s += 0.25 * spacing / radius;
        if (s >= total)
            break;
        const Vec3 point = std::cos(s) * p + std::sin(s) * w;
        const CellId c = grid.locate(point);
        add(c);
        if (c != nocell)
            current = c;
    }
    add(grid.locate(q));
    return result;
}

TangentFrame tangentframe(const Vec3& centre, double radius_m)
{
    const Vec3 origin = normalise(centre);
    const double latitude = latitudeof(origin);
    const double longitude = latitude >= 0.5 * pi || latitude <= -0.5 * pi ? 0.0 : longitudeof(origin);
    return { origin, eastof(latitude, longitude), northof(latitude, longitude), radius_m };
}

Vec3 fromtangent(const TangentFrame& frame, double east_km, double north_km)
{
    const double scale = 1000.0 / frame.radius_m;
    return normalise(frame.origin + (east_km * scale) * frame.east + (north_km * scale) * frame.north);
}

Vec3 uniformpoint(double u, double v)
{
    const double z = 2.0 * u - 1.0;
    const double longitude = 2.0 * pi * v - pi;
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    return { r * std::cos(longitude), r * std::sin(longitude), z };
}
}
