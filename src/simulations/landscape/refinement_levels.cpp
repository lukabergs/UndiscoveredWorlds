#include "simulations/landscape/refinement_levels.hpp"

#include <cmath>
#include <stdexcept>

// The cell size and the level assignment are step S0.4; the loop is step S2.3
// (docs/architecture/multiscale-plan.md).

namespace landscape
{
namespace
{
constexpr double pi = 3.14159265358979323846;
}

double levelcellsize_km(double radius_m, int worldwidth, int factor, int level)
{
    if (!std::isfinite(radius_m) || !(radius_m > 0.0))
        throw std::invalid_argument("levelcellsize_km: the radius must be positive and finite");
    if (worldwidth < 4 || worldwidth % 2 != 0)
        throw std::invalid_argument("levelcellsize_km: the world width must be even and at least 4");
    if (factor < 2 || level < 0)
        throw std::invalid_argument("levelcellsize_km: the factor must be at least 2 and the level not negative");
    const double rows = static_cast<double>(worldwidth / 2) * static_cast<double>(levelscale(factor, level));
    return pi * radius_m / rows / 1000.0;
}

int detaillevel(double size_km, double detailcells, double radius_m, int worldwidth, int factor, int levels)
{
    if (!std::isfinite(size_km) || !(size_km > 0.0))
        throw std::invalid_argument("detaillevel: the detail's size must be positive and finite");
    if (!std::isfinite(detailcells) || !(detailcells > 0.0))
        throw std::invalid_argument("detaillevel: the cells per detail must be positive and finite");
    if (levels < 1)
        throw std::invalid_argument("detaillevel: at least one level must run");
    const double cell_km = size_km / detailcells;
    for (int level = 1; level <= levels; level++)
        if (levelcellsize_km(radius_m, worldwidth, factor, level) <= cell_km)
            return level;
    return levels;
}

RefinementResult runrefinement(const RegionalGlobalInputs&, const LevelTarget&, const RegionalParameters&,
                               const RefinementParameters&)
{
    throw std::logic_error("not implemented: S2.3");
}
}
