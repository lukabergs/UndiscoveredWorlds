#include "simulations/landscape/refinement_projection.hpp"

#include <stdexcept>

// The projection, its error measure and the resample are step S2.2
// (docs/architecture/multiscale-plan.md); S0.4 wrote their contracts.

namespace landscape
{
std::vector<float> projectionweights(const LevelWindow&, const std::vector<std::uint8_t>&, const std::vector<std::uint8_t>&,
                                     const std::vector<std::uint8_t>&, double)
{
    throw std::logic_error("not implemented: S2.2");
}

ProjectionReport projectlevel(const LevelWindow&, std::vector<float>&, const std::vector<float>&, const LevelWindow&,
                              const std::vector<float>&)
{
    throw std::logic_error("not implemented: S2.2");
}

double levelparentmeanerror(const LevelWindow&, const std::vector<float>&, const LevelWindow&, const std::vector<float>&)
{
    throw std::logic_error("not implemented: S2.2");
}

ResampledField resampleto(const std::vector<LevelWindow>&, const std::vector<std::vector<float>>&, double, double)
{
    throw std::logic_error("not implemented: S2.2");
}
}
