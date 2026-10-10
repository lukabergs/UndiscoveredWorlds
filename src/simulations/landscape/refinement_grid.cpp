#include "simulations/landscape/refinement_grid.hpp"

#include <stdexcept>

// The accessors are inline in the header (step S0.4); the window construction, the
// snapping and the halo chain are step S2.1 (docs/architecture/multiscale-plan.md).

namespace landscape
{
LevelSpan snapoutwards(LevelSpan, int) { throw std::logic_error("not implemented: S2.1"); }

LevelWindow levelwindow(const RegionalWindow&, int, int) { throw std::logic_error("not implemented: S2.1"); }

std::vector<LevelWindow> halochain(const RegionalWindow&, int, int, int)
{
    throw std::logic_error("not implemented: S2.1");
}

HaloChainReport halochainreport(const RegionalWindow&, const std::vector<LevelWindow>&)
{
    throw std::logic_error("not implemented: S2.1");
}
}
