#include "simulations/landscape/regional_tributaries.hpp"

#include "simulations/landscape/regional_pass.hpp"

// Multiscale plan S3.0: the hook, inert (S3.3 grows the tributaries).

namespace landscape
{
void tributarieslevel(const LevelContext& context, LevelState& state)
{
    static_cast<void>(context);
    state.tributaries = 0;
    state.tributarycells = 0;
}
}
