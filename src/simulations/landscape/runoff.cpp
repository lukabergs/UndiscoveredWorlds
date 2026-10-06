#include "simulations/landscape/runoff.hpp"

#include <stdexcept>

namespace landscape
{
// Erosion plan C2 implements both functions.
double budykorunoff(double, double, double)
{
    throw std::logic_error("budykorunoff: not implemented (erosion plan C2)");
}

double accumulateflowwithlosses(const SphereMetric&, const std::vector<float>&, const std::vector<float>&,
                                const TransmissionLoss&, FlowGraph&)
{
    throw std::logic_error("accumulateflowwithlosses: not implemented (erosion plan C2)");
}
}
