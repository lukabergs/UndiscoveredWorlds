#include "simulations/landscape/glacial.hpp"

#include <stdexcept>

// Contracts of erosion plan step G0; the bodies are step G1.
namespace landscape
{
double equilibriumlinetemperature(double)
{
    throw std::logic_error("landscape::equilibriumlinetemperature: not implemented (erosion plan G1)");
}

void equilibriumline(const SphereMetric&, const std::vector<float>&, const std::vector<float>&,
                     const std::vector<float>&, double, double, std::vector<float>&)
{
    throw std::logic_error("landscape::equilibriumline: not implemented (erosion plan G1)");
}

void accumulateice(const SphereMetric&, const FlowGraph&, const std::vector<float>&, const std::vector<float>&,
                   const std::vector<float>&, double, std::vector<double>&)
{
    throw std::logic_error("landscape::accumulateice: not implemented (erosion plan G1)");
}

void icethickness(const std::vector<double>&, const GlacialParameters&, std::vector<float>&)
{
    throw std::logic_error("landscape::icethickness: not implemented (erosion plan G1)");
}
}
