#pragma once

#include <vector>

// Experimental pre-mode-separated thermal closures. Kept for controlled comparisons;
// production uses hypsometric local forcing. See README.md in this directory.
namespace climateatmosphere
{
float thermalSurfacePressureAnomalyHpa(
    float temperatureAnomalyK,
    float referencePressureHpa,
    float referenceTemperatureK,
    float massRedistributionEfficiency);
std::vector<float> nonlocalThermalResponse(
    int longitudeCells,
    int latitudeCells,
    const std::vector<float>& localTemperatureAnomalyK,
    float tropicalLatitudeDegrees,
    float longitudinalReachDegrees,
    float meridionalReachDegrees,
    float rotationDirection = 1.0f);
}
