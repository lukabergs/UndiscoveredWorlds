#include "ocean_feedback.hpp"
#include "climate_ocean_dynamics.hpp"
#include "climate_grid.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateocean::legacy
{
void applyEmpiricalSstWindFeedback(const climategrid::SphericalGrid& grid,
    const OceanForcing& forcing, const OceanConfig& config, float relaxation,
    OceanState& state, double& residual)
{
    const int columns = grid.columns;
    const int rows = grid.rows;
    const double dy = config.planetRadiusMetres * grid.latitudeSpacingRadians;
    const auto ocean = [&](std::size_t cell) { return !forcing.landMask[cell] && forcing.bathymetryMetres[cell] > 0.0f; };
    std::vector<float> targetPressure(state.sstC.size(), 0.0f);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto c = grid.index(x, y);
            if (!ocean(c)) continue;
            targetPressure[c] = -config.sstWindFeedbackMpsPerK * 10.0f *
                (state.surfaceSkinTemperatureC[c] - forcing.atmosphericTemperatureC[c]);
        }
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto c = grid.index(x, y);
            if (!ocean(c)) continue;
            const auto at = [&](int xx, int yy)
            {
                const auto other = grid.index(xx, std::clamp(yy, 0, rows - 1));
                return ocean(other) ? targetPressure[other] : targetPressure[c];
            };
            const double dx = grid.cellAreasSquareMetres[y] / grid.zonalFaceLengthsMetres[y];
            const double eastForce = -(at(x + 1, y) - at(x - 1, y)) * 100.0 / (2.0 * dx * config.airDensityKgM3);
            const double northForce = -(at(x, y - 1) - at(x, y + 1)) * 100.0 / (2.0 * dy * config.airDensityKgM3);
            const double f = 2.0 * config.rotationRatePerSecond * config.rotationDirection * std::sin(grid.latitudeCentresRadians[y]);
            const double drag = 1.0 / 43200.0;
            const double deltaU = forcing.eastWindMps[c] + (drag * eastForce + f * northForce) / (drag * drag + f * f) - state.coupledEastWindMps[c];
            const double deltaV = forcing.southWindMps[c] - (drag * northForce - f * eastForce) / (drag * drag + f * f) - state.coupledSouthWindMps[c];
            residual += grid.cellAreasSquareMetres[y] * (deltaU * deltaU + deltaV * deltaV); // 1 m/s reference
            state.coupledEastWindMps[c] += relaxation * static_cast<float>(deltaU);
            state.coupledSouthWindMps[c] += relaxation * static_cast<float>(deltaV);
            state.coupledPressureAnomalyHpa[c] += relaxation * (targetPressure[c] - state.coupledPressureAnomalyHpa[c]);
        }
}
}
