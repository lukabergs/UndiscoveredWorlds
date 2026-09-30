#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateocean::detail
{
// Backward Euler for liquid heat on a periodic latitude row. With h=E/C
// relative to freezing, theta=max(h,0): ice stores negative h, never negative
// liquid temperature. Zero face conductance closes coasts without special seams.
struct ZonalDiffusionWorkspace
{
    std::vector<double> initial, rate, capacityRatio, lower, diagonal, upper, rhs, theta, auxiliary, inverse, modifiedUpper;
    std::vector<unsigned char> liquid;

    explicit ZonalDiffusionWorkspace(int n)
        : initial(n), rate(n), capacityRatio(n), lower(n), diagonal(n), upper(n), rhs(n), theta(n), auxiliary(n),
          inverse(n), modifiedUpper(n), liquid(n) {}

    bool solveCyclic()
    {
        const int n = static_cast<int>(theta.size());
        const double alpha = lower[0], beta = upper[n - 1], gamma = -diagonal[0];
        // Sherman-Morrison reduction to two solves with the same tridiagonal
        // factorization. alpha is the top-right, beta the bottom-left entry.
        for (int i = 0; i < n; ++i)
        {
            double pivot = diagonal[i];
            if (i == 0) pivot -= gamma;
            if (i == n - 1) pivot -= alpha * beta / gamma;
            if (i > 0) pivot -= lower[i] * modifiedUpper[i - 1];
            if (!(pivot > 0.0) || !std::isfinite(pivot)) return false;
            inverse[i] = 1.0 / pivot;
            modifiedUpper[i] = i + 1 < n ? upper[i] * inverse[i] : 0.0;
            theta[i] = (rhs[i] - (i > 0 ? lower[i] * theta[i - 1] : 0.0)) * inverse[i];
            const double u = i == 0 ? gamma : (i == n - 1 ? beta : 0.0);
            auxiliary[i] = (u - (i > 0 ? lower[i] * auxiliary[i - 1] : 0.0)) * inverse[i];
        }
        for (int i = n - 2; i >= 0; --i)
        {
            theta[i] -= modifiedUpper[i] * theta[i + 1];
            auxiliary[i] -= modifiedUpper[i] * auxiliary[i + 1];
        }
        const double denominator = 1.0 + auxiliary[0] + alpha * auxiliary[n - 1] / gamma;
        if (!std::isfinite(denominator) || std::abs(denominator) < 1.0e-14) return false;
        const double correction = (theta[0] + alpha * theta[n - 1] / gamma) / denominator;
        for (int i = 0; i < n; ++i) theta[i] -= correction * auxiliary[i];
        return true;
    }

    bool solve(const double* enthalpy, const double* conductance, double heatCapacity,
        double dtOverArea, int& maximumIterations, double& maximumResidualK, const double* capacities = nullptr)
    {
        const int n = static_cast<int>(theta.size());
        bool hasLiquid = false, hasDiffusion = false;
        double scale = 1.0;
        for (int i = 0; i < n; ++i)
        {
            capacityRatio[i] = capacities ? heatCapacity / capacities[i] : 1.0;
            initial[i] = enthalpy[i] / (capacities ? capacities[i] : heatCapacity);
            rate[i] = dtOverArea * conductance[i];
            theta[i] = std::max(0.0, initial[i]);
            liquid[i] = initial[i] > 0.0;
            hasLiquid = hasLiquid || liquid[i];
            hasDiffusion = hasDiffusion || rate[i] > 0.0;
            scale = std::max(scale, std::abs(initial[i]));
        }
        if (!hasLiquid || !hasDiffusion) return true;
        // Diffusion alone cannot freeze an initially liquid cell: its coldest
        // donor is at freezing. Start with liquid cells, then activate ice cells
        // only when the received heat exceeds their latent deficit. This set
        // grows monotonically, terminating after at most n activations.
        for (int iteration = 1; iteration <= n + 1; ++iteration)
        {
            maximumIterations = std::max(maximumIterations, iteration);
            for (int i = 0; i < n; ++i)
            {
                const int west = (i + n - 1) % n;
                lower[i] = liquid[i] ? -capacityRatio[i] * rate[west] : 0.0;
                upper[i] = liquid[i] ? -capacityRatio[i] * rate[i] : 0.0;
                diagonal[i] = liquid[i] ? 1.0 + capacityRatio[i] * (rate[west] + rate[i]) : 1.0;
                rhs[i] = liquid[i] ? initial[i] : 0.0;
            }
            if (!solveCyclic()) return false;
            bool activated = false;
            double residual = 0.0;
            for (int i = 0; i < n; ++i)
            {
                const int west = (i + n - 1) % n, east = (i + 1) % n;
                const double h = initial[i] + capacityRatio[i] * (rate[west] * (theta[west] - theta[i]) +
                    rate[i] * (theta[east] - theta[i]));
                if (!std::isfinite(h) || !std::isfinite(theta[i])) return false;
                residual = std::max(residual, std::abs(std::max(0.0, h) - theta[i]));
                if (!liquid[i] && h > 0.0) { liquid[i] = true; activated = true; }
            }
            if (!activated)
            {
                maximumResidualK = std::max(maximumResidualK, residual);
                return residual <= 1.0e-9 * scale;
            }
        }
        return false;
    }
};
}
