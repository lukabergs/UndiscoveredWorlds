#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace climateocean::detail
{
// The upwind PV terms make this matrix nonsymmetric. BiCGSTAB therefore uses
// a forward/backward Gauss-Seidel preconditioner, rather than assuming SPD/CG.
// Algorithm: Barrett et al., Templates for the Solution of Linear Systems,
// section 2.3.8. Acceptance always uses the original, explicitly recomputed residual.
struct CirculationStencil
{
    int columns;
    const std::vector<bool>& active;
    const std::vector<double>& diagonal;
    const std::vector<double>& east;
    const std::vector<double>& west;
    const std::vector<double>& north;
    const std::vector<double>& south;

    std::size_t e(std::size_t i) const { return i % columns + 1 < static_cast<std::size_t>(columns) ? i + 1 : i + 1 - columns; }
    std::size_t w(std::size_t i) const { return i % columns ? i - 1 : i + columns - 1; }
    void multiply(const std::vector<double>& x, std::vector<double>& out) const
    {
        for (std::size_t i = 0; i < active.size(); ++i)
            out[i] = active[i] ? diagonal[i] * x[i] - east[i] * x[e(i)] - west[i] * x[w(i)] -
                north[i] * x[i - columns] - south[i] * x[i + columns] : 0.0;
    }
    void precondition(const std::vector<double>& b, std::vector<double>& x) const
    {
        // Include longitude-seam entries on their correct side of the global
        // triangular ordering. Inactive coast/pole vertices stay exactly zero.
        for (std::size_t i = 0; i < active.size(); ++i)
        {
            if (!active[i]) { x[i] = 0.0; continue; }
            const auto left = w(i), right = e(i);
            x[i] = (b[i] + north[i] * x[i - columns] +
                (left < i ? west[i] * x[left] : 0.0) + (right < i ? east[i] * x[right] : 0.0)) / diagonal[i];
        }
        for (std::size_t i = active.size(); i-- > 0;)
        {
            if (!active[i]) continue;
            const auto left = w(i), right = e(i);
            x[i] += (south[i] * x[i + columns] +
                (left > i ? west[i] * x[left] : 0.0) + (right > i ? east[i] * x[right] : 0.0)) / diagonal[i];
        }
    }
};

struct CirculationSolveResult { double relativeResidual = 0.0; int iterations = 0; };

inline CirculationSolveResult solveCirculationKrylov(const CirculationStencil& a,
    const std::vector<double>& b, std::vector<double>& x, int limit, double tolerance)
{
    const auto size = b.size();
    std::vector<double> r(size), shadow(size), p(size), v(size), s(size), t(size), phat(size), shat(size), ax(size);
    const auto dot = [](const std::vector<double>& u, const std::vector<double>& v)
    { double sum = 0.0; for (std::size_t i = 0; i < u.size(); ++i) sum += u[i] * v[i]; return sum; };
    const double norm = dot(b, b);
    CirculationSolveResult result;
    if (norm == 0.0) { std::fill(x.begin(), x.end(), 0.0); return result; }
    const auto measure = [&]
    {
        a.multiply(x, ax);
        double residual = 0.0;
        for (std::size_t i = 0; i < size; ++i) residual += (b[i] - ax[i]) * (b[i] - ax[i]);
        result.relativeResidual = std::sqrt(residual / norm);
        return result.relativeResidual <= tolerance;
    };
    if (measure()) return result;
    for (std::size_t i = 0; i < size; ++i) r[i] = shadow[i] = b[i] - ax[i];
    double oldRho = 1.0, alpha = 1.0, omega = 1.0;
    for (int iteration = 0; iteration < std::max(1, limit); ++iteration)
    {
        const double rho = dot(shadow, r);
        if (rho == 0.0 || !std::isfinite(rho) || omega == 0.0) break;
        const double beta = (rho / oldRho) * (alpha / omega);
        for (std::size_t i = 0; i < size; ++i) p[i] = iteration == 0 ? r[i] : r[i] + beta * (p[i] - omega * v[i]);
        a.precondition(p, phat); a.multiply(phat, v);
        const double denominator = dot(shadow, v);
        if (denominator == 0.0 || !std::isfinite(denominator)) break;
        alpha = rho / denominator;
        for (std::size_t i = 0; i < size; ++i) { s[i] = r[i] - alpha * v[i]; x[i] += alpha * phat[i]; }
        result.iterations = iteration + 1;
        if (measure()) return result;
        a.precondition(s, shat); a.multiply(shat, t);
        const double tt = dot(t, t);
        if (tt == 0.0 || !std::isfinite(tt)) break;
        omega = dot(t, s) / tt;
        if (!std::isfinite(omega)) break;
        for (std::size_t i = 0; i < size; ++i) { x[i] += omega * shat[i]; r[i] = s[i] - omega * t[i]; }
        if (measure()) return result;
        oldRho = rho;
    }
    measure();
    return result;
}
}
