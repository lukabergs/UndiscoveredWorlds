#pragma once

#include "core/deterministic_random.hpp"

#include <cmath>
#include <cstdint>

namespace mineralsystems::detail
{
// Counter-based stream; results depend only on the seed and call sequence.
class Random
{
public:
    explicit Random(std::uint64_t seed) : state_(seed) {}

    std::uint64_t next() { return deterministicmix64(state_++); }
    // Uniform in [0, 1).
    double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
    double uniform(double minimum, double maximum) { return minimum + (maximum - minimum) * uniform(); }
    double normal()
    {
        constexpr double tau = 6.28318530717958647692;
        const double radius = std::sqrt(-2.0 * std::log(1.0 - uniform()));
        return radius * std::cos(tau * uniform());
    }
    double lognormal(double median, double sigmaLn) { return median * std::exp(sigmaLn * normal()); }
    // Exact for any mean: sums Knuth samples over chunks of at most 16.
    std::uint64_t poisson(double mean)
    {
        std::uint64_t count = 0;
        while (mean > 0.0)
        {
            const double chunk = mean > 16.0 ? 16.0 : mean;
            mean -= chunk;
            const double limit = std::exp(-chunk);
            double product = uniform();
            while (product > limit)
            {
                ++count;
                product *= uniform();
            }
        }
        return count;
    }

private:
    std::uint64_t state_;
};
}
