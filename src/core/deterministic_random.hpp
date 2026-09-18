#pragma once

#include <algorithm>
#include <cstdint>

// Stable seed derivation shared by simulation stages. Changing these formulas
// changes generated worlds for existing seeds.
inline std::uint64_t deterministicmix64(std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ull;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
    return value ^ (value >> 31);
}

inline std::uint64_t deterministiccontextseed(long baseseed, int salt)
{
    std::uint64_t value = static_cast<std::uint64_t>(static_cast<std::int64_t>(baseseed));
    value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(salt)) * 0x9e3779b97f4a7c15ull;
    return deterministicmix64(value);
}

inline long deterministicfastseed(std::uint64_t seed)
{
    return static_cast<long>(seed & 0x7fffffffull);
}

template<typename... Ints>
inline std::uint64_t deterministiccombine(std::uint64_t seed, Ints... values)
{
    ((seed = deterministicmix64(seed ^ static_cast<std::uint64_t>(static_cast<std::int64_t>(values)))), ...);
    return seed;
}

template<typename... Ints>
inline int deterministicrandom(std::uint64_t seed, int minimum, int maximum, Ints... values)
{
    if (maximum < minimum)
        std::swap(minimum, maximum);

    const std::uint64_t combined = deterministiccombine(seed, values...);
    const std::uint64_t span = static_cast<std::uint64_t>(maximum - minimum) + 1ull;
    return minimum + static_cast<int>(combined % span);
}

template<typename... Ints>
inline int deterministicsignedrandom(std::uint64_t seed, int minimumabs, int maximumabs, Ints... values)
{
    const int magnitude = deterministicrandom(seed ^ 0x6a09e667f3bcc909ull, minimumabs, maximumabs, values...);
    return deterministicrandom(seed ^ 0xbb67ae8584caa73bull, 0, 1, values...) == 0 ? -magnitude : magnitude;
}

