#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace physical_layers::detail {

// Potential, not transported mass: preserve the strongest connected source.
// Functional drainage graphs are processed in O(cells), including lake/river cycles.
inline std::vector<uint8_t> route_ore_potential(const std::vector<int32_t>& receiver,
                                               const std::vector<uint8_t>& source) {
    const size_t count = receiver.size();
    if (source.size() != count) throw std::invalid_argument("drainage/source size mismatch");
    std::vector<uint8_t> potential = source;
    std::vector<uint32_t> incoming(count, 0);
    auto valid = [count](int32_t next) { return next >= 0 && size_t(next) < count; };
    for (auto next : receiver) if (valid(next)) ++incoming[size_t(next)];
    std::vector<size_t> queue;
    queue.reserve(count);
    for (size_t i = 0; i < count; ++i) if (incoming[i] == 0) queue.push_back(i);
    for (size_t head = 0; head < queue.size(); ++head) {
        const size_t i = queue[head];
        if (!valid(receiver[i])) continue;
        const size_t next = size_t(receiver[i]);
        potential[next] = std::max(potential[next], potential[i]);
        if (--incoming[next] == 0) queue.push_back(next);
    }
    // Remaining nodes are disjoint directed cycles, with tributaries already accumulated.
    for (size_t i = 0; i < count; ++i) {
        if (incoming[i] == 0) continue;
        uint8_t strongest = 0;
        size_t node = i;
        do {
            strongest = std::max(strongest, potential[node]);
            node = size_t(receiver[node]);
        } while (node != i);
        do {
            potential[node] = strongest;
            incoming[node] = 0;
            node = size_t(receiver[node]);
        } while (node != i);
    }
    return potential;
}

} // namespace physical_layers::detail
