#pragma once

#include <functional>
#include <vector>

class planet;

namespace terraindetail
{
// Eight-neighbour land distance, wrapping longitude and stopping at ocean cells.
bool buildlanddistancefield(planet& world, std::vector<std::vector<short>>& distances,
    int maximumdistance, const std::function<bool(int, int)>& issource);
}
